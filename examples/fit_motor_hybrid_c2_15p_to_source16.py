#!/usr/bin/env python3
"""Fit a 15-parameter structured C2 motion law to the 16-control 260 K champion.

NO thermodynamic integration is performed.

The goal is to test representation capacity before launching any expensive
thermodynamic optimization.

Structure
---------
Shared timing / extrema:
  * 3 parameters:
      small_max_deg
      small_down_duration_deg
      large_down_duration_deg

Extremum sharpness:
  * 4 positive curvature magnitudes:
      small_max_curvature
      small_min_curvature
      large_max_curvature
      large_min_curvature
  Velocity is zero at extrema, but acceleration is NOT forced to zero.
  The same extremum curvature is used by the two adjacent branches, preserving
  global C2 continuity.

BP-like exchange branches:
  * large maximum -> minimum
  * small minimum -> maximum
  Each branch has one midpoint-shape parameter q_mid.
  q_mid = 0.5 is the neutral middle shape; values above/below let the branch
  bow to either side.  This removes the unjustified "strictly linear BP" prior
  while structurally forbidding high-frequency wiggles.

Opposite / HP-side branches:
  * small maximum -> minimum
  * large minimum -> maximum
  Each has:
      kink_u
      kink_q
      kink_width_rel
  The kink is represented by a narrow C2 transition around (u, q).  Reducing
  kink_width_rel makes the change of slope progressively harder while the
  complete motion remains C2 and monotone.

Total: 3 + 4 + 2 + 6 = 15 parameters.

Outputs
-------
  outputs/motor_hybrid_c2_15p_fit_260k/fit.json
  outputs/motor_hybrid_c2_15p_fit_260k/motion.csv
  outputs/motor_hybrid_c2_15p_fit_260k/comparison.png
"""

from __future__ import annotations

from pathlib import Path
import argparse
import json
import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import BPoly
from scipy.optimize import least_squares

from compare_motor_motion_laws_stage7A5 import ROOT
from optimize_motor_temperature_point import base_geometry
from refine_motor_four_stage_hx9d_variable_gas import _load_basis
from refine_motor_fourier_c2_260k import build_design as fourier_build_design
from optimize_motor_free_spline_260k_v3 import (
    _make_kinematics,
    _normalized_motion,
    _source_best,
    _source_harmonics,
)

DEFAULT_SOURCE = ROOT / "outputs/motor_spline_from_linear_260k/report.json"
DEFAULT_FOURIER = ROOT / "outputs/motor_fourier_c2_8h_refine/report.json"
DEFAULT_OUTPUT = ROOT / "outputs/motor_hybrid_c2_15p_fit_260k"

MIN_BRANCH_DEG = 35.0
MAX_BRANCH_DEG = 325.0

CURVATURE_MIN = 0.1
CURVATURE_MAX = 4000.0

BP_MID_MIN = 0.20
BP_MID_MAX = 0.80

KINK_MIN = 0.035
KINK_MAX = 0.965
KINK_WIDTH_REL_MIN = 0.015
KINK_WIDTH_REL_MAX = 0.90

FIT_SAMPLES = 1440
VALIDATION_SAMPLES = 5760

# The position remains dominant.  Velocity is strongly represented; acceleration
# receives only a light weight because source16 acceleration contains local
# spline structure we do not necessarily want to reproduce literally.
VELOCITY_FIT_WEIGHT = 0.06
ACCELERATION_FIT_WEIGHT = 0.0015
MONOTONICITY_WEIGHT = 30.0

MULTISTART_WIDTHS = (0.04, 0.10, 0.22, 0.45)


def resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def _harmonic(a, b):
    a = max(float(a), 1e-12)
    b = max(float(b), 1e-12)
    return 2.0*a*b/(a+b)


def _bpoly_from_nodes(nodes):
    """nodes = [(u, q, dq/du, d2q/du2), ...]."""
    x = [float(n[0]) for n in nodes]
    yi = [
        [float(n[1]), float(n[2]), float(n[3])]
        for n in nodes
    ]
    return BPoly.from_derivatives(x, yi)


def _bp_branch(q_mid, acc_start, acc_end):
    """Increasing 0->1 branch with one global bowing parameter."""
    q_mid = float(q_mid)

    sec_left = q_mid / 0.5
    sec_right = (1.0-q_mid) / 0.5
    slope_mid = _harmonic(sec_left, sec_right)

    return _bpoly_from_nodes([
        (0.0, 0.0, 0.0, float(acc_start)),
        (0.5, q_mid, slope_mid, 0.0),
        (1.0, 1.0, 0.0, float(acc_end)),
    ])


def _kink_branch(kink_u, kink_q, width_rel, acc_start, acc_end):
    """Increasing 0->1 branch with a tunably hard internal C2 turn."""
    ku = float(kink_u)
    kq = float(kink_q)
    wr = float(width_rel)

    # Full transition width scales with the available distance to the nearer
    # endpoint.  Hence wr can approach zero without crossing an endpoint.
    full_w = wr * min(ku, 1.0-ku)
    half = 0.5*full_w
    ul = ku-half
    ur = ku+half

    # "Outer" slopes implied by the two sides of the kink.
    m1 = kq/ku
    m2 = (1.0-kq)/(1.0-ku)
    mk = _harmonic(m1, m2)

    # Side values lie on the two asymptotic secants, while the exact kink point
    # is retained at (ku,kq).  As width -> 0, the C2 turn becomes arbitrarily
    # concentrated around the kink.
    ql = kq - m1*half
    qr = kq + m2*half

    return _bpoly_from_nodes([
        (0.0, 0.0, 0.0, float(acc_start)),
        (ul, ql, m1, 0.0),
        (ku, kq, mk, 0.0),
        (ur, qr, m2, 0.0),
        (1.0, 1.0, 0.0, float(acc_end)),
    ])


class StructuredMotion15:
    def __init__(self, p):
        self.p = dict(p)

        sd = self.p["small_down_duration_deg"]/360.0
        su = 1.0-sd
        ld = self.p["large_down_duration_deg"]/360.0
        lu = 1.0-ld

        # Global normalized motor-time extrema curvatures:
        # q''(max) = -Amax ; q''(min) = +Amin.
        smax = self.p["small_max_curvature"]
        smin = self.p["small_min_curvature"]
        lmax = self.p["large_max_curvature"]
        lmin = self.p["large_min_curvature"]

        # Small down (non-BP), physical q: 1 -> 0.
        # Helper y=1-q increases 0->1.
        self.small_down = _kink_branch(
            self.p["small_down_kink_u"],
            self.p["small_down_kink_q"],
            self.p["small_down_kink_width_rel"],
            +smax*sd*sd,
            -smin*sd*sd,
        )

        # Small up (BP), physical q increases 0->1.
        self.small_up = _bp_branch(
            self.p["small_up_bp_mid_q"],
            +smin*su*su,
            -smax*su*su,
        )

        # Large down (BP), physical q: 1 -> 0, helper y=1-q.
        self.large_down = _bp_branch(
            self.p["large_down_bp_mid_q"],
            +lmax*ld*ld,
            -lmin*ld*ld,
        )

        # Large up (non-BP), physical q increases 0->1.
        self.large_up = _kink_branch(
            self.p["large_up_kink_u"],
            self.p["large_up_kink_q"],
            self.p["large_up_kink_width_rel"],
            +lmin*lu*lu,
            -lmax*lu*lu,
        )

    @staticmethod
    def _eval(poly, u):
        u = np.asarray(u, dtype=float)
        return (
            np.asarray(poly(u), dtype=float),
            np.asarray(poly.derivative(1)(u), dtype=float),
            np.asarray(poly.derivative(2)(u), dtype=float),
        )

    def normalized(self, t):
        t = np.asarray(t, dtype=float) % 1.0
        p = self.p

        # Small piston.
        smax_t = (p["small_max_deg"] % 360.0)/360.0
        sd = p["small_down_duration_deg"]/360.0
        from_smax = np.mod(t-smax_t, 1.0)
        sdown = from_smax <= sd

        sq = np.empty_like(t)
        sv = np.empty_like(t)
        sa = np.empty_like(t)

        if np.any(sdown):
            u = from_smax[sdown]/sd
            y, dy, ddy = self._eval(self.small_down, u)
            sq[sdown] = 1.0-y
            sv[sdown] = -dy/sd
            sa[sdown] = -ddy/(sd*sd)

        if np.any(~sdown):
            su = 1.0-sd
            u = (from_smax[~sdown]-sd)/su
            y, dy, ddy = self._eval(self.small_up, u)
            sq[~sdown] = y
            sv[~sdown] = dy/su
            sa[~sdown] = ddy/(su*su)

        # Large piston. Gauge: large maximum at t=0.
        ld = p["large_down_duration_deg"]/360.0
        ldown = t <= ld

        lq = np.empty_like(t)
        lv = np.empty_like(t)
        la = np.empty_like(t)

        if np.any(ldown):
            u = t[ldown]/ld
            y, dy, ddy = self._eval(self.large_down, u)
            lq[ldown] = 1.0-y
            lv[ldown] = -dy/ld
            la[ldown] = -ddy/(ld*ld)

        if np.any(~ldown):
            lu = 1.0-ld
            u = (t[~ldown]-ld)/lu
            y, dy, ddy = self._eval(self.large_up, u)
            lq[~ldown] = y
            lv[~ldown] = dy/lu
            la[~ldown] = ddy/(lu*lu)

        return (
            np.vstack((sq, lq)),
            np.vstack((sv, lv)),
            np.vstack((sa, la)),
        )

    def monotonicity_penalty(self, n=160):
        u = np.linspace(0.0, 1.0, n)
        derivatives = [
            self.small_down.derivative(1)(u),
            self.small_up.derivative(1)(u),
            self.large_down.derivative(1)(u),
            self.large_up.derivative(1)(u),
        ]
        return np.concatenate([
            np.minimum(np.asarray(v, dtype=float), 0.0)
            for v in derivatives
        ])


def _source_kinematics(best, limits):
    return _make_kinematics(
        limits,
        np.asarray(best["small_controls"], dtype=float),
        np.asarray(best["large_controls"], dtype=float),
        float(best["small_phase_deg"]),
        float(best["large_phase_deg"]),
    )


def _source_aligned(source, large_max_deg, t):
    t = np.asarray(t, dtype=float)
    tref = np.mod(t + large_max_deg/360.0, 1.0)
    q, dq = _normalized_motion(source, tref)

    theta = -2.0*math.pi*tref
    sl = source.small_volume_limits
    ll = source.large_volume_limits

    def evaluate_second(name):
        fn = getattr(source, name)
        try:
            values = np.asarray(fn(theta), dtype=float)
            if values.shape == theta.shape:
                return values
        except (TypeError, ValueError):
            pass
        return np.asarray([fn(float(x)) for x in theta], dtype=float)

    s2 = evaluate_second("small_cylinder_volume_second_derivative")
    l2 = evaluate_second("large_cylinder_volume_second_derivative")

    # d2q / d(motor_fraction)^2
    ddq = np.vstack((
        4.0*math.pi*math.pi*s2/sl.swept,
        4.0*math.pi*math.pi*l2/ll.swept,
    ))
    return q, dq, ddq


def _aligned_source_extrema(best):
    lmax = float(best["shape_diagnostics"]["large_max_deg"])
    return {
        "small_max": (
            float(best["shape_diagnostics"]["small_max_deg"])-lmax
        ) % 360.0,
        "small_min": (
            float(best["shape_diagnostics"]["small_min_deg"])-lmax
        ) % 360.0,
        "large_max": 0.0,
        "large_min": (
            float(best["shape_diagnostics"]["large_min_deg"])-lmax
        ) % 360.0,
    }


def _curvature_seed(source, best, extrema):
    lmax = float(best["shape_diagnostics"]["large_max_deg"])
    out = {}
    for name, deg in extrema.items():
        _, _, ddq = _source_aligned(
            source, lmax, np.array([deg/360.0], dtype=float)
        )
        row = 0 if name.startswith("small") else 1
        out[name] = max(abs(float(ddq[row, 0])), CURVATURE_MIN)
    return out


def _unpack(x):
    c = np.exp(np.asarray(x[3:7], dtype=float))
    return {
        "small_max_deg": float(x[0]) % 360.0,
        "small_down_duration_deg": float(x[1]),
        "large_down_duration_deg": float(x[2]),

        "small_max_curvature": float(c[0]),
        "small_min_curvature": float(c[1]),
        "large_max_curvature": float(c[2]),
        "large_min_curvature": float(c[3]),

        "small_up_bp_mid_q": float(x[7]),
        "large_down_bp_mid_q": float(x[8]),

        "small_down_kink_u": float(x[9]),
        "small_down_kink_q": float(x[10]),
        "small_down_kink_width_rel": float(x[11]),

        "large_up_kink_u": float(x[12]),
        "large_up_kink_q": float(x[13]),
        "large_up_kink_width_rel": float(x[14]),
    }


def _metrics(ref_q, ref_dq, ref_ddq, q, dq, ddq):
    return {
        "small_position_rms":
            float(np.sqrt(np.mean((q[0]-ref_q[0])**2))),
        "large_position_rms":
            float(np.sqrt(np.mean((q[1]-ref_q[1])**2))),
        "combined_position_rms":
            float(np.sqrt(np.mean((q-ref_q)**2))),
        "small_velocity_rms":
            float(np.sqrt(np.mean((dq[0]-ref_dq[0])**2))),
        "large_velocity_rms":
            float(np.sqrt(np.mean((dq[1]-ref_dq[1])**2))),
        "combined_velocity_rms":
            float(np.sqrt(np.mean((dq-ref_dq)**2))),
        "small_acceleration_rms":
            float(np.sqrt(np.mean((ddq[0]-ref_ddq[0])**2))),
        "large_acceleration_rms":
            float(np.sqrt(np.mean((ddq[1]-ref_ddq[1])**2))),
        "combined_acceleration_rms":
            float(np.sqrt(np.mean((ddq-ref_ddq)**2))),
        "maximum_position_difference":
            float(np.max(np.abs(q-ref_q))),
    }


def _extrema_count(d):
    x = np.asarray(d, dtype=float)
    threshold = max(1e-8, 1e-4*float(np.max(np.abs(x))))
    s = np.sign(x)
    s[np.abs(x) < threshold] = 0.0
    s = s[s != 0.0]
    if len(s) == 0:
        return 0
    return int(np.sum(s[1:] != s[:-1]) + (s[0] != s[-1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-report", type=Path, default=DEFAULT_SOURCE)
    ap.add_argument("--fourier-report", type=Path, default=DEFAULT_FOURIER)
    ap.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT)
    args = ap.parse_args()

    sp = resolve(args.source_report)
    fp = resolve(args.fourier_report)
    out = resolve(args.output_directory)
    out.mkdir(parents=True, exist_ok=True)

    sr = json.loads(sp.read_text())
    sb = sr["best_feasible"]

    fr = json.loads(fp.read_text())
    fb = _source_best(fr)
    h = _source_harmonics(fb)

    _, base = _load_basis()
    geometry = base_geometry(base)
    machine = fourier_build_design(
        base,
        geometry,
        fb["thermo"],
        np.asarray(fb["small_coefficients"]),
        np.asarray(fb["large_coefficients"]),
        h,
    )
    limits = machine.configuration.machine_volumes
    source = _source_kinematics(sb, limits)

    source_lmax = float(sb["shape_diagnostics"]["large_max_deg"])
    extrema = _aligned_source_extrema(sb)
    curv0 = _curvature_seed(source, sb, extrema)

    small_down = (extrema["small_min"]-extrema["small_max"]) % 360.0
    large_down = extrema["large_min"]

    lo = np.array([
        0.0,
        MIN_BRANCH_DEG,
        MIN_BRANCH_DEG,

        math.log(CURVATURE_MIN),
        math.log(CURVATURE_MIN),
        math.log(CURVATURE_MIN),
        math.log(CURVATURE_MIN),

        BP_MID_MIN,
        BP_MID_MIN,

        KINK_MIN,
        KINK_MIN,
        KINK_WIDTH_REL_MIN,

        KINK_MIN,
        KINK_MIN,
        KINK_WIDTH_REL_MIN,
    ], dtype=float)

    hi = np.array([
        360.0,
        MAX_BRANCH_DEG,
        MAX_BRANCH_DEG,

        math.log(CURVATURE_MAX),
        math.log(CURVATURE_MAX),
        math.log(CURVATURE_MAX),
        math.log(CURVATURE_MAX),

        BP_MID_MAX,
        BP_MID_MAX,

        KINK_MAX,
        KINK_MAX,
        KINK_WIDTH_REL_MAX,

        KINK_MAX,
        KINK_MAX,
        KINK_WIDTH_REL_MAX,
    ], dtype=float)

    tfit = np.linspace(0.0, 1.0, FIT_SAMPLES, endpoint=False)
    ref_q, ref_dq, ref_ddq = _source_aligned(source, source_lmax, tfit)

    def residual(x):
        try:
            motion = StructuredMotion15(_unpack(x))
            q, dq, ddq = motion.normalized(tfit)
        except (ValueError, FloatingPointError):
            return np.full(6*FIT_SAMPLES + 640, 1e3)

        mono = motion.monotonicity_penalty()
        return np.concatenate((
            (q-ref_q).ravel(),
            VELOCITY_FIT_WEIGHT*(dq-ref_dq).ravel(),
            ACCELERATION_FIT_WEIGHT*(ddq-ref_ddq).ravel(),
            MONOTONICITY_WEIGHT*mono,
        ))

    starts = []
    for width0 in MULTISTART_WIDTHS:
        starts.append(np.array([
            extrema["small_max"],
            small_down,
            large_down,

            math.log(curv0["small_max"]),
            math.log(curv0["small_min"]),
            math.log(curv0["large_max"]),
            math.log(curv0["large_min"]),

            0.50,
            0.50,

            0.33,
            0.59,
            width0,

            0.74,
            0.42,
            width0,
        ], dtype=float))

    # Two extra starts bias the BP branches slightly to opposite sides.
    for bp_s, bp_l in ((0.45, 0.55), (0.55, 0.45)):
        x = starts[1].copy()
        x[7] = bp_s
        x[8] = bp_l
        starts.append(x)

    best_fit = None
    best_x = None

    for i, x0 in enumerate(starts):
        fit = least_squares(
            residual,
            x0,
            bounds=(lo, hi),
            x_scale="jac",
            xtol=1e-11,
            ftol=1e-11,
            gtol=1e-11,
            max_nfev=5000,
            verbose=0,
        )
        print(json.dumps({
            "start": i,
            "success": bool(fit.success),
            "cost": float(fit.cost),
            "nfev": int(fit.nfev),
            "small_kink_width_rel": float(fit.x[11]),
            "large_kink_width_rel": float(fit.x[14]),
        }), flush=True)

        if best_fit is None or fit.cost < best_fit.cost:
            best_fit = fit
            best_x = fit.x.copy()

    params = _unpack(best_x)
    motion = StructuredMotion15(params)

    t = np.linspace(0.0, 1.0, VALIDATION_SAMPLES+1)
    ref_q, ref_dq, ref_ddq = _source_aligned(source, source_lmax, t)
    q, dq, ddq = motion.normalized(t)

    metrics = _metrics(ref_q, ref_dq, ref_ddq, q, dq, ddq)
    metrics["small_extrema_count"] = _extrema_count(dq[0])
    metrics["large_extrema_count"] = _extrema_count(dq[1])
    metrics["minimum_helper_branch_derivative"] = float(
        min(
            np.min(motion.small_down.derivative(1)(np.linspace(0,1,1000))),
            np.min(motion.small_up.derivative(1)(np.linspace(0,1,1000))),
            np.min(motion.large_down.derivative(1)(np.linspace(0,1,1000))),
            np.min(motion.large_up.derivative(1)(np.linspace(0,1,1000))),
        )
    )

    # Actual fitted accelerations at extrema, in normalized motor-time units.
    smax_t = params["small_max_deg"]/360.0
    smin_t = (
        params["small_max_deg"] + params["small_down_duration_deg"]
    ) % 360.0 / 360.0
    lmax_t = 0.0
    lmin_t = params["large_down_duration_deg"]/360.0

    extrema_t = {
        "small_max": smax_t,
        "small_min": smin_t,
        "large_max": lmax_t,
        "large_min": lmin_t,
    }
    fitted_extremum_acc = {}
    for name, et in extrema_t.items():
        _, _, aa = motion.normalized(np.array([et]))
        row = 0 if name.startswith("small") else 1
        fitted_extremum_acc[name] = float(aa[row, 0])

    result = {
        "model":
            "15-parameter structured C2 fit with BP bowing and hard-width HP turns",
        "thermodynamic_evaluations": 0,
        "success": bool(best_fit.success),
        "message": best_fit.message,
        "nfev": int(best_fit.nfev),
        "cost": float(best_fit.cost),
        "fit_weights": {
            "position": 1.0,
            "velocity": VELOCITY_FIT_WEIGHT,
            "acceleration": ACCELERATION_FIT_WEIGHT,
            "monotonicity": MONOTONICITY_WEIGHT,
        },
        "source_extrema_aligned_deg": extrema,
        "source_extremum_curvature_seed_abs_d2q_dt2": curv0,
        "parameters": params,
        "fitted_extremum_acceleration_d2q_dt2": fitted_extremum_acc,
        "metrics": metrics,
    }

    (out/"fit.json").write_text(json.dumps(result, indent=2)+"\n")

    np.savetxt(
        out/"motion.csv",
        np.column_stack((
            t, 360.0*t,
            ref_q[0], ref_q[1], q[0], q[1],
            ref_dq[0], ref_dq[1], dq[0], dq[1],
            ref_ddq[0], ref_ddq[1], ddq[0], ddq[1],
        )),
        delimiter=",",
        header=(
            "motor_time_fraction,motor_angle_deg,"
            "source_small,source_large,fit_small,fit_large,"
            "source_dsmall_dt,source_dlarge_dt,fit_dsmall_dt,fit_dlarge_dt,"
            "source_ddsmall_dt2,source_ddlarge_dt2,"
            "fit_ddsmall_dt2,fit_ddlarge_dt2"
        ),
        comments="",
    )

    fig, axes = plt.subplots(3, 2, figsize=(13, 11), sharex=True)
    labels = ("Small piston", "Large piston")

    for j in range(2):
        ax = axes[0, j]
        ax.plot(360*t, ref_q[j], "--", label="Source16 aligned")
        ax.plot(360*t, q[j], label="15p structured C2")
        ax.set_title(f"{labels[j]} — normalized position")
        ax.set_ylabel("q")
        ax.grid(True, alpha=0.25)
        ax.legend()

        ax = axes[1, j]
        ax.plot(360*t, ref_dq[j], "--", label="Source16 aligned")
        ax.plot(360*t, dq[j], label="15p structured C2")
        ax.axhline(0.0, linewidth=0.8)
        ax.set_title(f"{labels[j]} — normalized velocity")
        ax.set_ylabel("dq/dt")
        ax.grid(True, alpha=0.25)
        ax.legend()

        ax = axes[2, j]
        ax.plot(360*t, ref_ddq[j], "--", label="Source16 aligned")
        ax.plot(360*t, ddq[j], label="15p structured C2")
        ax.axhline(0.0, linewidth=0.8)
        ax.set_title(f"{labels[j]} — normalized acceleration")
        ax.set_ylabel("d²q/dt²")
        ax.set_xlabel("Motor angle [deg]")
        ax.grid(True, alpha=0.25)
        ax.legend()

    for ax in axes.ravel():
        ax.set_xlim(0, 360)
        ax.set_xticks(np.arange(0, 361, 45))

    fig.suptitle(
        "15-parameter structured C2 fit to source16 — no thermodynamics\n"
        f"position RMS={metrics['combined_position_rms']:.5f}   "
        f"velocity RMS={metrics['combined_velocity_rms']:.5f}   "
        f"acceleration RMS={metrics['combined_acceleration_rms']:.3f}\n"
        f"HP widths: S↓={params['small_down_kink_width_rel']:.4f}   "
        f"L↑={params['large_up_kink_width_rel']:.4f}   |   "
        f"BP mid-q: S↑={params['small_up_bp_mid_q']:.4f}   "
        f"L↓={params['large_down_bp_mid_q']:.4f}"
    )
    fig.tight_layout(rect=(0,0,1,0.93))
    fig.savefig(out/"comparison.png", dpi=180)

    print(json.dumps(result, indent=2))
    print(f"Saved {out/'fit.json'}")
    print(f"Saved {out/'motion.csv'}")
    print(f"Saved {out/'comparison.png'}")


if __name__ == "__main__":
    main()
