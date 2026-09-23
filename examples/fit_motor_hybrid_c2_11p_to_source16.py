#!/usr/bin/env python3
"""Fit an 11-parameter C2 compact motion law to the 16-control 260 K champion.

This script performs NO thermodynamic integration.

Purpose
-------
Before spending solver time on a new optimization campaign, test whether the
proposed compact family can reproduce the source16 kinematics in:
  * normalized position,
  * normalized velocity,
  * normalized acceleration.

Compared with the previous 9-parameter hybrid law, extrema no longer impose
zero acceleration.  Instead, each piston has one shared curvature magnitude at
its maximum and one at its minimum.  Those curvatures are continuous across the
two branches meeting at the extremum, so the complete periodic motion remains
C2 while velocity crosses zero with a finite slope.

11 parameters for BOTH pistons:
  3 extrema/timing parameters
  4 positive extremum-curvature magnitudes
  4 coordinates for the two non-BP interior kink points

BP branches:
  * large max -> min
  * small min -> max
are each represented by one quintic Hermite branch.  Their endpoint velocity is
zero and endpoint acceleration is supplied by the shared extremum curvatures.

Opposite branches:
  * small max -> min
  * large min -> max
use two quintic Hermite segments through one mobile (u_k, q_k).  The internal
knot uses a monotone harmonic-mean slope and zero acceleration.

Outputs
-------
  outputs/motor_hybrid_c2_11p_fit_260k/fit.json
  outputs/motor_hybrid_c2_11p_fit_260k/motion.csv
  outputs/motor_hybrid_c2_11p_fit_260k/comparison.png
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import argparse
import json
import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

from dada_solver.geometry import CylinderVolumeLimits

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
DEFAULT_OUTPUT = ROOT / "outputs/motor_hybrid_c2_11p_fit_260k"

MIN_BRANCH_DEG = 35.0
MAX_BRANCH_DEG = 325.0
KINK_MIN = 0.035
KINK_MAX = 0.965

CURVATURE_MIN = 0.25
CURVATURE_MAX = 2500.0

VELOCITY_FIT_WEIGHT = 0.05
MONOTONICITY_WEIGHT = 20.0
FIT_SAMPLES = 1440


def resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def _quintic_coefficients(y0, y1, width, slope0, slope1, acc0, acc1):
    """Quintic on local z=[0,1], derivatives specified wrt parent coordinate."""
    h = float(width)
    c0 = float(y0)
    c1 = h * float(slope0)
    c2 = 0.5 * h * h * float(acc0)

    r0 = float(y1) - c0 - c1 - c2
    r1 = h * float(slope1) - c1 - 2.0*c2
    r2 = h*h * float(acc1) - 2.0*c2

    c3 = 10.0*r0 - 4.0*r1 + 0.5*r2
    c4 = -15.0*r0 + 7.0*r1 - r2
    c5 = 6.0*r0 - 3.0*r1 + 0.5*r2
    return (c0, c1, c2, c3, c4, c5)


def _poly_eval(coeffs, z, width):
    c0, c1, c2, c3, c4, c5 = coeffs
    z = np.asarray(z, dtype=float)

    y = ((((c5*z + c4)*z + c3)*z + c2)*z + c1)*z + c0
    dy_dz = (((5*c5*z + 4*c4)*z + 3*c3)*z + 2*c2)*z + c1
    d2y_dz2 = ((20*c5*z + 12*c4)*z + 6*c3)*z + 2*c2

    h = float(width)
    return y, dy_dz/h, d2y_dz2/(h*h)


def _harmonic_slope(kink_u, kink_q):
    d1 = float(kink_q) / float(kink_u)
    d2 = (1.0-float(kink_q)) / (1.0-float(kink_u))
    return 2.0*d1*d2/(d1+d2)


@dataclass(frozen=True, slots=True)
class CompactC2FitKinematics:
    small_limits: CylinderVolumeLimits
    large_limits: CylinderVolumeLimits

    small_max_deg: float
    small_down_duration_deg: float
    large_down_duration_deg: float

    small_max_curvature: float
    small_min_curvature: float
    large_max_curvature: float
    large_min_curvature: float

    small_down_kink_u: float
    small_down_kink_q: float
    large_up_kink_u: float
    large_up_kink_q: float

    _small_down_left: tuple = field(init=False, repr=False)
    _small_down_right: tuple = field(init=False, repr=False)
    _small_up: tuple = field(init=False, repr=False)
    _large_down: tuple = field(init=False, repr=False)
    _large_up_left: tuple = field(init=False, repr=False)
    _large_up_right: tuple = field(init=False, repr=False)

    def __post_init__(self):
        for x in (self.small_down_duration_deg, self.large_down_duration_deg):
            if not MIN_BRANCH_DEG <= x <= MAX_BRANCH_DEG:
                raise ValueError("Branch duration outside admissible range.")

        for x in (
            self.small_max_curvature,
            self.small_min_curvature,
            self.large_max_curvature,
            self.large_min_curvature,
        ):
            if not CURVATURE_MIN <= x <= CURVATURE_MAX:
                raise ValueError("Extremum curvature outside admissible range.")

        for x in (
            self.small_down_kink_u,
            self.small_down_kink_q,
            self.large_up_kink_u,
            self.large_up_kink_q,
        ):
            if not KINK_MIN <= x <= KINK_MAX:
                raise ValueError("Kink coordinate outside admissible range.")

        # Durations expressed in motor-cycle fraction.
        sd = self.small_down_duration_deg / 360.0
        su = 1.0 - sd
        ld = self.large_down_duration_deg / 360.0
        lu = 1.0 - ld

        # All helper branches below increase y: 0 -> 1.
        # For a decreasing physical q branch, y = 1-q.
        #
        # Global q_tt curvature signs:
        #   maximum: -Amax
        #   minimum: +Amin
        #
        # Conversion to local branch-u curvature multiplies by duration^2.

        # Small down: q 1->0, y=1-q: max -> min.
        sdu0 = self.small_max_curvature * sd*sd
        sdu1 = -self.small_min_curvature * sd*sd
        sku = self.small_down_kink_u
        skq = self.small_down_kink_q
        sks = _harmonic_slope(sku, skq)
        left = _quintic_coefficients(
            0.0, skq, sku,
            0.0, sks,
            sdu0, 0.0,
        )
        right = _quintic_coefficients(
            skq, 1.0, 1.0-sku,
            sks, 0.0,
            0.0, sdu1,
        )
        object.__setattr__(self, "_small_down_left", left)
        object.__setattr__(self, "_small_down_right", right)

        # Small up BP: q 0->1: min -> max.
        sup0 = self.small_min_curvature * su*su
        sup1 = -self.small_max_curvature * su*su
        object.__setattr__(
            self,
            "_small_up",
            _quintic_coefficients(
                0.0, 1.0, 1.0,
                0.0, 0.0,
                sup0, sup1,
            ),
        )

        # Large down BP: q 1->0, y=1-q: max -> min.
        ldu0 = self.large_max_curvature * ld*ld
        ldu1 = -self.large_min_curvature * ld*ld
        object.__setattr__(
            self,
            "_large_down",
            _quintic_coefficients(
                0.0, 1.0, 1.0,
                0.0, 0.0,
                ldu0, ldu1,
            ),
        )

        # Large up: q 0->1: min -> max.
        luu0 = self.large_min_curvature * lu*lu
        luu1 = -self.large_max_curvature * lu*lu
        lku = self.large_up_kink_u
        lkq = self.large_up_kink_q
        lks = _harmonic_slope(lku, lkq)
        left = _quintic_coefficients(
            0.0, lkq, lku,
            0.0, lks,
            luu0, 0.0,
        )
        right = _quintic_coefficients(
            lkq, 1.0, 1.0-lku,
            lks, 0.0,
            0.0, luu1,
        )
        object.__setattr__(self, "_large_up_left", left)
        object.__setattr__(self, "_large_up_right", right)

    @property
    def small_volume_limits(self):
        return self.small_limits

    @property
    def large_volume_limits(self):
        return self.large_limits

    def _increasing_kink(self, u, ku, left, right):
        a = np.asarray(u, dtype=float)
        q = np.empty_like(a)
        dq = np.empty_like(a)
        ddq = np.empty_like(a)

        mask = a <= ku
        if np.any(mask):
            y, dy, ddy = _poly_eval(left, a[mask]/ku, ku)
            q[mask], dq[mask], ddq[mask] = y, dy, ddy
        if np.any(~mask):
            h = 1.0-ku
            y, dy, ddy = _poly_eval(right, (a[~mask]-ku)/h, h)
            q[~mask], dq[~mask], ddq[~mask] = y, dy, ddy

        return q, dq, ddq

    def _small_motor(self, motor_fraction):
        t = np.asarray(motor_fraction, dtype=float) % 1.0
        max_t = (self.small_max_deg % 360.0) / 360.0
        from_max = np.mod(t-max_t, 1.0)

        down_d = self.small_down_duration_deg/360.0
        on_down = from_max <= down_d

        q = np.empty_like(t)
        dqdt = np.empty_like(t)
        ddqdt2 = np.empty_like(t)

        if np.any(on_down):
            u = from_max[on_down]/down_d
            y, dy, ddy = self._increasing_kink(
                u,
                self.small_down_kink_u,
                self._small_down_left,
                self._small_down_right,
            )
            q[on_down] = 1.0-y
            dqdt[on_down] = -dy/down_d
            ddqdt2[on_down] = -ddy/(down_d*down_d)

        if np.any(~on_down):
            up_d = 1.0-down_d
            u = (from_max[~on_down]-down_d)/up_d
            y, dy, ddy = _poly_eval(self._small_up, u, 1.0)
            q[~on_down] = y
            dqdt[~on_down] = dy/up_d
            ddqdt2[~on_down] = ddy/(up_d*up_d)

        return q, dqdt, ddqdt2

    def _large_motor(self, motor_fraction):
        t = np.asarray(motor_fraction, dtype=float) % 1.0

        down_d = self.large_down_duration_deg/360.0
        on_down = t <= down_d

        q = np.empty_like(t)
        dqdt = np.empty_like(t)
        ddqdt2 = np.empty_like(t)

        if np.any(on_down):
            u = t[on_down]/down_d
            y, dy, ddy = _poly_eval(self._large_down, u, 1.0)
            q[on_down] = 1.0-y
            dqdt[on_down] = -dy/down_d
            ddqdt2[on_down] = -ddy/(down_d*down_d)

        if np.any(~on_down):
            up_d = 1.0-down_d
            u = (t[~on_down]-down_d)/up_d
            y, dy, ddy = self._increasing_kink(
                u,
                self.large_up_kink_u,
                self._large_up_left,
                self._large_up_right,
            )
            q[~on_down] = y
            dqdt[~on_down] = dy/up_d
            ddqdt2[~on_down] = ddy/(up_d*up_d)

        return q, dqdt, ddqdt2

    def normalized_motor(self, motor_fraction):
        sq, sd, sdd = self._small_motor(motor_fraction)
        lq, ld, ldd = self._large_motor(motor_fraction)
        return np.vstack((sq, lq)), np.vstack((sd, ld)), np.vstack((sdd, ldd))

    def monotonicity_penalty(self, samples=96):
        # Every helper branch is represented as an increasing 0->1 law.
        u = np.linspace(0.0, 1.0, samples)

        _, sd1, _ = self._increasing_kink(
            u, self.small_down_kink_u,
            self._small_down_left, self._small_down_right,
        )
        _, su1, _ = _poly_eval(self._small_up, u, 1.0)
        _, ld1, _ = _poly_eval(self._large_down, u, 1.0)
        _, lu1, _ = self._increasing_kink(
            u, self.large_up_kink_u,
            self._large_up_left, self._large_up_right,
        )
        return np.concatenate((
            np.minimum(sd1, 0.0),
            np.minimum(su1, 0.0),
            np.minimum(ld1, 0.0),
            np.minimum(lu1, 0.0),
        ))

    def _pair(self, theta, small):
        scalar = np.ndim(theta) == 0
        a = np.asarray(theta, dtype=float)
        t = np.mod(-a/(2.0*math.pi), 1.0)
        if small:
            q, dqdt, _ = self._small_motor(t)
            lim = self.small_limits
        else:
            q, dqdt, _ = self._large_motor(t)
            lim = self.large_limits

        v = lim.minimum + q*lim.swept
        dvdtheta = -(dqdt/(2.0*math.pi))*lim.swept
        if scalar:
            return float(v), float(dvdtheta)
        return v, dvdtheta

    def _second(self, theta, small):
        scalar = np.ndim(theta) == 0
        a = np.asarray(theta, dtype=float)
        t = np.mod(-a/(2.0*math.pi), 1.0)
        if small:
            _, _, ddqdt2 = self._small_motor(t)
            lim = self.small_limits
        else:
            _, _, ddqdt2 = self._large_motor(t)
            lim = self.large_limits
        out = ddqdt2*lim.swept/(4.0*math.pi*math.pi)
        return float(out) if scalar else out

    def small_cylinder_volume(self, theta):
        return self._pair(theta, True)[0]

    def large_cylinder_volume(self, theta):
        return self._pair(theta, False)[0]

    def small_cylinder_volume_derivative(self, theta):
        return self._pair(theta, True)[1]

    def large_cylinder_volume_derivative(self, theta):
        return self._pair(theta, False)[1]

    def small_cylinder_volume_second_derivative(self, theta):
        return self._second(theta, True)

    def large_cylinder_volume_second_derivative(self, theta):
        return self._second(theta, False)


def _source_kinematics(best, limits):
    return _make_kinematics(
        limits,
        np.asarray(best["small_controls"], dtype=float),
        np.asarray(best["large_controls"], dtype=float),
        float(best["small_phase_deg"]),
        float(best["large_phase_deg"]),
    )


def _source_aligned(source, large_max_deg, motor_fraction):
    t = np.asarray(motor_fraction, dtype=float)
    tref = np.mod(t + large_max_deg/360.0, 1.0)
    q, dq = _normalized_motion(source, tref)

    theta = -2.0*math.pi*tref
    sl = source.small_volume_limits
    ll = source.large_volume_limits

    def eval_second(name):
        fn = getattr(source, name)
        try:
            values = np.asarray(fn(theta), dtype=float)
            if values.shape == theta.shape:
                return values
        except (TypeError, ValueError):
            pass
        return np.asarray([fn(float(x)) for x in theta], dtype=float)

    s2 = eval_second("small_cylinder_volume_second_derivative")
    l2 = eval_second("large_cylinder_volume_second_derivative")
    ddq = np.vstack((
        4.0*math.pi*math.pi*s2/sl.swept,
        4.0*math.pi*math.pi*l2/ll.swept,
    ))
    return q, dq, ddq


def _extremum_curvature_estimates(source, best):
    lmax = float(best["shape_diagnostics"]["large_max_deg"])

    aligned = {
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

    values = {}
    for name, deg in aligned.items():
        _, _, a = _source_aligned(source, lmax, np.array([deg/360.0]))
        row = 0 if name.startswith("small") else 1
        values[name] = abs(float(a[row, 0]))

    return aligned, values


def _unpack(x):
    curv = np.exp(np.asarray(x[3:7], dtype=float))
    return {
        "small_max_deg": float(x[0]) % 360.0,
        "small_down_duration_deg": float(x[1]),
        "large_down_duration_deg": float(x[2]),
        "small_max_curvature": float(curv[0]),
        "small_min_curvature": float(curv[1]),
        "large_max_curvature": float(curv[2]),
        "large_min_curvature": float(curv[3]),
        "small_down_kink_u": float(x[7]),
        "small_down_kink_q": float(x[8]),
        "large_up_kink_u": float(x[9]),
        "large_up_kink_q": float(x[10]),
    }


def _make_kin(limits, p):
    return CompactC2FitKinematics(
        limits.small_cylinder,
        limits.large_cylinder,
        **p,
    )


def _metrics(ref_q, ref_dq, ref_ddq, q, dq, ddq):
    return {
        "position_rms": [
            float(np.sqrt(np.mean((q[i]-ref_q[i])**2))) for i in range(2)
        ],
        "combined_position_rms":
            float(np.sqrt(np.mean((q-ref_q)**2))),
        "velocity_rms": [
            float(np.sqrt(np.mean((dq[i]-ref_dq[i])**2))) for i in range(2)
        ],
        "combined_velocity_rms":
            float(np.sqrt(np.mean((dq-ref_dq)**2))),
        "acceleration_rms": [
            float(np.sqrt(np.mean((ddq[i]-ref_ddq[i])**2))) for i in range(2)
        ],
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
    aligned, curvature_seed = _extremum_curvature_estimates(source, sb)

    small_down = (aligned["small_min"]-aligned["small_max"]) % 360.0
    large_down = aligned["large_min"]

    # The kink coordinates are initialized from the excellent geometric fit of
    # the previous 9-parameter family, not from its thermodynamic champion.
    x0 = np.array([
        aligned["small_max"],
        small_down,
        large_down,
        math.log(max(curvature_seed["small_max"], CURVATURE_MIN)),
        math.log(max(curvature_seed["small_min"], CURVATURE_MIN)),
        math.log(max(curvature_seed["large_max"], CURVATURE_MIN)),
        math.log(max(curvature_seed["large_min"], CURVATURE_MIN)),
        0.32900632410553976,
        0.5908614005510023,
        0.7740272614547187,
        0.4245917638230884,
    ], dtype=float)

    lo = np.array([
        0.0,
        MIN_BRANCH_DEG,
        MIN_BRANCH_DEG,
        math.log(CURVATURE_MIN),
        math.log(CURVATURE_MIN),
        math.log(CURVATURE_MIN),
        math.log(CURVATURE_MIN),
        KINK_MIN, KINK_MIN, KINK_MIN, KINK_MIN,
    ])
    hi = np.array([
        360.0,
        MAX_BRANCH_DEG,
        MAX_BRANCH_DEG,
        math.log(CURVATURE_MAX),
        math.log(CURVATURE_MAX),
        math.log(CURVATURE_MAX),
        math.log(CURVATURE_MAX),
        KINK_MAX, KINK_MAX, KINK_MAX, KINK_MAX,
    ])

    tfit = np.linspace(0.0, 1.0, FIT_SAMPLES, endpoint=False)
    ref_q, ref_dq, ref_ddq = _source_aligned(source, source_lmax, tfit)

    def residual(x):
        p = _unpack(x)
        kin = _make_kin(limits, p)
        q, dq, _ = kin.normalized_motor(tfit)
        mono = kin.monotonicity_penalty()
        return np.concatenate((
            (q-ref_q).ravel(),
            VELOCITY_FIT_WEIGHT*(dq-ref_dq).ravel(),
            MONOTONICITY_WEIGHT*mono,
        ))

    fit = least_squares(
        residual,
        x0,
        bounds=(lo, hi),
        x_scale="jac",
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
        max_nfev=5000,
        verbose=1,
    )

    params = _unpack(fit.x)
    kin = _make_kin(limits, params)

    t = np.linspace(0.0, 1.0, 2881)
    ref_q, ref_dq, ref_ddq = _source_aligned(source, source_lmax, t)
    q, dq, ddq = kin.normalized_motor(t)

    metrics = _metrics(ref_q, ref_dq, ref_ddq, q, dq, ddq)
    metrics["small_extrema_count"] = _extrema_count(dq[0])
    metrics["large_extrema_count"] = _extrema_count(dq[1])
    metrics["minimum_increasing_branch_derivative"] = float(
        np.min(kin.monotonicity_penalty() + 0.0)
    )

    result = {
        "model": "11-parameter compact C2 fit with nonzero extremum acceleration",
        "thermodynamic_evaluations": 0,
        "success": bool(fit.success),
        "message": fit.message,
        "nfev": int(fit.nfev),
        "cost": float(fit.cost),
        "velocity_fit_weight": VELOCITY_FIT_WEIGHT,
        "source_extrema_aligned_deg": aligned,
        "source_extremum_curvature_estimates_d2q_dt2": curvature_seed,
        "parameters": params,
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
        ax.plot(360*t, q[j], label="11p C2 fit")
        ax.set_title(f"{labels[j]} — normalized position")
        ax.set_ylabel("q")
        ax.grid(True, alpha=0.25)
        ax.legend()

        ax = axes[1, j]
        ax.plot(360*t, ref_dq[j], "--", label="Source16 aligned")
        ax.plot(360*t, dq[j], label="11p C2 fit")
        ax.axhline(0.0, linewidth=0.8)
        ax.set_title(f"{labels[j]} — normalized velocity")
        ax.set_ylabel("dq/dt")
        ax.grid(True, alpha=0.25)
        ax.legend()

        ax = axes[2, j]
        ax.plot(360*t, ref_ddq[j], "--", label="Source16 aligned")
        ax.plot(360*t, ddq[j], label="11p C2 fit")
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
        "11-parameter compact C2 fit to source16 — no thermodynamics\n"
        f"position RMS={metrics['combined_position_rms']:.5f}   "
        f"velocity RMS={metrics['combined_velocity_rms']:.5f}   "
        f"acceleration RMS={metrics['combined_acceleration_rms']:.3f}"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out/"comparison.png", dpi=180)

    print(json.dumps(result, indent=2))
    print(f"Saved {out/'fit.json'}")
    print(f"Saved {out/'motion.csv'}")
    print(f"Saved {out/'comparison.png'}")


if __name__ == "__main__":
    main()
