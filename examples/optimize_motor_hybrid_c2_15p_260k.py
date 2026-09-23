#!/usr/bin/env python3
"""Thermodynamic optimization of the 15-parameter structured C2 motion at ΔT=260 K.

The motion family is exactly the one validated geometrically by
fit_motor_hybrid_c2_15p_to_source16.py:

  * 3 timing/extremum parameters
  * 4 nonzero extremum-curvature magnitudes
  * 2 BP bowing parameters
  * 2 × (kink_u, kink_q, kink_width_rel) on the opposite/HP-side branches

Hardware, exchanger geometry, total gas inventory and thermal conditions are
held fixed. Only kinematics is optimized.

The geometric source16 champion and the 15p fit are used only for initialization
and diagnostics. They do not contribute to the thermodynamic objective.

Storage policy:
  * history.jsonl       compact candidate records + 10-value periodic state
  * report.json         compact campaign/champion summary
  * best_diagnostics.json full solver result only for the current champion
  * best_motion.csv     source16-aligned vs current champion q, dq, ddq
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
import argparse
import hashlib
import json
import math
import time

import numpy as np
from scipy.interpolate import PPoly
from scipy.stats import qmc

from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
from dada_solver.integration import IntegrationInterrupted
from dada_solver.wall_backend import WallBackendSettings

from compare_motor_motion_laws_stage7A5 import ROOT, _evaluate
from optimize_motor_temperature_point import base_geometry, feasibility
from refine_motor_four_stage_hx9d_variable_gas import _load_basis
from refine_motor_fourier_c2_260k import build_design as fourier_build_design
from optimize_motor_free_spline_260k_v3 import (
    _fixed_inventory_design,
    _load_history,
    _make_kinematics,
    _normalized_motion,
    _safe_uniform_initial,
    _shape_diagnostics,
    _source_best,
    _source_harmonics,
)
from fit_motor_hybrid_c2_15p_to_source16 import (
    StructuredMotion15,
    _source_aligned,
    MIN_BRANCH_DEG,
    MAX_BRANCH_DEG,
    CURVATURE_MIN,
    CURVATURE_MAX,
    BP_MID_MIN,
    BP_MID_MAX,
    KINK_MIN,
    KINK_MAX,
    KINK_WIDTH_REL_MIN,
    KINK_WIDTH_REL_MAX,
)

POWER_FLOOR_W = 25.0

DEFAULT_SOURCE = ROOT / "outputs/motor_spline_from_linear_260k/report.json"
DEFAULT_FOURIER = ROOT / "outputs/motor_fourier_c2_8h_refine/report.json"
DEFAULT_FIT = ROOT / "outputs/motor_hybrid_c2_15p_fit_260k/fit.json"
DEFAULT_OUTPUT = ROOT / "outputs/motor_hybrid_c2_15p_260k"

ANGLE_RADII_DEG = (12.0, 8.0, 5.0, 3.0, 1.7, 0.9, 0.45, 0.20, 0.08)
LOG_CURV_RADII = (0.65, 0.48, 0.35, 0.24, 0.16, 0.10, 0.060, 0.035, 0.018)
BP_MID_RADII = (0.090, 0.070, 0.050, 0.035, 0.022, 0.014, 0.009, 0.005, 0.0025)
KINK_UQ_RADII = (0.100, 0.075, 0.055, 0.038, 0.025, 0.016, 0.010, 0.006, 0.003)
LOG_WIDTH_RADII = (0.80, 0.60, 0.45, 0.32, 0.22, 0.14, 0.09, 0.05, 0.025)

EVALS_PER_RADIUS = 112
WIDE_EVERY = 64

WIDE_ANGLE_DEG = 20.0
WIDE_LOG_CURV = 0.95
WIDE_BP_MID = 0.14
WIDE_KINK_UQ = 0.16
WIDE_LOG_WIDTH = 1.10

DIRECT_WARM_DISTANCE = 0.65
STRUCTURAL_SCREEN_SAMPLES = 384


def resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


class _FastPPoly:
    """Allocation-free scalar PPoly evaluator for the ODE hot path."""

    __slots__ = ("x", "c", "dc")

    def __init__(self, bpoly):
        pp = PPoly.from_bernstein_basis(bpoly)
        dp = pp.derivative(1)
        self.x = np.asarray(pp.x, dtype=float)
        self.c = np.asarray(pp.c, dtype=float)
        self.dc = np.asarray(dp.c, dtype=float)

    @staticmethod
    def _horner(coeffs, z):
        y = 0.0
        for a in coeffs:
            y = y*z + float(a)
        return y

    def pair(self, u):
        u = float(u)
        if u <= self.x[0]:
            i = 0
            z = 0.0
        elif u >= self.x[-1]:
            i = len(self.x)-2
            z = self.x[-1]-self.x[-2]
        else:
            i = int(np.searchsorted(self.x, u, side="right")-1)
            i = max(0, min(i, len(self.x)-2))
            z = u-self.x[i]

        return (
            self._horner(self.c[:, i], z),
            self._horner(self.dc[:, i], z),
        )


class StructuredKinematics15:
    """Solver-compatible adapter around StructuredMotion15.

    Vector paths use the validated BPoly representation.
    Scalar RHS paths use preconverted power-basis polynomials and manual Horner
    evaluation to avoid repeated scipy interpolation overhead.
    """

    __slots__ = (
        "small_limits", "large_limits", "params", "motion",
        "_small_max_t", "_sd", "_ld",
        "_small_down", "_small_up", "_large_down", "_large_up",
    )

    def __init__(self, limits, params):
        self.small_limits = limits.small_cylinder
        self.large_limits = limits.large_cylinder
        self.params = dict(params)
        self.motion = StructuredMotion15(self.params)

        self._small_max_t = (self.params["small_max_deg"] % 360.0)/360.0
        self._sd = self.params["small_down_duration_deg"]/360.0
        self._ld = self.params["large_down_duration_deg"]/360.0

        self._small_down = _FastPPoly(self.motion.small_down)
        self._small_up = _FastPPoly(self.motion.small_up)
        self._large_down = _FastPPoly(self.motion.large_down)
        self._large_up = _FastPPoly(self.motion.large_up)

    @property
    def small_volume_limits(self):
        return self.small_limits

    @property
    def large_volume_limits(self):
        return self.large_limits

    @property
    def small_physical_stroke(self):
        return None

    @property
    def large_physical_stroke(self):
        return None

    @property
    def diagnostics(self):
        return ()

    def require_feasible(self):
        return None

    def breakpoint_angles(self):
        # All joins are C2. No discontinuous solver breakpoint is needed.
        return ()

    def _small_scalar_q_dqdt(self, t):
        from_max = (float(t)-self._small_max_t) % 1.0

        if from_max <= self._sd:
            u = from_max/self._sd
            y, dy = self._small_down.pair(u)
            return 1.0-y, -dy/self._sd

        su = 1.0-self._sd
        u = (from_max-self._sd)/su
        y, dy = self._small_up.pair(u)
        return y, dy/su

    def _large_scalar_q_dqdt(self, t):
        t = float(t) % 1.0

        if t <= self._ld:
            u = t/self._ld
            y, dy = self._large_down.pair(u)
            return 1.0-y, -dy/self._ld

        lu = 1.0-self._ld
        u = (t-self._ld)/lu
        y, dy = self._large_up.pair(u)
        return y, dy/lu

    def _pair(self, theta, small):
        if np.ndim(theta) == 0:
            t = (-float(theta)/(2.0*math.pi)) % 1.0
            if small:
                q, dqdt = self._small_scalar_q_dqdt(t)
                lim = self.small_limits
            else:
                q, dqdt = self._large_scalar_q_dqdt(t)
                lim = self.large_limits

            v = lim.minimum + q*lim.swept
            dvdtheta = -(dqdt/(2.0*math.pi))*lim.swept
            return v, dvdtheta

        a = np.asarray(theta, dtype=float)
        t = np.mod(-a/(2.0*math.pi), 1.0)
        q, dqdt, _ = self.motion.normalized(t)
        row = 0 if small else 1
        lim = self.small_limits if small else self.large_limits
        return (
            lim.minimum + q[row]*lim.swept,
            -(dqdt[row]/(2.0*math.pi))*lim.swept,
        )

    def small_cylinder_volume(self, theta):
        return self._pair(theta, True)[0]

    def large_cylinder_volume(self, theta):
        return self._pair(theta, False)[0]

    def small_cylinder_volume_derivative(self, theta):
        return self._pair(theta, True)[1]

    def large_cylinder_volume_derivative(self, theta):
        return self._pair(theta, False)[1]

    def cylinder_volumes_and_derivatives(self, theta):
        sv, ds = self._pair(theta, True)
        lv, dl = self._pair(theta, False)
        return sv, lv, ds, dl


def _source_kinematics(best, limits):
    return _make_kinematics(
        limits,
        np.asarray(best["small_controls"], dtype=float),
        np.asarray(best["large_controls"], dtype=float),
        float(best["small_phase_deg"]),
        float(best["large_phase_deg"]),
    )


def _source_aligned_motion(source, large_max_deg, t):
    tref = np.mod(np.asarray(t, dtype=float) + large_max_deg/360.0, 1.0)
    return _normalized_motion(source, tref)


def _aligned_shape_distance(source, large_max_deg, candidate, n=4096):
    t = np.linspace(0.0, 1.0, n, endpoint=False)
    ref, dref = _source_aligned_motion(source, large_max_deg, t)
    cur, dcur = _normalized_motion(candidate, t)

    dv = cur-ref
    dd = dcur-dref

    return {
        "small_position_rms": float(np.sqrt(np.mean(dv[0]**2))),
        "large_position_rms": float(np.sqrt(np.mean(dv[1]**2))),
        "combined_position_rms": float(np.sqrt(np.mean(dv**2))),
        "small_velocity_rms": float(np.sqrt(np.mean(dd[0]**2))),
        "large_velocity_rms": float(np.sqrt(np.mean(dd[1]**2))),
        "combined_velocity_rms": float(np.sqrt(np.mean(dd**2))),
        "maximum_position_difference": float(np.max(np.abs(dv))),
    }


def _structurally_admissible(kin):
    u = np.linspace(0.0, 1.0, STRUCTURAL_SCREEN_SAMPLES)
    branches = (
        kin.motion.small_down,
        kin.motion.small_up,
        kin.motion.large_down,
        kin.motion.large_up,
    )
    for branch in branches:
        d = np.asarray(branch.derivative(1)(u), dtype=float)
        if not np.all(np.isfinite(d)):
            return False
        if float(np.min(d)) < -1e-9:
            return False
    return True


def _reflect(x, lo, hi):
    w = hi-lo
    y = (float(x)-lo) % (2.0*w)
    return lo + (2.0*w-y if y > w else y)


def _log_reflect(value, delta, lo, hi):
    return math.exp(
        _reflect(
            math.log(float(value)) + float(delta),
            math.log(lo),
            math.log(hi),
        )
    )


def _propose(center, u, angle_r, log_curv_r, bp_r, kink_r, log_width_r):
    u = np.asarray(u, dtype=float)
    p = dict(center)

    p["small_max_deg"] = (
        float(center["small_max_deg"]) + angle_r*(2*u[0]-1)
    ) % 360.0

    p["small_down_duration_deg"] = _reflect(
        float(center["small_down_duration_deg"]) + angle_r*(2*u[1]-1),
        MIN_BRANCH_DEG,
        MAX_BRANCH_DEG,
    )
    p["large_down_duration_deg"] = _reflect(
        float(center["large_down_duration_deg"]) + angle_r*(2*u[2]-1),
        MIN_BRANCH_DEG,
        MAX_BRANCH_DEG,
    )

    curvature_names = (
        "small_max_curvature",
        "small_min_curvature",
        "large_max_curvature",
        "large_min_curvature",
    )
    for j, name in enumerate(curvature_names):
        p[name] = _log_reflect(
            center[name],
            log_curv_r*(2*u[3+j]-1),
            CURVATURE_MIN,
            CURVATURE_MAX,
        )

    p["small_up_bp_mid_q"] = _reflect(
        float(center["small_up_bp_mid_q"]) + bp_r*(2*u[7]-1),
        BP_MID_MIN,
        BP_MID_MAX,
    )
    p["large_down_bp_mid_q"] = _reflect(
        float(center["large_down_bp_mid_q"]) + bp_r*(2*u[8]-1),
        BP_MID_MIN,
        BP_MID_MAX,
    )

    # Small non-BP/HP-side kink.
    p["small_down_kink_u"] = _reflect(
        float(center["small_down_kink_u"]) + kink_r*(2*u[9]-1),
        KINK_MIN,
        KINK_MAX,
    )
    p["small_down_kink_q"] = _reflect(
        float(center["small_down_kink_q"]) + kink_r*(2*u[10]-1),
        KINK_MIN,
        KINK_MAX,
    )
    p["small_down_kink_width_rel"] = _log_reflect(
        center["small_down_kink_width_rel"],
        log_width_r*(2*u[11]-1),
        KINK_WIDTH_REL_MIN,
        KINK_WIDTH_REL_MAX,
    )

    # Large non-BP/HP-side kink.
    p["large_up_kink_u"] = _reflect(
        float(center["large_up_kink_u"]) + kink_r*(2*u[12]-1),
        KINK_MIN,
        KINK_MAX,
    )
    p["large_up_kink_q"] = _reflect(
        float(center["large_up_kink_q"]) + kink_r*(2*u[13]-1),
        KINK_MIN,
        KINK_MAX,
    )
    p["large_up_kink_width_rel"] = _log_reflect(
        center["large_up_kink_width_rel"],
        log_width_r*(2*u[14]-1),
        KINK_WIDTH_REL_MIN,
        KINK_WIDTH_REL_MAX,
    )

    return p


def _candidate_id(p):
    return hashlib.sha256(
        json.dumps(p, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _phase_delta(a, b):
    return ((float(a)-float(b)+180.0) % 360.0)-180.0


def _param_distance(p, r):
    q = r["parameters"]

    terms = [
        _phase_delta(p["small_max_deg"], q["small_max_deg"])/25.0,
        (p["small_down_duration_deg"]-q["small_down_duration_deg"])/25.0,
        (p["large_down_duration_deg"]-q["large_down_duration_deg"])/25.0,
    ]

    for name in (
        "small_max_curvature",
        "small_min_curvature",
        "large_max_curvature",
        "large_min_curvature",
    ):
        terms.append(
            (math.log(p[name])-math.log(q[name]))/0.55
        )

    terms.extend((
        (p["small_up_bp_mid_q"]-q["small_up_bp_mid_q"])/0.10,
        (p["large_down_bp_mid_q"]-q["large_down_bp_mid_q"])/0.10,
    ))

    for name in (
        "small_down_kink_u",
        "small_down_kink_q",
        "large_up_kink_u",
        "large_up_kink_q",
    ):
        terms.append((p[name]-q[name])/0.14)

    for name in (
        "small_down_kink_width_rel",
        "large_up_kink_width_rel",
    ):
        terms.append(
            (math.log(p[name])-math.log(q[name]))/0.75
        )

    return float(np.linalg.norm(terms))


def _compact_result(r):
    c = r.get("convergence") or {}
    v = r.get("validity") or {}

    return {
        "status": r.get("status"),
        "indicated_thermal_efficiency": r.get("indicated_thermal_efficiency"),
        "indicated_power_w": r.get("indicated_power_w"),
        "heat_input_w": r.get("heat_input_w"),
        "heat_out_w": r.get("heat_out_w"),
        "maximum_tube_reynolds": r.get("maximum_tube_reynolds"),
        "maximum_tube_mach_number": r.get("maximum_tube_mach_number"),
        "maximum_absolute_mass_flow_kg_s":
            r.get("maximum_absolute_mass_flow_kg_s"),
        "maximum_pressure_pa": r.get("maximum_pressure_pa"),
        "maximum_temperature_k": r.get("maximum_temperature_k"),
        "total_mass_kg": r.get("total_mass_kg"),
        "validity_verdict": v.get("verdict"),
        "failed_criteria": v.get("failed_criteria"),
        "cycles_completed": c.get("cycles_completed"),
        "last_normalized_periodic_error":
            c.get("last_normalized_periodic_error"),
        "normalized_periodic_error_ratio":
            c.get("normalized_periodic_error_ratio"),
        "message":
            r.get("message") if r.get("status") != "converged" else None,
    }


def _eligible(r):
    return bool(
        r.get("feasible")
        and r.get("result", {}).get("status") == "converged"
        and r["result"].get("indicated_thermal_efficiency") is not None
    )


def _eta(r):
    return float(r["result"]["indicated_thermal_efficiency"])


def _backend_compact(s):
    if not s:
        return None

    cache = s.get("exact_kinematics_cache") or {}
    return {
        "actual_backend": s.get("actual_backend"),
        "fallback_calls": s.get("fallback_calls"),
        "exact_cache_hit_fraction": cache.get("hit_fraction"),
    }


def _save_motion(path, source, source_lmax, kin):
    t = np.linspace(0.0, 1.0, 2881)

    ref_q, ref_dq, ref_ddq = _source_aligned(
        source, source_lmax, t
    )
    q, dq, ddq = kin.motion.normalized(t)

    np.savetxt(
        path,
        np.column_stack((
            t,
            360.0*t,
            ref_q[0],
            ref_q[1],
            q[0],
            q[1],
            ref_dq[0],
            ref_dq[1],
            dq[0],
            dq[1],
            ref_ddq[0],
            ref_ddq[1],
            ddq[0],
            ddq[1],
        )),
        delimiter=",",
        header=(
            "motor_time_fraction,motor_angle_deg,"
            "source16_aligned_small,source16_aligned_large,"
            "structured15_small,structured15_large,"
            "source16_aligned_dsmall_dt,source16_aligned_dlarge_dt,"
            "structured15_dsmall_dt,structured15_dlarge_dt,"
            "source16_aligned_ddsmall_dt2,source16_aligned_ddlarge_dt2,"
            "structured15_ddsmall_dt2,structured15_ddlarge_dt2"
        ),
        comments="",
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-report", type=Path, default=DEFAULT_SOURCE)
    ap.add_argument("--fourier-report", type=Path, default=DEFAULT_FOURIER)
    ap.add_argument("--fit-report", type=Path, default=DEFAULT_FIT)
    ap.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT)

    ap.add_argument("--evaluations", type=int, default=1024)
    ap.add_argument("--evaluations-per-radius", type=int, default=EVALS_PER_RADIUS)
    ap.add_argument("--budget-seconds", type=float, default=21600)
    ap.add_argument("--candidate-seconds", type=float, default=150)
    ap.add_argument("--seed", type=int, default=260923151)

    args = ap.parse_args()

    sp = resolve(args.source_report)
    fp = resolve(args.fourier_report)
    fitp = resolve(args.fit_report)
    out = resolve(args.output_directory)
    out.mkdir(parents=True, exist_ok=True)

    sr = json.loads(sp.read_text())
    sb = sr.get("best_feasible")
    if sb is None:
        raise RuntimeError("Source report has no best_feasible.")

    fit_report = json.loads(fitp.read_text())
    if not fit_report.get("success"):
        raise RuntimeError("15p geometric fit did not succeed.")
    fit_params = {
        k: float(v)
        for k, v in fit_report["parameters"].items()
    }

    fr = json.loads(fp.read_text())
    fb = _source_best(fr)
    harmonics = _source_harmonics(fb)

    cd, base = _load_basis()
    geometry = base_geometry(base)

    machine = fourier_build_design(
        base,
        geometry,
        fb["thermo"],
        np.asarray(fb["small_coefficients"]),
        np.asarray(fb["large_coefficients"]),
        harmonics,
    )

    total_mass = float(sb["result"]["total_mass_kg"])
    machine = _fixed_inventory_design(machine, total_mass)
    limits = machine.configuration.machine_volumes

    source = _source_kinematics(sb, limits)
    source_lmax = float(sb["shape_diagnostics"]["large_max_deg"])

    fit_kin = StructuredKinematics15(limits, fit_params)
    if not _structurally_admissible(fit_kin):
        raise RuntimeError("Geometric 15p fit is unexpectedly non-monotone.")

    fit_distance = _aligned_shape_distance(
        source, source_lmax, fit_kin
    )

    print(json.dumps({
        "structured15_thermo_seed": {
            "parameters": fit_params,
            "fit_metrics": fit_report.get("metrics"),
            "recomputed_shape_distance": fit_distance,
        }
    }, separators=(",", ":")), flush=True)

    definition = SimpleNamespace(
        wall_numerical_settings=cd.wall_numerical_settings,
        wall_backend=WallBackendSettings("numba"),
        exact_kinematics_cache=True,
        shared_replay=True,
    )

    identity = {
        "study": "260 K structured C2 15-parameter thermodynamic optimization v1",
        "degrees_of_freedom_total": 15,
        "global_phase_gauge": "large maximum fixed at motor angle 0 deg",
        "source_report": str(sp),
        "source_sha256": hashlib.sha256(sp.read_bytes()).hexdigest(),
        "fourier_report": str(fp),
        "fourier_sha256": hashlib.sha256(fp.read_bytes()).hexdigest(),
        "fit_report": str(fitp),
        "fit_sha256": hashlib.sha256(fitp.read_bytes()).hexdigest(),
        "source_candidate_id": sb["candidate_id"],
        "source_large_max_deg_for_alignment": source_lmax,
        "fixed_total_mass_kg": total_mass,
        "fixed_thermo": fb["thermo"],
        "parameter_bounds": {
            "branch_deg": [MIN_BRANCH_DEG, MAX_BRANCH_DEG],
            "curvature": [CURVATURE_MIN, CURVATURE_MAX],
            "bp_mid_q": [BP_MID_MIN, BP_MID_MAX],
            "kink_u_q": [KINK_MIN, KINK_MAX],
            "kink_width_rel": [KINK_WIDTH_REL_MIN, KINK_WIDTH_REL_MAX],
        },
        "search_radii": {
            "angle_deg": list(ANGLE_RADII_DEG),
            "log_curvature": list(LOG_CURV_RADII),
            "bp_mid_q": list(BP_MID_RADII),
            "kink_u_q": list(KINK_UQ_RADII),
            "log_kink_width": list(LOG_WIDTH_RADII),
        },
        "evaluations_per_radius": args.evaluations_per_radius,
        "wide_probe_every": WIDE_EVERY,
        "wide_probe": {
            "angle_deg": WIDE_ANGLE_DEG,
            "log_curvature": WIDE_LOG_CURV,
            "bp_mid_q": WIDE_BP_MID,
            "kink_u_q": WIDE_KINK_UQ,
            "log_kink_width": WIDE_LOG_WIDTH,
        },
        "structural_screen_samples": STRUCTURAL_SCREEN_SAMPLES,
        "sobol_seed": args.seed,
        "storage_policy":
            "compact history/report; full solver result only for champion",
    }

    dp = out/"definition.json"
    hp = out/"history.jsonl"
    rp = out/"report.json"
    bp = out/"best_diagnostics.json"

    if dp.exists() and json.loads(dp.read_text()) != identity:
        raise RuntimeError(
            "Campaign definition changed; use a new output directory."
        )
    if not dp.exists():
        dp.write_text(json.dumps(identity, indent=2)+"\n")

    history = _load_history(hp)
    valid = [r for r in history if _eligible(r)]
    best = max(valid, key=_eta) if valid else None

    finished = {
        r["candidate_id"]
        for r in history
        if r.get("result", {}).get("status") != "interrupted"
    }

    points = qmc.Sobol(
        d=15, scramble=True, seed=args.seed
    ).random_base2(16)

    start = time.monotonic()
    deadline = start+args.budget_seconds
    completed = 0

    proposal_index = max(
        (r.get("proposal_index", -1) for r in history),
        default=-1,
    )+1

    skipped_structural = 0

    while completed < args.evaluations and time.monotonic() < deadline:
        idx = proposal_index
        proposal_index += 1

        if not history and idx == 0:
            params = dict(fit_params)
            kind = "source16_15p_geometric_fit_seed"
            radius_index = -1
            angle_r = log_curv_r = bp_r = kink_r = log_width_r = 0.0
        else:
            searched = sum(
                1
                for r in history
                if r.get("kind") != "source16_15p_geometric_fit_seed"
            )

            radius_index = min(
                searched//args.evaluations_per_radius,
                len(ANGLE_RADII_DEG)-1,
            )

            if searched > 0 and searched % WIDE_EVERY == WIDE_EVERY-1:
                center = fit_params
                angle_r = WIDE_ANGLE_DEG
                log_curv_r = WIDE_LOG_CURV
                bp_r = WIDE_BP_MID
                kink_r = WIDE_KINK_UQ
                log_width_r = WIDE_LOG_WIDTH
                kind = "fit_seed_wide_probe"
            else:
                center = (
                    best["parameters"]
                    if best is not None
                    else fit_params
                )
                angle_r = ANGLE_RADII_DEG[radius_index]
                log_curv_r = LOG_CURV_RADII[radius_index]
                bp_r = BP_MID_RADII[radius_index]
                kink_r = KINK_UQ_RADII[radius_index]
                log_width_r = LOG_WIDTH_RADII[radius_index]
                kind = "structured15_refinement"

            params = _propose(
                center,
                points[idx % len(points)],
                angle_r,
                log_curv_r,
                bp_r,
                kink_r,
                log_width_r,
            )

        cid = _candidate_id(params)
        if cid in finished:
            continue

        try:
            kin = StructuredKinematics15(limits, params)
        except (ValueError, FloatingPointError):
            skipped_structural += 1
            continue

        if not _structurally_admissible(kin):
            skipped_structural += 1
            continue

        design = replace(machine, kinematics=kin)

        before = time.monotonic()
        candidate_deadline = min(
            deadline, before+args.candidate_seconds
        )

        reusable = [
            r for r in history
            if r.get("result", {}).get("status") == "converged"
            and r.get("last_complete_state") is not None
        ]

        warm_distance = None
        if reusable:
            warm = min(
                reusable,
                key=lambda r: _param_distance(params, r),
            )
            warm_distance = _param_distance(params, warm)
            warm_source = warm["candidate_id"]

            if warm_distance <= DIRECT_WARM_DISTANCE:
                initial = np.asarray(
                    warm["last_complete_state"], dtype=float
                )
                warm_mode = "nearest_periodic_state"
            else:
                initial = _safe_uniform_initial(
                    design,
                    wall_source=np.asarray(
                        warm["last_complete_state"], dtype=float
                    ),
                )
                warm_mode = "safe_uniform_for_distant_shape"
        else:
            initial = _safe_uniform_initial(
                design,
                wall_source=np.asarray(
                    sb["last_complete_state"], dtype=float
                ),
            )
            warm_source = sb["candidate_id"]
            warm_mode = "safe_uniform_with_source16_wall_energy"

        def progress(_):
            now = time.monotonic()
            if now >= deadline:
                raise IntegrationInterrupted("Campaign budget exhausted.")
            if now >= candidate_deadline:
                raise IntegrationInterrupted("Candidate budget exhausted.")

        box = {}

        def observe(periodic):
            box["statistics"] = periodic.backend_statistics

        safe_retry = False

        try:
            try:
                full, state = _evaluate(
                    f"structured15_260k_{idx}",
                    design,
                    definition,
                    initial_state=initial,
                    progress_callback=progress,
                    periodic_observer=observe,
                )
            except MicrotubeDomainError:
                safe_retry = True
                safe = _safe_uniform_initial(
                    design, wall_source=initial
                )
                full, state = _evaluate(
                    f"structured15_260k_{idx}_safe",
                    design,
                    definition,
                    initial_state=safe,
                    progress_callback=progress,
                    periodic_observer=observe,
                )
                warm_mode = "safe_uniform_retry_after_domain_error"

        except MicrotubeDomainError as exc:
            full = {
                "status": "invalid_exchanger",
                "message": str(exc),
            }
            state = None

        except IntegrationInterrupted as exc:
            full = {
                "status": "interrupted",
                "message": str(exc),
            }
            state = None

        except (ValueError, RuntimeError, ArithmeticError) as exc:
            full = {
                "status": "integration_failure",
                "message": str(exc),
            }
            state = None

        feasible, reasons = (
            feasibility(full, POWER_FLOOR_W)
            if full.get("status") == "converged"
            else (False, [])
        )

        compact = _compact_result(full)
        shape = _shape_diagnostics(kin)
        distance = _aligned_shape_distance(
            source, source_lmax, kin, 2048
        )

        rec = {
            "index": len(history),
            "proposal_index": idx,
            "candidate_id": cid,
            "kind": kind,
            "radius_index": radius_index,
            "search_radii": {
                "angle_deg": angle_r,
                "log_curvature": log_curv_r,
                "bp_mid_q": bp_r,
                "kink_u_q": kink_r,
                "log_kink_width": log_width_r,
            },
            "parameters": {
                k: float(v) for k, v in params.items()
            },
            "shape_diagnostics": shape,
            "shape_distance_vs_source16_aligned": distance,
            "result": compact,
            "feasible": feasible,
            "physical_constraint_failures": reasons,
            "elapsed_seconds": time.monotonic()-before,
            "warm_start_source": warm_source,
            "warm_start_mode": warm_mode,
            "warm_start_distance": warm_distance,
            "backend": _backend_compact(box.get("statistics")),
            "safe_retry_used": safe_retry,
            "last_complete_state":
                state.tolist() if state is not None else None,
        }

        with hp.open("a") as f:
            f.write(json.dumps(rec, separators=(",", ":"))+"\n")
            f.flush()

        history.append(rec)
        completed += 1

        if compact.get("status") != "interrupted":
            finished.add(cid)

        moved = False
        if feasible and (best is None or _eta(rec) > _eta(best)):
            best = rec
            moved = True

            _save_motion(
                out/"best_motion.csv",
                source,
                source_lmax,
                kin,
            )

            bp.write_text(
                json.dumps({
                    "candidate_id": cid,
                    "parameters": rec["parameters"],
                    "shape_diagnostics": shape,
                    "shape_distance_vs_source16_aligned": distance,
                    "full_result": full,
                }, indent=2)+"\n"
            )

        report = {
            "best_feasible": best,
            "source16": {
                "candidate_id": sb["candidate_id"],
                "efficiency":
                    sb["result"]["indicated_thermal_efficiency"],
                "power_W":
                    sb["result"]["indicated_power_w"],
                "shape_diagnostics":
                    sb["shape_diagnostics"],
                "phase_alignment_large_max_deg":
                    source_lmax,
            },
            "geometric_fit_seed": {
                "parameters": fit_params,
                "fit_metrics": fit_report.get("metrics"),
                "shape_distance_recomputed":
                    fit_distance,
            },
            "attempted_total": len(history),
            "completed_this_run": completed,
            "skipped_structural_total_this_run":
                skipped_structural,
            "requested_duration_seconds":
                args.budget_seconds,
            "actual_duration_seconds":
                time.monotonic()-start,
            "definition_file": "definition.json",
            "best_diagnostics_file":
                "best_diagnostics.json",
        }

        rp.write_text(json.dumps(report, indent=2)+"\n")

        print(json.dumps({
            "index": rec["index"],
            "proposal_index": idx,
            "kind": kind,
            "radius_index": radius_index,
            "status": compact.get("status"),
            "feasible": feasible,
            "efficiency":
                compact.get("indicated_thermal_efficiency"),
            "power_W":
                compact.get("indicated_power_w"),
            "small_max_deg":
                params["small_max_deg"],
            "small_min_deg": (
                params["small_max_deg"]
                + params["small_down_duration_deg"]
            ) % 360.0,
            "large_max_deg": 0.0,
            "large_min_deg":
                params["large_down_duration_deg"],
            "bp_mid_q": {
                "small_up":
                    params["small_up_bp_mid_q"],
                "large_down":
                    params["large_down_bp_mid_q"],
            },
            "small_down_kink": [
                params["small_down_kink_u"],
                params["small_down_kink_q"],
                params["small_down_kink_width_rel"],
            ],
            "large_up_kink": [
                params["large_up_kink_u"],
                params["large_up_kink_q"],
                params["large_up_kink_width_rel"],
            ],
            "extremum_curvature": {
                "Smax": params["small_max_curvature"],
                "Smin": params["small_min_curvature"],
                "Lmax": params["large_max_curvature"],
                "Lmin": params["large_min_curvature"],
            },
            "rms_vs_source16":
                distance["combined_position_rms"],
            "center_moved": moved,
            "best_efficiency":
                _eta(best) if best else None,
            "actual_backend":
                (rec["backend"] or {}).get("actual_backend"),
            "elapsed_seconds":
                rec["elapsed_seconds"],
            "skipped_structural":
                skipped_structural,
        }, separators=(",", ":")), flush=True)

    print(f"Saved {rp}", flush=True)


if __name__ == "__main__":
    main()
