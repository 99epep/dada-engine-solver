#!/usr/bin/env python3
"""Compact 7-parameter power-sigmoid motor-motion search at DeltaT = 260 K.

Each piston has exactly one maximum and one minimum. The two monotone branches
between them use

    F_k(u) = u^k / (u^k + (1-u)^k),   0 <= u <= 1, k > 2

so velocity and acceleration vanish at both extrema and the periodic law is C2.

A common cycle rotation is physically redundant here, so the large-piston
maximum is fixed at motor angle 0 deg.

Seven optimized parameters for BOTH pistons:
    small_max_deg
    small_down_duration_deg
    large_down_duration_deg
    small_k_down
    small_k_up
    large_k_down
    large_k_up

Storage is intentionally compact:
  * history.jsonl: only 7 parameters, compact result metrics, warm-start state
  * report.json: compact champion and campaign summary
  * best_diagnostics.json: detailed solver result ONLY for the champion
  * best_motion.csv: phase-aligned source16 and power-sigmoid motion
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from types import SimpleNamespace
import argparse, hashlib, json, math, time

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import qmc

from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
from dada_solver.geometry import CylinderVolumeLimits
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

POWER_FLOOR_W = 25.0
DEFAULT_SOURCE = ROOT / "outputs/motor_spline_from_linear_260k/report.json"
DEFAULT_FOURIER = ROOT / "outputs/motor_fourier_c2_8h_refine/report.json"
DEFAULT_OUTPUT = ROOT / "outputs/motor_power_sigmoid_260k"

MIN_BRANCH_DEG = 35.0
MAX_BRANCH_DEG = 325.0
K_MIN = 2.05
K_MAX = 20.0

ANGLE_RADII_DEG = (18.0, 12.0, 8.0, 5.0, 3.0, 1.6, 0.8, 0.35, 0.15)
LOGK_RADII = (0.55, 0.40, 0.29, 0.20, 0.13, 0.085, 0.050, 0.028, 0.014)
EVALS_PER_RADIUS = 160
WIDE_EVERY = 32
WIDE_ANGLE_DEG = 28.0
WIDE_LOGK = 0.85
DIRECT_WARM_DISTANCE = 0.45


def resolve(path):
    return path if path.is_absolute() else ROOT / path


def _power_sigmoid_scalar(u, k):
    """Allocation-free scalar F_k(u), dF/du for the ODE hot path."""
    if u <= 0.0:
        return 0.0, 0.0
    if u >= 1.0:
        return 1.0, 0.0
    a = u ** k
    b = (1.0 - u) ** k
    den = a + b
    f = a / den
    fp = k * (u ** (k - 1.0)) * ((1.0 - u) ** (k - 1.0)) / (den * den)
    return f, fp


def _power_sigmoid(u, k):
    """Vectorized diagnostic/plotting path."""
    u = np.asarray(u, dtype=float)
    if u.ndim == 0:
        return _power_sigmoid_scalar(float(u), k)
    a = np.power(u, k)
    b = np.power(1.0 - u, k)
    den = a + b
    f = a / den
    fp = np.zeros_like(u)
    interior = (u > 0.0) & (u < 1.0)
    if np.any(interior):
        ui = u[interior]
        deni = np.power(ui, k) + np.power(1.0 - ui, k)
        fp[interior] = (
            k
            * np.power(ui, k - 1.0)
            * np.power(1.0 - ui, k - 1.0)
            / (deni * deni)
        )
    return f, fp


@dataclass(frozen=True, slots=True)
class PowerSigmoidKinematics:
    small_limits: CylinderVolumeLimits
    large_limits: CylinderVolumeLimits
    small_max_deg: float
    small_down_duration_deg: float
    large_down_duration_deg: float
    small_k_down: float
    small_k_up: float
    large_k_down: float
    large_k_up: float

    def __post_init__(self):
        for x in (self.small_down_duration_deg, self.large_down_duration_deg):
            if not MIN_BRANCH_DEG <= x <= MAX_BRANCH_DEG:
                raise ValueError("Branch duration outside admissible range.")
        for k in (
            self.small_k_down, self.small_k_up,
            self.large_k_down, self.large_k_up,
        ):
            if not K_MIN <= k <= K_MAX:
                raise ValueError("Hardness outside admissible range.")

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
        return ()

    @staticmethod
    def _normalized_motor(phi_deg, max_deg, down_deg, k_down, k_up):
        phi = np.asarray(phi_deg, dtype=float) % 360.0
        max_deg = float(max_deg) % 360.0
        min_deg = (max_deg + down_deg) % 360.0

        from_max = np.mod(phi - max_deg, 360.0)
        on_down = from_max <= down_deg

        q = np.empty_like(phi, dtype=float)
        dq_dphi_rad = np.empty_like(phi, dtype=float)

        if np.any(on_down):
            u = from_max[on_down] / down_deg
            f, fp = _power_sigmoid(u, k_down)
            q[on_down] = 1.0 - f
            dq_dphi_rad[on_down] = -fp * (180.0 / math.pi) / down_deg

        if np.any(~on_down):
            up_deg = 360.0 - down_deg
            from_min = np.mod(phi[~on_down] - min_deg, 360.0)
            u = from_min / up_deg
            f, fp = _power_sigmoid(u, k_up)
            q[~on_down] = f
            dq_dphi_rad[~on_down] = fp * (180.0 / math.pi) / up_deg

        if q.ndim == 0:
            return float(q), float(dq_dphi_rad)
        return q, dq_dphi_rad

    @staticmethod
    def _normalized_motor_scalar(phi_deg, max_deg, down_deg, k_down, k_up):
        """Allocation-free scalar motor-law evaluation for the ODE hot path."""
        phi = float(phi_deg) % 360.0
        max_deg = float(max_deg) % 360.0
        from_max = (phi - max_deg) % 360.0

        if from_max <= down_deg:
            u = from_max / down_deg
            f, fp = _power_sigmoid_scalar(u, k_down)
            q = 1.0 - f
            dq_dphi_rad = -fp * (180.0 / math.pi) / down_deg
            return q, dq_dphi_rad

        up_deg = 360.0 - down_deg
        from_min = from_max - down_deg
        u = from_min / up_deg
        f, fp = _power_sigmoid_scalar(u, k_up)
        q = f
        dq_dphi_rad = fp * (180.0 / math.pi) / up_deg
        return q, dq_dphi_rad

    def _pair(self, theta, small):
        if np.ndim(theta) == 0:
            phi = (-float(theta) * 180.0 / math.pi) % 360.0
            if small:
                q, dq = self._normalized_motor_scalar(
                    phi, self.small_max_deg, self.small_down_duration_deg,
                    self.small_k_down, self.small_k_up,
                )
                lim = self.small_limits
            else:
                q, dq = self._normalized_motor_scalar(
                    phi, 0.0, self.large_down_duration_deg,
                    self.large_k_down, self.large_k_up,
                )
                lim = self.large_limits
            return lim.minimum + q * lim.swept, -dq * lim.swept

        a = np.asarray(theta, dtype=float)
        phi = np.mod(-np.degrees(a), 360.0)
        if small:
            q, dq = self._normalized_motor(
                phi, self.small_max_deg, self.small_down_duration_deg,
                self.small_k_down, self.small_k_up,
            )
            lim = self.small_limits
        else:
            q, dq = self._normalized_motor(
                phi, 0.0, self.large_down_duration_deg,
                self.large_k_down, self.large_k_up,
            )
            lim = self.large_limits
        return lim.minimum + q * lim.swept, -dq * lim.swept

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


def params_to_kinematics(limits, p):
    return PowerSigmoidKinematics(
        limits.small_cylinder, limits.large_cylinder,
        float(p["small_max_deg"]) % 360.0,
        float(p["small_down_duration_deg"]),
        float(p["large_down_duration_deg"]),
        float(p["small_k_down"]),
        float(p["small_k_up"]),
        float(p["large_k_down"]),
        float(p["large_k_up"]),
    )


def _source_kinematics(best, limits):
    return _make_kinematics(
        limits,
        np.asarray(best["small_controls"], dtype=float),
        np.asarray(best["large_controls"], dtype=float),
        float(best["small_phase_deg"]),
        float(best["large_phase_deg"]),
    )


def _source_aligned_motion(source, large_max_deg, t):
    tref = np.mod(np.asarray(t, dtype=float) + large_max_deg / 360.0, 1.0)
    return _normalized_motion(source, tref)


def _aligned_shape_distance(source, large_max_deg, candidate, n=4096):
    t = np.linspace(0.0, 1.0, n, endpoint=False)
    ref, dref = _source_aligned_motion(source, large_max_deg, t)
    cur, dcur = _normalized_motion(candidate, t)
    dv = cur - ref
    dd = dcur - dref
    return {
        "small_position_rms": float(np.sqrt(np.mean(dv[0]**2))),
        "large_position_rms": float(np.sqrt(np.mean(dv[1]**2))),
        "combined_position_rms": float(np.sqrt(np.mean(dv**2))),
        "small_velocity_rms": float(np.sqrt(np.mean(dd[0]**2))),
        "large_velocity_rms": float(np.sqrt(np.mean(dd[1]**2))),
        "combined_velocity_rms": float(np.sqrt(np.mean(dd**2))),
        "maximum_position_difference": float(np.max(np.abs(dv))),
    }


def _fit_seed(source, best, limits):
    d = best["shape_diagnostics"]
    lmax = float(d["large_max_deg"])
    raw = {
        "small_max_deg": (float(d["small_max_deg"]) - lmax) % 360.0,
        "small_down_duration_deg": (
            float(d["small_min_deg"]) - float(d["small_max_deg"])
        ) % 360.0,
        "large_down_duration_deg": (
            float(d["large_min_deg"]) - lmax
        ) % 360.0,
    }

    t = np.linspace(0.0, 1.0, 1440, endpoint=False)
    ref, _ = _source_aligned_motion(source, lmax, t)

    h0 = math.log(1.0)
    x0 = np.array([
        raw["small_max_deg"],
        raw["small_down_duration_deg"],
        raw["large_down_duration_deg"],
        h0, h0, h0, h0,
    ])
    hmin = math.log(K_MIN - 2.0)
    hmax = math.log(K_MAX - 2.0)
    lo = np.array([0.0, MIN_BRANCH_DEG, MIN_BRANCH_DEG] + [hmin]*4)
    hi = np.array([360.0, MAX_BRANCH_DEG, MAX_BRANCH_DEG] + [hmax]*4)

    def unpack(x):
        k = 2.0 + np.exp(x[3:])
        return {
            "small_max_deg": float(x[0]) % 360.0,
            "small_down_duration_deg": float(x[1]),
            "large_down_duration_deg": float(x[2]),
            "small_k_down": float(k[0]),
            "small_k_up": float(k[1]),
            "large_k_down": float(k[2]),
            "large_k_up": float(k[3]),
        }

    def residual(x):
        kin = params_to_kinematics(limits, unpack(x))
        cur, _ = _normalized_motion(kin, t)
        return (cur - ref).ravel()

    fit = least_squares(
        residual, x0, bounds=(lo, hi),
        xtol=1e-11, ftol=1e-11, gtol=1e-11,
        max_nfev=2000,
    )
    p = unpack(fit.x)
    kin = params_to_kinematics(limits, p)
    return p, kin, {
        "success": bool(fit.success),
        "cost": float(fit.cost),
        "nfev": int(fit.nfev),
        "shape_distance_vs_source16_aligned":
            _aligned_shape_distance(source, lmax, kin),
    }


def _reflect(x, lo, hi):
    w = hi - lo
    y = (x - lo) % (2.0*w)
    return lo + (2.0*w - y if y > w else y)


def _h(k):
    return math.log(float(k) - 2.0)


def _k(h):
    return 2.0 + math.exp(float(h))


def _propose(center, u, arad, krad):
    p = dict(center)
    u = np.asarray(u, dtype=float)
    p["small_max_deg"] = (
        float(center["small_max_deg"]) + arad*(2*u[0]-1)
    ) % 360.0
    p["small_down_duration_deg"] = _reflect(
        float(center["small_down_duration_deg"]) + arad*(2*u[1]-1),
        MIN_BRANCH_DEG, MAX_BRANCH_DEG,
    )
    p["large_down_duration_deg"] = _reflect(
        float(center["large_down_duration_deg"]) + arad*(2*u[2]-1),
        MIN_BRANCH_DEG, MAX_BRANCH_DEG,
    )

    hmin = math.log(K_MIN - 2.0)
    hmax = math.log(K_MAX - 2.0)
    for j, name in enumerate(
        ("small_k_down", "small_k_up", "large_k_down", "large_k_up")
    ):
        hh = _reflect(
            _h(center[name]) + krad*(2*u[3+j]-1),
            hmin, hmax,
        )
        p[name] = _k(hh)
    return p


def _candidate_id(p):
    return hashlib.sha256(
        json.dumps(p, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _phase_delta(a, b):
    return ((float(a)-float(b)+180.0) % 360.0) - 180.0


def _param_distance(p, r):
    q = r["parameters"]
    terms = [
        _phase_delta(p["small_max_deg"], q["small_max_deg"]) / 30.0,
        (p["small_down_duration_deg"] - q["small_down_duration_deg"]) / 30.0,
        (p["large_down_duration_deg"] - q["large_down_duration_deg"]) / 30.0,
    ]
    for n in ("small_k_down", "small_k_up", "large_k_down", "large_k_up"):
        terms.append((_h(p[n]) - _h(q[n])) / 0.5)
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
        "maximum_absolute_mass_flow_kg_s": r.get("maximum_absolute_mass_flow_kg_s"),
        "maximum_pressure_pa": r.get("maximum_pressure_pa"),
        "maximum_temperature_k": r.get("maximum_temperature_k"),
        "total_mass_kg": r.get("total_mass_kg"),
        "validity_verdict": v.get("verdict"),
        "failed_criteria": v.get("failed_criteria"),
        "cycles_completed": c.get("cycles_completed"),
        "last_normalized_periodic_error": c.get("last_normalized_periodic_error"),
        "normalized_periodic_error_ratio": c.get("normalized_periodic_error_ratio"),
        "message": r.get("message") if r.get("status") != "converged" else None,
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


def _save_motion(path, source, lmax, candidate):
    t = np.linspace(0.0, 1.0, 2881)
    ref, dref = _source_aligned_motion(source, lmax, t)
    cur, dcur = _normalized_motion(candidate, t)
    np.savetxt(
        path,
        np.column_stack((
            t, 360*t,
            ref[0], ref[1], cur[0], cur[1],
            dref[0], dref[1], dcur[0], dcur[1],
        )),
        delimiter=",",
        header=(
            "motor_time_fraction,motor_angle_deg,"
            "source16_aligned_small,source16_aligned_large,"
            "power_sigmoid_small,power_sigmoid_large,"
            "source16_aligned_dsmall_dt,source16_aligned_dlarge_dt,"
            "power_sigmoid_dsmall_dt,power_sigmoid_dlarge_dt"
        ),
        comments="",
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-report", type=Path, default=DEFAULT_SOURCE)
    ap.add_argument("--fourier-report", type=Path, default=DEFAULT_FOURIER)
    ap.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--evaluations", type=int, default=2048)
    ap.add_argument("--budget-seconds", type=float, default=21600)
    ap.add_argument("--candidate-seconds", type=float, default=150)
    ap.add_argument("--evaluations-per-radius", type=int, default=EVALS_PER_RADIUS)
    ap.add_argument("--seed", type=int, default=26092371)
    args = ap.parse_args()

    sp, fp = resolve(args.source_report), resolve(args.fourier_report)
    out = resolve(args.output_directory)
    out.mkdir(parents=True, exist_ok=True)

    sr = json.loads(sp.read_text())
    sb = sr.get("best_feasible")
    if sb is None:
        raise RuntimeError("Source report has no best_feasible.")

    fr = json.loads(fp.read_text())
    fb = _source_best(fr)
    h = _source_harmonics(fb)

    cd, base = _load_basis()
    geometry = base_geometry(base)
    machine = fourier_build_design(
        base, geometry, fb["thermo"],
        np.asarray(fb["small_coefficients"]),
        np.asarray(fb["large_coefficients"]),
        h,
    )
    total_mass = float(sb["result"]["total_mass_kg"])
    machine = _fixed_inventory_design(machine, total_mass)
    limits = machine.configuration.machine_volumes

    source = _source_kinematics(sb, limits)
    source_lmax = float(sb["shape_diagnostics"]["large_max_deg"])
    fit_params, fit_kin, fit_info = _fit_seed(source, sb, limits)

    definition = SimpleNamespace(
        wall_numerical_settings=cd.wall_numerical_settings,
        wall_backend=WallBackendSettings("numba"),
        exact_kinematics_cache=True,
        shared_replay=True,
    )

    identity = {
        "study": "260 K compact power-sigmoid motion v1",
        "function": "F_k(u)=u^k/(u^k+(1-u)^k)",
        "degrees_of_freedom_total": 7,
        "global_phase_gauge": "large maximum fixed at motor angle 0 deg",
        "source_report": str(sp),
        "source_sha256": hashlib.sha256(sp.read_bytes()).hexdigest(),
        "fourier_report": str(fp),
        "fourier_sha256": hashlib.sha256(fp.read_bytes()).hexdigest(),
        "source_candidate_id": sb["candidate_id"],
        "source_large_max_deg_for_alignment": source_lmax,
        "fixed_total_mass_kg": total_mass,
        "fixed_thermo": fb["thermo"],
        "minimum_branch_deg": MIN_BRANCH_DEG,
        "k_bounds": [K_MIN, K_MAX],
        "angle_radii_deg": list(ANGLE_RADII_DEG),
        "logk_radii": list(LOGK_RADII),
        "evaluations_per_radius": args.evaluations_per_radius,
        "wide_probe_every": WIDE_EVERY,
        "wide_angle_deg": WIDE_ANGLE_DEG,
        "wide_logk": WIDE_LOGK,
        "sobol_seed": args.seed,
        "storage_policy": "compact history/report; full solver result only for champion",
    }

    dp, hp, rp = out/"definition.json", out/"history.jsonl", out/"report.json"
    bp = out/"best_diagnostics.json"

    if dp.exists() and json.loads(dp.read_text()) != identity:
        raise RuntimeError("Campaign definition changed; use a new output directory.")
    if not dp.exists():
        dp.write_text(json.dumps(identity, indent=2) + "\n")

    history = _load_history(hp)
    valid = [r for r in history if _eligible(r)]
    best = max(valid, key=_eta) if valid else None
    finished = {
        r["candidate_id"] for r in history
        if r.get("result", {}).get("status") != "interrupted"
    }

    points = qmc.Sobol(d=7, scramble=True, seed=args.seed).random_base2(16)
    start = time.monotonic()
    deadline = start + args.budget_seconds
    completed = 0
    pi = max((r.get("proposal_index", -1) for r in history), default=-1) + 1

    while completed < args.evaluations and time.monotonic() < deadline:
        idx = pi
        pi += 1

        if not history and idx == 0:
            p = dict(fit_params)
            kind = "source16_shape_fit_seed"
            ri = -1
            arad = krad = 0.0
        else:
            searched = sum(
                1 for r in history if r.get("kind") != "source16_shape_fit_seed"
            )
            ri = min(
                searched // args.evaluations_per_radius,
                len(ANGLE_RADII_DEG)-1,
            )
            if searched > 0 and searched % WIDE_EVERY == WIDE_EVERY - 1:
                center = fit_params
                arad, krad = WIDE_ANGLE_DEG, WIDE_LOGK
                kind = "fit_seed_wide_probe"
            else:
                center = best["parameters"] if best is not None else fit_params
                arad, krad = ANGLE_RADII_DEG[ri], LOGK_RADII[ri]
                kind = "power_sigmoid_refinement"
            p = _propose(center, points[idx % len(points)], arad, krad)

        cid = _candidate_id(p)
        if cid in finished:
            continue

        kin = params_to_kinematics(limits, p)
        design = replace(machine, kinematics=kin)
        before = time.monotonic()
        cand_deadline = min(deadline, before + args.candidate_seconds)

        reusable = [
            r for r in history
            if r.get("result", {}).get("status") == "converged"
            and r.get("last_complete_state") is not None
        ]
        warm_distance = None
        if reusable:
            warm = min(reusable, key=lambda r: _param_distance(p, r))
            warm_distance = _param_distance(p, warm)
            warm_source = warm["candidate_id"]
            if warm_distance <= DIRECT_WARM_DISTANCE:
                initial = np.asarray(warm["last_complete_state"])
                warm_mode = "nearest_periodic_state"
            else:
                initial = _safe_uniform_initial(
                    design, wall_source=np.asarray(warm["last_complete_state"])
                )
                warm_mode = "safe_uniform_for_distant_shape"
        else:
            initial = _safe_uniform_initial(
                design, wall_source=np.asarray(sb["last_complete_state"])
            )
            warm_source = sb["candidate_id"]
            warm_mode = "safe_uniform_with_source16_wall_energy"

        def progress(_):
            now = time.monotonic()
            if now >= deadline:
                raise IntegrationInterrupted("Campaign budget exhausted.")
            if now >= cand_deadline:
                raise IntegrationInterrupted("Candidate budget exhausted.")

        box = {}
        def observe(periodic):
            box["statistics"] = periodic.backend_statistics

        safe_retry = False
        try:
            try:
                full, state = _evaluate(
                    f"power_sigmoid_260k_{idx}",
                    design, definition,
                    initial_state=initial,
                    progress_callback=progress,
                    periodic_observer=observe,
                )
            except MicrotubeDomainError:
                safe_retry = True
                safe = _safe_uniform_initial(design, wall_source=initial)
                full, state = _evaluate(
                    f"power_sigmoid_260k_{idx}_safe",
                    design, definition,
                    initial_state=safe,
                    progress_callback=progress,
                    periodic_observer=observe,
                )
                warm_mode = "safe_uniform_retry_after_domain_error"
        except MicrotubeDomainError as exc:
            full, state = {"status":"invalid_exchanger","message":str(exc)}, None
        except IntegrationInterrupted as exc:
            full, state = {"status":"interrupted","message":str(exc)}, None
        except (ValueError, RuntimeError, ArithmeticError) as exc:
            full, state = {"status":"integration_failure","message":str(exc)}, None

        feasible, reasons = (
            feasibility(full, POWER_FLOOR_W)
            if full.get("status") == "converged" else (False, [])
        )
        compact = _compact_result(full)
        shape = _shape_diagnostics(kin)
        distance = _aligned_shape_distance(source, source_lmax, kin, 2048)

        rec = {
            "index": len(history),
            "proposal_index": idx,
            "candidate_id": cid,
            "kind": kind,
            "radius_index": ri,
            "angle_radius_deg": arad,
            "logk_radius": krad,
            "parameters": {k: float(v) for k, v in p.items()},
            "shape_diagnostics": shape,
            "shape_distance_vs_source16_aligned": distance,
            "result": compact,
            "feasible": feasible,
            "physical_constraint_failures": reasons,
            "elapsed_seconds": time.monotonic() - before,
            "warm_start_source": warm_source,
            "warm_start_mode": warm_mode,
            "warm_start_distance": warm_distance,
            "backend": _backend_compact(box.get("statistics")),
            "safe_retry_used": safe_retry,
            "last_complete_state": state.tolist() if state is not None else None,
        }

        with hp.open("a") as f:
            f.write(json.dumps(rec, separators=(",", ":")) + "\n")
            f.flush()

        history.append(rec)
        completed += 1
        if compact.get("status") != "interrupted":
            finished.add(cid)

        moved = False
        if feasible and (best is None or _eta(rec) > _eta(best)):
            best = rec
            moved = True
            _save_motion(out/"best_motion.csv", source, source_lmax, kin)
            bp.write_text(json.dumps({
                "candidate_id": cid,
                "parameters": rec["parameters"],
                "shape_diagnostics": shape,
                "shape_distance_vs_source16_aligned": distance,
                "full_result": full,
            }, indent=2) + "\n")

        report = {
            "best_feasible": best,
            "source16": {
                "candidate_id": sb["candidate_id"],
                "efficiency": sb["result"]["indicated_thermal_efficiency"],
                "power_W": sb["result"]["indicated_power_w"],
                "shape_diagnostics": sb["shape_diagnostics"],
                "phase_alignment_large_max_deg": source_lmax,
            },
            "fitted_seed": {
                "parameters": fit_params,
                "fit": fit_info,
                "shape_diagnostics": _shape_diagnostics(fit_kin),
            },
            "attempted_total": len(history),
            "completed_this_run": completed,
            "requested_duration_seconds": args.budget_seconds,
            "actual_duration_seconds": time.monotonic() - start,
            "definition_file": "definition.json",
            "best_diagnostics_file": "best_diagnostics.json",
        }
        rp.write_text(json.dumps(report, indent=2) + "\n")

        print(json.dumps({
            "index": rec["index"],
            "proposal_index": idx,
            "kind": kind,
            "angle_radius_deg": arad,
            "logk_radius": krad,
            "status": compact.get("status"),
            "feasible": feasible,
            "efficiency": compact.get("indicated_thermal_efficiency"),
            "power_W": compact.get("indicated_power_w"),
            "small_max_deg": p["small_max_deg"],
            "small_min_deg": (
                p["small_max_deg"] + p["small_down_duration_deg"]
            ) % 360.0,
            "large_max_deg": 0.0,
            "large_min_deg": p["large_down_duration_deg"],
            "k": {
                "S_down": p["small_k_down"],
                "S_up": p["small_k_up"],
                "L_down": p["large_k_down"],
                "L_up": p["large_k_up"],
            },
            "rms_vs_source16": distance["combined_position_rms"],
            "center_moved": moved,
            "best_efficiency": _eta(best) if best else None,
            "actual_backend": (rec["backend"] or {}).get("actual_backend"),
            "elapsed_seconds": rec["elapsed_seconds"],
        }, separators=(",", ":")), flush=True)

    print(f"Saved {rp}", flush=True)


if __name__ == "__main__":
    main()
