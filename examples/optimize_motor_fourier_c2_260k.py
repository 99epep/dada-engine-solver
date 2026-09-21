#!/usr/bin/env python3
"""Large smooth-motion search around the 260 K UU four-stage champion.

The four-stage piecewise-linear motion is projected onto a periodic Fourier
basis.  The resulting motion is C-infinity (therefore C2) and uses the full
stroke of each cylinder after exact numerical normalization.

Free variables:
  - 16 kinematic shape coefficients: 4 Fourier harmonics x sin/cos x 2 pistons
  - 5 tightly bounded thermodynamic/hardware variables:
      swept_ratio, n_i, length_i_m, n_o, length_o_m

The thermodynamic variables are hard-limited to +/-3 % around the retained
260 K UU champion.  Kinematics has much more freedom.

The search is resumable and append-only.

Typical run:
    PYTHONPATH=src python3 examples/optimize_motor_fourier_c2_260k.py \
        --budget-seconds 21600 --evaluations 512 --candidate-seconds 180

Outputs:
    outputs/motor_fourier_c2_260k/history.jsonl
    outputs/motor_fourier_c2_260k/report.json
    outputs/motor_fourier_c2_260k/best_motion.csv
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
from types import SimpleNamespace
import argparse
import hashlib
import json
import math
import time

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import qmc

from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.integration import IntegrationInterrupted
from dada_solver.wall_backend import WallBackendSettings

from compare_motor_motion_laws_stage7A5 import ROOT, _evaluate
from optimize_motor_temperature_point import (
    COLD_K,
    CHARGE_PA,
    base_geometry,
    build_design as _linear_build_design,
)
from optimize_motor_temperature_point import (
    feasibility,
    load_source,
    uniform_metadata,
    warm_start,
)
from refine_motor_four_stage_hx9d_variable_gas import _load_basis


SOURCE_REPORT = ROOT / "outputs" / "motor_valve_optimization_260K" / "UU" / "report.json"
OUTPUT_DIRECTORY = ROOT / "outputs" / "motor_fourier_c2_260k"

DELTA_T_K = 260.0
HOT_K = COLD_K + DELTA_T_K
POWER_FLOOR_W = 25.0
THERMO_RELATIVE_LIMIT = 0.03
HARMONICS = 4
SHAPE_DIM = 2 * HARMONICS * 2
THERMO_NAMES = ("swept_ratio", "n_i", "length_i_m", "n_o", "length_o_m")

# Broad-to-local coefficient radii.  Coefficient vectors are renormalized after
# every proposal, so these are angular-ish perturbation scales in shape space.
DEFAULT_SHAPE_RADII = (0.45, 0.28, 0.16, 0.09, 0.05, 0.025)
DEFAULT_EVALUATIONS_PER_RADIUS = 64
SOBOL_SEED = 26092101


@dataclass(frozen=True)
class FourierVolumeKinematics:
    """Periodic full-stroke Fourier kinematics, smooth to all orders."""

    small_volume_limits: CylinderVolumeLimits
    large_volume_limits: CylinderVolumeLimits
    small_coefficients: tuple[float, ...]
    large_coefficients: tuple[float, ...]
    harmonics: int = HARMONICS

    _small_offset: float = field(init=False, repr=False, compare=False)
    _small_scale: float = field(init=False, repr=False, compare=False)
    _large_offset: float = field(init=False, repr=False, compare=False)
    _large_scale: float = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        expected = 2 * self.harmonics
        if len(self.small_coefficients) != expected or len(self.large_coefficients) != expected:
            raise ValueError(f"Expected {expected} coefficients per piston.")
        for coeffs in (self.small_coefficients, self.large_coefficients):
            if not np.all(np.isfinite(np.asarray(coeffs, float))):
                raise ValueError("Fourier coefficients must be finite.")
            if np.linalg.norm(np.asarray(coeffs, float)) < 1e-12:
                raise ValueError("Fourier coefficient vector must be non-zero.")

        so, ss = self._normalization(self.small_coefficients)
        lo, ls = self._normalization(self.large_coefficients)
        object.__setattr__(self, "_small_offset", so)
        object.__setattr__(self, "_small_scale", ss)
        object.__setattr__(self, "_large_offset", lo)
        object.__setattr__(self, "_large_scale", ls)

    @property
    def small_physical_stroke(self):
        return None

    @property
    def large_physical_stroke(self):
        return None

    def breakpoint_angles(self):
        # No derivative discontinuity exists.
        return ()

    def _raw(self, t, coeffs):
        t = np.asarray(t, dtype=float)
        out = np.zeros_like(t)
        for k in range(1, self.harmonics + 1):
            a = coeffs[2 * (k - 1)]
            b = coeffs[2 * (k - 1) + 1]
            w = 2.0 * math.pi * k
            out = out + a * np.cos(w * t) + b * np.sin(w * t)
        return out

    def _raw_dt(self, t, coeffs):
        t = np.asarray(t, dtype=float)
        out = np.zeros_like(t)
        for k in range(1, self.harmonics + 1):
            a = coeffs[2 * (k - 1)]
            b = coeffs[2 * (k - 1) + 1]
            w = 2.0 * math.pi * k
            out = out + (-a * w * np.sin(w * t) + b * w * np.cos(w * t))
        return out

    def _normalization(self, coeffs):
        coeffs = tuple(float(x) for x in coeffs)
        grid = np.linspace(0.0, 1.0, 2049, endpoint=False)
        raw = self._raw(grid, coeffs)

        i_min = int(np.argmin(raw))
        i_max = int(np.argmax(raw))
        dx = 1.0 / len(grid)

        def wrapped(x):
            return float(self._raw(np.asarray([x % 1.0]), coeffs)[0])

        def local_extreme(index, sign):
            center = grid[index]
            result = minimize_scalar(
                lambda x: sign * wrapped(x),
                bounds=(center - 2.0 * dx, center + 2.0 * dx),
                method="bounded",
                options={"xatol": 1e-13},
            )
            return wrapped(result.x)

        minimum = local_extreme(i_min, +1.0)
        maximum = local_extreme(i_max, -1.0)
        span = maximum - minimum
        if not math.isfinite(span) or span <= 1e-9:
            raise ValueError("Degenerate Fourier motion.")
        return minimum, 1.0 / span

    def normalized_fractions(self, t):
        sr = self._raw(t, self.small_coefficients)
        lr = self._raw(t, self.large_coefficients)
        s = (sr - self._small_offset) * self._small_scale
        l = (lr - self._large_offset) * self._large_scale
        return s, l

    def normalized_fraction_derivatives(self, t):
        sd = self._raw_dt(t, self.small_coefficients) * self._small_scale
        ld = self._raw_dt(t, self.large_coefficients) * self._large_scale
        return sd, ld

    def cylinder_volumes_and_derivatives(self, theta):
        if not math.isfinite(theta):
            raise ValueError("Angle must be finite.")
        t = (-theta / (2.0 * math.pi)) % 1.0
        s, l = self.normalized_fractions(t)
        sd, ld = self.normalized_fraction_derivatives(t)

        s = float(s)
        l = float(l)
        sd = float(sd)
        ld = float(ld)

        sv = self.small_volume_limits.minimum + s * self.small_volume_limits.swept
        lv = self.large_volume_limits.minimum + l * self.large_volume_limits.swept
        # dt/dtheta = -1/(2*pi)
        svd = -sd * self.small_volume_limits.swept / (2.0 * math.pi)
        lvd = -ld * self.large_volume_limits.swept / (2.0 * math.pi)
        return sv, lv, svd, lvd

    def small_cylinder_volume(self, theta):
        return self.cylinder_volumes_and_derivatives(theta)[0]

    def large_cylinder_volume(self, theta):
        return self.cylinder_volumes_and_derivatives(theta)[1]

    def small_cylinder_volume_derivative(self, theta):
        return self.cylinder_volumes_and_derivatives(theta)[2]

    def large_cylinder_volume_derivative(self, theta):
        return self.cylinder_volumes_and_derivatives(theta)[3]


def _linear_fraction(t, p, piston):
    knots = np.asarray([0.0, p["t1"], p["t2"], p["t3"], 1.0])
    if piston == "small":
        values = np.asarray([p["b_s"], 1.0, p["a_s"], 0.0, p["b_s"]])
    else:
        values = np.asarray([1.0, p["b_l"], 0.0, p["a_l"], 1.0])
    return np.interp(np.asarray(t), knots, values)


def _fit_fourier_coefficients(p, piston):
    """Least-squares projection of the retained linear law."""
    t = np.linspace(0.0, 1.0, 4096, endpoint=False)
    y = _linear_fraction(t, p, piston)
    # Constant offset is irrelevant because every candidate is normalized.
    y = y - np.mean(y)
    columns = []
    for k in range(1, HARMONICS + 1):
        w = 2.0 * math.pi * k
        columns.append(np.cos(w * t))
        columns.append(np.sin(w * t))
    A = np.column_stack(columns)
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    n = float(np.linalg.norm(c))
    if n <= 0:
        raise RuntimeError("Could not project linear motion to Fourier basis.")
    return c / n


def _count_extrema(coeffs):
    """Count robust derivative sign changes over one period."""
    dummy_limits = CylinderVolumeLimits(minimum=1e-6, maximum=2e-6)
    kin = FourierVolumeKinematics(
        dummy_limits,
        dummy_limits,
        tuple(coeffs),
        tuple(coeffs),
    )
    t = np.linspace(0.0, 1.0, 4096, endpoint=False)
    d = np.asarray(kin.normalized_fraction_derivatives(t)[0])
    threshold = max(1e-8, 1e-4 * float(np.max(np.abs(d))))
    signs = np.sign(d)
    signs[np.abs(d) < threshold] = 0

    nz = np.flatnonzero(signs)
    if len(nz) == 0:
        return 0
    compact = signs[nz]
    changes = int(np.sum(compact[1:] != compact[:-1]))
    changes += int(compact[0] != compact[-1])
    return changes


def _shape_is_admissible(small, large):
    # One maximum and one minimum per revolution => two derivative sign changes.
    return _count_extrema(small) == 2 and _count_extrema(large) == 2


def _normalize_coefficients(v):
    v = np.asarray(v, float)
    n = float(np.linalg.norm(v))
    if not math.isfinite(n) or n < 1e-12:
        raise ValueError("Degenerate coefficient proposal.")
    return v / n


def _shape_distance(a, b):
    # Account for the irrelevant global sign: c and -c normalize to mirrored
    # raw extrema but produce q and 1-q.  They are not equivalent physically,
    # so keep ordinary Euclidean distance.
    return float(np.linalg.norm(np.asarray(a, float) - np.asarray(b, float)))


def _candidate_id(parameters, small, large):
    payload = {
        **{
            k: int(parameters[k]) if k in ("n_i", "n_o") else float(parameters[k])
            for k in THERMO_NAMES
        },
        "small_coefficients": [float(x) for x in small],
        "large_coefficients": [float(x) for x in large],
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _thermo_bounds(seed):
    bounds = {}
    for name in ("swept_ratio", "length_i_m", "length_o_m"):
        x = float(seed[name])
        bounds[name] = (x * (1.0 - THERMO_RELATIVE_LIMIT), x * (1.0 + THERMO_RELATIVE_LIMIT))
    for name in ("n_i", "n_o"):
        x = int(seed[name])
        lo = int(math.floor(x * (1.0 - THERMO_RELATIVE_LIMIT)))
        hi = int(math.ceil(x * (1.0 + THERMO_RELATIVE_LIMIT)))
        bounds[name] = (lo, hi)
    return bounds


def _thermo_valid(p, bounds):
    return all(bounds[name][0] <= p[name] <= bounds[name][1] for name in THERMO_NAMES)


def _propose_thermo(center, seed, bounds, u, global_fraction=0.0):
    """Small thermo motion; hard bounds remain +/-3 % from the UU seed."""
    p = dict(center)
    # Radius is intentionally much smaller than the hard envelope.
    rel = 0.008 if global_fraction == 0.0 else 0.015
    for i, name in enumerate(("swept_ratio", "length_i_m", "length_o_m")):
        value = float(center[name]) * (1.0 + (2.0 * u[i] - 1.0) * rel)
        p[name] = float(np.clip(value, *bounds[name]))
    for j, name in enumerate(("n_i", "n_o"), start=3):
        radius = max(8, int(round(seed[name] * (0.008 if global_fraction == 0.0 else 0.015))))
        value = int(round(center[name] + (2.0 * u[j] - 1.0) * radius))
        p[name] = int(np.clip(value, *bounds[name]))
    return p


def _build_design(base, geometry, thermo, small, large):
    # Reuse the production 260 K hardware construction, then replace only motion.
    placeholder = dict(
        thermo,
        t1=0.25,
        t2=0.50,
        t3=0.75,
        a_l=0.5,
        b_l=0.5,
        a_s=0.5,
        b_s=0.5,
    )
    design = _linear_build_design(base, geometry, placeholder, HOT_K)
    limits = design.configuration.machine_volumes
    kin = FourierVolumeKinematics(
        limits.small_cylinder,
        limits.large_cylinder,
        tuple(float(x) for x in small),
        tuple(float(x) for x in large),
    )
    return replace(
        design,
        configuration=replace(
            design.configuration,
            heat_in_valve_placement="upstream",
            heat_out_valve_placement="upstream",
        ),
        kinematics=kin,
    )


def _save_motion_csv(path, small, large):
    # Normalization does not depend on physical volume limits.
    lim = CylinderVolumeLimits(minimum=1e-6, maximum=2e-6)
    kin = FourierVolumeKinematics(lim, lim, tuple(small), tuple(large))
    t = np.linspace(0.0, 1.0, 1441)
    s, l = kin.normalized_fractions(t)
    sd, ld = kin.normalized_fraction_derivatives(t)
    arr = np.column_stack((t, 360.0 * t, s, l, sd, ld))
    np.savetxt(
        path,
        arr,
        delimiter=",",
        header="motor_time_fraction,motor_angle_deg,small_fraction_0_1,large_fraction_0_1,dsmall_dt,dlarge_dt",
        comments="",
    )


def _load_history(path):
    if not path.exists():
        return []
    return [json.loads(x) for x in path.read_text().splitlines() if x.strip()]


def _eligible(r):
    result = r.get("result", {})
    return bool(
        r.get("feasible")
        and result.get("status") == "converged"
        and result.get("indicated_thermal_efficiency") is not None
    )


def _efficiency(r):
    return float(r["result"]["indicated_thermal_efficiency"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--budget-seconds", type=float, default=21600.0)
    ap.add_argument("--evaluations", type=int, default=512)
    ap.add_argument("--candidate-seconds", type=float, default=180.0)
    ap.add_argument("--evaluations-per-radius", type=int, default=DEFAULT_EVALUATIONS_PER_RADIUS)
    ap.add_argument("--seed", type=int, default=SOBOL_SEED)
    ap.add_argument("--output-directory", type=Path, default=OUTPUT_DIRECTORY)
    args = ap.parse_args()

    if min(args.budget_seconds, args.evaluations, args.candidate_seconds, args.evaluations_per_radius) <= 0:
        raise ValueError("Budgets and evaluation counts must be positive.")
    if not SOURCE_REPORT.exists():
        raise FileNotFoundError(SOURCE_REPORT)

    campaign_definition, base = _load_basis()
    definition = SimpleNamespace(
        wall_numerical_settings=campaign_definition.wall_numerical_settings,
        wall_backend=WallBackendSettings("numba"),
        exact_kinematics_cache=True,
        shared_replay=True,
    )

    geometry = base_geometry(base)
    source_p, source_state, source_record = load_source(SOURCE_REPORT, base, geometry)
    thermo_seed = {name: source_p[name] for name in THERMO_NAMES}
    bounds = _thermo_bounds(thermo_seed)

    small_seed = _fit_fourier_coefficients(source_p, "small")
    large_seed = _fit_fourier_coefficients(source_p, "large")
    if not _shape_is_admissible(small_seed, large_seed):
        raise RuntimeError(
            "The 4-harmonic projection of the linear champion does not have "
            "exactly one maximum and one minimum per piston."
        )

    out = args.output_directory
    if not out.is_absolute():
        out = ROOT / out
    out.mkdir(parents=True, exist_ok=True)
    hp = out / "history.jsonl"
    rp = out / "report.json"
    dp = out / "definition.json"

    identity = {
        "study": "260 K UU smooth Fourier-motion search",
        "source_report": str(SOURCE_REPORT.relative_to(ROOT)),
        "source_sha256": hashlib.sha256(SOURCE_REPORT.read_bytes()).hexdigest(),
        "delta_t_k": DELTA_T_K,
        "topology": "UU",
        "harmonics": HARMONICS,
        "shape_dimensions": SHAPE_DIM,
        "thermo_dimensions": len(THERMO_NAMES),
        "thermo_hard_relative_limit": THERMO_RELATIVE_LIMIT,
        "thermo_bounds": {k: list(v) for k, v in bounds.items()},
        "shape_radii": list(DEFAULT_SHAPE_RADII),
        "evaluations_per_radius": args.evaluations_per_radius,
        "sobol_seed": args.seed,
        "linear_source_parameters": source_p,
        "small_seed_coefficients": small_seed.tolist(),
        "large_seed_coefficients": large_seed.tolist(),
        "kinematic_constraint": "exactly one maximum and one minimum per piston per revolution",
        "normalization": "each Fourier curve is numerically normalized to full [0,1] stroke",
        "strategy": (
            "Seed with least-squares Fourier projection of the retained linear UU champion; "
            "then alternating broad incumbent-centered and occasional global smooth-shape proposals. "
            "Thermodynamics remains within +/-3% of the retained UU champion."
        ),
    }

    if dp.exists() and json.loads(dp.read_text()) != identity:
        raise RuntimeError("Campaign definition changed; use a new output directory.")
    if not dp.exists():
        dp.write_text(json.dumps(identity, indent=2) + "\n")

    history = _load_history(hp)
    eligible = [r for r in history if _eligible(r)]
    best = max(eligible, key=_efficiency) if eligible else None
    existing = {
        r["candidate_id"]
        for r in history
        if r.get("result", {}).get("status") != "interrupted"
    }

    # 5 thermo + 16 shape Sobol dimensions.
    points = qmc.Sobol(d=5 + SHAPE_DIM, scramble=True, seed=args.seed).random_base2(15)

    start = time.monotonic()
    deadline = start + args.budget_seconds
    completed = 0
    proposal_index = max((r.get("proposal_index", -1) for r in history), default=-1) + 1

    source_hw = source_record.get("hardware")
    if source_hw is None:
        # Reconstruct exact seed hardware metadata if an older report omitted it.
        seed_design = _build_design(base, geometry, thermo_seed, small_seed, large_seed)
        _, source_hw, _ = uniform_metadata(seed_design)

    while completed < args.evaluations and time.monotonic() < deadline:
        idx = proposal_index
        proposal_index += 1

        if not history and idx == 0:
            thermo = dict(thermo_seed)
            small = small_seed.copy()
            large = large_seed.copy()
            kind = "linear_fourier_projection_seed"
            radius_index = -1
            radius = 0.0
        else:
            n_nonseed = sum(1 for r in history if r.get("kind") != "linear_fourier_projection_seed")
            radius_index = min(n_nonseed // args.evaluations_per_radius, len(DEFAULT_SHAPE_RADII) - 1)
            radius = DEFAULT_SHAPE_RADII[radius_index]
            u = points[idx % len(points)]

            if best is not None:
                center_thermo = dict(best["thermo"])
                center_small = np.asarray(best["small_coefficients"], float)
                center_large = np.asarray(best["large_coefficients"], float)
            else:
                center_thermo = dict(thermo_seed)
                center_small = small_seed
                center_large = large_seed

            # Every fifth non-seed proposal explores shape space much more widely.
            global_shape = (idx % 5 == 0)
            if global_shape:
                sv = (2.0 * u[5:5 + 2 * HARMONICS] - 1.0)
                lv = (2.0 * u[5 + 2 * HARMONICS:] - 1.0)
                # Suppress pathological high-frequency dominance while still
                # allowing a shape unrelated to the incumbent.
                weights = np.repeat(1.0 / (np.arange(1, HARMONICS + 1) ** 1.5), 2)
                small = _normalize_coefficients(sv * weights)
                large = _normalize_coefficients(lv * weights)
                thermo = _propose_thermo(center_thermo, thermo_seed, bounds, u[:5], global_fraction=1.0)
                kind = "global_smooth_shape"
            else:
                small = _normalize_coefficients(
                    center_small + radius * (2.0 * u[5:5 + 2 * HARMONICS] - 1.0)
                )
                large = _normalize_coefficients(
                    center_large + radius * (2.0 * u[5 + 2 * HARMONICS:] - 1.0)
                )
                thermo = _propose_thermo(center_thermo, thermo_seed, bounds, u[:5])
                kind = "incumbent_shape_neighborhood"

            if not _thermo_valid(thermo, bounds):
                continue
            if not _shape_is_admissible(small, large):
                continue

        cid = _candidate_id(thermo, small, large)
        if cid in existing:
            continue

        before = time.monotonic()
        design = _build_design(base, geometry, thermo, small, large)

        try:
            uniform, hardware, derived_geometry = uniform_metadata(design)

            reusable = [
                r for r in history
                if r.get("result", {}).get("status") == "converged"
                and r.get("last_complete_state") is not None
                and r.get("hardware") is not None
            ]

            warm_source = None
            initial = np.asarray(source_state, float)

            if reusable:
                def distance(r):
                    dt = 0.0
                    for name in ("swept_ratio", "length_i_m", "length_o_m"):
                        dt += ((float(thermo[name]) - float(r["thermo"][name])) / float(thermo_seed[name])) ** 2
                    for name in ("n_i", "n_o"):
                        dt += ((float(thermo[name]) - float(r["thermo"][name])) / float(thermo_seed[name])) ** 2
                    ds = _shape_distance(small, r["small_coefficients"])
                    dl = _shape_distance(large, r["large_coefficients"])
                    return math.sqrt(dt + ds * ds + dl * dl)

                warm_source = min(reusable, key=distance)
                initial = warm_start(
                    warm_source["last_complete_state"],
                    warm_source["hardware"],
                    uniform,
                    hardware,
                )
            else:
                initial = warm_start(source_state, source_hw, uniform, hardware)

            candidate_deadline = min(deadline, before + args.candidate_seconds)

            def progress(_):
                now = time.monotonic()
                if now >= deadline:
                    raise IntegrationInterrupted("Fourier C2 campaign wall-clock budget exhausted.")
                if now >= candidate_deadline:
                    raise IntegrationInterrupted("Candidate wall-clock budget exhausted.")

            try:
                result, state = _evaluate(
                    f"fourier_c2_260k_{idx}",
                    design,
                    definition,
                    initial_state=np.asarray(initial),
                    progress_callback=progress,
                )
            except IntegrationInterrupted as exc:
                result = {"status": "interrupted", "message": str(exc)}
                state = None
            except (ValueError, RuntimeError) as exc:
                result = {"status": "integration_failure", "message": str(exc)}
                state = None

        except (ValueError, RuntimeError) as exc:
            result = {"status": "setup_failure", "message": str(exc)}
            state = None
            hardware = None
            derived_geometry = None

        feasible, reasons = (
            feasibility(result, POWER_FLOOR_W)
            if result.get("status") == "converged"
            else (False, [])
        )

        record = {
            "index": len(history),
            "proposal_index": idx,
            "candidate_id": cid,
            "kind": kind,
            "radius_index": radius_index,
            "shape_radius": radius,
            "thermo": {
                k: int(thermo[k]) if k in ("n_i", "n_o") else float(thermo[k])
                for k in THERMO_NAMES
            },
            "small_coefficients": [float(x) for x in small],
            "large_coefficients": [float(x) for x in large],
            "small_extrema_count": _count_extrema(small),
            "large_extrema_count": _count_extrema(large),
            "hardware": hardware,
            "derived_geometry": derived_geometry,
            "result": result,
            "feasible": feasible,
            "physical_constraint_failures": reasons,
            "elapsed_seconds": time.monotonic() - before,
            "warm_start_source": (
                warm_source["candidate_id"] if warm_source is not None else "UU_260K_reference"
            ),
            "last_complete_state": state.tolist() if state is not None else None,
        }

        with hp.open("a") as stream:
            stream.write(json.dumps(record) + "\n")
            stream.flush()

        history.append(record)
        completed += 1
        if result.get("status") != "interrupted":
            existing.add(cid)

        moved = False
        if feasible and (best is None or _efficiency(record) > _efficiency(best)):
            best = record
            moved = True
            _save_motion_csv(out / "best_motion.csv", small, large)

        report = {
            "best_feasible": best,
            "seed_projection": next(
                (r for r in history if r.get("kind") == "linear_fourier_projection_seed"),
                None,
            ),
            "reference_linear_UU": {
                "efficiency": source_record.get("result", {}).get("indicated_thermal_efficiency"),
                "power_W": source_record.get("result", {}).get("indicated_power_w"),
                "parameters": source_p,
            },
            "attempted_total": len(history),
            "completed_this_run": completed,
            "requested_duration_seconds": args.budget_seconds,
            "actual_duration_seconds": time.monotonic() - start,
            "definition": identity,
        }
        rp.write_text(json.dumps(report, indent=2) + "\n")

        print(
            json.dumps(
                {
                    "index": record["index"],
                    "proposal_index": idx,
                    "kind": kind,
                    "radius": radius,
                    "status": result.get("status"),
                    "feasible": feasible,
                    "efficiency": result.get("indicated_thermal_efficiency"),
                    "power_W": result.get("indicated_power_w"),
                    "center_moved": moved,
                    "best_efficiency": _efficiency(best) if best is not None else None,
                    "thermo": record["thermo"],
                    "elapsed_seconds": record["elapsed_seconds"],
                }
            ),
            flush=True,
        )

    print(f"Saved {rp}", flush=True)


if __name__ == "__main__":
    main()
