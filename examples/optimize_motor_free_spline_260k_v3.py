#!/usr/bin/env python3
"""Independent C2 free-spline motion search seeded from the final 8H champion.

Scientific purpose
------------------
This campaign tests parameterization independence.  It freezes the complete
thermodynamic/hardware design and total gas inventory of the supplied Fourier
champion, then changes ONLY the imposed cylinder motion.

Each piston is represented by:
  * an independent periodic cubic spline (production FreeKinematics);
  * an independent continuous phase shift, so extrema are explicitly free;
  * an admissibility condition of exactly one minimum and one maximum per cycle.

The spline controls are local in angle and structurally unrelated to a Fourier
basis.  Candidates with extra extrema are rejected before any thermodynamic
integration.  No four-stage boundaries, BP/HP timing, compression timing or
expansion timing are imposed.

Typical run after the 8H campaign has finished:

    PYTHONPATH=src python3 examples/optimize_motor_free_spline_260k.py \
        --source-report outputs/motor_fourier_c2_8h_refine/report.json \
        --output-directory outputs/motor_free_spline_260k \
        --controls 16 \
        --evaluations 1024 \
        --budget-seconds 10800

The run is append-only and resumable.

Outputs
-------
    <output>/definition.json
    <output>/history.jsonl
    <output>/report.json
    <output>/best_motion.csv

Notes
-----
* Hardware, swept volumes, clearances, exchanger geometry, valve placement,
  frequency, source temperatures and TOTAL GAS MASS are frozen to the source
  Fourier champion.
* The source 8H law is sampled at the spline knots to create the initial seed.
* The production FreeKinematics normalization uses the true continuous spline
  extrema, not sampled clipping.
* Phase shifts act in the injected study-angle convention.  The CSV is written
  in forward motor-time convention to match the Fourier optimization reports.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from types import SimpleNamespace
import argparse
import hashlib
import json
import math
import time

import numpy as np
from scipy.stats import qmc

from dada_solver.configuration import ChargeConfiguration
from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
from dada_solver.free_kinematics import (
    FreeKinematics,
    FreeKinematicsConfiguration,
    FreeMotionDefinition,
)
from dada_solver.integration import IntegrationInterrupted
from dada_solver.wall_backend import WallBackendSettings

from compare_motor_motion_laws_stage7A5 import ROOT, _evaluate, _uniform_wall_initial
from optimize_motor_temperature_point import base_geometry, feasibility
from refine_motor_four_stage_hx9d_variable_gas import _load_basis
from refine_motor_fourier_c2_260k import build_design as fourier_build_design


DEFAULT_SOURCE = ROOT / "outputs" / "motor_fourier_c2_8h_refine" / "report.json"
DEFAULT_OUTPUT = ROOT / "outputs" / "motor_free_spline_260k"

POWER_FLOOR_W = 25.0
DEFAULT_CONTROLS = 16
DEFAULT_SEED = 26092231
EVALS_PER_RADIUS = 96

# Euclidean displacement of each canonical zero-mean/unit-norm control vector.
SHAPE_RADII = (0.18, 0.12, 0.080, 0.050, 0.032, 0.020, 0.012, 0.007, 0.0035)

# Independent phase-shift radius for each piston at the corresponding stage.
PHASE_RADII_DEG = (10.0, 7.0, 5.0, 3.5, 2.4, 1.5, 0.9, 0.5, 0.25)


@dataclass(frozen=True)
class PhaseShiftedFreeKinematics:
    """Production FreeKinematics with independent piston phase shifts.

    The generic FreeKinematics API remains authoritative.  The combined scalar
    evaluator below merely evaluates its already-built cubic coefficients
    directly, avoiding NumPy allocation on every LSODA RHS call.
    """

    base: FreeKinematics
    small_phase_rad: float = 0.0
    large_phase_rad: float = 0.0

    @property
    def small_volume_limits(self):
        return self.base.small_volume_limits

    @property
    def large_volume_limits(self):
        return self.base.large_volume_limits

    @property
    def small_physical_stroke(self):
        return None

    @property
    def large_physical_stroke(self):
        return None

    @property
    def diagnostics(self):
        return self.base.diagnostics

    def require_feasible(self):
        self.base.require_feasible()

    @staticmethod
    def _scalar_pair(motion, theta):
        """Return volume and dV/dtheta from production spline coefficients."""
        wrapped = float(theta) % (2.0 * math.pi)
        knots = motion._knots
        n = len(knots) - 1
        # Equally spaced knots by FreeKinematics construction.
        index = int(wrapped * n / (2.0 * math.pi))
        if index >= n:
            index = n - 1
        dx = wrapped - knots[index]
        coeff = motion._coefficients
        a = coeff[0][index]
        b = coeff[1][index]
        c = coeff[2][index]
        d = coeff[3][index]
        raw = ((a * dx + b) * dx + c) * dx + d
        draw = (3.0 * a * dx + 2.0 * b) * dx + c
        volume = (
            (raw - motion._minimum) * motion._scale
            + motion.definition.minimum_volume
        )
        derivative = draw * motion._scale
        return volume, derivative

    def small_cylinder_volume(self, theta):
        if np.ndim(theta):
            return self.base.small_cylinder_volume(
                np.asarray(theta) - self.small_phase_rad
            )
        return self._scalar_pair(
            self.base._small, float(theta) - self.small_phase_rad
        )[0]

    def large_cylinder_volume(self, theta):
        if np.ndim(theta):
            return self.base.large_cylinder_volume(
                np.asarray(theta) - self.large_phase_rad
            )
        return self._scalar_pair(
            self.base._large, float(theta) - self.large_phase_rad
        )[0]

    def small_cylinder_volume_derivative(self, theta):
        if np.ndim(theta):
            return self.base.small_cylinder_volume_derivative(
                np.asarray(theta) - self.small_phase_rad
            )
        return self._scalar_pair(
            self.base._small, float(theta) - self.small_phase_rad
        )[1]

    def large_cylinder_volume_derivative(self, theta):
        if np.ndim(theta):
            return self.base.large_cylinder_volume_derivative(
                np.asarray(theta) - self.large_phase_rad
            )
        return self._scalar_pair(
            self.base._large, float(theta) - self.large_phase_rad
        )[1]

    def small_cylinder_volume_second_derivative(self, theta):
        return self.base.small_cylinder_volume_second_derivative(
            np.asarray(theta) - self.small_phase_rad
        )

    def large_cylinder_volume_second_derivative(self, theta):
        return self.base.large_cylinder_volume_second_derivative(
            np.asarray(theta) - self.large_phase_rad
        )

    def cylinder_volumes_and_derivatives(self, theta):
        # LSODA / CompiledWallRHS uses scalar angles here.
        if np.ndim(theta):
            return (
                self.small_cylinder_volume(theta),
                self.large_cylinder_volume(theta),
                self.small_cylinder_volume_derivative(theta),
                self.large_cylinder_volume_derivative(theta),
            )
        sv, ds = self._scalar_pair(
            self.base._small, float(theta) - self.small_phase_rad
        )
        lv, dl = self._scalar_pair(
            self.base._large, float(theta) - self.large_phase_rad
        )
        return sv, lv, ds, dl

    def breakpoint_angles(self):
        return ()


def _source_best(report: dict) -> dict:
    best = report.get("best_feasible")
    if best is None:
        raise RuntimeError("Source report has no best_feasible.")
    if best.get("last_complete_state") is None:
        raise RuntimeError("Source champion has no saved periodic state.")
    if "small_coefficients" not in best or "large_coefficients" not in best:
        raise RuntimeError("Source champion is not a Fourier motion record.")
    return best


def _source_harmonics(best: dict) -> int:
    h = int(best.get("harmonics", len(best["small_coefficients"]) // 2))
    if len(best["small_coefficients"]) != 2 * h:
        raise RuntimeError("Invalid small Fourier coefficient count.")
    if len(best["large_coefficients"]) != 2 * h:
        raise RuntimeError("Invalid large Fourier coefficient count.")
    return h


def _fixed_inventory_design(design, total_mass_kg: float):
    charge = ChargeConfiguration(
        temperature=design.configuration.charge.temperature,
        total_mass=float(total_mass_kg),
    )
    return replace(
        design,
        configuration=replace(design.configuration, charge=charge),
    )


def _canonical(values):
    values = np.asarray(values, dtype=float)
    values = values - np.mean(values)
    norm = float(np.linalg.norm(values))
    if not math.isfinite(norm) or norm < 1e-12:
        raise ValueError("Degenerate spline control vector.")
    return values / norm


def _project_fourier_to_controls(source_kinematics, controls: int):
    """Sample both study-angle laws at equally spaced production spline knots."""
    theta = 2.0 * math.pi * np.arange(controls, dtype=float) / controls

    sl = source_kinematics.small_volume_limits
    ll = source_kinematics.large_volume_limits

    small = np.asarray(
        [
            (source_kinematics.small_cylinder_volume(float(a)) - sl.minimum)
            / sl.swept
            for a in theta
        ]
    )
    large = np.asarray(
        [
            (source_kinematics.large_cylinder_volume(float(a)) - ll.minimum)
            / ll.swept
            for a in theta
        ]
    )
    return _canonical(small), _canonical(large)


def _make_kinematics(limits, small, large, small_phase_deg, large_phase_deg):
    sm = FreeMotionDefinition(
        tuple(map(float, small)),
        limits.small_cylinder.minimum,
        limits.small_cylinder.maximum,
    )
    lg = FreeMotionDefinition(
        tuple(map(float, large)),
        limits.large_cylinder.minimum,
        limits.large_cylinder.maximum,
    )
    base = FreeKinematics(FreeKinematicsConfiguration(sm, lg))
    return PhaseShiftedFreeKinematics(
        base,
        math.radians(float(small_phase_deg)),
        math.radians(float(large_phase_deg)),
    )


def _build_candidate_design(source_design, small, large, small_phase_deg, large_phase_deg):
    kin = _make_kinematics(
        source_design.configuration.machine_volumes,
        small,
        large,
        small_phase_deg,
        large_phase_deg,
    )
    return replace(source_design, kinematics=kin)


def _normalized_motion(kinematics, motor_fraction):
    """Evaluate motion in forward motor convention.

    FreeKinematics accepts vector angles and uses the fast vectorized path.
    Historical FourierVolumeKinematics is scalar-only, so retain a compatible
    scalar fallback for source/reference motions.
    """
    t = np.asarray(motor_fraction, dtype=float)
    theta = -2.0 * math.pi * t

    sl = kinematics.small_volume_limits
    ll = kinematics.large_volume_limits

    def evaluate(name):
        fn = getattr(kinematics, name)
        try:
            values = np.asarray(fn(theta), dtype=float)
            if values.shape == theta.shape:
                return values
        except (TypeError, ValueError):
            pass
        return np.asarray([fn(float(a)) for a in theta], dtype=float)

    small_v = evaluate("small_cylinder_volume")
    large_v = evaluate("large_cylinder_volume")
    small_d = evaluate("small_cylinder_volume_derivative")
    large_d = evaluate("large_cylinder_volume_derivative")

    small = (small_v - sl.minimum) / sl.swept
    large = (large_v - ll.minimum) / ll.swept
    ds = -2.0 * math.pi * small_d / sl.swept
    dl = -2.0 * math.pi * large_d / ll.swept
    return np.vstack((small, large)), np.vstack((ds, dl))


def _robust_extrema_count(derivative):
    d = np.asarray(derivative, dtype=float)
    threshold = max(1e-8, 1e-4 * float(np.max(np.abs(d))))
    signs = np.sign(d)
    signs[np.abs(d) < threshold] = 0.0
    signs = signs[signs != 0.0]
    if len(signs) == 0:
        return 0
    return int(np.sum(signs[1:] != signs[:-1]) + (signs[0] != signs[-1]))


def _shape_diagnostics(kinematics, sample_count=8192):
    t = np.linspace(0.0, 1.0, sample_count, endpoint=False)
    q, dq = _normalized_motion(kinematics, t)
    return {
        "small_max_deg": 360.0 * float(t[int(np.argmax(q[0]))]),
        "small_min_deg": 360.0 * float(t[int(np.argmin(q[0]))]),
        "large_max_deg": 360.0 * float(t[int(np.argmax(q[1]))]),
        "large_min_deg": 360.0 * float(t[int(np.argmin(q[1]))]),
        "small_peak_abs_dq_dt": float(np.max(np.abs(dq[0]))),
        "large_peak_abs_dq_dt": float(np.max(np.abs(dq[1]))),
        "small_extrema_count": _robust_extrema_count(dq[0]),
        "large_extrema_count": _robust_extrema_count(dq[1]),
    }


def _admissible(kinematics):
    d = _shape_diagnostics(kinematics, 2048)
    return d["small_extrema_count"] == 2 and d["large_extrema_count"] == 2


def _shape_rms(reference_kinematics, test_kinematics, sample_count=4096):
    t = np.linspace(0.0, 1.0, sample_count, endpoint=False)
    ref, dref = _normalized_motion(reference_kinematics, t)
    test, dtest = _normalized_motion(test_kinematics, t)
    delta = test - ref
    ddelta = dtest - dref
    return {
        "small_position_rms": float(np.sqrt(np.mean(delta[0] ** 2))),
        "large_position_rms": float(np.sqrt(np.mean(delta[1] ** 2))),
        "combined_position_rms": float(np.sqrt(np.mean(delta ** 2))),
        "small_velocity_rms": float(np.sqrt(np.mean(ddelta[0] ** 2))),
        "large_velocity_rms": float(np.sqrt(np.mean(ddelta[1] ** 2))),
        "combined_velocity_rms": float(np.sqrt(np.mean(ddelta ** 2))),
        "maximum_position_difference": float(np.max(np.abs(delta))),
    }


def _phase_difference_deg(a, b):
    return ((float(a) - float(b) + 180.0) % 360.0) - 180.0


def _candidate_id(small, large, small_phase_deg, large_phase_deg):
    payload = {
        "small": list(map(float, small)),
        "large": list(map(float, large)),
        "small_phase_deg": float(small_phase_deg) % 360.0,
        "large_phase_deg": float(large_phase_deg) % 360.0,
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _record_distance(small, large, small_phase_deg, large_phase_deg, record):
    ds = np.asarray(small) - np.asarray(record["small_controls"])
    dl = np.asarray(large) - np.asarray(record["large_controls"])
    # 10 degrees ~= 0.1 in canonical-control Euclidean distance.
    ps = _phase_difference_deg(small_phase_deg, record["small_phase_deg"]) / 100.0
    pl = _phase_difference_deg(large_phase_deg, record["large_phase_deg"]) / 100.0
    return float(math.sqrt(ds @ ds + dl @ dl + ps * ps + pl * pl))


# Above this distance a previous periodic gas distribution is more likely to
# trigger an immediate hydraulic-domain failure than to accelerate convergence.
DIRECT_PERIODIC_WARM_DISTANCE = 0.055


def _direction(u, center):
    d = 2.0 * np.asarray(u, dtype=float) - 1.0
    d -= np.mean(d)
    # Stay in the tangent space of the canonical sphere to avoid wasting the
    # perturbation on the irrelevant scale direction.
    d -= center * float(d @ center)
    n = float(np.linalg.norm(d))
    if n < 1e-12:
        # Deterministic fallback; this is effectively impossible for Sobol.
        d = np.roll(center, 1) - center
        d -= np.mean(d)
        d -= center * float(d @ center)
        n = float(np.linalg.norm(d))
    return d / n


def _propose(center, u, radius, phase_radius_deg, territory_probe=False):
    small0 = np.asarray(center["small_controls"], dtype=float)
    large0 = np.asarray(center["large_controls"], dtype=float)
    n = len(small0)

    ds = _direction(u[:n], small0)
    dl = _direction(u[n : 2 * n], large0)

    # The final two Sobol dimensions vary proposal magnitude continuously.
    rs = radius * (0.25 + 0.75 * float(u[2 * n]))
    rl = radius * (0.25 + 0.75 * float(u[2 * n + 1]))

    small = _canonical(small0 + rs * ds)
    large = _canonical(large0 + rl * dl)

    # Use two deterministic mixtures of existing Sobol coordinates for phases
    # without increasing dimension again.
    us = (float(u[0]) + float(u[n])) % 1.0
    ul = (float(u[n - 1]) + float(u[2 * n - 1])) % 1.0
    multiplier = 1.5 if territory_probe else 1.0

    small_phase = (
        float(center["small_phase_deg"])
        + multiplier * phase_radius_deg * (2.0 * us - 1.0)
    ) % 360.0
    large_phase = (
        float(center["large_phase_deg"])
        + multiplier * phase_radius_deg * (2.0 * ul - 1.0)
    ) % 360.0

    return small, large, small_phase, large_phase


def _safe_uniform_initial(design, wall_source=None):
    wrapper = design.build()
    state = np.asarray(
        _uniform_wall_initial(design.configuration, wrapper),
        dtype=float,
    )
    if wall_source is not None:
        wall_source = np.asarray(wall_source, dtype=float)
        if wall_source.shape == (10,):
            # Identical hardware: reuse wall energies, but deliberately rebuild
            # the gas charge for the target theta=0 geometry.
            state[8:10] = wall_source[8:10]
    return state


def _load_history(path):
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def _eligible(record):
    result = record.get("result", {})
    return bool(
        record.get("feasible")
        and result.get("status") == "converged"
        and result.get("indicated_thermal_efficiency") is not None
    )


def _eta(record):
    return float(record["result"]["indicated_thermal_efficiency"])


def _save_motion(path, source_kinematics, candidate_kinematics):
    t = np.linspace(0.0, 1.0, 2881)
    src, dsrc = _normalized_motion(source_kinematics, t)
    cur, dcur = _normalized_motion(candidate_kinematics, t)

    data = np.column_stack(
        (
            t,
            360.0 * t,
            src[0],
            src[1],
            cur[0],
            cur[1],
            dsrc[0],
            dsrc[1],
            dcur[0],
            dcur[1],
        )
    )
    np.savetxt(
        path,
        data,
        delimiter=",",
        header=(
            "motor_time_fraction,motor_angle_deg,"
            "source_fourier_small,source_fourier_large,"
            "spline_small,spline_large,"
            "source_fourier_dsmall_dt,source_fourier_dlarge_dt,"
            "spline_dsmall_dt,spline_dlarge_dt"
        ),
        comments="",
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-report", type=Path, default=DEFAULT_SOURCE)
    ap.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--controls", type=int, default=DEFAULT_CONTROLS)
    ap.add_argument("--evaluations", type=int, default=1024)
    ap.add_argument("--budget-seconds", type=float, default=10800.0)
    ap.add_argument("--candidate-seconds", type=float, default=180.0)
    ap.add_argument("--evaluations-per-radius", type=int, default=EVALS_PER_RADIUS)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = ap.parse_args()

    if args.controls < 8:
        raise ValueError("--controls must be at least 8 for this scientific search.")
    if any(
        x <= 0
        for x in (
            args.evaluations,
            args.budget_seconds,
            args.candidate_seconds,
            args.evaluations_per_radius,
        )
    ):
        raise ValueError("Evaluation counts and time budgets must be positive.")

    source_path = (
        args.source_report
        if args.source_report.is_absolute()
        else ROOT / args.source_report
    )
    source_report = json.loads(source_path.read_text())
    source_best = _source_best(source_report)
    source_h = _source_harmonics(source_best)

    campaign_definition, base = _load_basis()
    geometry = base_geometry(base)

    source_design = fourier_build_design(
        base,
        geometry,
        source_best["thermo"],
        np.asarray(source_best["small_coefficients"], dtype=float),
        np.asarray(source_best["large_coefficients"], dtype=float),
        source_h,
    )

    source_state = np.asarray(source_best["last_complete_state"], dtype=float)
    total_mass_kg = float(
        source_best.get("result", {}).get(
            "total_mass_kg", np.sum(source_state[:8:2])
        )
    )
    source_design = _fixed_inventory_design(source_design, total_mass_kg)
    source_kinematics = source_design.kinematics

    small_seed, large_seed = _project_fourier_to_controls(
        source_kinematics, args.controls
    )
    seed_kinematics = _make_kinematics(
        source_design.configuration.machine_volumes,
        small_seed,
        large_seed,
        0.0,
        0.0,
    )
    seed_shape = _shape_diagnostics(seed_kinematics)
    if (
        seed_shape["small_extrema_count"] != 2
        or seed_shape["large_extrema_count"] != 2
    ):
        raise RuntimeError(
            "Projected Fourier seed gained extra spline extrema. "
            "Increase --controls or inspect the source motion."
        )

    seed_rms = _shape_rms(source_kinematics, seed_kinematics)

    definition = SimpleNamespace(
        wall_numerical_settings=campaign_definition.wall_numerical_settings,
        wall_backend=WallBackendSettings("numba"),
        exact_kinematics_cache=True,
        shared_replay=True,
    )

    out = (
        args.output_directory
        if args.output_directory.is_absolute()
        else ROOT / args.output_directory
    )
    out.mkdir(parents=True, exist_ok=True)
    history_path = out / "history.jsonl"
    report_path = out / "report.json"
    definition_path = out / "definition.json"

    identity = {
        "study": "260 K fixed-machine independent free-spline motion search v2-fast",
        "source_report": str(source_path),
        "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "source_harmonics": source_h,
        "source_candidate_id": source_best.get("candidate_id"),
        "source_efficiency": source_best["result"]["indicated_thermal_efficiency"],
        "source_power_W": source_best["result"]["indicated_power_w"],
        "fixed_total_mass_kg": total_mass_kg,
        "fixed_thermo": source_best["thermo"],
        "controls_per_piston": args.controls,
        "effective_shape_dof_per_piston": args.controls - 2,
        "independent_phase_dof": 2,
        "shape_radii": list(SHAPE_RADII),
        "phase_radii_deg": list(PHASE_RADII_DEG),
        "evaluations_per_radius": args.evaluations_per_radius,
        "sobol_seed": args.seed,
        "projected_seed_rms_vs_source": seed_rms,
        "constraints": {
            "exactly_one_minimum_and_one_maximum_per_piston": True,
            "C2_periodic_position": True,
            "hardware_frozen": True,
            "gas_inventory_frozen": True,
            "fast_scalar_spline_rhs": True,
            "four_stage_boundaries_imposed": False,
            "velocity_limit_imposed": False,
            "acceleration_limit_imposed": False,
        },
        "strategy": (
            "Periodic cubic FreeKinematics; local tangent-space Sobol refinement "
            "with independent piston phase shifts. Every eighth evaluated "
            "proposal probes around the original projected Fourier seed rather "
            "than the incumbent."
        ),
    }

    if definition_path.exists():
        existing = json.loads(definition_path.read_text())
        if existing != identity:
            raise RuntimeError(
                "Campaign definition changed. Use a new output directory."
            )
    else:
        definition_path.write_text(json.dumps(identity, indent=2) + "\n")

    history = _load_history(history_path)
    valid = [r for r in history if _eligible(r)]
    best = max(valid, key=_eta) if valid else None
    finished = {
        r["candidate_id"]
        for r in history
        if r.get("result", {}).get("status") != "interrupted"
    }

    # Enough deterministic Sobol points for long resumes.
    dimension = 2 * args.controls + 2
    sobol = qmc.Sobol(d=dimension, scramble=True, seed=args.seed)
    points = sobol.random_base2(16)

    start = time.monotonic()
    deadline = start + args.budget_seconds
    completed = 0
    proposal_index = (
        max((r.get("proposal_index", -1) for r in history), default=-1) + 1
    )

    seed_record = {
        "small_controls": list(map(float, small_seed)),
        "large_controls": list(map(float, large_seed)),
        "small_phase_deg": 0.0,
        "large_phase_deg": 0.0,
    }

    while completed < args.evaluations and time.monotonic() < deadline:
        this_proposal = proposal_index
        proposal_index += 1

        if not history and this_proposal == 0:
            small = small_seed.copy()
            large = large_seed.copy()
            small_phase = 0.0
            large_phase = 0.0
            kind = "projected_fourier_seed"
            radius_index = -1
            radius = 0.0
            phase_radius = 0.0
        else:
            searched = sum(
                1
                for r in history
                if r.get("kind") != "projected_fourier_seed"
            )
            radius_index = min(
                searched // args.evaluations_per_radius,
                len(SHAPE_RADII) - 1,
            )
            radius = SHAPE_RADII[radius_index]
            phase_radius = PHASE_RADII_DEG[radius_index]
            u = points[this_proposal % len(points)]

            territory_probe = searched % 8 == 7
            if territory_probe:
                center = seed_record
                kind = "source_territory_probe"
            elif best is not None:
                center = best
                kind = "incumbent_spline_refinement"
            else:
                center = seed_record
                kind = "seed_spline_refinement"

            small, large, small_phase, large_phase = _propose(
                center,
                u,
                radius,
                phase_radius,
                territory_probe=territory_probe,
            )

        cid = _candidate_id(small, large, small_phase, large_phase)
        if cid in finished:
            continue

        candidate_design = _build_candidate_design(
            source_design,
            small,
            large,
            small_phase,
            large_phase,
        )
        candidate_kinematics = candidate_design.kinematics

        # Cheap structural rejection before any solver work.
        if not _admissible(candidate_kinematics):
            continue

        before = time.monotonic()
        candidate_deadline = min(deadline, before + args.candidate_seconds)

        reusable = [
            r
            for r in history
            if r.get("result", {}).get("status") == "converged"
            and r.get("last_complete_state") is not None
        ]
        warm = None
        warm_distance = None
        if reusable:
            warm = min(
                reusable,
                key=lambda r: _record_distance(
                    small, large, small_phase, large_phase, r
                ),
            )
            warm_distance = _record_distance(
                small, large, small_phase, large_phase, warm
            )
            warm_source = warm["candidate_id"]
            if warm_distance <= DIRECT_PERIODIC_WARM_DISTANCE:
                initial = np.asarray(warm["last_complete_state"], dtype=float)
                warm_mode = "nearest_periodic_state"
            else:
                initial = _safe_uniform_initial(
                    candidate_design,
                    wall_source=np.asarray(warm["last_complete_state"], dtype=float),
                )
                warm_mode = "direct_safe_uniform_for_distant_shape"
        else:
            # The projected seed is extremely close to the Fourier source and
            # benefits from its converged state.  Later distant proposals use
            # the branch above.
            initial = source_state.copy()
            warm_source = source_best.get("candidate_id", "source_fourier")
            warm_mode = "source_fourier_periodic_state"

        def progress(_):
            now = time.monotonic()
            if now >= deadline:
                raise IntegrationInterrupted(
                    "Spline campaign wall-clock budget exhausted."
                )
            if now >= candidate_deadline:
                raise IntegrationInterrupted(
                    "Candidate wall-clock budget exhausted."
                )

        retried_safe = False
        try:
            try:
                backend_box = {}
                def observe_periodic(periodic):
                    backend_box["statistics"] = periodic.backend_statistics

                result, state = _evaluate(
                    f"free_spline_260k_{this_proposal}",
                    candidate_design,
                    definition,
                    initial_state=initial,
                    progress_callback=progress,
                    periodic_observer=observe_periodic,
                )
            except MicrotubeDomainError:
                # A distant motion can make a previous periodic state physically
                # incompatible at theta=0. Rebuild only the gas charge uniformly
                # for the target geometry, preserving same-hardware wall energy.
                retried_safe = True
                safe = _safe_uniform_initial(
                    candidate_design,
                    wall_source=initial,
                )
                result, state = _evaluate(
                    f"free_spline_260k_{this_proposal}_safe",
                    candidate_design,
                    definition,
                    initial_state=safe,
                    progress_callback=progress,
                    periodic_observer=observe_periodic,
                )
                warm_mode = "safe_uniform_retry_after_domain_error"
        except MicrotubeDomainError as exc:
            result = {
                "status": "invalid_exchanger",
                "message": str(exc),
            }
            state = None
        except IntegrationInterrupted as exc:
            result = {
                "status": "interrupted",
                "message": str(exc),
            }
            state = None
        except (ValueError, RuntimeError, ArithmeticError) as exc:
            result = {
                "status": "integration_failure",
                "message": str(exc),
            }
            state = None

        feasible, reasons = (
            feasibility(result, POWER_FLOOR_W)
            if result.get("status") == "converged"
            else (False, [])
        )

        shape = _shape_diagnostics(candidate_kinematics)
        source_distance = _shape_rms(
            source_kinematics,
            candidate_kinematics,
            sample_count=2048,
        )

        record = {
            "index": len(history),
            "proposal_index": this_proposal,
            "candidate_id": cid,
            "kind": kind,
            "controls_per_piston": args.controls,
            "radius_index": radius_index,
            "shape_radius": radius,
            "phase_radius_deg": phase_radius,
            "small_controls": list(map(float, small)),
            "large_controls": list(map(float, large)),
            "small_phase_deg": float(small_phase) % 360.0,
            "large_phase_deg": float(large_phase) % 360.0,
            "shape_diagnostics": shape,
            "shape_distance_vs_source_fourier": source_distance,
            "result": result,
            "feasible": feasible,
            "physical_constraint_failures": reasons,
            "elapsed_seconds": time.monotonic() - before,
            "warm_start_source": warm_source,
            "warm_start_mode": warm_mode,
            "warm_start_distance": warm_distance,
            "backend_statistics": backend_box.get("statistics"),
            "safe_retry_used": retried_safe,
            "last_complete_state": (
                state.tolist() if state is not None else None
            ),
        }

        with history_path.open("a") as stream:
            stream.write(json.dumps(record) + "\n")
            stream.flush()

        history.append(record)
        completed += 1

        if result.get("status") != "interrupted":
            finished.add(cid)

        moved = False
        if feasible and (best is None or _eta(record) > _eta(best)):
            best = record
            moved = True
            _save_motion(
                out / "best_motion.csv",
                source_kinematics,
                candidate_kinematics,
            )

        output_report = {
            "best_feasible": best,
            "source_fourier": {
                "candidate_id": source_best.get("candidate_id"),
                "harmonics": source_h,
                "efficiency": source_best["result"][
                    "indicated_thermal_efficiency"
                ],
                "power_W": source_best["result"]["indicated_power_w"],
                "total_mass_kg": total_mass_kg,
                "thermo": source_best["thermo"],
            },
            "projected_seed": {
                "controls_per_piston": args.controls,
                "shape_diagnostics": seed_shape,
                "shape_rms_vs_source_fourier": seed_rms,
            },
            "attempted_total": len(history),
            "completed_this_run": completed,
            "requested_duration_seconds": args.budget_seconds,
            "actual_duration_seconds": time.monotonic() - start,
            "definition": identity,
        }
        report_path.write_text(
            json.dumps(output_report, indent=2) + "\n"
        )

        print(
            json.dumps(
                {
                    "index": record["index"],
                    "proposal_index": this_proposal,
                    "kind": kind,
                    "radius": radius,
                    "phase_radius_deg": phase_radius,
                    "status": result.get("status"),
                    "feasible": feasible,
                    "efficiency": result.get(
                        "indicated_thermal_efficiency"
                    ),
                    "power_W": result.get("indicated_power_w"),
                    "position_rms_vs_8h": source_distance[
                        "combined_position_rms"
                    ],
                    "center_moved": moved,
                    "best_efficiency": _eta(best) if best else None,
                    "safe_retry_used": retried_safe,
                    "warm_start_mode": warm_mode,
                    "actual_backend": (
                        (backend_box.get("statistics") or {}).get("actual_backend")
                    ),
                    "elapsed_seconds": record["elapsed_seconds"],
                }
            ),
            flush=True,
        )

    print(f"Saved {report_path}", flush=True)


if __name__ == "__main__":
    main()
