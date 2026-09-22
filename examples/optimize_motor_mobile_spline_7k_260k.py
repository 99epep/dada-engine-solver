#!/usr/bin/env python3
"""Seven-mobile-knot periodic cubic spline optimization at fixed 260 K machine.

Purpose
-------
Replace the 16 equally-spaced spline controls by only 7 mobile knots per piston.
Each piston independently optimizes:
  * 7 spline control values (canonicalized: offset/amplitude are redundant),
  * 7 positive cyclic knot gaps whose sum is one revolution,
  * one global phase.

Knot ordering cannot change. A minimum 12 degree cyclic separation is enforced,
so knots cannot collapse onto each other. Hardware, gas inventory and thermal
conditions are frozen to the final 8H machine.

The campaign is seeded by projecting the current 16-control champion onto
7 initially uniform knots. The source 16-control motion is used only as the seed
and comparison reference; the 7-knot representation is then fully free.

Typical run:

    PYTHONPATH=src python3 examples/optimize_motor_mobile_spline_7k_260k.py \
      --source-report outputs/motor_spline_from_linear_260k/report.json \
      --fourier-report outputs/motor_fourier_c2_8h_refine/report.json \
      --output-directory outputs/motor_mobile_spline_7k_260k \
      --evaluations 2048 \
      --budget-seconds 21600

Outputs:
    definition.json
    history.jsonl
    report.json
    best_motion.csv
"""

from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass, replace
from pathlib import Path
from types import SimpleNamespace
import argparse
import hashlib
import json
import math
import time

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.stats import qmc

from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
from dada_solver.integration import IntegrationInterrupted
from dada_solver.wall_backend import WallBackendSettings

from compare_motor_motion_laws_stage7A5 import ROOT, _evaluate
from optimize_motor_temperature_point import base_geometry, feasibility
from refine_motor_four_stage_hx9d_variable_gas import _load_basis
from refine_motor_fourier_c2_260k import build_design as fourier_build_design

from optimize_motor_free_spline_260k_v3 import (
    _eta,
    _eligible,
    _fixed_inventory_design,
    _load_history,
    _make_kinematics,
    _normalized_motion,
    _safe_uniform_initial,
    _shape_diagnostics,
    _shape_rms,
    _source_best,
    _source_harmonics,
    _canonical,
)

TAU = 2.0 * math.pi
POWER_FLOOR_W = 25.0
KNOTS = 7
MIN_GAP_DEG = 12.0

DEFAULT_SOURCE = ROOT / "outputs/motor_spline_from_linear_260k/report.json"
DEFAULT_FOURIER = ROOT / "outputs/motor_fourier_c2_8h_refine/report.json"
DEFAULT_OUTPUT = ROOT / "outputs/motor_mobile_spline_7k_260k"

# Broad-to-fine radii. Shape radius acts in canonical control space.
# Knot radius acts on zero-mean log-gap coordinates.
SHAPE_RADII = (0.060, 0.038, 0.024, 0.014, 0.008, 0.0040, 0.0018)
KNOT_RADII = (0.40, 0.28, 0.18, 0.11, 0.065, 0.035, 0.016)
PHASE_RADII_DEG = (2.5, 1.7, 1.1, 0.70, 0.40, 0.22, 0.10)
EVALS_PER_RADIUS = 192

WIDE_PROBE_EVERY = 32
WIDE_SHAPE_RADIUS = 0.10
WIDE_KNOT_RADIUS = 0.75
WIDE_PHASE_RADIUS_DEG = 4.0

DIRECT_WARM_DISTANCE = 0.075


def _resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def _canonical_logits(values):
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) != KNOTS or not np.all(np.isfinite(x)):
        raise ValueError("Invalid knot logit vector.")
    return x - np.mean(x)


def _gap_fractions(logits):
    """Positive cyclic gaps with hard minimum separation and exact unit sum."""
    z = _canonical_logits(logits)
    z = z - np.max(z)
    w = np.exp(z)
    w /= np.sum(w)

    minimum = MIN_GAP_DEG / 360.0
    remainder = 1.0 - KNOTS * minimum
    if remainder <= 0:
        raise RuntimeError("MIN_GAP_DEG is too large.")
    return minimum + remainder * w


def _knots_from_logits(logits):
    gaps = _gap_fractions(logits)
    fractions = np.r_[0.0, np.cumsum(gaps)]
    fractions[-1] = 1.0
    return tuple(float(TAU * x) for x in fractions)


@dataclass(frozen=True)
class _MobilePeriodicMotion:
    control_values: tuple[float, ...]
    gap_logits: tuple[float, ...]
    minimum_volume: float
    maximum_volume: float

    def __post_init__(self):
        controls = _canonical(self.control_values)
        if len(controls) != KNOTS:
            raise ValueError(f"Require exactly {KNOTS} control values.")
        logits = _canonical_logits(self.gap_logits)
        knots = _knots_from_logits(logits)

        spline = CubicSpline(
            np.asarray(knots),
            np.r_[controls, controls[0]],
            bc_type="periodic",
        )

        roots = spline.derivative().roots(extrapolate=False)
        roots = roots[np.isfinite(roots) & (roots >= 0.0) & (roots <= TAU)]
        candidates = np.r_[np.asarray(knots), roots]
        values = spline(candidates)
        minimum = float(np.min(values))
        maximum = float(np.max(values))
        if not math.isfinite(minimum + maximum) or maximum - minimum <= 1e-12:
            raise ValueError("Degenerate mobile spline.")

        scale = (self.maximum_volume - self.minimum_volume) / (maximum - minimum)

        object.__setattr__(self, "control_values", tuple(float(x) for x in controls))
        object.__setattr__(self, "gap_logits", tuple(float(x) for x in logits))
        object.__setattr__(self, "_knots", knots)
        object.__setattr__(
            self, "_coefficients",
            tuple(tuple(float(x) for x in row) for row in spline.c),
        )
        object.__setattr__(self, "_minimum", minimum)
        object.__setattr__(self, "_scale", float(scale))

    @property
    def gap_fractions(self):
        return _gap_fractions(self.gap_logits)

    @property
    def knot_degrees(self):
        return np.degrees(np.asarray(self._knots[:-1]))

    def _scalar(self, theta, order=0):
        wrapped = float(theta) % TAU
        i = bisect_right(self._knots, wrapped) - 1
        if i >= KNOTS:
            i = KNOTS - 1
        dx = wrapped - self._knots[i]
        a, b, c, d = (row[i] for row in self._coefficients)

        if order == 0:
            raw = ((a * dx + b) * dx + c) * dx + d
            return (
                (raw - self._minimum) * self._scale
                + self.minimum_volume
            )
        if order == 1:
            return (3.0 * a * dx * dx + 2.0 * b * dx + c) * self._scale
        if order == 2:
            return (6.0 * a * dx + 2.0 * b) * self._scale
        raise ValueError("Only orders 0, 1 and 2 are supported.")

    def evaluate(self, theta, order=0):
        a = np.asarray(theta, dtype=float)
        if a.ndim == 0:
            return float(self._scalar(float(a), order))
        return np.asarray([self._scalar(float(x), order) for x in a], dtype=float)


@dataclass(frozen=True)
class MobileSplineKinematics:
    small: _MobilePeriodicMotion
    large: _MobilePeriodicMotion
    small_phase_deg: float
    large_phase_deg: float

    @property
    def small_volume_limits(self):
        from dada_solver.geometry import CylinderVolumeLimits
        return CylinderVolumeLimits(
            self.small.minimum_volume, self.small.maximum_volume
        )

    @property
    def large_volume_limits(self):
        from dada_solver.geometry import CylinderVolumeLimits
        return CylinderVolumeLimits(
            self.large.minimum_volume, self.large.maximum_volume
        )

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

    @property
    def small_phase_rad(self):
        return math.radians(float(self.small_phase_deg))

    @property
    def large_phase_rad(self):
        return math.radians(float(self.large_phase_deg))

    def small_cylinder_volume(self, theta):
        return self.small.evaluate(np.asarray(theta) - self.small_phase_rad, 0)

    def large_cylinder_volume(self, theta):
        return self.large.evaluate(np.asarray(theta) - self.large_phase_rad, 0)

    def small_cylinder_volume_derivative(self, theta):
        return self.small.evaluate(np.asarray(theta) - self.small_phase_rad, 1)

    def large_cylinder_volume_derivative(self, theta):
        return self.large.evaluate(np.asarray(theta) - self.large_phase_rad, 1)

    def small_cylinder_volume_second_derivative(self, theta):
        return self.small.evaluate(np.asarray(theta) - self.small_phase_rad, 2)

    def large_cylinder_volume_second_derivative(self, theta):
        return self.large.evaluate(np.asarray(theta) - self.large_phase_rad, 2)

    def cylinder_volumes_and_derivatives(self, theta):
        # Fast scalar path used by CompiledWallRHS.
        if np.ndim(theta):
            return (
                self.small_cylinder_volume(theta),
                self.large_cylinder_volume(theta),
                self.small_cylinder_volume_derivative(theta),
                self.large_cylinder_volume_derivative(theta),
            )
        x = float(theta)
        return (
            self.small._scalar(x - self.small_phase_rad, 0),
            self.large._scalar(x - self.large_phase_rad, 0),
            self.small._scalar(x - self.small_phase_rad, 1),
            self.large._scalar(x - self.large_phase_rad, 1),
        )

    def breakpoint_angles(self):
        return ()


def _make_mobile_kinematics(
    limits,
    small_controls,
    large_controls,
    small_logits,
    large_logits,
    small_phase_deg,
    large_phase_deg,
):
    s = _MobilePeriodicMotion(
        tuple(map(float, small_controls)),
        tuple(map(float, small_logits)),
        limits.small_cylinder.minimum,
        limits.small_cylinder.maximum,
    )
    l = _MobilePeriodicMotion(
        tuple(map(float, large_controls)),
        tuple(map(float, large_logits)),
        limits.large_cylinder.minimum,
        limits.large_cylinder.maximum,
    )
    return MobileSplineKinematics(
        s, l,
        float(small_phase_deg) % 360.0,
        float(large_phase_deg) % 360.0,
    )


def _source_kinematics(source_best, limits):
    return _make_kinematics(
        limits,
        np.asarray(source_best["small_controls"], dtype=float),
        np.asarray(source_best["large_controls"], dtype=float),
        float(source_best["small_phase_deg"]),
        float(source_best["large_phase_deg"]),
    )


def _project_source_to_uniform_mobile(source_kin, limits):
    """Sample the unshifted 16-control base at seven initially uniform knots."""
    base = source_kin.base
    theta = TAU * np.arange(KNOTS, dtype=float) / KNOTS

    sl = base.small_volume_limits
    ll = base.large_volume_limits
    small = np.asarray([
        (base.small_cylinder_volume(float(a)) - sl.minimum) / sl.swept
        for a in theta
    ])
    large = np.asarray([
        (base.large_cylinder_volume(float(a)) - ll.minimum) / ll.swept
        for a in theta
    ])

    return (
        _canonical(small),
        _canonical(large),
        np.zeros(KNOTS),
        np.zeros(KNOTS),
        math.degrees(float(source_kin.small_phase_rad)),
        math.degrees(float(source_kin.large_phase_rad)),
    )


def _robust_extrema_count(derivative):
    d = np.asarray(derivative, dtype=float)
    threshold = max(1e-8, 1e-4 * float(np.max(np.abs(d))))
    signs = np.sign(d)
    signs[np.abs(d) < threshold] = 0.0
    signs = signs[signs != 0.0]
    if len(signs) == 0:
        return 0
    return int(np.sum(signs[1:] != signs[:-1]) + (signs[0] != signs[-1]))


def _admissible(kinematics):
    t = np.linspace(0.0, 1.0, 4096, endpoint=False)
    _, dq = _normalized_motion(kinematics, t)
    return (
        _robust_extrema_count(dq[0]) == 2
        and _robust_extrema_count(dq[1]) == 2
    )


def _smooth_controls(values):
    v = np.asarray(values, dtype=float)
    return _canonical(
        (np.roll(v, 1) + 2.0 * v + np.roll(v, -1)) / 4.0
    )


def _make_admissible_seed(limits, seed):
    small, large, sl, ll, ps, pl = seed
    for passes in range(7):
        kin = _make_mobile_kinematics(
            limits, small, large, sl, ll, ps, pl
        )
        if _admissible(kin):
            return small, large, sl, ll, ps, pl, kin, passes
        small = _smooth_controls(small)
        large = _smooth_controls(large)
    raise RuntimeError("Could not create an admissible 7-knot seed.")


def _direction(u, center):
    d = 2.0 * np.asarray(u, dtype=float) - 1.0
    d -= np.mean(d)
    d -= center * float(d @ center)
    n = float(np.linalg.norm(d))
    if n < 1e-12:
        d = np.roll(center, 1) - center
        d -= np.mean(d)
        d -= center * float(d @ center)
        n = float(np.linalg.norm(d))
    return d / n


def _logit_direction(u):
    d = 2.0 * np.asarray(u, dtype=float) - 1.0
    d -= np.mean(d)
    n = float(np.linalg.norm(d))
    if n < 1e-12:
        d = np.arange(len(d), dtype=float)
        d -= np.mean(d)
        n = float(np.linalg.norm(d))
    return d / n


def _propose(center, u, shape_radius, knot_radius, phase_radius):
    n = KNOTS

    s0 = np.asarray(center["small_controls"], dtype=float)
    l0 = np.asarray(center["large_controls"], dtype=float)
    sg0 = np.asarray(center["small_gap_logits"], dtype=float)
    lg0 = np.asarray(center["large_gap_logits"], dtype=float)

    ds = _direction(u[0:n], s0)
    dl = _direction(u[n:2*n], l0)
    dsg = _logit_direction(u[2*n:3*n])
    dlg = _logit_direction(u[3*n:4*n])

    rs = shape_radius * (0.20 + 0.80 * float(u[4*n]))
    rl = shape_radius * (0.20 + 0.80 * float(u[4*n+1]))
    rsg = knot_radius * (0.20 + 0.80 * float(u[4*n+2]))
    rlg = knot_radius * (0.20 + 0.80 * float(u[4*n+3]))

    small = _canonical(s0 + rs * ds)
    large = _canonical(l0 + rl * dl)
    small_logits = _canonical_logits(sg0 + rsg * dsg)
    large_logits = _canonical_logits(lg0 + rlg * dlg)

    us = (float(u[0]) + float(u[2*n])) % 1.0
    ul = (float(u[n]) + float(u[3*n])) % 1.0
    ps = (float(center["small_phase_deg"]) + phase_radius*(2*us-1)) % 360.0
    pl = (float(center["large_phase_deg"]) + phase_radius*(2*ul-1)) % 360.0

    return small, large, small_logits, large_logits, ps, pl


def _candidate_id(s, l, sg, lg, ps, pl):
    payload = {
        "small_controls": [float(x) for x in s],
        "large_controls": [float(x) for x in l],
        "small_gap_logits": [float(x) for x in sg],
        "large_gap_logits": [float(x) for x in lg],
        "small_phase_deg": float(ps) % 360.0,
        "large_phase_deg": float(pl) % 360.0,
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _phase_difference_deg(a, b):
    return ((float(a)-float(b)+180.0) % 360.0) - 180.0


def _record_distance(s, l, sg, lg, ps, pl, r):
    ds = np.linalg.norm(np.asarray(s) - np.asarray(r["small_controls"]))
    dl = np.linalg.norm(np.asarray(l) - np.asarray(r["large_controls"]))
    dsg = np.linalg.norm(
        _canonical_logits(sg) - _canonical_logits(r["small_gap_logits"])
    )
    dlg = np.linalg.norm(
        _canonical_logits(lg) - _canonical_logits(r["large_gap_logits"])
    )
    dps = _phase_difference_deg(ps, r["small_phase_deg"]) / 30.0
    dpl = _phase_difference_deg(pl, r["large_phase_deg"]) / 30.0
    return float(math.sqrt(ds*ds + dl*dl + 0.15*dsg*dsg + 0.15*dlg*dlg + dps*dps + dpl*dpl))


def _motor_knot_degrees(motion, phase_deg):
    base = np.degrees(np.asarray(motion._knots[:-1]))
    return [float(x) for x in np.mod(-(base + float(phase_deg)), 360.0)]


def _knot_diagnostics(kin):
    return {
        "minimum_gap_deg": MIN_GAP_DEG,
        "small_gap_deg": [
            float(360.0*x) for x in kin.small.gap_fractions
        ],
        "large_gap_deg": [
            float(360.0*x) for x in kin.large.gap_fractions
        ],
        "small_motor_knot_deg": _motor_knot_degrees(
            kin.small, kin.small_phase_deg
        ),
        "large_motor_knot_deg": _motor_knot_degrees(
            kin.large, kin.large_phase_deg
        ),
    }


def _save_motion(path, source, candidate):
    t = np.linspace(0.0, 1.0, 2881)
    qs, ds = _normalized_motion(source, t)
    qc, dc = _normalized_motion(candidate, t)
    np.savetxt(
        path,
        np.column_stack((
            t, 360.0*t,
            qs[0], qs[1], qc[0], qc[1],
            ds[0], ds[1], dc[0], dc[1],
        )),
        delimiter=",",
        header=(
            "motor_time_fraction,motor_angle_deg,"
            "source16_small,source16_large,mobile7_small,mobile7_large,"
            "source16_dsmall_dt,source16_dlarge_dt,"
            "mobile7_dsmall_dt,mobile7_dlarge_dt"
        ),
        comments="",
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-report", type=Path, default=DEFAULT_SOURCE)
    ap.add_argument("--fourier-report", type=Path, default=DEFAULT_FOURIER)
    ap.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--evaluations", type=int, default=2048)
    ap.add_argument("--budget-seconds", type=float, default=21600.0)
    ap.add_argument("--candidate-seconds", type=float, default=150.0)
    ap.add_argument("--evaluations-per-radius", type=int, default=EVALS_PER_RADIUS)
    ap.add_argument("--seed", type=int, default=26092317)
    args = ap.parse_args()

    source_path = _resolve(args.source_report)
    fourier_path = _resolve(args.fourier_report)
    out = _resolve(args.output_directory)
    out.mkdir(parents=True, exist_ok=True)

    sr = json.loads(source_path.read_text())
    sb = sr.get("best_feasible")
    if sb is None:
        raise RuntimeError("Source report has no best_feasible.")

    fr = json.loads(fourier_path.read_text())
    fb = _source_best(fr)
    harmonics = _source_harmonics(fb)

    campaign_definition, base = _load_basis()
    geometry = base_geometry(base)

    machine = fourier_build_design(
        base, geometry, fb["thermo"],
        np.asarray(fb["small_coefficients"], dtype=float),
        np.asarray(fb["large_coefficients"], dtype=float),
        harmonics,
    )
    total_mass = float(sb["result"]["total_mass_kg"])
    machine = _fixed_inventory_design(machine, total_mass)
    limits = machine.configuration.machine_volumes

    source_kin = _source_kinematics(sb, limits)

    seed = _project_source_to_uniform_mobile(source_kin, limits)
    (
        seed_s, seed_l, seed_sg, seed_lg,
        seed_ps, seed_pl, seed_kin, smoothing_passes,
    ) = _make_admissible_seed(limits, seed)

    seed_distance = _shape_rms(source_kin, seed_kin, sample_count=4096)

    definition = SimpleNamespace(
        wall_numerical_settings=campaign_definition.wall_numerical_settings,
        wall_backend=WallBackendSettings("numba"),
        exact_kinematics_cache=True,
        shared_replay=True,
    )

    identity = {
        "study": "260 K seven-mobile-knot spline refinement v1",
        "source_report": str(source_path),
        "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "fourier_report": str(fourier_path),
        "fourier_sha256": hashlib.sha256(fourier_path.read_bytes()).hexdigest(),
        "source_candidate_id": sb["candidate_id"],
        "knots_per_piston": KNOTS,
        "mobile_knot_positions": True,
        "minimum_gap_deg": MIN_GAP_DEG,
        "projection_smoothing_passes": smoothing_passes,
        "projection_distance_vs_source16": seed_distance,
        "fixed_total_mass_kg": total_mass,
        "fixed_thermo": fb["thermo"],
        "shape_radii": list(SHAPE_RADII),
        "knot_radii": list(KNOT_RADII),
        "phase_radii_deg": list(PHASE_RADII_DEG),
        "evaluations_per_radius": args.evaluations_per_radius,
        "wide_probe_every": WIDE_PROBE_EVERY,
        "wide_shape_radius": WIDE_SHAPE_RADIUS,
        "wide_knot_radius": WIDE_KNOT_RADIUS,
        "wide_phase_radius_deg": WIDE_PHASE_RADIUS_DEG,
        "sobol_seed": args.seed,
        "constraints": {
            "exactly_one_minimum_and_one_maximum_per_piston": True,
            "C2_periodic_position": True,
            "knot_order_preserved": True,
            "minimum_knot_separation_deg": MIN_GAP_DEG,
            "hardware_frozen": True,
            "gas_inventory_frozen": True,
            "velocity_penalty": False,
            "acceleration_penalty": False,
        },
    }

    dp = out / "definition.json"
    hp = out / "history.jsonl"
    rp = out / "report.json"

    if dp.exists() and json.loads(dp.read_text()) != identity:
        raise RuntimeError("Campaign definition changed; use a new output directory.")
    if not dp.exists():
        dp.write_text(json.dumps(identity, indent=2) + "\n")

    history = _load_history(hp)
    valid = [r for r in history if _eligible(r)]
    best = max(valid, key=_eta) if valid else None
    finished = {
        r["candidate_id"]
        for r in history
        if r.get("result", {}).get("status") != "interrupted"
    }

    dimension = 4*KNOTS + 4
    points = qmc.Sobol(
        d=dimension, scramble=True, seed=args.seed
    ).random_base2(16)

    seed_record = {
        "small_controls": seed_s.tolist(),
        "large_controls": seed_l.tolist(),
        "small_gap_logits": seed_sg.tolist(),
        "large_gap_logits": seed_lg.tolist(),
        "small_phase_deg": seed_ps,
        "large_phase_deg": seed_pl,
    }

    start = time.monotonic()
    deadline = start + args.budget_seconds
    completed = 0
    proposal_index = max(
        (r.get("proposal_index", -1) for r in history), default=-1
    ) + 1

    while completed < args.evaluations and time.monotonic() < deadline:
        idx = proposal_index
        proposal_index += 1

        if not history and idx == 0:
            s, l = seed_s.copy(), seed_l.copy()
            sg, lg = seed_sg.copy(), seed_lg.copy()
            ps, pl = seed_ps, seed_pl
            kind = "uniform_7k_projection_seed"
            radius_index = -1
            sradius = kradius = pradius = 0.0
        else:
            searched = sum(
                1 for r in history
                if r.get("kind") != "uniform_7k_projection_seed"
            )
            radius_index = min(
                searched // args.evaluations_per_radius,
                len(SHAPE_RADII)-1,
            )

            if (
                searched > 0
                and searched % WIDE_PROBE_EVERY == WIDE_PROBE_EVERY - 1
            ):
                center = seed_record
                sradius = WIDE_SHAPE_RADIUS
                kradius = WIDE_KNOT_RADIUS
                pradius = WIDE_PHASE_RADIUS_DEG
                kind = "mobile7_wide_probe"
            else:
                center = best if best is not None else seed_record
                sradius = SHAPE_RADII[radius_index]
                kradius = KNOT_RADII[radius_index]
                pradius = PHASE_RADII_DEG[radius_index]
                kind = "mobile7_refinement"

            s, l, sg, lg, ps, pl = _propose(
                center,
                points[idx % len(points)],
                sradius, kradius, pradius,
            )

        cid = _candidate_id(s, l, sg, lg, ps, pl)
        if cid in finished:
            continue

        kin = _make_mobile_kinematics(
            limits, s, l, sg, lg, ps, pl
        )
        if not _admissible(kin):
            continue

        design = replace(machine, kinematics=kin)
        before = time.monotonic()
        candidate_deadline = min(
            deadline, before + args.candidate_seconds
        )

        reusable = [
            r for r in history
            if r.get("result", {}).get("status") == "converged"
            and r.get("last_complete_state") is not None
        ]

        if reusable:
            warm = min(
                reusable,
                key=lambda r: _record_distance(s, l, sg, lg, ps, pl, r),
            )
            warm_distance = _record_distance(
                s, l, sg, lg, ps, pl, warm
            )
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
                warm_mode = "direct_safe_uniform_for_distant_shape"
        else:
            initial = np.asarray(sb["last_complete_state"], dtype=float)
            warm_distance = None
            warm_source = sb["candidate_id"]
            warm_mode = "source16_periodic_state"

        def progress(_):
            now = time.monotonic()
            if now >= deadline:
                raise IntegrationInterrupted("Mobile-7 campaign budget exhausted.")
            if now >= candidate_deadline:
                raise IntegrationInterrupted("Candidate budget exhausted.")

        backend_box = {}
        def observe(periodic):
            backend_box["statistics"] = periodic.backend_statistics

        safe_retry = False
        try:
            try:
                result, state = _evaluate(
                    f"mobile7_260k_{idx}",
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
                result, state = _evaluate(
                    f"mobile7_260k_{idx}_safe",
                    design,
                    definition,
                    initial_state=safe,
                    progress_callback=progress,
                    periodic_observer=observe,
                )
                warm_mode = "safe_uniform_retry_after_domain_error"
        except MicrotubeDomainError as exc:
            result = {"status": "invalid_exchanger", "message": str(exc)}
            state = None
        except IntegrationInterrupted as exc:
            result = {"status": "interrupted", "message": str(exc)}
            state = None
        except (ValueError, RuntimeError, ArithmeticError) as exc:
            result = {"status": "integration_failure", "message": str(exc)}
            state = None

        feasible, reasons = (
            feasibility(result, POWER_FLOOR_W)
            if result.get("status") == "converged"
            else (False, [])
        )

        distance = _shape_rms(source_kin, kin, sample_count=2048)
        rec = {
            "index": len(history),
            "proposal_index": idx,
            "candidate_id": cid,
            "kind": kind,
            "knots_per_piston": KNOTS,
            "radius_index": radius_index,
            "shape_radius": sradius,
            "knot_radius": kradius,
            "phase_radius_deg": pradius,
            "small_controls": s.tolist(),
            "large_controls": l.tolist(),
            "small_gap_logits": sg.tolist(),
            "large_gap_logits": lg.tolist(),
            "small_phase_deg": float(ps) % 360.0,
            "large_phase_deg": float(pl) % 360.0,
            "knot_diagnostics": _knot_diagnostics(kin),
            "shape_diagnostics": _shape_diagnostics(kin),
            "shape_distance_vs_source16": distance,
            "result": result,
            "feasible": feasible,
            "physical_constraint_failures": reasons,
            "elapsed_seconds": time.monotonic() - before,
            "warm_start_source": warm_source,
            "warm_start_mode": warm_mode,
            "warm_start_distance": warm_distance,
            "backend_statistics": backend_box.get("statistics"),
            "safe_retry_used": safe_retry,
            "last_complete_state": (
                state.tolist() if state is not None else None
            ),
        }

        with hp.open("a") as f:
            f.write(json.dumps(rec) + "\n")
            f.flush()

        history.append(rec)
        completed += 1
        if result.get("status") != "interrupted":
            finished.add(cid)

        moved = False
        if feasible and (best is None or _eta(rec) > _eta(best)):
            best = rec
            moved = True
            _save_motion(out/"best_motion.csv", source_kin, kin)

        report = {
            "best_feasible": best,
            "source16": {
                "candidate_id": sb["candidate_id"],
                "efficiency": sb["result"]["indicated_thermal_efficiency"],
                "power_W": sb["result"]["indicated_power_w"],
                "shape_diagnostics": sb["shape_diagnostics"],
            },
            "projected_seed": {
                "knots_per_piston": KNOTS,
                "smoothing_passes": smoothing_passes,
                "shape_distance_vs_source16": seed_distance,
                "shape_diagnostics": _shape_diagnostics(seed_kin),
                "knot_diagnostics": _knot_diagnostics(seed_kin),
            },
            "attempted_total": len(history),
            "completed_this_run": completed,
            "requested_duration_seconds": args.budget_seconds,
            "actual_duration_seconds": time.monotonic() - start,
            "definition": identity,
        }
        rp.write_text(json.dumps(report, indent=2) + "\n")

        kd = rec["knot_diagnostics"]
        print(json.dumps({
            "index": rec["index"],
            "proposal_index": idx,
            "kind": kind,
            "shape_radius": sradius,
            "knot_radius": kradius,
            "phase_radius_deg": pradius,
            "status": result.get("status"),
            "feasible": feasible,
            "efficiency": result.get("indicated_thermal_efficiency"),
            "power_W": result.get("indicated_power_w"),
            "small_motor_knots_deg": [
                round(x, 2) for x in kd["small_motor_knot_deg"]
            ],
            "large_motor_knots_deg": [
                round(x, 2) for x in kd["large_motor_knot_deg"]
            ],
            "rms_vs_source16": distance["combined_position_rms"],
            "center_moved": moved,
            "best_efficiency": _eta(best) if best else None,
            "actual_backend": (
                (backend_box.get("statistics") or {}).get("actual_backend")
            ),
            "elapsed_seconds": rec["elapsed_seconds"],
        }), flush=True)

    print(f"Saved {rp}", flush=True)


if __name__ == "__main__":
    main()
