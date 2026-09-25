#!/usr/bin/env python3
"""Coupled thermodynamic optimization of paired six-bar mechanisms around candidate 3952.

Fixed in this campaign:
- the exact 260 K thermodynamic/hardware basis used by the structured-15 search
  that produced candidate 3952;
- exchanger geometry;
- gas inventory;
- source temperatures;
- charge;
- solver settings.

Free variables:
- 15 continuous dimensions of the SMALL six-bar;
- 15 continuous dimensions of the LARGE six-bar.

The two mechanisms are perturbed SIMULTANEOUSLY at every proposal. Their
primary/secondary assembly branches remain fixed to the selected family. The
initial S/L pair is the mirrored-and-polished pair in outputs/sixbar_pairs_3952/.

Candidate-3952 motion is NOT an optimization objective. Position RMS is computed
only as a diagnostic. The objective is actual indicated thermal efficiency of
the fixed 260 K machine.

Search strategy:
- one persistent local campaign per family (default 1,4,12,50);
- incumbent-centred 30-D scrambled-Sobol proposals;
- S and L both change in every non-seed proposal;
- cheap dense mechanical screening before thermodynamic integration;
- progressively shrinking trust radii;
- the exact current pair is evaluated first;
- nearest converged state from the same family is reused as warm start;
- families are interleaved round-robin.

Typical run:
    PYTHONPATH=src python3 examples/optimize_sixbar_pairs_thermo_3952.py \
      --evaluations-per-family 64 \
      --budget-seconds 21600

Rerunning resumes append-only histories.
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
from scipy.stats import qmc

from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
from dada_solver.integration import IntegrationInterrupted
from dada_solver.six_bar import (
    IndependentSixBarVolumeKinematics,
    SixBarCylinderMechanism,
)
from dada_solver.wall_backend import WallBackendSettings

import optimize_motor_hybrid_c2_15p_260k as hybrid15
import search_large_sixbar_stageL1 as stage
import search_large_second_dyad_v41_panel_3952 as v41
import search_motor_hybrid_3952_sixbar as motor3952


ROOT = Path.cwd()
PAIR_DIR = ROOT / "outputs" / "sixbar_pairs_3952"
OUTPUT = ROOT / "outputs" / "sixbar_thermo_coupled_3952_fine"
TARGET = ROOT / "outputs" / "motor_hybrid_c2_15p_260k" / "candidate_3952_motion_target.csv"
HISTORY_3952 = ROOT / "outputs" / "motor_hybrid_c2_15p_260k" / "history.jsonl"

NAMES = tuple(stage.NAMES)
PAIR_NAMES = tuple(f"S_{n}" for n in NAMES) + tuple(f"L_{n}" for n in NAMES)

POWER_FLOOR_W = 25.0
DEFAULT_RADIUS_MULTIPLIERS = (1.0, 0.75, 0.55, 0.40, 0.28, 0.20, 0.14)

# Fine trust region: at this stage the mechanism families are already good.
# We only allow small coordinated dimensional adjustments around the incumbent.
BASE_PRIMARY_LENGTH_FRAC = 0.015       # +/- 1.5 %
BASE_PRIMARY_POINT_FRAC = 0.025        # +/- 2.5 % of coordinate magnitude
BASE_PRIMARY_POINT_MIN = 0.05          # or +/- 0.05 crank radius minimum
BASE_PRIMARY_PHASE_DEG = 2.0           # +/- 2 deg
BASE_PIVOT_SPAN = 0.20                 # +/- 0.20 crank radius
BASE_SECONDARY_LENGTH_FRAC = 0.03      # +/- 3 %
BASE_H_RATIO_SPAN = 0.05               # +/- 0.05 EF-ratio
BASE_ROD_FRAC = 0.04                   # +/- 4 %
BASE_AXIS_OFFSET_SPAN = 0.20           # +/- 0.20 crank radius
BASE_AXIS_ANGLE_DEG = 2.0              # +/- 2 deg

# A Sobol direction can point outside the thin feasible mechanical domain,
# especially for rank 1 which starts almost on the H-lateral limits. Keep the
# direction (so S and L still move together in all 30 dimensions), but shorten
# the step until it is feasible.
DIRECTION_BACKOFFS = (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125)

# Mechanical feasibility envelope retained from synthesis.
STROKE_FLOOR = 1.0
STROKE_CEILING = 3.0
PRIMARY_SINE_FLOOR = 0.30
SECONDARY_SINE_FLOOR = 0.30
ROD_COS_FLOOR = 0.95
MAXIMUM_EH = 7.0
CRANK_CLEARANCE_FLOOR = 0.50
H_LATERAL_RMS_MAX = 0.25
H_LATERAL_SPAN_MAX = 0.65


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def append_jsonl(path: Path, record: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(record, separators=(",", ":")) + "\n")
        stream.flush()


def load_history(path: Path):
    if not path.exists():
        return []
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def wrap_pi(a: float) -> float:
    return (float(a) + math.pi) % (2.0 * math.pi) - math.pi


def candidate_from_loader_file(path: Path) -> dict:
    data = load_json(path)
    best = data.get("best")
    if not best or best.get("feasible") is not True:
        raise RuntimeError(f"No feasible best candidate in {path}")
    return best


def vector_from_candidate(best: dict) -> np.ndarray:
    p = best["parameters"]
    return np.array([float(p[name]) for name in NAMES], dtype=float)


def vector_to_parameters(v: np.ndarray) -> dict:
    return dict(zip(NAMES, map(float, v)))


def pair_vector(small: dict, large: dict) -> np.ndarray:
    return np.concatenate((vector_from_candidate(small), vector_from_candidate(large)))


def mechanism_from_vector(v: np.ndarray, primary_branch: int, second_branch: int):
    x = np.asarray(v, dtype=float)
    return SixBarCylinderMechanism(
        primary_ground=float(x[0]),
        primary_coupler=float(x[1]),
        primary_rocker=float(x[2]),
        primary_e_along=float(x[3]),
        primary_e_normal=float(x[4]),
        primary_phase=wrap_pi(x[5]),
        second_pivot_x=float(x[6]),
        second_pivot_y=float(x[7]),
        link_ef=float(x[8]),
        link_gf=float(x[9]),
        h_along_over_ef=float(x[10]),
        h_normal_over_ef=float(x[11]),
        piston_rod=float(x[12]),
        slider_axis_offset=float(x[13]),
        slider_axis_angle=wrap_pi(x[14]),
        primary_branch=int(primary_branch),
        second_branch=int(second_branch),
    )


def mechanism_candidate_dict(v, primary_branch, second_branch, metrics):
    out = dict(metrics)
    out["parameters"] = vector_to_parameters(v)
    out["primary"] = {
        "ground": float(v[0]),
        "coupler": float(v[1]),
        "rocker": float(v[2]),
        "assembly_branch": int(primary_branch),
        "E_along": float(v[3]),
        "E_normal": float(v[4]),
        "phase": wrap_pi(float(v[5])),
    }
    out["second_branch"] = int(second_branch)
    out["feasible"] = True
    return out


def load_3952_record():
    record, field = motor3952._load_candidate(HISTORY_3952, 3952)
    return record, field


def make_fixed_3952_machine():
    """Rebuild the same fixed hardware/thermo basis used by candidate 3952."""
    sp = hybrid15.resolve(hybrid15.DEFAULT_SOURCE)
    fp = hybrid15.resolve(hybrid15.DEFAULT_FOURIER)

    source_report = load_json(sp)
    source_best = source_report.get("best_feasible")
    if source_best is None:
        raise RuntimeError("Structured-15 source report has no best_feasible.")

    fourier_report = load_json(fp)
    fourier_best = hybrid15._source_best(fourier_report)
    harmonics = hybrid15._source_harmonics(fourier_best)

    campaign_definition, base = hybrid15._load_basis()
    geometry = hybrid15.base_geometry(base)

    machine = hybrid15.fourier_build_design(
        base,
        geometry,
        fourier_best["thermo"],
        np.asarray(fourier_best["small_coefficients"]),
        np.asarray(fourier_best["large_coefficients"]),
        harmonics,
    )

    total_mass = float(source_best["result"]["total_mass_kg"])
    machine = hybrid15._fixed_inventory_design(machine, total_mass)

    definition = SimpleNamespace(
        wall_numerical_settings=campaign_definition.wall_numerical_settings,
        wall_backend=WallBackendSettings("numba"),
        exact_kinematics_cache=True,
        shared_replay=True,
    )

    identity = {
        "source_report": str(sp),
        "source_sha256": hashlib.sha256(sp.read_bytes()).hexdigest(),
        "fourier_report": str(fp),
        "fourier_sha256": hashlib.sha256(fp.read_bytes()).hexdigest(),
        "fixed_total_mass_kg": total_mass,
        "fixed_thermo": fourier_best["thermo"],
    }
    return machine, definition, source_best, identity


def load_motion_target():
    raw = np.genfromtxt(TARGET, delimiter=",", names=True)
    if raw["theta_rad"][-1] >= 2.0 * math.pi - 1e-9:
        raw = raw[:-1]
    theta = np.asarray(raw["theta_rad"], float)
    bands = motor3952._detect_noise_bands(TARGET, 6)
    masks = {
        side: ~motor3952._theta_band_mask(theta, bands[side])
        for side in ("small", "large")
    }
    return {
        "theta": theta,
        "small_q": np.asarray(raw["small_fraction_0_1"], float),
        "small_dq": np.asarray(raw["small_dq_dtheta_per_rad"], float),
        "large_q": np.asarray(raw["large_fraction_0_1"], float),
        "large_dq": np.asarray(raw["large_dq_dtheta_per_rad"], float),
        "masks": masks,
        "bands": bands,
    }


MECH_ARGS = SimpleNamespace(
    stroke_floor=STROKE_FLOOR,
    stroke_ceiling=STROKE_CEILING,
    secondary_sine_floor=SECONDARY_SINE_FLOOR,
    rod_cos_floor=ROD_COS_FLOOR,
    maximum_EH=MAXIMUM_EH,
    crank_clearance_floor=CRANK_CLEARANCE_FLOOR,
    h_lateral_rms_max=H_LATERAL_RMS_MAX,
    h_lateral_span_max=H_LATERAL_SPAN_MAX,
)


def primary_closure_ok(v: np.ndarray) -> bool:
    g, c, r = map(float, v[:3])
    return bool(
        abs(g - 1.0) > abs(c - r) + 1e-10
        and c + r > g + 1.0 + 1e-10
    )


def side_mechanical_metrics(v, *, primary_branch, second_branch, side, motion):
    """Use the tested V4.1 evaluator only as geometry/diagnostic machinery."""
    if not primary_closure_ok(v):
        raise ValueError("primary closure")

    pv = np.asarray(v[:6], float)
    downstream = np.asarray(v[6:], float)
    theta = motion["theta"]
    E, Ed, psine, pdata = stage.primary(theta, pv, primary_branch)

    result = v41.evaluate_v41(
        downstream,
        E=E,
        Ed=Ed,
        psine=psine,
        pdata=pdata,
        pv=pv,
        second_branch=second_branch,
        tq=motion[f"{side}_q"],
        tdq=motion[f"{side}_dq"],
        mask=motion["masks"][side],
        velocity_weight=0.0005,
    )

    violation = v41.violation_v41(result, MECH_ARGS)
    violation += max(
        0.0,
        PRIMARY_SINE_FLOOR - result["minimum_primary_transmission_sine"],
    ) / PRIMARY_SINE_FLOOR
    if violation > 0.0:
        raise ValueError(f"mechanical constraints: {violation}")
    return result


def screen_pair(x, *, branches, motion):
    x = np.asarray(x, float)
    sv = x[:15].copy()
    lv = x[15:].copy()
    for v in (sv, lv):
        v[5] = wrap_pi(v[5])
        v[14] = wrap_pi(v[14])

    sm = side_mechanical_metrics(
        sv,
        primary_branch=branches["small_primary"],
        second_branch=branches["small_second"],
        side="small",
        motion=motion,
    )
    lm = side_mechanical_metrics(
        lv,
        primary_branch=branches["large_primary"],
        second_branch=branches["large_second"],
        side="large",
        motion=motion,
    )

    # Production kinematics constructor performs its own dense exact geometry check.
    small_mech = mechanism_from_vector(
        sv, branches["small_primary"], branches["small_second"]
    )
    large_mech = mechanism_from_vector(
        lv, branches["large_primary"], branches["large_second"]
    )
    return sv, lv, sm, lm, small_mech, large_mech


def perturb_mechanism(center, u, radius):
    c = np.asarray(center, float)
    z = 2.0 * np.asarray(u, float) - 1.0
    y = c.copy()

    for i in (0, 1, 2):
        y[i] = max(1e-4, c[i] * (1.0 + radius * BASE_PRIMARY_LENGTH_FRAC * z[i]))

    for i in (3, 4):
        span = max(abs(c[i]) * BASE_PRIMARY_POINT_FRAC, BASE_PRIMARY_POINT_MIN)
        y[i] = c[i] + radius * span * z[i]

    y[5] = wrap_pi(c[5] + math.radians(BASE_PRIMARY_PHASE_DEG) * radius * z[5])
    y[6] = c[6] + BASE_PIVOT_SPAN * radius * z[6]
    y[7] = c[7] + BASE_PIVOT_SPAN * radius * z[7]

    for i in (8, 9):
        y[i] = max(1e-4, c[i] * (1.0 + radius * BASE_SECONDARY_LENGTH_FRAC * z[i]))

    y[10] = c[10] + BASE_H_RATIO_SPAN * radius * z[10]
    y[11] = c[11] + BASE_H_RATIO_SPAN * radius * z[11]
    y[12] = max(1e-4, c[12] * (1.0 + radius * BASE_ROD_FRAC * z[12]))
    y[13] = c[13] + BASE_AXIS_OFFSET_SPAN * radius * z[13]
    y[14] = wrap_pi(c[14] + math.radians(BASE_AXIS_ANGLE_DEG) * radius * z[14])
    return y


def propose_pair(center, u, radius):
    """One 30-D proposal: S and L are both changed simultaneously."""
    s = perturb_mechanism(center[:15], u[:15], radius)
    l = perturb_mechanism(center[15:], u[15:], radius)
    return np.concatenate((s, l))


def candidate_id(rank, branches, x):
    payload = {
        "family_rank": int(rank),
        "branches": branches,
        "parameters": {name: float(value) for name, value in zip(PAIR_NAMES, x)},
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def compact_result(result: dict):
    return hybrid15._compact_result(result)


def eligible(record: dict) -> bool:
    return bool(
        record.get("feasible")
        and record.get("result", {}).get("status") == "converged"
        and record.get("result", {}).get("indicated_thermal_efficiency") is not None
    )


def eta(record: dict) -> float:
    return float(record["result"]["indicated_thermal_efficiency"])


def normalized_pair_distance(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    terms = []
    for offset in (0, 15):
        aa = a[offset:offset+15]
        bb = b[offset:offset+15]
        for i in (0, 1, 2, 3, 4, 8, 9, 12):
            terms.append((aa[i] - bb[i]) / max(abs(bb[i]), 0.5))
        for i in (6, 7, 10, 11, 13):
            terms.append((aa[i] - bb[i]) / 2.0)
        for i in (5, 14):
            terms.append(wrap_pi(aa[i] - bb[i]) / math.radians(30.0))
    return float(np.linalg.norm(terms))


def build_design(machine, small_mech, large_mech):
    limits = machine.configuration.machine_volumes
    kin = IndependentSixBarVolumeKinematics(
        small_mech, large_mech, limits.small_cylinder, limits.large_cylinder
    )
    return replace(machine, kinematics=kin)


def evaluate_thermo(*, label, design, definition, initial, deadline, candidate_seconds):
    started = time.monotonic()
    candidate_deadline = min(deadline, started + candidate_seconds)

    def progress(_):
        now = time.monotonic()
        if now >= deadline:
            raise IntegrationInterrupted("Campaign wall-clock budget exhausted.")
        if now >= candidate_deadline:
            raise IntegrationInterrupted("Candidate wall-clock budget exhausted.")

    box = {}

    def observe(periodic):
        box["statistics"] = periodic.backend_statistics

    safe_retry = False
    warm_mode = "direct_periodic_state"
    try:
        try:
            full, state = hybrid15._evaluate(
                label,
                design,
                definition,
                initial_state=np.asarray(initial, float),
                progress_callback=progress,
                periodic_observer=observe,
            )
        except MicrotubeDomainError:
            safe_retry = True
            safe = hybrid15._safe_uniform_initial(design, wall_source=np.asarray(initial, float))
            full, state = hybrid15._evaluate(
                label + "_safe",
                design,
                definition,
                initial_state=safe,
                progress_callback=progress,
                periodic_observer=observe,
            )
            warm_mode = "safe_uniform_retry_after_domain_error"
    except MicrotubeDomainError as exc:
        full = {"status": "invalid_exchanger", "message": str(exc)}
        state = None
    except IntegrationInterrupted as exc:
        full = {"status": "interrupted", "message": str(exc)}
        state = None
    except (ValueError, RuntimeError, ArithmeticError) as exc:
        full = {"status": "integration_failure", "message": str(exc)}
        state = None

    feasible, reasons = (
        hybrid15.feasibility(full, POWER_FLOOR_W)
        if full.get("status") == "converged"
        else (False, [])
    )
    return {
        "full": full,
        "compact": compact_result(full),
        "state": state,
        "feasible": bool(feasible),
        "reasons": reasons,
        "safe_retry": safe_retry,
        "warm_mode": warm_mode,
        "backend": hybrid15._backend_compact(box.get("statistics")),
        "elapsed_seconds": time.monotonic() - started,
    }


def family_files(rank, output):
    d = output / f"rank_{rank:02d}"
    return {
        "dir": d,
        "history": d / "history.jsonl",
        "report": d / "report.json",
        "best_pair": d / "best_pair.json",
    }


def source_pair(rank, pair_dir):
    sp = pair_dir / f"rank_{rank:02d}_small.json"
    lp = pair_dir / f"rank_{rank:02d}_large.json"
    small = candidate_from_loader_file(sp)
    large = candidate_from_loader_file(lp)
    branches = {
        "small_primary": int(small["primary"]["assembly_branch"]),
        "small_second": int(small["second_branch"]),
        "large_primary": int(large["primary"]["assembly_branch"]),
        "large_second": int(large["second_branch"]),
    }
    return {
        "small_path": sp,
        "large_path": lp,
        "small": small,
        "large": large,
        "branches": branches,
        "vector": pair_vector(small, large),
    }


def save_best_pair(path, record):
    payload = {
        "family_rank": record["family_rank"],
        "candidate_id": record["candidate_id"],
        "efficiency": record["result"].get("indicated_thermal_efficiency"),
        "power_W": record["result"].get("indicated_power_w"),
        "small": record["mechanical"]["small"],
        "large": record["mechanical"]["large"],
        "thermo_result": record["result"],
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def write_family_report(path, *, rank, history, source, identity):
    good = [r for r in history if eligible(r)]
    best = max(good, key=eta) if good else None
    report = {
        "family_rank": rank,
        "source_pair": {
            "small": str(source["small_path"]),
            "large": str(source["large_path"]),
            "branches": source["branches"],
        },
        "fixed_260K_basis": identity,
        "attempted_total": len(history),
        "feasible_total": len(good),
        "best_feasible": best,
    }
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return best


def parse_radii(text):
    values = tuple(float(x) for x in text.split(",") if x.strip())
    if not values or any(x <= 0.0 for x in values):
        raise argparse.ArgumentTypeError("All radius multipliers must be positive.")
    return values


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pair-directory", type=Path, default=PAIR_DIR)
    ap.add_argument("--output-directory", type=Path, default=OUTPUT)
    ap.add_argument("--families", default="1,4,12,50")
    ap.add_argument("--evaluations-per-family", type=int, default=64)
    ap.add_argument("--evaluations-per-radius", type=int, default=16)
    ap.add_argument("--radius-multipliers", type=parse_radii, default=DEFAULT_RADIUS_MULTIPLIERS)
    ap.add_argument("--budget-seconds", type=float, default=21600.0)
    ap.add_argument("--candidate-seconds", type=float, default=150.0)
    ap.add_argument("--seed", type=int, default=3959000)
    ap.add_argument("--maximum-structural-attempts-per-evaluation", type=int, default=500)
    args = ap.parse_args()

    if min(
        args.evaluations_per_family,
        args.evaluations_per_radius,
        args.budget_seconds,
        args.candidate_seconds,
        args.maximum_structural_attempts_per_evaluation,
    ) <= 0:
        raise ValueError("Budgets/counts must be positive.")

    ranks = [int(x) for x in args.families.split(",") if x.strip()]
    machine, definition, source_best, fixed_identity = make_fixed_3952_machine()
    record3952, selected_field = load_3952_record()
    motion = load_motion_target()

    fixed_identity = {
        **fixed_identity,
        "candidate_3952_selected_by": selected_field,
        "candidate_3952_index": record3952.get("index"),
        "candidate_3952_proposal_index": record3952.get("proposal_index"),
        "candidate_3952_efficiency": record3952.get("result", {}).get("indicated_thermal_efficiency"),
        "candidate_3952_power_W": record3952.get("result", {}).get("indicated_power_w"),
        "power_floor_W": POWER_FLOOR_W,
        "motion_target": str(TARGET),
        "objective": "maximize indicated thermal efficiency; 3952 motion RMS is diagnostic only",
        "free_dimensions": 30,
        "discrete_branches": "fixed per source family",
    }

    output = args.output_directory
    output.mkdir(parents=True, exist_ok=True)
    definition_path = output / "definition.json"
    definition_payload = {
        "study": "candidate-3952 fixed-thermo simultaneous 30-D six-bar pair optimization",
        "families": ranks,
        "fixed_260K_basis": fixed_identity,
        "mechanical_constraints": {
            "stroke_over_crank": [STROKE_FLOOR, STROKE_CEILING],
            "primary_sine_minimum": PRIMARY_SINE_FLOOR,
            "secondary_sine_minimum": SECONDARY_SINE_FLOOR,
            "rod_axis_cosine_minimum": ROD_COS_FLOOR,
            "EH_over_crank_maximum": MAXIMUM_EH,
            "crank_axis_clearance_minimum": CRANK_CLEARANCE_FLOOR,
            "H_axis_lateral_rms_over_stroke_maximum": H_LATERAL_RMS_MAX,
            "H_axis_lateral_span_over_stroke_maximum": H_LATERAL_SPAN_MAX,
        },
        "search": {
            "strategy": "round-robin fine incumbent-centred 30-D scrambled Sobol; S and L always perturbed simultaneously; infeasible directions use scalar step backoff",
            "evaluations_per_radius": args.evaluations_per_radius,
            "radius_multipliers": list(args.radius_multipliers),
            "base_primary_length_fraction": BASE_PRIMARY_LENGTH_FRAC,
            "base_primary_point_fraction": BASE_PRIMARY_POINT_FRAC,
            "base_primary_point_min_span": BASE_PRIMARY_POINT_MIN,
            "base_primary_phase_deg": BASE_PRIMARY_PHASE_DEG,
            "base_pivot_span": BASE_PIVOT_SPAN,
            "base_secondary_length_fraction": BASE_SECONDARY_LENGTH_FRAC,
            "base_H_ratio_span": BASE_H_RATIO_SPAN,
            "base_rod_fraction": BASE_ROD_FRAC,
            "base_axis_offset_span": BASE_AXIS_OFFSET_SPAN,
            "base_axis_angle_deg": BASE_AXIS_ANGLE_DEG,
            "sobol_seed": args.seed,
        },
    }

    if definition_path.exists():
        if load_json(definition_path) != definition_payload:
            raise ValueError("Search definition changed; use a new output directory.")
    else:
        definition_path.write_text(json.dumps(definition_payload, indent=2) + "\n", encoding="utf-8")

    families = {}
    for family_index, rank in enumerate(ranks):
        source = source_pair(rank, args.pair_directory)
        files = family_files(rank, output)
        files["dir"].mkdir(parents=True, exist_ok=True)
        history = load_history(files["history"])
        good = [r for r in history if eligible(r)]
        best = max(good, key=eta) if good else None
        points = qmc.Sobol(
            d=30,
            scramble=True,
            seed=args.seed + 10000 * family_index,
        ).random_base2(15)
        families[rank] = {
            "source": source,
            "files": files,
            "history": history,
            "best": best,
            "points": points,
            "skipped_structural": 0,
        }

    start = time.monotonic()
    deadline = start + args.budget_seconds
    state3952 = record3952.get("last_complete_state")
    if state3952 is not None:
        state3952 = np.asarray(state3952, float)
    else:
        state3952 = np.asarray(source_best["last_complete_state"], float)

    completed_this_run = {rank: 0 for rank in ranks}

    while time.monotonic() < deadline:
        made_progress = False
        for rank in ranks:
            if time.monotonic() >= deadline:
                break

            fam = families[rank]
            history = fam["history"]
            source = fam["source"]
            completed_nonseed = sum(1 for r in history if r.get("kind") != "source_pair_seed")
            if completed_nonseed >= args.evaluations_per_family:
                continue

            proposal_index = max((r.get("proposal_index", -1) for r in history), default=-1) + 1

            if not history:
                x = source["vector"].copy()
                kind = "source_pair_seed"
                radius_index = -1
                radius = 0.0
            else:
                radius_index = min(
                    completed_nonseed // args.evaluations_per_radius,
                    len(args.radius_multipliers) - 1,
                )
                radius = float(args.radius_multipliers[radius_index])
                kind = "simultaneous_30d_local_sobol"
                center = (
                    np.array([fam["best"]["parameters"][name] for name in PAIR_NAMES], dtype=float)
                    if fam["best"] is not None
                    else source["vector"]
                )

                found = None
                attempts_used = 0
                effective_radius = None

                for attempt in range(args.maximum_structural_attempts_per_evaluation):
                    point_index = (
                        proposal_index + fam["skipped_structural"] + attempt
                    ) % len(fam["points"])
                    direction = fam["points"][point_index]

                    # Keep one 30-D direction and shorten only its amplitude.
                    # This preserves simultaneous S/L motion instead of freezing
                    # troublesome coordinates one by one.
                    for backoff in DIRECTION_BACKOFFS:
                        trial_radius = radius * backoff
                        candidate = propose_pair(center, direction, trial_radius)
                        try:
                            screen_pair(
                                candidate,
                                branches=source["branches"],
                                motion=motion,
                            )
                        except (
                            ValueError,
                            FloatingPointError,
                            OverflowError,
                            np.linalg.LinAlgError,
                        ):
                            continue

                        found = candidate
                        effective_radius = trial_radius
                        attempts_used = attempt
                        break

                    if found is not None:
                        break

                fam["skipped_structural"] += attempts_used
                if found is None:
                    raise RuntimeError(
                        f"Family {rank}: no feasible fine 30-D direction found "
                        "even after step backoff."
                    )

                radius = float(effective_radius)
                x = found

            cid = candidate_id(rank, source["branches"], x)
            if any(
                r.get("candidate_id") == cid and r.get("result", {}).get("status") != "interrupted"
                for r in history
            ):
                fam["skipped_structural"] += 1
                continue

            try:
                sv, lv, sm, lm, small_mech, large_mech = screen_pair(
                    x, branches=source["branches"], motion=motion
                )
            except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
                fam["skipped_structural"] += 1
                continue

            design = build_design(machine, small_mech, large_mech)
            reusable = [
                r for r in history
                if r.get("result", {}).get("status") == "converged"
                and r.get("last_complete_state") is not None
            ]

            warm_distance = None
            if reusable:
                warm = min(
                    reusable,
                    key=lambda r: normalized_pair_distance(
                        x,
                        np.array([r["parameters"][name] for name in PAIR_NAMES], dtype=float),
                    ),
                )
                initial = np.asarray(warm["last_complete_state"], float)
                warm_source = warm["candidate_id"]
                warm_distance = normalized_pair_distance(
                    x,
                    np.array([warm["parameters"][name] for name in PAIR_NAMES], dtype=float),
                )
            else:
                initial = state3952
                warm_source = "candidate_3952_periodic_state"

            thermo = evaluate_thermo(
                label=f"sixbar_pair_3952_rank{rank}_{proposal_index}",
                design=design,
                definition=definition,
                initial=initial,
                deadline=deadline,
                candidate_seconds=args.candidate_seconds,
            )

            parameters = {name: float(value) for name, value in zip(PAIR_NAMES, x)}
            small_dict = mechanism_candidate_dict(
                sv, source["branches"]["small_primary"], source["branches"]["small_second"], sm
            )
            large_dict = mechanism_candidate_dict(
                lv, source["branches"]["large_primary"], source["branches"]["large_second"], lm
            )

            rec = {
                "index": len(history),
                "proposal_index": proposal_index,
                "family_rank": rank,
                "candidate_id": cid,
                "kind": kind,
                "radius_index": radius_index,
                "radius_multiplier": radius,
                "parameters": parameters,
                "branches": source["branches"],
                "mechanical": {
                    "small": small_dict,
                    "large": large_dict,
                    "combined_position_rms_vs_3952": float(
                        math.sqrt(0.5 * (
                            sm["position_rms_full_cycle"] ** 2
                            + lm["position_rms_full_cycle"] ** 2
                        ))
                    ),
                },
                "result": thermo["compact"],
                "feasible": thermo["feasible"],
                "physical_constraint_failures": thermo["reasons"],
                "elapsed_seconds": thermo["elapsed_seconds"],
                "warm_start_source": warm_source,
                "warm_start_distance": warm_distance,
                "warm_start_mode": thermo["warm_mode"],
                "safe_retry_used": thermo["safe_retry"],
                "backend": thermo["backend"],
                "last_complete_state": thermo["state"].tolist() if thermo["state"] is not None else None,
            }

            append_jsonl(fam["files"]["history"], rec)
            history.append(rec)
            completed_this_run[rank] += 1
            made_progress = True

            moved = False
            if eligible(rec) and (fam["best"] is None or eta(rec) > eta(fam["best"])):
                fam["best"] = rec
                moved = True
                save_best_pair(fam["files"]["best_pair"], rec)

            write_family_report(
                fam["files"]["report"],
                rank=rank,
                history=history,
                source=source,
                identity=fixed_identity,
            )

            print(json.dumps({
                "family_rank": rank,
                "index": rec["index"],
                "proposal_index": proposal_index,
                "kind": kind,
                "radius": radius,
                "status": rec["result"].get("status"),
                "feasible": rec["feasible"],
                "efficiency": rec["result"].get("indicated_thermal_efficiency"),
                "power_W": rec["result"].get("indicated_power_w"),
                "motion_rms_combined_pct": 100.0 * rec["mechanical"]["combined_position_rms_vs_3952"],
                "small_rms_pct": 100.0 * sm["position_rms_full_cycle"],
                "large_rms_pct": 100.0 * lm["position_rms_full_cycle"],
                "center_moved": moved,
                "best_efficiency": eta(fam["best"]) if fam["best"] is not None else None,
                "skipped_structural": fam["skipped_structural"],
                "elapsed_seconds": rec["elapsed_seconds"],
            }, separators=(",", ":")), flush=True)

        if not made_progress:
            break

    summary = {
        "description": "Fixed candidate-3952 260 K thermo; simultaneous 30-D optimization of paired six-bar mechanisms.",
        "fixed_260K_basis": fixed_identity,
        "requested_evaluations_per_family": args.evaluations_per_family,
        "completed_this_run": completed_this_run,
        "actual_duration_seconds": time.monotonic() - start,
        "families": {},
    }

    champions = []
    for rank in ranks:
        fam = families[rank]
        best = write_family_report(
            fam["files"]["report"],
            rank=rank,
            history=fam["history"],
            source=fam["source"],
            identity=fixed_identity,
        )
        summary["families"][str(rank)] = {
            "attempted_total": len(fam["history"]),
            "skipped_structural_total": fam["skipped_structural"],
            "best_feasible": best,
            "history": str(fam["files"]["history"]),
            "report": str(fam["files"]["report"]),
            "best_pair": str(fam["files"]["best_pair"]),
        }
        if best is not None:
            champions.append(best)

    summary["best_overall"] = max(champions, key=eta) if champions else None
    (output / "report.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Saved {output / 'report.json'}", flush=True)


if __name__ == "__main__":
    main()
