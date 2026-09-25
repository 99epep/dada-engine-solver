#!/usr/bin/env python3
"""Thermo-5D retuning of four fixed six-bar pairs around candidate 3952.

Purpose
-------
The mechanical pairs have already been synthesized and thermodynamically
polished in 30D.  This campaign freezes each pair completely and asks a
different question:

    How much efficiency is recovered if the 260 K machine is allowed to
    re-tune only its five established thermo/hardware coordinates?

Free variables, per fixed mechanism pair:
    swept_ratio = V_swept,S / V_swept,L
    n_i         = H_i tube count
    length_i_m  = H_i tube length
    n_o         = H_o tube count
    length_o_m  = H_o tube length

Fixed:
    - exact SMALL and LARGE six-bar geometry and branches;
    - ΔT = 260 K and the same hot/cold source conditions as candidate 3952;
    - total swept volume and cylinder clearance ratios;
    - total working-gas inventory of candidate 3952;
    - exchanger model, tube ID, valve scaling policy, solver settings;
    - all other machine parameters.

Families:
    ranks 1, 4, 12, 50 all use their final champions from
    outputs/sixbar_thermo_coupled_3952_fine_hlat25/.

Each family has its own persistent history.  Families are evaluated round-robin.

Typical run:
    PYTHONPATH=src python3 examples/optimize_sixbar_pairs_thermo5d_3952.py \
      --evaluations-per-family 96 \
      --budget-seconds 14400

Rerunning resumes the histories.
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
from optimize_motor_four_stage_thermo3d import _volume_limits
from refine_motor_four_stage_hx9d_variable_gas import _hardware_metrics


ROOT = Path.cwd()

OUTPUT = ROOT / "outputs" / "sixbar_thermo5d_3952"

PAIR_SOURCES = {
    rank: (
        ROOT / "outputs" / "sixbar_thermo_coupled_3952_fine_hlat25"
        / f"rank_{rank:02d}" / "best_pair.json",
        ROOT / "outputs" / "sixbar_thermo_coupled_3952_fine_hlat25"
        / f"rank_{rank:02d}" / "report.json",
    )
    for rank in (1, 4, 12, 50)
}

PARAMS = ("swept_ratio", "n_i", "length_i_m", "n_o", "length_o_m")

POWER_FLOOR_W = 25.0

# Same broad physically explored thermo domain used by the temperature campaign.
RATIO_BOUNDS = (0.50, 1.30)
COUNT_BOUNDS = (2400, 5200)
HI_LENGTH_BOUNDS = (0.03175, 0.06350)
HO_LENGTH_BOUNDS = (0.016933333333333334, 0.035983333333333335)

# Same local 5D radius hierarchy used by optimize_motor_temperature_point.py.
DEFAULT_RADII = (
    (0.08, 240, 0.10),
    (0.04, 120, 0.05),
    (0.02, 60, 0.025),
    (0.01, 30, 0.0125),
)


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def load_history(path: Path):
    if not path.exists():
        return []
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def append_record(path: Path, record: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(record, separators=(",", ":")) + "\n")
        stream.flush()


def mechanism_from_candidate(candidate: dict) -> SixBarCylinderMechanism:
    p = candidate["parameters"]
    primary = candidate["primary"]
    return SixBarCylinderMechanism(
        primary_ground=float(p["primary_ground"]),
        primary_coupler=float(p["primary_coupler"]),
        primary_rocker=float(p["primary_rocker"]),
        primary_e_along=float(p["primary_E_along"]),
        primary_e_normal=float(p["primary_E_normal"]),
        primary_phase=float(p["primary_phase"]),
        second_pivot_x=float(p["second_pivot_x"]),
        second_pivot_y=float(p["second_pivot_y"]),
        link_ef=float(p["link_EF"]),
        link_gf=float(p["link_GF"]),
        h_along_over_ef=float(p["H_along_over_EF"]),
        h_normal_over_ef=float(p["H_normal_over_EF"]),
        piston_rod=float(p["piston_rod"]),
        slider_axis_offset=float(p["slider_axis_offset"]),
        slider_axis_angle=float(p["slider_axis_angle"]),
        primary_branch=int(primary["assembly_branch"]),
        second_branch=int(candidate["second_branch"]),
    )


def make_3952_basis():
    """Rebuild the exact fixed 260 K hardware basis used by candidate 3952."""
    sp = hybrid15.resolve(hybrid15.DEFAULT_SOURCE)
    fp = hybrid15.resolve(hybrid15.DEFAULT_FOURIER)

    source_report = load_json(sp)
    source_best = source_report.get("best_feasible")
    if source_best is None:
        raise RuntimeError("Source report has no best_feasible.")

    fourier_report = load_json(fp)
    fb = hybrid15._source_best(fourier_report)
    harmonics = hybrid15._source_harmonics(fb)

    campaign_definition, base = hybrid15._load_basis()
    geometry = hybrid15.base_geometry(base)

    machine = hybrid15.fourier_build_design(
        base,
        geometry,
        fb["thermo"],
        np.asarray(fb["small_coefficients"]),
        np.asarray(fb["large_coefficients"]),
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

    limits = machine.configuration.machine_volumes
    s = limits.small_cylinder
    l = limits.large_cylinder

    geometry_info = {
        "total_swept_m3": float(s.swept + l.swept),
        "small_clearance_ratio": float(s.minimum / s.swept),
        "large_clearance_ratio": float(l.minimum / l.swept),
        "seed_ni": int(machine.heat_in.bank.tube_count),
        "seed_no": int(machine.heat_out.bank.tube_count),
        "seed_li": float(machine.heat_in.bank.tube_length_m),
        "seed_lo": float(machine.heat_out.bank.tube_length_m),
        "seed_hi_cda": float(machine.heat_in.outlet_valve_cda_m2),
        "seed_ho_cda": float(machine.heat_out.outlet_valve_cda_m2),
    }

    thermo_seed = {
        "swept_ratio": float(fb["thermo"]["swept_ratio"]),
        "n_i": int(fb["thermo"]["n_i"]),
        "length_i_m": float(fb["thermo"]["length_i_m"]),
        "n_o": int(fb["thermo"]["n_o"]),
        "length_o_m": float(fb["thermo"]["length_o_m"]),
    }

    for name, actual, expected in (
        ("n_i", geometry_info["seed_ni"], thermo_seed["n_i"]),
        ("n_o", geometry_info["seed_no"], thermo_seed["n_o"]),
    ):
        if actual != expected:
            raise RuntimeError(
                f"3952 basis mismatch for {name}: machine={actual}, thermo={expected}"
            )

    identity = {
        "source_report": str(sp),
        "source_sha256": hashlib.sha256(sp.read_bytes()).hexdigest(),
        "fourier_report": str(fp),
        "fourier_sha256": hashlib.sha256(fp.read_bytes()).hexdigest(),
        "fixed_total_mass_kg": total_mass,
        "thermo_seed_3952": thermo_seed,
        "total_swept_m3": geometry_info["total_swept_m3"],
        "clearance_ratios": [
            geometry_info["small_clearance_ratio"],
            geometry_info["large_clearance_ratio"],
        ],
    }

    return machine, definition, total_mass, geometry_info, thermo_seed, identity


def valid_params(p: dict) -> bool:
    return bool(
        RATIO_BOUNDS[0] <= float(p["swept_ratio"]) <= RATIO_BOUNDS[1]
        and COUNT_BOUNDS[0] <= int(p["n_i"]) <= COUNT_BOUNDS[1]
        and COUNT_BOUNDS[0] <= int(p["n_o"]) <= COUNT_BOUNDS[1]
        and HI_LENGTH_BOUNDS[0] <= float(p["length_i_m"]) <= HI_LENGTH_BOUNDS[1]
        and HO_LENGTH_BOUNDS[0] <= float(p["length_o_m"]) <= HO_LENGTH_BOUNDS[1]
    )


def build_design(
    seed_machine,
    geometry_info: dict,
    total_mass: float,
    params: dict,
    small_mech: SixBarCylinderMechanism,
    large_mech: SixBarCylinderMechanism,
):
    """Build one thermo candidate with the S/L mechanisms strictly frozen."""
    s_lim, l_lim = _volume_limits(
        geometry_info["total_swept_m3"],
        float(params["swept_ratio"]),
        geometry_info["small_clearance_ratio"],
        geometry_info["large_clearance_ratio"],
    )

    config = replace(
        seed_machine.configuration,
        machine_volumes=replace(
            seed_machine.configuration.machine_volumes,
            small_cylinder=s_lim,
            large_cylinder=l_lim,
        ),
    )

    ni = int(params["n_i"])
    no = int(params["n_o"])

    hi = replace(
        seed_machine.heat_in,
        bank=replace(
            seed_machine.heat_in.bank,
            tube_count=ni,
            tube_length_m=float(params["length_i_m"]),
        ),
        outlet_valve_cda_m2=(
            geometry_info["seed_hi_cda"]
            * ni / geometry_info["seed_ni"]
        ),
    )

    ho = replace(
        seed_machine.heat_out,
        bank=replace(
            seed_machine.heat_out.bank,
            tube_count=no,
            tube_length_m=float(params["length_o_m"]),
        ),
        outlet_valve_cda_m2=(
            geometry_info["seed_ho_cda"]
            * no / geometry_info["seed_no"]
        ),
    )

    design = replace(
        seed_machine,
        configuration=config,
        heat_in=hi,
        heat_out=ho,
    )

    # Preserve exactly the candidate-3952 working-gas inventory even though
    # exchanger dead volume and S/L swept split change.
    design = hybrid15._fixed_inventory_design(design, total_mass)

    limits = design.configuration.machine_volumes
    kin = IndependentSixBarVolumeKinematics(
        small_mech,
        large_mech,
        limits.small_cylinder,
        limits.large_cylinder,
    )

    return replace(design, kinematics=kin)


def hardware(design):
    return {
        "H_i": _hardware_metrics(design.heat_in),
        "H_o": _hardware_metrics(design.heat_out),
    }


def warm_state(source_state, source_hw, target_hw):
    """Same gas inventory; only scale wall-energy states with wall capacity."""
    x = np.asarray(source_state, dtype=float).copy()
    if x.shape != (10,):
        raise ValueError("Expected ten-state periodic state.")

    for i, side in ((8, "H_i"), (9, "H_o")):
        old = float(source_hw[side]["wall_capacity_j_k"])
        new = float(target_hw[side]["wall_capacity_j_k"])
        if old <= 0 or new <= 0:
            raise ValueError("Invalid exchanger wall capacity.")
        x[i] *= new / old

    return x


def thermo_distance(p: dict, record: dict, seed: dict):
    q = record["parameters"]
    a = np.asarray([
        float(p["swept_ratio"]),
        float(p["n_i"]) / seed["n_i"],
        float(p["length_i_m"]) / seed["length_i_m"],
        float(p["n_o"]) / seed["n_o"],
        float(p["length_o_m"]) / seed["length_o_m"],
    ])
    b = np.asarray([
        float(q["swept_ratio"]),
        float(q["n_i"]) / seed["n_i"],
        float(q["length_i_m"]) / seed["length_i_m"],
        float(q["n_o"]) / seed["n_o"],
        float(q["length_o_m"]) / seed["length_o_m"],
    ])
    return float(np.linalg.norm(a - b))


def candidate_id(rank: int, pair_id: str, params: dict):
    payload = {
        "family_rank": int(rank),
        "pair_candidate_id": pair_id,
        "parameters": {
            "swept_ratio": float(params["swept_ratio"]),
            "n_i": int(params["n_i"]),
            "length_i_m": float(params["length_i_m"]),
            "n_o": int(params["n_o"]),
            "length_o_m": float(params["length_o_m"]),
        },
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


def evaluate(
    *,
    label,
    design,
    definition,
    initial,
    deadline,
    candidate_seconds,
):
    started = time.monotonic()
    candidate_deadline = min(deadline, started + candidate_seconds)

    def progress(_):
        now = time.monotonic()
        if now >= deadline:
            raise IntegrationInterrupted("Thermo-5D campaign budget exhausted.")
        if now >= candidate_deadline:
            raise IntegrationInterrupted("Candidate budget exhausted.")

    box = {}

    def observe(periodic):
        box["statistics"] = periodic.backend_statistics

    safe_retry = False
    warm_mode = "scaled_wall_periodic_state"

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
            safe = hybrid15._safe_uniform_initial(
                design,
                wall_source=np.asarray(initial, float),
            )
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


def parse_radii(text: str):
    rows = []
    for item in text.split(","):
        if not item.strip():
            continue
        a, b, c = item.split(":")
        row = (float(a), int(b), float(c))
        if min(row) <= 0:
            raise argparse.ArgumentTypeError("All thermo radii must be positive.")
        rows.append(row)
    if not rows:
        raise argparse.ArgumentTypeError("At least one radius triplet is required.")
    return tuple(rows)


def load_pair(rank: int):
    pair_path, report_path = PAIR_SOURCES[rank]

    if not pair_path.exists():
        raise FileNotFoundError(pair_path)
    if not report_path.exists():
        raise FileNotFoundError(report_path)

    pair = load_json(pair_path)
    report = load_json(report_path)
    best = report.get("best_feasible")
    if best is None:
        raise RuntimeError(f"Family {rank} report has no best_feasible.")

    if pair.get("candidate_id") != best.get("candidate_id"):
        raise RuntimeError(
            f"Family {rank}: best_pair.json and report.json disagree on champion."
        )

    state = best.get("last_complete_state")
    if state is None:
        raise RuntimeError(f"Family {rank}: champion has no periodic state.")

    small = mechanism_from_candidate(pair["small"])
    large = mechanism_from_candidate(pair["large"])

    return {
        "pair_path": pair_path,
        "report_path": report_path,
        "pair_id": str(pair["candidate_id"]),
        "source_efficiency": float(pair["efficiency"]),
        "source_power_W": float(pair["power_W"]),
        "state": np.asarray(state, float),
        "small": small,
        "large": large,
    }


def family_paths(output: Path, rank: int):
    root = output / f"rank_{rank:02d}"
    return {
        "root": root,
        "history": root / "history.jsonl",
        "report": root / "report.json",
    }


def write_family_report(
    path: Path,
    rank: int,
    source: dict,
    history: list,
    definition_identity: dict,
    requested_nonseed: int,
):
    good = [r for r in history if eligible(r)]
    best = max(good, key=eta) if good else None
    n = sum(r.get("kind") != "fixed_pair_3952_thermo_seed" for r in history)

    report = {
        "family_rank": rank,
        "fixed_pair": {
            "path": str(source["pair_path"]),
            "candidate_id": source["pair_id"],
            "efficiency_before_thermo5d": source["source_efficiency"],
            "power_W_before_thermo5d": source["source_power_W"],
        },
        "definition": definition_identity,
        "attempted_total": len(history),
        "completed_nonseed": n,
        "requested_nonseed": requested_nonseed,
        "complete": n >= requested_nonseed,
        "best_feasible": best,
    }

    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--families", default="1,4,12,50")
    ap.add_argument("--output-directory", type=Path, default=OUTPUT)

    ap.add_argument("--evaluations-per-family", type=int, default=96)
    ap.add_argument("--evaluations-per-radius", type=int, default=24)
    ap.add_argument(
        "--radii",
        type=parse_radii,
        default=DEFAULT_RADII,
        help=(
            "comma-separated swept_ratio:tube_count:length_relative radii; "
            "default 0.08:240:0.10,0.04:120:0.05,0.02:60:0.025,0.01:30:0.0125"
        ),
    )

    ap.add_argument("--budget-seconds", type=float, default=14400.0)
    ap.add_argument("--candidate-seconds", type=float, default=180.0)
    ap.add_argument("--seed", type=int, default=3965000)
    args = ap.parse_args()

    if min(
        args.evaluations_per_family,
        args.evaluations_per_radius,
        args.budget_seconds,
        args.candidate_seconds,
    ) <= 0:
        raise ValueError("Budgets/counts must be positive.")

    ranks = [int(x) for x in args.families.split(",") if x.strip()]
    unknown = [x for x in ranks if x not in PAIR_SOURCES]
    if unknown:
        raise ValueError(f"Unknown families: {unknown}")

    (
        seed_machine,
        definition,
        total_mass,
        geometry_info,
        thermo_seed,
        basis_identity,
    ) = make_3952_basis()

    output = args.output_directory
    output.mkdir(parents=True, exist_ok=True)

    definition_identity = {
        "study": "fixed-sixbar-pair thermo5D retuning around candidate 3952",
        "fixed_delta_T_K": 260.0,
        "fixed_total_gas_mass_kg": total_mass,
        "free_parameters": list(PARAMS),
        "thermo_seed_3952": thermo_seed,
        "bounds": {
            "swept_ratio": list(RATIO_BOUNDS),
            "tube_count": list(COUNT_BOUNDS),
            "H_i_tube_length_m": list(HI_LENGTH_BOUNDS),
            "H_o_tube_length_m": list(HO_LENGTH_BOUNDS),
        },
        "radii": [list(x) for x in args.radii],
        "evaluations_per_radius": args.evaluations_per_radius,
        "power_floor_W": POWER_FLOOR_W,
        "strategy": (
            "persistent round-robin 5D incumbent-centred scrambled Sobol; "
            "each S/L six-bar pair is strictly fixed"
        ),
        "valve_policy": "outlet CdA scales linearly with tube count",
        "basis": basis_identity,
        "families": ranks,
        "seed": args.seed,
    }

    definition_path = output / "definition.json"
    if definition_path.exists():
        if load_json(definition_path) != definition_identity:
            raise ValueError(
                "Thermo-5D definition changed; use a new output directory."
            )
    else:
        definition_path.write_text(
            json.dumps(definition_identity, indent=2) + "\n",
            encoding="utf-8",
        )

    # Hardware metadata of the exact 3952 thermo point; used to scale wall
    # energy in the first warm start for each fixed pair.
    seed_hw = hardware(seed_machine)

    families = {}
    for j, rank in enumerate(ranks):
        source = load_pair(rank)
        paths = family_paths(output, rank)
        paths["root"].mkdir(parents=True, exist_ok=True)

        history = load_history(paths["history"])
        good = [r for r in history if eligible(r)]
        best = max(good, key=eta) if good else None

        points = qmc.Sobol(
            d=5,
            scramble=True,
            seed=args.seed + 10000 * j,
        ).random_base2(15)

        families[rank] = {
            "source": source,
            "paths": paths,
            "history": history,
            "best": best,
            "points": points,
        }

    started = time.monotonic()
    deadline = started + args.budget_seconds
    completed_this_run = {rank: 0 for rank in ranks}

    while time.monotonic() < deadline:
        made_progress = False

        for rank in ranks:
            if time.monotonic() >= deadline:
                break

            fam = families[rank]
            source = fam["source"]
            history = fam["history"]

            n = sum(
                r.get("kind") != "fixed_pair_3952_thermo_seed"
                for r in history
            )
            if n >= args.evaluations_per_family:
                continue

            proposal = max(
                (r.get("proposal_index", -1) for r in history),
                default=-1,
            ) + 1

            if not history:
                p = dict(thermo_seed)
                kind = "fixed_pair_3952_thermo_seed"
                ridx = -1
                rr = nr = lr = None
            else:
                ridx = min(
                    n // args.evaluations_per_radius,
                    len(args.radii) - 1,
                )
                rr, nr, lr = args.radii[ridx]

                center = (
                    fam["best"]["parameters"]
                    if fam["best"] is not None
                    else thermo_seed
                )
                u = fam["points"][(proposal - 1) % len(fam["points"])]

                p = {
                    "swept_ratio": float(
                        center["swept_ratio"] + (2*u[0]-1) * rr
                    ),
                    "n_i": int(round(
                        center["n_i"] + (2*u[1]-1) * nr
                    )),
                    "length_i_m": float(
                        center["length_i_m"] * (1 + (2*u[2]-1) * lr)
                    ),
                    "n_o": int(round(
                        center["n_o"] + (2*u[3]-1) * nr
                    )),
                    "length_o_m": float(
                        center["length_o_m"] * (1 + (2*u[4]-1) * lr)
                    ),
                }
                kind = "adaptive_local_thermo5d"

                if not valid_params(p):
                    # Out-of-domain Sobol point: advance proposal index without
                    # charging one requested evaluation.
                    dummy = {
                        "proposal_index": proposal,
                        "kind": "skipped_out_of_bounds",
                    }
                    # We do not append skipped points; shift the Sobol indexing
                    # deterministically by adding a tiny persistent marker file.
                    skip_path = fam["paths"]["root"] / "skipped.txt"
                    with skip_path.open("a", encoding="utf-8") as stream:
                        stream.write(f"{proposal}\n")
                    # Count past skipped indices when choosing the next point.
                    # History proposal indices alone cannot do this, so consume
                    # an in-bounds point below using a local forward scan.
                    found = None
                    for extra in range(1, 2048):
                        pi2 = proposal + extra
                        u2 = fam["points"][(pi2 - 1) % len(fam["points"])]
                        trial = {
                            "swept_ratio": float(
                                center["swept_ratio"] + (2*u2[0]-1) * rr
                            ),
                            "n_i": int(round(
                                center["n_i"] + (2*u2[1]-1) * nr
                            )),
                            "length_i_m": float(
                                center["length_i_m"]
                                * (1 + (2*u2[2]-1) * lr)
                            ),
                            "n_o": int(round(
                                center["n_o"] + (2*u2[3]-1) * nr
                            )),
                            "length_o_m": float(
                                center["length_o_m"]
                                * (1 + (2*u2[4]-1) * lr)
                            ),
                        }
                        if valid_params(trial):
                            proposal = pi2
                            p = trial
                            found = True
                            break
                    if not found:
                        raise RuntimeError(
                            f"Family {rank}: no in-bounds thermo proposal found."
                        )

            cid = candidate_id(rank, source["pair_id"], p)
            if any(
                r.get("candidate_id") == cid
                and r.get("result", {}).get("status") != "interrupted"
                for r in history
            ):
                continue

            design = build_design(
                seed_machine,
                geometry_info,
                total_mass,
                p,
                source["small"],
                source["large"],
            )
            target_hw = hardware(design)

            reusable = [
                r for r in history
                if r.get("result", {}).get("status") == "converged"
                and r.get("last_complete_state") is not None
                and r.get("hardware")
            ]

            if reusable:
                warm = min(
                    reusable,
                    key=lambda r: thermo_distance(p, r, thermo_seed),
                )
                initial = warm_state(
                    warm["last_complete_state"],
                    warm["hardware"],
                    target_hw,
                )
                warm_source = warm["candidate_id"]
                warm_distance = thermo_distance(p, warm, thermo_seed)
            else:
                initial = warm_state(
                    source["state"],
                    seed_hw,
                    target_hw,
                )
                warm_source = source["pair_id"]
                warm_distance = 0.0 if kind == "fixed_pair_3952_thermo_seed" else None

            ev = evaluate(
                label=f"sixbar_thermo5d_rank{rank}_{proposal}",
                design=design,
                definition=definition,
                initial=initial,
                deadline=deadline,
                candidate_seconds=args.candidate_seconds,
            )

            record = {
                "index": len(history),
                "proposal_index": proposal,
                "family_rank": rank,
                "candidate_id": cid,
                "pair_candidate_id": source["pair_id"],
                "kind": kind,
                "radius_index": ridx,
                "radii": {
                    "swept_ratio": rr,
                    "tube_count": nr,
                    "tube_length_relative": lr,
                },
                "parameters": {
                    "swept_ratio": float(p["swept_ratio"]),
                    "n_i": int(p["n_i"]),
                    "length_i_m": float(p["length_i_m"]),
                    "n_o": int(p["n_o"]),
                    "length_o_m": float(p["length_o_m"]),
                },
                "hardware": target_hw,
                "result": ev["compact"],
                "feasible": ev["feasible"],
                "physical_constraint_failures": ev["reasons"],
                "elapsed_seconds": ev["elapsed_seconds"],
                "warm_start_source": warm_source,
                "warm_start_distance": warm_distance,
                "warm_start_mode": ev["warm_mode"],
                "safe_retry_used": ev["safe_retry"],
                "backend": ev["backend"],
                "last_complete_state": (
                    ev["state"].tolist()
                    if ev["state"] is not None else None
                ),
            }

            append_record(fam["paths"]["history"], record)
            history.append(record)
            made_progress = True
            completed_this_run[rank] += 1

            moved = False
            if eligible(record) and (
                fam["best"] is None or eta(record) > eta(fam["best"])
            ):
                fam["best"] = record
                moved = True

            write_family_report(
                fam["paths"]["report"],
                rank,
                source,
                history,
                definition_identity,
                args.evaluations_per_family,
            )

            print(json.dumps({
                "family_rank": rank,
                "index": record["index"],
                "proposal_index": proposal,
                "kind": kind,
                "radius_index": ridx,
                "status": record["result"].get("status"),
                "feasible": record["feasible"],
                "efficiency":
                    record["result"].get("indicated_thermal_efficiency"),
                "power_W":
                    record["result"].get("indicated_power_w"),
                "swept_ratio": p["swept_ratio"],
                "n_i": p["n_i"],
                "length_i_mm": 1000.0 * p["length_i_m"],
                "n_o": p["n_o"],
                "length_o_mm": 1000.0 * p["length_o_m"],
                "center_moved": moved,
                "best_efficiency":
                    eta(fam["best"]) if fam["best"] is not None else None,
                "elapsed_seconds": record["elapsed_seconds"],
            }, separators=(",", ":")), flush=True)

        if not made_progress:
            break

    summary = {
        "description": (
            "Thermo-5D retuning of four fixed paired six-bar mechanisms "
            "around candidate 3952 at fixed ΔT=260 K."
        ),
        "definition": definition_identity,
        "requested_evaluations_per_family": args.evaluations_per_family,
        "completed_this_run": completed_this_run,
        "actual_duration_seconds": time.monotonic() - started,
        "families": {},
    }

    champions = []
    for rank in ranks:
        fam = families[rank]
        best = write_family_report(
            fam["paths"]["report"],
            rank,
            fam["source"],
            fam["history"],
            definition_identity,
            args.evaluations_per_family,
        )

        if best is not None:
            champions.append(best)

        summary["families"][str(rank)] = {
            "pair_candidate_id": fam["source"]["pair_id"],
            "efficiency_before_thermo5d":
                fam["source"]["source_efficiency"],
            "power_W_before_thermo5d":
                fam["source"]["source_power_W"],
            "attempted_total": len(fam["history"]),
            "best_feasible": best,
            "report": str(fam["paths"]["report"]),
            "history": str(fam["paths"]["history"]),
        }

    summary["best_overall"] = (
        max(champions, key=eta) if champions else None
    )

    (output / "report.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print(f"Saved {output / 'report.json'}", flush=True)


if __name__ == "__main__":
    main()
