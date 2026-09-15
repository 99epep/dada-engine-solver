"""Local refinement of the seven-variable four-stage K2 motor search.

This final refinement reads both the original four-stage global-search history
and the preceding local-refinement history.  Neither source history is modified.
New evaluations are written under outputs/motor_four_stage_k2_refine_final/.

The default fine-radius schedule is 0.0025 -> 0.001.
"""
from __future__ import annotations

from dataclasses import replace
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.stats import qmc

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.four_stage_kinematics import FourStageVolumeKinematics
from dada_solver.integration import IntegrationInterrupted

from compare_motor_motion_laws_stage7A5 import (
    A5_CAMPAIGN,
    ROOT,
    _candidate_design_and_mass,
    _evaluate,
    _same_inventory_design,
)
from optimize_motor_exchanger_asymmetry_stageK1 import _scaled_exchanger
from optimize_motor_four_stage_k2 import NAMES, _symmetry_diagnostics
from optimize_motor_piecewise_stageP3 import _feasibility


SOURCE_DIRECTORY = ROOT / "outputs" / "motor_four_stage_k2"
SOURCE_REFINEMENT_DIRECTORY = ROOT / "outputs" / "motor_four_stage_k2_refine"
DEFAULT_DIRECTORY = ROOT / "outputs" / "motor_four_stage_k2_refine_final"
MINIMUM_STAGE_FRACTION = 0.02
DEFAULT_RADII = (0.0025, 0.001)


def _load_history(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line]


def _parameter_vector(record: dict) -> np.ndarray:
    return np.array([record["parameters"][name] for name in NAMES], dtype=float)


def _efficiency(record: dict) -> float:
    return float(record["result"]["indicated_thermal_efficiency"])


def _eligible(record: dict) -> bool:
    return bool(
        record.get("feasible")
        and record.get("result", {}).get("status") == "converged"
        and record.get("result", {}).get("indicated_thermal_efficiency") is not None
    )


def _valid_parameters(parameters: np.ndarray) -> bool:
    if np.any(parameters[3:] < 0.0) or np.any(parameters[3:] > 1.0):
        return False
    stage_durations = np.diff(np.r_[0.0, parameters[:3], 1.0])
    return bool(np.min(stage_durations) >= MINIMUM_STAGE_FRACTION)


def _select_diverse_elite(
    records: list[dict],
    count: int,
    minimum_distance: float,
) -> list[dict]:
    """Select high-efficiency centers while retaining distinct local basins."""
    ranked = sorted(
        (record for record in records if _eligible(record)),
        key=_efficiency,
        reverse=True,
    )
    selected: list[dict] = []
    for record in ranked:
        vector = _parameter_vector(record)
        if all(
            np.linalg.norm(vector - _parameter_vector(other)) >= minimum_distance
            for other in selected
        ):
            selected.append(record)
            if len(selected) >= count:
                break

    # If the distance filter leaves too few centers, complete with the best
    # remaining records rather than silently reducing the requested pool.
    if len(selected) < count:
        selected_ids = {record["candidate_id"] for record in selected}
        for record in ranked:
            if record["candidate_id"] not in selected_ids:
                selected.append(record)
                selected_ids.add(record["candidate_id"])
                if len(selected) >= count:
                    break
    return selected


def _parse_radii(text: str) -> tuple[float, ...]:
    radii = tuple(float(value) for value in text.split(",") if value.strip())
    if not radii or any(radius <= 0.0 for radius in radii):
        raise argparse.ArgumentTypeError("radii must be positive comma-separated values")
    if any(second >= first for first, second in zip(radii, radii[1:])):
        raise argparse.ArgumentTypeError("radii must be strictly decreasing")
    return radii


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget-seconds", type=float, default=3600)
    parser.add_argument("--evaluations", type=int, default=128)
    parser.add_argument(
        "--candidate-seconds",
        type=float,
        default=90,
        help="Cooperative per-candidate wall-clock limit.",
    )
    parser.add_argument(
        "--radii",
        type=_parse_radii,
        default=DEFAULT_RADII,
        help="Decreasing local-radius schedule; default: 0.0025,0.001.",
    )
    parser.add_argument(
        "--evaluations-per-radius",
        type=int,
        default=64,
        help="Number of recorded refinement evaluations before reducing radius.",
    )
    parser.add_argument(
        "--elite-centers",
        type=int,
        default=4,
        help="Number of distinct high-efficiency centers maintained during refinement.",
    )
    parser.add_argument(
        "--elite-distance",
        type=float,
        default=0.003,
        help="Minimum Euclidean distance in the seven normalized parameters between elite centers.",
    )
    parser.add_argument("--seed", type=int, default=260916)
    parser.add_argument("--output-directory", type=Path, default=DEFAULT_DIRECTORY)
    args = parser.parse_args()

    if args.budget_seconds <= 0.0:
        raise ValueError("--budget-seconds must be positive.")
    if args.evaluations <= 0:
        raise ValueError("--evaluations must be positive.")
    if args.evaluations_per_radius <= 0:
        raise ValueError("--evaluations-per-radius must be positive.")
    if args.elite_centers <= 0:
        raise ValueError("--elite-centers must be positive.")
    if args.elite_distance < 0.0:
        raise ValueError("--elite-distance must be non-negative.")

    source_history_path = SOURCE_DIRECTORY / "history.jsonl"
    if not source_history_path.exists():
        raise FileNotFoundError(
            f"Missing global-search history: {source_history_path.relative_to(ROOT)}"
        )
    source_history = _load_history(source_history_path)
    source_eligible = [record for record in source_history if _eligible(record)]
    if not source_eligible:
        raise RuntimeError("The global four-stage search contains no feasible candidate.")

    source_refinement_history_path = SOURCE_REFINEMENT_DIRECTORY / "history.jsonl"
    if not source_refinement_history_path.exists():
        raise FileNotFoundError(
            "Missing preceding refinement history: "
            f"{source_refinement_history_path.relative_to(ROOT)}"
        )
    source_refinement_history = _load_history(source_refinement_history_path)
    source_refinement_eligible = [
        record for record in source_refinement_history if _eligible(record)
    ]
    if not source_refinement_eligible:
        raise RuntimeError("The preceding refinement contains no feasible candidate.")

    directory = args.output_directory
    directory.mkdir(parents=True, exist_ok=True)
    history_path = directory / "history.jsonl"
    refinement_history = _load_history(history_path)

    reference_path = ROOT / "outputs" / "motor_champion_sixbar_k2.json"
    reference = json.loads(reference_path.read_text())
    definition = CampaignDefinition(A5_CAMPAIGN)
    base, mass, _, _, _ = _candidate_design_and_mass(definition)
    basis = reference["thermodynamic_basis"]
    base = replace(
        base,
        heat_in=_scaled_exchanger(base.heat_in, basis["stage_k2_k_i"]),
        heat_out=_scaled_exchanger(base.heat_out, basis["stage_k2_k_o"]),
    )
    limits = base.configuration.machine_volumes

    identity = {
        "names": list(NAMES),
        "source_history_sha256": hashlib.sha256(source_history_path.read_bytes()).hexdigest(),
        "source_refinement_history_sha256": hashlib.sha256(
            source_refinement_history_path.read_bytes()
        ).hexdigest(),
        "source_reference_sha256": hashlib.sha256(reference_path.read_bytes()).hexdigest(),
        "minimum_stage_fraction": MINIMUM_STAGE_FRACTION,
        "radii": list(args.radii),
        "evaluations_per_radius": args.evaluations_per_radius,
        "elite_centers": args.elite_centers,
        "elite_distance": args.elite_distance,
        "sobol_seed": args.seed,
        "strategy": (
            "Final pure local Sobol refinement around elite points from the global search, "
            "preceding refinement and current final-refinement history; no global "
            "candidates and no objective penalty."
        ),
    }
    definition_path = directory / "definition.json"
    if definition_path.exists():
        saved_identity = json.loads(definition_path.read_text())
        if saved_identity != identity:
            raise ValueError(
                "Refinement definition changed; use a new output directory or restore "
                "the previous command-line settings."
            )
    else:
        definition_path.write_text(json.dumps(identity, indent=2) + "\n")

    source_all_history = source_history + source_refinement_history
    all_history = source_all_history + refinement_history
    base_best = max(source_eligible, key=_efficiency)
    previous_refinement_best = max(source_refinement_eligible, key=_efficiency)
    current_eligible = [record for record in all_history if _eligible(record)]
    best = max(current_eligible, key=_efficiency)
    phase_start_best = best

    sampler = qmc.Sobol(d=len(NAMES), scramble=True, seed=args.seed)
    points = sampler.random_base2(14)
    existing_ids = {record["candidate_id"] for record in all_history}

    next_proposal_index = max(
        (record.get("proposal_index", -1) for record in refinement_history),
        default=-1,
    ) + 1
    completed_this_run = 0
    start = time.monotonic()
    deadline = start + args.budget_seconds

    def progress(_: object) -> None:
        if time.monotonic() >= deadline:
            raise IntegrationInterrupted("Refinement wall-clock budget exhausted.")
        if time.monotonic() >= candidate_deadline:
            raise IntegrationInterrupted(
                "Candidate wall-clock budget exhausted; result remains unknown."
            )

    while completed_this_run < args.evaluations and time.monotonic() < deadline:
        refinement_count = len(refinement_history)
        radius_index = min(
            refinement_count // args.evaluations_per_radius,
            len(args.radii) - 1,
        )
        radius = args.radii[radius_index]

        elite = _select_diverse_elite(
            source_all_history + refinement_history,
            args.elite_centers,
            args.elite_distance,
        )
        if not elite:
            raise RuntimeError("No feasible elite center is available for refinement.")

        proposal_index = next_proposal_index
        next_proposal_index += 1
        if proposal_index >= len(points):
            raise RuntimeError(
                "Refinement exhausted its 16384 deterministic Sobol proposals; "
                "change --seed or start a new output directory."
            )
        center = elite[proposal_index % len(elite)]
        center_vector = _parameter_vector(center)
        unit = points[proposal_index]
        parameters = center_vector + (2.0 * unit - 1.0) * radius

        if not _valid_parameters(parameters):
            continue

        mapping = {
            name: float(value)
            for name, value in zip(NAMES, parameters, strict=True)
        }
        candidate_id = hashlib.sha256(
            json.dumps(mapping, sort_keys=True).encode()
        ).hexdigest()
        if candidate_id in existing_ids:
            continue

        kinematics = FourStageVolumeKinematics(
            limits.small_cylinder,
            limits.large_cylinder,
            **mapping,
        )
        design = _same_inventory_design(base, mass, kinematics)

        reusable = [
            record
            for record in source_all_history + refinement_history
            if record.get("result", {}).get("status") == "converged"
            and record.get("last_complete_state")
        ]
        source = (
            min(
                reusable,
                key=lambda record: np.linalg.norm(
                    parameters - _parameter_vector(record)
                ),
            )
            if reusable
            else None
        )
        initial_state = (
            source["last_complete_state"]
            if source is not None
            else reference["last_complete_state"]
        )

        before = time.monotonic()
        candidate_deadline = min(deadline, before + args.candidate_seconds)
        try:
            result, state = _evaluate(
                f"four_stage_refine_{len(refinement_history)}",
                design,
                definition,
                initial_state=np.array(initial_state),
                progress_callback=progress,
            )
        except (ValueError, RuntimeError) as exc:
            result = {"status": "integration_failure", "message": str(exc)}
            state = None

        feasible, reasons = (
            _feasibility(result)
            if result.get("status") == "converged"
            else (False, [])
        )
        record = {
            "index": len(refinement_history),
            "proposal_index": proposal_index,
            "candidate_id": candidate_id,
            "kind": "local_refinement_sobol",
            "radius": radius,
            "radius_index": radius_index,
            "center_candidate_id": center["candidate_id"],
            "center_efficiency": _efficiency(center),
            "parameters": mapping,
            "symmetry": _symmetry_diagnostics(mapping),
            "result": result,
            "feasible": feasible,
            "physical_constraint_failures": reasons,
            "elapsed_seconds": time.monotonic() - before,
            "warm_start_source": (
                source["candidate_id"] if source is not None else "sixbar_reference"
            ),
            "last_complete_state": state.tolist() if state is not None else None,
        }
        with history_path.open("a") as stream:
            stream.write(json.dumps(record) + "\n")
            stream.flush()

        refinement_history.append(record)
        existing_ids.add(candidate_id)
        completed_this_run += 1

        if feasible and _efficiency(record) > _efficiency(best):
            best = record

        refinement_eligible = [
            candidate for candidate in refinement_history if _eligible(candidate)
        ]
        best_refinement = (
            max(refinement_eligible, key=_efficiency)
            if refinement_eligible
            else None
        )
        report = {
            "best_overall": best,
            "best_refinement": best_refinement,
            "best_global_search": base_best,
            "best_previous_refinement": previous_refinement_best,
            "best_symmetry": _symmetry_diagnostics(best["parameters"]),
            "phase_start_best": phase_start_best,
            "refinement_attempted_total": len(refinement_history),
            "current_radius": radius,
            "current_radius_index": radius_index,
            "requested_duration_seconds": args.budget_seconds,
            "actual_duration_seconds": time.monotonic() - start,
            "reference_efficiency": basis["reference_efficiency"],
            "reference_sixbar_efficiency": reference["result"][
                "indicated_thermal_efficiency"
            ],
        }
        (directory / "report.json").write_text(json.dumps(report, indent=2) + "\n")

        print(
            json.dumps(
                {
                    "refinement_index": record["index"],
                    "proposal_index": proposal_index,
                    "radius": radius,
                    "center_efficiency": _efficiency(center),
                    "status": result.get("status"),
                    "feasible": feasible,
                    "efficiency": result.get("indicated_thermal_efficiency"),
                    "power_W": result.get("indicated_power_w"),
                    "best_efficiency": _efficiency(best),
                    "symmetry_rms": record["symmetry"]["rms"],
                    "elapsed_seconds": record["elapsed_seconds"],
                }
            ),
            flush=True,
        )

    print(f"Saved {directory / 'report.json'}", flush=True)


if __name__ == "__main__":
    main()
