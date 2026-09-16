"""Adaptive 7D local optimization starting only from 150/30/150/30.

Hardware is frozen to the retained K2 design:
    H_i: 3708 tubes
    H_o: 3708 tubes
Charge is 100 kPa absolute at theta=0 for every candidate.

The ONLY initial search center is:
    LP exchange 150 deg
    compression  30 deg
    HP exchange 150 deg
    expansion     30 deg

All seven kinematic coordinates vary:
    t1, t2, t3, a_l, b_l, a_s, b_s

No historical 7D/9D optimum is ever used as a search center or warm start.
Warm starts may only reuse states generated inside this campaign.

The search uses anisotropic incumbent-centered Sobol refinement. Timing is
allowed to move farther than the intermediate volume levels at first, so the
optimizer can walk toward a very different chronology if that is favorable.

Default radius schedule (timing, levels):
    (0.15,   0.06)
    (0.08,   0.04)
    (0.04,   0.02)
    (0.02,   0.01)
    (0.01,   0.005)
    (0.005,  0.0025)
    (0.0025, 0.00125)
    (0.001,  0.0005)

Run:
    PYTHONPATH=src python3 examples/refine_motor_four_stage_from_150_30.py \
        --budget-seconds 21600 --evaluations 512 --candidate-seconds 120

Outputs:
    outputs/motor_four_stage_from_150_30/history.jsonl
    outputs/motor_four_stage_from_150_30/report.json
    outputs/motor_four_stage_from_150_30/definition.json
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
from dada_solver.configuration import ChargeConfiguration
from dada_solver.four_stage_kinematics import FourStageVolumeKinematics
from dada_solver.integration import IntegrationInterrupted

from compare_motor_motion_laws_stage7A5 import (
    A5_CAMPAIGN,
    ROOT,
    _candidate_design_and_mass,
    _evaluate,
    _uniform_wall_initial,
)
from optimize_motor_exchanger_asymmetry_stageK1 import _scaled_exchanger
from optimize_motor_four_stage_k2 import NAMES, _symmetry_diagnostics
from optimize_motor_piecewise_stageP3 import _feasibility


THERMO3D_REPORT = ROOT / "outputs" / "motor_four_stage_thermo3d" / "report.json"
DEFAULT_DIRECTORY = ROOT / "outputs" / "motor_four_stage_from_150_30"

CHARGE_PRESSURE_PA = 100_000.0
MINIMUM_STAGE_FRACTION = 0.02
SOBOL_SEED = 26091615

DEFAULT_RADII = (
    (0.15, 0.06),
    (0.08, 0.04),
    (0.04, 0.02),
    (0.02, 0.01),
    (0.01, 0.005),
    (0.005, 0.0025),
    (0.0025, 0.00125),
    (0.001, 0.0005),
)


def _load_history(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def _eligible(record: dict) -> bool:
    return bool(
        record.get("feasible")
        and record.get("result", {}).get("status") == "converged"
        and record.get("result", {}).get("indicated_thermal_efficiency") is not None
    )


def _efficiency(record: dict) -> float:
    return float(record["result"]["indicated_thermal_efficiency"])


def _parameter_vector(record: dict) -> np.ndarray:
    return np.asarray([record["parameters"][name] for name in NAMES], dtype=float)


def _valid_parameters(vector: np.ndarray) -> bool:
    if np.any(~np.isfinite(vector)):
        return False
    if np.any(vector[3:] < 0.0) or np.any(vector[3:] > 1.0):
        return False
    durations = np.diff(np.r_[0.0, vector[:3], 1.0])
    return bool(np.min(durations) >= MINIMUM_STAGE_FRACTION)


def _phase_degrees(parameters: dict) -> dict:
    t1 = float(parameters["t1"])
    t2 = float(parameters["t2"])
    t3 = float(parameters["t3"])
    return {
        "low_pressure_exchange_deg": 360.0 * t1,
        "compression_deg": 360.0 * (t2 - t1),
        "high_pressure_exchange_deg": 360.0 * (t3 - t2),
        "expansion_deg": 360.0 * (1.0 - t3),
    }


def _candidate_id(parameters: dict) -> str:
    payload = json.dumps(
        {name: float(parameters[name]) for name in NAMES},
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def _load_basis():
    report = json.loads(THERMO3D_REPORT.read_text())
    control = report.get("control_1bar")
    if control is None or not control.get("feasible"):
        raise RuntimeError("Missing feasible 1-bar thermo-3D K2 control.")

    definition = CampaignDefinition(A5_CAMPAIGN)
    base, _legacy_mass, _saved, _best, _candidate = _candidate_design_and_mass(definition)

    heat_in = _scaled_exchanger(base.heat_in, float(control["parameters"]["k_i"]))
    heat_out = _scaled_exchanger(base.heat_out, float(control["parameters"]["k_o"]))

    if heat_in.bank.tube_count != 3708 or heat_out.bank.tube_count != 3708:
        raise RuntimeError("Expected retained K2 hardware to use 3708 tubes.")

    charge = ChargeConfiguration(
        temperature=base.configuration.charge.temperature,
        pressure=CHARGE_PRESSURE_PA,
    )
    base = replace(
        base,
        configuration=replace(base.configuration, charge=charge),
        heat_in=heat_in,
        heat_out=heat_out,
    )
    return definition, base, control


def _seed(control: dict) -> dict:
    retained = control["motion_parameters"]
    return {
        "t1": 150.0 / 360.0,
        "t2": 180.0 / 360.0,
        "t3": 330.0 / 360.0,
        "a_l": float(retained["a_l"]),
        "b_l": float(retained["b_l"]),
        "a_s": float(retained["a_s"]),
        "b_s": float(retained["b_s"]),
    }


def _build_design(base, parameters: dict):
    limits = base.configuration.machine_volumes
    kin = FourStageVolumeKinematics(
        limits.small_cylinder,
        limits.large_cylinder,
        **{name: float(parameters[name]) for name in NAMES},
    )
    return replace(base, kinematics=kin)


def _uniform_initial(design) -> np.ndarray:
    wrapper = design.build()
    return np.asarray(_uniform_wall_initial(design.configuration, wrapper), dtype=float)


def _rescaled_warm_start(source: dict, target_uniform: np.ndarray) -> np.ndarray:
    state = np.asarray(source["last_complete_state"], dtype=float).copy()
    if state.shape != (10,):
        raise ValueError("Expected ten-state dynamic-wall warm start.")
    source_mass = float(np.sum(state[:8:2]))
    target_mass = float(np.sum(target_uniform[:8:2]))
    if source_mass <= 0.0 or target_mass <= 0.0:
        raise ValueError("Warm-start inventories must be positive.")
    state[:8] *= target_mass / source_mass
    return state


def _parse_radii(text: str) -> tuple[tuple[float, float], ...]:
    result = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        timing, levels = item.split(":")
        pair = (float(timing), float(levels))
        if pair[0] <= 0 or pair[1] <= 0:
            raise argparse.ArgumentTypeError("Radii must be positive.")
        result.append(pair)
    if not result:
        raise argparse.ArgumentTypeError("At least one radius pair is required.")
    return tuple(result)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget-seconds", type=float, default=21600.0)
    parser.add_argument("--evaluations", type=int, default=512)
    parser.add_argument("--candidate-seconds", type=float, default=120.0)
    parser.add_argument("--evaluations-per-radius", type=int, default=64)
    parser.add_argument(
        "--radii",
        type=_parse_radii,
        default=DEFAULT_RADII,
        help="Comma-separated timing:level pairs.",
    )
    parser.add_argument("--seed", type=int, default=SOBOL_SEED)
    parser.add_argument("--output-directory", type=Path, default=DEFAULT_DIRECTORY)
    args = parser.parse_args()

    if args.budget_seconds <= 0 or args.evaluations <= 0 or args.candidate_seconds <= 0:
        raise ValueError("Budgets and evaluation counts must be positive.")
    if args.evaluations_per_radius <= 0:
        raise ValueError("--evaluations-per-radius must be positive.")

    definition, base, control = _load_basis()
    seed_parameters = _seed(control)

    directory = args.output_directory
    directory.mkdir(parents=True, exist_ok=True)
    history_path = directory / "history.jsonl"
    definition_path = directory / "definition.json"

    identity = {
        "names": list(NAMES),
        "only_search_seed": "150/30/150/30 with retained champion a/b levels",
        "hardware": {
            "H_i_tube_count": 3708,
            "H_o_tube_count": 3708,
            "H_i_tube_length_m": base.heat_in.bank.tube_length_m,
            "H_o_tube_length_m": base.heat_out.bank.tube_length_m,
            "tube_inner_diameter_m": base.heat_in.bank.inner_diameter_m,
        },
        "charge_pressure_pa": CHARGE_PRESSURE_PA,
        "minimum_stage_fraction": MINIMUM_STAGE_FRACTION,
        "radii": [list(pair) for pair in args.radii],
        "evaluations_per_radius": args.evaluations_per_radius,
        "sobol_seed": args.seed,
        "seed_parameters": seed_parameters,
        "strategy": (
            "Adaptive incumbent-centered local Sobol climb. The only initial "
            "center is the 150/30/150/30 seed; no historical optimum is used."
        ),
    }

    if definition_path.exists():
        if json.loads(definition_path.read_text()) != identity:
            raise ValueError("Definition changed; use a new output directory.")
    else:
        definition_path.write_text(json.dumps(identity, indent=2) + "\n")

    history = _load_history(history_path)
    eligible = [record for record in history if _eligible(record)]
    best = max(eligible, key=_efficiency) if eligible else None
    existing_finished = {
        record["candidate_id"]
        for record in history
        if record.get("result", {}).get("status") != "interrupted"
    }

    sobol = qmc.Sobol(d=7, scramble=True, seed=args.seed)
    points = sobol.random_base2(14)

    next_proposal = max(
        (record.get("proposal_index", -1) for record in history),
        default=-1,
    ) + 1

    start = time.monotonic()
    deadline = start + args.budget_seconds
    completed_this_run = 0

    def progress(_: object) -> None:
        if time.monotonic() >= deadline:
            raise IntegrationInterrupted("Campaign wall-clock budget exhausted.")
        if time.monotonic() >= candidate_deadline:
            raise IntegrationInterrupted("Candidate wall-clock budget exhausted.")

    while completed_this_run < args.evaluations and time.monotonic() < deadline:
        proposal_index = next_proposal
        next_proposal += 1

        if not history and proposal_index == 0:
            parameters = dict(seed_parameters)
            kind = "150_30_seed"
            radius_index = -1
            timing_radius = None
            level_radius = None
        else:
            nonseed_count = sum(1 for r in history if r["kind"] != "150_30_seed")
            radius_index = min(
                nonseed_count // args.evaluations_per_radius,
                len(args.radii) - 1,
            )
            timing_radius, level_radius = args.radii[radius_index]
            center = best["parameters"] if best is not None else seed_parameters
            center_vector = np.asarray([center[name] for name in NAMES], dtype=float)

            sobol_index = proposal_index - 1
            if sobol_index >= len(points):
                raise RuntimeError("Exhausted deterministic Sobol bank.")
            u = points[sobol_index]
            radii = np.asarray([timing_radius] * 3 + [level_radius] * 4)
            vector = center_vector + (2.0 * u - 1.0) * radii

            if not _valid_parameters(vector):
                continue

            parameters = {
                name: float(value)
                for name, value in zip(NAMES, vector, strict=True)
            }
            kind = "adaptive_local_sobol_from_150_30"

        candidate_id = _candidate_id(parameters)
        if candidate_id in existing_finished:
            continue

        design = _build_design(base, parameters)
        target_uniform = _uniform_initial(design)

        reusable = [
            record
            for record in history
            if record.get("result", {}).get("status") == "converged"
            and record.get("last_complete_state") is not None
        ]
        source = None
        initial_state = target_uniform
        if reusable:
            target_vector = np.asarray([parameters[name] for name in NAMES], dtype=float)
            source = min(
                reusable,
                key=lambda record: np.linalg.norm(
                    target_vector - _parameter_vector(record)
                ),
            )
            initial_state = _rescaled_warm_start(source, target_uniform)

        before = time.monotonic()
        candidate_deadline = min(deadline, before + args.candidate_seconds)

        try:
            result, state = _evaluate(
                f"four_stage_from_150_30_{proposal_index}",
                design,
                definition,
                initial_state=np.asarray(initial_state, dtype=float),
                progress_callback=progress,
            )
        except IntegrationInterrupted as exc:
            result = {"status": "interrupted", "message": str(exc)}
            state = None
        except (ValueError, RuntimeError) as exc:
            result = {"status": "integration_failure", "message": str(exc)}
            state = None

        feasible, reasons = (
            _feasibility(result)
            if result.get("status") == "converged"
            else (False, [])
        )

        record = {
            "index": len(history),
            "proposal_index": proposal_index,
            "candidate_id": candidate_id,
            "kind": kind,
            "radius_index": radius_index,
            "timing_radius": timing_radius,
            "level_radius": level_radius,
            "parameters": parameters,
            "phase_degrees": _phase_degrees(parameters),
            "symmetry": _symmetry_diagnostics(parameters),
            "result": result,
            "feasible": feasible,
            "physical_constraint_failures": reasons,
            "elapsed_seconds": time.monotonic() - before,
            "warm_start_source": source["candidate_id"] if source is not None else None,
            "last_complete_state": state.tolist() if state is not None else None,
        }

        with history_path.open("a") as stream:
            stream.write(json.dumps(record) + "\n")
            stream.flush()

        history.append(record)
        completed_this_run += 1
        if result.get("status") != "interrupted":
            existing_finished.add(candidate_id)

        moved = False
        if feasible and (best is None or _efficiency(record) > _efficiency(best)):
            best = record
            moved = True

        seed_record = next((r for r in history if r["kind"] == "150_30_seed"), None)

        report = {
            "best_feasible": best,
            "seed": seed_record,
            "attempted_total": len(history),
            "completed_this_run": completed_this_run,
            "requested_duration_seconds": args.budget_seconds,
            "actual_duration_seconds": time.monotonic() - start,
            "definition": identity,
            "current_radius_index": radius_index,
            "current_timing_radius": timing_radius,
            "current_level_radius": level_radius,
            "best_phase_degrees": best["phase_degrees"] if best is not None else None,
        }
        (directory / "report.json").write_text(json.dumps(report, indent=2) + "\n")

        print(
            json.dumps(
                {
                    "index": record["index"],
                    "radius": [timing_radius, level_radius],
                    "status": result.get("status"),
                    "feasible": feasible,
                    "efficiency": result.get("indicated_thermal_efficiency"),
                    "power_W": result.get("indicated_power_w"),
                    "candidate_phases_deg": record["phase_degrees"],
                    "center_moved": moved,
                    "best_efficiency": _efficiency(best) if best is not None else None,
                    "best_phases_deg": best["phase_degrees"] if best is not None else None,
                }
            ),
            flush=True,
        )

    print(f"Saved {directory / 'report.json'}", flush=True)


if __name__ == "__main__":
    main()
