"""Joint 9D optimization of four-stage motor kinematics and exchanger tube counts.

Optimization variables
----------------------
Kinematics (7):
    t1, t2, t3, a_l, b_l, a_s, b_s

Physical exchanger sizing (2):
    n_i : H_i tube count
    n_o : H_o tube count

The retained K2 tube lengths and 0.33 mm tube ID are held fixed. Varying tube
count therefore changes, through the existing physical exchanger model:
- gas-side and air-side area/conductance,
- tube flow area and hydraulic resistance,
- header frontal dimensions and header gas volume,
- tube gas hold-up / total exchanger dead volume,
- wall mass and heat capacity,
- external-air-side geometry.

Outlet-valve CdA is scaled in direct proportion to tube count relative to the
retained 3708-tube K2 hardware. This preserves the valve-area/tube-count ratio
instead of silently retaining a valve sized for another bank.

Every candidate is filled independently at 100 kPa absolute at theta=0.
Working-gas inventory is therefore NOT held fixed when motion or exchanger
dead volume changes.

Objective
---------
Maximize indicated thermal efficiency, subject to the repository's existing
motor feasibility constraints (including >=40 W indicated power). There is no
explicit size penalty: exchanger size is paid only through the physical model
(dead volume, wall capacity, hydraulics, etc.). Fan power remains excluded from
the objective per the current project decision.

The search is persistent and resumable. It begins with explicit controls:
1. retained K2 champion: 3708 / 3708 tubes;
2. retained levels, 120/60/120/60, H_o=6050;
3. retained levels, 150/30/150/30, H_o=4840;
4. retained levels, 150/30/150/30, current 3708 / 3708 hardware.

It then alternates broad 9D Sobol exploration and incumbent-centered Sobol
proposals.

Example:
    PYTHONPATH=src python3 examples/optimize_motor_four_stage_hx9d.py \
        --budget-seconds 10800 --evaluations 256 --candidate-seconds 120

Outputs:
    outputs/motor_four_stage_hx9d/history.jsonl
    outputs/motor_four_stage_hx9d/report.json
    outputs/motor_four_stage_hx9d/definition.json
"""

from __future__ import annotations

from dataclasses import asdict, replace
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
from optimize_motor_exchanger_asymmetry_stageK1 import (
    _hardware_metrics,
    _scaled_exchanger,
)
from optimize_motor_four_stage_k2 import _symmetry_diagnostics
from optimize_motor_piecewise_stageP3 import _feasibility


THERMO3D_REPORT = ROOT / "outputs" / "motor_four_stage_thermo3d" / "report.json"
DEFAULT_DIRECTORY = ROOT / "outputs" / "motor_four_stage_hx9d"

KIN_NAMES = ("t1", "t2", "t3", "a_l", "b_l", "a_s", "b_s")
PARAMETER_NAMES = (*KIN_NAMES, "n_i", "n_o")

CHARGE_PRESSURE_PA = 100_000.0
MINIMUM_STAGE_FRACTION = 0.02
DEFAULT_HI_COUNT_BOUNDS = (1600, 8000)
DEFAULT_HO_COUNT_BOUNDS = (1600, 8000)
SOBOL_SEED = 26091603


def _load_history(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return [
        json.loads(line)
        for line in path.read_text().splitlines()
        if line.strip()
    ]


def _eligible(record: dict) -> bool:
    return bool(
        record.get("feasible")
        and record.get("result", {}).get("status") == "converged"
        and record.get("result", {}).get("indicated_thermal_efficiency")
        is not None
    )


def _efficiency(record: dict) -> float:
    return float(record["result"]["indicated_thermal_efficiency"])


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


def _valid_parameters(parameters: dict) -> bool:
    times = np.asarray(
        [parameters["t1"], parameters["t2"], parameters["t3"]],
        dtype=float,
    )
    levels = np.asarray(
        [parameters["a_l"], parameters["b_l"],
         parameters["a_s"], parameters["b_s"]],
        dtype=float,
    )
    if np.any(~np.isfinite(times)) or np.any(~np.isfinite(levels)):
        return False
    if np.any(levels < 0.0) or np.any(levels > 1.0):
        return False
    durations = np.diff(np.r_[0.0, times, 1.0])
    return bool(np.min(durations) >= MINIMUM_STAGE_FRACTION)


def _candidate_id(parameters: dict) -> str:
    payload = {
        name: (
            int(parameters[name])
            if name in ("n_i", "n_o")
            else float(parameters[name])
        )
        for name in PARAMETER_NAMES
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _count_from_unit(u: float, bounds: tuple[int, int]) -> int:
    lo, hi = bounds
    value = lo + float(u) * (hi - lo)
    return int(np.clip(round(value), lo, hi))


def _count_to_unit(n: int, bounds: tuple[int, int]) -> float:
    lo, hi = bounds
    return (int(n) - lo) / (hi - lo)


def _ordered_times_from_unit(u: np.ndarray) -> tuple[float, float, float]:
    """Map 3 Sobol coordinates to four ordered stages with a 2% minimum."""
    remaining = 1.0 - 4.0 * MINIMUM_STAGE_FRACTION
    gaps = []
    for j in range(3):
        # Stick-breaking mapping: broad ordered durations without sorting,
        # retaining all four minimum durations.
        portion = remaining * (
            1.0 - (1.0 - float(u[j])) ** (1.0 / (3 - j))
        )
        gaps.append(MINIMUM_STAGE_FRACTION + portion)
        remaining -= portion
    return tuple(np.cumsum(gaps))


def _global_parameters(
    u: np.ndarray,
    hi_bounds: tuple[int, int],
    ho_bounds: tuple[int, int],
) -> dict:
    t1, t2, t3 = _ordered_times_from_unit(u[:3])
    return {
        "t1": float(t1),
        "t2": float(t2),
        "t3": float(t3),
        "a_l": float(u[3]),
        "b_l": float(u[4]),
        "a_s": float(u[5]),
        "b_s": float(u[6]),
        "n_i": _count_from_unit(u[7], hi_bounds),
        "n_o": _count_from_unit(u[8], ho_bounds),
    }


def _normalized_vector(
    parameters: dict,
    hi_bounds: tuple[int, int],
    ho_bounds: tuple[int, int],
) -> np.ndarray:
    return np.asarray(
        [
            parameters["t1"],
            parameters["t2"],
            parameters["t3"],
            parameters["a_l"],
            parameters["b_l"],
            parameters["a_s"],
            parameters["b_s"],
            _count_to_unit(parameters["n_i"], hi_bounds),
            _count_to_unit(parameters["n_o"], ho_bounds),
        ],
        dtype=float,
    )


def _local_parameters(
    center: dict,
    u: np.ndarray,
    radius: float,
    hi_bounds: tuple[int, int],
    ho_bounds: tuple[int, int],
) -> dict:
    vector = _normalized_vector(center, hi_bounds, ho_bounds)
    proposal = vector + (2.0 * u - 1.0) * radius

    # Kinematic coordinates stay literal fractions; exchanger coordinates are
    # normalized only for proposal generation.
    result = {
        name: float(proposal[i])
        for i, name in enumerate(KIN_NAMES)
    }
    result["n_i"] = _count_from_unit(
        float(np.clip(proposal[7], 0.0, 1.0)), hi_bounds
    )
    result["n_o"] = _count_from_unit(
        float(np.clip(proposal[8], 0.0, 1.0)), ho_bounds
    )
    return result


def _geometry_metrics(exchanger) -> dict:
    dimensions = exchanger.bank.dimensions()
    metrics = _hardware_metrics(exchanger)
    return {
        "tube_count": int(exchanger.bank.tube_count),
        "tube_length_m": float(exchanger.bank.tube_length_m),
        "inner_diameter_m": float(exchanger.bank.inner_diameter_m),
        "outlet_valve_cda_m2": float(exchanger.outlet_valve_cda_m2),
        "tube_flow_area_m2": float(dimensions["tube_flow_area_m2"]),
        "core_width_m": float(dimensions["core_width_m"]),
        "core_height_m": float(dimensions["core_height_m"]),
        **{key: float(value) for key, value in metrics.items()},
    }


def _load_control() -> dict:
    if not THERMO3D_REPORT.exists():
        raise FileNotFoundError(THERMO3D_REPORT)
    report = json.loads(THERMO3D_REPORT.read_text())
    control = report.get("control_1bar")
    if control is None or not control.get("feasible"):
        raise RuntimeError("Missing feasible 1-bar thermo-3D K2 control.")
    return control


def _build_basis():
    definition = CampaignDefinition(A5_CAMPAIGN)
    base, _legacy_mass, _saved, _best, _candidate = (
        _candidate_design_and_mass(definition)
    )
    control = _load_control()

    retained_hi = _scaled_exchanger(
        base.heat_in,
        float(control["parameters"]["k_i"]),
    )
    retained_ho = _scaled_exchanger(
        base.heat_out,
        float(control["parameters"]["k_o"]),
    )

    charge = ChargeConfiguration(
        temperature=base.configuration.charge.temperature,
        pressure=CHARGE_PRESSURE_PA,
    )
    base = replace(
        base,
        configuration=replace(base.configuration, charge=charge),
        heat_in=retained_hi,
        heat_out=retained_ho,
    )
    return definition, base, control


def _build_design(base, parameters: dict):
    retained_hi_count = int(base.heat_in.bank.tube_count)
    retained_ho_count = int(base.heat_out.bank.tube_count)

    hi_ratio = int(parameters["n_i"]) / retained_hi_count
    ho_ratio = int(parameters["n_o"]) / retained_ho_count

    heat_in = replace(
        base.heat_in,
        bank=replace(
            base.heat_in.bank,
            tube_count=int(parameters["n_i"]),
        ),
        outlet_valve_cda_m2=base.heat_in.outlet_valve_cda_m2 * hi_ratio,
    )
    heat_out = replace(
        base.heat_out,
        bank=replace(
            base.heat_out.bank,
            tube_count=int(parameters["n_o"]),
        ),
        outlet_valve_cda_m2=base.heat_out.outlet_valve_cda_m2 * ho_ratio,
    )

    limits = base.configuration.machine_volumes
    kin_mapping = {
        name: float(parameters[name])
        for name in KIN_NAMES
    }
    kinematics = FourStageVolumeKinematics(
        limits.small_cylinder,
        limits.large_cylinder,
        **kin_mapping,
    )

    return replace(
        base,
        heat_in=heat_in,
        heat_out=heat_out,
        kinematics=kinematics,
    )


def _target_uniform_state(design) -> tuple[np.ndarray, dict, dict]:
    wrapper = design.build()
    initial = _uniform_wall_initial(design.configuration, wrapper)
    hi = _geometry_metrics(design.heat_in)
    ho = _geometry_metrics(design.heat_out)
    return np.asarray(initial, dtype=float), hi, ho


def _warm_start_from_record(
    source: dict,
    target_uniform: np.ndarray,
    target_hi: dict,
    target_ho: dict,
) -> np.ndarray:
    """Rescale a periodic state to the new 1-bar candidate inventory.

    Gas masses AND internal energies are scaled by the same factor, preserving
    specific energy and temperature. Wall energies are scaled by wall capacity,
    preserving wall temperature.
    """
    state = np.asarray(source["last_complete_state"], dtype=float).copy()
    if state.shape != (10,):
        raise ValueError("Expected a 10-state dynamic-wall warm start.")

    source_mass = float(np.sum(state[:8:2]))
    target_mass = float(np.sum(target_uniform[:8:2]))
    if source_mass <= 0.0 or target_mass <= 0.0:
        raise ValueError("Warm-start gas inventory must be positive.")

    state[:8] *= target_mass / source_mass

    source_hi_capacity = float(source["hardware"]["H_i"]["wall_capacity_j_k"])
    source_ho_capacity = float(source["hardware"]["H_o"]["wall_capacity_j_k"])
    state[8] *= float(target_hi["wall_capacity_j_k"]) / source_hi_capacity
    state[9] *= float(target_ho["wall_capacity_j_k"]) / source_ho_capacity
    return state


def _seed_parameters(control: dict) -> list[tuple[str, dict]]:
    retained = control["motion_parameters"]

    def motion(t1, t2, t3):
        return {
            "t1": float(t1),
            "t2": float(t2),
            "t3": float(t3),
            "a_l": float(retained["a_l"]),
            "b_l": float(retained["b_l"]),
            "a_s": float(retained["a_s"]),
            "b_s": float(retained["b_s"]),
        }

    champion = {
        **{name: float(retained[name]) for name in KIN_NAMES},
        "n_i": 3708,
        "n_o": 3708,
    }
    sizing_120 = {
        **motion(120 / 360, 180 / 360, 300 / 360),
        "n_i": 3708,
        "n_o": 6050,
    }
    sizing_150_reduced = {
        **motion(150 / 360, 180 / 360, 330 / 360),
        "n_i": 3708,
        "n_o": 4840,
    }
    timing_150_current = {
        **motion(150 / 360, 180 / 360, 330 / 360),
        "n_i": 3708,
        "n_o": 3708,
    }

    return [
        ("retained_K2_champion_control", champion),
        ("120_60_120_60_Ho6050_control", sizing_120),
        ("150_30_150_30_Ho4840_control", sizing_150_reduced),
        ("150_30_150_30_current_hardware_control", timing_150_current),
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget-seconds", type=float, default=10800.0)
    parser.add_argument("--evaluations", type=int, default=256)
    parser.add_argument("--candidate-seconds", type=float, default=120.0)
    parser.add_argument(
        "--hi-count-min", type=int, default=DEFAULT_HI_COUNT_BOUNDS[0]
    )
    parser.add_argument(
        "--hi-count-max", type=int, default=DEFAULT_HI_COUNT_BOUNDS[1]
    )
    parser.add_argument(
        "--ho-count-min", type=int, default=DEFAULT_HO_COUNT_BOUNDS[0]
    )
    parser.add_argument(
        "--ho-count-max", type=int, default=DEFAULT_HO_COUNT_BOUNDS[1]
    )
    parser.add_argument("--seed", type=int, default=SOBOL_SEED)
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=DEFAULT_DIRECTORY,
    )
    args = parser.parse_args()

    if args.budget_seconds <= 0:
        raise ValueError("--budget-seconds must be positive.")
    if args.evaluations <= 0:
        raise ValueError("--evaluations must be positive.")
    if args.candidate_seconds <= 0:
        raise ValueError("--candidate-seconds must be positive.")

    hi_bounds = (args.hi_count_min, args.hi_count_max)
    ho_bounds = (args.ho_count_min, args.ho_count_max)
    if hi_bounds[0] < 1 or hi_bounds[0] >= hi_bounds[1]:
        raise ValueError("Invalid H_i tube-count bounds.")
    if ho_bounds[0] < 1 or ho_bounds[0] >= ho_bounds[1]:
        raise ValueError("Invalid H_o tube-count bounds.")

    definition, base, control = _build_basis()
    seeds = _seed_parameters(control)

    directory = args.output_directory
    directory.mkdir(parents=True, exist_ok=True)
    history_path = directory / "history.jsonl"
    definition_path = directory / "definition.json"

    identity = {
        "parameter_names": list(PARAMETER_NAMES),
        "minimum_stage_fraction": MINIMUM_STAGE_FRACTION,
        "H_i_tube_count_bounds": list(hi_bounds),
        "H_o_tube_count_bounds": list(ho_bounds),
        "charge_pressure_pa": CHARGE_PRESSURE_PA,
        "frequency_hz": abs(base.configuration.angular_speed) / (2 * math.pi),
        "retained_H_i_tube_length_m": base.heat_in.bank.tube_length_m,
        "retained_H_o_tube_length_m": base.heat_out.bank.tube_length_m,
        "tube_inner_diameter_m": base.heat_in.bank.inner_diameter_m,
        "retained_H_i_tube_count": base.heat_in.bank.tube_count,
        "retained_H_o_tube_count": base.heat_out.bank.tube_count,
        "sobol_seed": args.seed,
        "objective": (
            "Maximize indicated thermal efficiency among candidates satisfying "
            "the existing motor feasibility constraints, including >=40 W."
        ),
        "mass_policy": (
            "Candidate-specific uniform 100 kPa filling at theta=0 using actual "
            "cylinder and exchanger volumes."
        ),
        "valve_policy": (
            "Each exchanger outlet-valve CdA scales linearly with tube count "
            "relative to retained K2 hardware."
        ),
        "search_strategy": (
            "Four explicit controls, then alternating broad 9D Sobol proposals "
            "and incumbent-centered Sobol proposals. Local radius decreases "
            "0.12 -> 0.06 -> 0.03 in normalized coordinates."
        ),
        "fan_power_policy": (
            "Reported by exchanger metadata but excluded from the efficiency "
            "objective under the current project decision."
        ),
        "seeds": [
            {"kind": kind, "parameters": params}
            for kind, params in seeds
        ],
        "thermo3d_source_sha256": hashlib.sha256(
            THERMO3D_REPORT.read_bytes()
        ).hexdigest(),
    }

    if definition_path.exists():
        saved = json.loads(definition_path.read_text())
        if saved != identity:
            raise ValueError(
                "9D optimization definition changed; use a new output "
                "directory or restore the previous settings."
            )
    else:
        definition_path.write_text(json.dumps(identity, indent=2) + "\n")

    history = _load_history(history_path)
    existing_finished = {
        record["candidate_id"]
        for record in history
        if record.get("result", {}).get("status") != "interrupted"
    }
    eligible = [record for record in history if _eligible(record)]
    best = max(eligible, key=_efficiency) if eligible else None
    phase_start_best = best

    sobol = qmc.Sobol(d=9, scramble=True, seed=args.seed)
    points = sobol.random_base2(14)  # 16384 deterministic proposals

    next_proposal = max(
        (record.get("proposal_index", -1) for record in history),
        default=-1,
    ) + 1

    start = time.monotonic()
    deadline = start + args.budget_seconds
    completed_this_run = 0

    def progress(_: object) -> None:
        if time.monotonic() >= deadline:
            raise IntegrationInterrupted(
                "9D optimization wall-clock budget exhausted."
            )
        if time.monotonic() >= candidate_deadline:
            raise IntegrationInterrupted(
                "Candidate wall-clock budget exhausted; result remains unknown."
            )

    while (
        completed_this_run < args.evaluations
        and time.monotonic() < deadline
    ):
        proposal_index = next_proposal
        next_proposal += 1

        if proposal_index < len(seeds):
            kind, parameters = seeds[proposal_index]
            parameters = dict(parameters)
            radius = None
        else:
            sobol_index = proposal_index - len(seeds)
            if sobol_index >= len(points):
                raise RuntimeError("Exhausted deterministic Sobol proposal bank.")
            u = points[sobol_index]

            # Keep one broad proposal in four throughout the campaign.
            if best is None or sobol_index % 4 == 0:
                parameters = _global_parameters(u, hi_bounds, ho_bounds)
                kind = "global_sobol_9d"
                radius = None
            else:
                completed_after_seeds = max(0, len(history) - len(seeds))
                if completed_after_seeds < 64:
                    radius = 0.12
                elif completed_after_seeds < 128:
                    radius = 0.06
                else:
                    radius = 0.03
                parameters = _local_parameters(
                    best["parameters"],
                    u,
                    radius,
                    hi_bounds,
                    ho_bounds,
                )
                kind = "incumbent_neighborhood_sobol_9d"

        if not _valid_parameters(parameters):
            continue

        candidate_id = _candidate_id(parameters)
        if candidate_id in existing_finished:
            continue

        design = _build_design(base, parameters)
        target_uniform, hi_metrics, ho_metrics = _target_uniform_state(design)

        # Warm start from the nearest converged 9D candidate, with inventory and
        # wall-energy rescaling. Explicit control seeds are evaluated from their
        # own uniform filling state so they remain clean controls.
        reusable = [
            record
            for record in history
            if record.get("result", {}).get("status") == "converged"
            and record.get("last_complete_state") is not None
            and record.get("hardware")
        ]
        source = None
        initial_state = target_uniform
        if proposal_index >= len(seeds) and reusable:
            target_vector = _normalized_vector(
                parameters, hi_bounds, ho_bounds
            )
            source = min(
                reusable,
                key=lambda record: np.linalg.norm(
                    target_vector
                    - _normalized_vector(
                        record["parameters"], hi_bounds, ho_bounds
                    )
                ),
            )
            initial_state = _warm_start_from_record(
                source,
                target_uniform,
                hi_metrics,
                ho_metrics,
            )

        before = time.monotonic()
        candidate_deadline = min(
            deadline,
            before + args.candidate_seconds,
        )
        try:
            result, state = _evaluate(
                f"four_stage_hx9d_{proposal_index}",
                design,
                definition,
                initial_state=np.asarray(initial_state, dtype=float),
                progress_callback=progress,
            )
        except IntegrationInterrupted as exc:
            result = {
                "status": "interrupted",
                "message": str(exc),
            }
            state = None
        except (ValueError, RuntimeError) as exc:
            result = {
                "status": "integration_failure",
                "message": str(exc),
            }
            state = None

        feasible, reasons = (
            _feasibility(result)
            if result.get("status") == "converged"
            else (False, [])
        )

        total_mass = (
            result.get("total_mass_kg")
            if result.get("status") == "converged"
            else float(np.sum(target_uniform[:8:2]))
        )

        record = {
            "index": len(history),
            "proposal_index": proposal_index,
            "candidate_id": candidate_id,
            "kind": kind,
            "radius": radius,
            "parameters": {
                **{name: float(parameters[name]) for name in KIN_NAMES},
                "n_i": int(parameters["n_i"]),
                "n_o": int(parameters["n_o"]),
            },
            "phase_degrees": _phase_degrees(parameters),
            "symmetry": _symmetry_diagnostics(parameters),
            "hardware": {
                "H_i": hi_metrics,
                "H_o": ho_metrics,
            },
            "candidate_fill_mass_kg": total_mass,
            "result": result,
            "feasible": feasible,
            "physical_constraint_failures": reasons,
            "elapsed_seconds": time.monotonic() - before,
            "warm_start_source": (
                source["candidate_id"] if source is not None else None
            ),
            "last_complete_state": (
                state.tolist() if state is not None else None
            ),
        }

        with history_path.open("a") as stream:
            stream.write(json.dumps(record) + "\n")
            stream.flush()

        history.append(record)
        completed_this_run += 1
        if result.get("status") != "interrupted":
            existing_finished.add(candidate_id)

        if feasible and (
            best is None or _efficiency(record) > _efficiency(best)
        ):
            best = record

        seed_records = {
            record["kind"]: record
            for record in history
            if record["kind"].endswith("_control")
        }

        report = {
            "best_feasible": best,
            "phase_start_best": phase_start_best,
            "controls": seed_records,
            "attempted_total": len(history),
            "completed_this_run": completed_this_run,
            "requested_duration_seconds": args.budget_seconds,
            "actual_duration_seconds": time.monotonic() - start,
            "definition": identity,
        }
        (directory / "report.json").write_text(
            json.dumps(report, indent=2) + "\n"
        )

        print(
            json.dumps(
                {
                    "index": record["index"],
                    "proposal_index": proposal_index,
                    "kind": kind,
                    "status": result.get("status"),
                    "feasible": feasible,
                    "efficiency": result.get(
                        "indicated_thermal_efficiency"
                    ),
                    "power_W": result.get("indicated_power_w"),
                    "Qin_W": result.get("heat_input_w"),
                    "n_i": int(parameters["n_i"]),
                    "n_o": int(parameters["n_o"]),
                    "phases_deg": record["phase_degrees"],
                    "best_efficiency": (
                        _efficiency(best) if best is not None else None
                    ),
                    "elapsed_seconds": record["elapsed_seconds"],
                }
            ),
            flush=True,
        )

    print(f"Saved {directory / 'report.json'}", flush=True)


if __name__ == "__main__":
    main()
