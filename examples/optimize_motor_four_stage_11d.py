"""Persistent 11-variable timing-release search on the final K2 / 1-bar motor.

The final seven-variable four-stage law is embedded exactly as the seed, then
the two cylinder chronologies are released:

    t0_s
    t1_l, t2_l, t3_l
    t1_s, t2_s, t3_s
    a_l, b_l, a_s, b_s

The large-cylinder origin remains the global phase gauge. t0_s is the signed
cyclic origin of the small-cylinder law in [-0.5, 0.5).

Frozen physical basis:
- air;
- 100 kPa absolute filling pressure at global t=0;
- 2 Hz;
- Stage K2 exchanger multipliers;
- S/L swept-volume ratio 0.84 and the current total swept volume;
- one-percent clearance ratios;
- no mechanical-loss model.

Because t0_s changes V_S at the filling origin, total gas mass is recalculated
for every candidate so the configured filling pressure remains exactly 100 kPa.
Warm starts rescale gas masses and gas energies together, preserving specific
internal energies; wall energies are unchanged because exchanger hardware is
frozen.

All eleven coordinates vary simultaneously. The default local radius schedule
is 0.10 -> 0.05 -> 0.025, with a diverse elite pool retained.
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
from dada_solver.independent_four_stage_kinematics import (
    IndependentFourStageVolumeKinematics,
)
from dada_solver.integration import IntegrationInterrupted

from compare_motor_motion_laws_stage7A5 import (
    A5_CAMPAIGN,
    ROOT,
    _candidate_design_and_mass,
    _evaluate,
)
from optimize_motor_exchanger_asymmetry_stageK1 import _scaled_exchanger
from optimize_motor_piecewise_stageP3 import _feasibility


THERMO3D_REPORT = ROOT / "outputs" / "motor_four_stage_thermo3d" / "report.json"
DEFAULT_DIRECTORY = ROOT / "outputs" / "motor_four_stage_11d"

NAMES = (
    "t0_s",
    "t1_l",
    "t2_l",
    "t3_l",
    "t1_s",
    "t2_s",
    "t3_s",
    "a_l",
    "b_l",
    "a_s",
    "b_s",
)

MINIMUM_STAGE_FRACTION = 0.02
DEFAULT_RADII = (0.10, 0.05, 0.025)
SOBOL_SEED = 260918
CHARGE_PRESSURE_PA = 100_000.0


def _load_history(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return [
        json.loads(line)
        for line in path.read_text().splitlines()
        if line.strip()
    ]


def _parse_radii(text: str) -> tuple[float, ...]:
    values = tuple(float(item) for item in text.split(",") if item.strip())
    if not values or any(value <= 0.0 for value in values):
        raise argparse.ArgumentTypeError(
            "radii must be positive comma-separated values"
        )
    if any(second >= first for first, second in zip(values, values[1:])):
        raise argparse.ArgumentTypeError("radii must be strictly decreasing")
    return values


def _vector(mapping: dict[str, float]) -> np.ndarray:
    return np.asarray([mapping[name] for name in NAMES], dtype=float)


def _canonical_t0(value: float) -> float:
    return float(((value + 0.5) % 1.0) - 0.5)


def _distance(left: np.ndarray, right: np.ndarray) -> float:
    delta = np.asarray(left - right, dtype=float)
    delta[0] = ((delta[0] + 0.5) % 1.0) - 0.5
    return float(np.linalg.norm(delta))


def _efficiency(record: dict) -> float:
    return float(record["result"]["indicated_thermal_efficiency"])


def _eligible(record: dict) -> bool:
    return bool(
        record.get("feasible")
        and record.get("result", {}).get("status") == "converged"
        and record.get("result", {}).get("indicated_thermal_efficiency")
        is not None
    )


def _valid_parameters(vector: np.ndarray) -> bool:
    if not -0.5 <= vector[0] < 0.5:
        return False

    large_durations = np.diff(np.r_[0.0, vector[1:4], 1.0])
    small_durations = np.diff(np.r_[0.0, vector[4:7], 1.0])

    if min(np.min(large_durations), np.min(small_durations)) < MINIMUM_STAGE_FRACTION:
        return False

    return bool(np.all(vector[7:] >= 0.0) and np.all(vector[7:] <= 1.0))


def _candidate_id(mapping: dict[str, float]) -> str:
    return hashlib.sha256(
        json.dumps(mapping, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _select_diverse_elite(records, count, minimum_distance):
    ranked = sorted(
        (record for record in records if _eligible(record)),
        key=_efficiency,
        reverse=True,
    )
    selected = []

    for record in ranked:
        candidate = _vector(record["parameters"])
        if all(
            _distance(candidate, _vector(other["parameters"])) >= minimum_distance
            for other in selected
        ):
            selected.append(record)
            if len(selected) >= count:
                break

    if len(selected) < count:
        selected_ids = {record["candidate_id"] for record in selected}
        for record in ranked:
            if record["candidate_id"] not in selected_ids:
                selected.append(record)
                selected_ids.add(record["candidate_id"])
                if len(selected) >= count:
                    break

    return selected


def _timing_diagnostics(parameters):
    large = np.diff(
        np.asarray(
            [0.0, parameters["t1_l"], parameters["t2_l"], parameters["t3_l"], 1.0]
        )
    )
    small = np.diff(
        np.asarray(
            [0.0, parameters["t1_s"], parameters["t2_s"], parameters["t3_s"], 1.0]
        )
    )

    t0 = parameters["t0_s"]
    absolute_small = np.mod(
        np.asarray(
            [
                t0,
                t0 + parameters["t1_s"],
                t0 + parameters["t2_s"],
                t0 + parameters["t3_s"],
            ]
        ),
        1.0,
    )

    delta = small - large
    return {
        "t0_s": t0,
        "t0_s_degrees": 360.0 * t0,
        "large_stage_fractions": large.tolist(),
        "small_stage_fractions": small.tolist(),
        "large_stage_degrees": (360.0 * large).tolist(),
        "small_stage_degrees": (360.0 * small).tolist(),
        "small_absolute_event_fractions": absolute_small.tolist(),
        "small_absolute_event_degrees": (360.0 * absolute_small).tolist(),
        "stage_duration_delta_s_minus_l": delta.tolist(),
        "stage_duration_rms_difference": float(
            math.sqrt(float(np.mean(delta * delta)))
        ),
    }


def _level_symmetry(parameters):
    delta_a = float(parameters["a_l"] - parameters["a_s"])
    delta_b = float(parameters["b_l"] - parameters["b_s"])
    return {
        "delta_a": delta_a,
        "delta_b": delta_b,
        "rms": math.sqrt((delta_a * delta_a + delta_b * delta_b) / 2.0),
    }


def _target_filling_mass(design):
    """Return filling mass and connected gas volume at theta=0."""
    wrapper = design.build()
    model = getattr(wrapper, "model", wrapper)
    volumes = model.volumes(0.0)
    config = design.configuration
    mass = (
        CHARGE_PRESSURE_PA
        * volumes.total
        / (config.gas.gas_constant * config.charge.temperature)
    )
    return float(mass), float(volumes.total)


def _rescale_gas_inventory(state, target_mass):
    """Preserve gas specific energies while changing total inventory."""
    values = np.asarray(state, dtype=float).copy()
    if values.ndim != 1 or len(values) < 8:
        raise ValueError("Warm-start state has an unexpected layout.")

    source_mass = float(np.sum(values[:8:2]))
    if not math.isfinite(source_mass) or source_mass <= 0.0:
        raise ValueError("Warm-start state has invalid gas inventory.")

    values[:8] *= target_mass / source_mass
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget-seconds", type=float, default=10800.0)
    parser.add_argument("--evaluations", type=int, default=256)
    parser.add_argument(
        "--candidate-seconds",
        type=float,
        default=120.0,
        help="Cooperative per-candidate wall-clock limit.",
    )
    parser.add_argument(
        "--radii",
        type=_parse_radii,
        default=DEFAULT_RADII,
        help="Local radius schedule; default: 0.10,0.05,0.025.",
    )
    parser.add_argument(
        "--evaluations-per-radius",
        type=int,
        default=96,
        help="Completed local evaluations before reducing the radius.",
    )
    parser.add_argument(
        "--elite-centers",
        type=int,
        default=4,
        help="Number of distinct feasible centers retained.",
    )
    parser.add_argument(
        "--elite-distance",
        type=float,
        default=0.04,
        help="Minimum cyclic-aware 11D distance between elite centers.",
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=DEFAULT_DIRECTORY,
    )
    args = parser.parse_args()

    if args.budget_seconds <= 0.0:
        raise ValueError("--budget-seconds must be positive.")
    if args.evaluations <= 0:
        raise ValueError("--evaluations must be positive.")
    if args.candidate_seconds <= 0.0:
        raise ValueError("--candidate-seconds must be positive.")
    if args.evaluations_per_radius <= 0:
        raise ValueError("--evaluations-per-radius must be positive.")
    if args.elite_centers <= 0:
        raise ValueError("--elite-centers must be positive.")
    if args.elite_distance < 0.0:
        raise ValueError("--elite-distance must be non-negative.")

    if not THERMO3D_REPORT.exists():
        raise FileNotFoundError(
            f"Missing thermo-3D report: {THERMO3D_REPORT.relative_to(ROOT)}"
        )

    thermo = json.loads(THERMO3D_REPORT.read_text())
    control = thermo["control_1bar"]
    control_parameters = control["parameters"]
    motion = control["motion_parameters"]

    if abs(float(control_parameters["swept_ratio"]) - 0.84) > 1e-12:
        raise ValueError("Expected the retained K2 control at swept ratio 0.84.")
    if abs(float(control["charge_pressure_pa"]) - CHARGE_PRESSURE_PA) > 1e-9:
        raise ValueError("Expected the retained K2 control at 100 kPa.")

    seed_mapping = {
        "t0_s": 0.0,
        "t1_l": float(motion["t1"]),
        "t2_l": float(motion["t2"]),
        "t3_l": float(motion["t3"]),
        "t1_s": float(motion["t1"]),
        "t2_s": float(motion["t2"]),
        "t3_s": float(motion["t3"]),
        "a_l": float(motion["a_l"]),
        "b_l": float(motion["b_l"]),
        "a_s": float(motion["a_s"]),
        "b_s": float(motion["b_s"]),
    }
    seed_vector = _vector(seed_mapping)

    definition = CampaignDefinition(A5_CAMPAIGN)
    base, _legacy_mass, _saved, _a5, _candidate = (
        _candidate_design_and_mass(definition)
    )

    charge = ChargeConfiguration(
        temperature=base.configuration.charge.temperature,
        pressure=CHARGE_PRESSURE_PA,
    )
    base = replace(
        base,
        configuration=replace(base.configuration, charge=charge),
        heat_in=_scaled_exchanger(
            base.heat_in,
            float(control_parameters["k_i"]),
        ),
        heat_out=_scaled_exchanger(
            base.heat_out,
            float(control_parameters["k_o"]),
        ),
    )
    limits = base.configuration.machine_volumes

    directory = args.output_directory
    directory.mkdir(parents=True, exist_ok=True)
    history_path = directory / "history.jsonl"
    definition_path = directory / "definition.json"

    identity = {
        "names": list(NAMES),
        "seed": seed_mapping,
        "minimum_stage_fraction": MINIMUM_STAGE_FRACTION,
        "radii": list(args.radii),
        "evaluations_per_radius": args.evaluations_per_radius,
        "elite_centers": args.elite_centers,
        "elite_distance": args.elite_distance,
        "sobol_seed": SOBOL_SEED,
        "charge_pressure_pa": CHARGE_PRESSURE_PA,
        "frequency_hz": abs(base.configuration.angular_speed) / (2.0 * math.pi),
        "K2_k_i": float(control_parameters["k_i"]),
        "K2_k_o": float(control_parameters["k_o"]),
        "swept_ratio": float(control_parameters["swept_ratio"]),
        "source_thermo3d_sha256": hashlib.sha256(
            THERMO3D_REPORT.read_bytes()
        ).hexdigest(),
        "strategy": (
            "Exact shared-7D seed followed by cyclic-aware local scrambled "
            "Sobol refinement in all eleven coordinates simultaneously; "
            "diverse elite centers retained; no symmetry penalty."
        ),
        "mass_policy": (
            "Each candidate is filled at 100 kPa at global t=0. Warm-start "
            "gas masses and gas energies are rescaled together to the "
            "candidate filling inventory; wall energies remain unchanged."
        ),
    }

    if definition_path.exists():
        if json.loads(definition_path.read_text()) != identity:
            raise ValueError(
                "11D definition changed; use a new output directory or "
                "restore the original search settings."
            )
    else:
        definition_path.write_text(json.dumps(identity, indent=2) + "\n")

    history = _load_history(history_path)
    existing_ids = {
        record["candidate_id"]
        for record in history
        if record.get("result", {}).get("status") != "interrupted"
    }

    eligible = [record for record in history if _eligible(record)]
    best = max(eligible, key=_efficiency) if eligible else None
    phase_start_best = best

    sampler = qmc.Sobol(d=len(NAMES), scramble=True, seed=SOBOL_SEED)
    points = sampler.random_base2(14)

    next_proposal_index = max(
        (record.get("proposal_index", -1) for record in history),
        default=-1,
    ) + 1

    start = time.monotonic()
    deadline = start + args.budget_seconds
    completed_this_run = 0

    def progress(_: object) -> None:
        if time.monotonic() >= deadline:
            raise IntegrationInterrupted("11D phase wall-clock budget exhausted.")
        if time.monotonic() >= candidate_deadline:
            raise IntegrationInterrupted(
                "Candidate wall-clock budget exhausted; result remains unknown."
            )

    while (
        completed_this_run < args.evaluations
        and time.monotonic() < deadline
    ):
        proposal_index = next_proposal_index
        next_proposal_index += 1

        if proposal_index == 0 and not history:
            vector = seed_vector.copy()
            kind = "shared_7d_seed"
            center = None
            radius = 0.0
            radius_index = -1
        else:
            local_count = sum(
                record.get("kind") == "local_11d_sobol"
                for record in history
            )
            radius_index = min(
                local_count // args.evaluations_per_radius,
                len(args.radii) - 1,
            )
            radius = args.radii[radius_index]

            elite = _select_diverse_elite(
                history,
                args.elite_centers,
                args.elite_distance,
            )
            if not elite:
                raise RuntimeError(
                    "The exact 7D seed did not produce a feasible 11D baseline."
                )

            point_index = proposal_index - 1
            if point_index >= len(points):
                raise RuntimeError(
                    "11D search exhausted its 16384 deterministic Sobol proposals."
                )

            center = elite[point_index % len(elite)]
            center_vector = _vector(center["parameters"])
            vector = center_vector + (2.0 * points[point_index] - 1.0) * radius
            vector[0] = _canonical_t0(vector[0])
            kind = "local_11d_sobol"

        if not _valid_parameters(vector):
            continue

        mapping = {
            name: float(value)
            for name, value in zip(NAMES, vector, strict=True)
        }
        candidate_id = _candidate_id(mapping)

        if candidate_id in existing_ids:
            continue

        kinematics = IndependentFourStageVolumeKinematics(
            limits.small_cylinder,
            limits.large_cylinder,
            **mapping,
        )
        design = replace(base, kinematics=kinematics)

        target_mass, filling_volume = _target_filling_mass(design)

        reusable = [
            record
            for record in history
            if record.get("result", {}).get("status") == "converged"
            and record.get("last_complete_state")
        ]
        source = (
            min(
                reusable,
                key=lambda record: _distance(
                    vector,
                    _vector(record["parameters"]),
                ),
            )
            if reusable
            else None
        )

        if source is not None:
            initial_state = _rescale_gas_inventory(
                source["last_complete_state"],
                target_mass,
            )
            warm_start_source = source["candidate_id"]
        else:
            initial_state = _rescale_gas_inventory(
                control["last_complete_state"],
                target_mass,
            )
            warm_start_source = "thermo3d_K2_1bar_control"

        before = time.monotonic()
        candidate_deadline = min(
            deadline,
            before + args.candidate_seconds,
        )

        try:
            result, state = _evaluate(
                f"four_stage_11d_{len(history)}",
                design,
                definition,
                initial_state=initial_state,
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

        record = {
            "index": len(history),
            "proposal_index": proposal_index,
            "candidate_id": candidate_id,
            "kind": kind,
            "radius": radius,
            "radius_index": radius_index,
            "center_candidate_id": (
                center["candidate_id"] if center is not None else None
            ),
            "center_efficiency": (
                _efficiency(center) if center is not None else None
            ),
            "parameters": mapping,
            "timing_diagnostics": _timing_diagnostics(mapping),
            "level_symmetry": _level_symmetry(mapping),
            "filling_pressure_pa": CHARGE_PRESSURE_PA,
            "filling_total_volume_m3": filling_volume,
            "target_filling_mass_kg": target_mass,
            "result": result,
            "feasible": feasible,
            "physical_constraint_failures": reasons,
            "elapsed_seconds": time.monotonic() - before,
            "warm_start_source": warm_start_source,
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
            existing_ids.add(candidate_id)

        if feasible and (
            best is None or _efficiency(record) > _efficiency(best)
        ):
            best = record

        if kind == "shared_7d_seed" and result.get("status") == "converged":
            delta = (
                float(result["indicated_thermal_efficiency"])
                - float(control["result"]["indicated_thermal_efficiency"])
            )
            if abs(delta) > 1e-4:
                raise RuntimeError(
                    "The 11D shared seed does not reproduce the retained 7D "
                    f"control closely enough: efficiency delta={delta:+.6g}."
                )

        seed_record = next(
            (
                candidate
                for candidate in history
                if candidate.get("kind") == "shared_7d_seed"
                and _eligible(candidate)
            ),
            None,
        )

        report = {
            "best_feasible": best,
            "shared_7d_seed": seed_record,
            "external_7d_control": {
                "efficiency": control["result"][
                    "indicated_thermal_efficiency"
                ],
                "power_W": control["result"]["indicated_power_w"],
                "mass_kg": control["result"]["total_mass_kg"],
                "parameters": motion,
            },
            "best_timing_diagnostics": (
                _timing_diagnostics(best["parameters"])
                if best is not None
                else None
            ),
            "best_level_symmetry": (
                _level_symmetry(best["parameters"])
                if best is not None
                else None
            ),
            "phase_start_best": phase_start_best,
            "attempted_total": len(history),
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
                    "radius": radius,
                    "status": result.get("status"),
                    "feasible": feasible,
                    "efficiency": result.get(
                        "indicated_thermal_efficiency"
                    ),
                    "power_W": result.get("indicated_power_w"),
                    "filling_mass_kg": target_mass,
                    "t0_s_deg": 360.0 * mapping["t0_s"],
                    "timing_rms": record["timing_diagnostics"][
                        "stage_duration_rms_difference"
                    ],
                    "level_symmetry_rms": record["level_symmetry"]["rms"],
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
