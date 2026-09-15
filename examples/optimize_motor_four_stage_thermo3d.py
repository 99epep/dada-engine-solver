"""Thermodynamic retuning of the final seven-variable four-stage motor law.

The final four-stage motion shape is frozen.  This campaign varies only:

    k_i          H_i tube-length multiplier
    k_o          H_o tube-length multiplier
    swept_ratio  V_S,swept / V_L,swept

The total cylinder swept volume is held constant while swept_ratio changes.
Each cylinder keeps its reference clearance ratio.  The working gas is air and
the filling condition is fixed at 100 kPa absolute at theta=0 for every
candidate; total gas mass is therefore derived from the candidate's actual
cylinder and exchanger gas volumes and is NOT held fixed.

This is intentionally a separate persistent campaign.  It does not modify the
completed seven-variable kinematic histories.

Default broad bounds:
    k_i          [0.75, 1.50]
    k_o          [0.40, 0.85]
    swept_ratio  [0.50, 1.30]

The first recorded candidate is a 1-bar control using the old K2 exchanger
multipliers and the reference cylinder swept-volume ratio.  Subsequent points
use a deterministic scrambled Sobol sequence.

Run, for example:

    PYTHONPATH=src python3 examples/optimize_motor_four_stage_thermo3d.py \
        --budget-seconds 7200 --evaluations 256 --candidate-seconds 120

Rerunning with the same arguments resumes the append-only history.
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
from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.integration import IntegrationInterrupted

from compare_motor_motion_laws_stage7A5 import (
    A5_CAMPAIGN,
    ROOT,
    _candidate_design_and_mass,
    _evaluate,
)
from optimize_motor_exchanger_asymmetry_stageK1 import (
    _hardware_metrics,
    _scaled_exchanger,
)
from optimize_motor_piecewise_stageP3 import _feasibility


FINAL_MOTION_REPORT = (
    ROOT / "outputs" / "motor_four_stage_k2_refine_final" / "report.json"
)
K2_REPORT = ROOT / "outputs" / "motor_exchanger_asymmetry_stageK2.json"
DEFAULT_DIRECTORY = ROOT / "outputs" / "motor_four_stage_thermo3d"

PARAMETER_NAMES = ("k_i", "k_o", "swept_ratio")
BOUNDS = {
    "k_i": (0.75, 1.50),
    "k_o": (0.40, 0.85),
    "swept_ratio": (0.50, 1.30),
}
CHARGE_PRESSURE_PA = 100_000.0
SOBOL_SEED = 260917


def _load_history(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return [
        json.loads(line)
        for line in path.read_text().splitlines()
        if line.strip()
    ]


def _scale(unit: np.ndarray) -> dict[str, float]:
    return {
        name: lo + float(value) * (hi - lo)
        for value, (name, (lo, hi)) in zip(
            unit, BOUNDS.items(), strict=True
        )
    }


def _efficiency(record: dict) -> float:
    return float(record["result"]["indicated_thermal_efficiency"])


def _eligible(record: dict) -> bool:
    return bool(
        record.get("feasible")
        and record.get("result", {}).get("status") == "converged"
        and record.get("result", {}).get("indicated_thermal_efficiency")
        is not None
    )


def _candidate_id(parameters: dict[str, float]) -> str:
    payload = json.dumps(
        {name: float(parameters[name]) for name in PARAMETER_NAMES},
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def _volume_limits(
    total_swept_m3: float,
    swept_ratio: float,
    small_clearance_ratio: float,
    large_clearance_ratio: float,
) -> tuple[CylinderVolumeLimits, CylinderVolumeLimits]:
    """Return S/L limits at fixed total swept volume."""
    if not math.isfinite(swept_ratio) or swept_ratio <= 0.0:
        raise ValueError("Swept-volume ratio must be finite and positive.")

    large_swept = total_swept_m3 / (1.0 + swept_ratio)
    small_swept = total_swept_m3 - large_swept

    small_minimum = small_clearance_ratio * small_swept
    large_minimum = large_clearance_ratio * large_swept

    return (
        CylinderVolumeLimits(
            small_minimum,
            small_minimum + small_swept,
        ),
        CylinderVolumeLimits(
            large_minimum,
            large_minimum + large_swept,
        ),
    )


def _derived_geometry(
    small: CylinderVolumeLimits,
    large: CylinderVolumeLimits,
    hi,
    ho,
) -> dict:
    hi_metrics = _hardware_metrics(hi)
    ho_metrics = _hardware_metrics(ho)
    return {
        "small_swept_volume_m3": small.swept,
        "large_swept_volume_m3": large.swept,
        "total_swept_volume_m3": small.swept + large.swept,
        "small_minimum_volume_m3": small.minimum,
        "large_minimum_volume_m3": large.minimum,
        "small_maximum_volume_m3": small.maximum,
        "large_maximum_volume_m3": large.maximum,
        "small_clearance_ratio": small.minimum / small.swept,
        "large_clearance_ratio": large.minimum / large.swept,
        "H_i": hi_metrics,
        "H_o": ho_metrics,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget-seconds", type=float, default=7200.0)
    parser.add_argument("--evaluations", type=int, default=256)
    parser.add_argument(
        "--candidate-seconds",
        type=float,
        default=120.0,
        help="Cooperative per-candidate wall-clock limit.",
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

    if not FINAL_MOTION_REPORT.exists():
        raise FileNotFoundError(
            f"Missing final motion report: "
            f"{FINAL_MOTION_REPORT.relative_to(ROOT)}"
        )
    if not K2_REPORT.exists():
        raise FileNotFoundError(
            f"Missing K2 report: {K2_REPORT.relative_to(ROOT)}"
        )

    final_motion_report = json.loads(FINAL_MOTION_REPORT.read_text())
    final_motion = final_motion_report["best_overall"]
    motion_parameters = {
        name: float(final_motion["parameters"][name])
        for name in (
            "t1",
            "t2",
            "t3",
            "a_l",
            "b_l",
            "a_s",
            "b_s",
        )
    }

    k2_report = json.loads(K2_REPORT.read_text())
    k2 = k2_report["best_feasible"]

    definition = CampaignDefinition(A5_CAMPAIGN)
    base, _old_mass, _saved, _a5_best, _candidate = (
        _candidate_design_and_mass(definition)
    )

    base_small = base.configuration.machine_volumes.small_cylinder
    base_large = base.configuration.machine_volumes.large_cylinder
    total_swept_m3 = base_small.swept + base_large.swept
    reference_swept_ratio = base_small.swept / base_large.swept
    small_clearance_ratio = base_small.minimum / base_small.swept
    large_clearance_ratio = base_large.minimum / base_large.swept

    control_parameters = {
        "k_i": float(k2["k_i"]),
        "k_o": float(k2["k_o"]),
        "swept_ratio": float(reference_swept_ratio),
    }

    directory = args.output_directory
    directory.mkdir(parents=True, exist_ok=True)
    history_path = directory / "history.jsonl"
    definition_path = directory / "definition.json"

    identity = {
        "parameter_names": list(PARAMETER_NAMES),
        "bounds": {name: list(bounds) for name, bounds in BOUNDS.items()},
        "charge_pressure_pa": CHARGE_PRESSURE_PA,
        "charge_temperature_k": base.configuration.charge.temperature,
        "working_gas": "air_from_A5_configuration",
        "frequency_hz": abs(base.configuration.angular_speed)
        / (2.0 * math.pi),
        "total_swept_volume_m3": total_swept_m3,
        "reference_swept_ratio": reference_swept_ratio,
        "small_clearance_ratio": small_clearance_ratio,
        "large_clearance_ratio": large_clearance_ratio,
        "motion_parameters": motion_parameters,
        "motion_source_sha256": hashlib.sha256(
            FINAL_MOTION_REPORT.read_bytes()
        ).hexdigest(),
        "K2_source_sha256": hashlib.sha256(
            K2_REPORT.read_bytes()
        ).hexdigest(),
        "sobol_seed": SOBOL_SEED,
        "strategy": (
            "One atmospheric-pressure control followed by broad scrambled "
            "3D Sobol search; final seven-variable motion frozen; total "
            "cylinder swept volume fixed; S/L swept-volume ratio, H_i tube "
            "length and H_o tube length varied simultaneously."
        ),
        "mass_policy": (
            "Derived independently for every candidate from a uniform "
            "100 kPa filling state at theta=0 using actual connected "
            "exchanger gas volumes; mass is not an optimization coordinate."
        ),
    }

    if definition_path.exists():
        if json.loads(definition_path.read_text()) != identity:
            raise ValueError(
                "Thermo-3D definition changed; use a new output directory."
            )
    else:
        definition_path.write_text(
            json.dumps(identity, indent=2) + "\n"
        )

    history = _load_history(history_path)
    existing_ids = {record["candidate_id"] for record in history}
    feasible_history = [record for record in history if _eligible(record)]
    best = (
        max(feasible_history, key=_efficiency)
        if feasible_history
        else None
    )
    phase_start_best = best

    # Generate a deterministic oversized Sobol bank.  Evaluation count is a
    # campaign stopping rule, not a requirement that the consumed prefix be a
    # power of two.
    sobol = qmc.Sobol(d=3, scramble=True, seed=SOBOL_SEED)
    points = sobol.random_base2(12)  # 4096 proposals

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
                "Thermo-3D wall-clock budget exhausted."
            )
        if time.monotonic() >= candidate_deadline:
            raise IntegrationInterrupted(
                "Candidate wall-clock budget exhausted."
            )

    while (
        completed_this_run < args.evaluations
        and time.monotonic() < deadline
    ):
        proposal_index = next_proposal
        next_proposal += 1

        if proposal_index == 0:
            parameters = dict(control_parameters)
            kind = "K2_ratio_control_at_1bar"
        else:
            sobol_index = proposal_index - 1
            if sobol_index >= len(points):
                raise RuntimeError(
                    "Thermo-3D campaign exhausted its deterministic "
                    "4096-point Sobol bank."
                )
            parameters = _scale(points[sobol_index])
            kind = "global_sobol_3d"

        candidate_id = _candidate_id(parameters)
        if candidate_id in existing_ids:
            continue

        small_limits, large_limits = _volume_limits(
            total_swept_m3,
            parameters["swept_ratio"],
            small_clearance_ratio,
            large_clearance_ratio,
        )

        volumes = replace(
            base.configuration.machine_volumes,
            small_cylinder=small_limits,
            large_cylinder=large_limits,
        )
        charge = ChargeConfiguration(
            temperature=base.configuration.charge.temperature,
            pressure=CHARGE_PRESSURE_PA,
        )
        config = replace(
            base.configuration,
            machine_volumes=volumes,
            charge=charge,
        )

        kinematics = FourStageVolumeKinematics(
            small_limits,
            large_limits,
            **motion_parameters,
        )
        heat_in = _scaled_exchanger(
            base.heat_in, parameters["k_i"]
        )
        heat_out = _scaled_exchanger(
            base.heat_out, parameters["k_o"]
        )
        design = replace(
            base,
            configuration=config,
            heat_in=heat_in,
            heat_out=heat_out,
            kinematics=kinematics,
        )

        before = time.monotonic()
        candidate_deadline = min(
            deadline, before + args.candidate_seconds
        )
        try:
            result, state = _evaluate(
                f"four_stage_thermo3d_{proposal_index}",
                design,
                definition,
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
            "parameters": parameters,
            "motion_parameters": motion_parameters,
            "derived_geometry": _derived_geometry(
                small_limits,
                large_limits,
                heat_in,
                heat_out,
            ),
            "charge_pressure_pa": CHARGE_PRESSURE_PA,
            "result": result,
            "feasible": feasible,
            "physical_constraint_failures": reasons,
            "elapsed_seconds": time.monotonic() - before,
            "last_complete_state": (
                state.tolist() if state is not None else None
            ),
        }

        with history_path.open("a") as stream:
            stream.write(json.dumps(record) + "\n")
            stream.flush()

        history.append(record)
        existing_ids.add(candidate_id)
        completed_this_run += 1

        if feasible and (
            best is None or _efficiency(record) > _efficiency(best)
        ):
            best = record

        control = next(
            (
                candidate
                for candidate in history
                if candidate.get("kind")
                == "K2_ratio_control_at_1bar"
                and _eligible(candidate)
            ),
            None,
        )

        report = {
            "best_feasible": best,
            "control_1bar": control,
            "phase_start_best": phase_start_best,
            "attempted_total": len(history),
            "requested_duration_seconds": args.budget_seconds,
            "actual_duration_seconds": time.monotonic() - start,
            "final_7D_fixed_motion": {
                "source_candidate_id": final_motion["candidate_id"],
                "source_efficiency_with_old_fixed_mass": final_motion[
                    "result"
                ]["indicated_thermal_efficiency"],
                "parameters": motion_parameters,
            },
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
                    "mass_kg": result.get("total_mass_kg"),
                    "k_i": parameters["k_i"],
                    "k_o": parameters["k_o"],
                    "swept_ratio": parameters["swept_ratio"],
                    "best_efficiency": (
                        _efficiency(best)
                        if best is not None
                        else None
                    ),
                    "elapsed_seconds": record["elapsed_seconds"],
                }
            ),
            flush=True,
        )

    print(f"Saved {directory / 'report.json'}", flush=True)


if __name__ == "__main__":
    main()
