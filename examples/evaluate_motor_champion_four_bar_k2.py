"""Evaluate the synthesized shared-crank four-bar with the Stage K2 motor.

This is a composition/validation run.  It keeps the Stage K2 thermodynamic
design, total working-gas inventory, asymmetric microtube exchangers, dynamic
wall storage, and numerical convergence criterion.  Only the Stage F6 free
motion law is replaced by its projection-only four-bar synthesis.

Run from the repository root:

    PYTHONPATH=src python3 examples/evaluate_motor_champion_four_bar_k2.py
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import csv
import json
from pathlib import Path

import numpy as np

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.campaign.evaluator import rescale_wall_state
from dada_solver.coupler_projection import (
    load_shared_crank_coupler_projection_geometry,
)
from dada_solver.exchangers.air_wall import AirWallMotor
from compare_motor_motion_laws_stage7A5 import (
    A5_CAMPAIGN,
    ROOT,
    _candidate_design_and_mass,
    _evaluate,
    _same_inventory_design,
)
from optimize_motor_exchanger_asymmetry_stageK1 import (
    _hardware_metrics,
    _scaled_exchanger,
)


GEOMETRY = ROOT / "outputs" / "motor_champion_four_bar_projection_geometry.toml"
TARGET_CSV = ROOT / "outputs" / "motor_champion_motion_target.csv"
K2_REPORT = ROOT / "outputs" / "motor_exchanger_asymmetry_stageK2.json"
OUTPUT = ROOT / "outputs" / "motor_champion_four_bar_k2.json"


def _motion_fit(kinematics) -> dict[str, float]:
    with TARGET_CSV.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    theta = np.asarray([float(row["theta_rad"]) for row in rows])
    result: dict[str, float] = {}
    for side in ("small", "large"):
        limits = getattr(kinematics, f"{side}_volume_limits")
        volume = np.asarray(
            [getattr(kinematics, f"{side}_cylinder_volume")(x) for x in theta]
        )
        derivative = np.asarray(
            [
                getattr(kinematics, f"{side}_cylinder_volume_derivative")(x)
                for x in theta
            ]
        )
        target_volume = np.asarray(
            [float(row[f"{side}_volume_m3"]) for row in rows]
        )
        target_derivative = limits.swept * np.asarray(
            [float(row[f"{side}_dq_dtheta_per_rad"]) for row in rows]
        )
        result[f"{side}_normalized_volume_rms"] = float(
            np.sqrt(np.mean(((volume - target_volume) / limits.swept) ** 2))
        )
        result[f"{side}_normalized_derivative_rms"] = float(
            np.sqrt(np.mean(((derivative - target_derivative) / limits.swept) ** 2))
        )
    return result


def _rescaled_initial_state(saved, old_design, new_design) -> np.ndarray:
    old_wrapper = old_design.build()
    new_wrapper = new_design.build()
    if not isinstance(old_wrapper, AirWallMotor) or not isinstance(new_wrapper, AirWallMotor):
        raise TypeError("The K2 composition requires dynamic-wall microtube models.")
    old_capacities = [
        old_wrapper.heat_in.wall_capacity_j_k,
        old_wrapper.heat_out.wall_capacity_j_k,
    ]
    new_capacities = [
        new_wrapper.heat_in.wall_capacity_j_k,
        new_wrapper.heat_out.wall_capacity_j_k,
    ]
    total_mass = float(np.asarray(saved)[0:8:2].sum())
    return rescale_wall_state(saved, total_mass, old_capacities, new_capacities)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--finite-coupler", action="store_true")
    parser.add_argument("--geometry", type=Path, default=GEOMETRY)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--maximum-cycles", type=int)
    args = parser.parse_args()

    definition = CampaignDefinition(A5_CAMPAIGN)
    base, total_mass, saved, stage_a5, candidate_id = _candidate_design_and_mass(
        definition
    )
    k2 = json.loads(K2_REPORT.read_text())["best_feasible"]
    loader = load_shared_crank_coupler_projection_geometry
    if args.finite_coupler:
        from compact_coupler_geometry import load_compact_coupler
        loader = load_compact_coupler
    kinematics = loader(
        args.geometry,
        base.configuration.machine_volumes.small_cylinder,
        base.configuration.machine_volumes.large_cylinder,
    )
    inventory_design = _same_inventory_design(base, total_mass, kinematics)
    design = replace(
        inventory_design,
        heat_in=_scaled_exchanger(inventory_design.heat_in, float(k2["k_i"])),
        heat_out=_scaled_exchanger(inventory_design.heat_out, float(k2["k_o"])),
    )
    if args.maximum_cycles is not None:
        numerical = replace(
            design.configuration.numerical,
            maximum_cycles=args.maximum_cycles,
        )
        design = replace(
            design,
            configuration=replace(design.configuration, numerical=numerical),
        )

    initial_state = _rescaled_initial_state(saved, base, design)
    result, final_state = _evaluate(
        "motor_champion_four_bar_k2",
        design,
        definition,
        initial_state=initial_state,
    )
    report = {
        "experiment": ("Finite-rod coupler four-bar in the Stage K2 motor" if args.finite_coupler
                       else "Projection-only shared-crank four-bar in the Stage K2 motor"),
        "geometry_file": str(args.geometry.relative_to(ROOT)),
        "angle_convention": "study_angle_radians; motor reversal applied once by build_model",
        "thermodynamic_basis": {
            "stage_a5_candidate_id": candidate_id,
            "stage_a5_objective": stage_a5.get("objective"),
            "stage_k2_k_i": float(k2["k_i"]),
            "stage_k2_k_o": float(k2["k_o"]),
            "stage_k2_reference_efficiency": float(k2["indicated_thermal_efficiency"]),
            "stage_k2_reference_power_w": float(k2["indicated_power_w"]),
            "total_working_gas_mass_kg": total_mass,
            "heat_in": _hardware_metrics(design.heat_in),
            "heat_out": _hardware_metrics(design.heat_out),
            "wall_numerical_settings": asdict(definition.wall_numerical_settings),
        },
        "motion_fit_against_stage_f6_target": _motion_fit(kinematics),
        "result": result,
        "last_complete_state": final_state.tolist() if final_state is not None else None,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
