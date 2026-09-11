"""Stage P2: Sobol search of the historical ideal piecewise-linear motion law.

Purpose
-------
Test whether the poor historical piecewise-linear result came from its old
Lambda values rather than from the motion-law family itself.

Comparison basis:
- same Stage 7A5 total working-gas mass,
- same cylinder volume limits,
- same 2 Hz operation,
- same dynamic-wall microtube hardware,
- same physical validity limits,
- no mechanical-loss model.

Only three kinematic quantities vary:
    small_lambda_target
    large_lambda_target
    adiabatic_sector_fraction
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import qmc

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.kinematics import IdealPiecewiseLinearVolumeKinematics

from compare_motor_motion_laws_stage7A5 import (
    ROOT,
    A5_CAMPAIGN,
    _candidate_design_and_mass,
    _evaluate,
    _same_inventory_design,
)

DEFAULT_JSON = ROOT / "outputs" / "motor_piecewise_stageP3.json"
DEFAULT_CSV = ROOT / "outputs" / "motor_piecewise_stageP3.csv"

BOUNDS = {
    "small_lambda_target": (0.54, 0.72),
    "large_lambda_target": (0.50, 0.68),
    "adiabatic_sector_fraction": (0.02, 0.16),
}

MINIMUM_POWER_W = 40.0
MAXIMUM_PRESSURE_PA = 1_200_000.0
MAXIMUM_TEMPERATURE_K = 850.0
MAXIMUM_ABSOLUTE_MASS_FLOW_KG_S = 0.08
MAXIMUM_MACH = 0.2
MAXIMUM_REYNOLDS = 2300.0


def _scale(unit: np.ndarray) -> tuple[float, float, float]:
    values = []
    for x, name in zip(unit, tuple(BOUNDS), strict=True):
        lo, hi = BOUNDS[name]
        values.append(lo + float(x) * (hi - lo))
    return tuple(values)


def _feasibility(result: dict) -> tuple[bool, list[str]]:
    reasons: list[str] = []

    eta = result.get("indicated_thermal_efficiency")
    power = result.get("indicated_power_w")
    if eta is None or not math.isfinite(eta):
        reasons.append("motor_efficiency_unavailable")
    if power is None or power < MINIMUM_POWER_W:
        reasons.append("minimum_motor_power")

    diagnostics = result.get("diagnostics", {})
    pressure_extrema = diagnostics.get("pressure_extrema", {})
    temperature_extrema = diagnostics.get("temperature_extrema", {})

    max_pressure = (
        max(float(ext["maximum"]) for ext in pressure_extrema.values())
        if pressure_extrema else None
    )
    max_temperature = (
        max(float(ext["maximum"]) for ext in temperature_extrema.values())
        if temperature_extrema else None
    )

    if max_pressure is None or max_pressure > MAXIMUM_PRESSURE_PA:
        reasons.append("maximum_pressure")
    if max_temperature is None or max_temperature > MAXIMUM_TEMPERATURE_K:
        reasons.append("maximum_temperature")

    max_flow = result.get("maximum_absolute_mass_flow_kg_s")
    if max_flow is None or max_flow > MAXIMUM_ABSOLUTE_MASS_FLOW_KG_S:
        reasons.append("maximum_absolute_mass_flow")

    max_mach = result.get("maximum_tube_mach_number")
    if max_mach is None or max_mach > MAXIMUM_MACH:
        reasons.append("maximum_mach")

    max_re = result.get("maximum_tube_reynolds")
    if max_re is None or max_re >= MAXIMUM_REYNOLDS:
        reasons.append("maximum_reynolds")

    verdict = result.get("validity", {}).get("verdict")
    if verdict != "valid":
        reasons.append("valid_thermodynamic_model")

    result["maximum_pressure_pa"] = max_pressure
    result["maximum_temperature_k"] = max_temperature
    return not reasons, reasons


def _compact(index: int, params: dict, result: dict) -> dict:
    feasible, reasons = _feasibility(result)
    return {
        "index": index,
        **params,
        "status": result.get("status"),
        "feasible": feasible,
        "failure_reasons": ";".join(reasons),
        "indicated_thermal_efficiency": result.get(
            "indicated_thermal_efficiency"
        ),
        "indicated_power_w": result.get("indicated_power_w"),
        "heat_input_w": result.get("heat_input_w"),
        "heat_out_w": result.get("heat_out_w"),
        "maximum_pressure_pa": result.get("maximum_pressure_pa"),
        "maximum_temperature_k": result.get("maximum_temperature_k"),
        "maximum_absolute_mass_flow_kg_s": result.get(
            "maximum_absolute_mass_flow_kg_s"
        ),
        "maximum_tube_reynolds": result.get("maximum_tube_reynolds"),
        "maximum_tube_mach_number": result.get(
            "maximum_tube_mach_number"
        ),
        "cold_isothermality_error": result.get("validity", {}).get(
            "cold_isothermality_error"
        ),
        "hot_isothermality_error": result.get("validity", {}).get(
            "hot_isothermality_error"
        ),
        "pressure_equalization_error": result.get("validity", {}).get(
            "maximum_pressure_equalization_error"
        ),
        "convergence_cycles": result.get("convergence", {}).get(
            "cycles_completed"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sobol-power",
        type=int,
        default=5,
        help="Evaluate 2**m Sobol points; default m=5 gives 32.",
    )
    parser.add_argument("--seed", type=int, default=1301)
    parser.add_argument("--output", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    args = parser.parse_args()

    if args.sobol_power < 1:
        raise ValueError("--sobol-power must be positive.")

    definition = CampaignDefinition(A5_CAMPAIGN)
    base_design, total_mass, _saved_a5_state, a5_best, candidate_id = (
        _candidate_design_and_mass(definition)
    )

    small_limits = base_design.configuration.machine_volumes.small_cylinder
    large_limits = base_design.configuration.machine_volumes.large_cylinder

    sampler = qmc.Sobol(d=3, scramble=True, seed=args.seed)
    points = sampler.random_base2(args.sobol_power)

    parameter_sets: list[tuple[float, float, float]] = [
        (0.619035229459405, 0.5974640518426895, 0.08983197969384492),
        (0.5898326434195041, 0.5798693427443504, 0.12454577540047468),
        (0.5845900582149625, 0.6278566133975982, 0.15438122472725807),
    ]
    parameter_sets.extend(_scale(point) for point in points)

    full_results: list[dict] = []
    compact_results: list[dict] = []

    print(
        f"Stage P1: {len(parameter_sets)} evaluations "
        f"({2**args.sobol_power} Sobol + historical reference)"
    )
    print(f"Fixed total air mass: {total_mass:.12g} kg\n")

    for index, (lambda_s, lambda_l, adiabatic) in enumerate(parameter_sets):
        kinematics = IdealPiecewiseLinearVolumeKinematics(
            small_limits,
            large_limits,
            lambda_s,
            lambda_l,
            adiabatic,
        )
        design = _same_inventory_design(
            base_design,
            total_mass,
            kinematics,
        )
        result, _ = _evaluate(
            f"piecewise_{index}",
            design,
            definition,
        )
        params = {
            "small_lambda_target": lambda_s,
            "large_lambda_target": lambda_l,
            "adiabatic_sector_fraction": adiabatic,
        }
        result["parameters"] = params
        compact = _compact(index, params, result)
        full_results.append(result)
        compact_results.append(compact)

        eta = compact["indicated_thermal_efficiency"]
        eta_text = "n/a" if eta is None else f"{100*eta:8.4f}%"
        p = compact["indicated_power_w"]
        p_text = "n/a" if p is None else f"{p:7.2f} W"
        flag = "OK" if compact["feasible"] else "X"
        print(
            f"{index:2d} {flag:2s} "
            f"eta={eta_text:>9s}  P={p_text:>9s}  "
            f"Ls={lambda_s:.4f}  Ll={lambda_l:.4f}  "
            f"a={adiabatic:.4f}"
        )

    feasible = [
        row
        for row in compact_results
        if row["feasible"]
        and row["indicated_thermal_efficiency"] is not None
    ]
    feasible.sort(
        key=lambda row: row["indicated_thermal_efficiency"],
        reverse=True,
    )

    comparison_path = (
        ROOT / "outputs" / "motor_motion_law_comparison_stage7A5.json"
    )
    harmonic_reference = None
    if comparison_path.exists():
        comparison = json.loads(comparison_path.read_text())
        harmonic_reference = comparison.get("best_harmonic")

    output = {
        "experiment": "Stage P1 piecewise-linear Sobol screening",
        "comparison_basis": {
            "same_total_working_gas_mass": True,
            "total_mass_kg": total_mass,
            "same_cylinder_volume_limits": True,
            "same_operation": True,
            "same_microtube_hardware": True,
            "mechanical_losses_modeled": False,
            "stage7A5_candidate_id": candidate_id,
            "stage7A5_physical": a5_best["physical"],
        },
        "bounds": BOUNDS,
        "constraints": {
            "minimum_motor_power_w": MINIMUM_POWER_W,
            "maximum_pressure_pa": MAXIMUM_PRESSURE_PA,
            "maximum_temperature_k": MAXIMUM_TEMPERATURE_K,
            "maximum_absolute_mass_flow_kg_s": (
                MAXIMUM_ABSOLUTE_MASS_FLOW_KG_S
            ),
            "maximum_mach": MAXIMUM_MACH,
            "maximum_reynolds": MAXIMUM_REYNOLDS,
        },
        "sobol": {
            "power": args.sobol_power,
            "sample_count": 2**args.sobol_power,
            "seed": args.seed,
            "scramble": True,
        },
        "harmonic_reference": harmonic_reference,
        "best_feasible": feasible[0] if feasible else None,
        "ranking": feasible,
        "all_compact_results": compact_results,
        "all_full_results": full_results,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n")

    fieldnames = list(compact_results[0])
    with args.csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(compact_results)

    print("\nTOP FEASIBLE")
    for rank, row in enumerate(feasible[:10], 1):
        print(
            f"{rank:2d}  "
            f"eta={100*row['indicated_thermal_efficiency']:.5f}%  "
            f"P={row['indicated_power_w']:.3f} W  "
            f"Ls={row['small_lambda_target']:.5f}  "
            f"Ll={row['large_lambda_target']:.5f}  "
            f"a={row['adiabatic_sector_fraction']:.5f}"
        )

    if harmonic_reference:
        eta = harmonic_reference.get("indicated_thermal_efficiency")
        power = harmonic_reference.get("indicated_power_w")
        if eta is not None and power is not None:
            print(
                "\nHarmonic benchmark: "
                f"eta={100*eta:.5f}%  P={power:.3f} W"
            )

    print(f"\nWrote {args.output.relative_to(ROOT)}")
    print(f"Wrote {args.csv.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
