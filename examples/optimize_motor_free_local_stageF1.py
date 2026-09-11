"""Stage F1: local FreeKinematics perturbation of the 250-degree harmonic law.

Question
--------
Starting from the best exact harmonic pair found in the Stage 7A5 comparison,
do the volume laws spontaneously flatten near minimum volume when given a few
local shape degrees of freedom?

The experiment uses 12 equally spaced FreeKinematics controls per cylinder
(30-degree spacing). All controls are fixed to samples of the 250-degree
harmonic solution except four:

    small cylinder:  60 deg and  90 deg
    large cylinder: 150 deg and 210 deg

For the large cylinder the exact harmonic minimum is at 180 deg, so the two
varied controls are the immediate neighbours of the minimum control.

For the small cylinder the exact harmonic minimum is at 70 deg, so 60 and
90 deg are the two closest controls bracketing the minimum.

The four values are searched independently with Sobol. Bounds prevent them
from dropping below -1, so a move toward -1 means a direct tendency to flatten
the minimum rather than create a new deeper minimum away from the original one.

Comparison basis is identical to the Stage 7A5 harmonic experiment:
- same total working-gas mass,
- same cylinder volume limits,
- same 2 Hz operation,
- same dynamic-wall microtube hardware,
- no mechanical-loss model.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import qmc

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.free_kinematics import (
    FreeKinematics,
    FreeKinematicsConfiguration,
    FreeMotionDefinition,
)

from compare_motor_motion_laws_stage7A5 import (
    ROOT,
    A5_CAMPAIGN,
    _candidate_design_and_mass,
    _evaluate,
    _same_inventory_design,
)


PHASE_DEG = 250.0
CONTROL_COUNT = 12
ANGLES_DEG = tuple(30.0 * i for i in range(CONTROL_COUNT))

OUTPUT_JSON = ROOT / "outputs" / "motor_free_local_stageF1.json"
OUTPUT_CSV = ROOT / "outputs" / "motor_free_local_stageF1.csv"

# Search only the four controls nearest the harmonic minima.
PARAMETERS = (
    "small_60_deg",
    "small_90_deg",
    "large_150_deg",
    "large_210_deg",
)

BOUNDS = {
    # Harmonic baselines:
    # S60  = cos(60 - 250) = -0.984807753...
    # S90  = cos(90 - 250) = -0.939692621...
    # L150 = L210 = cos(150) = -0.866025404...
    #
    # Lower bound -1 means "flatten toward the minimum" without allowing
    # these neighbours to become a deeper minimum than the original one.
    "small_60_deg": (-1.0, -0.78),
    "small_90_deg": (-1.0, -0.72),
    "large_150_deg": (-1.0, -0.62),
    "large_210_deg": (-1.0, -0.62),
}

MINIMUM_POWER_W = 40.0


def _harmonic_controls(side: str) -> list[float]:
    values = []
    for angle_deg in ANGLES_DEG:
        theta = math.radians(angle_deg)
        if side == "large":
            values.append(math.cos(theta))
        elif side == "small":
            values.append(math.cos(theta - math.radians(PHASE_DEG)))
        else:
            raise ValueError(side)
    return values


BASE_SMALL = _harmonic_controls("small")
BASE_LARGE = _harmonic_controls("large")

INDEX = {
    "small_60_deg": ANGLES_DEG.index(60.0),
    "small_90_deg": ANGLES_DEG.index(90.0),
    "large_150_deg": ANGLES_DEG.index(150.0),
    "large_210_deg": ANGLES_DEG.index(210.0),
}

BASELINE = {
    "small_60_deg": BASE_SMALL[INDEX["small_60_deg"]],
    "small_90_deg": BASE_SMALL[INDEX["small_90_deg"]],
    "large_150_deg": BASE_LARGE[INDEX["large_150_deg"]],
    "large_210_deg": BASE_LARGE[INDEX["large_210_deg"]],
}


def _make_kinematics(small_limits, large_limits, values: dict[str, float]):
    small = list(BASE_SMALL)
    large = list(BASE_LARGE)

    small[INDEX["small_60_deg"]] = values["small_60_deg"]
    small[INDEX["small_90_deg"]] = values["small_90_deg"]
    large[INDEX["large_150_deg"]] = values["large_150_deg"]
    large[INDEX["large_210_deg"]] = values["large_210_deg"]

    return FreeKinematics(
        FreeKinematicsConfiguration(
            FreeMotionDefinition(
                tuple(small),
                small_limits.minimum,
                small_limits.maximum,
            ),
            FreeMotionDefinition(
                tuple(large),
                large_limits.minimum,
                large_limits.maximum,
            ),
        )
    )


def _scale(point: np.ndarray) -> dict[str, float]:
    result = {}
    for x, name in zip(point, PARAMETERS, strict=True):
        lo, hi = BOUNDS[name]
        result[name] = lo + float(x) * (hi - lo)
    return result


def _flattening_fraction(name: str, value: float) -> float:
    """0 = harmonic baseline, 1 = neighbour moved exactly to -1."""
    baseline = BASELINE[name]
    denominator = baseline - (-1.0)
    if denominator == 0:
        return 0.0
    return (baseline - value) / denominator


def _result_row(index: int, label: str, values: dict, result: dict, kin):
    eta = result.get("indicated_thermal_efficiency")
    power = result.get("indicated_power_w")
    validity = result.get("validity", {})
    feasible = (
        result.get("status") == "converged"
        and eta is not None
        and math.isfinite(eta)
        and power is not None
        and power >= MINIMUM_POWER_W
        and validity.get("verdict") == "valid"
    )

    diagnostics = kin.diagnostics
    row = {
        "index": index,
        "label": label,
        "feasible": feasible,
        **values,
        "small_60_flattening_fraction": _flattening_fraction(
            "small_60_deg", values["small_60_deg"]
        ),
        "small_90_flattening_fraction": _flattening_fraction(
            "small_90_deg", values["small_90_deg"]
        ),
        "large_150_flattening_fraction": _flattening_fraction(
            "large_150_deg", values["large_150_deg"]
        ),
        "large_210_flattening_fraction": _flattening_fraction(
            "large_210_deg", values["large_210_deg"]
        ),
        "indicated_thermal_efficiency": eta,
        "indicated_power_w": power,
        "heat_input_w": result.get("heat_input_w"),
        "heat_out_w": result.get("heat_out_w"),
        "maximum_tube_reynolds": result.get("maximum_tube_reynolds"),
        "maximum_tube_mach_number": result.get("maximum_tube_mach_number"),
        "maximum_absolute_mass_flow_kg_s": result.get(
            "maximum_absolute_mass_flow_kg_s"
        ),
        "cold_isothermality_error": validity.get(
            "cold_isothermality_error"
        ),
        "hot_isothermality_error": validity.get(
            "hot_isothermality_error"
        ),
        "pressure_equalization_error": validity.get(
            "maximum_pressure_equalization_error"
        ),
        "small_max_abs_dV_dtheta": (
            diagnostics[0].maximum_absolute_first_derivative
        ),
        "large_max_abs_dV_dtheta": (
            diagnostics[1].maximum_absolute_first_derivative
        ),
        "small_max_abs_d2V_dtheta2": (
            diagnostics[0].maximum_absolute_second_derivative
        ),
        "large_max_abs_d2V_dtheta2": (
            diagnostics[1].maximum_absolute_second_derivative
        ),
        "convergence_cycles": result.get("convergence", {}).get(
            "cycles_completed"
        ),
    }
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sobol-power",
        type=int,
        default=6,
        help="2**m Sobol points; default m=6 gives 64.",
    )
    parser.add_argument("--seed", type=int, default=1701)
    parser.add_argument("--output", type=Path, default=OUTPUT_JSON)
    parser.add_argument("--csv", type=Path, default=OUTPUT_CSV)
    args = parser.parse_args()

    if args.sobol_power < 1:
        raise ValueError("--sobol-power must be positive.")

    definition = CampaignDefinition(A5_CAMPAIGN)
    base_design, total_mass, _saved, _a5_best, candidate_id = (
        _candidate_design_and_mass(definition)
    )
    small_limits = base_design.configuration.machine_volumes.small_cylinder
    large_limits = base_design.configuration.machine_volumes.large_cylinder

    # Explicit controls before Sobol:
    # 0: harmonic samples represented through FreeKinematics.
    # 1: all four neighbours placed at -1: a strong local plateau probe.
    candidates: list[tuple[str, dict[str, float]]] = [
        ("free_harmonic_baseline", dict(BASELINE)),
        (
            "explicit_flattened_neighbours",
            {name: -1.0 for name in PARAMETERS},
        ),
    ]

    sampler = qmc.Sobol(d=4, scramble=True, seed=args.seed)
    for point in sampler.random_base2(args.sobol_power):
        candidates.append(("sobol", _scale(point)))

    rows = []
    full_results = []

    print("Stage F1: local FreeKinematics perturbation")
    print(f"250-degree harmonic reference; fixed gas mass {total_mass:.12g} kg")
    print(f"{len(candidates)} evaluations "
          f"({2**args.sobol_power} Sobol + 2 controls)\n")
    print("Harmonic raw control baselines:")
    for name in PARAMETERS:
        print(f"  {name:18s} {BASELINE[name]: .9f}")
    print()

    for index, (label, values) in enumerate(candidates):
        kin = _make_kinematics(
            small_limits,
            large_limits,
            values,
        )
        design = _same_inventory_design(
            base_design,
            total_mass,
            kin,
        )
        result, _ = _evaluate(
            f"free_local_{index}",
            design,
            definition,
        )
        result["control_values"] = values
        result["free_motion_diagnostics"] = [
            asdict(item) for item in kin.diagnostics
        ]
        row = _result_row(
            index,
            label,
            values,
            result,
            kin,
        )
        rows.append(row)
        full_results.append(result)

        eta = row["indicated_thermal_efficiency"]
        eta_text = "n/a" if eta is None else f"{100*eta:8.5f}%"
        power = row["indicated_power_w"]
        power_text = "n/a" if power is None else f"{power:7.3f} W"
        flag = "OK" if row["feasible"] else "X"
        print(
            f"{index:3d} {flag:2s} {eta_text:>10s} "
            f"P={power_text:>10s}  "
            f"S60={values['small_60_deg']:.4f} "
            f"S90={values['small_90_deg']:.4f} "
            f"L150={values['large_150_deg']:.4f} "
            f"L210={values['large_210_deg']:.4f}"
        )

    feasible = [
        row for row in rows
        if row["feasible"]
        and row["indicated_thermal_efficiency"] is not None
    ]
    feasible.sort(
        key=lambda row: row["indicated_thermal_efficiency"],
        reverse=True,
    )

    exact_harmonic = None
    comparison_path = (
        ROOT / "outputs" / "motor_motion_law_comparison_stage7A5.json"
    )
    if comparison_path.exists():
        exact_harmonic = json.loads(
            comparison_path.read_text()
        ).get("best_harmonic")

    output = {
        "experiment": (
            "Stage F1 local FreeKinematics perturbation of "
            "250-degree harmonic law"
        ),
        "comparison_basis": {
            "same_total_working_gas_mass": True,
            "total_mass_kg": total_mass,
            "same_cylinder_volume_limits": True,
            "same_operation": True,
            "same_microtube_hardware": True,
            "mechanical_losses_modeled": False,
            "stage7A5_candidate_id": candidate_id,
        },
        "phase_degrees": PHASE_DEG,
        "control_angles_degrees": ANGLES_DEG,
        "harmonic_baseline_controls": BASELINE,
        "search_bounds": BOUNDS,
        "exact_harmonic_reference": exact_harmonic,
        "sobol": {
            "power": args.sobol_power,
            "sample_count": 2 ** args.sobol_power,
            "seed": args.seed,
            "scramble": True,
        },
        "best_feasible": feasible[0] if feasible else None,
        "ranking": feasible,
        "all_rows": rows,
        "all_full_results": full_results,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n")

    with args.csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print("\nTOP FEASIBLE")
    for rank, row in enumerate(feasible[:12], 1):
        print(
            f"{rank:2d} "
            f"eta={100*row['indicated_thermal_efficiency']:.6f}% "
            f"P={row['indicated_power_w']:.3f} W  "
            f"S60={row['small_60_deg']:.5f} "
            f"S90={row['small_90_deg']:.5f} "
            f"L150={row['large_150_deg']:.5f} "
            f"L210={row['large_210_deg']:.5f}"
        )

    if exact_harmonic is not None:
        eta = exact_harmonic.get("indicated_thermal_efficiency")
        power = exact_harmonic.get("indicated_power_w")
        if eta is not None and power is not None:
            print(
                "\nExact harmonic benchmark: "
                f"eta={100*eta:.6f}% P={power:.3f} W"
            )

    print(f"\nWrote {args.output.relative_to(ROOT)}")
    print(f"Wrote {args.csv.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
