"""Stage K1: independent H_i / H_o exchanger-length multipliers at fixed F6 motion.

Purpose
-------
Test whether the Stage F6 champion benefits from asymmetric exchanger sizing
before releasing more kinematic degrees of freedom.

The two search coordinates are deliberately simple and physical:

    k_i = H_i tube-length multiplier
    k_o = H_o tube-length multiplier

Only tube length is scaled. Tube count, diameters, pitch, header depth,
additional dead volume, material properties and external-air mass flow remain
unchanged.

This means each k changes, consistently:
- gas-side and air-side heat-transfer area;
- tube gas hold-up;
- tube wall material and wall heat capacity;
- tube friction length;
- exchanger envelope length.

It does NOT scale:
- tube flow area (tube count and diameter stay fixed);
- header gas volume;
- additional internal volume;
- external-air mass flow.

The Stage F6 champion motion, cylinder volumes, frequency and total working-gas
mass are frozen. Thus this is an exchanger-asymmetry experiment, not a charge
or motion re-optimization.

Default search: 64 scrambled Sobol points, no controls.

Outputs
-------
    outputs/motor_exchanger_asymmetry_stageK1.csv
    outputs/motor_exchanger_asymmetry_stageK1.json
    outputs/motor_exchanger_asymmetry_stageK1_efficiency.png
    outputs/motor_exchanger_asymmetry_stageK1_power.png

Run
---
    PYTHONPATH=src python3 examples/optimize_motor_exchanger_asymmetry_stageK1.py

Resume
------
    PYTHONPATH=src python3 examples/optimize_motor_exchanger_asymmetry_stageK1.py --resume
"""

from __future__ import annotations

from dataclasses import replace
import argparse
import csv
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np
from scipy.stats import qmc

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.exchangers.microtube import MicrotubeExchanger

from compare_motor_motion_laws_stage7A5 import (
    ROOT,
    A5_CAMPAIGN,
    _candidate_design_and_mass,
    _evaluate,
    _same_inventory_design,
)
from optimize_motor_two_extrema_stageF6 import (
    PARAMETERS as F6_PARAMETERS,
    _make_kinematics,
)


F6_JSON = ROOT / "outputs" / "motor_two_extrema_stageF6.json"

OUTPUT_CSV = ROOT / "outputs" / "motor_exchanger_asymmetry_stageK1.csv"
OUTPUT_JSON = ROOT / "outputs" / "motor_exchanger_asymmetry_stageK1.json"
OUTPUT_ETA_PLOT = (
    ROOT / "outputs" / "motor_exchanger_asymmetry_stageK1_efficiency.png"
)
OUTPUT_POWER_PLOT = (
    ROOT / "outputs" / "motor_exchanger_asymmetry_stageK1_power.png"
)

# Broad first pass. Baseline tube length is 42.333... mm:
#   H_i: 21.2 .. 84.7 mm
#   H_o: 12.7 .. 63.5 mm
BOUNDS = {
    "k_i": (0.50, 2.00),
    "k_o": (0.30, 1.50),
}

MINIMUM_POWER_W = 40.0
MAXIMUM_PRESSURE_PA = 1_200_000.0
MAXIMUM_TEMPERATURE_K = 850.0
MAXIMUM_ABSOLUTE_MASS_FLOW_KG_S = 0.08


def _load_f6_champion() -> tuple[dict, dict[str, float]]:
    data = json.loads(F6_JSON.read_text())
    best = data.get("best_feasible")
    if best is None:
        raise RuntimeError("F6 contains no feasible champion.")

    values = {
        name: float(best[name])
        for name in F6_PARAMETERS
    }
    return best, values


def _scale(point: np.ndarray) -> dict[str, float]:
    return {
        name: lo + float(x) * (hi - lo)
        for x, (name, (lo, hi)) in zip(
            point,
            BOUNDS.items(),
            strict=True,
        )
    }


def _scaled_exchanger(
    exchanger: MicrotubeExchanger,
    multiplier: float,
) -> MicrotubeExchanger:
    if not math.isfinite(multiplier) or multiplier <= 0.0:
        raise ValueError("Exchanger multiplier must be positive.")
    bank = replace(
        exchanger.bank,
        tube_length_m=exchanger.bank.tube_length_m * multiplier,
    )
    return replace(exchanger, bank=bank)


def _hardware_metrics(exchanger: MicrotubeExchanger) -> dict[str, float]:
    built = exchanger.build()
    metadata = dict(built.metadata)
    return {
        "tube_length_m": exchanger.bank.tube_length_m,
        "working_gas_volume_m3": float(metadata["working_gas_volume_m3"]),
        "tube_gas_volume_m3": float(metadata["tube_gas_volume_m3"]),
        "header_gas_volume_m3": float(metadata["header_gas_volume_m3"]),
        "tube_internal_area_m2": float(metadata["tube_internal_area_m2"]),
        "wall_capacity_j_k": float(metadata["wall_capacity_j_k"]),
        "overall_static_conductance_w_k": float(
            metadata["overall_static_conductance_w_k"]
        ),
        "gas_film_resistance_k_w": float(
            metadata["gas_film_resistance_k_w"]
        ),
        "metal_resistance_k_w": float(metadata["metal_resistance_k_w"]),
        "air_film_resistance_k_w": float(
            metadata["air_film_resistance_k_w"]
        ),
        "fluid_envelope_length_m": float(
            metadata["fluid_envelope_length_m"]
        ),
    }


def _extrema_maximum(diagnostics: dict, group: str) -> float | None:
    values = diagnostics.get(group)
    if not isinstance(values, dict) or not values:
        return None

    maxima = []
    for item in values.values():
        if isinstance(item, dict) and item.get("maximum") is not None:
            maxima.append(float(item["maximum"]))
    return max(maxima) if maxima else None


def _result_row(
    index: int,
    values: dict[str, float],
    result: dict,
    hi_metrics: dict,
    ho_metrics: dict,
) -> dict:
    eta = result.get("indicated_thermal_efficiency")
    power = result.get("indicated_power_w")
    validity = result.get("validity", {})
    diagnostics = result.get("diagnostics", {})

    maximum_pressure = _extrema_maximum(
        diagnostics,
        "pressure_extrema",
    )
    maximum_temperature = _extrema_maximum(
        diagnostics,
        "temperature_extrema",
    )
    maximum_mass_flow = result.get("maximum_absolute_mass_flow_kg_s")

    constraints = {
        "power": (
            power is not None
            and math.isfinite(float(power))
            and float(power) >= MINIMUM_POWER_W
        ),
        "pressure": (
            maximum_pressure is not None
            and maximum_pressure <= MAXIMUM_PRESSURE_PA
        ),
        "temperature": (
            maximum_temperature is not None
            and maximum_temperature <= MAXIMUM_TEMPERATURE_K
        ),
        "mass_flow": (
            maximum_mass_flow is not None
            and maximum_mass_flow <= MAXIMUM_ABSOLUTE_MASS_FLOW_KG_S
        ),
        "validity": validity.get("verdict") == "valid",
    }

    feasible = (
        result.get("status") == "converged"
        and eta is not None
        and math.isfinite(float(eta))
        and all(constraints.values())
    )

    return {
        "index": index,
        "label": "sobol",
        "feasible": feasible,
        "k_i": values["k_i"],
        "k_o": values["k_o"],
        "indicated_thermal_efficiency": eta,
        "indicated_power_w": power,
        "heat_input_w": result.get("heat_input_w"),
        "heat_out_w": result.get("heat_out_w"),
        "maximum_pressure_pa": maximum_pressure,
        "maximum_temperature_k": maximum_temperature,
        "maximum_absolute_mass_flow_kg_s": maximum_mass_flow,
        "maximum_tube_reynolds": result.get("maximum_tube_reynolds"),
        "maximum_tube_mach_number": result.get("maximum_tube_mach_number"),
        "cold_isothermality_error": validity.get(
            "cold_isothermality_error"
        ),
        "hot_isothermality_error": validity.get(
            "hot_isothermality_error"
        ),
        "pressure_equalization_error": validity.get(
            "maximum_pressure_equalization_error"
        ),
        "convergence_cycles": result.get("convergence", {}).get(
            "cycles_completed"
        ),
        "constraint_power_ok": constraints["power"],
        "constraint_pressure_ok": constraints["pressure"],
        "constraint_temperature_ok": constraints["temperature"],
        "constraint_mass_flow_ok": constraints["mass_flow"],
        "constraint_validity_ok": constraints["validity"],
        "hi_tube_length_m": hi_metrics["tube_length_m"],
        "ho_tube_length_m": ho_metrics["tube_length_m"],
        "hi_working_gas_volume_m3": hi_metrics["working_gas_volume_m3"],
        "ho_working_gas_volume_m3": ho_metrics["working_gas_volume_m3"],
        "hi_wall_capacity_j_k": hi_metrics["wall_capacity_j_k"],
        "ho_wall_capacity_j_k": ho_metrics["wall_capacity_j_k"],
        "hi_static_conductance_w_k": hi_metrics[
            "overall_static_conductance_w_k"
        ],
        "ho_static_conductance_w_k": ho_metrics[
            "overall_static_conductance_w_k"
        ],
        "hi_internal_area_m2": hi_metrics["tube_internal_area_m2"],
        "ho_internal_area_m2": ho_metrics["tube_internal_area_m2"],
    }


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _load_checkpoint(path: Path) -> dict[int, dict]:
    if not path.exists():
        return {}
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        return {}
    required = {"index", "label", "k_i", "k_o"}
    if not required.issubset(rows[0]):
        return {}
    return {int(float(row["index"])): row for row in rows}


def _restore_cached_row(row: dict) -> dict:
    restored = dict(row)
    restored["index"] = int(float(restored["index"]))
    for key in list(restored):
        if key in ("label",):
            continue
        if key == "feasible" or key.startswith("constraint_"):
            restored[key] = (
                str(restored[key]).strip().lower() == "true"
            )
            continue
        if restored[key] in ("", None):
            restored[key] = None
            continue
        try:
            restored[key] = float(restored[key])
        except (TypeError, ValueError):
            pass
    return restored


def _setup_matplotlib():
    cache = Path(tempfile.gettempdir()) / "dada_solver_matplotlib"
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _plot_plane(
    path: Path,
    rows: list[dict],
    metric: str,
    title: str,
    colorbar_label: str,
) -> None:
    feasible = [
        row for row in rows
        if bool(row["feasible"])
        and row.get(metric) is not None
        and math.isfinite(float(row[metric]))
    ]
    if not feasible:
        return

    plt = _setup_matplotlib()
    fig, ax = plt.subplots(figsize=(8.5, 7.0))

    x = np.asarray([float(row["k_i"]) for row in feasible])
    y = np.asarray([float(row["k_o"]) for row in feasible])
    z = np.asarray([float(row[metric]) for row in feasible])

    if metric == "indicated_thermal_efficiency":
        z = 100.0 * z

    points = ax.scatter(x, y, c=z, s=55)
    colorbar = fig.colorbar(points, ax=ax)
    colorbar.set_label(colorbar_label)

    best = max(
        feasible,
        key=lambda row: float(row["indicated_thermal_efficiency"]),
    )
    ax.scatter(
        [float(best["k_i"])],
        [float(best["k_o"])],
        marker="x",
        s=130,
        label=(
            f"best eta: k_i={float(best['k_i']):.3f}, "
            f"k_o={float(best['k_o']):.3f}"
        ),
    )

    ax.axvline(1.0, linewidth=0.8, linestyle="--", alpha=0.45)
    ax.axhline(1.0, linewidth=0.8, linestyle="--", alpha=0.45)

    ax.set_xlabel("k_i = H_i tube-length multiplier")
    ax.set_ylabel("k_o = H_o tube-length multiplier")
    ax.set_xlim(BOUNDS["k_i"])
    ax.set_ylim(BOUNDS["k_o"])
    ax.grid(True, alpha=0.25)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.10),
        frameon=False,
    )
    ax.set_title(title)

    fig.tight_layout(rect=(0, 0.05, 1, 1))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sobol-power",
        type=int,
        default=6,
        help="2**m Sobol points; default m=6 gives 64.",
    )
    parser.add_argument("--seed", type=int, default=3251)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--csv", type=Path, default=OUTPUT_CSV)
    parser.add_argument("--output", type=Path, default=OUTPUT_JSON)
    args = parser.parse_args()

    if args.sobol_power < 1:
        raise ValueError("--sobol-power must be positive.")

    definition = CampaignDefinition(A5_CAMPAIGN)
    base_design, total_mass, _saved, _a5_best, candidate_id = (
        _candidate_design_and_mass(definition)
    )
    if not isinstance(base_design.heat_in, MicrotubeExchanger):
        raise TypeError("Expected microtube H_i design.")
    if not isinstance(base_design.heat_out, MicrotubeExchanger):
        raise TypeError("Expected microtube H_o design.")

    f6_best, f6_values = _load_f6_champion()
    small_limits = base_design.configuration.machine_volumes.small_cylinder
    large_limits = base_design.configuration.machine_volumes.large_cylinder
    kinematics = _make_kinematics(
        small_limits,
        large_limits,
        f6_values,
    )

    # Fix the F6 total gas inventory before changing exchanger hold-up.
    reference_design = _same_inventory_design(
        base_design,
        total_mass,
        kinematics,
    )

    baseline_hi = _hardware_metrics(reference_design.heat_in)
    baseline_ho = _hardware_metrics(reference_design.heat_out)

    sampler = qmc.Sobol(d=2, scramble=True, seed=args.seed)
    candidates = [
        _scale(point)
        for point in sampler.random_base2(args.sobol_power)
    ]

    print("Stage K1: independent exchanger tube-length multipliers")
    print(
        f"F6 motion fixed at index {int(f6_best['index'])}: "
        f"eta={100*float(f6_best['indicated_thermal_efficiency']):.6f}% "
        f"P={float(f6_best['indicated_power_w']):.3f} W"
    )
    print(f"Fixed total gas mass: {total_mass:.12g} kg")
    print(
        f"Baseline length: "
        f"H_i={1000*baseline_hi['tube_length_m']:.3f} mm, "
        f"H_o={1000*baseline_ho['tube_length_m']:.3f} mm"
    )
    print(
        f"{len(candidates)} Sobol evaluations, no controls\n"
        f"k_i in {BOUNDS['k_i']}, k_o in {BOUNDS['k_o']}"
    )
    print()

    checkpoint = _load_checkpoint(args.csv) if args.resume else {}
    rows: list[dict] = []

    for index, values in enumerate(candidates):
        cached = checkpoint.get(index)
        if cached is not None:
            if (
                cached.get("label") == "sobol"
                and abs(float(cached["k_i"]) - values["k_i"]) < 1e-12
                and abs(float(cached["k_o"]) - values["k_o"]) < 1e-12
            ):
                row = _restore_cached_row(cached)
                rows.append(row)
                eta = row.get("indicated_thermal_efficiency")
                eta_text = (
                    "n/a"
                    if eta is None
                    else f"{100*float(eta):8.5f}%"
                )
                print(
                    f"{index:3d} CK {eta_text:>10s} "
                    f"k_i={values['k_i']:.4f} "
                    f"k_o={values['k_o']:.4f}"
                )
                continue

        hi = _scaled_exchanger(
            reference_design.heat_in,
            values["k_i"],
        )
        ho = _scaled_exchanger(
            reference_design.heat_out,
            values["k_o"],
        )
        hi_metrics = _hardware_metrics(hi)
        ho_metrics = _hardware_metrics(ho)

        design = replace(
            reference_design,
            heat_in=hi,
            heat_out=ho,
        )

        result, _last_state = _evaluate(
            f"exchanger_asymmetry_k1_{index}",
            design,
            definition,
        )
        row = _result_row(
            index,
            values,
            result,
            hi_metrics,
            ho_metrics,
        )
        rows.append(row)
        _write_csv(args.csv, rows)

        eta = row["indicated_thermal_efficiency"]
        eta_text = (
            "n/a"
            if eta is None
            else f"{100*float(eta):8.5f}%"
        )
        power = row["indicated_power_w"]
        power_text = (
            "n/a"
            if power is None
            else f"{float(power):7.3f} W"
        )
        flag = "OK" if row["feasible"] else "X"

        print(
            f"{index:3d} {flag:2s} {eta_text:>10s} "
            f"P={power_text:>10s}  "
            f"k_i={values['k_i']:.4f} "
            f"k_o={values['k_o']:.4f}  "
            f"L_i={1000*hi_metrics['tube_length_m']:.2f} mm "
            f"L_o={1000*ho_metrics['tube_length_m']:.2f} mm"
        )

    rows.sort(key=lambda row: int(row["index"]))
    _write_csv(args.csv, rows)

    feasible = [
        row for row in rows
        if bool(row["feasible"])
        and row.get("indicated_thermal_efficiency") is not None
        and math.isfinite(float(row["indicated_thermal_efficiency"]))
    ]
    feasible.sort(
        key=lambda row: float(row["indicated_thermal_efficiency"]),
        reverse=True,
    )

    report = {
        "experiment": (
            "Stage K1 independent H_i/H_o tube-length multipliers "
            "with Stage F6 champion motion fixed"
        ),
        "definition": {
            "k_i": "H_i tube-length multiplier",
            "k_o": "H_o tube-length multiplier",
            "scaled_geometry": ["tube_length_m"],
            "unchanged_geometry": [
                "tube_count",
                "inner_diameter_m",
                "wall_thickness_m",
                "pitch_m",
                "header_depth_m",
                "additional_internal_volume_m3",
            ],
            "external_air_mass_flow_unchanged": True,
            "same_total_working_gas_mass": True,
            "same_F6_motion": True,
            "mechanical_losses_modeled": False,
            "external_aerodynamic_losses_in_objective": False,
        },
        "search": {
            "bounds": BOUNDS,
            "sobol_power": args.sobol_power,
            "sample_count": 2 ** args.sobol_power,
            "seed": args.seed,
            "scramble": True,
            "control_count": 0,
        },
        "constraints": {
            "minimum_indicated_power_w": MINIMUM_POWER_W,
            "maximum_pressure_pa": MAXIMUM_PRESSURE_PA,
            "maximum_temperature_k": MAXIMUM_TEMPERATURE_K,
            "maximum_absolute_mass_flow_kg_s":
                MAXIMUM_ABSOLUTE_MASS_FLOW_KG_S,
            "thermodynamic_validity_required": True,
        },
        "F6_reference": {
            "candidate_id": candidate_id,
            "index": int(f6_best["index"]),
            "parameters": f6_values,
            "indicated_thermal_efficiency": float(
                f6_best["indicated_thermal_efficiency"]
            ),
            "indicated_power_w": float(f6_best["indicated_power_w"]),
            "total_mass_kg": total_mass,
        },
        "baseline_hardware": {
            "H_i": baseline_hi,
            "H_o": baseline_ho,
        },
        "best_feasible": feasible[0] if feasible else None,
        "ranking": feasible,
        "all_rows": rows,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")

    _plot_plane(
        OUTPUT_ETA_PLOT,
        rows,
        "indicated_thermal_efficiency",
        "Stage K1 — efficiency over independent exchanger lengths",
        "Indicated thermal efficiency (%)",
    )
    _plot_plane(
        OUTPUT_POWER_PLOT,
        rows,
        "indicated_power_w",
        "Stage K1 — indicated power over independent exchanger lengths",
        "Indicated power (W)",
    )

    print("\nTOP FEASIBLE")
    for rank, row in enumerate(feasible[:12], 1):
        print(
            f"{rank:2d} "
            f"eta={100*float(row['indicated_thermal_efficiency']):.6f}% "
            f"P={float(row['indicated_power_w']):.3f} W  "
            f"k_i={float(row['k_i']):.4f} "
            f"k_o={float(row['k_o']):.4f}  "
            f"L_i={1000*float(row['hi_tube_length_m']):.2f} mm "
            f"L_o={1000*float(row['ho_tube_length_m']):.2f} mm"
        )

    print(f"\nWrote {args.csv.relative_to(ROOT)}")
    print(f"Wrote {args.output.relative_to(ROOT)}")
    if feasible:
        print(f"Wrote {OUTPUT_ETA_PLOT.relative_to(ROOT)}")
        print(f"Wrote {OUTPUT_POWER_PLOT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
