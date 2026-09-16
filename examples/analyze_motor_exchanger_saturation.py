"""Isolated thermal and hydraulic saturation sweeps for the retained motor.

This is a diagnostic experiment, not a physical exchanger redesign.

The retained K2 / 100 kPa / 2 Hz / shared-7D motor champion is held fixed.
Two independent artificial sweeps are performed:

THERMAL
    Multiply both gas-wall and air-wall conductances of H_i, H_o, or both,
    while leaving wall capacity, gas volume and hydraulics unchanged.

HYDRAULIC
    Multiply the hydraulic conductance of H_i, H_o, or both, while leaving
    thermal conductances, wall capacity and gas volume unchanged.

For the current TubeHalfLink closure, a hydraulic conductance factor C is
implemented by dividing both linear and quadratic pressure-drop coefficients
by C:
    core_loss_multiplier      <- core_loss_multiplier / C
    header_loss_coefficient   <- header_loss_coefficient / C
    valve CdA                 <- valve CdA * sqrt(C)
The compressible tube-area cap remains unchanged.

These variants are deliberately non-geometric. Their purpose is to answer:
- Is the current H_i or H_o thermally saturated?
- Is the current H_i or H_o hydraulically saturated?
- How much motor efficiency/power remains available if either limitation
  is relaxed without paying the geometric penalties of a real redesign?

Default factors:
    0.5, 1, 2, 4, 8

Outputs:
    outputs/motor_exchanger_saturation/history.jsonl
    outputs/motor_exchanger_saturation/report.json
    outputs/motor_exchanger_saturation/results.csv

Run:
    PYTHONPATH=src python3 examples/analyze_motor_exchanger_saturation.py

The output directory is resumable: completed scenario IDs are skipped.
"""

from __future__ import annotations

from dataclasses import replace
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.exchangers.air_wall import AirWallExchanger, AirWallMotor
from dada_solver.exchangers.hardware import TubeHalfLink
from dada_solver.exchangers.wall_cycle import (
    solve_periodic_wall_motor,
    wall_cycle_performance,
)

from compare_motor_motion_laws_stage7A5 import A5_CAMPAIGN, ROOT
from plot_motor_four_stage_7d_cycle_alt import _build_retained_design


DEFAULT_DIRECTORY = ROOT / "outputs" / "motor_exchanger_saturation"
DEFAULT_FACTORS = (0.5, 1.0, 2.0, 4.0, 8.0)
SIDES = ("H_i", "H_o", "both")


def _parse_factors(text: str) -> tuple[float, ...]:
    values = tuple(float(x) for x in text.split(",") if x.strip())
    if not values:
        raise argparse.ArgumentTypeError("At least one factor is required.")
    if any(not math.isfinite(x) or x <= 0.0 for x in values):
        raise argparse.ArgumentTypeError("Factors must be finite and positive.")
    return values


def _load_history(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return [
        json.loads(line)
        for line in path.read_text().splitlines()
        if line.strip()
    ]


def _scenario_id(kind: str, side: str, factor: float) -> str:
    payload = json.dumps(
        {"kind": kind, "side": side, "factor": float(factor)},
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def _scale_thermal_exchanger(exchanger, factor: float):
    if not isinstance(exchanger, AirWallExchanger):
        raise TypeError(
            "Expected AirWallExchanger for isolated thermal scaling, "
            f"got {type(exchanger).__name__}."
        )
    return replace(
        exchanger,
        gas_wall_conductance_w_k=exchanger.gas_wall_conductance_w_k * factor,
        air_wall_conductance_w_k=exchanger.air_wall_conductance_w_k * factor,
    )


def _thermal_variant(wrapper: AirWallMotor, side: str, factor: float) -> AirWallMotor:
    heat_in = wrapper.heat_in
    heat_out = wrapper.heat_out

    if side in ("H_i", "both"):
        heat_in = _scale_thermal_exchanger(heat_in, factor)
    if side in ("H_o", "both"):
        heat_out = _scale_thermal_exchanger(heat_out, factor)

    return replace(wrapper, heat_in=heat_in, heat_out=heat_out)


def _scale_tube_link_hydraulics(link, conductance_factor: float):
    """Scale pressure-drop coefficients inversely to conductance_factor."""
    if not isinstance(link, TubeHalfLink):
        raise TypeError(
            "Expected TubeHalfLink for isolated hydraulic scaling, "
            f"got {type(link).__name__}."
        )

    valve_cda = link.valve_cda_m2
    if valve_cda is not None:
        valve_cda = valve_cda * math.sqrt(conductance_factor)

    return replace(
        link,
        core_loss_multiplier=link.core_loss_multiplier / conductance_factor,
        header_loss_coefficient=(
            link.header_loss_coefficient / conductance_factor
        ),
        valve_cda_m2=valve_cda,
    )


def _hydraulic_variant(
    wrapper: AirWallMotor,
    side: str,
    conductance_factor: float,
) -> AirWallMotor:
    model = wrapper.model

    small_cold = model.small_cold_link
    large_hot = model.large_hot_link
    cold_large = model.cold_large_valve
    hot_small = model.hot_small_valve

    if side in ("H_i", "both"):
        small_cold = _scale_tube_link_hydraulics(
            small_cold,
            conductance_factor,
        )
        cold_large = replace(
            cold_large,
            flow_model=_scale_tube_link_hydraulics(
                cold_large.flow_model,
                conductance_factor,
            ),
        )

    if side in ("H_o", "both"):
        large_hot = _scale_tube_link_hydraulics(
            large_hot,
            conductance_factor,
        )
        hot_small = replace(
            hot_small,
            flow_model=_scale_tube_link_hydraulics(
                hot_small.flow_model,
                conductance_factor,
            ),
        )

    changed_model = replace(
        model,
        small_cold_link=small_cold,
        large_hot_link=large_hot,
        cold_large_valve=cold_large,
        hot_small_valve=hot_small,
    )
    return replace(wrapper, model=changed_model)


def _evaluate(wrapper, initial_state, settings, maximum_cycles):
    periodic = solve_periodic_wall_motor(
        wrapper,
        initial_state,
        maximum_cycles=maximum_cycles,
        settings=settings,
    )
    if not periodic.converged:
        return {
            "status": periodic.status,
            "message": periodic.message,
            "cycles_completed": len(periodic.history),
        }, periodic.last_complete_state

    assert periodic.trajectory is not None
    performance = wall_cycle_performance(wrapper, periodic.trajectory)

    return {
        "status": "converged",
        "message": periodic.message,
        "cycles_completed": len(periodic.history),
        "indicated_thermal_efficiency": performance.thermal_efficiency,
        "indicated_power_w": performance.gas_power,
        "heat_input_w": performance.heat_in_power,
        "heat_out_w": performance.heat_out_power,
        "gas_work_per_cycle_j": performance.gas_work_per_cycle,
        "heat_input_per_cycle_j": performance.heat_in_per_cycle,
        "heat_out_per_cycle_j": performance.heat_out_per_cycle,
    }, periodic.trajectory[:10, -1].copy()


def _scenario_sequence(factors):
    # Baseline is evaluated only once.
    yield ("baseline", "both", 1.0)

    nonunity = [factor for factor in factors if not math.isclose(factor, 1.0)]

    for side in SIDES:
        for factor in nonunity:
            yield ("thermal", side, factor)

    for side in SIDES:
        for factor in nonunity:
            yield ("hydraulic", side, factor)


def _write_csv(path: Path, records: list[dict]) -> None:
    rows = []
    for record in records:
        result = record.get("result", {})
        rows.append(
            {
                "kind": record["kind"],
                "side": record["side"],
                "factor": record["factor"],
                "status": result.get("status"),
                "cycles_completed": result.get("cycles_completed"),
                "efficiency_percent": (
                    None
                    if result.get("indicated_thermal_efficiency") is None
                    else 100.0 * result["indicated_thermal_efficiency"]
                ),
                "power_W": result.get("indicated_power_w"),
                "heat_input_W": result.get("heat_input_w"),
                "heat_out_W": result.get("heat_out_w"),
                "delta_efficiency_percentage_point_vs_baseline": record.get(
                    "delta_efficiency_percentage_point_vs_baseline"
                ),
                "delta_power_W_vs_baseline": record.get(
                    "delta_power_W_vs_baseline"
                ),
            }
        )

    if not rows:
        return

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _report(history: list[dict], definition: dict) -> dict:
    baseline = next(
        (
            record
            for record in history
            if record["kind"] == "baseline"
            and record.get("result", {}).get("status") == "converged"
        ),
        None,
    )

    if baseline is not None:
        eta0 = baseline["result"]["indicated_thermal_efficiency"]
        power0 = baseline["result"]["indicated_power_w"]
        for record in history:
            result = record.get("result", {})
            eta = result.get("indicated_thermal_efficiency")
            power = result.get("indicated_power_w")
            record["delta_efficiency_percentage_point_vs_baseline"] = (
                None if eta is None else 100.0 * (eta - eta0)
            )
            record["delta_power_W_vs_baseline"] = (
                None if power is None else power - power0
            )

    grouped = {}
    for kind in ("thermal", "hydraulic"):
        grouped[kind] = {}
        for side in SIDES:
            selected = [
                record
                for record in history
                if record["kind"] == kind and record["side"] == side
            ]
            selected.sort(key=lambda record: record["factor"])
            grouped[kind][side] = selected

    return {
        "definition": definition,
        "baseline": baseline,
        "sweeps": grouped,
        "all_records": history,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--factors",
        type=_parse_factors,
        default=DEFAULT_FACTORS,
        help="Comma-separated sweep factors; default: 0.5,1,2,4,8.",
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=DEFAULT_DIRECTORY,
    )
    args = parser.parse_args()

    control, baseline_wrapper, baseline_angles, baseline_trajectory, baseline_perf = (
        _build_retained_design()
    )
    if not isinstance(baseline_wrapper, AirWallMotor):
        raise TypeError("Expected retained design to build an AirWallMotor.")

    definition = CampaignDefinition(A5_CAMPAIGN)
    settings = definition.wall_numerical_settings
    maximum_cycles = definition.configuration.numerical.maximum_cycles

    directory = args.output_directory
    directory.mkdir(parents=True, exist_ok=True)

    history_path = directory / "history.jsonl"
    report_path = directory / "report.json"
    csv_path = directory / "results.csv"

    definition_payload = {
        "factors": list(args.factors),
        "thermal_definition": (
            "Multiply both gas-wall and air-wall conductances on selected "
            "exchanger(s); keep wall capacity, gas volume and hydraulics fixed."
        ),
        "hydraulic_definition": (
            "Multiply selected exchanger hydraulic conductance by dividing "
            "both linear and quadratic pressure-drop coefficients by the factor; "
            "thermal model, wall capacity and gas volume remain fixed."
        ),
        "kinematics": control["motion_parameters"],
        "charge_pressure_pa": 100000.0,
        "frequency_hz": abs(baseline_wrapper.model.signed_angular_speed)
        / (2.0 * math.pi),
        "baseline_efficiency": baseline_perf.thermal_efficiency,
        "baseline_power_W": baseline_perf.gas_power,
    }

    history = _load_history(history_path)
    existing = {record["scenario_id"] for record in history}

    baseline_state = np.asarray(baseline_trajectory[:10, -1], dtype=float)

    for kind, side, factor in _scenario_sequence(args.factors):
        scenario_id = _scenario_id(kind, side, factor)
        if scenario_id in existing:
            continue

        if kind == "baseline":
            wrapper = baseline_wrapper
        elif kind == "thermal":
            wrapper = _thermal_variant(baseline_wrapper, side, factor)
        elif kind == "hydraulic":
            wrapper = _hydraulic_variant(baseline_wrapper, side, factor)
        else:
            raise AssertionError(kind)

        result, final_state = _evaluate(
            wrapper,
            baseline_state,
            settings,
            maximum_cycles,
        )

        record = {
            "scenario_id": scenario_id,
            "kind": kind,
            "side": side,
            "factor": float(factor),
            "result": result,
            "last_complete_state": (
                None if final_state is None else np.asarray(final_state).tolist()
            ),
        }

        with history_path.open("a") as stream:
            stream.write(json.dumps(record) + "\n")
            stream.flush()

        history.append(record)
        existing.add(scenario_id)

        report = _report(history, definition_payload)
        report_path.write_text(json.dumps(report, indent=2) + "\n")
        _write_csv(csv_path, history)

        eta = result.get("indicated_thermal_efficiency")
        eta_text = "n/a" if eta is None else f"{100.0 * eta:.6f}%"
        power = result.get("indicated_power_w")
        power_text = "n/a" if power is None else f"{power:.4f} W"
        print(
            f"{kind:9s} {side:4s} x{factor:<5g} "
            f"{result.get('status'):>14s}  "
            f"eta={eta_text:>10s}  P={power_text}"
        )

    report = _report(history, definition_payload)
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    _write_csv(csv_path, history)

    print(f"\nWrote {report_path.relative_to(ROOT)}")
    print(f"Wrote {csv_path.relative_to(ROOT)}")
    print(f"Wrote {history_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
