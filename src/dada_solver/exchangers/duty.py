"""Time-resolved exchanger duty export for literature/manufacturer comparisons.

The intended through-flow is intermittent and unidirectional. Both ports are
exported separately because a lumped exchanger stores gas; inlet-port reflux
predicted by the existing bidirectional links must not be hidden by abs(flow).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from dada_solver.dynamics import ThermodynamicModel
from dada_solver.integration import CycleIntegrationResult
from dada_solver.state import ThermodynamicState


def _positive_integral(time, values) -> float:
    """Integrate the positive part of a piecewise-linear sampled signal."""
    a, b = values[:-1], values[1:]
    dt = np.diff(time)
    positive_a, positive_b = np.maximum(a, 0.0), np.maximum(b, 0.0)
    fraction = np.divide(
        positive_a + positive_b, np.abs(a) + np.abs(b),
        out=np.zeros_like(a), where=(np.abs(a) + np.abs(b)) > 0.0,
    )
    return float(np.sum(dt * 0.5 * (positive_a + positive_b) * fraction))


def _duration_above(time, values, threshold) -> float:
    a, b = values[:-1] - threshold, values[1:] - threshold
    fraction = np.where((a > 0.0) & (b > 0.0), 1.0, 0.0)
    crossing = (a > 0.0) != (b > 0.0)
    fraction[crossing] = (
        np.maximum(a[crossing], b[crossing])
        / (np.abs(a[crossing]) + np.abs(b[crossing]))
    )
    return float(np.sum(np.diff(time) * fraction))


def summarize_port(time, flow, activity_fraction: float = 0.01) -> dict:
    """Time-weighted port metrics, with an explicit reporting-only cutoff."""
    time, flow = np.asarray(time, dtype=float), np.asarray(flow, dtype=float)
    if (
        time.ndim != 1 or time.size < 2 or flow.shape != time.shape
        or not np.all(np.isfinite(time)) or not np.all(np.isfinite(flow))
        or np.any(np.diff(time) < 0.0) or time[-1] <= time[0]
    ):
        raise ValueError("Expected finite, chronological time/flow samples with positive duration.")
    if not math.isfinite(activity_fraction) or not 0.0 < activity_fraction < 1.0:
        raise ValueError("Activity fraction must lie strictly between zero and one.")
    period = time[-1] - time[0]
    peak = max(float(np.max(flow)), 0.0)
    forward = _positive_integral(time, flow)
    reverse = _positive_integral(time, -flow)
    active = _duration_above(time, flow, activity_fraction * peak)
    return {
        "minimum_mass_flow_kg_s": float(np.min(flow)),
        "peak_forward_mass_flow_kg_s": peak,
        "mean_forward_mass_flow_kg_s": forward / period,
        "forward_mass_per_cycle_kg": forward,
        "reverse_mass_per_cycle_kg": reverse,
        "reverse_to_forward_mass_ratio": reverse / forward if forward > 0.0 else None,
        "activity_threshold_kg_s": activity_fraction * peak,
        "duration_above_activity_threshold_s": active,
        "duty_fraction_above_activity_threshold": active / period,
        "reflux_above_activity_threshold_s": _duration_above(time, -flow, activity_fraction * peak),
    }


def extract_exchanger_duty(cycle: CycleIntegrationResult, model: ThermodynamicModel,
                           activity_fraction: float = 0.01) -> tuple[dict, dict]:
    if not cycle.completed:
        raise ValueError("Duty extraction requires a completed cycle.")
    time = (cycle.angles - cycle.angles[0]) / model.angular_speed
    period = float(time[-1])
    columns = {"time_s": time, "cycle_progress_rad": cycle.angles}
    for name in ("S", "L", "H_i", "H_o"):
        columns[f"{name}_pressure_Pa"] = np.empty(time.size)
        columns[f"{name}_temperature_K"] = np.empty(time.size)
    for name in ("S_to_H_i", "H_i_to_L", "L_to_H_o", "H_o_to_S"):
        columns[f"{name}_mass_flow_kg_s"] = np.empty(time.size)
    for name in ("H_i", "H_o"):
        columns[f"{name}_heat_received_W"] = np.empty(time.size)
    for i, (angle, values, topology) in enumerate(zip(cycle.angles, cycle.states.T, cycle.topologies, strict=True)):
        state = ThermodynamicState.from_array(values)
        temperatures, pressures = state.temperatures_and_pressures(model.gas, model.volumes(float(angle)))
        rates = model.evaluate(float(angle), state, topology)
        for j, name in enumerate(("S", "L", "H_i", "H_o")):
            columns[f"{name}_pressure_Pa"][i] = pressures[j]
            columns[f"{name}_temperature_K"][i] = temperatures[j]
        for name, flow in (
            ("S_to_H_i", rates.flows.small_to_cold),
            ("H_i_to_L", rates.flows.cold_to_large),
            ("L_to_H_o", rates.flows.large_to_hot),
            ("H_o_to_S", rates.flows.hot_to_small),
        ):
            columns[f"{name}_mass_flow_kg_s"][i] = flow
        columns["H_i_heat_received_W"][i] = rates.cold_heat_rate
        columns["H_o_heat_received_W"][i] = rates.hot_heat_rate

    summary = {
        "period_s": period,
        "cycle_frequency_Hz": 1.0 / period,
        "signed_crank_speed_rad_s": model.signed_angular_speed,
        "intended_flow": "intermittent_unidirectional",
        "activity_fraction_of_forward_peak": activity_fraction,
        "activity_threshold_is_reporting_only": True,
        "sample_count": time.size,
        "extrema_are_independent_not_a_simultaneous_operating_point": True,
        "gas": {"R_J_kg_K": model.gas.gas_constant, "Cp_J_kg_K": model.gas.heat_capacity_cp},
        "exchangers": {},
    }
    for name, inlet, outlet, thermal, volume, heat in (
        ("H_i", "S_to_H_i", "H_i_to_L", model.cold_heat_transfer, model.machine_volumes.cold_heat_exchanger, cycle.cold_heat),
        ("H_o", "L_to_H_o", "H_o_to_S", model.hot_heat_transfer, model.machine_volumes.hot_heat_exchanger, cycle.hot_heat),
    ):
        temperature = columns[f"{name}_temperature_K"]
        pressure = columns[f"{name}_pressure_Pa"]
        summary["exchangers"][name] = {
            "reservoir_temperature_K": thermal.reservoir_temperature,
            "gas_temperature_range_K": [float(np.min(temperature)), float(np.max(temperature))],
            "gas_pressure_range_Pa": [float(np.min(pressure)), float(np.max(pressure))],
            "assumed_UA_W_K": thermal.conductance,
            "assumed_gas_volume_m3": volume,
            "mean_heat_received_W": float((heat[-1] - heat[0]) / period),
            "inlet_port": summarize_port(time, columns[f"{inlet}_mass_flow_kg_s"], activity_fraction),
            "outlet_port": summarize_port(time, columns[f"{outlet}_mass_flow_kg_s"], activity_fraction),
        }
    return summary, columns


def main(arguments=None) -> int:
    from dada_solver.configuration import load_simulation_configuration
    from dada_solver.factory import build_initial_state, build_model, build_periodic_solver, initial_valve_topology
    from dada_solver.performance import calculate_cycle_performance
    from dada_solver.periodic import PeriodicStatus
    from dada_solver.topology import classify_cycle_topology
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("configuration", type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--activity-fraction", type=float, default=0.01,
                        help="Reporting cutoff relative to each port's forward peak; no physical threshold")
    options = parser.parse_args(arguments)
    config = load_simulation_configuration(options.configuration)
    model = build_model(config)
    result = build_periodic_solver(config, model).solve(build_initial_state(config, model), initial_valve_topology())
    if result.status is not PeriodicStatus.CONVERGED:
        print(f"periodic_status = {result.status.value}")
        return 2
    cycle = result.final_cycle
    summary, columns = extract_exchanger_duty(cycle, model, options.activity_fraction)
    performance = calculate_cycle_performance(cycle, ThermodynamicState.from_array(cycle.states[:, 0]), model.signed_angular_speed)
    summary.update({"source_configuration": str(options.configuration), "periodic_status": result.status.value,
                    "thermal_efficiency": performance.thermal_efficiency, "motor_power_W": performance.motor_power,
                    "cycle_topology": classify_cycle_topology(cycle.events).classification.value,
                    "cycles_completed": len(result.history)})
    options.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    json_path = Path(str(options.output_prefix) + ".json")
    csv_path = Path(str(options.output_prefix) + ".csv")
    json_path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(zip(*columns.values(), strict=True))
    print(f"duty_summary = {json_path}\nduty_history = {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
