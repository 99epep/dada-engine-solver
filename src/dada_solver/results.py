"""Scientific diagnostic extraction from a completed cycle."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from dada_solver.dynamics import NetworkFlows, ThermodynamicModel
from dada_solver.integration import CycleIntegrationResult
from dada_solver.topology import CycleTopologyDiagnostic, classify_cycle_topology


@dataclass(frozen=True, slots=True)
class Extrema:
    minimum: float
    maximum: float


@dataclass(frozen=True, slots=True)
class ValveEventDiagnostic:
    angle_radians: float
    angle_degrees: float
    valve_name: str
    transition: str
    small_cylinder_lambda: float
    large_cylinder_lambda: float


@dataclass(frozen=True, slots=True)
class CycleDiagnostics:
    pressure_extrema: dict[str, Extrema]
    temperature_extrema: dict[str, Extrema]
    mass_flow_extrema: dict[str, Extrema]
    cold_heat_rate_extrema: Extrema
    """Heat rate into working gas C/H_i from its exchanger boundary."""
    hot_heat_rate_extrema: Extrema
    """Heat rate into working gas H/H_o from its exchanger boundary."""
    valve_events: tuple[ValveEventDiagnostic, ...]
    topology: CycleTopologyDiagnostic
    orifice_pressure_regularization: float = 0.0
    regularization_sample_fraction: float = 0.0


def extract_cycle_diagnostics(
    cycle: CycleIntegrationResult,
    model: ThermodynamicModel,
    heat_rate_provider=None,
) -> CycleDiagnostics:
    """Re-evaluate sampled states to obtain diagnostic histories and extrema."""

    if not cycle.completed:
        raise ValueError("Diagnostics require a completed cycle integration.")
    sample_count = cycle.angles.size
    pressures = np.empty((4, sample_count), dtype=float)
    temperatures = np.empty((4, sample_count), dtype=float)
    flow_values = np.empty((4, sample_count), dtype=float)
    cold_heat_rates = np.empty(sample_count, dtype=float)
    hot_heat_rates = np.empty(sample_count, dtype=float)
    regularization_hits = np.zeros(sample_count, dtype=bool)
    regularization_width = max(
        float(getattr(flow_model, "pressure_regularization", 0.0))
        for flow_model in (
            model.large_hot_link,
            model.small_cold_link,
            model.hot_small_valve.flow_model,
            model.cold_large_valve.flow_model,
        )
    )

    from dada_solver.state import ThermodynamicState

    for index, (angle, topology) in enumerate(
        zip(cycle.angles, cycle.topologies, strict=True)
    ):
        state = ThermodynamicState.from_array(cycle.states[:, index])
        temperatures[:, index] = state.temperatures(model.gas)
        pressures[:, index] = state.pressures(model.gas, model.volumes(float(angle)))
        if regularization_width > 0.0:
            pressure_differences = (
                pressures[1, index] - pressures[3, index],
                pressures[0, index] - pressures[2, index],
                pressures[3, index] - pressures[0, index],
                pressures[2, index] - pressures[1, index],
            )
            regularization_hits[index] = any(
                abs(value) <= regularization_width
                for value in pressure_differences
            )
        rates = model.evaluate(float(angle), state, topology)
        flow_values[:, index] = _flow_array(rates.flows)
        if heat_rate_provider is None:
            cold_heat_rates[index] = rates.cold_heat_rate
            hot_heat_rates[index] = rates.hot_heat_rate
        else:
            cold_heat_rates[index], hot_heat_rates[index] = heat_rate_provider(
                index, float(angle), state)

    names = ("S", "L", "C", "H")
    flow_names = ("large_to_hot", "small_to_cold", "hot_to_small", "cold_to_large")
    event_diagnostics = []
    for event in cycle.events:
        volumes = model.volumes(event.angle)
        event_diagnostics.append(
            ValveEventDiagnostic(
                angle_radians=event.angle,
                angle_degrees=float(np.degrees(event.angle)),
                valve_name=event.valve_name,
                transition=event.current_state.value,
                small_cylinder_lambda=(
                    model.machine_volumes.small_cylinder.closure_fraction(
                        volumes.small_cylinder
                    )
                ),
                large_cylinder_lambda=(
                    model.machine_volumes.large_cylinder.closure_fraction(
                        volumes.large_cylinder
                    )
                ),
            )
        )

    return CycleDiagnostics(
        pressure_extrema={
            name: Extrema(float(np.min(values)), float(np.max(values)))
            for name, values in zip(names, pressures, strict=True)
        },
        temperature_extrema={
            name: Extrema(float(np.min(values)), float(np.max(values)))
            for name, values in zip(names, temperatures, strict=True)
        },
        mass_flow_extrema={
            name: Extrema(float(np.min(values)), float(np.max(values)))
            for name, values in zip(flow_names, flow_values, strict=True)
        },
        cold_heat_rate_extrema=Extrema(
            float(np.min(cold_heat_rates)), float(np.max(cold_heat_rates))
        ),
        hot_heat_rate_extrema=Extrema(
            float(np.min(hot_heat_rates)), float(np.max(hot_heat_rates))
        ),
        valve_events=tuple(event_diagnostics),
        topology=classify_cycle_topology(cycle.events),
        orifice_pressure_regularization=regularization_width,
        regularization_sample_fraction=float(np.mean(regularization_hits)),
    )


def _flow_array(flows: NetworkFlows) -> NDArray[np.float64]:
    return np.array(
        [
            flows.large_to_hot,
            flows.small_to_cold,
            flows.hot_to_small,
            flows.cold_to_large,
        ]
    )
