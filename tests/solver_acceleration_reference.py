"""Frozen pre-Stage-2 reference algebra, used only by equivalence tests.

Copied before refactoring on Debian 13. Do not update this oracle when changing
production code. There is only one conservative assembly in production.
"""
from dataclasses import replace
import math
import numpy as np
from dada_solver.dynamics import ModelRates, NetworkFlows, ValveTopology
from dada_solver.geometry import InstantaneousVolumes
from dada_solver.state import ThermodynamicState
from dada_solver.heat_transfer import PrescribedHeatRate
from dada_solver.valves import ValveState

def evaluate(
    self,
    theta: float,
    state: ThermodynamicState,
    topology: ValveTopology,
) -> ModelRates:
    """Evaluate conservative rates without imposing a nominal phase sequence."""

    combined_provider = getattr(
        self.kinematics, "cylinder_volumes_and_derivatives", None
    )
    if combined_provider is None:
        volumes = self.volumes(theta)
        volume_rate_small, volume_rate_large = self.cylinder_volume_rates(theta)
    else:
        small, large, small_derivative, large_derivative = combined_provider(theta)
        volumes = InstantaneousVolumes(
            small,
            large,
            self.machine_volumes.cold_heat_exchanger,
            self.machine_volumes.hot_heat_exchanger,
        )
        volume_rate_small = self.angular_speed * small_derivative
        volume_rate_large = self.angular_speed * large_derivative
    temperatures, pressures = state.temperatures_and_pressures(self.gas, volumes)

    large_hot = self.large_hot_link.bidirectional_flow(
        pressures[1], pressures[3], temperatures[1], temperatures[3], self.gas
    )
    small_cold = self.small_cold_link.bidirectional_flow(
        pressures[0], pressures[2], temperatures[0], temperatures[2], self.gas
    )
    selected_topology = effective_topology(self, theta, state, topology)
    hot_small = self.hot_small_valve.flow(
        selected_topology.hot_to_small,
        pressures[3],
        pressures[0],
        temperatures[3],
        self.gas,
    )
    cold_large = self.cold_large_valve.flow(
        selected_topology.cold_to_large,
        pressures[2],
        pressures[1],
        temperatures[2],
        self.gas,
    )

    mass_rates = np.zeros(4, dtype=float)
    energy_rates = np.zeros(4, dtype=float)
    self._apply_bidirectional_link(
        mass_rates, energy_rates, 1, 3, large_hot, temperatures
    )
    self._apply_bidirectional_link(
        mass_rates, energy_rates, 0, 2, small_cold, temperatures
    )
    self._apply_directed_link(
        mass_rates, energy_rates, 3, 0, hot_small, temperatures[3]
    )
    self._apply_directed_link(
        mass_rates, energy_rates, 2, 1, cold_large, temperatures[2]
    )

    cold_heat_rate = self.cold_heat_transfer.heat_rate(temperatures[2])
    hot_heat_rate = self.hot_heat_transfer.heat_rate(temperatures[3])
    energy_rates[2] += cold_heat_rate
    energy_rates[3] += hot_heat_rate

    energy_rates[0] -= pressures[0] * volume_rate_small
    energy_rates[1] -= pressures[1] * volume_rate_large
    gas_work_rate = (
        pressures[0] * volume_rate_small + pressures[1] * volume_rate_large
    )

    derivative = np.empty(8, dtype=float)
    derivative[0::2] = mass_rates
    derivative[1::2] = energy_rates
    mass_residual_rate = float(np.sum(mass_rates))
    energy_residual_rate = float(
        np.sum(energy_rates)
        - cold_heat_rate
        - hot_heat_rate
        + gas_work_rate
    )

    return ModelRates(
        state_derivative=derivative,
        flows=NetworkFlows(
            large_to_hot=large_hot.mass_flow_rate,
            small_to_cold=small_cold.mass_flow_rate,
            hot_to_small=hot_small.mass_flow_rate,
            cold_to_large=cold_large.mass_flow_rate,
            large_to_hot_choked=large_hot.is_choked,
            small_to_cold_choked=small_cold.is_choked,
            hot_to_small_choked=hot_small.is_choked,
            cold_to_large_choked=cold_large.is_choked,
        ),
        cold_heat_rate=cold_heat_rate,
        hot_heat_rate=hot_heat_rate,
        gas_work_rate=gas_work_rate,
        mass_residual_rate=mass_residual_rate,
        energy_residual_rate=energy_residual_rate,
    )

def effective_topology(
    self,
    theta: float,
    state: ThermodynamicState,
    stored_topology: ValveTopology,
) -> ValveTopology:
    """Return instantaneous diode states or the stored hysteretic states."""

    if not self.continuous_ideal_diodes:
        return stored_topology
    pressures = state.pressures(self.gas, self.volumes(theta))
    return ValveTopology(
        ValveState.OPEN if pressures[3] > pressures[0] else ValveState.CLOSED,
        ValveState.OPEN if pressures[2] > pressures[1] else ValveState.CLOSED,
    )

def flow_contexts(self, angle, values):
    """Geometry-independent pressure/flow state for contextual wall closures."""
    gas = ThermodynamicState.from_array(np.asarray(values[:8]))
    pressures = gas.pressures(self.model.gas,self.model.volumes(angle))
    flows = evaluate(self.model,angle,gas,ValveTopology(ValveState.CLOSED,ValveState.CLOSED)).flows
    frequency = self.model.angular_speed/(2*math.pi)
    contexts = (
        dict(frequency_hz=frequency,passages=((flows.small_to_cold,pressures[0],pressures[2]),
                                             (flows.cold_to_large,pressures[2],pressures[1]))),
        dict(frequency_hz=frequency,passages=((flows.large_to_hot,pressures[1],pressures[3]),
                                             (flows.hot_to_small,pressures[3],pressures[0]))))
    for index,context in zip((2,3),contexts):
        context['port_pressure_ratios'] = tuple(max(p1,p2)/min(p1,p2) for _,p1,p2 in context['passages'])
        # A closed valve supports a pressure jump; it is not an axial tube
        # pressure gradient at zero flow. Retain the port ratio separately.
        context['passages'] = tuple((flow,p1,p2) if flow else
            (flow,pressures[index],pressures[index]) for flow,p1,p2 in context['passages'])
    return contexts

def thermal_rates(self, angle, values):
    """Single diagnostic/integration access path; never substitutes static UA."""
    gas = ThermodynamicState.from_array(np.asarray(values[:8]))
    temperatures = gas.temperatures(self.model.gas)
    contextual = any(getattr(x,'requires_flow_context',False) for x in (self.heat_in,self.heat_out))
    contexts = flow_contexts(self,angle,values) if contextual else (None,None)
    return tuple(exchanger.rates(temperatures[index],values[index+6],context=context)
                 if getattr(exchanger,'requires_flow_context',False)
                 else exchanger.rates(temperatures[index],values[index+6])
                 for index,exchanger,context in zip((2,3),(self.heat_in,self.heat_out),contexts))

def derivative(self, angle, values):
    gas = ThermodynamicState.from_array(np.asarray(values[:8]))
    incoming, outgoing = thermal_rates(self,angle, values)
    # Reuse all existing conservative transport and passive-valve equations.
    instantaneous = replace(self.model,
        cold_heat_transfer=PrescribedHeatRate(incoming['gas_heat_w']),
        hot_heat_transfer=PrescribedHeatRate(outgoing['gas_heat_w']))
    rates = evaluate(instantaneous, angle, gas, ValveTopology(ValveState.CLOSED, ValveState.CLOSED))
    return np.r_[rates.state_derivative,
                 incoming['wall_energy_rate_w'], outgoing['wall_energy_rate_w'],
                 incoming['air_heat_w'], outgoing['air_heat_w'],
                 incoming['gas_heat_w'], outgoing['gas_heat_w'], rates.gas_work_rate] / self.model.angular_speed
