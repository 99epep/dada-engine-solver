"""Conservative four-control-volume thermodynamic dynamics."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from numpy.typing import NDArray

from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.geometry import InstantaneousVolumes, MachineVolumes
from dada_solver.heat_transfer import ReservoirHeatTransfer
from dada_solver.hydraulics import HydraulicFlowModel, FlowResult
from dada_solver.kinematics import VolumeKinematics
from dada_solver.state import ThermodynamicState
from dada_solver.valves import PassiveCheckValve, ValveState


@dataclass(frozen=True, slots=True)
class ValveTopology:
    """Current hydraulic state of both independently passive valves."""

    hot_to_small: ValveState
    cold_to_large: ValveState


@dataclass(frozen=True, slots=True)
class NetworkFlows:
    """Signed link flows in their namesake directions, in kg/s."""

    large_to_hot: float
    small_to_cold: float
    hot_to_small: float
    cold_to_large: float
    large_to_hot_choked: bool
    small_to_cold_choked: bool
    hot_to_small_choked: bool
    cold_to_large_choked: bool


@dataclass(frozen=True, slots=True)
class ModelRates:
    """State derivative and diagnostic rates at one instant."""

    state_derivative: NDArray[np.float64]
    flows: NetworkFlows
    cold_heat_rate: float
    hot_heat_rate: float
    gas_work_rate: float
    mass_residual_rate: float
    energy_residual_rate: float


@dataclass(frozen=True, slots=True)
class ThermodynamicModel:
    """First-level conservative model with four independent pressures."""

    gas: CaloricallyPerfectGas
    machine_volumes: MachineVolumes
    kinematics: VolumeKinematics
    angular_speed: float
    cold_heat_transfer: ReservoirHeatTransfer
    hot_heat_transfer: ReservoirHeatTransfer
    large_hot_link: HydraulicFlowModel
    small_cold_link: HydraulicFlowModel
    hot_small_valve: PassiveCheckValve
    cold_large_valve: PassiveCheckValve
    continuous_ideal_diodes: bool = False
    study_crank_direction: int = 1

    @property
    def signed_angular_speed(self) -> float:
        """Study speed; integration uses increasing progress and its magnitude."""
        return self.study_crank_direction * self.angular_speed

    def __post_init__(self) -> None:
        if not math.isfinite(self.angular_speed) or self.angular_speed <= 0.0:
            raise ValueError("Angular speed must be finite and positive.")
        if self.study_crank_direction not in (-1, 1):
            raise ValueError("Study crank direction must be -1 or 1.")

    def volumes(self, theta: float) -> InstantaneousVolumes:
        """Return all four control-volume values at an angular position."""

        return InstantaneousVolumes(
            small_cylinder=self.kinematics.small_cylinder_volume(theta),
            large_cylinder=self.kinematics.large_cylinder_volume(theta),
            cold_heat_exchanger=self.machine_volumes.cold_heat_exchanger,
            hot_heat_exchanger=self.machine_volumes.hot_heat_exchanger,
        )

    def cylinder_volume_rates(self, theta: float) -> tuple[float, float]:
        """Return signed dV/dt values for S and L in m^3/s."""

        return (
            self.angular_speed
            * self.kinematics.small_cylinder_volume_derivative(theta),
            self.angular_speed
            * self.kinematics.large_cylinder_volume_derivative(theta),
        )

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
        effective_topology = self.effective_topology(theta, state, topology)
        hot_small = self.hot_small_valve.flow(
            effective_topology.hot_to_small,
            pressures[3],
            pressures[0],
            temperatures[3],
            self.gas,
        )
        cold_large = self.cold_large_valve.flow(
            effective_topology.cold_to_large,
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

    def valve_transitions(
        self,
        theta: float,
        state: ThermodynamicState,
        topology: ValveTopology,
    ) -> ValveTopology:
        """Apply both valves independently, including simultaneous transitions."""

        pressures = state.pressures(self.gas, self.volumes(theta))
        hot_small = self.hot_small_valve.updated_state(
            topology.hot_to_small, pressures[3], pressures[0]
        )
        cold_large = self.cold_large_valve.updated_state(
            topology.cold_to_large, pressures[2], pressures[1]
        )
        return ValveTopology(hot_small.current_state, cold_large.current_state)

    def _apply_bidirectional_link(
        self,
        mass_rates: NDArray[np.float64],
        energy_rates: NDArray[np.float64],
        first: int,
        second: int,
        flow: FlowResult,
        temperatures: NDArray[np.float64],
    ) -> None:
        signed_flow = flow.mass_flow_rate
        if signed_flow >= 0.0:
            directed = FlowResult(signed_flow, flow.is_choked)
            self._apply_directed_link(
                mass_rates,
                energy_rates,
                first,
                second,
                directed,
                temperatures[first],
            )
        else:
            directed = FlowResult(-signed_flow, flow.is_choked)
            self._apply_directed_link(
                mass_rates,
                energy_rates,
                second,
                first,
                directed,
                temperatures[second],
            )

    def _apply_directed_link(
        self,
        mass_rates: NDArray[np.float64],
        energy_rates: NDArray[np.float64],
        source: int,
        destination: int,
        flow: FlowResult,
        upstream_temperature: float,
    ) -> None:
        mass_flow = flow.mass_flow_rate
        enthalpy_flow = mass_flow * self.gas.enthalpy(upstream_temperature)
        mass_rates[source] -= mass_flow
        mass_rates[destination] += mass_flow
        energy_rates[source] -= enthalpy_flow
        energy_rates[destination] += enthalpy_flow
