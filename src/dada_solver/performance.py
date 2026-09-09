"""Cycle-integrated thermodynamic performance and conservation checks."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math

from dada_solver.integration import CycleIntegrationResult
from dada_solver.state import ThermodynamicState


class OperatingMode(Enum):
    REFRIGERATION = "refrigeration"
    MOTOR = "motor"
    NON_REFRIGERATION = "non_refrigeration"


@dataclass(frozen=True, slots=True)
class ConservationReport:
    absolute_mass_residual: float
    relative_mass_residual: float
    absolute_energy_residual: float
    relative_energy_residual: float


@dataclass(frozen=True, slots=True)
class CyclePerformance:
    """Signed branch quantities; legacy cold/hot names mean H_i/H_o.

    Those branches stay fixed when the external reservoirs are exchanged.
    In motor operation their powers must not be interpreted as refrigeration.
    """
    cold_heat_per_cycle: float
    hot_heat_per_cycle: float
    gas_work_per_cycle: float
    mechanical_input_per_cycle: float
    cooling_cop: float | None
    heating_cop: float | None
    cooling_power: float
    heating_power: float
    mechanical_input_power: float
    operating_mode: OperatingMode
    conservation: ConservationReport

    @property
    def heat_in_per_cycle(self) -> float:
        return self.cold_heat_per_cycle

    @property
    def heat_out_per_cycle(self) -> float:
        """Heat received by H_o (negative for intended heat rejection)."""
        return self.hot_heat_per_cycle

    @property
    def heat_in_power(self) -> float:
        return self.cooling_power

    @property
    def heat_out_power(self) -> float:
        """Signed heat received by H_o."""
        return -self.heating_power

    @property
    def gas_power(self) -> float:
        return -self.mechanical_input_power

    @property
    def thermal_efficiency(self) -> float | None:
        if self.operating_mode is not OperatingMode.MOTOR:
            return None
        return self.gas_work_per_cycle / self.heat_in_per_cycle

    @property
    def motor_power(self) -> float | None:
        return self.gas_power if self.operating_mode is OperatingMode.MOTOR else None


def calculate_cycle_performance(
    cycle: CycleIntegrationResult,
    initial_state: ThermodynamicState,
    angular_speed: float,
) -> CyclePerformance:
    """Calculate signed quantities using the study's signed crank speed.

    Negative speed denotes the motor reservoir assignment. Merely consuming
    work with that assignment does not constitute refrigeration.
    """

    if not cycle.completed:
        raise ValueError("Performance requires a completed cycle integration.")
    if angular_speed == 0.0 or not math.isfinite(angular_speed):
        raise ValueError("Angular speed must be finite and non-zero.")

    cold_heat = float(cycle.cold_heat[-1] - cycle.cold_heat[0])
    hot_heat = float(cycle.hot_heat[-1] - cycle.hot_heat[0])
    gas_work = float(cycle.gas_work[-1] - cycle.gas_work[0])
    mechanical_input = -gas_work
    frequency = abs(angular_speed) / (2.0 * math.pi)
    refrigeration = angular_speed > 0.0 and cold_heat > 0.0 and hot_heat < 0.0 and gas_work < 0.0
    motor = angular_speed < 0.0 and cold_heat > 0.0 and hot_heat < 0.0 and gas_work > 0.0
    cooling_cop = (
        cold_heat / mechanical_input
        if angular_speed > 0.0 and mechanical_input > 0.0 and cold_heat > 0.0
        else None
    )
    heating_cop = (
        -hot_heat / mechanical_input
        if angular_speed > 0.0 and mechanical_input > 0.0 and hot_heat < 0.0
        else None
    )

    final_state = cycle.final_state
    mass_residual = final_state.total_mass - initial_state.total_mass
    relative_mass_residual = mass_residual / initial_state.total_mass
    internal_energy_change = (
        final_state.total_internal_energy - initial_state.total_internal_energy
    )
    energy_residual = internal_energy_change - cold_heat - hot_heat + gas_work
    energy_scale = max(
        abs(initial_state.total_internal_energy),
        abs(cold_heat) + abs(hot_heat) + abs(gas_work),
    )

    return CyclePerformance(
        cold_heat_per_cycle=cold_heat,
        hot_heat_per_cycle=hot_heat,
        gas_work_per_cycle=gas_work,
        mechanical_input_per_cycle=mechanical_input,
        cooling_cop=cooling_cop,
        heating_cop=heating_cop,
        cooling_power=cold_heat * frequency,
        heating_power=-hot_heat * frequency,
        mechanical_input_power=mechanical_input * frequency,
        operating_mode=(
            OperatingMode.REFRIGERATION
            if refrigeration
            else OperatingMode.MOTOR if motor else OperatingMode.NON_REFRIGERATION
        ),
        conservation=ConservationReport(
            absolute_mass_residual=mass_residual,
            relative_mass_residual=relative_mass_residual,
            absolute_energy_residual=energy_residual,
            relative_energy_residual=energy_residual / energy_scale,
        ),
    )
