"""External lumped cooling loads, independent of the DADA cycle model."""

from __future__ import annotations

from dataclasses import dataclass
import math

from dada_solver.performance import CyclePerformance, OperatingMode


@dataclass(frozen=True, slots=True)
class LumpedPhaseChangeLoad:
    """Sensible cooling followed by an isothermal phase change."""

    material_name: str
    mass: float
    initial_temperature: float
    target_temperature: float
    initial_phase_specific_heat: float
    phase_change_specific_energy: float
    phase_change_fraction: float = 1.0

    def __post_init__(self) -> None:
        if not self.material_name.strip():
            raise ValueError("Material name must not be empty.")
        for name, value in (
            ("mass", self.mass),
            ("initial temperature", self.initial_temperature),
            ("target temperature", self.target_temperature),
            ("specific heat", self.initial_phase_specific_heat),
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"Load {name} must be finite and positive.")
        if self.initial_temperature < self.target_temperature:
            raise ValueError("Initial load temperature must not be below target temperature.")
        if (
            not math.isfinite(self.phase_change_specific_energy)
            or self.phase_change_specific_energy < 0.0
        ):
            raise ValueError("Phase-change specific energy must be non-negative.")
        if (
            not math.isfinite(self.phase_change_fraction)
            or not 0.0 <= self.phase_change_fraction <= 1.0
        ):
            raise ValueError("Phase-change fraction must lie in [0, 1].")

    @property
    def sensible_energy(self) -> float:
        return (
            self.mass
            * self.initial_phase_specific_heat
            * (self.initial_temperature - self.target_temperature)
        )

    @property
    def phase_change_energy(self) -> float:
        return (
            self.mass
            * self.phase_change_specific_energy
            * self.phase_change_fraction
        )

    @property
    def ideal_cooling_energy(self) -> float:
        """Return the ideal load only, excluding container and parasitic gains."""

        return self.sensible_energy + self.phase_change_energy


@dataclass(frozen=True, slots=True)
class CoolingTask:
    load: LumpedPhaseChangeLoad
    allowed_time: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.allowed_time) or self.allowed_time <= 0.0:
            raise ValueError("Allowed cooling time must be finite and positive.")

    @property
    def required_average_cooling_power(self) -> float:
        return self.load.ideal_cooling_energy / self.allowed_time

    def required_cop(self, mechanical_input_power: float) -> float:
        if not math.isfinite(mechanical_input_power) or mechanical_input_power <= 0.0:
            raise ValueError("Mechanical input power must be finite and positive.")
        return self.required_average_cooling_power / mechanical_input_power


@dataclass(frozen=True, slots=True)
class HumanPowerEnvelope:
    minimum_power: float
    nominal_power: float
    maximum_power: float

    def __post_init__(self) -> None:
        values = self.minimum_power, self.nominal_power, self.maximum_power
        if not all(math.isfinite(value) and value > 0.0 for value in values):
            raise ValueError("Human power values must be finite and positive.")
        if not self.minimum_power <= self.nominal_power <= self.maximum_power:
            raise ValueError("Nominal human power must lie within its range.")


@dataclass(frozen=True, slots=True)
class CoolingTaskAssessment:
    ideal_load_energy: float
    cooling_energy_per_cycle: float
    average_cooling_power: float
    required_mechanical_input_power: float
    cooling_cop: float | None
    ideal_estimated_completion_time: float | None
    required_pedaling_power_for_target_time: float | None
    operable_at_nominal_human_power: bool
    target_met_at_operating_point: bool


def assess_cooling_task(
    performance: CyclePerformance,
    task: CoolingTask,
    human_power: HumanPowerEnvelope,
) -> CoolingTaskAssessment:
    """Compare one fixed machine operating point with an external ideal load."""

    cooling_available = (
        performance.operating_mode is OperatingMode.REFRIGERATION
        and performance.cooling_power > 0.0
    )
    completion_time = (
        task.load.ideal_cooling_energy / performance.cooling_power
        if cooling_available
        else None
    )
    cop = performance.cooling_cop
    target_pedaling_power = (
        task.required_average_cooling_power / cop
        if cop is not None and cop > 0.0
        else None
    )
    operable = (
        cooling_available
        and performance.mechanical_input_power > 0.0
        and performance.mechanical_input_power <= human_power.nominal_power
    )
    target_met = (
        operable
        and completion_time is not None
        and completion_time <= task.allowed_time
    )
    return CoolingTaskAssessment(
        ideal_load_energy=task.load.ideal_cooling_energy,
        cooling_energy_per_cycle=performance.cold_heat_per_cycle,
        average_cooling_power=performance.cooling_power,
        required_mechanical_input_power=performance.mechanical_input_power,
        cooling_cop=cop,
        ideal_estimated_completion_time=completion_time,
        required_pedaling_power_for_target_time=target_pedaling_power,
        operable_at_nominal_human_power=operable,
        target_met_at_operating_point=target_met,
    )
