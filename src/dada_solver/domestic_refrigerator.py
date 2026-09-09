"""External domestic-refrigerator comparison boundaries.

This module does not insert cabinet heat leaks into the DADA cycle. It keeps
standard test temperatures, cabinet loads and auxiliary electricity separate
so machine COP is not confused with appliance energy consumption.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True, slots=True)
class RefrigeratorTemperaturePoint:
    """One ambient-to-compartment temperature lift in kelvin."""

    name: str
    ambient_temperature: float
    compartment_temperature: float

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("Temperature-point name must not be empty.")
        for name, value in (
            ("Ambient temperature", self.ambient_temperature),
            ("Compartment temperature", self.compartment_temperature),
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")
        if self.ambient_temperature <= self.compartment_temperature:
            raise ValueError("Ambient temperature must exceed compartment temperature.")

    @property
    def temperature_lift(self) -> float:
        return self.ambient_temperature - self.compartment_temperature

    @property
    def carnot_cooling_cop(self) -> float:
        """Return the reversible bound at the declared boundary temperatures."""

        return self.compartment_temperature / self.temperature_lift


@dataclass(frozen=True, slots=True)
class SteadyCabinetLoad:
    """Simple external cabinet load, independent of refrigerating technology."""

    thermal_conductance: float
    internal_heat_load: float = 0.0

    def __post_init__(self) -> None:
        if not math.isfinite(self.thermal_conductance) or self.thermal_conductance <= 0.0:
            raise ValueError("Cabinet thermal conductance must be finite and positive.")
        if not math.isfinite(self.internal_heat_load) or self.internal_heat_load < 0.0:
            raise ValueError("Internal heat load must be finite and non-negative.")

    def cooling_power(self, point: RefrigeratorTemperaturePoint) -> float:
        return self.thermal_conductance * point.temperature_lift + self.internal_heat_load


@dataclass(frozen=True, slots=True)
class AppliancePowerAssessment:
    cabinet_cooling_power: float
    thermodynamic_input_power: float
    shaft_input_power: float
    total_appliance_input_power: float


@dataclass(frozen=True, slots=True)
class AnnualEnergyBudget:
    """Declared appliance electricity budget over one conventional year."""

    annual_energy_kwh: float
    hours_per_year: float = 8760.0

    def __post_init__(self) -> None:
        for name, value in (
            ("Annual energy", self.annual_energy_kwh),
            ("Hours per year", self.hours_per_year),
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")

    @property
    def average_input_power(self) -> float:
        """Return the continuous-equivalent electrical power in watts."""

        return 1000.0 * self.annual_energy_kwh / self.hours_per_year

    def maximum_on_input_power(
        self,
        duty_fraction: float,
        off_cycle_power: float = 0.0,
    ) -> float:
        """Return the allowed total input while running for a declared duty cycle."""

        _validate_duty_fraction(duty_fraction)
        if not math.isfinite(off_cycle_power) or off_cycle_power < 0.0:
            raise ValueError("Off-cycle power must be finite and non-negative.")
        available = self.average_input_power - (1.0 - duty_fraction) * off_cycle_power
        return available / duty_fraction


@dataclass(frozen=True, slots=True)
class CyclingApplianceAssessment:
    """Annualized result for an appliance that cycles between on and off states."""

    duty_fraction: float
    on_cycle_input_power: float
    off_cycle_input_power: float
    average_input_power: float
    annual_energy_kwh: float


def assess_appliance_power(
    point: RefrigeratorTemperaturePoint,
    cabinet: SteadyCabinetLoad,
    machine_cooling_cop: float,
    transmission_efficiency: float,
    auxiliary_power: float,
) -> AppliancePowerAssessment:
    """Propagate an explicit cabinet load through machine and appliance boundaries."""

    for name, value in (
        ("Machine cooling COP", machine_cooling_cop),
        ("Transmission efficiency", transmission_efficiency),
    ):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive.")
    if transmission_efficiency > 1.0:
        raise ValueError("Transmission efficiency must not exceed one.")
    if not math.isfinite(auxiliary_power) or auxiliary_power < 0.0:
        raise ValueError("Auxiliary power must be finite and non-negative.")
    cooling = cabinet.cooling_power(point)
    thermodynamic = cooling / machine_cooling_cop
    shaft = thermodynamic / transmission_efficiency
    return AppliancePowerAssessment(
        cabinet_cooling_power=cooling,
        thermodynamic_input_power=thermodynamic,
        shaft_input_power=shaft,
        total_appliance_input_power=shaft + auxiliary_power,
    )


def annualize_cycling_appliance(
    on_cycle: AppliancePowerAssessment,
    duty_fraction: float,
    off_cycle_input_power: float = 0.0,
    hours_per_year: float = 8760.0,
) -> CyclingApplianceAssessment:
    """Convert explicit on/off powers and duty fraction to annual electricity."""

    _validate_duty_fraction(duty_fraction)
    if not math.isfinite(off_cycle_input_power) or off_cycle_input_power < 0.0:
        raise ValueError("Off-cycle input power must be finite and non-negative.")
    if not math.isfinite(hours_per_year) or hours_per_year <= 0.0:
        raise ValueError("Hours per year must be finite and positive.")
    average = (
        duty_fraction * on_cycle.total_appliance_input_power
        + (1.0 - duty_fraction) * off_cycle_input_power
    )
    return CyclingApplianceAssessment(
        duty_fraction=duty_fraction,
        on_cycle_input_power=on_cycle.total_appliance_input_power,
        off_cycle_input_power=off_cycle_input_power,
        average_input_power=average,
        annual_energy_kwh=average * hours_per_year / 1000.0,
    )


def _validate_duty_fraction(value: float) -> None:
    if not math.isfinite(value) or not 0.0 < value <= 1.0:
        raise ValueError("Duty fraction must lie in (0, 1].")
