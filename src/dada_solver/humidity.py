"""Moist-air phase-change screening outside the dry-gas cycle model."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math

import numpy as np

from dada_solver.dynamics import ThermodynamicModel
from dada_solver.integration import CycleIntegrationResult
from dada_solver.state import ThermodynamicState


class MoistureScreeningVerdict(Enum):
    """Availability and phase-change risk of a dry-gas simulation."""

    CLEAR = "clear"
    CONDENSATION_OR_FROST_RISK = "condensation_or_frost_risk"
    UNAVAILABLE = "unavailable"


@dataclass(frozen=True, slots=True)
class HumidityScreeningConfiguration:
    """Humidity of the uniformly charged gas before any phase separation."""

    relative_humidity: float

    def __post_init__(self) -> None:
        if (
            not math.isfinite(self.relative_humidity)
            or not 0.0 <= self.relative_humidity <= 1.0
        ):
            raise ValueError("Relative humidity must lie in [0, 1].")


@dataclass(frozen=True, slots=True)
class ControlVolumeMoistureDiagnostic:
    """Maximum ideal-mixture saturation ratio in one control volume."""

    maximum_saturation_ratio: float
    angle_degrees_at_maximum: float
    temperature_at_maximum: float
    pressure_at_maximum: float
    first_saturated_angle_degrees: float | None
    first_saturated_temperature: float | None
    first_saturated_pressure: float | None
    first_predicted_phase: str | None


@dataclass(frozen=True, slots=True)
class MoistureScreeningReport:
    """Screening result; it does not model condensate mass or latent heat."""

    verdict: MoistureScreeningVerdict
    initial_water_mole_fraction: float | None
    control_volumes: dict[str, ControlVolumeMoistureDiagnostic]
    assumptions: tuple[str, ...]
    unavailable_reason: str | None = None


def saturation_vapor_pressure(temperature: float, *, over_ice: bool) -> float:
    """Return water saturation pressure in Pa using Murphy--Koop equations.

    The ice expression is used below the triple point for equilibrium frost
    screening. The liquid-water expression is used above it and for the
    initial relative-humidity definition at ordinary ambient temperature.
    """

    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("Temperature must be finite and positive.")
    if over_ice:
        if not 110.0 <= temperature <= 273.16:
            raise ValueError("Ice saturation correlation requires 110 <= T <= 273.16 K.")
        logarithm = (
            9.550426
            - 5723.265 / temperature
            + 3.53068 * math.log(temperature)
            - 0.00728332 * temperature
        )
    else:
        if not 123.0 <= temperature <= 332.0:
            raise ValueError("Liquid saturation correlation requires 123 <= T <= 332 K.")
        logarithm = (
            54.842763
            - 6763.22 / temperature
            - 4.210 * math.log(temperature)
            + 0.000367 * temperature
            + math.tanh(0.0415 * (temperature - 218.8))
            * (
                53.878
                - 1331.22 / temperature
                - 9.44523 * math.log(temperature)
                + 0.014025 * temperature
            )
        )
    return math.exp(logarithm)


def assess_moisture_phase_change_risk(
    cycle: CycleIntegrationResult,
    model: ThermodynamicModel,
    configuration: HumidityScreeningConfiguration | None,
    charge_temperature: float,
    charge_pressure: float,
) -> MoistureScreeningReport:
    """Screen a completed dry-gas cycle for condensation or frost onset.

    Before phase change, a closed, initially uniform ideal mixture retains one
    water mole fraction. If saturation is exceeded anywhere, the dry-gas model
    has crossed its validity boundary; no condensate removal or latent heat is
    calculated beyond that point.
    """

    assumptions = (
        "uniform initial water mole fraction",
        "ideal gas mixture before phase change",
        "no condensation, frost, drainage, blockage, or latent heat model",
        "equilibrium ice saturation below the water triple point",
    )
    if configuration is None:
        return MoistureScreeningReport(
            verdict=MoistureScreeningVerdict.UNAVAILABLE,
            initial_water_mole_fraction=None,
            control_volumes={},
            assumptions=assumptions,
            unavailable_reason="initial relative humidity was not configured",
        )
    if not cycle.completed:
        raise ValueError("Moisture screening requires a completed cycle.")
    initial_saturation_pressure = saturation_vapor_pressure(
        charge_temperature, over_ice=charge_temperature < 273.16
    )
    water_mole_fraction = (
        configuration.relative_humidity
        * initial_saturation_pressure
        / charge_pressure
    )
    if water_mole_fraction >= 1.0:
        raise ValueError("Initial humidity implies an invalid water mole fraction.")

    ratios = np.full((4, cycle.angles.size), np.nan, dtype=float)
    temperatures = np.empty_like(ratios)
    pressures = np.empty_like(ratios)
    unavailable_samples: list[str] = []
    for index, angle in enumerate(cycle.angles):
        state = ThermodynamicState.from_array(cycle.states[:, index])
        temperatures[:, index] = state.temperatures(model.gas)
        pressures[:, index] = state.pressures(model.gas, model.volumes(float(angle)))
        for volume_index in range(4):
            temperature = float(temperatures[volume_index, index])
            try:
                saturation_pressure = saturation_vapor_pressure(
                    temperature, over_ice=temperature < 273.16
                )
            except ValueError:
                unavailable_samples.append(
                    f"{('S', 'L', 'C', 'H')[volume_index]} at {temperature:.6g} K"
                )
                continue
            ratios[volume_index, index] = (
                water_mole_fraction
                * pressures[volume_index, index]
                / saturation_pressure
            )

    diagnostics: dict[str, ControlVolumeMoistureDiagnostic] = {}
    maximum_ratio = 0.0
    for volume_index, name in enumerate(("S", "L", "C", "H")):
        if np.all(np.isnan(ratios[volume_index])):
            continue
        sample_index = int(np.nanargmax(ratios[volume_index]))
        ratio = float(ratios[volume_index, sample_index])
        maximum_ratio = max(maximum_ratio, ratio)
        saturated_indices = np.flatnonzero(ratios[volume_index] >= 1.0)
        first_index = (
            None if saturated_indices.size == 0 else int(saturated_indices[0])
        )
        first_temperature = (
            None if first_index is None else float(temperatures[volume_index, first_index])
        )
        diagnostics[name] = ControlVolumeMoistureDiagnostic(
            maximum_saturation_ratio=ratio,
            angle_degrees_at_maximum=math.degrees(float(cycle.angles[sample_index])),
            temperature_at_maximum=float(temperatures[volume_index, sample_index]),
            pressure_at_maximum=float(pressures[volume_index, sample_index]),
            first_saturated_angle_degrees=(
                None
                if first_index is None
                else math.degrees(float(cycle.angles[first_index]))
            ),
            first_saturated_temperature=first_temperature,
            first_saturated_pressure=(
                None if first_index is None else float(pressures[volume_index, first_index])
            ),
            first_predicted_phase=(
                None
                if first_temperature is None
                else ("frost" if first_temperature < 273.16 else "liquid_condensation")
            ),
        )
    if maximum_ratio >= 1.0:
        verdict = MoistureScreeningVerdict.CONDENSATION_OR_FROST_RISK
        unavailable_reason = None
    elif unavailable_samples:
        verdict = MoistureScreeningVerdict.UNAVAILABLE
        unavailable_reason = (
            "sampled states outside the selected saturation-pressure "
            f"correlation range, including {unavailable_samples[0]}"
        )
    else:
        verdict = MoistureScreeningVerdict.CLEAR
        unavailable_reason = None
    return MoistureScreeningReport(
        verdict=verdict,
        initial_water_mole_fraction=water_mole_fraction,
        control_volumes=diagnostics,
        assumptions=assumptions,
        unavailable_reason=unavailable_reason,
    )


def moisture_saturation_ratio_history(
    cycle: CycleIntegrationResult,
    model: ThermodynamicModel,
    configuration: HumidityScreeningConfiguration,
    charge_temperature: float,
    charge_pressure: float,
) -> np.ndarray:
    """Return S/L/C/H saturation-ratio histories, with NaN when unavailable."""

    initial_pressure = saturation_vapor_pressure(
        charge_temperature, over_ice=charge_temperature < 273.16
    )
    water_mole_fraction = (
        configuration.relative_humidity * initial_pressure / charge_pressure
    )
    ratios = np.full((4, cycle.angles.size), np.nan, dtype=float)
    for sample_index, angle in enumerate(cycle.angles):
        state = ThermodynamicState.from_array(cycle.states[:, sample_index])
        temperatures = state.temperatures(model.gas)
        pressures = state.pressures(model.gas, model.volumes(float(angle)))
        for volume_index, (temperature, pressure) in enumerate(
            zip(temperatures, pressures, strict=True)
        ):
            try:
                saturation_pressure = saturation_vapor_pressure(
                    float(temperature), over_ice=float(temperature) < 273.16
                )
            except ValueError:
                continue
            ratios[volume_index, sample_index] = (
                water_mole_fraction * float(pressure) / saturation_pressure
            )
    return ratios
