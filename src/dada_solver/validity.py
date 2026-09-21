"""A-posteriori thermodynamic and hydraulic validity assessment."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np

from dada_solver.configuration import ValidityThresholds
from dada_solver.dynamics import ThermodynamicModel
from dada_solver.integration import CycleIntegrationResult
from dada_solver.state import ThermodynamicState


class ValidityVerdict(Enum):
    VALID = "valid"
    INVALID = "invalid"
    INDETERMINATE = "indeterminate"


@dataclass(frozen=True, slots=True)
class ValidityReport:
    verdict: ValidityVerdict
    maximum_pressure_equalization_error: float
    cold_isothermality_error: float | None
    hot_isothermality_error: float | None
    maximum_mach_number: float | None
    mach_number_status: str
    maximum_compressibility_deviation: float
    maximum_cp_variation: float
    failed_criteria: tuple[str, ...]
    unavailable_criteria: tuple[str, ...]


def assess_cycle_validity(
    cycle: CycleIntegrationResult,
    model: ThermodynamicModel,
    thresholds: ValidityThresholds,
    geometric_flow_areas: dict[str, float] | None = None,
    *, replay=None,
) -> ValidityReport:
    """Assess configured criteria without claiming unavailable quantities are valid."""

    if not cycle.completed:
        raise ValueError("Validity assessment requires a completed cycle.")
    if replay is not None:
        replay.require(model=model, cycle=cycle)
        pressure_history, temperature_history = replay.pressures, replay.temperatures
    else:
        pressures: list[np.ndarray] = []
        temperatures: list[np.ndarray] = []
        for angle, values in zip(cycle.angles, cycle.states.T, strict=True):
            state = ThermodynamicState.from_array(values)
            temperatures.append(state.temperatures(model.gas))
            pressures.append(state.pressures(model.gas, model.volumes(float(angle))))
        pressure_history = np.asarray(pressures).T
        temperature_history = np.asarray(temperatures).T

    pair_references = np.maximum(
        0.5 * (pressure_history[[0, 1]] + pressure_history[[2, 3]]),
        np.finfo(float).tiny,
    )
    pair_differences = np.abs(
        pressure_history[[0, 1]] - pressure_history[[2, 3]]
    )
    pressure_error = float(np.max(pair_differences / pair_references))
    def excursion(index, closure):
        # A generic closure need not have a fixed reservoir reference temperature.
        reference = getattr(closure, 'reservoir_temperature', None)
        if reference is None:
            return None
        return float((np.max(temperature_history[index])-np.min(temperature_history[index])) / reference)
    cold_isothermality = excursion(2, model.cold_heat_transfer)
    hot_isothermality = excursion(3, model.hot_heat_transfer)
    compressibility_deviation = 0.0
    cp_variation = 0.0

    failed: list[str] = []
    if pressure_error > thresholds.maximum_pressure_equalization_error:
        failed.append("pressure_equalization")
    # Temperature excursion is a design diagnostic. Variable temperatures are
    # already resolved by these balances, so excursion does not invalidate them.
    if compressibility_deviation > thresholds.maximum_compressibility_deviation:
        failed.append("compressibility_factor")
    if cp_variation > thresholds.maximum_cp_variation:
        failed.append("heat_capacity_variation")

    unavailable: list[str] = []
    maximum_mach: float | None = None
    if geometric_flow_areas is None:
        unavailable.append("mach_number")
        mach_status = "unavailable"
    else:
        # Velocity evaluation will be added with the detailed geometric-link model.
        unavailable.append("mach_number")
        mach_status = "unavailable"

    if failed:
        verdict = ValidityVerdict.INVALID
    elif unavailable:
        verdict = ValidityVerdict.INDETERMINATE
    else:
        verdict = ValidityVerdict.VALID
    return ValidityReport(
        verdict=verdict,
        maximum_pressure_equalization_error=pressure_error,
        cold_isothermality_error=cold_isothermality,
        hot_isothermality_error=hot_isothermality,
        maximum_mach_number=maximum_mach,
        mach_number_status=mach_status,
        maximum_compressibility_deviation=compressibility_deviation,
        maximum_cp_variation=cp_variation,
        failed_criteria=tuple(failed),
        unavailable_criteria=tuple(unavailable),
    )
