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
        pressure_history = replay.pressures
    else:
        pressures: list[np.ndarray] = []
        for angle, values in zip(cycle.angles, cycle.states.T, strict=True):
            state = ThermodynamicState.from_array(values)
            pressures.append(state.pressures(model.gas, model.volumes(float(angle))))
        pressure_history = np.asarray(pressures).T

    pair_references = np.maximum(
        0.5 * (pressure_history[[0, 1]] + pressure_history[[2, 3]]),
        np.finfo(float).tiny,
    )
    pair_differences = np.abs(
        pressure_history[[0, 1]] - pressure_history[[2, 3]]
    )
    pressure_error = float(np.max(pair_differences / pair_references))
    compressibility_deviation = 0.0
    cp_variation = 0.0

    from dada_solver.fluids import CaloricallyPerfectGas
    if type(model.gas) is not CaloricallyPerfectGas:
        states=[p for angle,values in zip(cycle.angles,cycle.states.T,strict=True)
                for p in ThermodynamicState.from_array(values).fluid_states(model.gas,model.volumes(float(angle)))]
        if any(p.compressibility_factor is None or p.cp is None for p in states):
            raise ValueError('This validity policy requires compressibility and Cp diagnostics.')
        compressibility_deviation=max(abs(p.compressibility_factor-1) for p in states)
        capacities=[p.cp for p in states]
        cp_variation=(max(capacities)-min(capacities))/min(capacities)
    failed: list[str] = []
    # Pair pressure differences are resolved by the network, not a validity veto.
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
        maximum_mach_number=maximum_mach,
        mach_number_status=mach_status,
        maximum_compressibility_deviation=compressibility_deviation,
        maximum_cp_variation=cp_variation,
        failed_criteria=tuple(failed),
        unavailable_criteria=tuple(unavailable),
    )
