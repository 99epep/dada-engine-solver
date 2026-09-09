"""Thermodynamic evaluation of one sizing design point."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from dada_solver.configuration import SimulationConfiguration
from dada_solver.dynamics import ThermodynamicModel, ValveTopology
from dada_solver.factory import (
    build_initial_state,
    build_model,
    build_periodic_solver,
    initial_valve_topology,
)
from dada_solver.performance import CyclePerformance, calculate_cycle_performance
from dada_solver.periodic import PeriodicResult, PeriodicStatus
from dada_solver.results import CycleDiagnostics, extract_cycle_diagnostics
from dada_solver.integration import CycleIntegrationResult
from dada_solver.integration import SegmentIntegrationProgress
from typing import Callable
from dada_solver.sizing.design import DesignParameter, DesignPoint, apply_design_point
from dada_solver.state import ThermodynamicState
from dada_solver.validity import ValidityReport, assess_cycle_validity


class EvaluationStatus(Enum):
    CONVERGED = "converged"
    NOT_CONVERGED = "not_converged"
    INVALID_PHYSICAL_STATE = "invalid_physical_state"
    NUMERICAL_FAILURE = "numerical_failure"


@dataclass(frozen=True, slots=True)
class DesignEvaluation:
    point: DesignPoint
    configuration: SimulationConfiguration
    status: EvaluationStatus
    periodic: PeriodicResult
    performance: CyclePerformance | None
    diagnostics: CycleDiagnostics | None
    validity: ValidityReport | None
    model: ThermodynamicModel | None
    cycle: CycleIntegrationResult | None

    @property
    def usable(self) -> bool:
        return self.status is EvaluationStatus.CONVERGED


class ThermodynamicSizingEvaluator:
    """Evaluate and cache expensive periodic simulations by design point."""

    def __init__(self, base_configuration: SimulationConfiguration) -> None:
        self.base_configuration = base_configuration
        self._cache: dict[tuple[tuple[str, float], ...], DesignEvaluation] = {}

    def evaluate(self, point: DesignPoint) -> DesignEvaluation:
        key = tuple(sorted((item.value, float(value)) for item, value in point.values.items()))
        if key in self._cache:
            return self._cache[key]
        configuration = apply_design_point(self.base_configuration, point)
        model = build_model(configuration)
        initial_state = build_initial_state(configuration, model)
        periodic = build_periodic_solver(configuration, model).solve(
            initial_state, initial_valve_topology()
        )
        status = _evaluation_status(periodic.status)
        if periodic.status is not PeriodicStatus.CONVERGED or periodic.final_cycle is None:
            evaluation = DesignEvaluation(
                point, configuration, status, periodic, None, None, None, None, None
            )
        else:
            cycle = periodic.final_cycle
            cycle_initial_state = ThermodynamicState.from_array(cycle.states[:, 0])
            evaluation = DesignEvaluation(
                point=point,
                configuration=configuration,
                status=status,
                periodic=periodic,
                performance=calculate_cycle_performance(
                    cycle, cycle_initial_state, model.signed_angular_speed
                ),
                diagnostics=extract_cycle_diagnostics(cycle, model),
                validity=assess_cycle_validity(cycle, model, configuration.validity),
                model=model,
                cycle=cycle,
            )
        self._cache[key] = evaluation
        return evaluation

    @property
    def cached_evaluation_count(self) -> int:
        return len(self._cache)


def evaluate_configuration(
    configuration: SimulationConfiguration,
    initial_state: ThermodynamicState | None = None,
    initial_topology: ValveTopology | None = None,
    integration_progress_callback: Callable[[SegmentIntegrationProgress], None] | None = None,
) -> DesignEvaluation:
    """Evaluate one configuration with an optional periodic-state warm start."""

    model = build_model(configuration)
    filling_state = build_initial_state(configuration, model)
    if initial_state is None:
        state = filling_state
    else:
        state = _rescale_state_inventory(initial_state, filling_state.total_mass)
    topology = initial_topology or initial_valve_topology()
    periodic = build_periodic_solver(
        configuration, model, integration_progress_callback
    ).solve(state, topology)
    status = _evaluation_status(periodic.status)
    if periodic.status is not PeriodicStatus.CONVERGED or periodic.final_cycle is None:
        return DesignEvaluation(
            DesignPoint({}), configuration, status, periodic, None, None, None, None, None
        )
    cycle = periodic.final_cycle
    cycle_initial_state = ThermodynamicState.from_array(cycle.states[:, 0])
    return DesignEvaluation(
        point=DesignPoint({}),
        configuration=configuration,
        status=status,
        periodic=periodic,
        performance=calculate_cycle_performance(
            cycle, cycle_initial_state, model.signed_angular_speed
        ),
        diagnostics=extract_cycle_diagnostics(cycle, model),
        validity=assess_cycle_validity(cycle, model, configuration.validity),
        model=model,
        cycle=cycle,
    )


def _rescale_state_inventory(
    state: ThermodynamicState, target_total_mass: float
) -> ThermodynamicState:
    """Preserve composition and specific energies while matching charge inventory."""

    factor = target_total_mass / state.total_mass
    values = state.as_array() * factor
    return ThermodynamicState.from_array(values)


def _evaluation_status(status: PeriodicStatus) -> EvaluationStatus:
    if status is PeriodicStatus.CONVERGED:
        return EvaluationStatus.CONVERGED
    if status is PeriodicStatus.INVALID_PHYSICAL_STATE:
        return EvaluationStatus.INVALID_PHYSICAL_STATE
    if status is PeriodicStatus.NUMERICAL_INTEGRATION_FAILURE:
        return EvaluationStatus.NUMERICAL_FAILURE
    return EvaluationStatus.NOT_CONVERGED
