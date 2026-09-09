"""Periodic steady-state search strategies."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math

import numpy as np

from dada_solver.dynamics import ValveTopology
from dada_solver.integration import (
    CycleIntegrationResult,
    CycleIntegrator,
    IntegrationStatistics,
)
from dada_solver.state import ThermodynamicState


class PeriodicStatus(Enum):
    CONVERGED = "converged"
    NOT_CONVERGED = "not_converged"
    MAXIMUM_CYCLE_COUNT_REACHED = "maximum_cycle_count_reached"
    INVALID_PHYSICAL_STATE = "invalid_physical_state"
    NUMERICAL_INTEGRATION_FAILURE = "numerical_integration_failure"


@dataclass(frozen=True, slots=True)
class CycleConvergence:
    cycle_number: int
    normalized_state_error: float
    topology_recovered: bool
    absolute_mass_residual: float
    relative_mass_residual: float
    integration_statistics: IntegrationStatistics


@dataclass(frozen=True, slots=True)
class PeriodicResult:
    status: PeriodicStatus
    message: str
    history: tuple[CycleConvergence, ...]
    final_cycle: CycleIntegrationResult | None


@dataclass(frozen=True, slots=True)
class SuccessiveCycleSolver:
    """Robust fixed-point iteration by repeated complete cycles."""

    integrator: CycleIntegrator
    maximum_cycles: int = 500
    relative_tolerance: float = 1.0e-7
    mass_absolute_tolerance: float = 1.0e-12
    energy_absolute_tolerance: float = 1.0e-6

    def __post_init__(self) -> None:
        if self.maximum_cycles <= 0:
            raise ValueError("Maximum cycle count must be positive.")
        for value in (
            self.relative_tolerance,
            self.mass_absolute_tolerance,
            self.energy_absolute_tolerance,
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError("Periodic convergence tolerances must be positive.")

    def solve(
        self,
        initial_state: ThermodynamicState,
        initial_topology: ValveTopology,
    ) -> PeriodicResult:
        state = initial_state
        topology = initial_topology
        reference_mass = initial_state.total_mass
        history: list[CycleConvergence] = []
        final_cycle: CycleIntegrationResult | None = None

        for cycle_number in range(1, self.maximum_cycles + 1):
            final_cycle = self.integrator.integrate_cycle(state, topology)
            if not final_cycle.completed:
                status = (
                    PeriodicStatus.INVALID_PHYSICAL_STATE
                    if final_cycle.status.value == "invalid_physical_state"
                    else PeriodicStatus.NUMERICAL_INTEGRATION_FAILURE
                )
                return PeriodicResult(
                    status,
                    final_cycle.message,
                    tuple(history),
                    final_cycle,
                )

            end_state = final_cycle.final_state
            error = self._normalized_error(state, end_state)
            topology_recovered = (
                self.integrator.model.continuous_ideal_diodes
                or final_cycle.final_topology == topology
            )
            absolute_mass_residual = end_state.total_mass - reference_mass
            relative_mass_residual = absolute_mass_residual / reference_mass
            history.append(
                CycleConvergence(
                    cycle_number,
                    error,
                    topology_recovered,
                    absolute_mass_residual,
                    relative_mass_residual,
                    final_cycle.statistics,
                )
            )
            if error <= 1.0 and topology_recovered:
                return PeriodicResult(
                    PeriodicStatus.CONVERGED,
                    "Periodic steady state converged.",
                    tuple(history),
                    final_cycle,
                )
            state = end_state
            topology = final_cycle.final_topology

        return PeriodicResult(
            PeriodicStatus.MAXIMUM_CYCLE_COUNT_REACHED,
            "Maximum cycle count reached before periodic convergence.",
            tuple(history),
            final_cycle,
        )

    def _normalized_error(
        self,
        start: ThermodynamicState,
        end: ThermodynamicState,
    ) -> float:
        start_values = start.as_array()
        end_values = end.as_array()
        absolute_scales = np.empty(8, dtype=float)
        absolute_scales[0::2] = self.mass_absolute_tolerance
        absolute_scales[1::2] = self.energy_absolute_tolerance
        scales = absolute_scales + self.relative_tolerance * np.maximum(
            np.abs(start_values), np.abs(end_values)
        )
        return float(np.sqrt(np.mean(((end_values - start_values) / scales) ** 2)))
