"""Passive pressure-driven check valves with physical hysteresis."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math

from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.hydraulics import HydraulicFlowModel, FlowResult


class ValveState(Enum):
    CLOSED = "closed"
    OPEN = "open"


@dataclass(frozen=True, slots=True)
class ValveTransition:
    previous_state: ValveState
    current_state: ValveState

    @property
    def changed(self) -> bool:
        return self.previous_state is not self.current_state


@dataclass(frozen=True, slots=True)
class PassiveCheckValve:
    """One-way valve controlled only by its pressure difference and hysteresis."""

    flow_model: HydraulicFlowModel
    opening_pressure_difference: float
    closing_pressure_difference: float

    def __post_init__(self) -> None:
        thresholds = (
            self.opening_pressure_difference,
            self.closing_pressure_difference,
        )
        if not all(math.isfinite(value) for value in thresholds):
            raise ValueError("Valve pressure thresholds must be finite.")
        if self.closing_pressure_difference > self.opening_pressure_difference:
            raise ValueError("Valve thresholds must satisfy delta_P_close <= delta_P_open.")

    def updated_state(
        self,
        state: ValveState,
        upstream_pressure: float,
        downstream_pressure: float,
    ) -> ValveTransition:
        """Apply the published passive opening and closing conditions."""

        pressure_difference = upstream_pressure - downstream_pressure
        new_state = state
        if state is ValveState.CLOSED:
            zero_width = (
                self.opening_pressure_difference
                == self.closing_pressure_difference
            )
            opens = (
                pressure_difference > self.opening_pressure_difference
                if zero_width
                else pressure_difference >= self.opening_pressure_difference
            )
            if opens:
                new_state = ValveState.OPEN
        elif pressure_difference <= self.closing_pressure_difference:
            new_state = ValveState.CLOSED
        return ValveTransition(previous_state=state, current_state=new_state)

    def flow(
        self,
        state: ValveState,
        upstream_pressure: float,
        downstream_pressure: float,
        upstream_temperature: float,
        gas: CaloricallyPerfectGas,
    ) -> FlowResult:
        """Return forward flow only; an open valve never permits reverse flow."""

        if state is ValveState.CLOSED:
            return FlowResult(mass_flow_rate=0.0, is_choked=False)
        return self.flow_model.directed_flow(
            upstream_pressure,
            downstream_pressure,
            upstream_temperature,
            gas,
        )
