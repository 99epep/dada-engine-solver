"""Interchangeable scalar sizing objectives."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from dada_solver.factory import build_model
from dada_solver.sizing.evaluator import DesignEvaluation


@dataclass(frozen=True, slots=True)
class ObjectiveValue:
    name: str
    value: float | None
    available: bool


class SizingObjective(Protocol):
    name: str

    def evaluate(self, evaluation: DesignEvaluation) -> ObjectiveValue:
        """Return a value to minimize."""


@dataclass(frozen=True, slots=True)
class MaximizeCoolingCop:
    name: str = "maximize_cooling_cop"

    def evaluate(self, evaluation: DesignEvaluation) -> ObjectiveValue:
        cop = evaluation.performance.cooling_cop if evaluation.performance else None
        return ObjectiveValue(self.name, -cop if cop is not None else None, cop is not None)


@dataclass(frozen=True, slots=True)
class MaximizeThermalEfficiency:
    name: str = "maximize_thermal_efficiency"

    def evaluate(self, evaluation: DesignEvaluation) -> ObjectiveValue:
        efficiency = evaluation.performance.thermal_efficiency if evaluation.performance else None
        return ObjectiveValue(
            self.name, -efficiency if efficiency is not None else None,
            efficiency is not None,
        )


@dataclass(frozen=True, slots=True)
class MaximizeMotorPower:
    name: str = "maximize_motor_power"

    def evaluate(self, evaluation: DesignEvaluation) -> ObjectiveValue:
        power = evaluation.performance.motor_power if evaluation.performance else None
        return ObjectiveValue(self.name, -power if power is not None else None, power is not None)


@dataclass(frozen=True, slots=True)
class MinimizeTotalSweptVolume:
    name: str = "minimize_total_swept_volume"

    def evaluate(self, evaluation: DesignEvaluation) -> ObjectiveValue:
        volumes = evaluation.configuration.machine_volumes
        value = volumes.small_cylinder.swept + volumes.large_cylinder.swept
        return ObjectiveValue(self.name, value, True)


@dataclass(frozen=True, slots=True)
class MinimizeChargingPressure:
    name: str = "minimize_charging_pressure"

    def evaluate(self, evaluation: DesignEvaluation) -> ObjectiveValue:
        configuration = evaluation.configuration
        filling_volume = build_model(configuration).volumes(0.0).total
        pressure = configuration.charge.resolved_pressure(
            configuration.gas,
            filling_volume,
        )
        return ObjectiveValue(self.name, pressure, True)


@dataclass(frozen=True, slots=True)
class MinimizeTotalUa:
    name: str = "minimize_total_ua"

    def evaluate(self, evaluation: DesignEvaluation) -> ObjectiveValue:
        value = (
            evaluation.configuration.cold_thermal_conductance
            + evaluation.configuration.hot_thermal_conductance
        )
        return ObjectiveValue(self.name, value, True)
