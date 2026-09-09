"""Physical sizing constraints expressed as non-negative feasibility margins."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Protocol

from dada_solver.sizing.evaluator import DesignEvaluation
from dada_solver.mechanics import PistonFaceAreas, extract_generalized_gas_torque
from dada_solver.topology import CycleTopologyClassification
from dada_solver.validity import ValidityVerdict
from dada_solver.thermal_loads import CoolingTask
from dada_solver.performance import OperatingMode


@dataclass(frozen=True, slots=True)
class ConstraintValue:
    name: str
    margin: float | None
    satisfied: bool
    available: bool


class SizingConstraint(Protocol):
    name: str

    def evaluate(self, evaluation: DesignEvaluation) -> ConstraintValue:
        """Return margin >= 0 when the constraint is satisfied."""


@dataclass(frozen=True, slots=True)
class MinimumCoolingPower:
    required_power: float
    name: str = "minimum_cooling_power"

    def __post_init__(self) -> None:
        _require_positive_limit(self.name, self.required_power)

    def evaluate(self, evaluation: DesignEvaluation) -> ConstraintValue:
        if evaluation.performance is None or evaluation.performance.operating_mode is not OperatingMode.REFRIGERATION:
            return _unavailable(self.name)
        margin = evaluation.performance.cooling_power - self.required_power
        return ConstraintValue(self.name, margin, margin >= 0.0, True)


@dataclass(frozen=True, slots=True)
class CompleteCoolingTaskWithinTime:
    task: CoolingTask
    name: str = "complete_cooling_task_within_time"

    def evaluate(self, evaluation: DesignEvaluation) -> ConstraintValue:
        if evaluation.performance is None or evaluation.performance.operating_mode is not OperatingMode.REFRIGERATION:
            return _unavailable(self.name)
        margin = (
            evaluation.performance.cooling_power
            - self.task.required_average_cooling_power
        )
        return ConstraintValue(self.name, margin, margin >= 0.0, True)


@dataclass(frozen=True, slots=True)
class MinimumMotorPower:
    required_power: float
    name: str = "minimum_motor_power"

    def __post_init__(self) -> None:
        _require_positive_limit(self.name, self.required_power)

    def evaluate(self, evaluation: DesignEvaluation) -> ConstraintValue:
        power = evaluation.performance.motor_power if evaluation.performance else None
        if power is None:
            return _unavailable(self.name)
        margin = power - self.required_power
        return ConstraintValue(self.name, margin, margin >= 0.0, True)


@dataclass(frozen=True, slots=True)
class MaximumMechanicalInputPower:
    limit: float
    name: str = "maximum_mechanical_input_power"

    def __post_init__(self) -> None:
        _require_positive_limit(self.name, self.limit)

    def evaluate(self, evaluation: DesignEvaluation) -> ConstraintValue:
        if evaluation.performance is None:
            return _unavailable(self.name)
        margin = self.limit - evaluation.performance.mechanical_input_power
        return ConstraintValue(self.name, margin, margin >= 0.0, True)


@dataclass(frozen=True, slots=True)
class MaximumPressure:
    limit: float
    name: str = "maximum_pressure"

    def __post_init__(self) -> None:
        _require_positive_limit(self.name, self.limit)

    def evaluate(self, evaluation: DesignEvaluation) -> ConstraintValue:
        if evaluation.diagnostics is None:
            return _unavailable(self.name)
        maximum = max(item.maximum for item in evaluation.diagnostics.pressure_extrema.values())
        margin = self.limit - maximum
        return ConstraintValue(self.name, margin, margin >= 0.0, True)


@dataclass(frozen=True, slots=True)
class MaximumTemperature:
    limit: float
    name: str = "maximum_temperature"

    def __post_init__(self) -> None:
        _require_positive_limit(self.name, self.limit)

    def evaluate(self, evaluation: DesignEvaluation) -> ConstraintValue:
        if evaluation.diagnostics is None:
            return _unavailable(self.name)
        maximum = max(
            item.maximum for item in evaluation.diagnostics.temperature_extrema.values()
        )
        margin = self.limit - maximum
        return ConstraintValue(self.name, margin, margin >= 0.0, True)


@dataclass(frozen=True, slots=True)
class MaximumAbsoluteMassFlow:
    limit: float
    name: str = "maximum_absolute_mass_flow"

    def __post_init__(self) -> None:
        _require_positive_limit(self.name, self.limit)

    def evaluate(self, evaluation: DesignEvaluation) -> ConstraintValue:
        if evaluation.diagnostics is None:
            return _unavailable(self.name)
        maximum = max(
            max(abs(item.minimum), abs(item.maximum))
            for item in evaluation.diagnostics.mass_flow_extrema.values()
        )
        margin = self.limit - maximum
        return ConstraintValue(self.name, margin, margin >= 0.0, True)


@dataclass(frozen=True, slots=True)
class MaximumMachNumber:
    limit: float
    name: str = "maximum_mach_number"

    def __post_init__(self) -> None:
        _require_positive_limit(self.name, self.limit)

    def evaluate(self, evaluation: DesignEvaluation) -> ConstraintValue:
        if evaluation.validity is None or evaluation.validity.maximum_mach_number is None:
            return _unavailable(self.name)
        margin = self.limit - evaluation.validity.maximum_mach_number
        return ConstraintValue(self.name, margin, margin >= 0.0, True)


@dataclass(frozen=True, slots=True)
class MaximumPressureEqualizationError:
    limit: float
    name: str = "maximum_pressure_equalization_error"

    def __post_init__(self) -> None:
        _require_positive_limit(self.name, self.limit)

    def evaluate(self, evaluation: DesignEvaluation) -> ConstraintValue:
        if evaluation.validity is None:
            return _unavailable(self.name)
        margin = self.limit - evaluation.validity.maximum_pressure_equalization_error
        return ConstraintValue(self.name, margin, margin >= 0.0, True)


@dataclass(frozen=True, slots=True)
class MaximumIsothermalityError:
    limit: float
    name: str = "maximum_isothermality_error"

    def __post_init__(self) -> None:
        _require_positive_limit(self.name, self.limit)

    def evaluate(self, evaluation: DesignEvaluation) -> ConstraintValue:
        if evaluation.validity is None:
            return _unavailable(self.name)
        maximum = max(
            evaluation.validity.cold_isothermality_error,
            evaluation.validity.hot_isothermality_error,
        )
        margin = self.limit - maximum
        return ConstraintValue(self.name, margin, margin >= 0.0, True)


@dataclass(frozen=True, slots=True)
class RequireNominalCycleTopology:
    name: str = "nominal_cycle_topology"

    def evaluate(self, evaluation: DesignEvaluation) -> ConstraintValue:
        if evaluation.diagnostics is None:
            return _unavailable(self.name)
        satisfied = (
            evaluation.diagnostics.topology.classification
            is CycleTopologyClassification.NOMINAL
        )
        return ConstraintValue(self.name, 1.0 if satisfied else -1.0, satisfied, True)


@dataclass(frozen=True, slots=True)
class RequireValidThermodynamicModel:
    name: str = "valid_thermodynamic_model"

    def evaluate(self, evaluation: DesignEvaluation) -> ConstraintValue:
        if evaluation.validity is None:
            return _unavailable(self.name)
        satisfied = evaluation.validity.verdict is ValidityVerdict.VALID
        return ConstraintValue(self.name, 1.0 if satisfied else -1.0, satisfied, True)


@dataclass(frozen=True, slots=True)
class MaximumPistonGasForce:
    limit: float
    piston_areas: PistonFaceAreas
    name: str = "maximum_piston_gas_force"

    def __post_init__(self) -> None:
        _require_positive_limit(self.name, self.limit)

    def evaluate(self, evaluation: DesignEvaluation) -> ConstraintValue:
        if evaluation.diagnostics is None:
            return _unavailable(self.name)
        maximum_small = (
            evaluation.diagnostics.pressure_extrema["S"].maximum
            * self.piston_areas.small
        )
        maximum_large = (
            evaluation.diagnostics.pressure_extrema["L"].maximum
            * self.piston_areas.large
        )
        maximum = max(maximum_small, maximum_large)
        margin = self.limit - maximum
        return ConstraintValue(self.name, margin, margin >= 0.0, True)


@dataclass(frozen=True, slots=True)
class MaximumAbsoluteGeneralizedGasTorque:
    limit: float
    name: str = "maximum_absolute_generalized_gas_torque"

    def __post_init__(self) -> None:
        _require_positive_limit(self.name, self.limit)

    def evaluate(self, evaluation: DesignEvaluation) -> ConstraintValue:
        if evaluation.model is None or evaluation.cycle is None:
            return _unavailable(self.name)
        torque = extract_generalized_gas_torque(evaluation.cycle, evaluation.model)
        maximum = float(abs(torque).max())
        margin = self.limit - maximum
        return ConstraintValue(self.name, margin, margin >= 0.0, True)


@dataclass(frozen=True, slots=True)
class RequirePeriodicConvergence:
    name: str = "periodic_convergence"

    def evaluate(self, evaluation: DesignEvaluation) -> ConstraintValue:
        margin = 1.0 if evaluation.usable else -1.0
        return ConstraintValue(self.name, margin, evaluation.usable, True)


def _unavailable(name: str) -> ConstraintValue:
    return ConstraintValue(name, None, False, False)


def _require_positive_limit(name: str, value: float) -> None:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"Constraint {name} limit must be finite and positive.")
