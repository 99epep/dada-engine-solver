"""Normalized constrained optimization adapter for sizing problems."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import OptimizeResult, minimize

from dada_solver.sizing.design import DesignPoint
from dada_solver.sizing.problem import SizingAssessment, SizingProblem


@dataclass(frozen=True, slots=True)
class OptimizationSettings:
    """Explicit numerical scales and termination controls for SLSQP."""

    objective_scale: float
    constraint_scales: Mapping[str, float]
    unavailable_objective_penalty: float
    unavailable_constraint_margin: float
    maximum_iterations: int
    function_tolerance: float

    def __post_init__(self) -> None:
        for name, value in (
            ("objective scale", self.objective_scale),
            ("unavailable objective penalty", self.unavailable_objective_penalty),
            ("function tolerance", self.function_tolerance),
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")
        if (
            not math.isfinite(self.unavailable_constraint_margin)
            or self.unavailable_constraint_margin >= 0.0
        ):
            raise ValueError("Unavailable constraint margin must be finite and negative.")
        if self.maximum_iterations <= 0:
            raise ValueError("Maximum iteration count must be positive.")
        for name, value in self.constraint_scales.items():
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"Constraint scale {name} must be finite and positive.")


@dataclass(frozen=True, slots=True)
class SizingOptimizationResult:
    success: bool
    message: str
    iteration_count: int
    function_evaluation_count: int
    best_point: DesignPoint
    best_assessment: SizingAssessment
    raw_result: OptimizeResult


class SlsqpSizingOptimizer:
    """Run SLSQP without selecting a physical objective or hidden scale."""

    def __init__(self, problem: SizingProblem, settings: OptimizationSettings) -> None:
        self.problem = problem
        self.settings = settings
        names = {constraint.name for constraint in problem.constraints}
        missing = names - set(settings.constraint_scales)
        extra = set(settings.constraint_scales) - names
        if missing:
            raise ValueError(
                "Missing explicit scales for constraints: " + ", ".join(sorted(missing))
            )
        if extra:
            raise ValueError(
                "Constraint scales do not match the problem: "
                + ", ".join(sorted(extra))
            )

    def optimize(self) -> SizingOptimizationResult:
        variables = self.problem.variables
        if not variables:
            raise ValueError("At least one design variable is required for optimization.")
        initial = np.array(
            [
                (variable.initial_value - variable.lower_bound)
                / (variable.upper_bound - variable.lower_bound)
                for variable in variables
            ],
            dtype=float,
        )

        def assess(normalized: NDArray[np.float64]) -> SizingAssessment:
            return self.problem.assess(self._decode(normalized))

        def objective(normalized: NDArray[np.float64]) -> float:
            value = assess(normalized).objective
            if not value.available or value.value is None:
                return self.settings.unavailable_objective_penalty
            return value.value / self.settings.objective_scale

        scipy_constraints = []
        for index, constraint in enumerate(self.problem.constraints):
            scale = self.settings.constraint_scales[constraint.name]

            def constraint_function(
                normalized: NDArray[np.float64],
                constraint_index: int = index,
                constraint_scale: float = scale,
            ) -> float:
                value = assess(normalized).constraints[constraint_index]
                if not value.available or value.margin is None:
                    return self.settings.unavailable_constraint_margin
                return value.margin / constraint_scale

            scipy_constraints.append({"type": "ineq", "fun": constraint_function})

        raw = minimize(
            objective,
            initial,
            method="SLSQP",
            bounds=[(0.0, 1.0)] * len(variables),
            constraints=scipy_constraints,
            options={
                "maxiter": self.settings.maximum_iterations,
                "ftol": self.settings.function_tolerance,
                "disp": False,
            },
        )
        best_point = self._decode(np.asarray(raw.x, dtype=float))
        best_assessment = self.problem.assess(best_point)
        success = bool(raw.success) and best_assessment.feasible
        message = str(raw.message)
        if raw.success and not best_assessment.feasible:
            message = "Optimizer terminated successfully but the final point is infeasible."
        return SizingOptimizationResult(
            success=success,
            message=message,
            iteration_count=int(raw.nit),
            function_evaluation_count=int(raw.nfev),
            best_point=best_point,
            best_assessment=best_assessment,
            raw_result=raw,
        )

    def _decode(self, normalized: NDArray[np.float64]) -> DesignPoint:
        if normalized.shape != (len(self.problem.variables),):
            raise ValueError("Normalized design vector has an invalid shape.")
        values = {}
        for coordinate, variable in zip(
            np.clip(normalized, 0.0, 1.0), self.problem.variables, strict=True
        ):
            values[variable.parameter] = variable.lower_bound + coordinate * (
                variable.upper_bound - variable.lower_bound
            )
        return DesignPoint(values)
