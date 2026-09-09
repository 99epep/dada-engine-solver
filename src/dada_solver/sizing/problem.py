"""Composition of design variables, objectives and physical constraints."""

from __future__ import annotations

from dataclasses import dataclass

from dada_solver.sizing.constraints import ConstraintValue, SizingConstraint
from dada_solver.sizing.design import DesignPoint, DesignVariable
from dada_solver.sizing.evaluator import DesignEvaluation, ThermodynamicSizingEvaluator
from dada_solver.sizing.objectives import ObjectiveValue, SizingObjective


@dataclass(frozen=True, slots=True)
class SizingAssessment:
    evaluation: DesignEvaluation
    objective: ObjectiveValue
    constraints: tuple[ConstraintValue, ...]

    @property
    def feasible(self) -> bool:
        return self.evaluation.usable and all(
            item.available and item.satisfied for item in self.constraints
        )


@dataclass(frozen=True, slots=True)
class SizingProblem:
    variables: tuple[DesignVariable, ...]
    objective: SizingObjective
    constraints: tuple[SizingConstraint, ...]
    evaluator: ThermodynamicSizingEvaluator

    def __post_init__(self) -> None:
        parameters = [variable.parameter for variable in self.variables]
        if len(parameters) != len(set(parameters)):
            raise ValueError("Each design parameter may appear only once.")

    @property
    def initial_point(self) -> DesignPoint:
        return DesignPoint(
            {variable.parameter: variable.initial_value for variable in self.variables}
        )

    def assess(self, point: DesignPoint) -> SizingAssessment:
        evaluation = self.evaluator.evaluate(point)
        return SizingAssessment(
            evaluation=evaluation,
            objective=self.objective.evaluate(evaluation),
            constraints=tuple(
                constraint.evaluate(evaluation) for constraint in self.constraints
            ),
        )
