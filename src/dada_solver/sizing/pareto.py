"""Objective comparison and non-dominated filtering without hidden weights."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from dada_solver.sizing.design import DesignPoint
from dada_solver.sizing.evaluator import ThermodynamicSizingEvaluator
from dada_solver.sizing.objectives import ObjectiveValue, SizingObjective
from dada_solver.sizing.problem import SizingAssessment
from dada_solver.sizing.constraints import SizingConstraint


@dataclass(frozen=True, slots=True)
class MultiObjectiveAssessment:
    point: DesignPoint
    assessment: SizingAssessment
    objectives: tuple[ObjectiveValue, ...]

    @property
    def eligible(self) -> bool:
        return self.assessment.feasible and all(item.available for item in self.objectives)


def assess_objective_set(
    points: Iterable[DesignPoint],
    evaluator: ThermodynamicSizingEvaluator,
    objectives: tuple[SizingObjective, ...],
    constraints: tuple[SizingConstraint, ...],
) -> tuple[MultiObjectiveAssessment, ...]:
    """Evaluate supplied points once and report every requested objective."""

    if not objectives:
        raise ValueError("At least one objective is required.")
    from dada_solver.sizing.problem import SizingProblem

    primary = objectives[0]
    problem = SizingProblem((), primary, constraints, evaluator)
    results = []
    for point in points:
        assessment = problem.assess(point)
        results.append(
            MultiObjectiveAssessment(
                point=point,
                assessment=assessment,
                objectives=tuple(
                    objective.evaluate(assessment.evaluation)
                    for objective in objectives
                ),
            )
        )
    return tuple(results)


def pareto_front(
    assessments: Iterable[MultiObjectiveAssessment],
) -> tuple[MultiObjectiveAssessment, ...]:
    """Return feasible non-dominated points for minimization objectives."""

    eligible = tuple(item for item in assessments if item.eligible)
    front = []
    for candidate in eligible:
        if not any(
            _dominates(other, candidate)
            for other in eligible
            if other is not candidate
        ):
            front.append(candidate)
    return tuple(front)


def _dominates(
    first: MultiObjectiveAssessment,
    second: MultiObjectiveAssessment,
) -> bool:
    first_values = tuple(item.value for item in first.objectives)
    second_values = tuple(item.value for item in second.objectives)
    assert all(value is not None for value in first_values + second_values)
    return all(
        first_value <= second_value
        for first_value, second_value in zip(first_values, second_values, strict=True)
    ) and any(
        first_value < second_value
        for first_value, second_value in zip(first_values, second_values, strict=True)
    )

