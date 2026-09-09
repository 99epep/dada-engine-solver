"""Reproducible feasibility exploration before local optimization."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import math
from typing import Callable

import numpy as np
from scipy.stats import qmc

from dada_solver.sizing.design import DesignPoint
from dada_solver.sizing.problem import SizingAssessment, SizingProblem


@dataclass(frozen=True, slots=True)
class FeasibilityStudySettings:
    sample_count: int
    random_seed: int
    include_initial_point: bool = True

    def __post_init__(self) -> None:
        if self.sample_count <= 0:
            raise ValueError("Feasibility sample count must be positive.")
        if self.random_seed < 0:
            raise ValueError("Feasibility random seed must be non-negative.")


@dataclass(frozen=True, slots=True)
class FeasibilityStudyResult:
    assessments: tuple[SizingAssessment, ...]
    feasible_assessments: tuple[SizingAssessment, ...]
    evaluation_status_counts: dict[str, int]
    constraint_rejection_counts: dict[str, int]

    @property
    def feasible_count(self) -> int:
        return len(self.feasible_assessments)


ProgressCallback = Callable[[int, int, SizingAssessment], None]


def run_feasibility_study(
    problem: SizingProblem,
    settings: FeasibilityStudySettings,
    progress_callback: ProgressCallback | None = None,
) -> FeasibilityStudyResult:
    """Evaluate a Latin-hypercube design without selecting a preferred point."""

    if not problem.variables:
        raise ValueError("Feasibility exploration requires design variables.")
    sampler = qmc.LatinHypercube(
        d=len(problem.variables),
        seed=settings.random_seed,
    )
    normalized = sampler.random(settings.sample_count)
    points = [_decode(problem, row) for row in normalized]
    if settings.include_initial_point:
        points.insert(0, problem.initial_point)

    assessments = []
    total = len(points)
    for index, point in enumerate(points, start=1):
        assessment = problem.assess(point)
        assessments.append(assessment)
        if progress_callback is not None:
            progress_callback(index, total, assessment)

    feasible = tuple(item for item in assessments if item.feasible)
    statuses = Counter(item.evaluation.status.value for item in assessments)
    rejections: Counter[str] = Counter()
    for assessment in assessments:
        for constraint in assessment.constraints:
            if not constraint.available or not constraint.satisfied:
                rejections[constraint.name] += 1
    return FeasibilityStudyResult(
        assessments=tuple(assessments),
        feasible_assessments=feasible,
        evaluation_status_counts=dict(statuses),
        constraint_rejection_counts=dict(rejections),
    )


def _decode(problem: SizingProblem, normalized: np.ndarray) -> DesignPoint:
    values = {}
    for coordinate, variable in zip(normalized, problem.variables, strict=True):
        if not math.isfinite(float(coordinate)) or not 0.0 <= coordinate <= 1.0:
            raise ValueError("Normalized feasibility coordinate must lie in [0, 1].")
        values[variable.parameter] = variable.lower_bound + float(coordinate) * (
            variable.upper_bound - variable.lower_bound
        )
    return DesignPoint(values)

