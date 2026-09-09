"""One-at-a-time sensitivity studies around a configured design point."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np

from dada_solver.sizing.design import DesignParameter, DesignPoint
from dada_solver.sizing.problem import SizingAssessment, SizingProblem


class SamplingScale(Enum):
    LINEAR = "linear"
    LOGARITHMIC = "logarithmic"


@dataclass(frozen=True, slots=True)
class SensitivityStudyResult:
    parameter: DesignParameter
    scale: SamplingScale
    assessments: tuple[SizingAssessment, ...]


def run_one_at_a_time_sensitivity(
    problem: SizingProblem,
    parameter: DesignParameter,
    point_count: int,
    scale: SamplingScale = SamplingScale.LINEAR,
    lower_bound: float | None = None,
    upper_bound: float | None = None,
) -> SensitivityStudyResult:
    """Sweep one configured variable while holding all others at initial values."""

    if point_count < 2:
        raise ValueError("Sensitivity point count must be at least two.")
    try:
        variable = next(item for item in problem.variables if item.parameter is parameter)
    except StopIteration as error:
        raise ValueError(f"Parameter {parameter.value} is not a configured variable.") from error
    if (lower_bound is None) != (upper_bound is None):
        raise ValueError("Local sensitivity bounds must be provided together.")
    sweep_lower = variable.lower_bound if lower_bound is None else lower_bound
    sweep_upper = variable.upper_bound if upper_bound is None else upper_bound
    if sweep_lower < variable.lower_bound or sweep_upper > variable.upper_bound:
        raise ValueError("Local sensitivity bounds must lie within the configured bounds.")
    if sweep_lower >= sweep_upper:
        raise ValueError("Sensitivity lower bound must be less than its upper bound.")
    if scale is SamplingScale.LOGARITHMIC:
        if sweep_lower <= 0.0:
            raise ValueError("Logarithmic sensitivity requires a positive lower bound.")
        values = np.geomspace(sweep_lower, sweep_upper, point_count)
    else:
        values = np.linspace(sweep_lower, sweep_upper, point_count)
    initial_values = dict(problem.initial_point.values)
    assessments = []
    for value in values:
        point_values = dict(initial_values)
        point_values[parameter] = float(value)
        assessments.append(problem.assess(DesignPoint(point_values)))
    return SensitivityStudyResult(parameter, scale, tuple(assessments))
