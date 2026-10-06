"""Objectives requiring explicit manufacturing context from a built design."""
from dataclasses import dataclass
import math

from dada_solver.performance import OperatingMode
from dada_solver.sizing.objectives import ObjectiveValue
from dada_solver.sizing.evaluator import DesignEvaluation


def cooling_power_per_total_microtube(
    cooling_power_w: float | None,
    heat_in_microtube_count: int,
    heat_out_microtube_count: int,
) -> float | None:
    """Return cooling watts per physical microtube across both exchangers."""
    if (cooling_power_w is None or not math.isfinite(cooling_power_w)
            or heat_in_microtube_count <= 0 or heat_out_microtube_count <= 0):
        return None
    return cooling_power_w / (heat_in_microtube_count + heat_out_microtube_count)


def cooling_cop_times_power_per_total_microtube(
    cooling_cop: float | None,
    cooling_power_per_total_microtube_w: float | None,
) -> float | None:
    """Return COP times cooling watts per physical microtube."""
    if cooling_cop is None or cooling_power_per_total_microtube_w is None:
        return None
    return cooling_cop * cooling_power_per_total_microtube_w


@dataclass(frozen=True, slots=True)
class MaximizeCoolingPowerPerTotalMicrotube:
    name: str = 'maximize_cooling_power_per_total_microtube'

    def evaluate(self, evaluation: DesignEvaluation, *,
                 heat_in_microtube_count: int | None = None,
                 heat_out_microtube_count: int | None = None) -> ObjectiveValue:
        """Minimize the negative score; missing manufacturing context is unavailable."""
        p = evaluation.performance
        value = None
        if (p is not None and p.operating_mode is OperatingMode.REFRIGERATION
                and heat_in_microtube_count is not None
                and heat_out_microtube_count is not None):
            value = cooling_power_per_total_microtube(
                p.cooling_power, heat_in_microtube_count, heat_out_microtube_count)
        return ObjectiveValue(self.name, -value if value is not None else None,
                              value is not None)


@dataclass(frozen=True, slots=True)
class MaximizeCoolingCopTimesPowerPerTotalMicrotube(
        MaximizeCoolingPowerPerTotalMicrotube):
    name: str = 'maximize_cooling_cop_times_power_per_total_microtube'

    def evaluate(self, evaluation: DesignEvaluation, *,
                 heat_in_microtube_count: int | None = None,
                 heat_out_microtube_count: int | None = None) -> ObjectiveValue:
        """Minimize the negative COP-times-productivity score."""
        p = evaluation.performance
        value = None
        if (p is not None and p.operating_mode is OperatingMode.REFRIGERATION
                and heat_in_microtube_count is not None
                and heat_out_microtube_count is not None):
            productivity = cooling_power_per_total_microtube(
                p.cooling_power, heat_in_microtube_count, heat_out_microtube_count)
            value = cooling_cop_times_power_per_total_microtube(
                p.cooling_cop, productivity)
        return ObjectiveValue(self.name, -value if value is not None else None,
                              value is not None)
