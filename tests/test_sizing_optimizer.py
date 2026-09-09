from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import pytest

from dada_solver.configuration import load_simulation_configuration
from dada_solver.sizing.constraints import ConstraintValue
from dada_solver.sizing.design import DesignParameter, DesignPoint, DesignVariable, apply_design_point
from dada_solver.sizing.objectives import ObjectiveValue
from dada_solver.sizing.optimizer import OptimizationSettings, SlsqpSizingOptimizer
from dada_solver.sizing.problem import SizingProblem


EXAMPLE_PATH = (
    Path(__file__).parents[1] / "examples" / "harmonic_controlled_example.toml"
)


class AlgebraicEvaluator:
    def __init__(self) -> None:
        self.base = load_simulation_configuration(EXAMPLE_PATH)

    def evaluate(self, point: DesignPoint):
        return SimpleNamespace(
            point=point,
            configuration=apply_design_point(self.base, point),
            usable=True,
        )


@dataclass(frozen=True)
class TargetUaObjective:
    name: str = "target_cold_ua"

    def evaluate(self, evaluation) -> ObjectiveValue:
        difference = evaluation.configuration.cold_thermal_conductance - 7.0
        return ObjectiveValue(self.name, difference * difference, True)


@dataclass(frozen=True)
class MinimumUaConstraint:
    name: str = "minimum_cold_ua"

    def evaluate(self, evaluation) -> ConstraintValue:
        margin = evaluation.configuration.cold_thermal_conductance - 6.0
        return ConstraintValue(self.name, margin, margin >= 0.0, True)


def test_slsqp_uses_explicit_normalization_and_finds_feasible_point() -> None:
    problem = SizingProblem(
        variables=(
            DesignVariable(DesignParameter.COLD_UA, 1.0, 20.0, 10.0),
        ),
        objective=TargetUaObjective(),
        constraints=(MinimumUaConstraint(),),
        evaluator=AlgebraicEvaluator(),  # type: ignore[arg-type]
    )
    settings = OptimizationSettings(
        objective_scale=1.0,
        constraint_scales={"minimum_cold_ua": 1.0},
        unavailable_objective_penalty=1.0e6,
        unavailable_constraint_margin=-1.0,
        maximum_iterations=100,
        function_tolerance=1.0e-12,
    )

    result = SlsqpSizingOptimizer(problem, settings).optimize()

    assert result.success
    assert result.best_assessment.feasible
    assert result.best_point.values[DesignParameter.COLD_UA] == pytest.approx(
        7.0, abs=1.0e-5
    )


def test_optimizer_requires_a_scale_for_every_constraint() -> None:
    problem = SizingProblem(
        variables=(DesignVariable(DesignParameter.COLD_UA, 1.0, 20.0, 10.0),),
        objective=TargetUaObjective(),
        constraints=(MinimumUaConstraint(),),
        evaluator=AlgebraicEvaluator(),  # type: ignore[arg-type]
    )
    settings = OptimizationSettings(
        objective_scale=1.0,
        constraint_scales={},
        unavailable_objective_penalty=1.0e6,
        unavailable_constraint_margin=-1.0,
        maximum_iterations=100,
        function_tolerance=1.0e-9,
    )

    with pytest.raises(ValueError, match="Missing explicit scales"):
        SlsqpSizingOptimizer(problem, settings)
