from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

from dada_solver.configuration import load_simulation_configuration
from dada_solver.sizing.constraints import ConstraintValue
from dada_solver.sizing.design import DesignParameter, DesignPoint, DesignVariable, apply_design_point
from dada_solver.sizing.export import write_feasibility_csv
from dada_solver.sizing.feasibility import FeasibilityStudySettings, run_feasibility_study
from dada_solver.sizing.objectives import ObjectiveValue
from dada_solver.sizing.problem import SizingProblem
from dada_solver.sizing.sensitivity import SamplingScale, run_one_at_a_time_sensitivity


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
            status=SimpleNamespace(value="converged"),
        )


@dataclass(frozen=True)
class MinimizeUa:
    name: str = "minimize_cold_ua"

    def evaluate(self, evaluation) -> ObjectiveValue:
        return ObjectiveValue(
            self.name,
            evaluation.configuration.cold_thermal_conductance,
            True,
        )


@dataclass(frozen=True)
class MinimumUa:
    name: str = "minimum_cold_ua"

    def evaluate(self, evaluation) -> ConstraintValue:
        margin = evaluation.configuration.cold_thermal_conductance - 5.0
        return ConstraintValue(self.name, margin, margin >= 0.0, True)


def create_problem() -> SizingProblem:
    return SizingProblem(
        variables=(DesignVariable(DesignParameter.COLD_UA, 0.0, 10.0, 6.0),),
        objective=MinimizeUa(),
        constraints=(MinimumUa(),),
        evaluator=AlgebraicEvaluator(),  # type: ignore[arg-type]
    )


def test_latin_hypercube_feasibility_is_reproducible_and_bounded() -> None:
    settings = FeasibilityStudySettings(8, random_seed=42, include_initial_point=True)

    first = run_feasibility_study(create_problem(), settings)
    second = run_feasibility_study(create_problem(), settings)

    first_values = [
        item.evaluation.point.values[DesignParameter.COLD_UA]
        for item in first.assessments
    ]
    second_values = [
        item.evaluation.point.values[DesignParameter.COLD_UA]
        for item in second.assessments
    ]
    assert first_values == second_values
    assert first_values[0] == 6.0
    assert all(0.0 <= value <= 10.0 for value in first_values)
    assert len(first.assessments) == 9
    assert first.evaluation_status_counts == {"converged": 9}
    assert first.constraint_rejection_counts["minimum_cold_ua"] > 0


def test_feasibility_csv_preserves_margins_and_status(tmp_path: Path) -> None:
    result = run_feasibility_study(
        create_problem(),
        FeasibilityStudySettings(2, random_seed=1, include_initial_point=False),
    )
    path = tmp_path / "feasibility.csv"

    write_feasibility_csv(result, path)
    text = path.read_text(encoding="utf-8")

    assert "evaluation_status" in text
    assert "constraint_minimum_cold_ua_margin" in text
    assert "converged" in text


def test_feasibility_settings_reject_invalid_sample_count() -> None:
    import pytest

    with pytest.raises(ValueError, match="positive"):
        FeasibilityStudySettings(0, random_seed=0)


def test_one_at_a_time_sensitivity_includes_bounds() -> None:
    result = run_one_at_a_time_sensitivity(
        create_problem(), DesignParameter.COLD_UA, 3, SamplingScale.LINEAR
    )
    values = [
        assessment.evaluation.point.values[DesignParameter.COLD_UA]
        for assessment in result.assessments
    ]
    assert values == [0.0, 5.0, 10.0]


def test_one_at_a_time_sensitivity_accepts_local_logarithmic_bounds() -> None:
    result = run_one_at_a_time_sensitivity(
        create_problem(),
        DesignParameter.COLD_UA,
        4,
        SamplingScale.LOGARITHMIC,
        lower_bound=1.0,
        upper_bound=8.0,
    )
    values = [
        assessment.evaluation.point.values[DesignParameter.COLD_UA]
        for assessment in result.assessments
    ]
    assert values == [1.0, 2.0, 4.0, 8.0]


def test_one_at_a_time_sensitivity_rejects_local_bounds_outside_design_space() -> None:
    import pytest

    with pytest.raises(ValueError, match="configured bounds"):
        run_one_at_a_time_sensitivity(
            create_problem(),
            DesignParameter.COLD_UA,
            3,
            lower_bound=1.0,
            upper_bound=11.0,
        )
