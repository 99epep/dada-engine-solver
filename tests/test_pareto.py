from types import SimpleNamespace

from dada_solver.sizing.design import DesignParameter, DesignPoint
from dada_solver.sizing.objectives import ObjectiveValue
from dada_solver.sizing.pareto import MultiObjectiveAssessment, pareto_front


def item(name: str, x: float, y: float, feasible: bool = True):
    point = DesignPoint({DesignParameter.COLD_UA: x + 1.0})
    return MultiObjectiveAssessment(
        point=point,
        assessment=SimpleNamespace(feasible=feasible),
        objectives=(
            ObjectiveValue("first", x, True),
            ObjectiveValue("second", y, True),
        ),
    )


def test_pareto_front_keeps_tradeoffs_and_removes_dominated_or_infeasible_points() -> None:
    first = item("a", 1.0, 4.0)
    second = item("b", 2.0, 2.0)
    third = item("c", 4.0, 1.0)
    dominated = item("d", 3.0, 3.0)
    infeasible = item("e", 0.5, 0.5, feasible=False)

    front = pareto_front((first, second, third, dominated, infeasible))

    assert front == (first, second, third)


def test_unavailable_objective_excludes_point_from_pareto_front() -> None:
    candidate = MultiObjectiveAssessment(
        point=DesignPoint({DesignParameter.COLD_UA: 5.0}),
        assessment=SimpleNamespace(feasible=True),
        objectives=(ObjectiveValue("first", None, False),),
    )

    assert pareto_front((candidate,)) == ()
