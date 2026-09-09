"""Plain-text reports for constrained sizing results."""

from __future__ import annotations

from dada_solver.sizing.optimizer import SizingOptimizationResult
from dada_solver.sizing.problem import SizingAssessment
from dada_solver.sizing.feasibility import FeasibilityStudyResult


def format_optimization_report(result: SizingOptimizationResult) -> str:
    lines = [
        "DADA constrained sizing",
        f"optimization_success = {str(result.success).lower()}",
        f"message = {result.message}",
        f"iterations = {result.iteration_count}",
        f"function_evaluations = {result.function_evaluation_count}",
        f"feasible = {str(result.best_assessment.feasible).lower()}",
        f"objective_name = {result.best_assessment.objective.name}",
        f"objective_value = {_optional(result.best_assessment.objective.value)}",
        "",
        "Design point",
    ]
    for parameter, value in sorted(
        result.best_point.values.items(), key=lambda item: item[0].value
    ):
        lines.append(f"{parameter.value} = {value:.12e}")
    lines.extend(("", "Constraints"))
    for constraint in result.best_assessment.constraints:
        lines.append(
            f"{constraint.name}: available={str(constraint.available).lower()}, "
            f"satisfied={str(constraint.satisfied).lower()}, "
            f"margin={_optional(constraint.margin)}"
        )
    return "\n".join(lines)


def format_assessment_report(assessment: SizingAssessment) -> str:
    """Report one design point before an optimization is attempted."""

    lines = [
        "DADA sizing point assessment",
        f"evaluation_status = {assessment.evaluation.status.value}",
        f"feasible = {str(assessment.feasible).lower()}",
        f"objective_name = {assessment.objective.name}",
        f"objective_available = {str(assessment.objective.available).lower()}",
        f"objective_value = {_optional(assessment.objective.value)}",
        "",
        "Design point",
    ]
    for parameter, value in sorted(
        assessment.evaluation.point.values.items(), key=lambda item: item[0].value
    ):
        lines.append(f"{parameter.value} = {value:.12e}")
    lines.extend(("", "Constraints"))
    for constraint in assessment.constraints:
        lines.append(
            f"{constraint.name}: available={str(constraint.available).lower()}, "
            f"satisfied={str(constraint.satisfied).lower()}, "
            f"margin={_optional(constraint.margin)}"
        )
    return "\n".join(lines)


def format_feasibility_report(result: FeasibilityStudyResult) -> str:
    lines = [
        "DADA sizing feasibility study",
        f"evaluated_points = {len(result.assessments)}",
        f"feasible_points = {result.feasible_count}",
        "",
        "Evaluation statuses",
    ]
    for name, count in sorted(result.evaluation_status_counts.items()):
        lines.append(f"{name} = {count}")
    lines.extend(("", "Constraint rejection counts"))
    if not result.constraint_rejection_counts:
        lines.append("none")
    for name, count in sorted(result.constraint_rejection_counts.items()):
        lines.append(f"{name} = {count}")
    return "\n".join(lines)


def _optional(value: float | None) -> str:
    return "unavailable" if value is None else f"{value:.12e}"
