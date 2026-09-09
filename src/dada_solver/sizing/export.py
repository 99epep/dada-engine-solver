"""Machine-readable export of sizing feasibility assessments."""

from __future__ import annotations

import csv
from pathlib import Path

from dada_solver.sizing.feasibility import FeasibilityStudyResult


def write_feasibility_csv(
    result: FeasibilityStudyResult,
    path: str | Path,
) -> None:
    """Write design values, statuses, objective values and constraint margins."""

    parameters = sorted(
        {
            parameter
            for assessment in result.assessments
            for parameter in assessment.evaluation.point.values
        },
        key=lambda item: item.value,
    )
    constraint_names = sorted(
        {
            constraint.name
            for assessment in result.assessments
            for constraint in assessment.constraints
        }
    )
    fieldnames = (
        [parameter.value for parameter in parameters]
        + ["evaluation_status", "feasible", "objective_name", "objective_value"]
        + [f"constraint_{name}_available" for name in constraint_names]
        + [f"constraint_{name}_satisfied" for name in constraint_names]
        + [f"constraint_{name}_margin" for name in constraint_names]
    )
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for assessment in result.assessments:
            row: dict[str, object] = {
                parameter.value: assessment.evaluation.point.values.get(parameter)
                for parameter in parameters
            }
            row.update(
                {
                    "evaluation_status": assessment.evaluation.status.value,
                    "feasible": assessment.feasible,
                    "objective_name": assessment.objective.name,
                    "objective_value": assessment.objective.value,
                }
            )
            by_name = {item.name: item for item in assessment.constraints}
            for name in constraint_names:
                constraint = by_name[name]
                row[f"constraint_{name}_available"] = constraint.available
                row[f"constraint_{name}_satisfied"] = constraint.satisfied
                row[f"constraint_{name}_margin"] = constraint.margin
            writer.writerow(row)
