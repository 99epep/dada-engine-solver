"""Text reports for ideal external cooling-cell requirements."""

from __future__ import annotations

from dada_solver.thermal_load_configuration import CoolingCellScenario
from dada_solver.thermal_loads import CoolingTaskAssessment


def format_cooling_cell_requirements(scenario: CoolingCellScenario) -> str:
    load = scenario.load
    reference = scenario.reference_task
    ambitious = scenario.ambitious_task
    human = scenario.human_power
    return "\n".join(
        (
            "Cooling cell ideal thermal-load requirements",
            f"material = {load.material_name}",
            f"mass_kg = {load.mass:.12e}",
            f"sensible_energy_J = {load.sensible_energy:.12e}",
            f"phase_change_energy_J = {load.phase_change_energy:.12e}",
            f"ideal_total_energy_J = {load.ideal_cooling_energy:.12e}",
            f"reference_time_s = {reference.allowed_time:.12e}",
            f"reference_cooling_power_W = {reference.required_average_cooling_power:.12e}",
            f"reference_COP_at_nominal_human_power = {reference.required_cop(human.nominal_power):.12e}",
            f"ambitious_time_s = {ambitious.allowed_time:.12e}",
            f"ambitious_cooling_power_W = {ambitious.required_average_cooling_power:.12e}",
            f"ambitious_COP_at_nominal_human_power = {ambitious.required_cop(human.nominal_power):.12e}",
            "load_scope = ideal_only_excludes_container_and_parasitic_heat_gains",
        )
    )


def format_cooling_task_assessment(assessment: CoolingTaskAssessment) -> str:
    return "\n".join(
        (
            "Cooling task at fixed machine operating point",
            f"ideal_load_energy_J = {assessment.ideal_load_energy:.12e}",
            f"cooling_energy_J_per_cycle = {assessment.cooling_energy_per_cycle:.12e}",
            f"average_cooling_power_W = {assessment.average_cooling_power:.12e}",
            f"required_mechanical_input_power_W = {assessment.required_mechanical_input_power:.12e}",
            f"cooling_COP = {_optional(assessment.cooling_cop)}",
            f"ideal_estimated_completion_time_s = {_optional(assessment.ideal_estimated_completion_time)}",
            f"required_pedaling_power_for_target_time_W = {_optional(assessment.required_pedaling_power_for_target_time)}",
            f"operable_at_nominal_human_power = {str(assessment.operable_at_nominal_human_power).lower()}",
            f"target_met_at_operating_point = {str(assessment.target_met_at_operating_point).lower()}",
        )
    )


def _optional(value: float | None) -> str:
    return "unavailable" if value is None else f"{value:.12e}"

