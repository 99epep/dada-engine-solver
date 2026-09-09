"""Plain-text scientific simulation reports."""

from __future__ import annotations

import math

from dada_solver.performance import CyclePerformance, OperatingMode
from dada_solver.periodic import PeriodicResult
from dada_solver.results import CycleDiagnostics
from dada_solver.validity import ValidityReport
from dada_solver.humidity import MoistureScreeningReport
from dada_solver.nomenclature import study_name


def format_simulation_report(
    periodic: PeriodicResult,
    performance: CyclePerformance,
    diagnostics: CycleDiagnostics,
    validity: ValidityReport,
    total_gas_inventory: float,
    include_integration_segments: bool = False,
    moisture: MoistureScreeningReport | None = None,
) -> str:
    """Return a compact English diagnostic report with explicit SI units."""

    last_convergence = periodic.history[-1]
    total_elapsed = sum(
        item.integration_statistics.elapsed_seconds for item in periodic.history
    )
    total_rhs = sum(
        item.integration_statistics.right_hand_side_evaluations
        for item in periodic.history
    )
    total_events = sum(
        item.integration_statistics.event_function_evaluations
        for item in periodic.history
    )
    worst_segment = max(
        last_convergence.integration_statistics.segments,
        key=lambda item: item.right_hand_side_evaluations,
    )
    lines = [
        "DADA thermodynamic simulation",
        f"periodic_status = {periodic.status.value}",
        f"cycles_completed = {last_convergence.cycle_number}",
        f"periodic_normalized_error = {last_convergence.normalized_state_error:.6e}",
        f"integration_elapsed_seconds = {total_elapsed:.6e}",
        f"right_hand_side_evaluations = {total_rhs}",
        f"event_function_evaluations = {total_events}",
        "final_cycle_solve_segments = "
        f"{last_convergence.integration_statistics.solve_segments}",
        "final_cycle_jacobian_evaluations = "
        f"{last_convergence.integration_statistics.jacobian_evaluations}",
        "final_cycle_linear_decompositions = "
        f"{last_convergence.integration_statistics.linear_decompositions}",
        "final_cycle_worst_segment_angle_deg = "
        f"[{math.degrees(worst_segment.start_angle):.6f}, "
        f"{math.degrees(worst_segment.end_angle):.6f}]",
        "final_cycle_worst_segment_rhs_evaluations = "
        f"{worst_segment.right_hand_side_evaluations}",
        f"cycle_topology = {diagnostics.topology.classification.value}",
        f"validity_verdict = {validity.verdict.value}",
        f"operating_mode = {performance.operating_mode.value}",
        f"total_gas_inventory_kg = {total_gas_inventory:.12e}",
        "",
        "Cycle-integrated quantities",
        "exchangers = H_i:heat_in, H_o:heat_out",
        "angle_coordinate = increasing_cycle_progress_phi",
        f"heat_in_J_per_cycle = {performance.heat_in_per_cycle:.12e}",
        f"heat_out_J_per_cycle = {performance.heat_out_per_cycle:.12e}",
        f"heat_in_power_W = {performance.heat_in_power:.12e}",
        f"heat_out_power_W = {performance.heat_out_power:.12e}",
        f"gas_work_J_per_cycle = {performance.gas_work_per_cycle:.12e}",
        f"gas_power_W = {performance.gas_power:.12e}",
        f"mechanical_input_J_per_cycle = {performance.mechanical_input_per_cycle:.12e}",
        f"thermal_efficiency = {_optional_number(performance.thermal_efficiency)}",
        f"motor_power_W = {_optional_number(performance.motor_power)}",
        f"cooling_COP = {_optional_number(performance.cooling_cop)}",
        f"heating_COP = {_optional_number(performance.heating_cop)}",
        *(
            [] if performance.operating_mode is OperatingMode.MOTOR else [
                f"cooling_power_W = {performance.cooling_power:.12e}",
                f"heating_power_W = {performance.heating_power:.12e}",
            ]
        ),
        f"mechanical_input_power_W = {performance.mechanical_input_power:.12e}",
        "",
        "Conservation",
        f"absolute_mass_residual_kg = {performance.conservation.absolute_mass_residual:.12e}",
        f"relative_mass_residual = {performance.conservation.relative_mass_residual:.12e}",
        f"absolute_energy_residual_J = {performance.conservation.absolute_energy_residual:.12e}",
        f"relative_energy_residual = {performance.conservation.relative_energy_residual:.12e}",
        "",
        "Validity indicators",
        f"maximum_pressure_equalization_error = {validity.maximum_pressure_equalization_error:.12e}",
        f"H_i_isothermality_error = {validity.cold_isothermality_error:.12e}",
        f"H_o_isothermality_error = {validity.hot_isothermality_error:.12e}",
        f"maximum_mach_number = {_optional_number(validity.maximum_mach_number)}",
        f"mach_number_status = {validity.mach_number_status}",
        f"maximum_compressibility_deviation = {validity.maximum_compressibility_deviation:.12e}",
        f"maximum_cp_variation = {validity.maximum_cp_variation:.12e}",
        "orifice_pressure_regularization_Pa = "
        f"{diagnostics.orifice_pressure_regularization:.12e}",
        "regularization_sample_fraction = "
        f"{diagnostics.regularization_sample_fraction:.12e}",
    ]
    if diagnostics.topology.reasons:
        lines.append("topology_reasons = " + ", ".join(diagnostics.topology.reasons))
    if validity.failed_criteria:
        lines.append("failed_validity_criteria = " + ", ".join(validity.failed_criteria))
    if validity.unavailable_criteria:
        lines.append(
            "unavailable_validity_criteria = "
            + ", ".join(validity.unavailable_criteria)
        )
    if moisture is not None:
        lines.extend(("", "Moisture screening"))
        lines.append(f"moisture_phase_change_status = {moisture.verdict.value}")
        lines.append(
            "initial_water_mole_fraction = "
            + _optional_number(moisture.initial_water_mole_fraction)
        )
        if moisture.unavailable_reason is not None:
            lines.append(f"moisture_unavailable_reason = {moisture.unavailable_reason}")
        for name, item in moisture.control_volumes.items():
            name = study_name(name)
            lines.append(
                f"{name}_maximum_saturation_ratio = "
                f"{item.maximum_saturation_ratio:.12e}"
            )
            lines.append(
                f"{name}_maximum_saturation_ratio_angle_deg = "
                f"{item.angle_degrees_at_maximum:.9f}"
            )
            lines.append(
                f"{name}_first_saturated_angle_deg = "
                + _optional_number(item.first_saturated_angle_degrees)
            )
            lines.append(
                f"{name}_first_predicted_condensed_phase = "
                + (item.first_predicted_phase or "none")
            )

    lines.extend(("", "Pressure extrema"))
    for name, extrema in diagnostics.pressure_extrema.items():
        name = study_name(name)
        lines.append(
            f"{name}_pressure_Pa = [{extrema.minimum:.12e}, {extrema.maximum:.12e}]"
        )
    lines.extend(("", "Temperature extrema"))
    for name, extrema in diagnostics.temperature_extrema.items():
        name = study_name(name)
        lines.append(
            f"{name}_temperature_K = [{extrema.minimum:.12e}, {extrema.maximum:.12e}]"
        )
    lines.extend(("", "Mass-flow extrema"))
    for name, extrema in diagnostics.mass_flow_extrema.items():
        name = study_name(name)
        lines.append(
            f"{name}_mass_flow_kg_per_s = [{extrema.minimum:.12e}, {extrema.maximum:.12e}]"
        )
    lines.extend(("", "Valve events"))
    if not diagnostics.valve_events:
        lines.append("none")
    for event in diagnostics.valve_events:
        lines.append(
            f"angle_deg={event.angle_degrees:.9f}, valve={study_name(event.valve_name)}, "
            f"transition={event.transition}, Lambda_S={event.small_cylinder_lambda:.9f}, "
            f"Lambda_L={event.large_cylinder_lambda:.9f}"
        )
    if include_integration_segments:
        lines.extend(("", "Final-cycle integration segments"))
        for segment in last_convergence.integration_statistics.segments:
            lines.append(
                f"angle_deg=[{math.degrees(segment.start_angle):.9f}, "
                f"{math.degrees(segment.end_angle):.9f}], "
                f"H_o_to_S={segment.hot_to_small_state.value}, "
                f"H_i_to_L={segment.cold_to_large_state.value}, "
                f"ended_at_event={str(segment.ended_at_event).lower()}, "
                f"elapsed_s={segment.elapsed_seconds:.6e}, "
                f"rhs={segment.right_hand_side_evaluations}, "
                f"jac={segment.jacobian_evaluations}, "
                f"lu={segment.linear_decompositions}, "
                f"samples={segment.stored_samples}"
            )
    return "\n".join(lines)


def _optional_number(value: float | None) -> str:
    return "unavailable" if value is None else f"{value:.12e}"
