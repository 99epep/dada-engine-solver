"""Plain-text report for coupled cycle/exchanger sizing."""

from dada_solver.exchangers.coupling import CoupledExchangerSizingResult
from dada_solver.nomenclature import study_name


def format_coupled_exchanger_report(result: CoupledExchangerSizingResult) -> str:
    lines = [
        "DADA coupled cycle/exchanger sizing",
        f"status = {study_name(result.status.value)}",
        f"message = {result.message}",
        f"iterations = {len(result.history)}",
        "",
        "Iteration history",
    ]
    if not result.history:
        lines.append("none")
    for item in result.history:
        lines.append(
            f"iteration={item.iteration}, normalized_change={item.normalized_change:.6e}, "
            f"H_i_flow={item.cold_screening_mass_flow:.6e}, "
            f"H_o_flow={item.hot_screening_mass_flow:.6e}, "
            f"H_i_UA={item.cold_ua:.6e}, H_o_UA={item.hot_ua:.6e}, "
            f"H_i_volume={item.cold_gas_volume:.6e}, "
            f"H_o_volume={item.hot_gas_volume:.6e}"
        )
    lines.extend(
        (
            "",
            "Final cycle inputs",
            f"H_i_UA_W_per_K = {result.configuration.cold_thermal_conductance:.12e}",
            f"H_o_UA_W_per_K = {result.configuration.hot_thermal_conductance:.12e}",
            "H_i_heat_exchanger_volume_m3 = "
            f"{result.configuration.machine_volumes.cold_heat_exchanger:.12e}",
            "H_o_heat_exchanger_volume_m3 = "
            f"{result.configuration.machine_volumes.hot_heat_exchanger:.12e}",
        )
    )
    for side, optimization in (
        ("H_i", result.cold_optimization),
        ("H_o", result.hot_optimization),
    ):
        lines.extend(("", f"{side} exchanger geometry"))
        if optimization is None or optimization.best is None:
            lines.append("unavailable")
            continue
        performance = optimization.best.performance
        geometry = performance.geometry
        lines.extend(
            (
                f"channel_count = {geometry.channel_count}",
                f"channel_width_m = {geometry.channel_width:.12e}",
                f"channel_height_m = {geometry.channel_height:.12e}",
                f"channel_length_m = {geometry.channel_length:.12e}",
                f"total_flow_area_m2 = {geometry.total_flow_area:.12e}",
                f"gas_side_area_m2 = {geometry.gas_side_area:.12e}",
                f"pressure_drop_Pa = {_optional(performance.pressure_drop)}",
                f"pressure_drop_fraction = {_optional(performance.pressure_drop_fraction)}",
                f"mach_number = {performance.mach_number:.12e}",
                f"regime = {performance.regime.value}",
            )
        )
    lines.extend(("", "Final cycle performance"))
    performance = None if result.evaluation is None else result.evaluation.performance
    if performance is None:
        lines.append("unavailable")
    else:
        lines.extend(
            (
                f"H_i_heat_received_power_W = {performance.heat_in_power:.12e}",
                f"H_o_heat_received_power_W = {performance.heat_out_power:.12e}",
                f"gas_power_W = {performance.gas_power:.12e}",
                f"motor_power_W = {_optional(performance.motor_power)}",
                f"thermal_efficiency = {_optional(performance.thermal_efficiency)}",
                f"mechanical_input_power_W = {performance.mechanical_input_power:.12e}",
                f"cooling_COP = {_optional(performance.cooling_cop)}",
            )
        )
    reference = (
        None
        if result.reference_evaluation is None
        else result.reference_evaluation.performance
    )
    lines.extend(("", "Reference CdA comparison"))
    if reference is None or performance is None or not result.history:
        lines.append("unavailable")
    else:
        lines.extend(
            (
                f"reference_H_i_heat_received_power_W = {reference.heat_in_power:.12e}",
                f"geometric_H_i_heat_received_power_W = {performance.heat_in_power:.12e}",
                "H_i_heat_received_power_change_W = "
                f"{performance.heat_in_power - reference.heat_in_power:.12e}",
                "reference_mechanical_input_power_W = "
                f"{reference.mechanical_input_power:.12e}",
                "geometric_mechanical_input_power_W = "
                f"{performance.mechanical_input_power:.12e}",
            )
        )
    return "\n".join(lines)


def _optional(value: float | None) -> str:
    return "unavailable" if value is None else f"{value:.12e}"
