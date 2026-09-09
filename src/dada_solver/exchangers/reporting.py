"""Plain-text reporting for geometric heat-exchanger screening."""

from __future__ import annotations

from dada_solver.exchangers.sizing import ExchangerSizingResult
from dada_solver.exchangers.optimization import ExchangerOptimizationResult


def format_exchanger_screening_report(
    results: tuple[ExchangerSizingResult, ...],
) -> str:
    feasible = [item for item in results if item.feasible]
    unavailable = [item for item in results if item.unavailable_reasons]
    lines = [
        "DADA heat-exchanger screening",
        f"evaluated_designs = {len(results)}",
        f"feasible_designs = {len(feasible)}",
        f"designs_with_unavailable_results = {len(unavailable)}",
    ]
    if not feasible:
        lines.extend(("", "No feasible design was found in the supplied grid."))
        return "\n".join(lines)
    lines.extend(("", "Feasible designs"))
    for index, result in enumerate(
        sorted(feasible, key=lambda item: item.performance.geometry.gas_volume), 1
    ):
        item = result.performance
        geometry = item.geometry
        lines.extend(
            (
                f"design_{index}",
                f"  channel_count = {geometry.channel_count}",
                f"  channel_width = {geometry.channel_width:.12e}",
                f"  channel_height = {geometry.channel_height:.12e}",
                f"  channel_length = {geometry.channel_length:.12e}",
                f"  gas_volume = {geometry.gas_volume:.12e}",
                f"  total_flow_area = {geometry.total_flow_area:.12e}",
                f"  overall_conductance = {_optional(item.overall_conductance)}",
                f"  pressure_drop = {_optional(item.pressure_drop)}",
                f"  mach_number = {item.mach_number:.12e}",
                f"  reynolds_number = {item.reynolds_number:.12e}",
                f"  regime = {item.regime.value}",
            )
        )
    return "\n".join(lines)


def _optional(value: float | None) -> str:
    return "unavailable" if value is None else f"{value:.12e}"


def format_exchanger_optimization_report(
    result: ExchangerOptimizationResult,
) -> str:
    lines = [
        "DADA heat-exchanger geometry optimization",
        f"success = {str(result.success).lower()}",
        f"objective = {result.objective.value}",
        f"message = {result.message}",
        f"screened_designs = {result.screened_design_count}",
        f"local_optimizations = {result.local_optimization_count}",
        f"pareto_designs = {len(result.pareto_candidates)}",
    ]
    if result.best is None:
        return "\n".join(lines)
    performance = result.best.performance
    geometry = performance.geometry
    lines.extend(
        (
            "",
            "Best feasible geometry",
            f"channel_count = {geometry.channel_count}",
            f"channel_width_m = {geometry.channel_width:.12e}",
            f"channel_height_m = {geometry.channel_height:.12e}",
            f"channel_length_m = {geometry.channel_length:.12e}",
            f"gas_volume_m3 = {geometry.gas_volume:.12e}",
            f"gas_side_area_m2 = {geometry.gas_side_area:.12e}",
            f"total_flow_area_m2 = {geometry.total_flow_area:.12e}",
            f"overall_conductance_W_per_K = {_optional(performance.overall_conductance)}",
            f"pressure_drop_Pa = {_optional(performance.pressure_drop)}",
            f"pressure_drop_fraction = {_optional(performance.pressure_drop_fraction)}",
            f"mach_number = {performance.mach_number:.12e}",
            f"reynolds_number = {performance.reynolds_number:.12e}",
            f"regime = {performance.regime.value}",
        )
    )
    return "\n".join(lines)
