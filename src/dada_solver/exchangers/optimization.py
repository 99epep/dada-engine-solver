"""Local constrained refinement of screened parallel-channel geometries."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math
from typing import Callable, Iterable

import numpy as np
from scipy.optimize import minimize

from dada_solver.exchangers.models import (
    ExchangerOperatingPoint,
    ParallelRectangularChannels,
    ThermalResistanceModel,
    TransportProperties,
    evaluate_parallel_channels,
)
from dada_solver.exchangers.sizing import (
    ExchangerRequirements,
    ExchangerSizingResult,
    assess_exchanger,
    pareto_exchanger_designs,
    screen_parallel_channel_designs,
)
from dada_solver.fluids import CaloricallyPerfectGas


class ExchangerObjective(Enum):
    """One physical objective; constraints are never folded into it."""

    MINIMUM_GAS_VOLUME = "minimum_gas_volume"
    MINIMUM_PRESSURE_DROP = "minimum_pressure_drop"
    MINIMUM_GAS_SIDE_AREA = "minimum_gas_side_area"


@dataclass(frozen=True, slots=True)
class ContinuousGeometryBounds:
    """Positive manufacturing bounds for continuous channel dimensions."""

    minimum_width: float
    maximum_width: float
    minimum_height: float
    maximum_height: float
    minimum_length: float
    maximum_length: float

    def __post_init__(self) -> None:
        for minimum, maximum, name in (
            (self.minimum_width, self.maximum_width, "width"),
            (self.minimum_height, self.maximum_height, "height"),
            (self.minimum_length, self.maximum_length, "length"),
        ):
            if not math.isfinite(minimum) or minimum <= 0.0:
                raise ValueError(f"Minimum channel {name} must be finite and positive.")
            if not math.isfinite(maximum) or maximum <= minimum:
                raise ValueError(f"Maximum channel {name} must exceed its minimum.")


@dataclass(frozen=True, slots=True)
class ExchangerOptimizationResult:
    objective: ExchangerObjective
    success: bool
    message: str
    best: ExchangerSizingResult | None
    feasible_candidates: tuple[ExchangerSizingResult, ...]
    pareto_candidates: tuple[ExchangerSizingResult, ...]
    screened_design_count: int
    local_optimization_count: int


def optimize_parallel_channel_geometry(
    *,
    channel_counts: Iterable[int],
    bounds: ContinuousGeometryBounds,
    seed_widths: Iterable[float],
    seed_heights: Iterable[float],
    seed_lengths: Iterable[float],
    operating_point: ExchangerOperatingPoint,
    gas: CaloricallyPerfectGas,
    transport: TransportProperties,
    thermal_resistances: ThermalResistanceModel,
    requirements: ExchangerRequirements,
    objective: ExchangerObjective = ExchangerObjective.MINIMUM_GAS_VOLUME,
    surface_roughness: float = 0.0,
    minor_loss_coefficient: float = 0.0,
    maximum_iterations: int = 200,
    function_tolerance: float = 1.0e-9,
) -> ExchangerOptimizationResult:
    """Screen a deterministic grid, then refine feasible seeds with SLSQP.

    Each integer channel count is handled independently. Only feasible grid
    points seed local optimization, preventing an unavailable correlation from
    being disguised by an arbitrary numerical penalty. This is not a proof of
    global optimality.
    """

    counts = tuple(dict.fromkeys(channel_counts))
    if not counts:
        raise ValueError("At least one channel count is required.")
    if maximum_iterations <= 0:
        raise ValueError("Maximum iterations must be positive.")
    if not math.isfinite(function_tolerance) or function_tolerance <= 0.0:
        raise ValueError("Function tolerance must be finite and positive.")

    screened = screen_parallel_channel_designs(
        channel_counts=counts,
        channel_widths=seed_widths,
        channel_heights=seed_heights,
        channel_lengths=seed_lengths,
        operating_point=operating_point,
        gas=gas,
        transport=transport,
        thermal_resistances=thermal_resistances,
        requirements=requirements,
        surface_roughness=surface_roughness,
        minor_loss_coefficient=minor_loss_coefficient,
    )
    feasible_seeds = [item for item in screened if item.feasible]
    if not feasible_seeds:
        return ExchangerOptimizationResult(
            objective=objective,
            success=False,
            message="No feasible deterministic seed was found within the supplied bounds.",
            best=None,
            feasible_candidates=(),
            pareto_candidates=(),
            screened_design_count=len(screened),
            local_optimization_count=0,
        )

    candidates = list(feasible_seeds)
    local_count = 0
    for count in counts:
        count_seeds = [
            item for item in feasible_seeds
            if item.performance.geometry.channel_count == count
        ]
        if not count_seeds:
            continue
        seed = min(count_seeds, key=lambda item: _objective_value(item, objective))
        refined = _refine_one_count(
            count=count,
            seed=seed,
            bounds=bounds,
            operating_point=operating_point,
            gas=gas,
            transport=transport,
            thermal_resistances=thermal_resistances,
            requirements=requirements,
            objective=objective,
            surface_roughness=surface_roughness,
            minor_loss_coefficient=minor_loss_coefficient,
            maximum_iterations=maximum_iterations,
            function_tolerance=function_tolerance,
        )
        local_count += 1
        if refined is not None and refined.feasible:
            candidates.append(refined)

    unique = _unique_results(candidates)
    best = min(unique, key=lambda item: _objective_value(item, objective))
    return ExchangerOptimizationResult(
        objective=objective,
        success=True,
        message=(
            "Best feasible result from deterministic screening and local refinement; "
            "global optimality is not established."
        ),
        best=best,
        feasible_candidates=tuple(unique),
        pareto_candidates=pareto_exchanger_designs(unique),
        screened_design_count=len(screened),
        local_optimization_count=local_count,
    )


def _refine_one_count(
    *,
    count: int,
    seed: ExchangerSizingResult,
    bounds: ContinuousGeometryBounds,
    operating_point: ExchangerOperatingPoint,
    gas: CaloricallyPerfectGas,
    transport: TransportProperties,
    thermal_resistances: ThermalResistanceModel,
    requirements: ExchangerRequirements,
    objective: ExchangerObjective,
    surface_roughness: float,
    minor_loss_coefficient: float,
    maximum_iterations: int,
    function_tolerance: float,
) -> ExchangerSizingResult | None:
    lower = np.log(
        [bounds.minimum_width, bounds.minimum_height, bounds.minimum_length]
    )
    upper = np.log(
        [bounds.maximum_width, bounds.maximum_height, bounds.maximum_length]
    )
    span = upper - lower

    def decode(normalized: np.ndarray) -> ParallelRectangularChannels:
        width, height, length = np.exp(lower + normalized * span)
        return ParallelRectangularChannels(
            channel_count=count,
            channel_width=float(width),
            channel_height=float(height),
            channel_length=float(length),
            surface_roughness=surface_roughness,
            minor_loss_coefficient=minor_loss_coefficient,
        )

    def evaluate(normalized: np.ndarray) -> ExchangerSizingResult:
        performance = evaluate_parallel_channels(
            decode(normalized),
            operating_point,
            gas,
            transport,
            thermal_resistances,
        )
        return assess_exchanger(performance, requirements)

    seed_geometry = seed.performance.geometry
    seed_log = np.log(
        [
            seed_geometry.channel_width,
            seed_geometry.channel_height,
            seed_geometry.channel_length,
        ]
    )
    initial = np.clip((seed_log - lower) / span, 0.0, 1.0)
    objective_scale = _objective_value(seed, objective)

    def objective_function(normalized: np.ndarray) -> float:
        return _objective_value(evaluate(normalized), objective) / objective_scale

    constraints: list[dict[str, object]] = []
    for margin in _constraint_margin_functions(requirements):
        constraints.append(
            {
                "type": "ineq",
                "fun": lambda normalized, function=margin: function(
                    evaluate(normalized)
                ),
            }
        )
    result = minimize(
        objective_function,
        initial,
        method="SLSQP",
        bounds=((0.0, 1.0),) * 3,
        constraints=constraints,
        options={"maxiter": maximum_iterations, "ftol": function_tolerance},
    )
    refined = evaluate(np.asarray(result.x, dtype=float))
    return refined if refined.feasible else None


def _constraint_margin_functions(
    requirements: ExchangerRequirements,
) -> tuple[Callable[[ExchangerSizingResult], float], ...]:
    def available(value: float | None, margin: Callable[[float], float]) -> float:
        return -1.0 if value is None else margin(value)

    margins: list[Callable[[ExchangerSizingResult], float]] = [
        lambda item: available(
            item.performance.overall_conductance,
            lambda value: value / requirements.minimum_overall_conductance - 1.0,
        ),
        lambda item: available(
            item.performance.pressure_drop,
            lambda value: 1.0 - value / requirements.maximum_pressure_drop,
        ),
        lambda item: available(
            item.performance.pressure_drop_fraction,
            lambda value: 1.0 - value / requirements.maximum_pressure_drop_fraction,
        ),
        lambda item: 1.0
        - item.performance.mach_number / requirements.maximum_mach_number,
        lambda item: 1.0
        - item.performance.geometry.gas_volume / requirements.maximum_gas_volume,
    ]
    if requirements.minimum_effectiveness is not None:
        margins.append(
            lambda item: available(
                item.performance.reservoir_effectiveness,
                lambda value: value / requirements.minimum_effectiveness - 1.0,
            )
        )
    return tuple(margins)


def _objective_value(
    result: ExchangerSizingResult,
    objective: ExchangerObjective,
) -> float:
    performance = result.performance
    if objective is ExchangerObjective.MINIMUM_GAS_VOLUME:
        return performance.geometry.gas_volume
    if objective is ExchangerObjective.MINIMUM_PRESSURE_DROP:
        return math.inf if performance.pressure_drop is None else performance.pressure_drop
    return performance.geometry.gas_side_area


def _unique_results(
    results: Iterable[ExchangerSizingResult],
) -> list[ExchangerSizingResult]:
    unique: dict[tuple[int, float, float, float], ExchangerSizingResult] = {}
    for item in results:
        geometry = item.performance.geometry
        key = (
            geometry.channel_count,
            round(geometry.channel_width, 14),
            round(geometry.channel_height, 14),
            round(geometry.channel_length, 14),
        )
        unique[key] = item
    return list(unique.values())
