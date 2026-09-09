"""Constraint-based screening of geometric heat-exchanger candidates."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import math
from typing import Iterable

from dada_solver.exchangers.models import (
    ExchangerOperatingPoint,
    ExchangerPerformance,
    ParallelRectangularChannels,
    ThermalResistanceModel,
    TransportProperties,
    evaluate_parallel_channels,
)
from dada_solver.fluids import CaloricallyPerfectGas


@dataclass(frozen=True, slots=True)
class ExchangerRequirements:
    """Hard constraints for exchanger screening, all in SI units."""

    minimum_overall_conductance: float
    maximum_pressure_drop: float
    maximum_mach_number: float
    maximum_gas_volume: float
    maximum_pressure_drop_fraction: float
    minimum_effectiveness: float | None = None

    def __post_init__(self) -> None:
        for name, value in (
            ("Minimum overall conductance", self.minimum_overall_conductance),
            ("Maximum pressure drop", self.maximum_pressure_drop),
            ("Maximum Mach number", self.maximum_mach_number),
            ("Maximum gas volume", self.maximum_gas_volume),
            ("Maximum pressure-drop fraction", self.maximum_pressure_drop_fraction),
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")
        if self.minimum_effectiveness is not None and not (
            0.0 < self.minimum_effectiveness <= 1.0
        ):
            raise ValueError("Minimum effectiveness must lie in (0, 1].")


@dataclass(frozen=True, slots=True)
class ExchangerSizingResult:
    """One evaluated geometry with individual constraint diagnostics."""

    performance: ExchangerPerformance
    feasible: bool
    unavailable_reasons: tuple[str, ...]
    violated_constraints: tuple[str, ...]


def assess_exchanger(
    performance: ExchangerPerformance,
    requirements: ExchangerRequirements,
) -> ExchangerSizingResult:
    unavailable = []
    violated = []
    if performance.overall_conductance is None:
        unavailable.append("overall_conductance")
    elif performance.overall_conductance < requirements.minimum_overall_conductance:
        violated.append("minimum_overall_conductance")
    if performance.pressure_drop is None:
        unavailable.append("pressure_drop")
    elif performance.pressure_drop > requirements.maximum_pressure_drop:
        violated.append("maximum_pressure_drop")
    if performance.pressure_drop_fraction is None:
        if "pressure_drop" not in unavailable:
            unavailable.append("pressure_drop_fraction")
    elif (
        performance.pressure_drop_fraction
        > requirements.maximum_pressure_drop_fraction
    ):
        violated.append("maximum_pressure_drop_fraction")
    if performance.mach_number > requirements.maximum_mach_number:
        violated.append("maximum_mach_number")
    if performance.geometry.gas_volume > requirements.maximum_gas_volume:
        violated.append("maximum_gas_volume")
    if requirements.minimum_effectiveness is not None:
        if performance.reservoir_effectiveness is None:
            unavailable.append("reservoir_effectiveness")
        elif performance.reservoir_effectiveness < requirements.minimum_effectiveness:
            violated.append("minimum_effectiveness")
    return ExchangerSizingResult(
        performance=performance,
        feasible=not unavailable and not violated,
        unavailable_reasons=tuple(unavailable),
        violated_constraints=tuple(violated),
    )


def screen_parallel_channel_designs(
    *,
    channel_counts: Iterable[int],
    channel_widths: Iterable[float],
    channel_heights: Iterable[float],
    channel_lengths: Iterable[float],
    operating_point: ExchangerOperatingPoint,
    gas: CaloricallyPerfectGas,
    transport: TransportProperties,
    thermal_resistances: ThermalResistanceModel,
    requirements: ExchangerRequirements,
    surface_roughness: float = 0.0,
    minor_loss_coefficient: float = 0.0,
) -> tuple[ExchangerSizingResult, ...]:
    """Evaluate a finite caller-defined design grid without hidden bounds."""

    results = []
    for count, width, height, length in product(
        channel_counts, channel_widths, channel_heights, channel_lengths
    ):
        geometry = ParallelRectangularChannels(
            channel_count=count,
            channel_width=width,
            channel_height=height,
            channel_length=length,
            surface_roughness=surface_roughness,
            minor_loss_coefficient=minor_loss_coefficient,
        )
        performance = evaluate_parallel_channels(
            geometry, operating_point, gas, transport, thermal_resistances
        )
        results.append(assess_exchanger(performance, requirements))
    return tuple(results)


def pareto_exchanger_designs(
    results: Iterable[ExchangerSizingResult],
) -> tuple[ExchangerSizingResult, ...]:
    """Return feasible designs non-dominated in volume, pressure drop, and UA."""

    feasible = [item for item in results if item.feasible]
    front = []
    for candidate in feasible:
        cp = candidate.performance
        assert cp.pressure_drop is not None and cp.overall_conductance is not None
        dominated = False
        for other in feasible:
            if other is candidate:
                continue
            op = other.performance
            assert op.pressure_drop is not None and op.overall_conductance is not None
            no_worse = (
                op.geometry.gas_volume <= cp.geometry.gas_volume
                and op.pressure_drop <= cp.pressure_drop
                and op.overall_conductance >= cp.overall_conductance
            )
            strictly_better = (
                op.geometry.gas_volume < cp.geometry.gas_volume
                or op.pressure_drop < cp.pressure_drop
                or op.overall_conductance > cp.overall_conductance
            )
            if no_worse and strictly_better:
                dominated = True
                break
        if not dominated:
            front.append(candidate)
    return tuple(front)
