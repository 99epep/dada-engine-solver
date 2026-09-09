"""Pointwise comparison of reference CdA and geometric hydraulic closures."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.hydraulics import HydraulicFlowModel


@dataclass(frozen=True, slots=True)
class HydraulicComparisonPoint:
    downstream_pressure: float
    pressure_difference: float
    reference_mass_flow: float
    geometric_mass_flow: float
    relative_flow_difference: float | None
    reference_choked: bool
    geometric_choked: bool


def compare_hydraulic_models(
    reference: HydraulicFlowModel,
    geometric: HydraulicFlowModel,
    upstream_pressure: float,
    downstream_pressures: Iterable[float],
    upstream_temperature: float,
    gas: CaloricallyPerfectGas,
) -> tuple[HydraulicComparisonPoint, ...]:
    """Compare closures over caller-selected pressure ratios without fitting."""

    points = []
    for downstream in downstream_pressures:
        reference_result = reference.directed_flow(
            upstream_pressure, downstream, upstream_temperature, gas
        )
        geometric_result = geometric.directed_flow(
            upstream_pressure, downstream, upstream_temperature, gas
        )
        relative = None
        if reference_result.mass_flow_rate > 0.0:
            relative = (
                geometric_result.mass_flow_rate - reference_result.mass_flow_rate
            ) / reference_result.mass_flow_rate
        points.append(
            HydraulicComparisonPoint(
                downstream_pressure=downstream,
                pressure_difference=upstream_pressure - downstream,
                reference_mass_flow=reference_result.mass_flow_rate,
                geometric_mass_flow=geometric_result.mass_flow_rate,
                relative_flow_difference=relative,
                reference_choked=reference_result.is_choked,
                geometric_choked=geometric_result.is_choked,
            )
        )
    return tuple(points)
