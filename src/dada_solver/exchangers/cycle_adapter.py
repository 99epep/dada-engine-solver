"""Explicit bridge from cycle diagnostics to exchanger screening inputs."""

from __future__ import annotations

from enum import Enum

from dada_solver.exchangers.models import ExchangerOperatingPoint
from dada_solver.results import CycleDiagnostics


class ExchangerSide(Enum):
    COLD = "cold"
    HOT = "hot"


def conservative_port_flow_operating_point(
    diagnostics: CycleDiagnostics,
    side: ExchangerSide,
    representative_pressure: float,
    representative_temperature: float,
) -> ExchangerOperatingPoint:
    """Use the largest adjacent-port flow as a screening mass flow.

    A lumped exchanger can accumulate mass, so its two port flows need not be
    equal instantaneously. This adapter deliberately uses the larger absolute
    port-flow extremum; it does not claim to reconstruct an internal velocity
    field from the zero-dimensional model.
    """

    names = (
        ("small_to_cold", "cold_to_large")
        if side is ExchangerSide.COLD
        else ("large_to_hot", "hot_to_small")
    )
    maximum_flow = max(
        max(
            abs(diagnostics.mass_flow_extrema[name].minimum),
            abs(diagnostics.mass_flow_extrema[name].maximum),
        )
        for name in names
    )
    return ExchangerOperatingPoint(
        absolute_mass_flow_rate=maximum_flow,
        pressure=representative_pressure,
        temperature=representative_temperature,
    )
