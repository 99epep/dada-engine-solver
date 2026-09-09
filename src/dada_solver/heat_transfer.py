"""Zero-dimensional heat-transfer closures."""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True, slots=True)
class ReservoirHeatTransfer:
    """Finite-UA heat transfer between a gas volume and a fixed reservoir."""

    conductance: float
    reservoir_temperature: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.conductance) or self.conductance < 0.0:
            raise ValueError("Thermal conductance must be finite and non-negative.")
        if (
            not math.isfinite(self.reservoir_temperature)
            or self.reservoir_temperature <= 0.0
        ):
            raise ValueError("Reservoir temperature must be finite and positive.")

    def heat_rate(self, gas_temperature: float) -> float:
        """Return heat received by the gas in W."""

        if not math.isfinite(gas_temperature) or gas_temperature <= 0.0:
            raise ValueError("Gas temperature must be finite and positive.")
        return self.conductance * (self.reservoir_temperature - gas_temperature)

