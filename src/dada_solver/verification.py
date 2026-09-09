"""Independent analytical residuals derived from the published appendix."""

from __future__ import annotations

import math


def closed_adiabatic_invariant(pressure: float, volume: float, gamma: float) -> float:
    """Return P*V**gamma for a closed adiabatic ideal-gas volume."""

    _require_positive("pressure", pressure)
    _require_positive("volume", volume)
    _require_positive("gamma", gamma)
    return pressure * volume**gamma


def closed_adiabatic_log_rate_residual(
    pressure: float,
    pressure_rate: float,
    volume: float,
    volume_rate: float,
    gamma: float,
) -> float:
    """Return d(ln(P V**gamma))/dt, which must vanish."""

    _require_positive("pressure", pressure)
    _require_positive("volume", volume)
    _require_positive("gamma", gamma)
    return pressure_rate / pressure + gamma * volume_rate / volume


def adiabatic_donor_log_rate_residual(
    mass: float,
    mass_rate: float,
    temperature: float,
    temperature_rate: float,
    volume: float,
    volume_rate: float,
    gamma: float,
) -> float:
    """Return the appendix A.5 differential residual for an outflow-only donor."""

    _require_positive("mass", mass)
    _require_positive("temperature", temperature)
    _require_positive("volume", volume)
    _require_positive("gamma", gamma)
    return temperature_rate / temperature - (gamma - 1.0) * (
        mass_rate / mass - volume_rate / volume
    )


def _require_positive(name: str, value: float) -> None:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive.")

