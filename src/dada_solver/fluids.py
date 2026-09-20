"""Working-fluid models and consistency checks."""

from __future__ import annotations

from dataclasses import dataclass
import math
from dada_solver import numerical_primitives as numeric


@dataclass(frozen=True, slots=True)
class CaloricallyPerfectGas:
    """Single-phase ideal gas with constant caloric properties in SI units."""

    gas_constant: float
    heat_capacity_cp: float
    heat_capacity_cv: float
    consistency_tolerance: float = 1.0e-10

    def __post_init__(self) -> None:
        values = (
            self.gas_constant,
            self.heat_capacity_cp,
            self.heat_capacity_cv,
        )
        if not all(math.isfinite(value) and value > 0.0 for value in values):
            raise ValueError("Gas properties must be finite and strictly positive.")
        if self.heat_capacity_cp <= self.heat_capacity_cv:
            raise ValueError("Cp must be greater than Cv.")
        if not math.isclose(
            self.gas_constant,
            self.heat_capacity_cp - self.heat_capacity_cv,
            rel_tol=self.consistency_tolerance,
            abs_tol=0.0,
        ):
            raise ValueError("Gas properties must satisfy R = Cp - Cv.")

    @property
    def heat_capacity_ratio(self) -> float:
        """Return gamma = Cp / Cv."""

        return self.heat_capacity_cp / self.heat_capacity_cv

    def temperature(self, mass: float, internal_energy: float) -> float:
        """Reconstruct temperature from mass and internal energy."""

        if not math.isfinite(mass) or mass <= 0.0:
            raise ValueError("Mass must be finite and strictly positive.")
        if not math.isfinite(internal_energy) or internal_energy <= 0.0:
            raise ValueError("Internal energy must be finite and strictly positive.")
        return numeric.temperature(mass,internal_energy,self.heat_capacity_cv)

    def pressure(self, mass: float, internal_energy: float, volume: float) -> float:
        """Reconstruct pressure from the conservative state and volume."""

        if not math.isfinite(volume) or volume <= 0.0:
            raise ValueError("Volume must be finite and strictly positive.")
        temperature = self.temperature(mass, internal_energy)
        return numeric.pressure(mass,self.gas_constant,temperature,volume)

    def mass(self, pressure: float, temperature: float, volume: float) -> float:
        """Return ideal-gas mass for a prescribed equilibrium state."""

        if not math.isfinite(pressure) or pressure <= 0.0:
            raise ValueError("Pressure must be finite and strictly positive.")
        if not math.isfinite(temperature) or temperature <= 0.0:
            raise ValueError("Temperature must be finite and strictly positive.")
        if not math.isfinite(volume) or volume <= 0.0:
            raise ValueError("Volume must be finite and strictly positive.")
        return pressure * volume / (self.gas_constant * temperature)

    def internal_energy(self, mass: float, temperature: float) -> float:
        """Return internal energy for a prescribed mass and temperature."""

        if not math.isfinite(mass) or mass <= 0.0:
            raise ValueError("Mass must be finite and strictly positive.")
        if not math.isfinite(temperature) or temperature <= 0.0:
            raise ValueError("Temperature must be finite and strictly positive.")
        return mass * self.heat_capacity_cv * temperature

    def enthalpy(self, temperature: float) -> float:
        """Return specific enthalpy at the given temperature."""

        if not math.isfinite(temperature) or temperature <= 0.0:
            raise ValueError("Temperature must be finite and strictly positive.")
        return numeric.enthalpy(self.heat_capacity_cp, temperature)

    def compressibility_factor(self, pressure: float, temperature: float) -> float:
        """Return the first-level ideal-gas compressibility factor."""

        if pressure <= 0.0 or temperature <= 0.0:
            raise ValueError("Pressure and temperature must be strictly positive.")
        return 1.0

