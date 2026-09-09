"""Conservative thermodynamic state and ideal-gas reconstruction."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from numpy.typing import NDArray

from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.geometry import InstantaneousVolumes


STATE_NAMES = ("m_S", "U_S", "m_L", "U_L", "m_C", "U_C", "m_H", "U_H")


@dataclass(frozen=True, slots=True)
class ThermodynamicState:
    """Conservative state of the four independent control volumes."""

    mass_small: float
    energy_small: float
    mass_large: float
    energy_large: float
    mass_cold: float
    energy_cold: float
    mass_hot: float
    energy_hot: float

    def __post_init__(self) -> None:
        values = (
            self.mass_small,
            self.energy_small,
            self.mass_large,
            self.energy_large,
            self.mass_cold,
            self.energy_cold,
            self.mass_hot,
            self.energy_hot,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("Every state component must be finite.")
        if any(value <= 0.0 for value in values):
            raise ValueError("Every mass and internal energy must be strictly positive.")

    def as_array(self) -> NDArray[np.float64]:
        """Return state using the published ordering."""

        return np.array(
            [
                self.mass_small,
                self.energy_small,
                self.mass_large,
                self.energy_large,
                self.mass_cold,
                self.energy_cold,
                self.mass_hot,
                self.energy_hot,
            ],
            dtype=float,
        )

    @classmethod
    def from_array(cls, values: NDArray[np.float64]) -> "ThermodynamicState":
        """Build a checked state from the published eight-component ordering."""

        array = np.asarray(values, dtype=float)
        if array.shape != (8,):
            raise ValueError("Thermodynamic state must contain exactly eight values.")
        return cls(*array)

    @property
    def masses(self) -> NDArray[np.float64]:
        return self.as_array()[0::2]

    @property
    def energies(self) -> NDArray[np.float64]:
        return self.as_array()[1::2]

    @property
    def total_mass(self) -> float:
        return float(np.sum(self.masses))

    @property
    def total_internal_energy(self) -> float:
        return float(np.sum(self.energies))

    def temperatures(self, gas: CaloricallyPerfectGas) -> NDArray[np.float64]:
        """Return temperatures ordered as S, L, C, H in K."""

        values = self.as_array()
        return values[1::2] / (values[0::2] * gas.heat_capacity_cv)

    def pressures(
        self,
        gas: CaloricallyPerfectGas,
        volumes: InstantaneousVolumes,
    ) -> NDArray[np.float64]:
        """Return pressures ordered as S, L, C, H in Pa."""

        _, pressures = self.temperatures_and_pressures(gas, volumes)
        return pressures

    def temperatures_and_pressures(
        self,
        gas: CaloricallyPerfectGas,
        volumes: InstantaneousVolumes,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Reconstruct both primitive arrays from one conservative-state copy."""

        values = self.as_array()
        masses = values[0::2]
        temperatures = values[1::2] / (masses * gas.heat_capacity_cv)
        volume_array = np.fromiter(volumes.as_tuple(), dtype=float, count=4)
        pressures = masses * gas.gas_constant * temperatures / volume_array
        return temperatures, pressures


@dataclass(frozen=True, slots=True)
class UniformCharge:
    """Uniform equilibrium filling state used only as an initial condition."""

    pressure: float
    temperature: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.pressure) or self.pressure <= 0.0:
            raise ValueError("Charge pressure must be finite and strictly positive.")
        if not math.isfinite(self.temperature) or self.temperature <= 0.0:
            raise ValueError("Charge temperature must be finite and strictly positive.")

    def create_state(
        self,
        gas: CaloricallyPerfectGas,
        volumes: InstantaneousVolumes,
    ) -> ThermodynamicState:
        """Distribute charge mass according to each volume at uniform P and T."""

        masses = [
            gas.mass(self.pressure, self.temperature, volume)
            for volume in volumes.as_tuple()
        ]
        energies = [gas.internal_energy(mass, self.temperature) for mass in masses]
        return ThermodynamicState(
            masses[0],
            energies[0],
            masses[1],
            energies[1],
            masses[2],
            energies[2],
            masses[3],
            energies[3],
        )

    def total_mass(
        self,
        gas: CaloricallyPerfectGas,
        volumes: InstantaneousVolumes,
    ) -> float:
        """Return the gas inventory defined by the filling configuration."""

        return self.pressure * volumes.total / (
            gas.gas_constant * self.temperature
        )
