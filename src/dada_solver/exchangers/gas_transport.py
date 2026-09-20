"""Explicit dilute-gas transport models; provenance in MICROTUBE_GAS_MODEL.md."""
from dataclasses import dataclass
from typing import Protocol
import math
import numpy as np
from dada_solver import numerical_primitives as numeric


class GasTransportModel(Protocol):
    gas_constant: float
    provenance: str
    mean_free_path_convention: str
    def viscosity(self, temperature: float) -> float: ...
    def conductivity(self, temperature: float) -> float: ...
    def cp(self, temperature: float) -> float: ...
    def mean_free_path(self, pressure: float, temperature: float) -> float: ...


@dataclass(frozen=True)
class DiluteGasTransport:
    """Air/N2/Ar Sutherland; He NIST interpolation; ideal-gas Shomate cp.

    Air cp uses an explicitly approximate 79/21 mole N2/O2 mixture. R stays
    compatible with the existing calorically perfect thermodynamic gas. cp(T)
    here informs transport only, never silently replaces caloric state energy.
    """
    species: str = 'air'
    minimum_temperature: float = 200.
    maximum_temperature: float = 1000.
    provenance: str = 'COMSOL Sutherland tables 5-2/5-3; NIST Shomate; NIST helium table'

    @property
    def mean_free_path_convention(self):
        return 'viscosity_hard_sphere'

    def __post_init__(self):
        if self.species not in ('air','nitrogen','argon','helium'):
            raise ValueError('Unsupported transport species.')
        if not 200 <= self.minimum_temperature < self.maximum_temperature <= 1000:
            raise ValueError('Transport implementation is limited to 200-1000 K.')

    @property
    def gas_constant(self):
        return numeric.gas_constant(numeric.SPECIES.index(self.species))

    def _check(self, t):
        if not math.isfinite(t) or not self.minimum_temperature <= t <= self.maximum_temperature:
            raise ValueError(f'Transport temperature {t:g} K outside declared domain.')

    def viscosity(self, t):
        self._check(t)
        return numeric.viscosity(t, numeric.SPECIES.index(self.species))

    def conductivity(self, t):
        self._check(t)
        return numeric.conductivity(t, numeric.SPECIES.index(self.species))

    _helium = staticmethod(numeric.helium_property)
    _molar_cp = staticmethod(numeric.molar_cp)

    def cp(self,t):
        self._check(t)
        return numeric.transport_cp(t, numeric.SPECIES.index(self.species))

    def density(self,p,t):
        self._check(t)
        if not math.isfinite(p) or p<=0: raise ValueError('Pressure must be positive.')
        return p/(self.gas_constant*t)

    def mean_free_path(self,p,t):
        """Viscosity-based hard-sphere convention: mu/p * sqrt(pi*R*T/2)."""
        self.density(p,t)
        return numeric.mean_free_path(self.viscosity(t),p,self.gas_constant,t)


@dataclass(frozen=True)
class ConstantGasTransport:
    """Explicit legacy/screening transport; not the production property law."""
    gas_constant: float
    dynamic_viscosity: float
    thermal_conductivity: float
    heat_capacity_cp: float
    provenance: str = 'Explicit caller-supplied legacy/screening constants'

    @property
    def mean_free_path_convention(self):
        return 'viscosity_hard_sphere'

    def __post_init__(self):
        if any(not math.isfinite(v) or v<=0 for v in (self.gas_constant,self.dynamic_viscosity,self.thermal_conductivity,self.heat_capacity_cp)):
            raise ValueError('Positive finite legacy properties required.')
        if self.heat_capacity_cp<=self.gas_constant: raise ValueError('cp must exceed R.')
    def viscosity(self,t): return self.dynamic_viscosity
    def conductivity(self,t): return self.thermal_conductivity
    def cp(self,t): return self.heat_capacity_cp
    def mean_free_path(self,p,t):
        if p<=0 or t<=0: raise ValueError('Positive state required.')
        return numeric.mean_free_path(self.dynamic_viscosity,p,self.gas_constant,t)
