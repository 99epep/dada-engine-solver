"""Explicit dilute-gas transport models; provenance in MICROTUBE_GAS_MODEL.md."""
from dataclasses import dataclass
from typing import Protocol
import math
import numpy as np
from dada_solver import numerical_primitives as numeric


PROPERTY_TEMPERATURE_DOMAINS = {
    'air': (100.,1000.), 'nitrogen': (200.,1000.),
    'argon': (200.,1000.), 'helium': (50.,1000.),
}


class TransportDomainError(ValueError):
    """A property-temperature failure, not a phase/EOS validity verdict."""
    def __init__(self, temperature, species, minimum, maximum):
        self.category = ('transport_temperature_nonfinite' if not math.isfinite(temperature)
                         else 'transport_temperature_below_domain' if temperature<minimum
                         else 'transport_temperature_above_domain')
        self.diagnostics = dict(category=self.category, species=species,
            temperature_k=temperature if math.isfinite(temperature) else None,
            minimum_temperature_k=minimum, maximum_temperature_k=maximum,
            domain_kind='property_temperature_domain')
        super().__init__(f'Transport temperature {temperature:g} K outside declared domain '
                         f'[{minimum:g}, {maximum:g}] K for {species}: {self.category}')


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
    """Dilute air/He correlations, N2/Ar Sutherland; ideal-gas Shomate cp.

    Air cp uses an explicitly approximate 79/21 mole N2/O2 mixture. R stays
    compatible with the existing calorically perfect thermodynamic gas. cp(T)
    here informs transport only, never silently replaces caloric state energy.
    """
    species: str = 'air'
    minimum_temperature: float | None = None
    maximum_temperature: float = 1000.
    provenance: str = 'Air: Lemmon/Jacobsen 2004; He: Arp/McCarty/Friend and Hands/Arp; N2/Ar: COMSOL Sutherland; NIST Shomate cp'
    correlation_version: str = 'dilute_species_v2'

    @property
    def mean_free_path_convention(self):
        return 'viscosity_hard_sphere'

    def __post_init__(self):
        if self.species not in ('air','nitrogen','argon','helium'):
            raise ValueError('Unsupported transport species.')
        lower,upper=PROPERTY_TEMPERATURE_DOMAINS[self.species]
        if self.minimum_temperature is None:
            object.__setattr__(self,'minimum_temperature',lower)
        if not lower <= self.minimum_temperature < self.maximum_temperature <= upper:
            raise ValueError(f'{self.species} transport is limited to {lower:g}-{upper:g} K; user bounds may only restrict it.')
        if self.correlation_version!='dilute_species_v2':
            raise ValueError('Unsupported transport correlation version.')

    @property
    def gas_constant(self):
        return numeric.gas_constant(numeric.SPECIES.index(self.species))

    def _check(self, t):
        if not math.isfinite(t) or not self.minimum_temperature <= t <= self.maximum_temperature:
            raise TransportDomainError(t,self.species,self.minimum_temperature,self.maximum_temperature)

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
