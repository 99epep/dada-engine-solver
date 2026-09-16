"""Explicit dilute-gas transport models; provenance in MICROTUBE_GAS_MODEL.md."""
from dataclasses import dataclass
from typing import Protocol
import math
import numpy as np


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
        return {'air':287.05,'nitrogen':296.803,'argon':208.132,'helium':2077.1}[self.species]

    def _check(self, t):
        if not math.isfinite(t) or not self.minimum_temperature <= t <= self.maximum_temperature:
            raise ValueError(f'Transport temperature {t:g} K outside declared domain.')

    def viscosity(self, t):
        self._check(t)
        if self.species == 'helium': return self._helium(t, False)
        mu, s = {'air':(1.716e-5,111.),'nitrogen':(1.663e-5,107.),'argon':(2.125e-5,114.)}[self.species]
        return mu*(t/273.)**1.5*(273.+s)/(t+s)

    def conductivity(self, t):
        self._check(t)
        if self.species == 'helium': return self._helium(t, True)
        k, s = {'air':(.0241,194.),'nitrogen':(.0242,150.),'argon':(.0163,170.)}[self.species]
        return k*(t/273.)**1.5*(273.+s)/(t+s)

    @staticmethod
    def _helium(t, conductivity):
        temperatures=(200,225,250,273.15,275,300,325,350,375,400,450,500,600,700,800,900,1000)
        values=((118.5,128.2,137.7,146.2,146.9,155.9,164.6,173.2,181.6,189.9,206.,221.7,251.9,280.9,308.9,336.1,362.6)
                if conductivity else (15.1,16.4,17.6,18.7,18.8,19.9,21.,22.1,23.2,24.3,26.3,28.3,32.2,35.9,39.5,43.,46.4))
        return float(np.interp(t,temperatures,values))*(1e-3 if conductivity else 1e-6)

    @staticmethod
    def _molar_cp(t, oxygen=False):
        if oxygen:
            a,b,c,d,e=((31.32234,-20.23531,57.86644,-36.50624,-.007374) if t<700 else
                       (30.03235,8.772972,-3.988133,.788313,-.741599))
        else:
            a,b,c,d,e=((28.98641,1.853978,-9.647459,16.63537,.000117) if t<500 else
                       (19.50583,19.88705,-8.598535,1.369784,.527601))
        x=t/1000
        return a+b*x+c*x*x+d*x**3+e/x**2

    def cp(self,t):
        self._check(t)
        if self.species in ('argon','helium'): return 2.5*self.gas_constant
        if self.species=='nitrogen': return self._molar_cp(t)/.0280134
        return (.79*self._molar_cp(t)+.21*self._molar_cp(t,True))/(.79*.0280134+.21*.0319988)

    def density(self,p,t):
        self._check(t)
        if not math.isfinite(p) or p<=0: raise ValueError('Pressure must be positive.')
        return p/(self.gas_constant*t)

    def mean_free_path(self,p,t):
        """Viscosity-based hard-sphere convention: mu/p * sqrt(pi*R*T/2)."""
        self.density(p,t)
        return self.viscosity(t)/p*math.sqrt(math.pi*self.gas_constant*t/2)


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
        return self.dynamic_viscosity/p*math.sqrt(math.pi*self.gas_constant*t/2)
