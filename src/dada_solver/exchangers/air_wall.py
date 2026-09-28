"""Conservative lumped wall storage with a finite-capacity external air stream.

The air passage is quasi-steady against a uniform wall. Wall conductances and
capacity are explicit inputs; Doty's overall UA does not identify their split.
"""
from dataclasses import dataclass
import math
from dada_solver import numerical_primitives as numeric
import time

import numpy as np
from scipy.integrate import solve_ivp

from dada_solver.dynamics import ThermodynamicModel
from dada_solver.exchangers.base import LumpedWallThermalModel
from dada_solver.state import ThermodynamicState
from dada_solver.dynamics import ValveTopology
from dada_solver.valves import ValveState
from dada_solver.integration import IntegrationInterrupted


@dataclass(frozen=True, slots=True)
class WallThermalPoint:
    """Typed numerical facts; dictionaries are only the public rates boundary."""
    gas_heat_w: float
    air_heat_w: float
    wall_energy_rate_w: float
    air_outlet_temperature_k: float | None
    wall_temperature_k: float
    film_diagnostics: tuple = ()

    external_heat_w = property(lambda s: s.air_heat_w)
    external_outlet_temperature_k = property(lambda s: s.air_outlet_temperature_k)

    def rates(self):
        return dict(gas_heat_w=self.gas_heat_w, air_heat_w=self.air_heat_w,
                    wall_energy_rate_w=self.wall_energy_rate_w,
                    air_outlet_temperature_k=self.air_outlet_temperature_k)


@dataclass(frozen=True)
class AirWallExchanger:
    gas_wall_conductance_w_k: float
    air_wall_conductance_w_k: float
    wall_capacity_j_k: float
    air_mass_flow_kg_s: float
    air_cp_j_kg_k: float
    air_inlet_temperature_k: float
    gas_film: object | None = None

    def __post_init__(self):
        for value in (self.wall_capacity_j_k, self.air_cp_j_kg_k, self.air_inlet_temperature_k):
            if not math.isfinite(value) or value <= 0:
                raise ValueError('Wall capacity, air cp and inlet temperature must be positive.')
        for value in (self.gas_wall_conductance_w_k, self.air_wall_conductance_w_k, self.air_mass_flow_kg_s):
            if not math.isfinite(value) or value < 0:
                raise ValueError('Conductances and air flow must be finite and nonnegative.')
        capacity_rate = self.air_mass_flow_kg_s*self.air_cp_j_kg_k
        effective = (capacity_rate * -math.expm1(-self.air_wall_conductance_w_k/capacity_rate)
                     if capacity_rate > 0 else 0.0)
        object.__setattr__(self, '_air_capacity_rate', capacity_rate)
        object.__setattr__(self, '_effective_air_conductance', effective)

    @property
    def external_stream(self):
        from .external_stream import ExternalFluidStream
        return ExternalFluidStream(self.air_inlet_temperature_k, self.air_mass_flow_kg_s,
            self.air_cp_j_kg_k, self.air_wall_conductance_w_k, 'air (legacy)')

    external_inlet_temperature_k = property(lambda s: s.air_inlet_temperature_k)
    external_mass_flow_kg_s = property(lambda s: s.air_mass_flow_kg_s)
    external_cp_j_kg_k = property(lambda s: s.air_cp_j_kg_k)
    effective_external_conductance_w_k = property(lambda s: s._effective_air_conductance)

    @property
    def requires_flow_context(self):
        return self.gas_film is not None

    def rates(self, gas_temperature_k, wall_energy_j, *, context=None):
        return self.thermal_point(gas_temperature_k, wall_energy_j, context=context).rates()

    def thermal_point(self, gas_temperature_k, wall_energy_j, *, context=None):
        """Evaluate shared wall/film physics once and retain its diagnostic facts."""
        wall_temperature = wall_energy_j / self.wall_capacity_j_k
        if any(not math.isfinite(v) or v <= 0 for v in (wall_temperature, gas_temperature_k)):
            raise ValueError('Gas and wall temperatures must be positive.')
        capacity_rate = self._air_capacity_rate
        effective = self._effective_air_conductance
        conductance = self.gas_wall_conductance_w_k
        diagnostics = ()
        if self.gas_film is not None:
            if context is None: raise ValueError('Variable gas film requires instantaneous flow context.')
            conductance, diagnostics = self.gas_film.evaluate(gas_temperature_k,wall_temperature,context)
        air_heat,gas_heat,wall_rate = numeric.wall_heat_rates(effective,self.air_inlet_temperature_k,
            wall_temperature,conductance,gas_temperature_k)
        return WallThermalPoint(gas_heat, air_heat, wall_rate,
            self.air_inlet_temperature_k-air_heat/capacity_rate if capacity_rate > 0 else None,
            wall_temperature, diagnostics)


# Historical name retained for callers and stored study schemas.
from .external_stream import ExternalStreamWallMachine as AirWallMotor
