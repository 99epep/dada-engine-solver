"""Versioned single-phase (rho,u) tables; no extrapolation or property callbacks."""
from dataclasses import dataclass, asdict
import hashlib
import json
import math
import numpy as np
from .fluids import FluidState, CaloricallyPerfectGas

PROPERTIES = ('temperature', 'pressure', 'specific_enthalpy', 'compressibility_factor', 'cp', 'cv', 'sound_speed')


class FluidDomainError(ValueError):
    """Unavailable thermodynamic state, including invalid single-phase cells."""


def interpolate_state(rho, u, rho_axis, u_axis, properties, valid_cells, limits):
    """Bilinear rho/u interpolation, shared by Python and Numba.

    Closed axis endpoints are supported; no clamping/extrapolation is permitted.
    A mask marks complete valid cells, not just sampled valid nodes.
    """
    result = np.zeros(7)
    if not (math.isfinite(rho) and math.isfinite(u)) or rho < rho_axis[0] or rho > rho_axis[-1] or u < u_axis[0] or u > u_axis[-1]:
        return False, result
    i = min(np.searchsorted(rho_axis, rho, side='right')-1, len(rho_axis)-2)
    j = min(np.searchsorted(u_axis, u, side='right')-1, len(u_axis)-2)
    if not valid_cells[i,j]: return False, result
    x = (rho-rho_axis[i])/(rho_axis[i+1]-rho_axis[i])
    y = (u-u_axis[j])/(u_axis[j+1]-u_axis[j])
    for k in range(7):
        a = properties[i,j,k]*(1-y)+properties[i,j+1,k]*y
        b = properties[i+1,j,k]*(1-y)+properties[i+1,j+1,k]*y
        result[k] = a*(1-x)+b*x
    if result[0] < limits[0] or result[0] > limits[1] or result[1] < limits[2] or result[1] > limits[3]:
        return False, result
    return True, result


def _frozen_array(value, dtype=float):
    a = np.asarray(value, dtype=dtype)
    return np.frombuffer(a.tobytes(), dtype=a.dtype).reshape(a.shape)


@dataclass(frozen=True, eq=False)
class TabulatedFluid:
    """Runtime table with immutable storage and a canonical scientific artifact.

    Schema 1 uses positive internal-energy reference and complete-cell validity. This
    is not a phase detector: a data provider must certify its single-phase cells.
    Optional caloric/sound fields are presently required for this prototype.
    """
    rho_axis: object
    u_axis: object
    properties: object
    valid_cells: object
    limits: object  # Tmin, Tmax, Pmin, Pmax
    provenance: str
    ideal_reference: CaloricallyPerfectGas | None = None

    def __post_init__(self):
        if np.asarray(self.valid_cells).dtype.kind != 'b': raise ValueError('Single-phase cell mask requires explicit booleans.')
        for name in ('rho_axis','u_axis','properties','limits','valid_cells'):
            object.__setattr__(self, name, _frozen_array(getattr(self,name), bool if name=='valid_cells' else float))
        for axis in (self.rho_axis,self.u_axis):
            if axis.ndim!=1 or len(axis)<2 or not np.all(np.isfinite(axis)) or axis[0]<=0 or np.any(np.diff(axis)<=0):
                raise ValueError('Table axes must be finite, positive and strictly increasing.')
        if self.properties.shape!=(len(self.rho_axis),len(self.u_axis),7) or not np.all(np.isfinite(self.properties)) or np.any(self.properties<=0):
            raise ValueError('Table requires seven positive finite property fields.')
        if self.valid_cells.shape!=(len(self.rho_axis)-1,len(self.u_axis)-1): raise ValueError('Invalid single-phase cell mask.')
        if self.limits.shape!=(4,) or not np.all(np.isfinite(self.limits)) or min(self.limits)<=0 or self.limits[0]>=self.limits[1] or self.limits[2]>=self.limits[3]:
            raise ValueError('Invalid T/P validity domain.')
        if not isinstance(self.provenance,str) or not self.provenance.strip(): raise ValueError('Table provenance is required.')
        if self.ideal_reference is not None:
            if type(self.ideal_reference) is not CaloricallyPerfectGas: raise ValueError('Invalid ideal reference.')
            expected = _ideal_properties(self.ideal_reference,self.rho_axis,self.u_axis)
            if not np.allclose(self.properties,expected,rtol=2e-14,atol=0):
                raise ValueError('Table does not reproduce its declared analytic ideal reference; ideal hydraulics forbidden.')

    def to_data(self):
        return dict(schema_version=1,model='tabulated_single_phase_rho_u',interpolation='bilinear_rho_u_v1',
            fields=list(PROPERTIES),rho_axis=self.rho_axis.tolist(),u_axis=self.u_axis.tolist(),
            properties=self.properties.tolist(),valid_cells=self.valid_cells.tolist(),
            limits=self.limits.tolist(),provenance=self.provenance,
            ideal_reference=asdict(self.ideal_reference) if self.ideal_reference else None)

    @property
    def content_hash(self):
        return hashlib.sha256(json.dumps(self.to_data(),sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

    @property
    def identity(self):
        return dict(model='tabulated_single_phase_rho_u',version=1,table_sha256=self.content_hash,
            interpolation='bilinear_rho_u_v1',rho_domain=self.rho_axis[[0,-1]].tolist(),
            u_domain=self.u_axis[[0,-1]].tolist(),temperature_pressure_limits=self.limits.tolist(),
            validity='explicit_single_phase_cell_mask; no extrapolation')

    @classmethod
    def from_data(cls,data):
        d=dict(data)
        expected={'schema_version','model','interpolation','fields','rho_axis','u_axis','properties','valid_cells','limits','provenance','ideal_reference'}
        if set(d)!=expected or d.pop('schema_version')!=1 or d.pop('model')!='tabulated_single_phase_rho_u' or d.pop('interpolation')!='bilinear_rho_u_v1' or d.pop('fields')!=list(PROPERTIES):
            raise ValueError('Unsupported fluid table schema.')
        if d['ideal_reference'] is not None: d['ideal_reference']=CaloricallyPerfectGas(**d['ideal_reference'])
        return cls(**d)

    def state_from_rho_u(self,density,specific_energy):
        ok,p=interpolate_state(density,specific_energy,self.rho_axis,self.u_axis,self.properties,self.valid_cells,self.limits)
        if not ok: raise FluidDomainError('State outside declared single-phase rho/u/T/P table domain.')
        return FluidState(density,specific_energy,*p)

    def state_from_rho_t(self,density,temperature):
        # Initialization only. Runtime RHS uses direct rho/u lookup.
        if self.ideal_reference:
            return self.state_from_rho_u(density,self.ideal_reference.heat_capacity_cv*temperature)
        raise NotImplementedError('Supply a validated rho/T inversion for initial filling of this table.')

    def density_from_pt(self,pressure,temperature):
        if self.ideal_reference:
            rho=self.ideal_reference.density_from_pt(pressure,temperature)
            self.state_from_rho_t(rho,temperature)
            return rho
        raise NotImplementedError('Supply a validated P/T inversion for initial filling of this table.')

    # These are ideal-reference properties for explicitly ideal-only diagnostics.
    # Access on a non-ideal table is a clear rejection, never an approximation.
    def _ideal(self):
        from .fluids import require_ideal_hydraulics
        return require_ideal_hydraulics(self)
    gas_constant=property(lambda s:s._ideal().gas_constant)
    heat_capacity_cp=property(lambda s:s._ideal().heat_capacity_cp)
    heat_capacity_cv=property(lambda s:s._ideal().heat_capacity_cv)
    heat_capacity_ratio=property(lambda s:s._ideal().heat_capacity_ratio)
    def compressibility_factor(self,p,t):
        return self.state_from_rho_t(self.density_from_pt(p,t),t).compressibility_factor


def _ideal_properties(gas,rhos,energies):
    return np.array([[[getattr(gas.state_from_rho_u(float(r),float(u)),key) for key in PROPERTIES]
                     for u in energies] for r in rhos])


def ideal_validation_table(gas,rho_axis=(.1,1.,5.,20.),temperature_axis=(200.,300.,450.,650.,1000.)):
    r=np.asarray(rho_axis,float);u=gas.heat_capacity_cv*np.asarray(temperature_axis,float)
    p=_ideal_properties(gas,r,u)
    return TabulatedFluid(r,u,p,np.ones((len(r)-1,len(u)-1),bool),
        (p[:,:,0].min(),p[:,:,0].max(),p[:,:,1].min(),p[:,:,1].max()),
        'Analytic CaloricallyPerfectGas validation data; not a real-gas or cryogenic dataset.',gas)
