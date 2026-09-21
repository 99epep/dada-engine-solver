"""Optional wall and mass-conserving initial guesses for periodic fixed points."""
from dataclasses import dataclass
import math
import numpy as np


def extrapolate_wall_energies(first, second, third, capacities):
    """Componentwise Aitken extrapolation, bounded to a 30 K guess change.

    Applied only between cycles, never as a physical state change within a cycle.
    A subsequent unmodified cycle must satisfy the original periodic criterion.
    """
    first, second, third, capacities = (np.asarray(v,dtype=float) for v in (first,second,third,capacities))
    if any(v.shape != (2,) for v in (first,second,third,capacities)):
        raise ValueError('Expected two wall values per vector.')
    if any(not np.all(np.isfinite(v)) for v in (first,second,third,capacities)) or np.any(capacities<=0):
        raise ValueError('Expected finite values and positive capacities.')
    step=third-second
    previous=second-first
    denominator=step-previous
    safe=(np.abs(denominator)>1e-10)&(step*previous>0)&(np.abs(step)<np.abs(previous))
    correction=np.divide(-step*step,denominator,out=np.zeros(2),where=safe)
    correction=np.clip(correction,-30*capacities,30*capacities)
    proposal=third+correction
    if not np.any(safe) or np.any(proposal<=0):
        return None
    return proposal


@dataclass(frozen=True)
class AdaptiveWallAccelerationSettings:
    """Opt-in numerical heuristics; never physical changes or looser tolerances.

    Four retained endpoints supply three increments and two contraction ratios.
    Fixed-cycle acceleration remains the default when this option is absent.
    """
    damping: float = 0.99
    minimum_confidence: float = 0.95
    maximum_contraction_ratio: float = 0.98
    relative_ratio_drift_limit: float = 0.5
    maximum_correction_kelvin: float = 30.0
    maximum_recent_step_multiple: float = 20.0
    minimum_increment_kelvin: float = 1e-7
    maximum_residual_growth: float = 2.0
    recovery_cycles: int = 4

    def __post_init__(self):
        values = (self.damping, self.minimum_confidence, self.maximum_contraction_ratio,
            self.relative_ratio_drift_limit, self.maximum_correction_kelvin,
            self.maximum_recent_step_multiple, self.minimum_increment_kelvin,
            self.maximum_residual_growth)
        if any(not math.isfinite(v) or v <= 0 for v in values):
            raise ValueError('Adaptive wall controls must be finite and positive.')
        if self.damping > 1 or self.minimum_confidence > 1 or self.maximum_contraction_ratio >= 1:
            raise ValueError('Damping and confidence must not exceed one; contraction must be below one.')
        if self.relative_ratio_drift_limit > 1 or self.maximum_residual_growth < 1:
            raise ValueError('Ratio drift must not exceed one and residual growth must be at least one.')
        if type(self.recovery_cycles) is not int or self.recovery_cycles < 4:
            raise ValueError('Recovery requires at least four ordinary physical cycles.')


@dataclass(frozen=True)
class AdaptiveWallProposal:
    wall_energies: tuple[float, float]
    correction_kelvin: tuple[float, float]
    previous_contraction_ratio: tuple[float, float]
    latest_contraction_ratio: tuple[float, float]
    confidence: tuple[float, float]
    eligible: tuple[bool, bool]


def adaptive_wall_proposal(endpoints, capacities, *,
        settings=AdaptiveWallAccelerationSettings()):
    """Damped Aitken proposal sensitive to changing contraction and recent speed.

    A trend drifting by a large fraction of its remaining contraction margin
    is rejected. Nearly constant, oscillating, growing or negligible increments
    are not extrapolated. Every jump is bounded by recent motion AND 30 K by
    default. The caller must validate the proposal with an ordinary ODE cycle.
    """
    energies=np.asarray(endpoints,dtype=float); capacities=np.asarray(capacities,dtype=float)
    if energies.shape!=(4,2) or capacities.shape!=(2,):
        raise ValueError('Expected four two-wall endpoints and two capacities.')
    if (not np.all(np.isfinite(energies)) or not np.all(np.isfinite(capacities))
            or np.any(energies<=0) or np.any(capacities<=0)):
        raise ValueError('Wall energies and capacities must be finite and positive.')
    temperatures=energies/capacities
    increments=np.diff(temperatures,axis=0)
    first,previous,latest=increments
    significant=np.all(np.abs(increments)>settings.minimum_increment_kelvin,axis=0)
    r0=np.divide(previous,first,out=np.zeros(2),where=first!=0)
    r1=np.divide(latest,previous,out=np.zeros(2),where=previous!=0)
    eligible=(significant & (r0>0) & (r1>0)
        & (r0<settings.maximum_contraction_ratio) & (r1<settings.maximum_contraction_ratio))
    margin=settings.relative_ratio_drift_limit*np.maximum(1-r1,np.finfo(float).eps)
    confidence=np.clip(1-np.abs(r1-r0)/margin,0,1)
    eligible &= confidence>0
    # A weak estimate on an otherwise extrapolatable wall can distort the
    # coupled thermal state. Wait for a coherent reliable window instead.
    if np.any(eligible) and np.min(confidence[eligible]) < settings.minimum_confidence:
        return None
    correction=np.divide(latest*r1,1-r1,out=np.zeros(2),where=eligible)
    correction *= settings.damping*confidence
    bound=np.minimum(settings.maximum_correction_kelvin,
                     settings.maximum_recent_step_multiple*np.abs(latest))
    correction=np.clip(correction,-bound,bound)
    correction=np.where(eligible,correction,0.)
    proposed=energies[-1]+correction*capacities
    eligible &= (proposed>0) & (np.abs(correction)>settings.minimum_increment_kelvin)
    if not np.any(eligible):
        return None
    correction=np.where(eligible,correction,0.)
    proposed=energies[-1]+correction*capacities
    return AdaptiveWallProposal(tuple(float(x) for x in proposed),
        tuple(float(x) for x in correction),tuple(float(x) for x in r0),
        tuple(float(x) for x in r1),tuple(float(x) for x in confidence),
        tuple(bool(x) for x in eligible))


MASS_INDICES = np.array([0, 2, 4, 6])
MASS_ROUNDOFF_TOLERANCE = 512*np.finfo(float).eps


def _physical_state(state):
    values=np.array(state,dtype=float,copy=True)
    if values.shape!=(10,) or not np.all(np.isfinite(values)) or np.any(values<=0):
        raise ValueError('Expected ten finite positive physical wall states.')
    return values


def mass_compatible(state, total_mass):
    return abs(float(np.sum(state[MASS_INDICES]))-total_mass)<=MASS_ROUNDOFF_TOLERANCE*total_mass


def tangent_coordinates(anchor):
    """The diagnostic's original SVD basis, in relative physical-state units."""
    scales=np.abs(_physical_state(anchor))
    normal=np.zeros(10);normal[MASS_INDICES]=scales[MASS_INDICES]
    _,_,vh=np.linalg.svd(normal[None,:]/np.linalg.norm(normal),full_matrices=True)
    return scales,vh[1:].T


@dataclass(frozen=True)
class MassConservingCoordinates:
    anchor: np.ndarray
    scales: np.ndarray
    basis: np.ndarray
    total_mass: float

    @classmethod
    def from_anchor(cls,state):
        anchor=_physical_state(state);scales,basis=tangent_coordinates(anchor)
        for array in (anchor,scales,basis):array.flags.writeable=False
        return cls(anchor,scales,basis,float(anchor[MASS_INDICES].sum()))

    def project(self,state):
        values=_physical_state(state)
        if not mass_compatible(values,self.total_mass):
            raise ValueError('Historical state is incompatible with the fixed gas charge.')
        return self.basis.T@((values-self.anchor)/self.scales)

    def reconstruct(self,z):
        z=np.asarray(z,dtype=float)
        if z.shape!=(9,) or not np.all(np.isfinite(z)):
            raise ValueError('Expected nine finite reduced coordinates.')
        return self.anchor+self.scales*(self.basis@z)


@dataclass(frozen=True)
class PeriodicMapPair:
    initial_state: np.ndarray
    final_state: np.ndarray
    normalized_error: float

    def __post_init__(self):
        for name in ('initial_state','final_state'):
            values=_physical_state(getattr(self,name));values.flags.writeable=False
            object.__setattr__(self,name,values)
        if not math.isfinite(self.normalized_error) or self.normalized_error<0:
            raise ValueError('Periodic error must be finite and nonnegative.')


@dataclass(frozen=True)
class AndersonAccelerationSettings:
    """Experimental map-pair count and safeguards; no change to physical Phi."""
    memory: int = 3
    damping: float = 0.8
    maximum_coefficient_l1: float = 10.
    maximum_relative_state_correction: float = 0.5
    maximum_residual_growth: float = 2.
    minimum_history: int = 2

    def __post_init__(self):
        if type(self.memory) is not int or not 2<=self.memory<=10:
            raise ValueError('Anderson memory must be an integer from two to ten.')
        if type(self.minimum_history) is not int or not 2<=self.minimum_history<=self.memory:
            raise ValueError('Minimum history must be between two and memory.')
        for name in ('damping','maximum_coefficient_l1','maximum_relative_state_correction','maximum_residual_growth'):
            value=getattr(self,name)
            if isinstance(value,bool) or not isinstance(value,(float,int)) or not math.isfinite(value) or value<=0:
                raise ValueError('Anderson controls must be finite positive numbers.')
        if self.damping>1 or self.maximum_coefficient_l1<1 or self.maximum_residual_growth<1:
            raise ValueError('Damping must not exceed one; coefficient/growth bounds must be at least one.')


@dataclass(frozen=True)
class AndersonProposal:
    state: np.ndarray | None
    statistics: dict


def anderson_proposal(pairs, *, settings=AndersonAccelerationSettings()):
    """Constrained residual least squares, rebased at the latest physical end.

    Memory counts physical map pairs. Rank/conditioning failure drops the oldest
    pair. A sqrt(machine epsilon) singular-value ratio rejects numerically
    useless fits; this fixed linear-algebra safeguard is not a tuning parameter.
    """
    stats=dict(status='skipped_insufficient_history',requested_memory=settings.memory,
        effective_memory=min(len(pairs),settings.memory),damping=settings.damping)
    if len(pairs)<settings.minimum_history:return AndersonProposal(None,stats)
    selected=list(pairs[-settings.memory:])
    try:
        coordinates=MassConservingCoordinates.from_anchor(selected[-1].final_state)
        z=np.column_stack([coordinates.project(p.initial_state) for p in selected])
        g=np.column_stack([coordinates.project(p.final_state) for p in selected])
    except ValueError as error:
        return AndersonProposal(None,dict(stats,status='rejected_preflight',reason=str(error)))
    f=g-z
    while f.shape[1]>=2:
        d=f[:,:-1]-f[:,-1,None]
        try:gamma,_,rank,singular=np.linalg.lstsq(d,-f[:,-1],rcond=None)
        except np.linalg.LinAlgError:
            rank=0;singular=np.array([])
        ratio=float(singular[-1]/singular[0]) if len(singular) and singular[0]>0 else 0.
        stats.update(effective_memory=f.shape[1],least_squares_rank=int(rank),
            singular_values=singular.tolist(),singular_value_ratio=ratio)
        if rank==d.shape[1] and ratio>=np.sqrt(np.finfo(float).eps):break
        f=f[:,1:];z=z[:,1:];g=g[:,1:]
    else:return AndersonProposal(None,dict(stats,status='skipped_rank'))
    alpha=np.r_[gamma,1.-gamma.sum()];l1=float(np.sum(np.abs(alpha)))
    stats.update(coefficients=alpha.tolist(),coefficient_l1=l1)
    if not np.all(np.isfinite(alpha)) or l1>settings.maximum_coefficient_l1:
        return AndersonProposal(None,dict(stats,status='skipped_coefficients'))
    candidate=((1-settings.damping)*z+settings.damping*g)@alpha
    try:
        raw=coordinates.reconstruct(candidate)
        norm=float(np.max(np.abs(raw-coordinates.anchor)/coordinates.scales))
        factor=min(1.,settings.maximum_relative_state_correction/norm) if norm else 1.
        proposed=coordinates.reconstruct(factor*candidate)
        applied=float(np.max(np.abs(proposed-coordinates.anchor)/coordinates.scales))
        stats.update(raw_relative_correction=norm,applied_relative_correction=applied,
            correction_scale=factor,reduced_correction_norm=float(np.linalg.norm(factor*candidate)))
        _physical_state(proposed)
        if not mass_compatible(proposed,coordinates.total_mass):raise ValueError('Proposal changed gas charge.')
        if applied>settings.maximum_relative_state_correction*(1+64*np.finfo(float).eps):
            raise ValueError('Proposal correction exceeded its bound.')
    except (ValueError,ArithmeticError) as error:
        return AndersonProposal(None,dict(stats,status='rejected_preflight',reason=str(error)))
    proposed.flags.writeable=False
    return AndersonProposal(proposed,dict(stats,status='proposed'))
