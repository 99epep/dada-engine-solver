"""Optional wall-only initial-guess acceleration for periodic fixed points."""
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
