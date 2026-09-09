"""Optional wall-only initial-guess acceleration for periodic fixed points."""
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
