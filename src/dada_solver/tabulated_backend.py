"""Compiled table reconstruction feeding the same conservative wall balances."""
from functools import lru_cache
import inspect
import math
import numpy as np
from . import numerical_primitives as numeric
from .tabulated_fluid import interpolate_state


def wall_kernel_tabulated(values,volumes,volume_rates,gas,links,walls,sources,destinations,one_way,ports,
                          rho_axis,u_axis,properties,valid_cells,limits):
    t=np.empty(4);pressures=np.empty(4);h=np.empty(4)
    for j in range(4):
        m,energy=values[2*j],values[2*j+1]
        if not math.isfinite(m) or m<=0 or not math.isfinite(volumes[j]) or volumes[j]<=0:
            return False,np.zeros(15)
        ok,state=interpolate_state(m/volumes[j],energy/m,rho_axis,u_axis,properties,valid_cells,limits)
        if not ok: return False,np.zeros(15)
        t[j],pressures[j],h[j]=state[0],state[1],state[2]
    return numeric.wall_balance_kernel(values,volume_rates,gas,links,walls,sources,destinations,
        one_way,ports,t,pressures,h)


@lru_cache(maxsize=1)
def compiled_tabulated_dispatcher():
    from numba import njit
    from numba.extending import register_jitable
    register_jitable(interpolate_state)
    for name,function in vars(numeric).items():
        if inspect.isfunction(function) and function.__module__==numeric.__name__:
            register_jitable(function)
    # No disk cache in this prototype: never risk transitive source invalidation.
    return njit(fastmath=False,parallel=False)(wall_kernel_tabulated)
