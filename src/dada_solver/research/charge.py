"""Geometry-derived ideal-gas filling inventory, before integration."""
from dataclasses import replace
import math

import numpy as np
from scipy.optimize import brentq, minimize_scalar

POLICY = 'reference_pressure_at_maximum_total_volume_v1'
METHOD = 'periodic_total_volume_grid8192_breakpoints_derivative_roots_local_refinement_v1'


def maximum_total_volume(model):
    """Numerical full-cycle maximum of the assembled thermodynamic volumes.

    Both cylinders share the same solver angle. Stationary roots, local sampled
    peaks and declared breakpoints are considered. This deterministic numerical
    search is not an analytic global-optimality proof for arbitrary motion laws.
    """
    period = 2*math.pi
    k = model.kinematics
    breaks = getattr(k, 'breakpoint_angles', lambda: ())()
    grid = np.unique(np.r_[np.linspace(0., period, 8193),
                           [float(a) % period for a in breaks]])
    def volume(angle):
        value = float(model.volumes(float(angle)).total)
        if not math.isfinite(value) or value <= 0:
            raise ValueError('Reference gas volume must be finite and positive.')
        return value
    def slope(angle):
        return float(k.small_cylinder_volume_derivative(angle) +
                     k.large_cylinder_volume_derivative(angle))
    volumes = np.array([volume(a) for a in grid])
    if np.ptp(volumes) <= 2e-14*float(np.max(volumes)):
        return float(volumes[0]), 0.
    slopes = np.array([slope(a) for a in grid])
    if not np.all(np.isfinite(slopes)):
        raise ValueError('Reference volume derivatives must be finite.')
    candidates = [(float(a) % period, float(v)) for a,v in zip(grid,volumes)]
    for i in range(len(grid)-1):
        if slopes[i] > 0 and slopes[i+1] < 0:
            angle = brentq(slope, grid[i], grid[i+1], xtol=1e-13)
            candidates.append((angle % period, volume(angle)))
    for i in range(1,len(grid)-1):
        if volumes[i] >= volumes[i-1] and volumes[i] >= volumes[i+1] and (
                volumes[i] > volumes[i-1] or volumes[i] > volumes[i+1]):
            optimum = minimize_scalar(lambda a: -volume(a), bounds=(grid[i-1],grid[i+1]),
                                      method='bounded', options={'xatol':1e-13})
            if not optimum.success: raise ValueError('Reference volume local refinement failed.')
            candidates.append((float(optimum.x) % period, volume(optimum.x)))
    # Deterministic first-angle tie break for flat/symmetric maxima. Return the
    # volume at that angle, not a different candidate's (slightly larger) value.
    maximum = max(v for _,v in candidates)
    angle,value = min((a,v) for a,v in candidates if maximum-v <= 2e-14*maximum)
    return value, angle


def apply_reference_charge(design, reference):
    """Use the production connection path once; no independent HX volume math."""
    assembled = design.build()
    model = getattr(assembled, 'model', assembled)
    volume,angle = maximum_total_volume(model)
    pressure,temperature = reference['pressure_pa'],reference['temperature_k']
    mass = pressure*volume/(model.gas.gas_constant*temperature)
    charge = replace(design.configuration.charge, pressure=None,
                     total_mass=mass, temperature=temperature)
    diagnostics = dict(policy=POLICY, reference_pressure_pa=pressure,
        reference_temperature_k=temperature, reference_total_gas_volume_m3=volume,
        reference_angle_rad=angle, derived_total_mass_kg=mass,
        reference_volume_method=METHOD)
    return replace(design,configuration=replace(design.configuration,charge=charge),
                   charge_diagnostics=diagnostics)
