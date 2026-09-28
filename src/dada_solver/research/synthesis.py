"""Declarative synthesis handoff, not a new optimizer or a generic curve fitter.

See MECHANISM_SYNTHESIS_SEARCH_THEORY.md. Thermodynamic releases execute through
OptimizationCampaign; geometric fitting and fresh-island saturation are deferred.
"""
from dataclasses import dataclass
import math
from .families import PRIMARY_COORDINATES, DOWNSTREAM_COORDINATES, SIXBAR_CONTINUOUS
from .margins import validate_mechanical_constraint

STAGES=('primary_discovery','downstream_fit','full_local_polish','mirror_initialization',
        'opposite_local_adaptation','paired_thermodynamic','hardware_retuning')


def release_coordinates(stage, sides=('large',)):
    """Six-bar coordinate groups only; no target-fitting search is executed."""
    if stage not in STAGES: raise ValueError('Unknown hierarchical synthesis stage.')
    if not sides or len(set(sides))!=len(sides) or set(sides)-{'small','large'}: raise ValueError('Choose distinct cylinder sides.')
    if stage=='paired_thermodynamic' and set(sides)!={'small','large'}: raise ValueError('Paired release requires both independent sides.')
    names=PRIMARY_COORDINATES if stage=='primary_discovery' else DOWNSTREAM_COORDINATES if stage=='downstream_fit' else SIXBAR_CONTINUOUS
    if stage in ('mirror_initialization','hardware_retuning'): return ()
    return tuple(f'kinematics.{side}.{name}' for side in sides for name in names)


@dataclass(frozen=True)
class FreshIslandPolicy:
    """Explicit future saturation inputs, separate from seeded design searches."""
    bounds: dict
    independent_seeds: tuple
    search_policy: str
    primary_branches: tuple = (-1,1)

    def __post_init__(self):
        if not self.bounds or not self.search_policy: raise ValueError('Saturation needs fixed bounds and an explicit search policy.')
        for interval in self.bounds.values():
            if len(interval)!=2 or any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) for v in interval) or interval[0]>=interval[1]:
                raise ValueError('Invalid fixed saturation bounds.')
        if not self.independent_seeds or len(set(self.independent_seeds))!=len(self.independent_seeds) or any(type(s) is not int or s<0 for s in self.independent_seeds):
            raise ValueError('Use distinct deterministic nonnegative island seeds.')
        if len(self.primary_branches)!=2 or any(type(b) is not int for b in self.primary_branches) or set(self.primary_branches)!={-1,1}:
            raise ValueError('Fresh primary islands must balance the two discrete branches.')


@dataclass(frozen=True)
class SynthesisRequest:
    target_study_id: str
    mechanism_family: str
    protocol: str
    stages: tuple
    retained_family_ids: tuple = ()
    position_role: str = 'primary_synthesis_reference'
    acceleration_role: str = 'diagnostic_only'
    mirror_policy: str = 'initialization_only'
    objective_after_pairing: str = 'study_thermodynamic_objective'
    mechanical_constraints: tuple = ()
    saturation_policy: FreshIslandPolicy | None = None

    def __post_init__(self):
        if not self.target_study_id: raise ValueError('A target scientific identity is required.')
        if self.mechanism_family not in ('four_bar','six_bar'): raise ValueError('Unsupported synthesis family.')
        if self.protocol not in ('design_exploitation','fresh_island_saturation'): raise ValueError('Choose design or saturation explicitly.')
        if not self.stages or any(s not in STAGES for s in self.stages) or list(self.stages)!=sorted(set(self.stages),key=STAGES.index):
            raise ValueError('Stages must follow the documented hierarchy without duplication.')
        if self.protocol=='fresh_island_saturation' and (self.retained_family_ids or self.stages!=('primary_discovery',)):
            raise ValueError('Saturation uses fresh starts only and primary discovery, never seeded exploitation.')
        if (self.protocol=='fresh_island_saturation')!=(self.saturation_policy is not None):
            raise ValueError('Only fresh-island saturation requires its explicit bounds, seeds and policy.')
        for row in self.mechanical_constraints:
            validate_mechanical_constraint(row,self.mechanism_family,scoped=False)
        if (self.position_role,self.acceleration_role,self.mirror_policy,self.objective_after_pairing)!=('primary_synthesis_reference','diagnostic_only','initialization_only','study_thermodynamic_objective'):
            raise ValueError('Preserve position/topology guidance, diagnostic acceleration, initialization-only mirroring and thermodynamic final assessment.')
