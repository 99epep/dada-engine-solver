"""Family-owned synthesis protocols and separate candidate evidence.

Geometric search operators are not implemented here. Thermodynamic adaptation
uses the ordinary Research evaluator, never a kinematic proxy score.
"""
from dataclasses import dataclass
import math
from types import MappingProxyType

from .families import (PRIMARY_COORDINATES, DOWNSTREAM_COORDINATES,
                       SIXBAR_CONTINUOUS, parameter_specs, side_metrics, build_side)
from .margins import validate_mechanical_constraint

STAGES = ('primary_discovery', 'downstream_fit', 'full_local_polish',
          'mirror_initialization', 'opposite_local_adaptation',
          'paired_thermodynamic', 'hardware_retuning')


@dataclass(frozen=True)
class SynthesisProtocol:
    """Ordered stage/coordinate ownership; categories are never released."""
    family: str
    releases: tuple

    @property
    def stages(self):
        return tuple(stage for stage, _ in self.releases)

    def coordinates(self, stage):
        for name, coordinates in self.releases:
            if name == stage:
                return coordinates
        raise ValueError(f'Stage {stage!r} is not supported by {self.family}.')


def _continuous(family):
    settings = dict(family=family)
    if family == 'four_bar':
        settings['output'] = 'rocker'
    return tuple(name for name, spec in parameter_specs(settings, 'small').items()
                 if spec.kind == 'continuous')


def _direct_protocol(family):
    coordinates = _continuous(family)
    return SynthesisProtocol(family, (
        ('global_discovery', coordinates), ('full_local_polish', coordinates),
        ('paired_thermodynamic', coordinates), ('hardware_retuning', ())))


PROTOCOLS = MappingProxyType({
    'slider_crank': _direct_protocol('slider_crank'),
    'four_bar': _direct_protocol('four_bar'),
    'six_bar': SynthesisProtocol('six_bar', (
        ('primary_discovery', PRIMARY_COORDINATES),
        ('downstream_fit', DOWNSTREAM_COORDINATES),
        ('full_local_polish', SIXBAR_CONTINUOUS), ('mirror_initialization', ()),
        ('opposite_local_adaptation', SIXBAR_CONTINUOUS),
        ('paired_thermodynamic', SIXBAR_CONTINUOUS), ('hardware_retuning', ()))),
})


def synthesis_protocol(family):
    try:
        return PROTOCOLS[family]
    except KeyError:
        raise ValueError(f'Unsupported synthesis family: {family!r}.') from None


def release_coordinates(stage, sides=('large',), *, family='six_bar'):
    """Release a family's continuous coordinates; default retains six-bar calls."""
    names = synthesis_protocol(family).coordinates(stage)
    if not sides or len(set(sides)) != len(sides) or set(sides)-{'small', 'large'}:
        raise ValueError('Choose distinct cylinder sides.')
    if stage == 'paired_thermodynamic' and set(sides) != {'small', 'large'}:
        raise ValueError('Paired release requires both independent sides.')
    return tuple(f'kinematics.{side}.{name}' for side in sides for name in names)


@dataclass(frozen=True)
class FreshIslandPolicy:
    """Explicit primary-loop saturation inputs, separate from seeded design."""
    bounds: dict
    independent_seeds: tuple
    search_policy: str
    primary_branches: tuple = (-1, 1)

    def __post_init__(self):
        if not self.bounds or not self.search_policy:
            raise ValueError('Saturation needs fixed bounds and an explicit search policy.')
        for interval in self.bounds.values():
            if len(interval) != 2 or any(isinstance(v, bool) or not isinstance(v, (int, float))
                                        or not math.isfinite(v) for v in interval) or interval[0] >= interval[1]:
                raise ValueError('Invalid fixed saturation bounds.')
        if (not self.independent_seeds or len(set(self.independent_seeds)) != len(self.independent_seeds)
                or any(type(s) is not int or s < 0 for s in self.independent_seeds)):
            raise ValueError('Use distinct deterministic nonnegative island seeds.')
        if (len(self.primary_branches) != 2 or any(type(b) is not int for b in self.primary_branches)
                or set(self.primary_branches) != {-1, 1}):
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
        if not self.target_study_id:
            raise ValueError('A target scientific identity is required.')
        owner = synthesis_protocol(self.mechanism_family)
        if self.protocol not in ('design_exploitation', 'fresh_island_saturation'):
            raise ValueError('Choose design or saturation explicitly.')
        if (not self.stages or any(s not in owner.stages for s in self.stages)
                or list(self.stages) != sorted(set(self.stages), key=owner.stages.index)):
            raise ValueError('Stages must follow the family protocol without duplication.')
        if self.protocol == 'fresh_island_saturation':
            if (self.mechanism_family == 'slider_crank' or self.retained_family_ids
                    or self.stages != (owner.stages[0],)):
                raise ValueError('Primary-loop saturation uses fresh starts and discovery only.')
        if (self.protocol == 'fresh_island_saturation') != (self.saturation_policy is not None):
            raise ValueError('Only fresh-island saturation requires its explicit bounds, seeds and policy.')
        for row in self.mechanical_constraints:
            validate_mechanical_constraint(row, self.mechanism_family, scoped=False)
        if (self.position_role, self.acceleration_role, self.mirror_policy, self.objective_after_pairing) != (
                'primary_synthesis_reference', 'diagnostic_only', 'initialization_only', 'study_thermodynamic_objective'):
            raise ValueError('Preserve position/topology guidance, diagnostic acceleration, initialization-only mirroring and thermodynamic final assessment.')


@dataclass(frozen=True)
class KinematicFit:
    """Dimensionless position and optional per-radian velocity errors, not power."""
    position_rms: float
    position_maximum_error: float
    velocity_rms_per_rad: float | None
    target_hash: str
    side: str
    method: str = 'trapezoidal_periodic_angle'


def assess_mechanism(artifact, target, side, *, mechanical_samples=1440, volume_limits=None):
    """Return distinct fit/mechanical evidence using production motion and geometry."""
    import numpy as np
    from dataclasses import asdict
    from dada_solver.geometry import CylinderVolumeLimits

    if side not in ('small', 'large'):
        raise ValueError('Choose SMALL or LARGE.')
    if type(mechanical_samples) is not int or mechanical_samples < 360:
        raise ValueError('Mechanical screening requires at least 360 samples.')
    raw = artifact.scientific
    limits = volume_limits if volume_limits is not None else CylinderVolumeLimits(1., 2.)
    law, geometry = build_side(raw['settings'], raw['geometry'], side, limits)
    angles = np.asarray(target.scientific['angles_rad'])
    reference = target.scientific['sides'][side]
    error = np.array([(law.value(float(t))-limits.minimum)/limits.swept for t in angles])-reference['position']
    def rms(values):
        return float(np.sqrt(np.sum(np.diff(angles)*(values[:-1]**2+values[1:]**2)/2)/(2*math.pi)))
    velocity_error = None
    if reference['first_derivative'] is not None:
        velocity_error = rms(np.array([law.value(float(t), 1)/limits.swept for t in angles])-reference['first_derivative'])
    fit = KinematicFit(rms(error), float(np.max(np.abs(error))), velocity_error, target.content_hash, side)
    mechanical = side_metrics(raw['settings'], law, geometry, mechanical_samples)
    if volume_limits is None:
        # No machine volume scale is implied by a portable normalized mechanism.
        for name in ('maximum_absolute_first_derivative', 'maximum_absolute_second_derivative'):
            if name in mechanical:
                mechanical[name] = None
    from .margins import margin_record
    constraints = [margin_record(row['metric'], mechanical.get(row['metric']), row['limit'],
                                 row['relation'], row['unit'], method='production mechanical screen')
                   for row in raw['constraints']]
    return dict(fit=asdict(fit), mechanical=dict(metrics=mechanical, constraints=constraints,
                samples=mechanical_samples, volume_limits_m3=None if volume_limits is None else [limits.minimum, limits.maximum]), thermodynamic=None)


def synthesis_member(family_id, mechanisms, target, *, provenance=None, thermodynamic=None, volume_limits=None):
    """Build a library member retaining per-side evidence, without selecting a winner.

    Thermodynamic evidence, when supplied, is an existing solver record; it is
    not calculated here. Its identity and metrics remain separate from the fit.
    """
    if not mechanisms or set(mechanisms)-{'small', 'large'}:
        raise ValueError('Provide mechanisms keyed by small and/or large.')
    volume_limits = volume_limits or {}
    if set(volume_limits)-set(mechanisms):
        raise ValueError('Volume limits must belong to a submitted piston mechanism.')
    evidence = {side: assess_mechanism(artifact, target, side, volume_limits=volume_limits.get(side))
                for side, artifact in mechanisms.items()}
    if thermodynamic is not None:
        if not {'candidate_id', 'metrics', 'status'} <= set(thermodynamic):
            raise ValueError('Thermodynamic evidence requires an identified solver record.')
    return dict(family_id=family_id, mechanisms={s: a.data for s, a in mechanisms.items()},
                metadata=dict(target_hash=target.content_hash, evidence=evidence,
                              provenance=provenance or {}, thermodynamic=thermodynamic))


def mechanism_catalogue(library):
    """One row per member and piston; input order is retained, never ranked."""
    rows = []
    for member in library.members:
        metadata = member['metadata']
        for side, data in member['mechanisms'].items():
            evidence = metadata.get('evidence', {}).get(side, {})
            rows.append(dict(family_id=member['family_id'], side=side,
                mechanism_family=data['scientific']['settings']['family'], artifact_hash=data['content_hash'],
                fit=evidence.get('fit'), mechanical=evidence.get('mechanical'),
                provenance=dict(family=metadata.get('provenance', {}), mechanism=data['provenance']),
                thermodynamic=metadata.get('thermodynamic')))
    return rows


@dataclass(frozen=True)
class SynthesisPlan:
    """Bind a family protocol to a target and retain submitted physical candidates.

    This boundary validates and records externally obtained geometries; it does
    not invent geometric search operators or run thermodynamics.
    """
    target: object
    request: SynthesisRequest

    def __post_init__(self):
        from .motion_target import MotionTarget
        if not isinstance(self.target, MotionTarget):
            raise ValueError('A synthesis plan requires a MotionTarget.')
        source = self.target.scientific['source']
        if self.request.target_study_id not in (self.target.content_hash, source.get('study_id')):
            raise ValueError('Synthesis request and target scientific identity disagree.')

    def artifact(self, geometry, *, settings=None, provenance=None):
        from .artifacts import MechanismArtifact
        return MechanismArtifact.create(self.request.mechanism_family, geometry,
            settings=settings, constraints=self.request.mechanical_constraints,
            provenance=dict(provenance or {}, target_hash=self.target.content_hash))

    def member(self, family_id, mechanisms, *, provenance=None, thermodynamic=None, volume_limits=None):
        if any(a.scientific['settings']['family'] != self.request.mechanism_family
               for a in mechanisms.values()):
            raise ValueError('Submitted mechanisms must belong to the requested family.')
        return synthesis_member(family_id, mechanisms, self.target,
                                provenance=provenance, thermodynamic=thermodynamic, volume_limits=volume_limits)

    def execute(self):
        raise NotImplementedError('Geometric synthesis search operators are not implemented; validate a plan or submit physical geometries.')
