"""Immutable periodic motion targets in study angle, independent of source family."""
from dataclasses import dataclass
from contextlib import contextmanager
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline, CubicHermiteSpline

from dada_solver.campaign.candidate import canonical_json, content_hash
from dada_solver.campaign.history import atomic_json

ANGLE_DOMAIN = 'study_angle_before_operation_transform'
PERIOD = 2*math.pi


@dataclass(frozen=True)
class MotionEvent:
    """Extensible named feature; a wrapped end angle describes a directed interval.

    Metadata can carry exact branch boundaries or source-defined descriptors.
    No default kink-neighborhood width or acceleration fit is inferred.
    """
    kind: str
    angle_rad: float
    end_angle_rad: float | None = None
    metadata_json: str = '{}'

    def to_data(self):
        return dict(kind=self.kind, angle_rad=self.angle_rad,
                    end_angle_rad=self.end_angle_rad, metadata=json.loads(self.metadata_json))


@dataclass(frozen=True)
class MotionTarget:
    payload_json: str

    @property
    def data(self):
        return json.loads(self.payload_json)

    @property
    def scientific(self):
        return self.data['scientific']

    @property
    def content_hash(self):
        return self.data['content_hash']

    @classmethod
    def create(cls, angles_rad, sides, *, source, provenance=None):
        angles = np.asarray(angles_rad, dtype=float)
        if (angles.ndim != 1 or len(angles) < 3 or not np.all(np.isfinite(angles))
                or np.any(np.diff(angles) <= 0) or angles[0] != 0. or angles[-1] != PERIOD):
            raise ValueError('Target angles must increase from 0 through 2*pi, including the periodic endpoint.')
        if set(sides) != {'small', 'large'} or not isinstance(source, dict) or not source:
            raise ValueError('A target needs SMALL/LARGE and explicit scientific source identity.')
        validated = {}
        for side, raw in sides.items():
            if set(raw) != {'position', 'first_derivative', 'second_derivative', 'events'}:
                raise ValueError('Target side must explicitly declare position, optional derivatives and events.')
            out = {}
            for key in ('position', 'first_derivative', 'second_derivative'):
                values = raw[key]
                if values is None and key != 'position':
                    out[key] = None
                    continue
                a = np.asarray(values, dtype=float)
                if a.shape != angles.shape or not np.all(np.isfinite(a)):
                    raise ValueError('Target samples must be finite and match the angle grid.')
                if not math.isclose(float(a[0]), float(a[-1]), rel_tol=1e-9, abs_tol=1e-9):
                    raise ValueError('Target position and available derivatives must be periodic.')
                if key == 'position' and (a.min() < -1e-9 or a.max() > 1+1e-9):
                    raise ValueError('Target position must be normalized to declared stroke; no clipping is performed.')
                out[key] = a.tolist()
            events = []
            for event in raw['events']:
                event = event.to_data() if isinstance(event, MotionEvent) else event
                if (set(event) != {'kind', 'angle_rad', 'end_angle_rad', 'metadata'}
                        or not isinstance(event['kind'], str) or not event['kind']
                        or not isinstance(event['metadata'], dict)):
                    raise ValueError('Invalid structured motion event.')
                for value in (event['angle_rad'], event['end_angle_rad']):
                    if value is not None and (isinstance(value, bool) or not isinstance(value, (int, float))
                            or not math.isfinite(value) or not 0 <= value < PERIOD):
                        raise ValueError('Event angles must lie in [0, 2*pi).')
                if event['angle_rad'] is None:
                    raise ValueError('A motion event requires an angle.')
                events.append(event)
            out['events'] = events
            validated[side] = out
        scientific = dict(angle_domain=ANGLE_DOMAIN, angle_unit='rad', derivative_variable='study_angle_radians',
                          period_rad=PERIOD, normalization='(volume-minimum_volume)/swept_volume',
                          angles_rad=angles.tolist(), sides=validated, source=source)
        return cls(canonical_json(dict(schema_version=1, artifact_type='motion_target',
            scientific=scientific, content_hash=content_hash(scientific), provenance=provenance or {})))

    @classmethod
    def from_data(cls, data):
        if (set(data) != {'schema_version', 'artifact_type', 'scientific', 'content_hash', 'provenance'}
                or type(data['schema_version']) is not int or data['schema_version'] != 1
                or data['artifact_type'] != 'motion_target'):
            raise ValueError('Unsupported motion target schema.')
        raw = data['scientific']
        result = cls.create(raw['angles_rad'], raw['sides'], source=raw['source'], provenance=data['provenance'])
        if result.scientific != raw or result.content_hash != data['content_hash']:
            raise ValueError('Motion target convention or content hash mismatch.')
        return result

    @classmethod
    def load(cls, path):
        return cls.from_data(json.loads(Path(path).read_text()))

    def save(self, path):
        path = Path(path)
        if path.exists():
            raise ValueError('Motion target already exists; choose a new path.')
        path.parent.mkdir(parents=True, exist_ok=True)
        atomic_json(path, self.data)

    @classmethod
    def from_kinematics(cls, model, *, source, samples=721, events=None, provenance=None):
        """Sample exact production values; absent derivatives are not differenced."""
        if type(samples) is not int or samples < 3:
            raise ValueError('At least three samples are required.')
        angles = np.linspace(0., PERIOD, samples)
        sides = {}
        for side in ('small', 'large'):
            limits = getattr(model, side+'_volume_limits')
            out = {}
            for order, key in enumerate(('position', 'first_derivative', 'second_derivative')):
                suffix = ('', '_derivative', '_second_derivative')[order]
                method = getattr(model, side+'_cylinder_volume'+suffix, None)
                try:
                    values = None if method is None else np.array([method(float(t)) for t in angles])
                except (AttributeError, NotImplementedError):
                    values = None
                if order == 0 and values is None:
                    raise ValueError('A motion target requires an available position law.')
                out[key] = None if values is None else ((values-limits.minimum if order == 0 else values)/limits.swept)
            features = [MotionEvent('periodic_seam', 0.)]
            if out['first_derivative'] is not None:
                method = getattr(model, side+'_cylinder_volume_derivative')
                roots = []
                velocity = out['first_derivative']
                for i in range(len(angles)-1):
                    left, right = velocity[i:i+2]
                    if left*right < 0:
                        roots.append(brentq(method, float(angles[i]), float(angles[i+1]), xtol=1e-13) % PERIOD)
                    elif left == 0 and velocity[i-1 if i else -2]*right < 0:
                        roots.append(float(angles[i]) % PERIOD)
                for angle in roots:
                    features.append(MotionEvent('turnaround', angle, metadata_json=canonical_json(
                        dict(method='production_derivative_root_on_declared_grid'))))
            features.extend((events or {}).get(side, ()))
            out['events'] = features
            sides[side] = out
        return cls.create(angles, sides, source=source, provenance=provenance)


def _hybrid_events(model, side):
    """Exact parameter-defined features, not a refit or sampled cadence heuristic."""
    maximum = model.small_max_deg if side == 'small' else 0.
    down = getattr(model, side+'_down_duration_deg')
    minimum = maximum+down
    events = [MotionEvent('maximum', math.radians(-maximum) % PERIOD),
              MotionEvent('minimum', math.radians(-minimum) % PERIOD)]
    if side == 'small':
        start, width, rounding = minimum, 360.-down, model.small_up_rounding
        kink_start, kink_width, kink = maximum, down, model.small_down_kink_u
    else:
        start, width, rounding = maximum, down, model.large_down_rounding
        kink_start, kink_width, kink = minimum, 360.-down, model.large_up_kink_u
    angle = lambda degree: math.radians(-degree) % PERIOD
    # Intervals follow increasing study angle, opposite the model's motor-angle parameter.
    events.extend((MotionEvent('rounding', angle(start+width*rounding), angle(start)),
                   MotionEvent('rounding', angle(start+width), angle(start+width*(1-rounding))),
                   MotionEvent('cadence_change', angle(start+width*rounding)),
                   MotionEvent('cadence_change', angle(start+width*(1-rounding))),
                   MotionEvent('kink', angle(kink_start+kink_width*kink), metadata_json=canonical_json(
                       dict(branch_start_rad=angle(kink_start+kink_width), branch_end_rad=angle(kink_start),
                            neighborhood_policy=None)))))
    limits = getattr(model, side+'_volume_limits')
    position = getattr(model, side+'_cylinder_volume')
    velocity = getattr(model, side+'_cylinder_volume_derivative')
    enriched = []
    for event in events:
        metadata = json.loads(event.metadata_json)
        metadata.update(normalized_position=(position(event.angle_rad)-limits.minimum)/limits.swept,
                        first_derivative_per_rad=velocity(event.angle_rad)/limits.swept)
        enriched.append(MotionEvent(event.kind, event.angle_rad, event.end_angle_rad, canonical_json(metadata)))
    return enriched


def target_from_study(study, active_values=None, *, samples=721):
    from .schema import compile_study, candidate_for_values
    from dada_solver.hybrid_compact_kinematics import HybridCompactKinematics

    definition = compile_study(study)
    values = {p.name: p.initial for p in study.space.parameters} if active_values is None else dict(active_values)
    candidate = candidate_for_values(definition, values)
    design = definition.adapter.build(dict(definition.fixed_parameters, **values))
    events = {}
    for side in ('small', 'large'):
        law = getattr(design.kinematics, side, None)
        backend = getattr(law, 'model', None)
        if isinstance(backend, HybridCompactKinematics):
            events[side] = _hybrid_events(backend, side)
    return MotionTarget.from_kinematics(design.kinematics, samples=samples, events=events,
        source=dict(study_id=study.study_id, candidate_id=candidate.candidate_id,
                    definition_id=definition.definition_id, families=study.settings),
        provenance=dict(extraction='production study-angle kinematics; no thermodynamic integration'))


def load_motion_target(source, *, candidate=None, samples=721):
    """Resolve a target, current study, or identified stored Research result."""
    source = Path(source)
    if source.suffix == '.toml':
        if candidate is not None:
            raise ValueError('A study TOML uses its configured values; candidate selectors require a stored result.')
        from .schema import load_study
        return target_from_study(load_study(source), samples=samples)
    if source.is_file():
        data = json.loads(source.read_text())
        if data.get('artifact_type') == 'motion_target':
            if candidate is not None:
                raise ValueError('A motion target does not accept a candidate selector.')
            return MotionTarget.from_data(data)
    with research_motion_source(source, candidate=candidate) as (study, values, identity):
        target = target_from_study(study, values, samples=samples)
    raw = target.scientific
    raw['source'].update(identity)
    return MotionTarget.create(raw['angles_rad'], raw['sides'], source=raw['source'], provenance=target.data['provenance'])


@contextmanager
def research_motion_source(source, *, candidate=None):
    """Reconstruct a complete current Research input, never invent a machine.

    MotionTarget alone deliberately does not provide this context. Returned
    coordinates are exact configured/selected values, not decoded approximations.
    """
    source = Path(source)
    if source.suffix == '.toml':
        if candidate is not None:
            raise ValueError('Candidate selectors require a stored Research result.')
        from .schema import load_study, compile_study, candidate_for_values
        study = load_study(source)
        values = {p.name:p.initial for p in study.space.parameters}
        definition = compile_study(study)
        chosen = candidate_for_values(definition,values)
        yield study, values, dict(study_id=study.study_id,candidate_id=chosen.candidate_id,
                                 source_definition_id=definition.definition_id)
        return
    if source.is_file() and json.loads(source.read_text()).get('artifact_type') == 'motion_target':
        raise ValueError('A MotionTarget alone cannot generate a machine study; use --report only or a Research source.')
    from .report import inspect, select_records
    from .snapshots import stored_study
    data = inspect(source)
    if candidate is None:
        if source.is_dir() or len(data['records']) != 1:
            raise ValueError('A stored campaign requires an explicit candidate selector.')
        record = data['records'][0]
    else:
        record, = select_records(data,[candidate])
    with stored_study(data) as study:
        study.space.encode(record['physical'])
        yield study, record['physical'], dict(study_id=data['study_id'],candidate_id=record['candidate_id'],
                                             source_definition_id=data['definition_id'])


class PeriodicTargetSide:
    """Periodic interpolation of frozen samples, enriched by exact source events.

    Hermite interpolation uses available source first derivatives. With position
    only, a periodic cubic interpolant is used, but no target velocity is claimed.
    Target acceleration is never read.
    """
    def __init__(self, target, side):
        raw = target.scientific
        self.data = raw['sides'][side]
        self.events = self.data['events']
        self.has_velocity = self.data['first_derivative'] is not None
        samples = {float(t): [float(q), None if not self.has_velocity else float(v)]
                   for t, q, v in zip(raw['angles_rad'], self.data['position'],
                                      self.data['first_derivative'] or [None]*len(raw['angles_rad']))}
        for event in self.events:
            m = event['metadata']
            if 'normalized_position' in m and (not self.has_velocity or 'first_derivative_per_rad' in m):
                t = float(event['angle_rad'])
                near = next((old for old in samples if abs(old-t) < 1e-12), None)
                if near is not None: del samples[near]
                samples[t] = [float(m['normalized_position']), m.get('first_derivative_per_rad')]
        angles = sorted(samples)
        q = [samples[t][0] for t in angles]
        if self.has_velocity:
            self.interpolant = CubicHermiteSpline(angles, q, [samples[t][1] for t in angles])
        else:
            self.interpolant = CubicSpline(angles, q, bc_type='periodic')

    def position(self, angles):
        return self.interpolant(np.asarray(angles) % PERIOD)

    def velocity(self, angles):
        if not self.has_velocity:
            return None
        return self.interpolant(np.asarray(angles) % PERIOD, 1)

    def extrema_angles(self):
        known = {e['kind']:e['angle_rad'] for e in self.events if e['kind'] in ('maximum','minimum')}
        if set(known) == {'maximum','minimum'}: return known
        roots = self.interpolant.derivative().roots(extrapolate=False)
        roots = roots[np.isfinite(roots) & (roots >= 0) & (roots < PERIOD)]
        if not len(roots): raise ValueError('Target has no resolved extrema.')
        values = self.position(roots)
        return dict(maximum=float(roots[np.argmax(values)]), minimum=float(roots[np.argmin(values)]))
