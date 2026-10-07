"""Deterministic feature-aware refit of periodic targets to production splines.

This is geometric initialization only. It never evaluates thermodynamics.
"""
from dataclasses import dataclass, field, asdict
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from dada_solver.free_kinematics import FreeMotionDefinition, _PeriodicMotion
from dada_solver.campaign.candidate import canonical_json, content_hash
from dada_solver.campaign.history import atomic_json
from .motion_target import MotionTarget, PERIOD, PeriodicTargetSide

COUNT = 15
STEP = PERIOD/COUNT
POLICY_VERSION = 'uniform15_feature_position_phase_v1'


@dataclass(frozen=True)
class RefitPolicy:
    """Numerical fit policy, not a physical domain or thermodynamic objective."""
    phase_samples: int = 32
    dense_samples: int = 1440
    feature_samples: int = 33
    maximum_evaluations: int = 120
    velocity_weight: float = 1e-4
    position_tie_mse: float = 1e-12

    def __post_init__(self):
        for name, minimum in (('phase_samples',8), ('dense_samples',360), ('feature_samples',9), ('maximum_evaluations',30)):
            value = getattr(self, name)
            if type(value) is not int or value < minimum:
                raise ValueError(f'{name} must be an integer >= {minimum}.')
        if not math.isfinite(self.velocity_weight) or not 0 <= self.velocity_weight <= 1e-3:
            raise ValueError('Velocity tie-break weight must be finite and between zero and 1e-3.')
        if not math.isfinite(self.position_tie_mse) or not 0 <= self.position_tie_mse <= 1e-10:
            raise ValueError('Position tie tolerance must be finite and at most 1e-10.')



def _feature_mesh(source, policy):
    background = np.linspace(0., PERIOD, policy.dense_samples, endpoint=False)
    blocks = []
    masses = dict(maximum=2., minimum=2., turnaround=1., rounding=.5,
                  cadence_change=.25, kink=1., periodic_seam=.1)
    for index, event in enumerate(source.events):
        kind, angle = event['kind'], event['angle_rad']
        radius = STEP/4
        rule = 'one quarter of spline spacing'
        if kind == 'rounding' and event['end_angle_rad'] is not None:
            width = (event['end_angle_rad']-angle) % PERIOD
            if width == 0: raise ValueError('A rounding interval must have positive width.')
            radius = min(STEP/8, width/4)
            points = np.r_[angle+np.linspace(0.,width,policy.feature_samples),
                           angle+np.linspace(-radius,0.,5), angle+width+np.linspace(0.,radius,5)]
            rule = 'full directed interval; exterior radius min(spacing/8, interval/4)'
        else:
            if kind == 'kink':
                m = event['metadata']
                width = (m['branch_end_rad']-m['branch_start_rad']) % PERIOD if {'branch_start_rad','branch_end_rad'} <= set(m) else PERIOD
                radius = min(STEP/2, width/8)
                rule = 'min(one half spline spacing, one eighth directed branch length)'
            elif kind == 'cadence_change':
                radius = STEP/8
                rule = 'one eighth of spline spacing'
            points = angle+np.linspace(-radius,radius,policy.feature_samples)
        points = np.r_[points, angle] % PERIOD
        weights = np.full(len(points), masses.get(kind,.1)/len(points))
        # Exact turnarounds receive half of this feature's position weight.
        if kind in ('maximum','minimum','turnaround'):
            weights *= .5
            weights[-1] += masses[kind]/2
        blocks.append(dict(event_index=index, kind=kind, angle_rad=angle,
                          samples=points, weights=weights, radius_rad=radius, rule=rule))
    angles = np.r_[background, *[b['samples'] for b in blocks]]
    weights = np.r_[np.full(len(background),1/len(background)), *[b['weights'] for b in blocks]]
    return angles, weights/weights.sum(), blocks


class _FitMotion:
    def __init__(self, definition):
        self.definition = definition
        self.production = _PeriodicMotion(definition)
        self.diagnostics = self.production.diagnostics

    def evaluate(self, angle, order=0):
        values = self.production.evaluate(angle,order)
        return values-1. if order == 0 else values

    def stationary_points(self):
        return tuple(dict(p,volume=p['volume']-1.) for p in self.production.stationary_points())


def _motion(coordinates):
    return _FitMotion(FreeMotionDefinition.from_shape_coordinates(coordinates,1.,2.))


def _simple_topology(motion):
    try:
        points = motion.stationary_points()
    except ValueError:
        return False
    return len(points) == 2 and {p['kind'] for p in points} == {'maximum','minimum'}


def _terms(motion, phase, source, angles, weights, policy):
    error = motion.evaluate(angles-phase)-source.position(angles)
    position = float(weights @ (error**2))
    velocity = None
    if source.has_velocity:
        dv = STEP*(motion.evaluate(angles-phase,1)-source.velocity(angles))
        velocity = float(weights @ (dv**2))
    return dict(weighted_position_mse=position, normalized_velocity_mse=velocity,
                velocity_term=0. if velocity is None else policy.velocity_weight*velocity,
                selection_score=position+(0. if velocity is None else policy.velocity_weight*velocity))


def _diagnostics(motion, phase, source, policy, angles, weights, blocks):
    dense = np.linspace(0.,PERIOD,4*policy.dense_samples,endpoint=False)
    error = motion.evaluate(dense-phase)-source.position(dense)
    dv = None if not source.has_velocity else motion.evaluate(dense-phase,1)-source.velocity(dense)
    points = [dict(angle_rad=(p['angle_rad']+phase) % PERIOD, kind=p['kind'], position=p['volume'])
              for p in motion.stationary_points()]
    expected = source.extrema_angles()
    angular = {p['kind']: float((p['angle_rad']-expected[p['kind']]+math.pi) % PERIOD-math.pi)
               for p in points if p['kind'] in expected}
    features = []
    for block in blocks:
        local = block['samples']
        e = motion.evaluate(local-phase)-source.position(local)
        value = motion.evaluate(block['angle_rad']-phase)-source.position(block['angle_rad'])
        features.append(dict(event_index=block['event_index'], kind=block['kind'], angle_rad=block['angle_rad'],
            radius_rad=block['radius_rad'], sampling_rule=block['rule'], sample_count=len(local),
            position_weight_mass=float(block['weights'].sum()), position_rms=float(np.sqrt(np.mean(e**2))),
            maximum_absolute_position_error=float(np.max(abs(e))), position_error_at_event=float(value),
            velocity_error_at_event=None if not source.has_velocity else float(motion.evaluate(block['angle_rad']-phase,1)-source.velocity(block['angle_rad']))))
    feature_errors = np.concatenate([motion.evaluate(b['samples']-phase)-source.position(b['samples']) for b in blocks]) if blocks else None
    return dict(position_rms=float(np.sqrt(np.mean(error**2))), maximum_absolute_position_error=float(np.max(abs(error))),
                feature_position_rms=None if feature_errors is None else float(np.sqrt(np.mean(feature_errors**2))),
                velocity_rms=None if dv is None else float(np.sqrt(np.mean(dv**2))), extrema=points,
                extrema_count=len(points), extrema_angular_errors_rad=angular,
                maximum_absolute_first_derivative=motion.diagnostics.maximum_absolute_first_derivative,
                maximum_absolute_second_derivative=motion.diagnostics.maximum_absolute_second_derivative,
                events=features, terms=_terms(motion,phase,source,angles,weights,policy),
                error_evaluation=dict(method='uniform periodic samples of frozen target interpolant',samples=len(dense)),
                target_acceleration_used=False)


def _fit_side(target, side, policy):
    source = PeriodicTargetSide(target,side)
    angles, weights, blocks = _feature_mesh(source,policy)
    reference = source.position(angles)
    sqrt_weights = np.sqrt(weights)
    naive_definition = FreeMotionDefinition(tuple(source.position(np.arange(COUNT)*STEP)),1.,2.)
    naive_motion = _FitMotion(naive_definition)
    candidates = []
    exploration = []
    for phase in np.linspace(0.,STEP,policy.phase_samples,endpoint=False):
        seed = FreeMotionDefinition(tuple(source.position(phase+np.arange(COUNT)*STEP)),1.,2.).to_shape_coordinates()
        def residual(z): return sqrt_weights*(_motion(z).evaluate(angles-phase)-reference)
        fit = least_squares(residual,seed,method='lm',ftol=1e-11,xtol=1e-11,gtol=1e-11,max_nfev=policy.maximum_evaluations)
        motion = _motion(fit.x)
        valid = _simple_topology(motion)
        terms = _terms(motion,float(phase),source,angles,weights,policy)
        exploration.append(dict(phase_rad=float(phase),topology_valid=valid,converged=bool(fit.success),terms=terms))
        if valid: candidates.append((motion,float(phase),fit.x))
    if not candidates:
        raise ValueError(f'No topology-preserving 15-control fit found for {side}; no study was produced.')
    def choose(rows):
        scores = [_terms(m,p,source,angles,weights,policy) for m,p,_ in rows]
        best = min(s['weighted_position_mse'] for s in scores)
        eligible = [i for i,s in enumerate(scores) if s['weighted_position_mse'] <= best+policy.position_tie_mse]
        return rows[min(eligible,key=lambda i:(scores[i]['selection_score'],i))]
    selected = choose(candidates)
    motion, phase, z = selected
    spacing = STEP/policy.phase_samples
    def residual_pair(x): return sqrt_weights*(_motion(x[:-1]).evaluate(angles-x[-1])-reference)
    polish = least_squares(residual_pair, np.r_[z,phase],
        bounds=(np.r_[np.full(COUNT-2,-np.inf),phase-spacing],np.r_[np.full(COUNT-2,np.inf),phase+spacing]),
        method='dogbox',ftol=1e-12,xtol=1e-12,gtol=1e-12,max_nfev=policy.maximum_evaluations)
    polished = _motion(polish.x[:-1])
    if _simple_topology(polished): candidates.append((polished,float(polish.x[-1]),polish.x[:-1]))
    motion, phase, _ = choose(candidates)
    # Moving by whole nodes changes only the control indexing, not the motion.
    shift = math.floor(phase/STEP)
    controls = np.roll(motion.definition.control_values,shift)
    phase -= shift*STEP
    definition = FreeMotionDefinition(tuple(controls),1.,2.)
    motion = _FitMotion(definition)
    coordinates = definition.to_shape_coordinates()
    reconstructed = FreeMotionDefinition.from_shape_coordinates(coordinates,1.,2.)
    if not np.allclose(reconstructed.control_values,definition.control_values,rtol=2e-13,atol=2e-13):
        raise ValueError('Fitted control chart round-trip failed.')
    final = _diagnostics(motion,phase,source,policy,angles,weights,blocks)
    naive = _diagnostics(naive_motion,0.,source,policy,angles,weights,blocks)
    return dict(controls=list(definition.control_values),shape_coordinates=list(coordinates),phase_rad=phase,
                node_angles_rad=((phase+np.arange(COUNT)*STEP) % PERIOD).tolist(), diagnostics=final,
                naive=dict(controls=list(naive_definition.control_values),phase_rad=0.,diagnostics=naive),
                phase_exploration=exploration,polish_converged=bool(polish.success))


@dataclass(frozen=True)
class MotionRefitResult:
    payload_json: str

    @property
    def data(self): return json.loads(self.payload_json)
    @property
    def content_hash(self): return self.data['content_hash']

    def save(self,path):
        path = Path(path)
        if path.exists(): raise ValueError('Refit report already exists; choose a new path.')
        path.parent.mkdir(parents=True,exist_ok=True)
        atomic_json(path,self.data)

    @classmethod
    def load(cls,path):
        return cls.from_data(json.loads(Path(path).read_text()))

    @classmethod
    def from_data(cls,data):
        if (set(data) != {'artifact_type','schema_version','scientific','content_hash'}
                or data['artifact_type'] != 'motion_refit' or type(data['schema_version']) is not int or data['schema_version'] != 1
                or content_hash(data['scientific']) != data['content_hash']):
            raise ValueError('Invalid refit report or content hash.')
        return cls(canonical_json(data))


@dataclass(frozen=True)
class MotionRefitRequest:
    target: MotionTarget
    destination_family: str = 'free_spline'
    points_per_piston: int = COUNT
    position_role: str = 'primary_synthesis_reference'
    acceleration_role: str = 'diagnostic_only'
    policy: RefitPolicy = field(default_factory=RefitPolicy)

    def __post_init__(self):
        if not isinstance(self.target,MotionTarget): raise ValueError('Refit requires a MotionTarget.')
        if self.destination_family != 'free_spline' or type(self.points_per_piston) is not int or self.points_per_piston != COUNT:
            raise ValueError('Refit requires free_spline with 15 points per piston.')
        if (self.position_role,self.acceleration_role) != ('primary_synthesis_reference','diagnostic_only'):
            raise ValueError('Refit uses position as reference and acceleration as diagnostic only.')
        if not isinstance(self.policy,RefitPolicy): raise ValueError('Refit requires an explicit valid numerical policy.')

    def execute(self):
        sides = {side:_fit_side(self.target,side,self.policy) for side in ('small','large')}
        scientific = dict(target_hash=self.target.content_hash,target=self.target.data,
            destination_family='free_spline',representation='shape_coordinates',count=COUNT,
            angle_domain=self.target.scientific['angle_domain'],policy_version=POLICY_VERSION,
            policy=asdict(self.policy),position_role=self.position_role,acceleration_role=self.acceleration_role,
            selection_policy='minimum weighted position MSE; velocity score only within numerical position tie',sides=sides,
            combined=dict(position_rms=math.sqrt(sum(s['diagnostics']['position_rms']**2 for s in sides.values())/2),
                          maximum_absolute_position_error=max(s['diagnostics']['maximum_absolute_position_error'] for s in sides.values()),
                          topology_valid=all(s['diagnostics']['extrema_count']==2 for s in sides.values())),
            thermodynamic_evaluation=False)
        return MotionRefitResult(canonical_json(dict(artifact_type='motion_refit',schema_version=1,
                scientific=scientific,content_hash=content_hash(scientific))))


def plot_refit(result):
    """Two-piston static inspection with source events and shifted uniform nodes."""
    import matplotlib.pyplot as plt
    raw = result.data['scientific']
    target = MotionTarget.from_data(raw['target'])
    angles = np.linspace(0.,PERIOD,1441)
    figure, axes = plt.subplots(2,3,figsize=(15,8),sharex=True)
    colors = dict(maximum='tab:red',minimum='tab:purple',rounding='tab:green',
                  cadence_change='tab:orange',kink='tab:brown',periodic_seam='gray')
    for row,side in enumerate(('small','large')):
        fitted = raw['sides'][side]
        phase = fitted['phase_rad']
        source = PeriodicTargetSide(target,side)
        motion = _motion(fitted['shape_coordinates'])
        naive = _FitMotion(FreeMotionDefinition(tuple(fitted['naive']['controls']),1.,2.))
        q = motion.evaluate(angles-phase)
        q0 = naive.evaluate(angles)
        reference = source.position(angles)
        axes[row,0].plot(angles,reference,'k--',label='Target')
        axes[row,0].plot(angles,q,label='Refit')
        axes[row,0].plot(angles,q0,color='gray',alpha=.4,label='Naive phase zero')
        nodes = np.asarray(fitted['node_angles_rad'])
        axes[row,0].scatter(nodes,motion.evaluate(nodes-phase),s=24,marker='o',facecolors='none',edgecolors='tab:blue',label='15 spline nodes',zorder=5)
        points = fitted['diagnostics']['extrema']
        axes[row,0].scatter([p['angle_rad'] for p in points],[p['position'] for p in points],marker='x',c='tab:red',label='Spline extrema',zorder=6)
        axes[row,1].plot(angles,q-reference,label='Refit error')
        axes[row,1].plot(angles,q0-reference,color='gray',alpha=.6,label='Naive error')
        axes[row,2].plot(angles,motion.evaluate(angles-phase,1),label='Refit velocity')
        if source.has_velocity: axes[row,2].plot(angles,source.velocity(angles),'k--',label='Target velocity')
        else: axes[row,2].text(.03,.95,'Target velocity unavailable',transform=axes[row,2].transAxes,va='top')
        event_labels = {}
        for event in source.events:
            if event['kind'] not in colors: continue
            event_labels.setdefault(round(event['angle_rad'],10), []).append(event['kind'])
            for ax in axes[row]:
                ax.axvline(event['angle_rad'],color=colors[event['kind']],alpha=.25,lw=.8)
                if event['kind'] == 'rounding' and event['end_angle_rad'] is not None:
                    start,end = event['angle_rad'],event['end_angle_rad']
                    intervals = [(start,end)] if end >= start else [(start,PERIOD),(0.,end)]
                    for left,right in intervals: ax.axvspan(left,right,color='tab:green',alpha=.05)
        for angle,kinds in event_labels.items():
            axes[row,0].annotate(', '.join(dict.fromkeys(kinds)),xy=(angle,1.02),rotation=90,
                                fontsize=6,va='bottom',color='dimgray')
        for ax,label in zip(axes[row],('Position / stroke','Position error / stroke','Velocity / rad')):
            ax.set_xlim(0.,PERIOD); ax.set_ylabel(label); ax.grid(alpha=.2); ax.legend(fontsize=8,loc='best')
        d,n = fitted['diagnostics'],fitted['naive']['diagnostics']
        axes[row,0].set_title(f"{side.upper()}: RMS {n['position_rms']:.4g} → {d['position_rms']:.4g}; phase {math.degrees(phase):.3f}°",pad=65)
    for ax in axes[-1]: ax.set_xlabel('Study angle (rad)')
    figure.suptitle('Uniform 15-control refit — events guide initialization; no thermodynamic evaluation',fontsize=11)
    figure.tight_layout(rect=(0,0,1,.94))
    return figure
