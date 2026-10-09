"""Six-bar stage adapters for the shared GeometrySearch island/polish engine.

No differential-evolution implementation or thermodynamic evaluator lives here.
Intermediate primary evidence is a projection opportunity, never a piston fit.
"""
from dataclasses import asdict, replace, dataclass
import math
import json
import numpy as np

from dada_solver.six_bar import SixBarPrimaryMechanism, SixBarCylinderMechanism, velocity_stationary_points
from .artifacts import MechanismArtifact, MechanismLibrary
from .families import PRIMARY_COORDINATES, DOWNSTREAM_COORDINATES, SIXBAR_CONTINUOUS
from .margins import margin_record
from .motion_target import PERIOD, PeriodicTargetSide
from .synthesis_search import GeometrySearch, SearchPolicy, SearchStopped, BOUNDS, CATEGORY_VALUES

from .mechanical_screen import FINAL_MECHANICAL_SAMPLES, inherit_parent_plan, screen_provenance, PrimaryMechanicalPolicy

STAGE_POLICY = 'hierarchical_six_bar_geometry_v1'
PRIMARY_POLICY = 'primary_topology_cadence_v2'
STAGES = ('primary_discovery','downstream_fit','full_local_polish','mirror_initialization','opposite_local_adaptation')


def primary_cycle_data(primary, samples):
    """Sign-independent full-cycle metrics from production primary joints."""
    angles=np.linspace(0.,PERIOD,samples,endpoint=False)
    state=primary.joint_state(angles)
    points=np.column_stack(state['joints']['E'])
    metrics=dict(minimum_primary_transmission_sine=float(np.min(state['primary_transmission_sine'])),
                 E_span_over_crank=float(np.max(np.ptp(points,axis=0))))
    return state,points,metrics


def primary_rejection_reason(records):
    names={r['name'] for r in records if not r['satisfied']}
    if 'minimum_primary_transmission_sine' in names: return 'primary_transmission_below_minimum'
    if 'E_span_over_crank' in names: return 'primary_e_span_below_minimum'
    return None


class PrimaryProjection:
    """PCA projection of E for inspectable chronology, not a cylinder law."""
    def __init__(self, primary, *, axis=None, sign=1., root_samples=360, metric_samples=721, cycle=None):
        self.geometry=primary
        state,points,base_metrics=cycle if cycle is not None else primary_cycle_data(primary,metric_samples)
        if axis is None:
            _,vectors=np.linalg.eigh(np.cov(points.T))
            axis=vectors[:,-1]
            if axis[np.argmax(abs(axis))]<0: axis=-axis
        self.axis=np.asarray(axis)*sign
        roots=velocity_stationary_points(self.coordinate,self.derivative,root_samples)
        if len(roots)<2:
            raise ValueError('Projected E has no useful chronology.')
        values=[r['position'] for r in roots]
        self.low,self.span=min(values),max(values)-min(values)
        if self.span<=1e-10: raise ValueError('Degenerate primary E projection.')
        self.roots=roots
        self.lateral_rms=float(np.std(points@np.array([-self.axis[1],self.axis[0]])))/self.span
        self.metrics=dict(base_metrics,
                          E_projection_span_over_crank=self.span,E_lateral_rms_over_projection_span=self.lateral_rms,
                          E_projection_extrema_count=len(roots),
                          E_axis_variance_fraction=float(np.var(points@self.axis)/max(np.var(points,axis=0).sum(),1e-30)))

    def coordinate(self,angles):
        state=self.geometry.joint_state(angles)
        return state['joints']['E'][0]*self.axis[0]+state['joints']['E'][1]*self.axis[1]

    def derivative(self,angles):
        state=self.geometry.joint_state(angles)
        return state['E_derivative'][0]*self.axis[0]+state['E_derivative'][1]*self.axis[1]

    def value(self,angles,order=0):
        if order==0: return 1+(self.coordinate(angles)-self.low)/self.span
        if order==1: return self.derivative(angles)/self.span
        raise NotImplementedError('Primary projection acceleration is not a synthesis target.')


@dataclass(frozen=True)
class PrimaryCadencePolicy:
    """Method 6.4 numerical preferences, never universal mechanical limits."""
    monotonicity_weight: float = 2.
    endpoint_zero_weight: float = .40
    turning_weight: float = .65
    turning_scale_rad: float = math.pi/12
    extra_crossing_weight: float = .50
    speed_ratio_weight: float = .55
    displacement_fraction_weight: float = .80
    long_symmetry_weight: float = .45
    transmission_preference_weight: float = .12
    directionality_preference_weight: float = .03
    preferred_transmission_sine: float = .35
    preferred_axis_variance_fraction: float = .70
    sign_guard_fraction: float = 3/360
    turnaround_guard_fraction: float = 8/360
    transition_guard_fraction: float = 6/360
    minimum_subphase_fraction: float = .15
    minimum_speed_ratio: float = 1.5
    minimum_split_explained_variance: float = .5
    missing_cadence: str = 'disable'
    reject_extra_turnarounds: bool = False

    def __post_init__(self):
        for name,value in asdict(self).items():
            if name in ('missing_cadence','reject_extra_turnarounds'): continue
            if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0:
                raise ValueError(f'primary_cadence.{name} must be finite and nonnegative.')
        if self.turning_scale_rad<=0 or self.minimum_speed_ratio<=1:
            raise ValueError('Primary timing scale must be positive and minimum_speed_ratio greater than one.')
        for name in ('sign_guard_fraction','turnaround_guard_fraction','transition_guard_fraction'):
            if not 0<=getattr(self,name)<.125: raise ValueError(f'primary_cadence.{name} must be in [0, 0.125).')
        if not 0<self.minimum_subphase_fraction<.5: raise ValueError('minimum_subphase_fraction must be in (0, 0.5).')
        for name in ('preferred_transmission_sine','preferred_axis_variance_fraction','minimum_split_explained_variance'):
            if not 0<=getattr(self,name)<=1: raise ValueError(f'primary_cadence.{name} must be in [0, 1].')
        if self.missing_cadence not in ('disable','error') or type(self.reject_extra_turnarounds) is not bool:
            raise ValueError('Use missing_cadence=disable/error and a boolean reject_extra_turnarounds.')


def _cyclic_error(a,b):
    return (a-b+math.pi)%PERIOD-math.pi


def primary_target_features(source,policy):
    """Target-owned turns and ordered cadence; absent cadence is explicit."""
    extrema=source.extrema_angles()
    maximum,minimum=extrema['maximum'],extrema['minimum']
    down=(minimum-maximum)%PERIOD
    ascending=down>math.pi
    start=minimum if ascending else maximum
    span=PERIOD-down if ascending else down
    result=dict(maximum_rad=maximum,minimum_rad=minimum,short_start_rad=start,
        short_span_rad=span,long_span_rad=PERIOD-span,short_direction=1 if ascending else -1,
        cadence=dict(available=False,reason='target_velocity_unavailable'))
    if source.has_velocity:
        guard=min(PERIOD*policy.turnaround_guard_fraction,.1*span)
        u=np.linspace(guard,span-guard,257)
        speed=np.abs(source.velocity(start+u))
        total=float(np.sum((speed-speed.mean())**2))
        event_splits=[]
        for event in source.events:
            offset=(event['angle_rad']-start)%PERIOD
            if event['kind'] in ('kink','cadence_change') and policy.minimum_subphase_fraction*span<offset<(1-policy.minimum_subphase_fraction)*span:
                event_splits.append((offset,'event:'+event['kind']))
        splits=event_splits or [(float(x),'two_level_velocity_estimate') for x in np.linspace(
            policy.minimum_subphase_fraction*span,(1-policy.minimum_subphase_fraction)*span,65)]
        candidates=[]
        for split,origin in splits:
            transition=min(PERIOD*policy.transition_guard_fraction,.2*min(split,span-split))
            first=(u<split-transition); second=(u>split+transition)
            if first.sum()<3 or second.sum()<3: continue
            means=[float(np.mean(np.abs(source.velocity(start+np.linspace(lo,hi,65)))))
                   for lo,hi in ((guard,split-transition),(split+transition,span-guard))]
            ratio=max(means)/max(min(means),1e-12)
            # Numerical two-level detection must improve on a single cadence;
            # event-supported transitions still require a resolved speed contrast.
            left=u<=split; right=~left
            sse=float(sum(np.sum((speed[m]-speed[m].mean())**2) for m in (left,right)))
            explained=max(0.,1-sse/max(total,1e-30))
            if ratio<policy.minimum_speed_ratio or (not event_splits and explained<policy.minimum_split_explained_variance): continue
            first_fast=means[0]>=means[1]
            q=source.position(np.array([start,start+split,start+span]))
            displacement=np.abs(np.diff(q));fraction=float(displacement[0 if first_fast else 1]/max(displacement.sum(),1e-12))
            candidates.append((sse,dict(available=True,source=origin,split_offset_rad=split,
                split_angle_rad=(start+split)%PERIOD,first_is_fast=first_fast,fast_slow_speed_ratio=ratio,
                fast_displacement_fraction=fraction,explained_variance=explained,
                first_core_rad=[guard,split-transition],second_core_rad=[split+transition,span-guard])))
        result['cadence']=min(candidates,key=lambda c:c[0])[1] if candidates else dict(
            available=False,reason='no_resolved_two_cadence_transition')
    if not result['cadence']['available'] and policy.missing_cadence=='error':
        raise ValueError('Primary target cadence unavailable: '+result['cadence']['reason'])
    return result


def primary_cadence_score(projection,target_features,policy,samples=721):
    """Signed topology/cadence score; no target position or velocity-profile fit."""
    angles=np.linspace(0.,PERIOD,samples,endpoint=False)
    velocity=np.asarray(projection.value(angles,1)); scale=max(float(np.sqrt(np.mean(velocity**2))),1e-12)
    start=target_features['short_start_rad'];span=target_features['short_span_rad'];direction=target_features['short_direction']
    guard=min(PERIOD*policy.sign_guard_fraction,.1*span,.1*(PERIOD-span))
    u=(angles-start)%PERIOD
    short=(u>=guard)&(u<=span-guard); long=(u>=span+guard)&(u<=PERIOD-guard)
    wrong=np.r_[np.maximum(-direction*velocity[short],0.),np.maximum(direction*velocity[long],0.)]
    monotonicity=float(np.sqrt(np.mean(wrong**2))/scale)
    roots=projection.roots
    high=[r for r in roots if r['kind']=='maximum'];low=[r for r in roots if r['kind']=='minimum']
    if not high or not low: raise ValueError('Projected E lacks a maximum/minimum pair.')
    high_turn=min(high,key=lambda r:abs(_cyclic_error(r['angle_rad'],target_features['maximum_rad'])))['angle_rad']
    low_turn=min(low,key=lambda r:abs(_cyclic_error(r['angle_rad'],target_features['minimum_rad'])))['angle_rad']
    turning=float(np.sqrt(np.mean([_cyclic_error(high_turn,target_features['maximum_rad'])**2,
                                  _cyclic_error(low_turn,target_features['minimum_rad'])**2])))
    endpoints=np.asarray(projection.value(np.array([target_features['maximum_rad'],target_features['minimum_rad']]),1))
    zero=float(np.sqrt(np.mean(endpoints**2))/scale)
    extra=max(0,len(high)+len(low)-2)
    if policy.reject_extra_turnarounds and extra: raise ValueError('Additional primary projection turnarounds are excluded by policy.')
    cadence=target_features['cadence'];ratio=fraction=None;ratio_error=fraction_error=0.
    if cadence['available']:
        means=[]
        for lower,upper in (cadence['first_core_rad'],cadence['second_core_rad']):
            means.append(float(np.mean(np.abs(projection.value(start+np.linspace(lower,upper,65),1)))))
        fast=0 if cadence['first_is_fast'] else 1
        ratio=means[fast]/max(means[1-fast],1e-12)
        ratio_error=abs(math.log(max(ratio,1e-12)/cadence['fast_slow_speed_ratio']))
        q=projection.value(start+np.array([0.,cadence['split_offset_rad'],span]))
        displacement=np.abs(np.diff(q));fraction=float(displacement[fast]/max(displacement.sum(),1e-12))
        fraction_error=abs(fraction-cadence['fast_displacement_fraction'])
    own_down_span=(low_turn-high_turn)%PERIOD
    own_start=low_turn if own_down_span<=math.pi else high_turn
    own_end=high_turn if own_down_span<=math.pi else low_turn
    own_span=(own_end-own_start)%PERIOD
    distance=np.linspace(min(PERIOD*policy.turnaround_guard_fraction,.20*own_span),.5*own_span,64)
    a=np.asarray(projection.value(own_start+distance,1));b=np.asarray(projection.value(own_end-distance,1))
    asymmetry=float(np.sqrt(np.mean((a-b)**2))/max(float(np.sqrt(np.mean(.5*(a*a+b*b)))),1e-12))
    metrics=projection.metrics
    transmission=max(0.,policy.preferred_transmission_sine-metrics['minimum_primary_transmission_sine'])/max(policy.preferred_transmission_sine,1e-12)
    directionality=max(0.,policy.preferred_axis_variance_fraction-metrics['E_axis_variance_fraction'])/max(policy.preferred_axis_variance_fraction,1e-12)
    terms=dict(monotonicity=policy.monotonicity_weight*monotonicity,
        endpoint_zero=policy.endpoint_zero_weight*zero,turnaround=policy.turning_weight*turning/policy.turning_scale_rad,
        extra_crossing=policy.extra_crossing_weight*extra,speed_ratio=policy.speed_ratio_weight*ratio_error,
        displacement_fraction=policy.displacement_fraction_weight*fraction_error,
        long_mirror_asymmetry=policy.long_symmetry_weight*asymmetry,
        transmission_preference=policy.transmission_preference_weight*transmission,
        directionality_preference=policy.directionality_preference_weight*directionality)
    return dict(policy_version=PRIMARY_POLICY,score=float(sum(terms.values())),terms=terms,
        target=target_features,projection_axis=projection.axis.tolist(),
        monotonicity_wrong_sign_rms=monotonicity,endpoint_zero_rms=zero,turning_rms_rad=turning,
        matched_turnarounds_rad=dict(maximum=high_turn,minimum=low_turn),
        zero_crossing_count=len(high)+len(low),stationary_point_count=len(roots),extra_crossings=extra,
        fast_slow_speed_ratio=ratio,fast_displacement_fraction=fraction,
        speed_ratio_log_error=ratio_error,displacement_fraction_error=fraction_error,
        long_mirror_asymmetry_rms=asymmetry,candidate_long_span_rad=own_span,
        candidate_long_start_rad=own_start,candidate_long_end_rad=own_end,long_symmetry_enabled=policy.long_symmetry_weight>0,
        mechanical_preferences=dict(transmission=transmission,directionality=directionality))


def primary_evidence(artifact,target,side,samples=1440,*,cadence_policy=None,mechanical_policy=None):
    source=PeriodicTargetSide(target,side); angles=np.linspace(0.,PERIOD,samples,endpoint=False)
    provenance=artifact.data['provenance']
    saved=provenance.get('primary_cadence_policy',{})
    score_samples=provenance.get('primary_score_samples',samples)
    root_samples=provenance.get('primary_root_samples',samples)
    policy=cadence_policy or PrimaryCadencePolicy(**saved)
    features=primary_target_features(source,policy)
    mechanical_policy=mechanical_policy or PrimaryMechanicalPolicy(**provenance.get('primary_mechanical_policy',{}))
    primary=artifact.reconstruct()
    cycle=primary_cycle_data(primary,samples)
    primary_constraints=mechanical_policy.constraints(cycle[2])
    candidates=[]
    for sign in (1.,-1.):
        projection=PrimaryProjection(primary,sign=sign,root_samples=root_samples,cycle=cycle)
        try: cadence=primary_cadence_score(projection,features,policy,score_samples)
        except ValueError: continue
        candidates.append((cadence['score'],projection,cadence))
    if not candidates: raise ValueError('No admissible primary projection sign.')
    _,projection,cadence=min(candidates,key=lambda row:row[0])
    error=projection.value(angles)-1-source.position(angles)
    dv=None if not source.has_velocity else float(np.sqrt(np.mean((projection.value(angles,1)-source.velocity(angles))**2)))
    constraints=[]; deferred=[]
    for row in artifact.scientific['constraints']:
        if row['metric'] in projection.metrics:
            constraints.append(margin_record(row['metric'],projection.metrics[row['metric']],row['limit'],row['relation'],row['unit']))
        else: deferred.append(row)
    return dict(fit=None,primary_cadence=cadence,primary_projection=dict(position_rms=float(np.sqrt(np.mean(error**2))),
        position_maximum_error=float(np.max(abs(error))),velocity_rms_per_rad=dv,
        axis=projection.axis.tolist(),meaning='Diagnostic only: E is not P; position resemblance is excluded from primary search'),
        mechanical=dict(metrics=projection.metrics,constraints=constraints,primary_constraints=primary_constraints,
            primary_mechanical_screen=mechanical_policy.profile(),deferred_constraints=deferred,samples=samples),thermodynamic=None)


def relative_downstream_bounds(primary,policy):
    points=np.column_stack(primary.joint_state(np.linspace(0.,PERIOD,721))['joints']['E'])
    center=(points.min(axis=0)+points.max(axis=0))/2
    extent=max(float(np.linalg.norm(np.ptp(points,axis=0))),policy.minimum_trajectory_extent)
    radius=policy.pivot_envelope_radius*extent
    bounds=dict(second_pivot_x=(center[0]-radius,center[0]+radius),second_pivot_y=(center[1]-radius,center[1]+radius),
        link_ef=(policy.link_extent_minimum*extent,policy.link_extent_maximum*extent),
        link_gf=(policy.link_extent_minimum*extent,policy.link_extent_maximum*extent),
        h_along_over_ef=(-1.,2.),h_normal_over_ef=(-1.,1.),
        piston_rod=(policy.link_extent_minimum*extent,policy.rod_extent_maximum*extent),
        slider_axis_offset=(-float(np.linalg.norm(center))-radius,float(np.linalg.norm(center))+radius),
        slider_axis_angle=(-math.pi,math.pi))
    bounds.update({n:v for n,v in policy.bounds.items() if n in bounds})
    return bounds


class SixBarStageSearch(GeometrySearch):
    """Only stage-owned geometry callbacks; DE, budgets and polish are inherited."""
    def __init__(self,plan,policy,sides,volume_limits,parent=None):
        self.stage=plan.request.stages[0]; self.parent=parent
        self.parent_artifact=None if parent is None else MechanismArtifact.from_data(parent['mechanisms'][sides[0]])
        self.fixed={} if parent is None else dict(self.parent_artifact.scientific['geometry'])
        if parent is not None: plan=inherit_parent_plan(plan,parent,sides[0])
        super().__init__(plan,policy,sides,volume_limits)
        self.cadence_policy=PrimaryCadencePolicy(**policy.primary_cadence)
        self.primary_mechanical_policy=PrimaryMechanicalPolicy(**policy.primary_mechanical)
        self.target_features={side:primary_target_features(self.targets[side],self.cadence_policy) for side in sides} if self.stage=='primary_discovery' else {}

    def coordinate_bounds(self):
        full=self.policy.coordinates('six_bar')
        names=PRIMARY_COORDINATES if self.stage=='primary_discovery' else DOWNSTREAM_COORDINATES if self.stage=='downstream_fit' else SIXBAR_CONTINUOUS
        if self.stage not in ('primary_discovery','downstream_fit'):
            search=self.parent['metadata'].get('search',{})
            full.update(search.get('primary_bounds',{}))
            full.update(search.get('effective_bounds',{}))
            full.update(self.policy.bounds)
        if self.stage=='downstream_fit':
            primary=SixBarPrimaryMechanism(**{n:self.fixed[n] for n in (*PRIMARY_COORDINATES,'primary_branch')})
            full.update(relative_downstream_bounds(primary,self.policy))
        return {n:full[n] for n in names}

    def category_choices(self):
        choices=self.policy.choices('six_bar')
        if self.stage=='primary_discovery': return tuple(dict(primary_branch=b) for b in dict.fromkeys(c['primary_branch'] for c in choices))
        if self.stage=='downstream_fit':
            if self.fixed['primary_branch'] not in [c['primary_branch'] for c in choices]:
                raise ValueError('Parent primary branch is excluded by the search policy.')
            return tuple(dict(primary_branch=self.fixed['primary_branch'],second_branch=b)
                         for b in dict.fromkeys(c['second_branch'] for c in choices))
        categories={n:self.fixed[n] for n in CATEGORY_VALUES['six_bar']}
        if categories not in choices: raise ValueError('Parent branches are excluded by the search policy.')
        return (categories,)

    def decode(self,x,categories):
        p={} if self.stage=='primary_discovery' else dict(self.fixed)
        p.update(super().decode(x,categories))
        return p

    def evaluate(self,x,side,categories,origin,*,retain=True):
        origin=dict(origin)
        if self.parent is not None:
            if self.stage!='downstream_fit': origin['source_settings']=self.parent_artifact.scientific['settings']
            origin.update(parent_family_id=self.parent['family_id'],parent_artifact_hash=self.parent_artifact.content_hash,
                          primary_family_id=self.parent['metadata'].get('provenance',{}).get('primary_family_id',self.parent['family_id']))
        if self.stage!='primary_discovery':
            # The inherited closure -> topology -> constraints -> dense-fit chain.
            return super().evaluate(x,side,categories,origin,retain=retain)
        self.check_budget(); self.evaluations+=1
        p=self.decode(x,categories)
        row=dict(valid=False,side=side,geometry=p,categories=categories,settings=dict(family='six_bar',component='primary'),origin=origin)
        reason='primary_closure'; penalty=40.
        try:
            primary=SixBarPrimaryMechanism(**p)
            cycle=primary_cycle_data(primary,self.policy.mechanical_samples)
            records=self.primary_mechanical_policy.constraints(cycle[2])
            row.update(metrics=cycle[2],primary_constraints=records)
            failed=primary_rejection_reason(records)
            if failed:
                reason=failed; penalty=1000.
                self.record_constraint_violations('primary_mechanical_constraints',records)
                for record in records:
                    if not record['satisfied']:
                        named='primary_transmission_below_minimum' if record['name']=='minimum_primary_transmission_sine' else 'primary_e_span_below_minimum'
                        if named!=reason: self.rejections[named]+=1
                raise ValueError('Primary design screen fails.')
            reason='primary_projection_topology'; penalty=20.
            candidates=[PrimaryProjection(primary,sign=sign,root_samples=self.policy.root_samples,cycle=cycle) for sign in (1.,-1.)]
            for projection in candidates:
                reason='primary_constraints'; penalty=1000.
                constraints=[margin_record(r['metric'],projection.metrics[r['metric']],r['limit'],r['relation'],r['unit'])
                             for r in self.plan.request.mechanical_constraints if r['metric'] in projection.metrics]
                if any(not r['satisfied'] for r in constraints): raise ValueError('Primary constraint fails.')
                reason='primary_projection_topology'
                try: cadence=primary_cadence_score(projection,self.target_features[side],self.cadence_policy,self.policy.samples)
                except ValueError: continue
                score=cadence['score']
                if not row['valid'] or score<row['score']:
                    row.update(valid=True,score=score,residual=np.sqrt(list(cadence['terms'].values())),
                        metrics=projection.metrics,reason=None,terms=cadence['terms'],cadence=cadence,
                        projection_axis=projection.axis.tolist())
            if not row['valid']: raise ValueError('No admissible projection sign.')
            if retain: self.archive.add(row)
        except (ValueError,ArithmeticError,np.linalg.LinAlgError):
            row.update(score=penalty,reason=reason); self.rejections[reason]+=1
        return row

    def initial(self,rng,categories,index,side):
        if self.stage=='primary_discovery':
            p=dict(primary_ground=3.,primary_coupler=3.5,primary_rocker=2.5,
                   primary_e_along=2.5,primary_e_normal=.3,primary_phase=0.)
            x=np.clip(self.encode(p),0.,1.) if index==0 else rng.uniform(size=len(self.names))
            try:
                projection=PrimaryProjection(SixBarPrimaryMechanism(**self.decode(x,categories)))
                current=max(projection.roots,key=lambda r:r['position'])['angle_rad']
                desired=self.targets[side].extrema_angles()['maximum']
                p=self.decode(x,categories); p['primary_phase']+=current-desired
                x=self.encode(p); i=self.names.index('primary_phase')
                x[i]=x[i]%1 if math.isclose(self.span[i],PERIOD) else np.clip(x[i],0.,1.)
            except ValueError: pass
            return x
        if self.stage!='downstream_fit': return self.encode(self.fixed)
        # Pivot around E envelope; link sum/difference chosen against the full
        # sampled E-G distance range before constructing the second dyad.
        primary=SixBarPrimaryMechanism(**{n:self.fixed[n] for n in (*PRIMARY_COORDINATES,'primary_branch')})
        from dada_solver.six_bar import _rr, _rigid_point
        state=primary.joint_state(np.linspace(0.,PERIOD,721))
        points=np.column_stack(state['joints']['E']); center=points.mean(axis=0)
        extent=max(float(np.linalg.norm(np.ptp(points,axis=0))),self.policy.minimum_trajectory_extent)
        angle=0. if index==0 else rng.uniform(-math.pi,math.pi)
        pivot=center+extent*np.array([math.cos(angle),math.sin(angle)])
        distances=np.linalg.norm(points-pivot,axis=1)
        link=(float(distances.max())/2+extent*.25)*(1. if index==0 else rng.uniform(1.,1.5))
        p=dict(second_pivot_x=pivot[0],second_pivot_y=pivot[1],link_ef=link,link_gf=link,
               h_along_over_ef=.5 if index==0 else rng.uniform(0.,1.),h_normal_over_ef=0. if index==0 else rng.uniform(-.25,.25),
               piston_rod=3*extent,slider_axis_offset=0.,slider_axis_angle=0.)
        try:
            f,fd=_rr(state['joints']['E'],state['E_derivative'],pivot,link,link,categories['second_branch'])
            h,_=_rigid_point(state['joints']['E'],state['E_derivative'],f,fd,p['h_along_over_ef'],p['h_normal_over_ef'])
            h=np.column_stack(h); _,vectors=np.linalg.eigh(np.cov(h.T)); axis=vectors[:,-1]
            if axis[np.argmax(abs(axis))]<0: axis=-axis
            p['slider_axis_angle']=math.atan2(axis[1],axis[0]); normal=np.array([-axis[1],axis[0]])
            p['slider_axis_offset']=float(h.mean(axis=0)@normal)
            p['piston_rod']=max(2*extent,float(np.max(abs(h@normal-p['slider_axis_offset'])))+extent)
        except ValueError: pass
        # Both slider-axis directions are initialization alternatives, not new
        # scientific categories. The positive rod closure remains unchanged.
        from types import SimpleNamespace
        seed_angles=np.linspace(0.,PERIOD,361,endpoint=False)
        reference=self.targets[side].position(seed_angles)
        seeds=[]
        for reverse in (False,True):
            candidate=dict(p)
            if reverse:
                candidate['slider_axis_angle']=(candidate['slider_axis_angle']+2*math.pi)%(2*math.pi)-math.pi
                candidate['slider_axis_offset']=-candidate['slider_axis_offset']
            x=np.clip(self.encode(candidate),0.,1.)
            try:
                slider,_=SixBarCylinderMechanism._evaluate(SimpleNamespace(**self.decode(x,categories)),seed_angles)
                span=float(np.ptp(slider))
                if span<=1e-10: continue
                q=1-(slider-float(np.min(slider)))/span
                seeds.append((float(np.mean((q-reference)**2)),x))
            except ValueError: pass
        return min(seeds,key=lambda entry:entry[0])[1] if seeds else np.clip(self.encode(p),0.,1.)

    def final_library(self):
        from .synthesis import assess_mechanism
        members=[]
        for row in self.archive.rows:
            origin=row['origin']; stage=self.stage
            parent_id=origin.get('parent_family_id')
            family_id=(f"primary-{row['side'][0].upper()}-{len(members)+1:03d}" if stage=='primary_discovery' else
                       f"{parent_id}/{stage}-{len(members)+1:03d}")
            provenance=dict(target_hash=self.plan.target.content_hash,stage=stage,protocol=self.plan.request.protocol,
                            origin=origin,categories=row['categories'],primary_family_id=origin.get('primary_family_id',family_id))
            provenance['mechanical_screen']=screen_provenance(self.plan,self.parent,row['side'])
            if stage=='primary_discovery':
                provenance.update(search_policy=PRIMARY_POLICY,primary_cadence_policy=asdict(self.cadence_policy),
                                  primary_score_samples=self.policy.samples,primary_root_samples=self.policy.root_samples,
                                  primary_mechanical_policy=asdict(self.primary_mechanical_policy),
                                  primary_mechanical_screen=self.primary_mechanical_policy.profile())
            artifact=MechanismArtifact.create('six_bar',row['geometry'],settings=row['settings'],
                constraints=self.plan.request.mechanical_constraints,provenance=provenance)
            if stage=='primary_discovery':
                try:
                    evidence=primary_evidence(artifact,self.plan.target,row['side'],FINAL_MECHANICAL_SAMPLES,
                        cadence_policy=self.cadence_policy,mechanical_policy=self.primary_mechanical_policy)
                except (ValueError,ArithmeticError):
                    self.rejections['final_primary_closure_or_topology']+=1; continue
                records=evidence['mechanical']['primary_constraints']
                if self.record_constraint_violations('final_primary_mechanical_constraints',records):
                    for record in records:
                        if not record['satisfied']:
                            named='primary_transmission_below_minimum' if record['name']=='minimum_primary_transmission_sine' else 'primary_e_span_below_minimum'
                            self.rejections[named]+=1; self.rejections['final_'+named]+=1
                    continue
            else:
                try:
                    roots=artifact.reconstruct().stationary_points(samples=self.policy.root_samples)
                except (ValueError,ArithmeticError):
                    self.rejections['final_refined_closure']+=1; continue
                if len(roots)!=2 or {r['kind'] for r in roots}!={'maximum','minimum'}:
                    self.rejections['final_refined_topology']+=1; continue
                try:
                    evidence=assess_mechanism(artifact,self.plan.target,row['side'],mechanical_samples=FINAL_MECHANICAL_SAMPLES,
                                              volume_limits=self.volume_limits.get(row['side']))
                except (ValueError,ArithmeticError):
                    self.rejections['final_refined_closure']+=1; continue
                evidence['mechanical']['closure']=dict(satisfied=True,method='production reconstruction and 1440-sample geometry screen')
                evidence['mechanical']['topology']=dict(extrema=roots,method='analytic velocity roots; doubled-grid and tangency screen, not a certified proof')
            if self.record_constraint_violations('final_mechanical_constraints',evidence['mechanical']['constraints']):
                self.rejections['final_mechanical_constraints']+=1; continue
            if stage!='primary_discovery' and self.policy.maximum_position_rms is not None and evidence['fit']['position_rms']>self.policy.maximum_position_rms:
                self.rejections['final_maximum_position_rms']+=1; continue
            mechanisms={} if self.parent is None else dict(self.parent['mechanisms'])
            mechanisms[row['side']]=artifact.data
            previous={} if self.parent is None else dict(self.parent['metadata'].get('evidence',{}))
            previous[row['side']]=evidence
            rejected_parent=False
            for other,raw in mechanisms.items():
                if other==row['side'] or raw['scientific']['settings'].get('component'): continue
                try:
                    parent_artifact=MechanismArtifact.from_data(raw)
                    roots=parent_artifact.reconstruct().stationary_points(samples=self.policy.root_samples)
                    if len(roots)!=2 or {r['kind'] for r in roots}!={'maximum','minimum'}:
                        self.rejections['final_refined_topology']+=1; rejected_parent=True; continue
                    previous[other]=assess_mechanism(parent_artifact,self.plan.target,other,mechanical_samples=FINAL_MECHANICAL_SAMPLES,volume_limits=self.volume_limits.get(other))
                    previous[other]['mechanical']['topology']=dict(extrema=roots,method='analytic velocity roots; doubled-grid and tangency screen')
                    previous[other]['mechanical']['closure']=dict(satisfied=True,method='production reconstruction and 1440-sample geometry screen')
                    if self.record_constraint_violations('final_mechanical_constraints',previous[other]['mechanical']['constraints']):
                        self.rejections['final_mechanical_constraints']+=1; rejected_parent=True
                    if self.policy.maximum_position_rms is not None and previous[other]['fit']['position_rms']>self.policy.maximum_position_rms:
                        self.rejections['final_maximum_position_rms']+=1; rejected_parent=True
                except (ValueError,ArithmeticError):
                    self.rejections['final_refined_closure']+=1; rejected_parent=True
            if rejected_parent: continue
            metadata=dict(target_hash=self.plan.target.content_hash,evidence=previous,thermodynamic=None,provenance=provenance,
                search=dict(score=row['score'],terms=row['terms'],policy_version=PRIMARY_POLICY if stage=='primary_discovery' else STAGE_POLICY,policy=asdict(self.policy),
                    effective_bounds=self.bounds,primary_bounds=({n:self.bounds[n] for n in PRIMARY_COORDINATES} if stage=='primary_discovery' else
                        self.parent['metadata'].get('search',{}).get('primary_bounds',
                            {n:BOUNDS['six_bar'][n] for n in PRIMARY_COORDINATES})),
                    released_coordinates=self.names,stopping_reason=self.stop_reason,
                    evaluations=self.evaluations,cache_hits=self.cache_hits,rejections=dict(self.rejections),
                    islands=[{k:t[k] for k in ('side','categories','island','seed','generation','evaluations')} for t in self.tasks]))
            members.append(dict(family_id=family_id,mechanisms=mechanisms,metadata=metadata))
            if stage in ('full_local_polish','opposite_local_adaptation'): break
        if not members: raise ValueError('No admissible six-bar basin found under this stage, bounds and budget.')
        for member in members:
            member['metadata']['search'].update(rejections=dict(self.rejections),rejection_evidence=self.rejection_evidence,
                rejection_evidence_limit=64,mechanical_screen=member['metadata']['provenance']['mechanical_screen'])
        return MechanismLibrary(tuple(json.loads(json.dumps(members))))


def mirror_geometry(geometry):
    """Reflect local y and reverse crank time: q_mirror(theta)=q_source(-theta)."""
    p=dict(geometry)
    for name in ('primary_e_normal','primary_phase','second_pivot_y','h_normal_over_ef',
                 'slider_axis_offset','slider_axis_angle','primary_branch','second_branch'):
        p[name]=-p[name]
    p['primary_phase']%=PERIOD
    return p


def mirror_member(plan,parent,side,policy,volume_limits):
    from .synthesis import assess_mechanism
    other='small' if side=='large' else 'large'
    if other not in parent['mechanisms']: raise ValueError('Mirror needs a complete mechanism on the opposite side.')
    if side in parent['mechanisms']: raise ValueError('Mirror initialization will not overwrite an existing opposite mechanism.')
    artifact=MechanismArtifact.from_data(parent['mechanisms'][other]); raw=artifact.scientific
    if raw['settings'].get('family')!='six_bar' or raw['settings'].get('component'): raise ValueError('Mirror requires a complete six-bar mechanism.')
    plan=inherit_parent_plan(plan,parent,other)
    p=mirror_geometry(raw['geometry']); mechanism=SixBarCylinderMechanism(**p)
    desired=PeriodicTargetSide(plan.target,side).extrema_angles()['maximum']
    p['primary_phase']=(p['primary_phase']+mechanism.minimum_angle-desired)%PERIOD
    provenance=dict(stage='mirror_initialization',target_hash=plan.target.content_hash,
        parent_family_id=parent['family_id'],parent_artifact_hash=artifact.content_hash,
        source_side=other,destination_side=side,mirror_policy='initialization_only',
        angle_rule='q_seed(theta)=q_source(-theta-delta); delta aligns the destination maximum',
        primary_family_id=parent['metadata'].get('provenance',{}).get('primary_family_id',parent['family_id']),
        categories={n:p[n] for n in ('primary_branch','second_branch')})
    provenance['mechanical_screen']=screen_provenance(plan,parent,other)
    result=MechanismArtifact.create('six_bar',p,settings=raw['settings'],constraints=plan.request.mechanical_constraints,provenance=provenance)
    roots=result.reconstruct().stationary_points(samples=policy.root_samples)
    if len(roots)!=2 or {r['kind'] for r in roots}!={'maximum','minimum'}: raise ValueError('Mirror topology fails.')
    evidence=assess_mechanism(result,plan.target,side,volume_limits=(volume_limits or {}).get(side),mechanical_samples=FINAL_MECHANICAL_SAMPLES)
    if any(not c['satisfied'] for c in evidence['mechanical']['constraints']): raise ValueError('Mirror violates declared constraints.')
    evidence['mechanical']['topology']=dict(extrema=roots,method='reflection/time reversal followed by analytic-root screen')
    evidence['mechanical']['closure']=dict(satisfied=True,method='production reconstruction and 1440-sample geometry screen')
    if policy.maximum_position_rms is not None and evidence['fit']['position_rms']>policy.maximum_position_rms: raise ValueError('Mirror exceeds maximum_position_rms.')
    mechanisms=dict(parent['mechanisms']); mechanisms[side]=result.data
    previous=dict(parent['metadata']['evidence']); previous[side]=evidence
    source_roots=artifact.reconstruct().stationary_points(samples=policy.root_samples)
    if len(source_roots)!=2 or {r['kind'] for r in source_roots}!={'maximum','minimum'}: raise ValueError('Mirror source topology fails.')
    previous[other]=assess_mechanism(artifact,plan.target,other,mechanical_samples=FINAL_MECHANICAL_SAMPLES,volume_limits=(volume_limits or {}).get(other))
    previous[other]['mechanical']['topology']=dict(extrema=source_roots,method='analytic velocity roots; doubled-grid and tangency screen')
    previous[other]['mechanical']['closure']=dict(satisfied=True,method='production reconstruction and 1440-sample geometry screen')
    if any(not c['satisfied'] for c in previous[other]['mechanical']['constraints']): raise ValueError('Mirror source violates inherited constraints.')
    if policy.maximum_position_rms is not None and previous[other]['fit']['position_rms']>policy.maximum_position_rms: raise ValueError('Mirror source exceeds maximum_position_rms.')
    bounds=dict(BOUNDS['six_bar'])
    bounds.update(parent['metadata'].get('search',{}).get('primary_bounds',{}))
    bounds.update(parent['metadata'].get('search',{}).get('effective_bounds',{}))
    for name in ('primary_e_normal','second_pivot_y','h_normal_over_ef','slider_axis_offset','slider_axis_angle'):
        low,high=bounds[name];bounds[name]=(-high,-low)
    width=bounds['primary_phase'][1]-bounds['primary_phase'][0]
    bounds['primary_phase']=(0.,PERIOD) if math.isclose(width,PERIOD) else (p['primary_phase']-width/2,p['primary_phase']+width/2)
    return dict(family_id=parent['family_id']+'/mirror-'+side,mechanisms=mechanisms,
                metadata=dict(target_hash=plan.target.content_hash,evidence=previous,thermodynamic=None,provenance=provenance,
                              search=dict(policy_version=STAGE_POLICY,policy=asdict(policy),mechanical_screen=provenance['mechanical_screen'],effective_bounds=bounds,
                                          primary_bounds={n:bounds[n] for n in PRIMARY_COORDINATES},released_coordinates=[])))


def execute_six_bar(plan,*,policy=None,sides=('large',),library=None,volume_limits=None):
    if len(plan.request.stages)!=1 or plan.request.stages[0] not in STAGES:
        raise NotImplementedError('Execute one six_bar geometric stage; use mechanism adapt or mechanism retune for thermodynamic study generation.')
    if plan.request.protocol!='design_exploitation': raise NotImplementedError('Fresh-island saturation is a separate operator.')
    if not sides or len(set(sides))!=len(sides) or set(sides)-{'small','large'}: raise ValueError('Choose distinct cylinder sides.')
    policy=policy or SearchPolicy()
    if not isinstance(policy,SearchPolicy): raise ValueError('Synthesis requires a SearchPolicy.')
    policy.coordinates('six_bar');policy.choices('six_bar')
    stage=plan.request.stages[0]
    if stage=='primary_discovery':
        if library is not None or plan.request.retained_family_ids: raise ValueError('Primary discovery uses fresh islands, not complete source mechanisms.')
        engine=SixBarStageSearch(plan,policy,sides,volume_limits);engine.discover();return engine.final_library()
    if not isinstance(library,MechanismLibrary) or not plan.request.retained_family_ids:
        raise ValueError('This six-bar stage requires a library and explicit retained family IDs.')
    if stage in ('mirror_initialization','opposite_local_adaptation') and len(sides)!=1:
        raise ValueError('Choose exactly one destination side for opposite-piston stages.')
    results=[]; failed=[]; remaining=policy.max_evaluations
    import time
    deadline=None if policy.budget_seconds is None else time.monotonic()+policy.budget_seconds
    pending=len(plan.request.retained_family_ids)*len(sides)
    if remaining is not None and remaining<pending: raise ValueError('Evaluation budget must cover every selected parent and side.')
    for identifier in plan.request.retained_family_ids:
        parent=library.member(identifier)
        if parent['metadata'].get('target_hash')!=plan.target.content_hash: raise ValueError('Library and target identities disagree.')
        for side in sides:
            if stage=='mirror_initialization':
                results.append(mirror_member(plan,parent,side,policy,volume_limits));continue
            if side not in parent['mechanisms']: raise ValueError('Selected family has no mechanism on this side.')
            artifact=MechanismArtifact.from_data(parent['mechanisms'][side]);raw=artifact.scientific
            if raw['settings']['family']!='six_bar': raise ValueError('Selected artifact is not six_bar.')
            primary=raw['settings'].get('component')=='primary'
            if primary!=(stage=='downstream_fit'): raise ValueError('Downstream needs a primary; local stages need a complete six-bar.')
            if stage=='opposite_local_adaptation':
                other='small' if side=='large' else 'large'
                if other not in parent['mechanisms'] or artifact.data['provenance'].get('stage')!='mirror_initialization':
                    raise ValueError('Opposite adaptation requires an independent pair with a mirrored seed on the selected side.')
            if remaining is not None and remaining<=0: break
            seconds=None if deadline is None else deadline-time.monotonic()
            if seconds is not None and seconds<=0: break
            effective=replace(policy,max_evaluations=None if remaining is None else max(1,remaining//pending),
                              budget_seconds=None if seconds is None else seconds/pending)
            pending-=1
            engine=SixBarStageSearch(plan,effective,(side,),volume_limits,parent)
            if stage=='downstream_fit': engine.discover()
            else:
                x=engine.encode(raw['geometry'])
                if np.any(x<0) or np.any(x>1): raise ValueError('Parent lies outside configured bounds.')
                try:
                    row=engine.evaluate(x,side,engine.choices[0],{})
                    if not row['valid']: raise ValueError('Parent fails current constraints/topology.')
                    engine.fallback_parent=row
                    engine.polish(row,policy.polish_evaluations)
                except SearchStopped: pass
                fallback=[engine.fallback_parent] if hasattr(engine,'fallback_parent') else []
                engine.archive.rows=sorted([*engine.archive.rows,*fallback],key=lambda r:r['score'])
            if remaining is not None: remaining-=engine.evaluations
            if not engine.archive.rows:
                failed.append(dict(parent_family_id=identifier,side=side,evaluations=engine.evaluations,
                                   rejections=dict(engine.rejections),rejection_evidence=engine.rejection_evidence,reason='no admissible coarse basin'))
                continue
            try:
                results.extend(engine.final_library().members)
            except ValueError as error:
                if str(error)!='No admissible six-bar basin found under this stage, bounds and budget.': raise
                failed.append(dict(parent_family_id=identifier,side=side,evaluations=engine.evaluations,
                                   rejections=dict(engine.rejections),rejection_evidence=engine.rejection_evidence,reason='no admissible refined basin'))
    if not results:
        summary=[{k:v for k,v in row.items() if k!='rejection_evidence'} for row in failed]
        raise ValueError(f'No admissible six-bar descendants from the selected parents: {summary}.')
    for member in results:
        member['metadata']['search']['failed_parent_searches']=failed
    return MechanismLibrary(tuple(json.loads(json.dumps(results))))
