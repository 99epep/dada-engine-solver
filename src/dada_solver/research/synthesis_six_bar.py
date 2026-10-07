"""Six-bar stage adapters for the shared GeometrySearch island/polish engine.

No differential-evolution implementation or thermodynamic evaluator lives here.
Intermediate primary evidence is a projection opportunity, never a piston fit.
"""
from dataclasses import asdict, replace
import math
import json
import numpy as np

from dada_solver.six_bar import SixBarPrimaryMechanism, SixBarCylinderMechanism, velocity_stationary_points
from .artifacts import MechanismArtifact, MechanismLibrary
from .families import PRIMARY_COORDINATES, DOWNSTREAM_COORDINATES, SIXBAR_CONTINUOUS
from .margins import margin_record
from .motion_target import PERIOD, PeriodicTargetSide
from .synthesis_search import GeometrySearch, SearchPolicy, SearchStopped, BOUNDS, CATEGORY_VALUES

STAGE_POLICY = 'hierarchical_six_bar_geometry_v1'
STAGES = ('primary_discovery','downstream_fit','full_local_polish','mirror_initialization','opposite_local_adaptation')


class PrimaryProjection:
    """PCA projection of E for inspectable chronology, not a cylinder law."""
    def __init__(self, primary, *, axis=None, sign=1.):
        self.geometry=primary
        angles=np.linspace(0.,PERIOD,721,endpoint=False)
        state=primary.joint_state(angles)
        points=np.column_stack(state['joints']['E'])
        if axis is None:
            _,vectors=np.linalg.eigh(np.cov(points.T))
            axis=vectors[:,-1]
            if axis[np.argmax(abs(axis))]<0: axis=-axis
        self.axis=np.asarray(axis)*sign
        roots=velocity_stationary_points(self.coordinate,self.derivative,360)
        if len(roots)<2:
            raise ValueError('Projected E has no useful chronology.')
        values=[r['position'] for r in roots]
        self.low,self.span=min(values),max(values)-min(values)
        if self.span<=1e-10: raise ValueError('Degenerate primary E projection.')
        self.roots=roots
        self.lateral_rms=float(np.std(points@np.array([-self.axis[1],self.axis[0]])))/self.span
        self.metrics=dict(minimum_primary_transmission_sine=float(np.min(state['primary_transmission_sine'])),
                          E_projection_span_over_crank=self.span,E_lateral_rms_over_projection_span=self.lateral_rms,
                          E_projection_extrema_count=len(roots))

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


def primary_evidence(artifact,target,side,samples=1440):
    source=PeriodicTargetSide(target,side); angles=np.linspace(0.,PERIOD,samples,endpoint=False)
    projections=[PrimaryProjection(artifact.reconstruct(),sign=sign) for sign in (1.,-1.)]
    projection=min(projections,key=lambda p:float(np.mean((p.value(angles)-1-source.position(angles))**2)))
    error=projection.value(angles)-1-source.position(angles)
    dv=None if not source.has_velocity else float(np.sqrt(np.mean((projection.value(angles,1)-source.velocity(angles))**2)))
    constraints=[]; deferred=[]
    for row in artifact.scientific['constraints']:
        if row['metric'] in projection.metrics:
            constraints.append(margin_record(row['metric'],projection.metrics[row['metric']],row['limit'],row['relation'],row['unit']))
        else: deferred.append(row)
    return dict(fit=None,primary_projection=dict(position_rms=float(np.sqrt(np.mean(error**2))),
        position_maximum_error=float(np.max(abs(error))),velocity_rms_per_rad=dv,
        axis=projection.axis.tolist(),meaning='E projection chronology opportunity; E is not P'),
        mechanical=dict(metrics=projection.metrics,constraints=constraints,deferred_constraints=deferred,samples=samples),thermodynamic=None)


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
        super().__init__(plan,policy,sides,volume_limits)

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
            reason='primary_projection_topology'; penalty=20.
            candidates=[PrimaryProjection(primary,sign=sign) for sign in (1.,-1.)]
            angles,weights=self.mesh[side]; target,velocity=self.references[side]
            for projection in candidates:
                reason='primary_constraints'; penalty=10.
                constraints=[margin_record(r['metric'],projection.metrics[r['metric']],r['limit'],r['relation'],r['unit'])
                             for r in self.plan.request.mechanical_constraints if r['metric'] in projection.metrics]
                if any(not r['satisfied'] for r in constraints): raise ValueError('Primary constraint fails.')
                error=projection.value(angles)-1-target; position=float(weights@(error**2))
                vmse=None; vt=0.; residual=np.sqrt(weights)*error
                if velocity is not None:
                    scale=max(float(np.sqrt(weights@(velocity**2))),1/PERIOD)
                    dv=(projection.value(angles,1)-velocity)/scale; vmse=float(weights@(dv**2))
                    vt=self.policy.velocity_weight*vmse/(1+vmse)
                    residual=np.r_[residual,math.sqrt(self.policy.velocity_weight/(1+vmse))*np.sqrt(weights)*dv]
                # Projection is a guide, not an equality constraint E=P. Lateral
                # potential is exposed; no universal near-linearity floor exists.
                lateral_term=self.policy.primary_linearity_weight*projection.lateral_rms**2/(1+projection.lateral_rms**2)
                topology_term=self.policy.primary_topology_weight*max(0.,len(projection.roots)-2)
                residual=np.r_[residual,math.sqrt(lateral_term),math.sqrt(topology_term)]
                score=position+vt+lateral_term+topology_term
                if not row['valid'] or score<row['score']:
                    row.update(valid=True,score=score,residual=residual,metrics=projection.metrics,reason=None,
                        terms=dict(projected_chronology_mse=position,normalized_velocity_mse=vmse,velocity_term=vt,
                                   projection_linearity_preference=lateral_term,projection_topology_preference=topology_term),
                        projection_axis=projection.axis.tolist())
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
            artifact=MechanismArtifact.create('six_bar',row['geometry'],settings=row['settings'],
                constraints=self.plan.request.mechanical_constraints,provenance=provenance)
            if stage=='primary_discovery':
                evidence=primary_evidence(artifact,self.plan.target,row['side'],self.policy.root_samples)
            else:
                try:
                    roots=artifact.reconstruct().stationary_points(samples=self.policy.root_samples)
                except (ValueError,ArithmeticError):
                    self.rejections['final_refined_closure']+=1; continue
                if len(roots)!=2 or {r['kind'] for r in roots}!={'maximum','minimum'}:
                    self.rejections['final_refined_topology']+=1; continue
                evidence=assess_mechanism(artifact,self.plan.target,row['side'],mechanical_samples=self.policy.mechanical_samples,
                                          volume_limits=self.volume_limits.get(row['side']))
                evidence['mechanical']['topology']=dict(extrema=roots,method='analytic velocity roots; doubled-grid and tangency screen, not a certified proof')
            if any(not r['satisfied'] for r in evidence['mechanical']['constraints']): continue
            mechanisms={} if self.parent is None else dict(self.parent['mechanisms'])
            mechanisms[row['side']]=artifact.data
            previous={} if self.parent is None else dict(self.parent['metadata'].get('evidence',{}))
            previous[row['side']]=evidence
            metadata=dict(target_hash=self.plan.target.content_hash,evidence=previous,thermodynamic=None,provenance=provenance,
                search=dict(score=row['score'],terms=row['terms'],policy_version=STAGE_POLICY,policy=asdict(self.policy),
                    effective_bounds=self.bounds,primary_bounds=({n:self.bounds[n] for n in PRIMARY_COORDINATES} if stage=='primary_discovery' else
                        self.parent['metadata'].get('search',{}).get('primary_bounds',
                            {n:BOUNDS['six_bar'][n] for n in PRIMARY_COORDINATES})),
                    released_coordinates=self.names,stopping_reason=self.stop_reason,
                    evaluations=self.evaluations,cache_hits=self.cache_hits,rejections=dict(self.rejections),
                    islands=[{k:t[k] for k in ('side','categories','island','seed','generation','evaluations')} for t in self.tasks]))
            members.append(dict(family_id=family_id,mechanisms=mechanisms,metadata=metadata))
            if stage in ('full_local_polish','opposite_local_adaptation'): break
        if not members: raise ValueError('No admissible six-bar basin found under this stage, bounds and budget.')
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
    if any(c not in plan.request.mechanical_constraints for c in raw['constraints']):
        raise ValueError('Retain all source artifact constraints.')
    p=mirror_geometry(raw['geometry']); mechanism=SixBarCylinderMechanism(**p)
    desired=PeriodicTargetSide(plan.target,side).extrema_angles()['maximum']
    p['primary_phase']=(p['primary_phase']+mechanism.minimum_angle-desired)%PERIOD
    provenance=dict(stage='mirror_initialization',target_hash=plan.target.content_hash,
        parent_family_id=parent['family_id'],parent_artifact_hash=artifact.content_hash,
        source_side=other,destination_side=side,mirror_policy='initialization_only',
        angle_rule='q_seed(theta)=q_source(-theta-delta); delta aligns the destination maximum',
        primary_family_id=parent['metadata'].get('provenance',{}).get('primary_family_id',parent['family_id']),
        categories={n:p[n] for n in ('primary_branch','second_branch')})
    result=MechanismArtifact.create('six_bar',p,settings=raw['settings'],constraints=plan.request.mechanical_constraints,provenance=provenance)
    roots=result.reconstruct().stationary_points(samples=policy.root_samples)
    if len(roots)!=2 or {r['kind'] for r in roots}!={'maximum','minimum'}: raise ValueError('Mirror topology fails.')
    evidence=assess_mechanism(result,plan.target,side,volume_limits=(volume_limits or {}).get(side),mechanical_samples=policy.mechanical_samples)
    if any(not c['satisfied'] for c in evidence['mechanical']['constraints']): raise ValueError('Mirror violates declared constraints.')
    evidence['mechanical']['topology']=dict(extrema=roots,method='reflection/time reversal followed by analytic-root screen')
    mechanisms=dict(parent['mechanisms']); mechanisms[side]=result.data
    previous=dict(parent['metadata']['evidence']); previous[side]=evidence
    bounds=dict(BOUNDS['six_bar'])
    bounds.update(parent['metadata'].get('search',{}).get('primary_bounds',{}))
    bounds.update(parent['metadata'].get('search',{}).get('effective_bounds',{}))
    for name in ('primary_e_normal','second_pivot_y','h_normal_over_ef','slider_axis_offset','slider_axis_angle'):
        low,high=bounds[name];bounds[name]=(-high,-low)
    width=bounds['primary_phase'][1]-bounds['primary_phase'][0]
    bounds['primary_phase']=(0.,PERIOD) if math.isclose(width,PERIOD) else (p['primary_phase']-width/2,p['primary_phase']+width/2)
    return dict(family_id=parent['family_id']+'/mirror-'+side,mechanisms=mechanisms,
                metadata=dict(target_hash=plan.target.content_hash,evidence=previous,thermodynamic=None,provenance=provenance,
                              search=dict(policy_version=STAGE_POLICY,policy=asdict(policy),effective_bounds=bounds,
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
            if any(c not in plan.request.mechanical_constraints for c in raw['constraints']): raise ValueError('Retain all source artifact constraints.')
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
                                   rejections=dict(engine.rejections),reason='no admissible coarse basin'))
                continue
            try:
                results.extend(engine.final_library().members)
            except ValueError as error:
                if str(error)!='No admissible six-bar basin found under this stage, bounds and budget.': raise
                failed.append(dict(parent_family_id=identifier,side=side,evaluations=engine.evaluations,
                                   rejections=dict(engine.rejections),reason='no admissible refined basin'))
    if not results: raise ValueError(f'No admissible six-bar descendants from the selected parents: {failed}.')
    for member in results:
        member['metadata']['search']['failed_parent_searches']=failed
    return MechanismLibrary(tuple(json.loads(json.dumps(results))))
