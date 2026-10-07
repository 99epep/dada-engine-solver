"""Deterministic multi-island geometry discovery and category-preserving polish.

Search scores organize a diverse archive; they are never thermodynamic evidence.
The engine does not import, construct or call a thermodynamic evaluator.
"""
from collections import Counter, OrderedDict
from dataclasses import dataclass, field, asdict
import itertools
import json
import math
import time

import numpy as np
from scipy.optimize import least_squares

from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.four_bar import FourBarLoop, FourBarSliderAssembly, SliderConstraint, RockerOutputPoint, CouplerOutputPoint
from dada_solver.mechanism_diagnostics import zero_crossings
from .families import parameter_specs, build_side, side_metrics
from .margins import margin_record
from .motion_target import PERIOD, PeriodicTargetSide
from .artifacts import MechanismLibrary, MechanismArtifact

POLICY_VERSION = 'direct_geometry_islands_v1'
BOUNDS_VERSION = 'crank_normalized_direct_bounds_v1'
BOUNDS = {
    'slider_crank': dict(rod_over_crank=(1.2,12.),offset_over_crank=(-4.,4.),phase_rad=(0.,PERIOD)),
    'six_bar': dict(primary_ground=(1.2,10.),primary_coupler=(1.2,12.),primary_rocker=(1.2,10.),
        primary_e_along=(-8.,12.),primary_e_normal=(-8.,8.),primary_phase=(0.,PERIOD),
        second_pivot_x=(-20.,20.),second_pivot_y=(-20.,20.),link_ef=(.1,30.),link_gf=(.1,30.),
        h_along_over_ef=(-1.,2.),h_normal_over_ef=(-1.,1.),piston_rod=(.1,60.),
        slider_axis_offset=(-30.,30.),slider_axis_angle=(-math.pi,math.pi)),
    'four_bar': dict(coupler=(1.2,8.),rocker=(1.2,8.),ground_x=(-6.,6.),ground_y=(-6.,6.),
        output_along=(-6.,6.),output_normal=(-4.,4.),rod_length=(1.5,16.),
        slider_origin_x=(-8.,8.),slider_origin_y=(-8.,8.),axis_angle=(-math.pi,math.pi),phase_rad=(0.,PERIOD)),
}
CATEGORY_VALUES = {
    'slider_crank': dict(volume_increases_with_coordinate=(True,False)),
    'six_bar': dict(primary_branch=(-1,1),second_branch=(-1,1)),
    'four_bar': dict(output=('rocker','coupler'),loop_branch=(-1,1),slider_branch=(-1,1),
                     crank_direction=(-1,1),volume_increases_with_coordinate=(True,False)),
}
PERIODIC_COORDINATES = {'phase_rad','axis_angle','primary_phase','slider_axis_angle'}


@dataclass(frozen=True)
class SearchPolicy:
    """Versioned numerical search policy; bounds are not physical domains."""
    seed: int = 1234
    islands: int = 4
    population: int = 16
    generations: int = 24
    max_evaluations: int | None = None
    budget_seconds: float | None = None
    samples: int = 721
    mechanical_samples: int = 720
    root_samples: int = 1440
    cluster_distance: float = .05
    retain_per_side: int = 64
    velocity_weight: float = .001
    feature_weight: float = .1
    event_radius_fraction: float = 1/36
    local_radius: float = .1
    polish_evaluations: int = 80
    discovery_polish_evaluations: int = 24
    bounds: dict = field(default_factory=dict)
    categories: dict = field(default_factory=dict)
    primary_linearity_weight: float = .001
    primary_topology_weight: float = .001
    minimum_trajectory_extent: float = .1
    pivot_envelope_radius: float = 2.
    link_extent_minimum: float = .25
    link_extent_maximum: float = 4.
    rod_extent_maximum: float = 8.

    def validate_relative_bounds(self):
        for name in ('minimum_trajectory_extent','pivot_envelope_radius','link_extent_minimum','link_extent_maximum','rod_extent_maximum'):
            value=getattr(self,name)
            if isinstance(value,bool) or not math.isfinite(value) or value<=0:
                raise ValueError(f'{name} must be positive and finite.')
        if self.link_extent_minimum >= self.link_extent_maximum:
            raise ValueError('Relative link bounds must be increasing.')

    def __post_init__(self):
        self.validate_relative_bounds()
        for name in ('primary_linearity_weight','primary_topology_weight'):
            value=getattr(self,name)
            if isinstance(value,bool) or not math.isfinite(value) or not 0<=value<=.1:
                raise ValueError(f'{name} must lie in [0, 0.1].')
        for name,minimum in (('seed',0),('islands',1),('population',4),('generations',0),
                ('samples',360),('mechanical_samples',360),('root_samples',360),('retain_per_side',2),
                ('polish_evaluations',1),('discovery_polish_evaluations',0)):
            if type(getattr(self,name)) is not int or getattr(self,name)<minimum:
                raise ValueError(f'{name} must be an integer >= {minimum}.')
        if self.max_evaluations is not None and (type(self.max_evaluations) is not int or self.max_evaluations<1):
            raise ValueError('max_evaluations must be a positive integer.')
        if self.budget_seconds is not None and (isinstance(self.budget_seconds,bool) or not math.isfinite(self.budget_seconds) or self.budget_seconds<=0):
            raise ValueError('Geometry budget must be finite and positive.')
        for name in ('cluster_distance','local_radius','event_radius_fraction'):
            value=getattr(self,name)
            if isinstance(value,bool) or not math.isfinite(value) or not 0<value<1:
                raise ValueError(f'{name} must be finite and lie in (0, 1).')
        for name in ('velocity_weight','feature_weight'):
            value=getattr(self,name)
            if isinstance(value,bool) or not math.isfinite(value) or not 0<=value<=.1:
                raise ValueError(f'{name} must lie in [0, 0.1].')

    def coordinates(self,family):
        if family not in BOUNDS: raise NotImplementedError('Unsupported geometric synthesis family.')
        bounds=dict(BOUNDS[family])
        if set(self.bounds)-set(bounds): raise ValueError('Unknown synthesis search-bound coordinate.')
        bounds.update(self.bounds)
        specs=parameter_specs(dict(family=family,**({'output':'rocker'} if family=='four_bar' else {})),'small')
        for name,interval in bounds.items():
            if (len(interval)!=2 or any(isinstance(v,bool) or not isinstance(v,(float,int)) or not math.isfinite(v) for v in interval)
                    or interval[0]>=interval[1]): raise ValueError(f'Invalid search bounds for {name}.')
            for value in interval: specs[name].validate(value)
            if name in PERIODIC_COORDINATES and interval[1]-interval[0]>PERIOD+1e-12:
                raise ValueError('Angular search bounds cannot exceed one full period.')
        return bounds

    def choices(self,family):
        choices=dict(CATEGORY_VALUES[family])
        if set(self.categories)-set(choices): raise ValueError('Unknown synthesis category.')
        for name,values in self.categories.items():
            if not isinstance(values,(tuple,list)) or not values or len(values)!=len(set(values)):
                raise ValueError('Categories require distinct explicit choices.')
            for value in values:
                if not any(type(value) is type(v) and value==v for v in choices[name]):
                    raise ValueError(f'Invalid scientific category: {name}.')
            choices[name]=tuple(values)
        return tuple(dict(zip(choices,values)) for values in itertools.product(*choices.values()))


def geometric_descriptor(family,geometry,bounds):
    """Crank-unit geometry with rigid drawing rotation/axis-origin gauge removed.

    Categories are compared separately. Angular distance is circular, never the
    raw difference between coordinates at opposite sides of a periodic seam.
    """
    p=dict(geometry)
    if family=='four_bar':
        c,s=math.cos(p['axis_angle']),math.sin(p['axis_angle'])
        for x,y in (('ground_x','ground_y'),('slider_origin_x','slider_origin_y')):
            p[x],p[y]=c*p[x]+s*p[y],-s*p[x]+c*p[y]
        p['phase_rad']=(p['phase_rad']-p['axis_angle'])%PERIOD
        names=[n for n in bounds if n not in ('axis_angle','slider_origin_x')]
    else: names=list(bounds)
    return tuple((p[name],PERIOD if name in PERIODIC_COORDINATES else bounds[name][1]-bounds[name][0]) for name in names)


def geometric_distance(family,left,right,bounds):
    a,b=geometric_descriptor(family,left,bounds),geometric_descriptor(family,right,bounds)
    differences=[]
    names=[n for n in bounds if family!='four_bar' or n not in ('axis_angle','slider_origin_x')]
    for name,(x,scale),(y,_) in zip(names,a,b):
        delta=(x-y+math.pi)%PERIOD-math.pi if name in PERIODIC_COORDINATES else x-y
        differences.append((delta/scale)**2)
    if family=='six_bar':
        from .families import PRIMARY_COORDINATES
        groups=[[d for name,d in zip(names,differences) if (name in PRIMARY_COORDINATES)==primary]
                for primary in (True,False)]
        groups=[g for g in groups if g]
        return math.sqrt(sum(sum(g)/len(g) for g in groups)/len(groups))
    return math.sqrt(sum(differences)/len(differences))


class BasinArchive:
    """Greedy score-ordered geometric clustering, category-stratified retention."""
    def __init__(self,family,bounds,policy):
        self.family,self.bounds,self.policy=family,bounds,policy
        self.rows=[]

    def add(self,row):
        if not row['valid']: return
        rows=sorted([*self.rows,row],key=lambda r:(r['score'],tuple(r['geometry'].values())))
        retained=[]
        for candidate in rows:
            if any(candidate['side']==old['side'] and candidate['categories']==old['categories'] and
                   candidate['origin'].get('parent_family_id')==old['origin'].get('parent_family_id') and
                   geometric_distance(self.family,candidate['geometry'],old['geometry'],self.bounds)<self.policy.cluster_distance
                   for old in retained): continue
            retained.append(candidate)
        # Round-robin categories: no one category can occupy the entire archive.
        grouped={}
        for candidate in retained:
            key=(candidate['side'],tuple(candidate['categories'].items()),candidate['origin'].get('parent_family_id'))
            grouped.setdefault(key,[]).append(candidate)
        self.rows=[]
        for side in ('small','large'):
            groups=[rows for (s,_,_),rows in grouped.items() if s==side]
            groups.sort(key=lambda rows:tuple(rows[0]['categories'].items()))
            for level in range(max([len(rows) for rows in groups],default=0)):
                for group in groups:
                    if level<len(group) and sum(r['side']==side for r in self.rows)<self.policy.retain_per_side:
                        self.rows.append(group[level])


class SearchStopped(Exception):
    pass


class GeometrySearch:
    """One shared evaluation engine, independent island populations and archive."""
    def __init__(self,plan,policy,sides,volume_limits):
        self.plan,self.policy=plan,policy
        self.family=plan.request.mechanism_family
        self.bounds=self.coordinate_bounds()
        self.names=list(self.bounds)
        self.lower=np.array([self.bounds[n][0] for n in self.names])
        self.span=np.array([self.bounds[n][1]-self.bounds[n][0] for n in self.names])
        self.choices=self.category_choices()
        if policy.retain_per_side<len(self.choices): raise ValueError('Retention must allow at least one member per explored category.')
        self.sides=tuple(sides)
        self.volume_limits=volume_limits or {}
        self.targets={side:PeriodicTargetSide(plan.target,side) for side in sides}
        self.mesh={}; self.references={}
        for side,source in self.targets.items():
            # A compatible target has two direction changes, not an oscillatory topology.
            angles=np.linspace(0.,PERIOD,policy.samples,endpoint=False)
            derivative=source.velocity(angles) if source.has_velocity else source.interpolant(angles,1)
            if zero_crossings(derivative)!=2: raise ValueError('Synthesis target must have one maximum and one minimum.')
            blocks=[]
            for event in source.events:
                if event['kind'] not in ('maximum','minimum','turnaround','cadence_change','rounding','kink'): continue
                angle=event['angle_rad']; width=PERIOD*policy.event_radius_fraction
                if event['end_angle_rad'] is not None:
                    points=angle+np.linspace(0.,(event['end_angle_rad']-angle)%PERIOD,9)
                else: points=angle+np.linspace(-width,width,9)
                blocks.append(points%PERIOD)
            feature=np.concatenate(blocks) if blocks else np.array([])
            grid=np.r_[angles,feature]
            weights=np.r_[np.full(len(angles),1/len(angles)),np.full(len(feature),policy.feature_weight/len(feature)) if len(feature) else []]
            self.mesh[side]=(grid,weights/weights.sum())
            self.references[side]=(source.position(grid),source.velocity(grid))
            if side not in self.volume_limits and any(r['metric'].startswith('maximum_absolute_') for r in plan.request.mechanical_constraints):
                raise ValueError('Dimensional derivative constraints require explicit cylinder volume limits.')
        self.archive=BasinArchive(self.family,self.bounds,policy)
        self.evaluations=0; self.cache_hits=0; self.rejections=Counter(); self.cache=OrderedDict()
        self.deadline=None if policy.budget_seconds is None else time.monotonic()+policy.budget_seconds
        self.tasks=[]; self.stop_reason='completed_generations'

    def coordinate_bounds(self):
        return self.policy.coordinates(self.family)

    def category_choices(self):
        return self.policy.choices(self.family)

    def check_budget(self):
        if self.policy.max_evaluations is not None and self.evaluations>=self.policy.max_evaluations:
            self.stop_reason='evaluation_limit'; raise SearchStopped
        if self.deadline is not None and time.monotonic()>=self.deadline:
            self.stop_reason='cooperative_wall_deadline'; raise SearchStopped

    def decode(self,x,categories):
        values=self.lower+self.span*np.asarray(x)
        p={name:float(v) for name,v in zip(self.names,values)}
        for name in PERIODIC_COORDINATES & set(p):
            low,high=self.bounds[name]
            if math.isclose(high-low,PERIOD,abs_tol=1e-12): p[name]=low+(p[name]-low)%PERIOD
        p.update({n:v for n,v in categories.items() if n!='output'})
        return p

    def encode(self,geometry): return (np.array([geometry[n] for n in self.names])-self.lower)/self.span

    def evaluate(self,x,side,categories,origin,*,retain=True):
        self.check_budget()
        geometry=self.decode(x,categories)
        key=(side,tuple(categories.items()),tuple(geometry.values()),tuple(sorted(origin.get('source_settings',{}).items())))
        if key in self.cache:
            self.cache_hits+=1
            row=dict(self.cache[key],origin=origin)
            if retain and row['valid']: self.archive.add(row)
            return row
        self.evaluations+=1
        settings=origin.get('source_settings') or dict(family=self.family,**({'output':categories['output']} if self.family=='four_bar' else {}))
        row=dict(valid=False,side=side,geometry=geometry,categories=categories,settings=settings,origin=origin)
        stage='loop_closure'; penalty=40.
        try:
            if self.family=='four_bar':
                p=geometry; ground=math.hypot(p['ground_x'],p['ground_y'])
                margin=min(p['coupler']+p['rocker']-ground-1,abs(ground-1)-abs(p['coupler']-p['rocker']),abs(ground-1))
                if margin<=1e-10*max(p['coupler'],p['rocker'],ground+1):
                    penalty+=max(0.,-margin)/(1+ground+p['coupler']+p['rocker'])
                    raise ValueError('Full-revolution loop closure fails.')
                assembly=FourBarSliderAssembly(FourBarLoop(p['coupler'],p['rocker'],p['ground_x'],p['ground_y'],p['loop_branch']),
                    (RockerOutputPoint if settings['output']=='rocker' else CouplerOutputPoint)(p['output_along'],p['output_normal']),
                    SliderConstraint(p['slider_origin_x'],p['slider_origin_y'],p['axis_angle'],p['rod_length'],p['slider_branch']))
                angles=np.linspace(0.,PERIOD,361,endpoint=False)
                crank=np.column_stack((np.cos(angles),np.sin(angles)))
                derivative=np.column_stack((-np.sin(angles),np.cos(angles)))
                primary=assembly.evaluate_many(crank,derivative,primary_only=True)
                point=np.column_stack(primary['output_point'])
                normal=np.array([-math.sin(p['axis_angle']),math.cos(p['axis_angle'])])
                transverse=(point-np.array([p['slider_origin_x'],p['slider_origin_y']]))@normal
                stage='rod_closure'; penalty=30.+max(0.,float(np.max(abs(transverse)))/p['rod_length']-1)
                if np.max(abs(transverse))>=p['rod_length']: raise ValueError('Piston rod cannot close.')
            elif self.family=='six_bar':
                from dada_solver.six_bar import SixBarPrimaryMechanism, SixBarCylinderMechanism, _rr
                from .families import PRIMARY_COORDINATES
                from types import SimpleNamespace
                primary=SixBarPrimaryMechanism(**{n:geometry[n] for n in (*PRIMARY_COORDINATES,'primary_branch')})
                angles=np.linspace(0.,PERIOD,361,endpoint=False)
                state=primary.joint_state(angles)
                stage='secondary_closure'; penalty=35.
                _rr(state['joints']['E'],state['E_derivative'],
                    (geometry['second_pivot_x'],geometry['second_pivot_y']),geometry['link_ef'],geometry['link_gf'],geometry['second_branch'])
                stage='rod_closure'; penalty=30.
                SixBarCylinderMechanism._evaluate(SimpleNamespace(**geometry),angles)
            elif self.family=='slider_crank':
                gap=geometry['rod_over_crank']-1-abs(geometry['offset_over_crank'])
                if gap<=0:
                    penalty+=abs(gap)/geometry['rod_over_crank']; raise ValueError('Slider-crank cannot close.')
            stage='stroke_or_singularity'; penalty=25.
            limits=self.volume_limits.get(side,CylinderVolumeLimits(1.,2.))
            law,backend=build_side(settings,geometry,side,limits)
            stage='piston_topology'; penalty=20.
            q=law.value(np.linspace(0.,PERIOD,721,endpoint=False),1)
            count=zero_crossings(q)
            if count!=2:
                penalty+=abs(count-2); raise ValueError('Additional piston reversals.')
            stage='mechanical_constraints'; penalty=10.
            metrics=side_metrics(settings,law,backend,self.policy.mechanical_samples)
            constraints=[margin_record(r['metric'],metrics.get(r['metric']),r['limit'],r['relation'],r['unit']) for r in self.plan.request.mechanical_constraints]
            violations=[r for r in constraints if not r['satisfied']]
            if violations:
                penalty+=sum(1. if not r['available'] else min(1.,abs(r['margin'])/max(abs(r['limit']),1.)) for r in violations)
                raise ValueError('Declared mechanical constraint fails.')
            stage='dense_fit'; penalty=5.
            angles,weights=self.mesh[side]; target,velocity=self.references[side]
            error=(law.value(angles)-limits.minimum)/limits.swept-target
            position=float(weights@(error**2)); velocity_mse=None
            if velocity is not None:
                scale=max(float(np.sqrt(weights@(velocity**2))),1/PERIOD)
                dv=(law.value(angles,1)/limits.swept-velocity)/scale
                velocity_mse=float(weights@(dv**2))
            velocity_term=0. if velocity_mse is None else self.policy.velocity_weight*velocity_mse/(1+velocity_mse)
            residual=np.sqrt(weights)*error
            if velocity_mse is not None:
                residual=np.r_[residual,math.sqrt(self.policy.velocity_weight/(1+velocity_mse))*np.sqrt(weights)*dv]
            score=position+velocity_term
            row.update(valid=True,score=score,residual=residual,terms=dict(weighted_position_mse=position,normalized_velocity_mse=velocity_mse,
                velocity_term=velocity_term),metrics=metrics,reason=None)
            if retain: self.archive.add(row)
        except (ValueError,ArithmeticError,np.linalg.LinAlgError) as error:
            row.update(score=penalty,reason=stage)
            self.rejections[stage]+=1
        if len(self.cache)>=1024: self.cache.popitem(last=False)
        self.cache[key]=row
        return row

    def initial(self,rng,categories,index,side):
        x=rng.uniform(0.,1.,len(self.names))
        if index==0:
            p=dict(rod_over_crank=4.,offset_over_crank=0.,phase_rad=0.) if self.family=='slider_crank' else dict(
                coupler=3.5,rocker=2.5,ground_x=3.,ground_y=0.,output_along=2.5,output_normal=.3,
                rod_length=8.,slider_origin_x=0.,slider_origin_y=0.,axis_angle=0.,phase_rad=0.)
            x=np.clip(self.encode(p),0.,1.)
            # Align the generic valid seed's maximum; no source geometry is used.
            seed=self.decode(x,categories)
            settings=dict(family=self.family,**({'output':categories['output']} if self.family=='four_bar' else {}))
            try:
                law,_=build_side(settings,seed,side,CylinderVolumeLimits(1.,2.))
                angles=np.linspace(0.,PERIOD,181,endpoint=False)
                current=angles[np.argmax(law.value(angles))]
                desired=self.targets[side].extrema_angles()['maximum']
                direction=categories.get('crank_direction',1)
                seed['phase_rad']+=direction*(current-desired)
                x=self.encode(seed)
                phase=self.names.index('phase_rad')
                if math.isclose(self.span[phase],PERIOD): x[phase]%=1.
                else: x[phase]=np.clip(x[phase],0.,1.)
            except ValueError: pass
        return x

    def discover(self):
        for category_index,categories in enumerate(self.choices):
            for island in range(self.policy.islands):
                for side in self.sides:
                    seed=int(np.random.SeedSequence([self.policy.seed,('small','large').index(side),category_index,island]).generate_state(1)[0])
                    rng=np.random.default_rng(seed)
                    self.tasks.append(dict(side=side,categories=categories,island=island,seed=seed,rng=rng,
                        population=[],scores=[],generation=0,evaluations=0))
        try:
            # Interleave populations across sides/categories/islands before DE.
            for index in range(self.policy.population):
                for task in self.tasks:
                    x=self.initial(task['rng'],task['categories'],index,task['side'])
                    row=self.evaluate(x,task['side'],task['categories'],dict(island=task['island'],seed=task['seed'],generation=-1))
                    task['population'].append(x); task['scores'].append(row['score']); task['evaluations']+=1
            for generation in range(self.policy.generations):
                for index in range(self.policy.population):
                    for task in self.tasks:
                        population=np.array(task['population'])
                        pool=[i for i in range(self.policy.population) if i!=index]
                        a,b,c=task['rng'].choice(pool,3,replace=False)
                        mutant=population[a]+task['rng'].uniform(.5,.9)*(population[b]-population[c])
                        mutant=np.where((mutant<0)|(mutant>1),task['rng'].uniform(size=len(self.names)),mutant)
                        mask=task['rng'].uniform(size=len(self.names))<.8; mask[task['rng'].integers(len(self.names))]=True
                        proposal=np.where(mask,mutant,population[index])
                        row=self.evaluate(proposal,task['side'],task['categories'],dict(island=task['island'],seed=task['seed'],generation=generation))
                        task['evaluations']+=1; task['generation']=generation+1
                        if row['score']<=task['scores'][index]:
                            task['population'][index]=proposal; task['scores'][index]=row['score']
            if self.policy.discovery_polish_evaluations:
                for row in list(self.archive.rows):
                    self.polish(row,self.policy.discovery_polish_evaluations)
        except SearchStopped: pass

    def polish(self,parent,maximum_evaluations):
        x0=self.encode(parent['geometry'])
        lower=np.maximum(0.,x0-self.policy.local_radius); upper=np.minimum(1.,x0+self.policy.local_radius)
        for i,name in enumerate(self.names):
            if name in PERIODIC_COORDINATES and math.isclose(self.span[i],PERIOD):
                lower[i]=x0[i]-self.policy.local_radius; upper[i]=x0[i]+self.policy.local_radius
        if not np.all(lower<x0) or not np.all(x0<upper):
            # Active global bounds remain valid starting points for bounded LS.
            x0=np.clip(x0,lower+1e-12,upper-1e-12)
        side,categories=parent['side'],parent['categories']
        angles,weights=self.mesh[side]
        residual_length=len(parent['residual'])
        origin=dict(parent.get('origin',{}),local_polish=True)
        def residual(x):
            row=self.evaluate(x,side,categories,origin)
            if not row['valid']: return np.full(residual_length,math.sqrt(row['score']/residual_length))
            return row['residual']
        try:
            least_squares(residual,x0,bounds=(lower,upper),method='trf',max_nfev=maximum_evaluations,
                          ftol=1e-10,xtol=1e-10,gtol=1e-10)
        except SearchStopped: raise

    def final_library(self):
        members=[]; from .synthesis import assess_mechanism
        for row in self.archive.rows:
            settings,geometry=row['settings'],row['geometry']
            artifact=self.plan.artifact(geometry,settings=settings,provenance=dict(target_hash=self.plan.target.content_hash,
                protocol=self.plan.request.protocol,stage=self.plan.request.stages[-1],search_policy=POLICY_VERSION,
                categories=row['categories'],origin=row['origin']))
            raw=artifact.scientific
            limits=self.volume_limits.get(row['side'],CylinderVolumeLimits(1.,2.))
            law,backend=build_side(settings,geometry,row['side'],limits)
            points=backend.stationary_points(samples=self.policy.root_samples) if self.family=='six_bar' else backend.stationary_points() if self.family=='slider_crank' else law.model.stationary_points(row['side'],samples=self.policy.root_samples)
            if len(points)!=2 or {p['kind'] for p in points}!={'maximum','minimum'}:
                self.rejections['final_refined_topology']+=1; continue
            evidence=assess_mechanism(artifact,self.plan.target,row['side'],mechanical_samples=max(self.policy.mechanical_samples,self.policy.root_samples),
                volume_limits=self.volume_limits.get(row['side']))
            if any(not r['satisfied'] for r in evidence['mechanical']['constraints']):
                self.rejections['final_mechanical_constraints']+=1; continue
            topology=dict(extrema=points,method='exact collinear dead centers' if self.family=='slider_crank' else
                'analytic velocity roots, doubled-grid and tangency screen; not a certified continuous proof',initial_samples=self.policy.root_samples)
            evidence['mechanical']['topology']=topology
            family_id=f"{row['side']}-family-{len(members)+1:03d}"
            metadata=dict(target_hash=self.plan.target.content_hash,evidence={row['side']:evidence},thermodynamic=None,
                search=dict(score=row['score'],terms=row['terms'],policy_version=POLICY_VERSION,bounds_version=BOUNDS_VERSION,
                    policy=asdict(self.policy),effective_bounds=self.bounds,stopping_reason=self.stop_reason,
                    evaluations=self.evaluations,cache_hits=self.cache_hits,rejections=dict(self.rejections),
                    islands=[dict(side=t['side'],categories=t['categories'],island=t['island'],seed=t['seed'],
                        generations=t['generation'],evaluations=t['evaluations']) for t in self.tasks]),
                provenance=dict(protocol=self.plan.request.protocol,stage=self.plan.request.stages[-1],origin=row['origin'],categories=row['categories']))
            members.append(dict(family_id=family_id,mechanisms={row['side']:artifact.data},metadata=metadata))
        if not members: raise ValueError('No mechanically admissible basin found under these bounds, constraints and search budget.')
        for member in members:
            member['metadata']['search']['rejections']=dict(self.rejections)
        return MechanismLibrary(tuple(json.loads(json.dumps(members))))


def execute_synthesis(plan,*,policy=None,sides=('large',),library=None,volume_limits=None):
    if plan.request.mechanism_family=='six_bar':
        from .synthesis_six_bar import execute_six_bar
        return execute_six_bar(plan,policy=policy,sides=sides,library=library,volume_limits=volume_limits)
    if plan.request.stages not in (('global_discovery',),('full_local_polish',)):
        raise NotImplementedError('Execute one geometric discovery or polish stage; thermodynamic stages require the actual Research solver.')
    if not sides or len(set(sides))!=len(sides) or set(sides)-{'small','large'}:
        raise ValueError('Choose distinct cylinder sides.')
    policy=policy or SearchPolicy()
    if not isinstance(policy,SearchPolicy): raise ValueError('Synthesis requires a SearchPolicy.')
    engine=GeometrySearch(plan,policy,sides,volume_limits)
    if plan.request.protocol!='design_exploitation':
        raise NotImplementedError('Fresh primary-loop saturation remains a distinct unimplemented operator.')
    if plan.request.stages==('global_discovery',):
        if library is not None or plan.request.retained_family_ids: raise ValueError('Discovery does not load retained candidates.')
        engine.discover()
    else:
        if not isinstance(library,MechanismLibrary) or not plan.request.retained_family_ids:
            raise ValueError('Local polish requires a library and explicit retained family IDs.')
        selected=[library.member(identifier) for identifier in plan.request.retained_family_ids]
        if any(sum(side in member['mechanisms'] for member in selected)>policy.retain_per_side for side in sides):
            raise ValueError('Polish retention must cover every selected parent on each piston side.')
        try:
            for family_id in plan.request.retained_family_ids:
                member=library.member(family_id)
                if member['metadata'].get('target_hash')!=plan.target.content_hash:
                    raise ValueError('Polish library and target identities disagree.')
                if not any(side in member['mechanisms'] for side in sides):
                    raise ValueError('Selected family has no mechanism on the requested piston sides.')
                for side in sides:
                    if side not in member['mechanisms']: continue
                    artifact=MechanismArtifact.from_data(member['mechanisms'][side])
                    raw=artifact.scientific
                    if raw['settings']['family']!=plan.request.mechanism_family: raise ValueError('Polish cannot change mechanism family.')
                    categories={name:raw['settings']['output'] if name=='output' else raw['geometry'][name] for name in CATEGORY_VALUES[engine.family]}
                    if categories not in engine.choices: raise ValueError('Polish categories are excluded by this search policy.')
                    if any(c not in plan.request.mechanical_constraints for c in raw['constraints']):
                        raise ValueError('Polish must retain all source artifact constraints.')
                    x=engine.encode(raw['geometry'])
                    if np.any(x<0) or np.any(x>1): raise ValueError('Polish geometry lies outside configured search bounds.')
                    parent=engine.evaluate(x,side,categories,dict(parent_family_id=family_id,parent_artifact_hash=artifact.content_hash,source_settings=raw['settings']))
                    if not parent['valid']: raise ValueError('Selected parent fails current synthesis constraints/topology.')
                    engine.polish(parent,policy.polish_evaluations)
        except SearchStopped: pass
    if plan.request.stages==('full_local_polish',):
        best={}
        for row in engine.archive.rows:
            key=(row['side'],row['origin']['parent_family_id'])
            if key not in best or row['score']<best[key]['score']: best[key]=row
        engine.archive.rows=list(best.values())
    return engine.final_library()
