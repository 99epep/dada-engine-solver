"""Deterministic position-led initialization of the structured C2 pair."""
from dataclasses import dataclass, field, asdict
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.interpolate import PPoly
from scipy.optimize import least_squares

from dada_solver.structured_kinematics import StructuredMotion15, StructuredKinematics15
from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.campaign.candidate import canonical_json, content_hash
from dada_solver.campaign.history import atomic_json
from .families import STRUCTURED_DEFAULTS
from .motion_target import MotionTarget, PERIOD, PeriodicTargetSide

NAMES = tuple(STRUCTURED_DEFAULTS)
POLICY_VERSION = 'structured_c2_position_multistart_v1'


@dataclass(frozen=True)
class RefitPolicy:
    """Numerical search settings, not physical domains."""
    dense_samples: int = 1440
    maximum_evaluations: int = 200
    velocity_weight: float = 1e-6
    kink_width_starts: tuple = (.25, .65)
    curvature_scales: tuple = (.7, 1., 1.4)

    def __post_init__(self):
        for name, minimum in (('dense_samples',360), ('maximum_evaluations',30)):
            if type(getattr(self,name)) is not int or getattr(self,name)<minimum:
                raise ValueError(f'{name} must be an integer >= {minimum}.')
        if not math.isfinite(self.velocity_weight) or not 0<=self.velocity_weight<=1e-4:
            raise ValueError('Velocity weight must lie in [0, 1e-4].')
        if not self.kink_width_starts or any(not 0<w<2 for w in self.kink_width_starts):
            raise ValueError('Kink-width starts must lie in (0, 2).')
        if not self.curvature_scales or any(not math.isfinite(s) or s<=0 for s in self.curvature_scales):
            raise ValueError('Curvature scales must be finite and positive.')


def fit_bounds(seed):
    """Finite numerical fitting boxes inside the existing structural domains."""
    lo=[]; hi=[]
    for n in NAMES:
        if n=='small_max_deg': a,b=seed[n]-180.,seed[n]+180.
        elif n.endswith('duration_deg'): a,b=.1,359.9
        elif n.endswith('curvature'): a,b=1e-5,max(5000.,seed[n]*10.)
        elif n.endswith('width_rel'): a,b=.005,1.995
        else: a,b=.001,.999
        lo.append(a);hi.append(b)
    return np.array(lo),np.array(hi)


def structured_model(parameters):
    limits=CylinderVolumeLimits(1.,2.)
    return StructuredKinematics15(SimpleNamespace(small_cylinder=limits,large_cylinder=limits),parameters)


def monotonicity(motion):
    """Minimum branch derivative at all polynomial critical points and endpoints.

    This checks each complete continuous branch, rather than a velocity grid.
    Each increasing helper has unit net rise and positive endpoint curvature;
    nonnegative derivatives guarantee exactly the two physical reversals.
    """
    minima={}
    for name in ('small_down','small_up','large_down','large_up'):
        derivative=PPoly.from_bernstein_basis(getattr(motion,name)).derivative()
        roots=derivative.derivative().roots(extrapolate=False)
        points=np.r_[derivative.x,roots[np.isfinite(roots)&(roots>=0)&(roots<=1)]]
        minima[name]=float(np.min(derivative(points)))
    return dict(valid=all(v>=-1e-9 for v in minima.values()),minimum_branch_derivatives=minima,
                method='piecewise polynomial derivative minima at endpoints and all critical roots')


def _seed(target, sources, supplied):
    p=dict(STRUCTURED_DEFAULTS)
    for side in ('small','large'):
        source=sources[side];ext=source.extrema_angles()
        maximum_t=(-ext['maximum']/PERIOD)%1.
        duration=((ext['maximum']-ext['minimum'])%PERIOD)/PERIOD
        p[side+'_down_duration_deg']=360*duration
        if side=='small':p['small_max_deg']=360*maximum_t
        # Curvature is a local position-shape seed only, never an acceleration
        # target in the objective. Time here is normalized motor time.
        h=PERIOD/512
        for kind,sign in (('maximum',-1.),('minimum',1.)):
            a=ext[kind]
            value=sign*(source.position(a+h)+source.position(a-h)-2*source.position(a))/(h/PERIOD)**2
            p[side+('_max_curvature' if kind=='maximum' else '_min_curvature')]=max(.1,float(value))
        if side=='small':
            midpoint=ext['minimum']-PERIOD*(1-duration)/2
            p['small_up_bp_mid_q']=float(source.position(midpoint))
        else:
            p['large_down_bp_mid_q']=float(1-source.position(-PERIOD*duration/2))
        kink=next((e for e in source.events if e['kind']=='kink'),None)
        if kink is not None:
            start=ext['maximum'] if side=='small' else ext['minimum']
            width=PERIOD*(duration if side=='small' else 1-duration)
            prefix='small_down' if side=='small' else 'large_up'
            p[prefix+'_kink_u']=float(((start-kink['angle_rad'])%PERIOD)/width)
            p[prefix+'_kink_q']=float(1-source.position(kink['angle_rad']) if side=='small' else source.position(kink['angle_rad']))
    for n,v in (supplied or {}).items():
        if n in p:p[n]=float(v)
    lo,hi=fit_bounds(p)
    return dict(zip(NAMES,np.clip([p[n] for n in NAMES],lo+1e-8,hi-1e-8)))


def _mesh(source,policy):
    base=np.linspace(0,PERIOD,policy.dense_samples,endpoint=False)
    features=[]
    for e in source.events:
        # Nine points within 1/128 cycle of a feature; rounding also contributes
        # its directed interior. Background always carries 90% of position mass.
        points=e['angle_rad']+np.linspace(-PERIOD/128,PERIOD/128,9)
        if e['kind']=='rounding' and e['end_angle_rad'] is not None:
            points=np.r_[points,e['angle_rad']+np.linspace(0,(e['end_angle_rad']-e['angle_rad'])%PERIOD,9)]
        features.extend(points%PERIOD)
    angles=np.r_[base,features]
    weights=np.r_[np.full(len(base),(.9 if features else 1.)/len(base)),
                  np.full(len(features),.1/len(features)) if features else []]
    return angles,np.sqrt(weights)


def _diagnostics(model,side,source,policy):
    theta=np.linspace(0,PERIOD,4*policy.dense_samples,endpoint=False)
    q=getattr(model,side+'_cylinder_volume')(theta)-1
    v=getattr(model,side+'_cylinder_volume_derivative')(theta)
    error=q-source.position(theta)
    p=model.params;maxangle=(-math.radians(p['small_max_deg']))%PERIOD if side=='small' else 0.
    minangle=(maxangle-math.radians(p[side+'_down_duration_deg']))%PERIOD
    expected=source.extrema_angles()
    extrema=[dict(kind=k,angle_rad=a,position=1. if k=='maximum' else 0.) for k,a in (('maximum',maxangle),('minimum',minangle))]
    events=[dict(kind=e['kind'],angle_rad=e['angle_rad'],position_error=float(getattr(model,side+'_cylinder_volume')(e['angle_rad'])-1-source.position(e['angle_rad']))) for e in source.events]
    return dict(position_rms=float(np.sqrt(np.mean(error**2))),maximum_absolute_position_error=float(np.max(abs(error))),
        velocity_rms=None if not source.has_velocity else float(np.sqrt(np.mean((v-source.velocity(theta))**2))),
        extrema=extrema,extrema_count=2 if monotonicity(model.motion)['valid'] else None,extrema_angular_errors_rad={k:float((a-expected[k]+math.pi)%PERIOD-math.pi) for k,a in (('maximum',maxangle),('minimum',minangle))},
        events=events,target_acceleration_used=False)


@dataclass(frozen=True)
class MotionRefitResult:
    payload_json: str
    @property
    def data(self):return json.loads(self.payload_json)
    @property
    def content_hash(self):return self.data['content_hash']
    def save(self,path):
        path=Path(path)
        if path.exists():raise ValueError('Refit report already exists; choose a new path.')
        path.parent.mkdir(parents=True,exist_ok=True);atomic_json(path,self.data)
    @classmethod
    def load(cls,path):return cls.from_data(json.loads(Path(path).read_text()))
    @classmethod
    def from_data(cls,data):
        if (set(data)!={'artifact_type','schema_version','scientific','content_hash'} or data['artifact_type']!='motion_refit'
            or type(data['schema_version']) is not int or data['schema_version']!=2
            or data['scientific']['destination_family']!='structured_c2_15p' or content_hash(data['scientific'])!=data['content_hash']):
            raise ValueError('Invalid structured refit report or content hash.')
        return cls(canonical_json(data))


@dataclass(frozen=True)
class MotionRefitRequest:
    target: MotionTarget
    destination_family: str = 'structured_c2_15p'
    position_role: str = 'primary_synthesis_reference'
    acceleration_role: str = 'diagnostic_only'
    policy: RefitPolicy = field(default_factory=RefitPolicy)
    initial_parameters: dict | None = None

    def __post_init__(self):
        if not isinstance(self.target,MotionTarget):raise ValueError('Refit requires a MotionTarget.')
        if self.destination_family!='structured_c2_15p':raise ValueError('Refit destination is structured_c2_15p.')
        if (self.position_role,self.acceleration_role)!=('primary_synthesis_reference','diagnostic_only'):
            raise ValueError('Refit uses position as reference and acceleration as diagnostic only.')
        if not isinstance(self.policy,RefitPolicy):raise ValueError('Refit requires a valid numerical policy.')

    def execute(self):
        sources={s:PeriodicTargetSide(self.target,s) for s in ('small','large')}
        seed=_seed(self.target,sources,self.initial_parameters)
        low,high=fit_bounds(seed);meshes={s:_mesh(source,self.policy) for s,source in sources.items()}
        def residual(x):
            model=structured_model(dict(zip(NAMES,x)));parts=[]
            for side,source in sources.items():
                angles,w=meshes[side]
                parts.append(w*(getattr(model,side+'_cylinder_volume')(angles)-1-source.position(angles)))
                if source.has_velocity and self.policy.velocity_weight:
                    parts.append(math.sqrt(self.policy.velocity_weight)*w*PERIOD*(getattr(model,side+'_cylinder_volume_derivative')(angles)-source.velocity(angles)))
            parts.append(10*model.motion.monotonicity_penalty(n=48)/math.sqrt(192))
            return np.concatenate(parts)
        starts=[];accepted=[]
        for width in self.policy.kink_width_starts:
            for scale in self.policy.curvature_scales:
                p=dict(seed)
                for n in NAMES:
                    if n.endswith('width_rel'):p[n]=width
                    if n.endswith('curvature'):p[n]*=scale
                x=np.clip([p[n] for n in NAMES],low+1e-8,high-1e-8)
                fit=least_squares(residual,x,bounds=(low,high),x_scale='jac',ftol=1e-10,xtol=1e-10,gtol=1e-10,max_nfev=self.policy.maximum_evaluations)
                parameters=dict(zip(NAMES,map(float,fit.x)));model=structured_model(parameters);mono=monotonicity(model.motion)
                score=float(np.dot(residual(fit.x),residual(fit.x)))
                starts.append(dict(initial_parameters=dict(zip(NAMES,map(float,x))),parameters=parameters,converged=bool(fit.success),
                    status=int(fit.status),message=fit.message,evaluations=int(fit.nfev),monotonicity=mono,least_squares_cost=score))
                if mono['valid']:accepted.append((score,len(starts)-1,parameters))
        if not accepted:raise ValueError('No monotone structured C2 fit found; no study was produced.')
        _,chosen,parameters=min(accepted,key=lambda r:(r[0],r[1]));model=structured_model(parameters)
        sides={s:dict(diagnostics=_diagnostics(model,s,source,self.policy),initial_diagnostics=_diagnostics(structured_model(seed),s,source,self.policy)) for s,source in sources.items()}
        x=np.array([parameters[n] for n in NAMES])
        scientific=dict(target_hash=self.target.content_hash,target=self.target.data,destination_family=self.destination_family,
            parameter_count=15,initial_parameters=seed,parameters=parameters,starts=starts,selected_start=chosen,
            convergence=starts[chosen]['converged'],monotonicity=monotonicity(model.motion),sides=sides,
            combined=dict(position_rms=math.sqrt(sum(s['diagnostics']['position_rms']**2 for s in sides.values())/2)),
            near_bounds=[n for i,n in enumerate(NAMES) if min(x[i]-low[i],high[i]-x[i])/(high[i]-low[i])<1e-3],
            fit_bounds={n:[float(low[i]),float(high[i])] for i,n in enumerate(NAMES)},
            angle_domain=self.target.scientific['angle_domain'],policy_version=POLICY_VERSION,policy=asdict(self.policy),
            position_role=self.position_role,acceleration_role=self.acceleration_role,thermodynamic_evaluation=False)
        return MotionRefitResult(canonical_json(dict(artifact_type='motion_refit',schema_version=2,scientific=scientific,content_hash=content_hash(scientific))))


def plot_refit(result):
    """Position, error and optional velocity in the current study-angle convention."""
    import matplotlib.pyplot as plt
    raw=result.data['scientific'];target=MotionTarget.from_data(raw['target']);model=structured_model(raw['parameters'])
    angles=np.linspace(0,PERIOD,1441);fig,axes=plt.subplots(2,3,figsize=(15,8),sharex=True)
    for row,side in enumerate(('small','large')):
        source=PeriodicTargetSide(target,side);q=getattr(model,side+'_cylinder_volume')(angles)-1
        axes[row,0].plot(angles,source.position(angles),'k--',label='Target');axes[row,0].plot(angles,q,label='Structured refit')
        axes[row,1].plot(angles,q-source.position(angles),label='Position error')
        if source.has_velocity:axes[row,2].plot(angles,source.velocity(angles),'k--',label='Target velocity')
        axes[row,2].plot(angles,getattr(model,side+'_cylinder_volume_derivative')(angles),label='Structured velocity')
        for ax in axes[row]:
            seen=set()
            for event in source.events:
                label=event['kind'] if ax is axes[row,0] and event['kind'] not in seen else None
                ax.axvline(event['angle_rad'],color='gray',alpha=.25,label=label)
                seen.add(event['kind'])
            for kind,angle in source.extrema_angles().items():ax.axvline(angle,color='black',ls=':',alpha=.4)
            for e in raw['sides'][side]['diagnostics']['extrema']:ax.axvline(e['angle_rad'],color='tab:red',ls='--',alpha=.4)
            ax.legend();ax.set_xlabel('Study angle [rad]');ax.grid(alpha=.2)
        axes[row,0].set_ylabel(side.upper()+' normalized position')
    fig.suptitle('structured_c2_15p refit — no thermodynamic evaluation');fig.tight_layout()
    return fig
