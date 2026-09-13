"""Independent small-cylinder E-output synthesis with finite slider rod.

Crank radius is the length unit. Ground orientation fixes the rotational gauge;
slider orientation and phase remain independent. No thermodynamic integration.
"""
import argparse
from dataclasses import asdict
import json
import math
from pathlib import Path
import time
import numpy as np
from scipy.optimize import differential_evolution
from dada_solver.four_bar import FourBarLoop, RockerOutputPoint, CouplerOutputPoint, SliderConstraint, FourBarSliderAssembly
from search_compact_motor_motion import ROOT, evaluate

NAMES=['ground','coupler','rocker','output_radius','output_angle','rod','axis_offset','axis_angle','phase']
BOUNDS=[(1.5,5.),(1.2,6.),(1.2,6.),(.5,6.),(-math.pi,math.pi),(2.,12.),(-4.,4.),(-math.pi,math.pi),(-math.pi,math.pi)]
FLOORS={'stroke_to_envelope':.12384,'minimum_transmission_sine':.48733,'minimum_rod_axis_cosine':.94927}


def build(v, branch, family="E"):
    if family not in ("E", "F"):
        raise ValueError("Output family must be E or F")
    g,c,r,out,beta,rod,h,axis,phase=v
    assembly=FourBarSliderAssembly(FourBarLoop(c,r,g,0.,branch),
        (RockerOutputPoint if family == "E" else CouplerOutputPoint)(out*math.cos(beta),out*math.sin(beta)),
        SliderConstraint(-h*math.sin(axis),h*math.cos(axis),axis,rod,1))
    return assembly,phase


def assess(v, branch, sign, theta, target, derivative, family="E"):
    assembly,phase=build(v,branch,family)
    # Exact full-revolution loop closure before sampling slider motion.
    ground=v[0]
    if abs(v[1]-v[2]) >= ground-1 or v[1]+v[2] <= ground+1:
        raise ValueError('Loop cannot close throughout the revolution')
    q,dq,metrics=evaluate(assembly,1.,1,phase,theta)
    metrics['stroke_over_crank']=metrics.pop('stroke_m')
    metrics['envelope_diagonal_over_crank']=metrics.pop('envelope_diagonal_m')
    if sign<0:q=1-q;dq=-dq
    pos=float(np.sqrt(np.mean((q-target)**2)));vel=float(np.sqrt(np.mean((dq-derivative)**2)))
    margins={k:metrics[k]-floor for k,floor in FLOORS.items()}
    return dict(position_rms=pos,derivative_rms=vel,score=math.sqrt(pos*pos+.1*vel*vel),
                metrics=metrics,margins=margins),q,dq


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--iterations',type=int,default=160)
    parser.add_argument('--family',choices=['E','F'],default='E')
    parser.add_argument('--restarts',type=int,default=1)
    parser.add_argument('--seed',type=int,default=127)
    parser.add_argument('--warm-start',type=Path)
    parser.add_argument('--output',type=Path,default=ROOT/'outputs/small_rocker_e_search.json')
    args=parser.parse_args()
    if args.iterations<1 or args.restarts<1: parser.error('Positive iteration and restart counts required')
    target_path=ROOT/'outputs/motor_champion_motion_target.csv'
    raw=np.genfromtxt(target_path,delimiter=',',names=True)[:-1]
    coarse=raw[::12];dense=raw
    output=args.output
    if args.family == 'F' and output.name == 'small_rocker_e_search.json':
        output=output.with_name('small_coupler_f_search.json')
    previous=json.loads(args.warm_start.read_text()) if args.warm_start else None
    if previous and previous['family'] != ('rocker_E' if args.family == 'E' else 'coupler_F'):
        parser.error('Warm-start geometry must belong to the selected output family')
    output.parent.mkdir(parents=True,exist_ok=True)
    report=dict(family=('rocker_E' if args.family == 'E' else 'coupler_F'),side='small',length_unit='crank_radius',bounds=dict(zip(NAMES,BOUNDS)),
        floors=FLOORS,derivative_weight=.1,iterations=args.iterations,restarts=args.restarts,warm_start=previous,runs=[],
        limitations=['Sampled slider clearance and transmission screens','Cylinder collisions not evaluated','No thermodynamic evaluation'])
    start=time.monotonic()
    for index,(branch,sign) in enumerate([(1,1),(1,-1),(-1,1),(-1,-1)]*args.restarts):
        archive=[]
        def objective(v):
            try:r,_,_=assess(v,branch,sign,coarse['theta_rad'],coarse['small_fraction_0_1'],coarse['small_dq_dtheta_per_rad'],args.family)
            except ValueError:return 1e4
            violation=sum(max(0,-m) for m in r['margins'].values())
            if violation:return 10+violation
            if not archive or r['score']<archive[-1]['score']:
                archive.append(dict(r,parameters=list(map(float,v))))
            return r['score']
        options={}
        if previous:
            matches=[r for r in previous['runs'] if r['branch']==branch and r['volume_sign']==sign and r['best']]
            if matches:
                seed_point=min(matches,key=lambda r:r['best']['score'])['best']['parameters']
                objective(seed_point)
                rng=np.random.default_rng(args.seed+index)
                lo=np.array([b[0] for b in BOUNDS]);hi=np.array([b[1] for b in BOUNDS])
                population=np.clip(np.array(seed_point)+rng.normal(0,.08,(72,9))*(hi-lo),lo,hi)
                population[0]=seed_point
                population[-18:]=rng.uniform(lo,hi,(18,9))
                options['init']=population
        result=differential_evolution(objective,BOUNDS,seed=args.seed+index,popsize=8,maxiter=args.iterations,polish=False,tol=1e-8,**options)
        validated=[]
        for candidate in sorted(archive,key=lambda r:r['score'])[:12]:
            try:r,_,_=assess(candidate['parameters'],branch,sign,dense['theta_rad'],dense['small_fraction_0_1'],dense['small_dq_dtheta_per_rad'],args.family)
            except ValueError:continue
            if all(m>=0 for m in r['margins'].values()):validated.append(dict(r,parameters=candidate['parameters']))
        best=min(validated,key=lambda r:r['score']) if validated else None
        run=dict(seed=args.seed+index,branch=branch,volume_sign=sign,nfev=result.nfev,best=best)
        if best:
            assembly,phase=build(best['parameters'],branch,args.family)
            run.update(assembly=asdict(assembly),phase=phase)
        report['runs'].append(run)
        winners=[r for r in report['runs'] if r['best']]
        report['best']=min(winners,key=lambda r:r['best']['score']) if winners else None
        report['elapsed_seconds']=time.monotonic()-start
        output.write_text(json.dumps(report,indent=2)+'\n')
        print(index,'best',best,flush=True)
    if report['best']:
        plot(report,raw,output.with_suffix(".png"))


def plot(report,raw,output=None):
    import os,tempfile
    os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'dada_solver_matplotlib'))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
    for run in report['runs']:
        if not run['best']:continue
        _,q,dq=assess(run['best']['parameters'],run['branch'],run['volume_sign'],raw['theta_rad'],raw['small_fraction_0_1'],raw['small_dq_dtheta_per_rad'], 'F' if report['family']=='coupler_F' else 'E')
        label=f"{report['family']} fit, branch {run['branch']:+d}, volume sign {run['volume_sign']:+d}"
        axes[0].plot(raw['theta_deg'],q,label=label);axes[1].plot(raw['theta_deg'],dq)
    for ax,key in zip(axes,['small_fraction_0_1','small_dq_dtheta_per_rad']):
        ax.plot(raw['theta_deg'],raw[key],'k--',label='Free-motion target',lw=2);ax.grid(alpha=.2)
    axes[0].legend(fontsize=8);axes[0].set_ylabel('Normalized volume')
    axes[1].set_ylabel('Derivative (rad⁻¹)');axes[1].set_xlabel('Study angle (deg)')
    fig.suptitle(f"Independent small-cylinder {report['family']} synthesis — finite piston rod")
    fig.tight_layout();fig.savefig(output or ROOT/'outputs/small_rocker_e_search.png',dpi=140)

if __name__=='__main__':main()
