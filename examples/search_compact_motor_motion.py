"""Fit complete reference slider mechanisms to the motor target, without physics.

Independent side searches preserve the reference branch and output family.
Compactness and transmission floors are explicit, reference-relative screens.
They are sampled geometric screens, not collision or structural certification.
"""
import argparse
from dataclasses import asdict, replace
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import differential_evolution

from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_model
from dada_solver.four_bar import CouplerOutputPoint

ROOT = Path(__file__).resolve().parents[1]


def evaluate(assembly, radius, direction, offset, angles):
    axis = np.array([math.cos(assembly.slider.axis_angle), math.sin(assembly.slider.axis_angle)])
    origin = np.array([assembly.slider.axis_origin_x, assembly.slider.axis_origin_y])
    states = [assembly.evaluate(radius*np.array([math.cos(t), math.sin(t)]),
              direction*radius*np.array([-math.sin(t), math.cos(t)]))
              for t in direction*angles+offset]
    x = np.array([s.coordinate for s in states])
    stroke = np.ptp(x)
    if stroke <= 1e-12:
        raise ValueError('Zero stroke')
    points = [[0., 0.], [assembly.loop.rocker_pivot_x, assembly.loop.rocker_pivot_y]]
    for s in states:
        points.extend([s.crank_pin, s.coupler_joint, s.output_point, origin+axis*s.coordinate])
    envelope = float(np.linalg.norm(np.ptp(points, axis=0)))
    metrics = dict(stroke_m=float(stroke), envelope_diagonal_m=envelope,
        stroke_to_envelope=float(stroke/envelope),
        minimum_transmission_sine=min(abs(s.four_bar_cross_product) for s in states)/(assembly.loop.coupler_length*assembly.loop.rocker_length),
        minimum_rod_axis_cosine=min(s.connecting_rod_transverse_margin for s in states)/assembly.slider.connecting_rod_length)
    return (x-x.min())/stroke, np.array([s.coordinate_derivative for s in states])/stroke, metrics


def candidate(base, vector, family="rocker"):
    if family not in ("rocker", "coupler"):
        raise ValueError("Unknown output family")
    a,b,c,d,e,f,phase = vector
    output = (CouplerOutputPoint(c*base.loop.coupler_length, d*base.loop.coupler_length)
              if family == 'coupler' else replace(base.output,
              along_rocker=base.output.along_rocker*c, normal_to_rocker=base.output.normal_to_rocker*d))
    return replace(base,
        loop=replace(base.loop,coupler_length=base.loop.coupler_length*a,rocker_length=base.loop.rocker_length*b),
        output=output,
        slider=replace(base.slider,connecting_rod_length=base.slider.connecting_rod_length*e,
                       axis_origin_y=base.slider.axis_origin_y+f)), phase


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--family', choices=['rocker','coupler'], default='rocker')
    ap.add_argument('--iterations',type=int,default=25)
    ap.add_argument('--seed',type=int,default=27)
    ap.add_argument('--output',type=Path,default=ROOT/'outputs/compact_motor_motion_search.json')
    args=ap.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    start=time.monotonic()
    reference=ROOT/'examples/cooling_cell_x2_reference.toml'
    kin=build_model(load_simulation_configuration(reference)).kinematics
    target_path=ROOT/'outputs/motor_champion_motion_target.csv'
    target=np.genfromtxt(target_path,delimiter=',',names=True)[::12]
    angles=target['theta_rad']
    # Drop the duplicated endpoint for uniform least-squares weighting.
    if angles[-1] >= 2*math.pi-1e-9:
        target=target[:-1]; angles=angles[:-1]
    report=dict(reference=str(reference.relative_to(ROOT)),target=str(target_path.relative_to(ROOT)),
        family=args.family,seed=args.seed,iterations=args.iterations,derivative_weight=.1,
        parameter_names=['coupler_multiplier','rocker_multiplier','output_along_multiplier',
                         'output_normal_multiplier','rod_multiplier','slider_axis_y_shift_m','phase_shift_rad'],
        limitations=['Sampled geometric constraints','No thickness, collisions, bearing loads or thermodynamic evaluation'],sides={})
    for name in ['small','large']:
        base=getattr(kin,name+'_assembly')
        if args.family == 'coupler':
            report['parameter_names'][2:4] = ['output_along_over_reference_coupler','output_normal_over_reference_coupler']
        sign=1 if getattr(kin,name+'_volume_increases_with_coordinate') else -1
        q0,dq0,m0=evaluate(base,kin.crank_radius,kin.crank_direction,kin.crank_angle_offset,angles)
        tq=target[name+'_fraction_0_1']; tdq=target[name+'_dq_dtheta_per_rad']
        floors={key:.95*m0[key] for key in ['stroke_to_envelope','minimum_transmission_sine','minimum_rod_axis_cosine']}
        archive=[]
        def score(v, keep=True):
            assembly,phase=candidate(base,v,args.family)
            try:
                q,dq,m=evaluate(assembly,kin.crank_radius,kin.crank_direction,kin.crank_angle_offset+phase,angles)
            except ValueError:
                return 1e6
            if any(m[k]<value for k,value in floors.items()):
                return 1e3+sum(max(0,value-m[k]) for k,value in floors.items())
            if sign<0: q=1-q; dq=-dq
            pos=float(np.sqrt(np.mean((q-tq)**2))); der=float(np.sqrt(np.mean((dq-tdq)**2)))
            value=math.sqrt(pos**2+.1*der**2)
            if keep and (not archive or value<archive[-1]['score']):
                archive.append(dict(score=value,position_rms=pos,derivative_rms=der,parameters=list(map(float,v)),metrics=m))
            return value
        bounds=[(.75,1.25)]*4+[(.65,1.15),(-.25*kin.crank_radius,.25*kin.crank_radius),(-math.pi,math.pi)]
        initial=[1,1,1,1,1,0,0]
        if args.family == 'coupler':
            bounds[2:4] = [(-1.5,1.5),(-1.5,1.5)]
            initial=[1,1,.5,0,1,0,0]
        report['bounds']=bounds
        score(initial)
        differential_evolution(score,bounds,seed=args.seed,popsize=5,maxiter=args.iterations,polish=False,x0=initial)
        if not archive:
            report['sides'][name] = dict(status='no_feasible_candidate', constraint_floors=floors)
            args.output.write_text(json.dumps(report,indent=2)+'\n')
            continue
        best=min(archive,key=lambda x:x['score'])
        assembly,phase=candidate(base,best['parameters'],args.family)
        _,_,dense=evaluate(assembly,kin.crank_radius,kin.crank_direction,kin.crank_angle_offset+phase,np.linspace(0,2*math.pi,2880,endpoint=False))
        report['sides'][name]=dict(reference_metrics=m0,constraint_floors=floors,best=best,
            dense_validation_metrics=dense,dense_screen_passed=all(dense[k]>=v for k,v in floors.items()),
            assembly=asdict(assembly),crank_radius=kin.crank_radius,crank_direction=kin.crank_direction,
            crank_angle_offset=kin.crank_angle_offset+phase,volume_increases_with_coordinate=sign>0,
            improvements=archive)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(name,best,flush=True)
    report['elapsed_seconds']=time.monotonic()-start
    args.output.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
