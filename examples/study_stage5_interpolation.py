"""Bounded geometry and full thermodynamic gates after exact Stage-5 changes."""
import argparse
from dataclasses import replace, asdict
import json
from pathlib import Path
from types import SimpleNamespace
import time
import numpy as np
from scipy.optimize import brentq
from benchmark_solver_acceleration import ROOT, DEFAULT_MANIFEST, digest, environment
from refine_motor_four_stage_hx9d_variable_gas import _load_basis,_build_design
from compare_motor_motion_laws_stage7A5 import _evaluate
from dada_solver.exchangers.wall_cycle import WallCycleNumericalSettings,wall_cycle_performance
from dada_solver.wall_backend import WallRHS, WallBackendSettings
from dada_solver.campaign.evaluator import json_values
from dada_solver.state import ThermodynamicState
from experimental_periodic_interpolation import PeriodicVolumeInterpolant,TAU
from stage5_mechanisms import six_bar


def compare_fields(a,b,path='',differences=None):
    """Strict original 1e-6/1e-12 scale for all exposed numerical report fields."""
    if differences is None:differences=[]
    if isinstance(a,dict) and isinstance(b,dict):
        if a.keys()!=b.keys():differences.append(dict(path=path,reason='keys'))
        for k in a.keys()&b.keys():compare_fields(a[k],b[k],path+'.'+k,differences)
    elif isinstance(a,(list,tuple)) and isinstance(b,(list,tuple)):
        if len(a)!=len(b):differences.append(dict(path=path,reason='length'))
        for i,(x,y) in enumerate(zip(a,b)):compare_fields(x,y,f'{path}[{i}]',differences)
    elif isinstance(a,(float,int)) and not isinstance(a,bool) and isinstance(b,(float,int)):
        if a==b:return differences
        distance=abs(a-b)/(1e-12+1e-6*max(abs(a),abs(b)))
        if not np.isfinite(distance) or distance>1:differences.append(dict(path=path,reference=a,actual=b,scaled_distance=distance))
    elif a!=b:differences.append(dict(path=path,reference=a,actual=b))
    return differences


def geometry_study(exact,resolutions):
    # Prime-sized half-offset grid independent of all power-of-two construction grids.
    dense=(np.arange(32771)+.5)*TAU/32771
    roots=[]
    knots=np.linspace(0,TAU,1441)
    velocities=np.array([exact.cylinder_volumes_and_derivatives(float(a))[2:] for a in knots])
    for j in (0,1):
        for a,b,da,db in zip(knots[:-1],knots[1:],velocities[:-1,j],velocities[1:,j]):
            if da*db<0:roots.append(brentq(lambda t:exact.cylinder_volumes_and_derivatives(t)[2+j],a,b))
    boundary=np.array([0.,1e-12,1e-8,TAU-1e-8,TAU-1e-12,TAU])
    dense_values=np.array([exact.cylinder_volumes_and_derivatives(float(a)) for a in dense])
    curvature=np.max(np.abs(np.gradient(dense_values[:,2:],dense,axis=0)),axis=1)
    focus=dense[np.argsort(curvature)[-64:]]
    # Four-bar's public closure margins locate branch-near audit neighborhoods.
    branch_angles=[];branch_margins=None
    if hasattr(exact,'slider_states'):
        margins=[]
        for a in knots:
            states=exact.slider_states(float(a))
            margins.append(min(min(abs(s.four_bar_cross_product),s.connecting_rod_transverse_margin) for s in states))
        branch_angles=knots[np.argsort(margins)[:16]].tolist()
        branch_margins=dict(minimum_public_margin=float(min(margins)),note='Cross-product and rod margins have different units; used only to select additional audit angles.')
    audit=np.unique(np.r_[dense,boundary,roots,focus,branch_angles,
        [float(a)+eps for a in [*roots,*branch_angles,*focus] for eps in (-1e-7,1e-7)]])
    reference=np.array([exact.cylinder_volumes_and_derivatives(float(a)) for a in audit])
    scale=np.r_[np.ptp(reference[:,:2],axis=0),np.max(np.abs(reference[:,2:]),axis=0)]
    rows=[]
    for n in resolutions:
        start=time.perf_counter();spline=PeriodicVolumeInterpolant(exact,n);preparation=time.perf_counter()-start
        actual=np.array([spline.cylinder_volumes_and_derivatives(float(a)) for a in audit])
        error=np.abs(actual-reference)
        rows.append(dict(intervals=n,preparation_seconds=preparation,maximum_absolute_error=error.max(axis=0).tolist(),
            maximum_normalized_error=(error/scale).max(axis=0).tolist(),normalization=scale.tolist(),
            worst_angles=audit[np.argmax(error,axis=0)].tolist(),
            boundary_absolute_error=np.abs(np.array([spline.cylinder_volumes_and_derivatives(float(a)) for a in boundary])-
                np.array([exact.cylinder_volumes_and_derivatives(float(a)) for a in boundary])).max(axis=0).tolist(),
            periodic_closure_exact=spline.cylinder_volumes_and_derivatives(0.)==spline.cylinder_volumes_and_derivatives(TAU)))
    return dict(samples=len(audit),extremum_angles=roots,high_curvature_angles=focus.tolist(),
        branch_angles=branch_angles,branch_margins=branch_margins,
        note='Six-bar uses fixed validated branches; independent dense/high-curvature checks are not a global proof of clearance from singularity.',
        resolutions=rows)


def valve_behavior(wrapper,periodic):
    angles=periodic.angles;trajectory=periodic.trajectory
    pressures=np.array([ThermodynamicState.from_array(values[:8]).pressures(wrapper.model.gas,wrapper.model.volumes(float(a)))
        for a,values in zip(angles,trajectory.T)])
    result={}
    for name,u,d in (('hot_to_small',3,0),('cold_to_large',2,1)):
        delta=pressures[:,u]-pressures[:,d];opened=delta>0
        changes=np.flatnonzero(opened[1:]!=opened[:-1]);events=[]
        for i in changes:
            a,b=delta[i:i+2]
            events.append(float(angles[i]+(angles[i+1]-angles[i])*(-a)/(b-a)))
        result[name]=dict(open_time_fraction=float(np.trapz(opened.astype(float),angles)/(angles[-1]-angles[0])),
            transition_count=len(events),estimated_crossing_angles=events)
    return result


def run(output,resolutions):
    if output.exists():raise FileExistsError('Use a new interpolation output directory.')
    output.mkdir(parents=True)
    m=json.loads(DEFAULT_MANIFEST.read_text());_,base=_load_basis()
    raw=dict(m['numerical_settings']);raw['integration_absolute_tolerances']=tuple(raw['integration_absolute_tolerances'])
    definition=SimpleNamespace(wall_numerical_settings=WallCycleNumericalSettings(**raw),wall_backend=WallBackendSettings('numba'))
    warm=next(c for c in m['cases'] if c['name']=='production_warm')
    rhs=WallRHS(_build_design(base,warm['parameters']).build(),definition.wall_backend)
    start=time.perf_counter();rhs(0,np.r_[warm['initial_state'],np.zeros(5)])
    metadata=dict(environment=environment(),excluded_compilation_seconds=time.perf_counter()-start,
        source_hashes={str(p.relative_to(ROOT)):digest(p) for p in sorted((ROOT/'src').rglob('*.py'))},resolutions=resolutions)
    (output/'environment.json').write_text(json.dumps(metadata,indent=2)+'\n')
    extra=json.loads((ROOT/'outputs/solver_acceleration_stage5/manifest_sixbar.json').read_text())['cases'][0]
    for case in [next(c for c in m['cases'] if c['name']=='smooth_four_bar'),extra]:
        design=_build_design(base,case['parameters'])
        exact=base.kinematics if case['name']=='smooth_four_bar' else six_bar(design.configuration.machine_volumes)
        design=replace(design,kinematics=exact)
        study=geometry_study(exact,resolutions)
        (output/(case['name']+'_geometry.json')).write_text(json.dumps(study,indent=2)+'\n')
        reference=None;physical_rows=[]
        selected=[0]+[r['intervals'] for r in study['resolutions'] if max(r['maximum_normalized_error'])<1e-6]
        for n in selected:
            started=time.perf_counter()
            candidate=design if n==0 else replace(design,kinematics=PeriodicVolumeInterpolant(exact,n))
            preparation=time.perf_counter()-started;records=[];capture=[]
            try:
                result,state=_evaluate(case['name'],candidate,definition,initial_state=np.array(case['initial_state']),
                    statistics_callback=records.append,periodic_observer=capture.append)
                exception=None
            except (ValueError,RuntimeError,ArithmeticError) as error:
                result=dict(status=type(error).__name__,message=str(error));state=None;exception=type(error).__name__
            elapsed=time.perf_counter()-started
            performance=json_values(asdict(wall_cycle_performance(candidate.build(),capture[0].trajectory))) if capture and capture[0].converged else None
            row=dict(case=case['name'],intervals=n,result=result,final_state=None if state is None else state.tolist(),
                performance=performance,exception=exception,elapsed_seconds=elapsed,interpolation_preparation_seconds=preparation,
                solver_statistics=records,backend_statistics=capture[0].backend_statistics if capture else {})
            row['valve_behavior']=valve_behavior(candidate.build(),capture[0]) if capture and capture[0].converged else None
            if reference is None:reference=row
            else:
                # Convergence history is recorded, not an accuracy observable.
                compared={k:v for k,v in reference['result'].items() if k not in ('convergence','label')}
                actual={k:v for k,v in result.items() if k not in ('convergence','label')}
                differences=compare_fields(compared,actual)
                differences+=compare_fields(reference['final_state'],row['final_state'],'final_state')
                # Conservation residuals are near-zero solver roundoff, assessed separately.
                if reference['performance'] and performance:
                    differences+=compare_fields({k:v for k,v in reference['performance'].items() if k!='conservation'},
                        {k:v for k,v in performance.items() if k!='conservation'},'performance')
                if reference['performance'] and performance:
                    differences+=compare_fields({k:v for k,v in reference['performance']['conservation'].items() if k.startswith('relative_')},
                        {k:v for k,v in performance['conservation'].items() if k.startswith('relative_')},'conservation')
                differences+=compare_fields(reference['valve_behavior'],row['valve_behavior'],'valve_behavior')
                row['comparison']=dict(passed=not differences,failed_fields=differences,
                    reference_cycles=reference['result'].get('convergence',{}).get('cycles_completed'),
                    actual_cycles=result.get('convergence',{}).get('cycles_completed'))
            if capture and capture[0].trajectory is not None:
                np.savez_compressed(output/f'{case["name"]}_{n}.npz',angles=capture[0].angles,trajectory=capture[0].trajectory)
            physical_rows.append(row)
            (output/f'{case["name"]}_{n}.json').write_text(json.dumps(row,indent=2)+'\n')
            print(case['name'],n,result['status'],f'{elapsed:.3f}s','gate',row.get('comparison',{}).get('passed'),flush=True)
        (output/(case['name']+'_summary.json')).write_text(json.dumps([dict(intervals=r['intervals'],seconds=r['elapsed_seconds'],
            passed=r.get('comparison',{}).get('passed'),failed_fields=len(r.get('comparison',{}).get('failed_fields',[]))) for r in physical_rows],indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--resolutions',type=int,nargs='+',default=[128,256,512,1024,2048])
    a=p.parse_args();run(a.output,a.resolutions)
