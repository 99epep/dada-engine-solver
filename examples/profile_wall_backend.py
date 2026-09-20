"""Bounded native/dispatcher microbenchmarks supplement full section timings.

The native loop removes Python crossings but includes native call/allocation and
loop cost. Samples are final-cycle states, not the full cold transient. Differences
are estimates, not additive exact attribution of an instrumented evaluation.
"""
from dataclasses import replace
import argparse
import json
from pathlib import Path
import time
import numpy as np
from numba import njit
from benchmark_solver_acceleration import ROOT,DEFAULT_MANIFEST
from refine_motor_four_stage_hx9d_variable_gas import _load_basis,_build_design
from dada_solver.wall_backend import WallRHS,WallBackendSettings


@njit
def native_batch(function,states,volumes,rates,gas,links,walls,sources,destinations,one_way,ports,passes):
    checksum=0.
    for _ in range(passes):
        for i in range(len(states)):
            ok,result=function(states[i],volumes[i],rates[i],gas,links,walls,sources,destinations,one_way,ports)
            checksum+=result[0]
    return checksum


def measure(output,repeats=5,passes=64):
    if output.exists(): raise FileExistsError('Use a new profile artifact.')
    manifest=json.loads(DEFAULT_MANIFEST.read_text());_,base=_load_basis();results={}
    for case in manifest['cases']:
        if case['name'] not in ('production_warm','production_cold','smooth_four_bar'): continue
        d=_build_design(base,case['parameters'])
        if case.get('motion')=='base':d=replace(d,kinematics=base.kinematics)
        w=d.build();rhs=WallRHS(w,WallBackendSettings('numba'));p=rhs.implementation
        with np.load(ROOT/'outputs/solver_acceleration_stage2/p0_baseline'/(case['name']+'_0.npz')) as saved:
            ids=np.linspace(0,len(saved['angles'])-1,128,dtype=int)
            angles=saved['angles'][ids];states=np.ascontiguousarray(saved['trajectory'][:,ids].T)
        volumes=[];rates=[]
        for a in angles:
            point=w.model.volumes(float(a));volumes.append(point.as_tuple());rates.append(w.model.cylinder_volume_rates(float(a)))
        volumes=np.array(volumes);rates=np.array(rates)
        args=(p.gas,p.links,p.walls,p.sources,p.destinations,p.one_way,p.ports)
        rhs(float(angles[0]),states[0])
        native_batch(p.dispatcher,states,volumes,rates,*args,1)
        native=[];dispatch=[];kinematics=[];checksum=None
        for _ in range(repeats):
            start=time.perf_counter();expected=native_batch(p.dispatcher,states,volumes,rates,*args,passes)
            native.append(time.perf_counter()-start)
            start=time.perf_counter();checksum=0.
            for _ in range(passes):
                for i in range(len(states)):
                    ok,result=p.dispatcher(states[i],volumes[i],rates[i],*args);checksum+=result[0]
            dispatch.append(time.perf_counter()-start)
            np.testing.assert_allclose(checksum,expected,rtol=1e-12,atol=1e-12)
            provider=getattr(w.model.kinematics,'cylinder_volumes_and_derivatives',None)
            start=time.perf_counter()
            for _ in range(passes):
                for a in angles:
                    if provider is not None:provider(float(a))
                    else:
                        w.model.volumes(float(a));w.model.cylinder_volume_rates(float(a))
            kinematics.append(time.perf_counter()-start)
        calls=len(states)*passes
        results[case['name']]=dict(calls_per_repeat=calls,repeats=repeats,
            native_batch_seconds=native,python_dispatch_batch_seconds=dispatch,kinematics_batch_seconds=kinematics,
            native_seconds_per_call=float(np.median(native))/calls,
            dispatcher_seconds_per_call=float(np.median(dispatch))/calls,
            kinematics_seconds_per_call=float(np.median(kinematics))/calls)
    output.write_text(json.dumps(dict(cases=results,limitations=__doc__),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path)
    a=p.parse_args();measure(a.output)
