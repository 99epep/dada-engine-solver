"""Bounded same-state RHS benchmark; no optimization or helium accuracy claim.

PYTHONPATH=src python3 examples/benchmark_research_v3_backends.py --output result.json
"""
import argparse
from dataclasses import replace
from importlib.resources import files
import json
from pathlib import Path
import platform
import time
import numpy as np
from dada_solver.research.schema import load_study, compile_study
from dada_solver.research.sixbar import PARAMETERS
from dada_solver.wall_backend import WallRHS, WallBackendSettings


def representative():
    d=compile_study(load_study(files('dada_solver.research').joinpath('data/sixbar-thermo5d.toml')))
    p={name:d.study.basis.data['reference_parameters'][spec[2]] for name,spec in PARAMETERS.items()}
    w=d.adapter.build(dict(d.fixed_parameters,**p)).build()
    x=np.array(d.initial_wall_state['values'])
    x[8:]*=np.array([w.heat_in.wall_capacity_j_k,w.heat_out.wall_capacity_j_k])/d.initial_wall_state['wall_capacities_j_k']
    return w,x


def benchmark(calls=10000):
    from dada_solver.tabulated_fluid import ideal_validation_table
    if calls<1: raise ValueError('Positive call count required.')
    w,x=representative()
    table=ideal_validation_table(w.model.gas,rho_axis=(.01,.1,.5,1.,2.,5.,10.,20.,50.))
    tabulated=replace(w,model=replace(w.model,gas=table))
    measurements={};values={}
    for name,wrapper,backend in (('ideal_compiled',w,'numba'),('tabulated_compiled',tabulated,'numba'),('python',w,'python')):
        r=WallRHS(wrapper,WallBackendSettings(backend))
        start=time.perf_counter();values[name]=r(0.,x);first=time.perf_counter()-start
        count=calls if backend=='numba' else min(calls,1000)
        durations=[]
        for repeat in range(3):
            start=time.perf_counter()
            for _ in range(count): r(0.,x)
            durations.append((time.perf_counter()-start)/count)
        steady=float(np.median(durations))
        measurements[name]=dict(first_call_seconds=first,calls_per_repeat=count,repeats=3,
            steady_seconds_per_call=steady,calls_per_second=1/steady,backend=r.snapshot())
    diff=values['tabulated_compiled']-values['ideal_compiled']
    return dict(scope='Same rank01 warm state and angle, one local machine; JIT and steady time separated. No helium accuracy claim.',
        platform=platform.platform(),processor=platform.processor(),fluid=table.identity,
        measurements=measurements,tabulated_to_ideal_ratio=measurements['tabulated_compiled']['steady_seconds_per_call']/measurements['ideal_compiled']['steady_seconds_per_call'],
        rhs_maximum_absolute_error=float(np.max(np.abs(diff))),
        rhs_maximum_scaled_error=float(np.max(np.abs(diff)/np.maximum(1.,abs(values['ideal_compiled'])))))


def benchmark_ideal_only(calls=30000):
    """Also runnable against exported pre-V3 sources through PYTHONPATH."""
    if calls<1: raise ValueError('Positive call count required.')
    w,x=representative();r=WallRHS(w,WallBackendSettings('numba'))
    start=time.perf_counter();r(0.,x);first=time.perf_counter()-start
    samples=[]
    for _ in range(5):
        start=time.perf_counter()
        for i in range(calls): r(0.,x)
        samples.append((time.perf_counter()-start)/calls)
    return dict(first_call_seconds=first,steady_seconds_per_call=float(np.median(samples)),
        repeats=5,calls_per_repeat=calls,samples=samples,backend=r.snapshot())


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--calls',type=int,default=10000)
    parser.add_argument('--ideal-only',action='store_true',help='Benchmark only the ideal path, including with an exported pre-V3 PYTHONPATH.')
    args=parser.parse_args();result=(benchmark_ideal_only if args.ideal_only else benchmark)(args.calls)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='measurements'},indent=2))
