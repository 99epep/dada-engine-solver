"""Supplementary fixed-state RHS timing; never substitutes for cycle benchmarks."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from benchmark_solver_acceleration import DEFAULT_MANIFEST, ROOT, environment
from refine_motor_four_stage_hx9d_variable_gas import _load_basis, _build_design


def run(output,repeats=7,passes=8):
    if output.exists(): raise FileExistsError('Use a new RHS timing output.')
    manifest=json.loads(DEFAULT_MANIFEST.read_text())
    case=next(x for x in manifest['cases'] if x['name']=='production_warm')
    _,base=_load_basis();wrapper=_build_design(base,case['parameters']).build()
    with np.load(ROOT/'outputs/solver_acceleration_stage2/p0_baseline/production_warm_0.npz') as saved:
        indices=np.linspace(0,len(saved['angles'])-1,256,dtype=int)
        points=[(float(saved['angles'][i]),saved['trajectory'][:,i].copy()) for i in indices]
    for angle,state in points: wrapper.derivative(angle,state)
    timings=[]
    for repeat in range(repeats):
        start=time.perf_counter()
        for _ in range(passes):
            for angle,state in points: wrapper.derivative(angle,state)
        timings.append(time.perf_counter()-start)
    result=dict(environment=environment(),rhs_calls_per_repeat=passes*len(points),
        elapsed_seconds=timings,median_seconds=float(np.median(timings)),maximum_seconds=max(timings))
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path);args=parser.parse_args();run(args.output)
