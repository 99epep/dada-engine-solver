"""Check the compiled prototype on frozen real states and their extrema."""
from dataclasses import replace
import json
import numpy as np
from benchmark_solver_acceleration import DEFAULT_MANIFEST,ROOT
from refine_motor_four_stage_hx9d_variable_gas import _load_basis,_build_design
from numba_wall_prototype import PreparedRHS
from dada_solver.exchangers.air_wall import AirWallMotor


def check():
    manifest=json.loads(DEFAULT_MANIFEST.read_text());_,base=_load_basis();results={}
    for case in manifest['cases']:
        if case['name'] not in ('production_warm','production_cold','nearby_warm','smooth_four_bar'): continue
        design=_build_design(base,case['parameters'])
        if case.get('motion')=='base': design=replace(design,kinematics=base.kinematics)
        wrapper=design.build();prepared=PreparedRHS(wrapper,AirWallMotor.derivative)
        with np.load(ROOT/'outputs/solver_acceleration_stage2/p0_baseline'/(case['name']+'_0.npz')) as saved:
            angles=saved['angles'];trajectory=saved['trajectory']
        selected=set(np.linspace(0,len(angles)-1,256,dtype=int))
        for row in trajectory[:10]:
            selected.update((int(np.argmin(row)),int(np.argmax(row))))
        for edge in getattr(wrapper.model.kinematics,'breakpoint_angles',lambda:())():
            i=int(np.searchsorted(angles,edge));selected.update(j for j in (i-1,i,i+1) if 0<=j<len(angles))
        maximum=0.
        for i in sorted(selected):
            expected=wrapper.derivative(float(angles[i]),trajectory[:,i])
            actual=prepared(float(angles[i]),trajectory[:,i])
            scale=1e-11+2e-11*np.abs(expected)
            maximum=max(maximum,float(np.max(np.abs(actual-expected)/scale)))
            np.testing.assert_allclose(actual,expected,rtol=2e-11,atol=1e-11)
        results[case['name']]=dict(points=len(selected),fallbacks=prepared.fallbacks,
            maximum_scaled_rhs_difference=maximum,first_call_seconds=prepared.first_call_seconds)
    print(json.dumps(results,indent=2))


if __name__=='__main__': check()
