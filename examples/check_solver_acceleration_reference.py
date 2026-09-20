"""Check saved real trajectory samples against the frozen test-only RHS oracle."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np

from benchmark_solver_acceleration import DEFAULT_MANIFEST, ROOT, legacy_design
from refine_motor_four_stage_hx9d_variable_gas import _load_basis, _build_design
from dada_solver.factory import initial_valve_topology
from dada_solver.state import ThermodynamicState
sys.path.insert(0,str(ROOT))
from tests import solver_acceleration_reference as reference
from tests.test_solver_acceleration import assert_rates_equal


def check(directory):
    manifest=json.loads(DEFAULT_MANIFEST.read_text());_,base=_load_basis()
    results={}
    for case in manifest['cases']:
        path=directory/(case['name']+'_0.npz')
        if not path.exists(): continue
        design=legacy_design() if case.get('legacy') else _build_design(base,case['parameters'])
        if case.get('motion')=='base': design=replace(design,kinematics=base.kinematics)
        wrapper=design.build()
        with np.load(path) as saved:
            angles=saved['angles'];trajectory=saved['trajectory']
        selected=set(np.linspace(0,len(angles)-1,64,dtype=int))
        temperatures=trajectory[1:8:2]/(trajectory[:8:2]*wrapper.model.gas.heat_capacity_cv)
        for row in temperatures:
            selected.update((int(np.argmin(row)),int(np.argmax(row))))
        for edge in getattr(wrapper.model.kinematics,'breakpoint_angles',lambda:())():
            index=np.searchsorted(angles,edge)
            selected.update(i for i in (index-1,index,index+1) if 0<=i<len(angles))
        for index in sorted(selected):
            angle=float(angles[index]);values=trajectory[:,index]
            gas=ThermodynamicState.from_array(values[:8]);topology=initial_valve_topology()
            expected=reference.evaluate(wrapper.model,angle,gas,topology)
            assert_rates_equal(wrapper.model.evaluate(angle,gas,topology),expected)
            assert wrapper.flow_contexts(angle,values)==reference.flow_contexts(wrapper,angle,values)
            assert wrapper.thermal_rates(angle,values)==reference.thermal_rates(wrapper,angle,values)
            np.testing.assert_array_equal(wrapper.derivative(angle,values),reference.derivative(wrapper,angle,values))
        results[case['name']]=len(selected)
    print(json.dumps(dict(exact_saved_state_checks=results),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline_directory',type=Path)
    args=parser.parse_args();check(args.baseline_directory)
