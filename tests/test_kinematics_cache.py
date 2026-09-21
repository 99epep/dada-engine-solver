"""Exact provider reuse, isolation and interruption equivalence."""
import json
import sys
from pathlib import Path
import numpy as np
import pytest
from dada_solver.kinematics_cache import prepare_exact_kinematics
from dada_solver.exchangers.wall_cycle import solve_periodic_wall_motor
from dada_solver.integration import IntegrationInterrupted
from dada_solver.wall_backend import WallBackendSettings


def test_exact_adjacent_angles_and_failed_call_do_not_pollute_cache():
    class Provider:
        calls=0
        def cylinder_volumes_and_derivatives(self,a):
            self.calls+=1
            if a<0: raise ValueError('Invalid test angle.')
            return [a,a*2,1.,2.]
    p=Provider();a=prepare_exact_kinematics(p);b=prepare_exact_kinematics(p)
    x=.5;y=np.nextafter(x,1.)
    for angle in [x,x,y,y,x]: a.cylinder_volumes_and_derivatives(angle)
    assert a.snapshot()['hits']==2 and p.calls==3
    with pytest.raises(ValueError): a.cylinder_volumes_and_derivatives(-1.)
    assert a.cylinder_volumes_and_derivatives(x)==(x,2*x,1.,2.)
    b.cylinder_volumes_and_derivatives(x)
    assert b.snapshot()['hits']==0 and p.calls==5
    assert prepare_exact_kinematics(a).snapshot()['calls']==0


def test_no_combined_provider_is_not_wrapped():
    p=object()
    assert prepare_exact_kinematics(p) is p


@pytest.mark.parametrize('backend',['python','numba'])
def test_cycle_and_retained_endpoint_after_interruption_are_exact(backend):
    if backend=='numba': pytest.importorskip('numba')
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'examples'))
    from benchmark_solver_acceleration import DEFAULT_MANIFEST
    from refine_motor_four_stage_hx9d_variable_gas import _load_basis,_build_design
    case=next(c for c in json.loads(DEFAULT_MANIFEST.read_text())['cases'] if c['name']=='production_cold')
    _,base=_load_basis();w=_build_design(base,case['parameters']).build()
    initial=np.asarray(case['initial_state'])
    results=[]
    for enabled in (False,True):
        completed=[]
        def progress(info):
            if completed: raise IntegrationInterrupted('Test deadline after one cycle.')
        result=solve_periodic_wall_motor(w,initial,maximum_cycles=3,
            progress_callback=progress,cycle_callback=lambda *_:completed.append(True),
            backend=WallBackendSettings(backend),exact_kinematics_cache=enabled)
        results.append(result)
    a,b=results
    assert a.status==b.status=='interrupted' and a.message==b.message
    assert a.history==b.history
    np.testing.assert_array_equal(a.angles,b.angles)
    np.testing.assert_array_equal(a.trajectory,b.trajectory)
    np.testing.assert_array_equal(a.last_complete_state,b.last_complete_state)
    assert b.backend_statistics['exact_kinematics_cache']['hits']>0


def test_runtime_protocol_and_scalar_methods_remain_available():
    from dada_solver.kinematics import KinematicsModel
    from tests.test_solver_acceleration import variable_wrapper
    original=variable_wrapper().model.kinematics
    prepared=prepare_exact_kinematics(original)
    assert isinstance(prepared,KinematicsModel)
    assert prepared.small_volume_limits==original.small_volume_limits
    assert prepared.large_volume_limits==original.large_volume_limits
    for method in ('small_cylinder_volume','large_cylinder_volume',
                   'small_cylinder_volume_derivative','large_cylinder_volume_derivative'):
        for angle in (0.,.314,-.314):
            assert getattr(prepared,method)(angle)==getattr(original,method)(angle)
    assert prepared.breakpoint_angles()==original.breakpoint_angles()
