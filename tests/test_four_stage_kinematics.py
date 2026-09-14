"""Four-stage user chronology, angle adaptation and exact linear derivatives."""
import math
from pathlib import Path
import numpy as np
import pytest
from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_model
from dada_solver.four_stage_kinematics import FourStageVolumeKinematics


@pytest.fixture
def model():
    config=load_simulation_configuration(Path(__file__).resolve().parents[1]/'examples/motor_demonstrator_original_325c.toml')
    limits=config.machine_volumes
    k=FourStageVolumeKinematics(limits.small_cylinder,limits.large_cylinder,.3,.55,.8,.4,.7,.6,.2)
    return build_model(config,kinematics=k)


def test_forward_motor_chronology_and_breakpoints(model):
    k=model.kinematics
    for t,qs,ql in [(0,.2,1),(.3,1,.7),(.55,.6,0),(.8,0,.4),(1,.2,1)]:
        assert k.small_cylinder_volume(t*2*math.pi)==pytest.approx(k.small_volume_limits.minimum+qs*k.small_volume_limits.swept)
        assert k.large_cylinder_volume(t*2*math.pi)==pytest.approx(k.large_volume_limits.minimum+ql*k.large_volume_limits.swept)
    np.testing.assert_allclose(k.breakpoint_angles(),2*math.pi*np.array([.3,.55,.8]))


def test_linear_slopes_wrap_and_bounded_volumes(model):
    k=model.kinematics
    for side in ['small','large']:
        v=getattr(k,f'{side}_cylinder_volume');dv=getattr(k,f'{side}_cylinder_volume_derivative')
        limits=getattr(k,f'{side}_volume_limits')
        for t in [.1,.4,.6,.9]:
            a=t*2*math.pi;h=1e-6
            assert dv(a)==pytest.approx((v(a+h)-v(a-h))/(2*h),rel=1e-8)
            assert v(a)==pytest.approx(v(a+2*math.pi),abs=1e-15)
        values=np.array([v(a) for a in np.linspace(0,2*math.pi,1001)])
        assert values.min()>=limits.minimum-1e-15
        assert values.max()<=limits.maximum+1e-15
        for t in [0,.3,.55,.8,1]:
            assert v(t*2*math.pi-1e-10)==pytest.approx(v(t*2*math.pi+1e-10),abs=1e-12)


@pytest.mark.parametrize('updates',[{'t1':0},{'t2':.1},{'t3':1},{'a_s':-1},{'b_l':float('nan')}])
def test_rejects_invalid_parameterization(model,updates):
    from dataclasses import replace
    with pytest.raises(ValueError): replace(model.kinematics.forward,**updates)
