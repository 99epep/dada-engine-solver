"""Check exact offset-slider derivatives and inversion used in the comparison."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

spec=importlib.util.spec_from_file_location('slider_comparison',Path(__file__).resolve().parents[1]/'examples/compare_slider_crank_motion.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

@pytest.mark.parametrize('offset',[0.,.7,-.7])
def test_offset_slider_derivative_wrap_and_inversion(offset):
    theta=np.linspace(0,2*np.pi,100,endpoint=False)
    q,dq=module.motion(theta,3.,offset,.3)
    h=1e-5
    qp,_=module.motion(theta+h,3.,offset,.3)
    qm,_=module.motion(theta-h,3.,offset,.3)
    assert np.allclose(dq,(qp-qm)/(2*h),atol=1e-9)
    inv,di=module.motion(theta,3.,offset,.3,True)
    assert np.allclose(inv,1-q)
    assert np.allclose(di,-dq)
    wrap,dwrap=module.motion(theta+2*np.pi,3.,offset,.3)
    assert np.allclose(wrap,q)
    assert np.allclose(dwrap,dq)

def test_invalid_full_rotation_is_rejected():
    with pytest.raises(ValueError):module.motion(np.array([0.]),2.,1.,0.)
