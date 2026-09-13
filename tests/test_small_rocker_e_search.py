"""Check the independent E parameterization closes a finite slider mechanism."""
import importlib.util
from pathlib import Path
import sys
import numpy as np
import pytest

EXAMPLES=Path(__file__).resolve().parents[1]/'examples'
sys.path.insert(0,str(EXAMPLES))
try:
    from search_small_rocker_e import build, assess
finally:
    sys.path.pop(0)


def test_independent_e_finite_rod_derivatives():
    v=[2.5,2.5,2.925,2.925,-.97,10.,0.,0.,.3]
    assembly,phase=build(v,1)
    t=.7;h=1e-6
    def state(angle):
        a=angle+phase
        return assembly.evaluate(np.array([np.cos(a),np.sin(a)]),np.array([-np.sin(a),np.cos(a)]))
    derivative=(state(t+h).coordinate-state(t-h).coordinate)/(2*h)
    assert state(t).coordinate_derivative == pytest.approx(derivative,abs=1e-8)


def test_full_revolution_closure_is_checked():
    theta=np.linspace(0,2*np.pi,20)
    with pytest.raises(ValueError,match='throughout'):
        assess([5,1.2,1.2,1,0,10,0,0,0],1,1,theta,theta,theta)


def test_f_output_uses_coupler_frame_with_finite_rod():
    from dada_solver.four_bar import CouplerOutputPoint
    v=[2.5,2.5,2.925,1.2,.3,10.,0.,0.,.3]
    assembly,phase=build(v,1,'F')
    assert isinstance(assembly.output,CouplerOutputPoint)
    assert assembly.slider.connecting_rod_length == 10.
    def state(t):
        a=t+phase
        return assembly.evaluate(np.array([np.cos(a),np.sin(a)]),np.array([-np.sin(a),np.cos(a)]))
    t=.7;h=1e-6
    numerical=(state(t+h).coordinate-state(t-h).coordinate)/(2*h)
    assert state(t).coordinate_derivative == pytest.approx(numerical,abs=1e-8)
