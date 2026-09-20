"""Numerical and fallback gates for the optional research-only compiled RHS."""
from dataclasses import replace
from pathlib import Path
import sys
import math
import numpy as np
import pytest

pytest.importorskip('numba')
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'examples'))
from numba_wall_prototype import PreparedRHS, kernel, properties, prototype_backend
from dada_solver.exchangers.air_wall import AirWallMotor
from dada_solver.exchangers.gas_transport import DiluteGasTransport
from dada_solver.state import UniformCharge
from tests.test_solver_acceleration import variable_wrapper


def state(wrapper,angle=0.,temperature=350.,pressure=2e5):
    gas=UniformCharge(pressure,temperature).create_state(wrapper.model.gas,wrapper.model.volumes(angle))
    return np.r_[gas.as_array(),wrapper.heat_in.wall_capacity_j_k*450,
        wrapper.heat_out.wall_capacity_j_k*330,np.zeros(5)]


@pytest.mark.parametrize('temperature',[200.,298.15,499.999,500.,598.15,699.999,700.,1000.])
def test_transport_piecewise_boundaries(temperature):
    tr=DiluteGasTransport()
    np.testing.assert_allclose(properties(temperature),
        (tr.viscosity(temperature),tr.conductivity(temperature),tr.cp(temperature)),rtol=5e-15,atol=0.)


@pytest.mark.parametrize('delta',[-1e-5,-1e-12,0.,1e-12,1e-5])
def test_compiled_rhs_valves_reflux_breakpoints_and_conservation(delta):
    w=variable_wrapper();p=PreparedRHS(w,AirWallMotor.derivative)
    for edge in (0.,*w.model.kinematics.breakpoint_angles(),2*math.pi):
        for a in (edge-1e-10,edge,edge+1e-10):
            x=state(w,a);x[1]*=1+delta;x[7]*=1-delta
            actual=p(a,x);expected=w.derivative(a,x)
            np.testing.assert_allclose(actual,expected,rtol=2e-11,atol=1e-11)
            assert abs(sum(actual[:8:2]))<1e-15
            # Total gas + wall storage = external heat - work, per radian.
            assert abs(sum(actual[1:8:2])+sum(actual[8:10])-sum(actual[10:12])+actual[14])<1e-9
    assert p.fallbacks==0
    assert kernel.nopython_signatures
    assert not kernel.targetoptions.get('fastmath',False)


@pytest.mark.parametrize('temperature,pressure',[(2000.,2e5),(350.,100.),(350.,2e5)])
def test_invalid_states_keep_reference_error(temperature,pressure):
    w=variable_wrapper();p=PreparedRHS(w,AirWallMotor.derivative)
    x=state(w,temperature=temperature,pressure=pressure)
    if pressure==2e5 and temperature==350: x[0]=-1.
    with pytest.raises(ValueError) as expected: w.derivative(0.,x)
    with pytest.raises(type(expected.value)) as actual: p(0.,x)
    assert str(actual.value)==str(expected.value)
    assert p.fallbacks==1


def test_candidate_constants_and_unsupported_family_fallback():
    w=variable_wrapper();x=state(w)
    p=PreparedRHS(w,AirWallMotor.derivative)
    altered=replace(w,heat_in=replace(w.heat_in,air_mass_flow_kg_s=2*w.heat_in.air_mass_flow_kg_s))
    other=PreparedRHS(altered,AirWallMotor.derivative)
    assert not np.array_equal(p.walls,other.walls)
    assert not p.walls.flags.writeable
    legacy=replace(w,heat_in=replace(w.heat_in,gas_film=None))
    with pytest.raises(TypeError): PreparedRHS(legacy,AirWallMotor.derivative)
    expected=legacy.derivative(0.,x)
    with prototype_backend() as stats:
        np.testing.assert_array_equal(legacy.derivative(0.,x),expected)
    assert stats['unsupported_calls']==1


def test_domain_fallback_calls_authoritative_reference_once():
    w=variable_wrapper();calls=[]
    def reference(wrapper,angle,values):
        calls.append((wrapper,angle));return np.arange(15.)
    p=PreparedRHS(w,reference)
    x=state(w,pressure=100.)
    np.testing.assert_array_equal(p(0.,x),np.arange(15.))
    assert len(calls)==1 and p.fallbacks==1


@pytest.mark.parametrize('ratio',[1.+1e-12,1.01,1.2,2.,10.])
def test_large_pressure_gradients_preserve_result_or_failure(ratio):
    w=variable_wrapper();p=PreparedRHS(w,AirWallMotor.derivative)
    x=state(w)
    # Scale mass and energy together: change pressure without changing T.
    x[:2]*=ratio
    try:
        expected=w.derivative(0.,x)
    except ValueError as error:
        with pytest.raises(type(error)) as actual: p(0.,x)
        assert str(actual.value)==str(error)
        assert p.fallbacks==1
    else:
        np.testing.assert_allclose(p(0.,x),expected,rtol=2e-11,atol=1e-11)


def test_short_array_delegates_instead_of_entering_unchecked_jit_memory():
    w=variable_wrapper();calls=[]
    def reference(wrapper,angle,values):
        calls.append(values.shape)
        raise ValueError('Reference shape validation')
    p=PreparedRHS(w,reference)
    with pytest.raises(ValueError,match='Reference shape validation'): p(0.,np.ones(3))
    assert calls==[(3,)]


def test_backend_scope_restores_reference_after_failure():
    original=AirWallMotor.derivative
    with pytest.raises(RuntimeError,match='Scope failure'):
        with prototype_backend():
            raise RuntimeError('Scope failure')
    assert AirWallMotor.derivative is original
