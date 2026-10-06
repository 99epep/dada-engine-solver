"""Explicit transition closure, unchanged guards, and compiled backend parity."""
from dataclasses import replace
import math
from types import SimpleNamespace

import numpy as np
import pytest

from dada_solver import numerical_primitives as numeric
from dada_solver.exchangers.gas_correlations import (
    MicrotubeGasModel, MicrotubeDomainError, transition_nusselt, transition_darcy,
    laminar_entry_nusselt, gnielinski, darcy_smooth)
from dada_solver.exchangers.hardware import TubeHalfLink
from dada_solver.exchangers.microtube_geometry import MicrotubeBank
from dada_solver.fluids import CaloricallyPerfectGas

BANK = MicrotubeBank(1000, .8, .000770, .0001524, .00125, .001)
MODEL = MicrotubeGasModel()


def diagnostic(re, **kwargs):
    flow = re*BANK.tube_flow_area_m2*MODEL.transport.viscosity(350)/BANK.inner_diameter_m
    return MODEL.diagnose(BANK, flow, kwargs.pop('p1', 1e6), kwargs.pop('p2', 1e6), 350, **kwargs)


@pytest.mark.parametrize('re,regime', [(2299.999,'laminar'),(2300,'transition'),
    (2500,'transition'),(3000,'transition'),(3999.999,'transition'),(4000,'turbulent'),(4500,'turbulent')])
def test_regimes_and_domain(re, regime):
    d = diagnostic(re)
    assert d.flow_regime == regime
    assert 'reynolds_outside_correlation_domain' not in d.issues
    MODEL.require(d)
    assert d.transition_model_used == (regime == 'transition')
    if regime == 'transition': assert d.correlation_id == 'bennett_gnielinski_transition_interpolation'


def test_thermal_endpoints():
    pr = .7
    low = numeric.bennett_mean_nusselt(2300,pr,BANK.inner_diameter_m/BANK.tube_length_m)
    assert transition_nusselt(2300,pr,BANK.inner_diameter_m/BANK.tube_length_m) == low
    assert transition_nusselt(4000,pr,BANK.inner_diameter_m/BANK.tube_length_m) == pytest.approx(gnielinski(4000,pr))
    for edge in (2300,4000):
        assert diagnostic(edge-1e-5).nusselt == pytest.approx(diagnostic(edge+1e-5).nusselt,rel=1e-7)
    assert transition_nusselt(2300,pr,.01,False) == 3.66


@pytest.mark.parametrize('re', [2299.999,2300,2500,3000,3999.999,4000,4500])
def test_hydraulic_flow_and_friction_continuity(re):
    gas = CaloricallyPerfectGas(287.05,1005.,717.95)
    mu = MODEL.transport.viscosity(350)
    flow = re*BANK.tube_flow_area_m2*mu/BANK.inner_diameter_m
    f = 64/re if re<2300 else transition_darcy(re) if re<4000 else darcy_smooth(re)
    # Include entrance loss and its continuous transition endpoint.
    if re<2300:
        f += numeric.shah_entry_excess(BANK.tube_length_m/(2*BANK.inner_diameter_m*re))*BANK.inner_diameter_m/(BANK.tube_length_m/2)
    elif re<4000:
        extra=numeric.shah_entry_excess(BANK.tube_length_m/(2*BANK.inner_diameter_m*2300))*BANK.inner_diameter_m/(BANK.tube_length_m/2)
        f=numeric.transition_friction(re,1.,extra)
    pout = 1e6
    c = f*(BANK.tube_length_m/2)/BANK.inner_diameter_m * flow**2 * gas.gas_constant*350/BANK.tube_flow_area_m2**2
    pin = math.sqrt(pout*pout+c)
    link = TubeHalfLink(BANK,mu,1.,0.,gas_model=MODEL)
    actual = link.directed_flow(pin,pout,350,gas).mass_flow_rate
    assert actual == pytest.approx(flow,rel=2e-10)
    if re in (2300,4000):
        left = link.directed_flow(pin-1e-5,pout,350,gas).mass_flow_rate
        right = link.directed_flow(pin+1e-5,pout,350,gas).mass_flow_rate
        assert left == pytest.approx(right,rel=1e-7)


def test_other_guards_remain():
    assert 'large_relative_pressure_drop' in diagnostic(3000,p1=1.5e6).issues
    assert 'high_mach' in diagnostic(3000,p1=1e4,p2=1e4).issues
    assert 'thermal_slip_not_implemented' in diagnostic(3000,p1=100,p2=100).issues
    assert 'beyond_continuum_model' in diagnostic(3000,p1=1,p2=1).issues
    assert 'reynolds_outside_correlation_domain' in diagnostic(5e6+1).issues
    with pytest.raises(ValueError): MODEL.transport.viscosity(99.9)
    with pytest.raises(MicrotubeDomainError): MODEL.require(diagnostic(3000,p1=1.5e6))
    assert numeric.thermal_kind(3000,.7) == 3
    assert not numeric.laminar_diagnostics(3000*BANK.tube_flow_area_m2*MODEL.transport.viscosity(350)/BANK.inner_diameter_m,1e6,1e6,350,.00077,.8,
        BANK.tube_flow_area_m2,.3,.2,True,0)[0]  # The laminar-only diagnostic rejects transition states.


def test_cycle_transition_weights():
    from dada_solver.exchangers.gas_diagnostics import cycle_microtube_diagnostics
    rows = [diagnostic(re) for re in (1000,3000,4500)]
    walls = [SimpleNamespace(film_diagnostics=(d,d),gas_heat_w=q) for d,q in zip(rows,(1.,2.,3.))]
    replay = SimpleNamespace(require=lambda **kwargs: None, samples=[
        SimpleNamespace(walls=(w,w),hydraulic_diagnostics=(d,d,d,d)) for w,d in zip(walls,rows)])
    wrapper = SimpleNamespace(heat_in=SimpleNamespace(requires_flow_context=True),heat_out=None,
                              model=SimpleNamespace(angular_speed=1.))
    report = cycle_microtube_diagnostics(wrapper,np.array([0.,1.,2.]),np.zeros((10,3)),replay=replay)
    assert report['model_validity'] == 'valid'
    assert report['transition_model_used']
    assert report['transition_time_fraction'] == .5
    assert report['transition_absolute_heat_fraction'] == .5
    p = report['passages']['Hi.inlet']
    assert p['transition_reynolds_range']['minimum'] == pytest.approx(3000)
    assert p['domains']['flow_regime']['laminar']['time_fraction'] == .25
    assert p['domains']['flow_regime']['turbulent']['time_fraction'] == .25


def test_transition_rhs_python_numba_without_fallback():
    pytest.importorskip('numba')
    from tests.test_solver_acceleration import variable_wrapper
    from tests.test_wall_backend import values
    from dada_solver.wall_backend import WallRHS, WallBackendSettings
    w = variable_wrapper()
    def link(x): return replace(x,bank=BANK)
    model = replace(w.model,small_cold_link=link(w.model.small_cold_link),
        large_hot_link=link(w.model.large_hot_link),
        cold_large_valve=replace(w.model.cold_large_valve,flow_model=link(w.model.cold_large_valve.flow_model)),
        hot_small_valve=replace(w.model.hot_small_valve,flow_model=link(w.model.hot_small_valve.flow_model)))
    def wall(x): return replace(x,gas_film=replace(x.gas_film,bank=BANK))
    w = replace(w,model=model,heat_in=wall(w.heat_in),heat_out=wall(w.heat_out))
    rhs = WallRHS(w,WallBackendSettings('numba'))
    for delta in np.linspace(.001,.2,60):
        x = values(w,pressure=1e6)
        x[:2] *= 1+delta
        contexts = w.flow_contexts(0.,x)
        _, diagnostics = w.heat_in.gas_film.evaluate(350.,350.,contexts[0])
        if any(d.transition_model_used for d in diagnostics):
            expected = w.derivative(0.,x)
            np.testing.assert_allclose(rhs(0.,x),expected,rtol=2e-11,atol=1e-11)
            assert rhs.snapshot()['fallback_calls'] == 0
            assert rhs.snapshot()['fallback_reason'] is None
            break
    else: pytest.fail('Fixture never entered transition')
