"""Compiled continuum flow/film parity; rejection remains the Python authority."""
from dataclasses import replace
import math

import numpy as np
import pytest

from dada_solver import numerical_primitives as numeric
from dada_solver.exchangers.gas_correlations import MicrotubeGasModel, MicrotubeDomainError
from dada_solver.exchangers.gas_transport import DiluteGasTransport
from dada_solver.exchangers.hardware import TubeHalfLink
from dada_solver.exchangers.microtube_geometry import MicrotubeBank
from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.wall_backend import WallRHS, WallBackendSettings, compiled_dispatcher
from tests.test_solver_acceleration import variable_wrapper
from tests.test_wall_backend import values


@pytest.fixture(scope='module')
def compiled():
    numba = pytest.importorskip('numba')
    compiled_dispatcher()  # Registers the same transitive helpers as production.
    return numba.njit(numeric.directed_flow), numba.njit(numeric.continuum_diagnostics)


def flow_case(re, species='air', cda_ratio=0., header=0.):
    tr = DiluteGasTransport(species)
    model = MicrotubeGasModel(transport=tr)
    bank = MicrotubeBank(1000, .8, .000770, .0001524, .00125, .001)
    cv_ratio = 1.5 if species in ('helium', 'argon') else 2.5
    gas = CaloricallyPerfectGas(tr.gas_constant, (cv_ratio+1)*tr.gas_constant, cv_ratio*tr.gas_constant)
    t = 350.; area = bank.tube_flow_area_m2; d = bank.inner_diameter_m
    mu = tr.viscosity(t); flow = re*area*mu/d
    f = 64/re if re<2300 else numeric.transition_friction(re) if re<4000 else numeric.turbulent_darcy(re)
    cda = area*cda_ratio
    # Invert tube + header + optional legacy valve losses at mean ideal density.
    coefficient = f*(bank.tube_length_m/2)/d/(2*area**2) + header/(4*area**2)
    if cda: coefficient += 1/(2*cda**2)
    pout = 1e7
    pin = math.sqrt(pout**2+2*gas.gas_constant*t*coefficient*flow**2)
    link = TubeHalfLink(bank, mu, 1., header, cda or None, model)
    parameters = np.array([d, bank.tube_length_m, bank.tube_count, area, 1., header, cda,
        .3, .2, 1., tr.minimum_temperature, tr.maximum_temperature, numeric.SPECIES.index(species)])
    return bank, model, gas, link, parameters, pin, pout, t, flow


@pytest.mark.parametrize('re', [2299.999,2300.,2300.001,2500.,3000.,3999.999,4000.,4500.,10000.])
@pytest.mark.parametrize('species', numeric.SPECIES)
def test_compiled_flow_and_film_endpoints(compiled, re, species):
    flow_fn, film_fn = compiled
    bank, model, gas, link, p, pin, pout, t, target = flow_case(re,species)
    expected = link.directed_flow(pin,pout,t,gas).mass_flow_rate
    ok, actual = flow_fn(pin,pout,t,p,gas.gas_constant,gas.heat_capacity_ratio)
    assert ok
    assert actual == pytest.approx(expected,rel=2e-11,abs=1e-15)
    assert actual == pytest.approx(target,rel=2e-10)
    for entry in (True,False):
        selected = replace(model,thermal_entry=entry)
        diag = selected.diagnose(bank,expected,pin,pout,t)
        selected.require(diag)
        ok, nu = film_fn(expected,pin,pout,t,p[0],p[1],p[3],.3,.2,entry,int(p[12]))
        assert ok
        assert nu == pytest.approx(diag.nusselt,rel=2e-13,abs=1e-13)


@pytest.mark.parametrize('header,cda',[(0.,0.),(1.2,0.),(1.2,1.),(3.,.7)])
def test_header_and_legacy_valve_losses(compiled, header, cda):
    fn,_=compiled
    _,_,gas,link,p,pin,pout,t,_=flow_case(3000.,cda_ratio=cda,header=header)
    ok, actual = fn(pin,pout,t,p,gas.gas_constant,gas.heat_capacity_ratio)
    assert ok
    assert actual == pytest.approx(link.directed_flow(pin,pout,t,gas).mass_flow_rate,rel=2e-11,abs=1e-15)
    assert fn(pout,pin,t,p,gas.gas_constant,gas.heat_capacity_ratio)==(True,0.)


def transition_wrapper(placement_in='downstream',placement_out='downstream'):
    w=variable_wrapper(placement_in,placement_out)
    bank=MicrotubeBank(1000,.8,.000770,.0001524,.00125,.001)
    def link(x):return replace(x,bank=bank)
    m=w.model
    m=replace(m,small_cold_link=link(m.small_cold_link),large_hot_link=link(m.large_hot_link),
        cold_large_valve=replace(m.cold_large_valve,flow_model=link(m.cold_large_valve.flow_model)),
        hot_small_valve=replace(m.hot_small_valve,flow_model=link(m.hot_small_valve.flow_model)))
    def wall(x):return replace(x,gas_film=replace(x.gas_film,bank=bank))
    return replace(w,model=m,heat_in=wall(w.heat_in),heat_out=wall(w.heat_out))


@pytest.mark.parametrize('hi,ho',[('downstream','downstream'),('upstream','downstream'),
    ('downstream','upstream'),('upstream','upstream')])
@pytest.mark.parametrize('tabulated', [False, True])
def test_full_rhs_placements_without_python(compiled, hi, ho, tabulated):
    w=transition_wrapper(hi,ho)
    if tabulated:
        from dada_solver.tabulated_fluid import ideal_validation_table
        w=replace(w,model=replace(w.model,gas=ideal_validation_table(w.model.gas)))
    rhs=WallRHS(w,WallBackendSettings('numba'))
    exercised=False
    for delta in np.linspace(.001,.2,60):
        x=values(w,pressure=1e6);x[:2]*=1+delta
        try: expected=w.derivative(0.,x)
        except MicrotubeDomainError:continue
        contexts=w.flow_contexts(0.,x)
        diagnostics=w.heat_in.gas_film.evaluate(350.,350.,contexts[0])[1]
        exercised |= any(d.transition_model_used for d in diagnostics)
        # Prove the kernel does not merely produce parity by falling back.
        rhs.implementation.reference=lambda *_:pytest.fail('Unexpected Python RHS fallback')
        actual=rhs(0.,x)
        np.testing.assert_allclose(actual,expected,rtol=2e-11,atol=1e-11)
        assert abs(sum(actual[:8:2]))<1e-15
        assert abs(sum(actual[1:8:2])+sum(actual[8:10])-sum(actual[10:12])+actual[14])<1e-9
        if exercised: break
    assert exercised
    assert rhs.snapshot()['fallback_calls']==0


@pytest.mark.parametrize('guard',['mach','drop','entry','knudsen','temperature'])
def test_physical_rejection_keeps_python_message(compiled, guard):
    w=transition_wrapper()
    if guard in ('mach','drop'):
        field='maximum_mach' if guard=='mach' else 'maximum_relative_pressure_drop'
        def wall(x):return replace(x,gas_film=replace(x.gas_film,model=replace(x.gas_film.model,**{field:1e-5})))
        w=replace(w,heat_in=wall(w.heat_in),heat_out=wall(w.heat_out))
    elif guard=='entry':
        def wall(x):return replace(x,gas_film=replace(x.gas_film,bank=replace(x.gas_film.bank,tube_length_m=.001)))
        w=replace(w,heat_in=wall(w.heat_in),heat_out=wall(w.heat_out))
    x=values(w,pressure=100. if guard=='knudsen' else 1e6,temperature=99.9 if guard=='temperature' else 350.)
    x[:2]*=1.05
    with pytest.raises(ValueError) as expected:w.derivative(0.,x)
    rhs=WallRHS(w,WallBackendSettings('numba'))
    with pytest.raises(type(expected.value)) as actual:rhs(0.,x)
    assert str(actual.value)==str(expected.value)
    assert rhs.snapshot()['fallback_calls']==1


@pytest.mark.parametrize('ratio', [-.01, float('nan'), float('inf')])
def test_public_transition_preserves_graetz_validation(ratio):
    from dada_solver.exchangers.gas_correlations import transition_nusselt
    with pytest.raises(ValueError, match='Graetz must be nonnegative'):
        transition_nusselt(2500., .7, ratio)
