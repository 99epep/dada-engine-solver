"""Independent HeatLib oracle, Shah Eq.192 and physical segment ownership."""
from dataclasses import replace
from pathlib import Path
import json
import math
import numpy as np
import pytest
from dada_solver import numerical_primitives as n
from dada_solver.exchangers.gas_correlations import MicrotubeGasModel
from dada_solver.exchangers.hardware import TubeHalfLink
from dada_solver.exchangers.microtube_geometry import MicrotubeBank
from dada_solver.fluids import CaloricallyPerfectGas

REFERENCE=json.loads((Path(__file__).parent/'data/bennett_heatlib_reference.json').read_text())
BANK=MicrotubeBank(1000,.02,.00077,.0001524,.00125,.001)
MODEL=MicrotubeGasModel()
GAS=CaloricallyPerfectGas(287.05,1005.,717.95)

@pytest.mark.parametrize('row',REFERENCE['values'])
def test_bennett_author_octave_reference(row):
    pr,z=row['prandtl'],row['inverse_graetz']
    assert n.bennett_mean_nusselt(1000,pr,1/(1000*pr*z))==pytest.approx(row['nusselt'],rel=3e-14)


def test_limits_and_entry_boundary():
    assert n.bennett_mean_nusselt(0,.7,1)==3.66
    assert n.bennett_mean_nusselt(1000,.7,1e-12)==pytest.approx(3.66,rel=1e-10)
    assert n.bennett_mean_nusselt(1000,.7,.1)>n.bennett_mean_nusselt(1000,.7,.01)>3.66
    sides=[n.bennett_mean_nusselt(1000,.7,1/(50*(1+eps))) for eps in (-1e-8,1e-8)]
    assert sides[0]==pytest.approx(sides[1],rel=2e-8)

@pytest.mark.parametrize('xp',[1e-8,1e-5,.001,.01,.05,.1,1.,100.,1e8])
def test_shah_original_eq192_and_darcy(xp):
    re=1200.
    f_re=3.44/math.sqrt(xp)+(1.25/(4*xp)+16-3.44/math.sqrt(xp))/(1+.00021/xp**2)
    assert n.shah_apparent_darcy(re,xp)==pytest.approx(4*f_re/re,rel=2e-14)
    if xp==1e8: assert n.shah_apparent_darcy(re,xp)==pytest.approx(64/re,rel=3e-10)

@pytest.mark.parametrize('re',[.01,10,600,2200])
def test_segment_additivity_and_no_midpoint_restart(re):
    mu=MODEL.transport.viscosity(300);d=BANK.inner_diameter_m;a=BANK.tube_flow_area_m2
    m=re*a*mu/d;length=BANK.tube_length_m
    loss=lambda x,y:n.shah_segment_pressure_loss(m,d,a,mu,1.2,x,y)
    assert loss(0,length)==pytest.approx(loss(0,length/2)+loss(length/2,length),rel=2e-15)
    assert loss(0,length)<2*loss(0,length/2)
    assert n.shah_entry_excess(0)==0
    # In the long-tube limit the entrance penalty is negligible relative to Poiseuille.
    ratio=loss(0,1e6)/(32*mu*1e6*m/(1.2*a*d*d))
    assert ratio<1e-6


def test_developing_laminar_valid_but_transition_and_turbulent_entry_rejected():
    mu=MODEL.transport.viscosity(300);a=BANK.tube_flow_area_m2;d=BANK.inner_diameter_m
    diag=lambda re,bank=BANK:MODEL.diagnose(bank,re*a*mu/d,1e6,1e6,300)
    result=diag(1000)
    assert result.hydrodynamic_developing and not result.issues
    assert result.correlation_id=='bennett_2020_combined_entry_constant_wall'
    MODEL.require(result)
    assert 'hydrodynamic_entry_unresolved' in diag(3000).issues
    assert 'turbulent_entry_unresolved' in diag(4500,replace(BANK,tube_length_m=.005)).issues
    assert 'inverse_graetz_outside_domain' in diag(1000,replace(BANK,tube_length_m=1e-9)).issues


def test_reverse_mirrors_segment_and_does_not_restart():
    first=TubeHalfLink(BANK,MODEL.transport.viscosity(300),1,0,gas_model=MODEL,axial_half=0)
    second=replace(first,axial_half=1)
    p1,p2=100010.,100000.
    forward=first.directed_flow(p1,p2,300,GAS).mass_flow_rate
    downstream=second.directed_flow(p1,p2,300,GAS).mass_flow_rate
    assert downstream>forward
    assert -second.bidirectional_flow(p2,p1,300,300,GAS).mass_flow_rate==forward
    assert -first.bidirectional_flow(p2,p1,300,300,GAS).mass_flow_rate==downstream


def test_production_assembly_assigns_physical_halves():
    from tests.test_solver_acceleration import variable_wrapper
    w=variable_wrapper();m=w.model
    assert m.small_cold_link.axial_half==m.large_hot_link.axial_half==0
    assert m.cold_large_valve.flow_model.axial_half==m.hot_small_valve.flow_model.axial_half==1


def test_developing_primitives_compiled_parity():
    numba=pytest.importorskip('numba')
    from dada_solver.wall_backend import compiled_dispatcher
    compiled_dispatcher()
    heat=numba.njit(n.bennett_mean_nusselt);loss=numba.njit(n.shah_segment_pressure_loss)
    for row in REFERENCE['values']:
        pr,z=row['prandtl'],row['inverse_graetz']
        assert heat(1000.,pr,1/(1000*pr*z))==pytest.approx(row['nusselt'],rel=3e-14)
    for x1,x2 in ((0.,.01),(.01,.02),(0.,.02)):
        args=(.001,.00077,BANK.tube_flow_area_m2,2e-5,1.2,x1,x2)
        assert loss(*args)==pytest.approx(n.shah_segment_pressure_loss(*args),rel=2e-14)


@pytest.mark.parametrize('half', [0, 1])
@pytest.mark.parametrize('reverse', [False, True])
def test_compiled_segment_direction(half, reverse):
    numba=pytest.importorskip('numba')
    from dada_solver.wall_backend import compiled_dispatcher
    compiled_dispatcher()
    link=TubeHalfLink(BANK,MODEL.transport.viscosity(300),1,0,gas_model=MODEL,axial_half=half)
    p=np.array([BANK.inner_diameter_m,BANK.tube_length_m,BANK.tube_count,
        BANK.tube_flow_area_m2,1,0,0,.3,.2,1,100,1000,0,half],dtype=float)
    ok,flow=numba.njit(n.directed_flow)(100010.,100000.,300.,p,GAS.gas_constant,GAS.heat_capacity_ratio,reverse)
    assert ok
    assert flow==pytest.approx(link.directed_flow(100010.,100000.,300.,GAS,reverse=reverse).mass_flow_rate,rel=2e-11)


@pytest.mark.parametrize('hi,ho',[('downstream','downstream'),('upstream','downstream'),
    ('downstream','upstream'),('upstream','upstream')])
def test_developing_full_rhs_without_fallback(hi,ho):
    pytest.importorskip('numba')
    from tests.test_compiled_transition import transition_wrapper
    from tests.test_wall_backend import values
    from dada_solver.wall_backend import WallRHS, WallBackendSettings
    w=transition_wrapper(hi,ho);m=w.model
    def link(x): return replace(x,bank=replace(x.bank,tube_length_m=.02))
    m=replace(m,small_cold_link=link(m.small_cold_link),large_hot_link=link(m.large_hot_link),
        cold_large_valve=replace(m.cold_large_valve,flow_model=link(m.cold_large_valve.flow_model)),
        hot_small_valve=replace(m.hot_small_valve,flow_model=link(m.hot_small_valve.flow_model)))
    def wall(x): return replace(x,gas_film=replace(x.gas_film,bank=replace(x.gas_film.bank,tube_length_m=.02)))
    w=replace(w,model=m,heat_in=wall(w.heat_in),heat_out=wall(w.heat_out))
    rhs=WallRHS(w,WallBackendSettings('numba'))
    for delta in (-.0001,.0001):
        x=values(w,pressure=1e6);x[4:6]*=1+delta
        expected=w.derivative(0.,x)
        ds=[d for wall,context in zip((w.heat_in,w.heat_out),w.flow_contexts(0.,x))
            for d in wall.gas_film.evaluate(350.,350.,context)[1]]
        assert any(d.hydrodynamic_developing and d.reynolds<2300 for d in ds)
        rhs.implementation.reference=lambda *_:pytest.fail('Unexpected Python fallback')
        np.testing.assert_allclose(rhs(0.,x),expected,rtol=2e-11,atol=1e-11)
    assert rhs.snapshot()['fallback_calls']==0


def test_nearly_equal_pressures_do_not_lose_root_bracket():
    link=TubeHalfLink(BANK,MODEL.transport.viscosity(300),1,0,gas_model=MODEL,axial_half=1)
    flow=link.directed_flow(2e5,math.nextafter(2e5,0),300,GAS).mass_flow_rate
    assert 0<flow<1e-12


@pytest.mark.parametrize('placement', ['downstream','upstream'])
def test_external_stream_production_assigns_physical_halves(tmp_path,placement):
    from dada_solver.research.presets import initialize_v3
    from dada_solver.research.schema import compile_study,load_study
    d=compile_study(load_study(initialize_v3(tmp_path/'study.toml')))
    design=d.adapter.build(dict(d.fixed_parameters))
    for exchanger in (design.heat_in,design.heat_out):
        parts=replace(exchanger,valve_placement=placement).build()
        assert parts.inlet.axial_half==0
        assert parts.outlet.axial_half==1
        metadata=dict(parts.metadata)
        assert metadata['laminar_hydraulic_correlation']=='compressible_poiseuille_plus_shah_london_1978_eq192'
        assert 'no_midpoint_restart' in metadata['entrance_segmentation']


@pytest.mark.parametrize('half',[0,1])
@pytest.mark.parametrize('length_over_diameter',[115.,200.,1000.])
def test_transition_segment_pressure_loss_monotonic_and_continuous(half,length_over_diameter):
    d=.001;a=1e-3;mu=2e-5;rho=2.;length=length_over_diameter*d/2
    linear=32*mu*length/(rho*a*d*d)
    def pressure(re):
        return n.tube_network_residual(re*a*mu/d,d,a,mu,linear,0.,1.,length,rho,0.,half*length,True)
    points=sorted([*np.linspace(2200,4100,301),2300-1e-6,2300.,2300+1e-6,4000-1e-6,4000.,4000+1e-6])
    assert np.all(np.diff([pressure(re) for re in points])>0)
    for edge in (2300,4000):
        assert pressure(edge-1e-6)==pytest.approx(pressure(edge+1e-6),rel=3e-9)


@pytest.mark.parametrize('re',[1000,2300,3000,4000,5000])
def test_transition_preserves_axial_segment_additivity(re):
    d=.001;a=1e-3;mu=2e-5;rho=2.;length=.2;flow=re*a*mu/d
    def loss(start,end):
        linear=32*mu*(end-start)/(rho*a*d*d)
        return n.tube_network_residual(flow,d,a,mu,linear,0.,1.,end-start,rho,0.,start,True)
    assert loss(0,length)==pytest.approx(loss(0,length/2)+loss(length/2,length),rel=2e-15)
