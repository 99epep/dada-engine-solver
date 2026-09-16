"""Analytic limits, source regressions and connected variable-film checks."""
from dataclasses import replace
from pathlib import Path
import math
import numpy as np
import pytest
from dada_solver.exchangers.gas_transport import DiluteGasTransport
from dada_solver.exchangers.gas_correlations import (
    MicrotubeGasModel, MicrotubeDomainError, GasSurfaceAccommodation, SecondOrderSlip,
    compressible_poiseuille, knudsen_regime, graur_silica_fit, laminar_entry_nusselt, gnielinski)
from dada_solver.exchangers.microtube_geometry import MicrotubeBank
from dada_solver.exchangers.hardware import load_hardware_definition, build_exchanger, TubeHalfLink
from dada_solver.exchangers.doty import validate_gas_hydraulics
from dada_solver.fluids import CaloricallyPerfectGas

ROOT=Path(__file__).resolve().parents[1]
BANK=MicrotubeBank(309,.127,.00033,.0001524,.00125,.001)


def test_compressible_poiseuille_recovers_incompressible_limit():
    tr=DiluteGasTransport();p=2e5;t=350;dp=1e-3;mu=tr.viscosity(t)
    actual=compressible_poiseuille(p+dp/2,p-dp/2,t,mu,.127,.00033,309,tr.gas_constant)
    expected=dp*(p/(tr.gas_constant*t))*309*math.pi*.00033**4/(128*mu*.127)
    assert actual==pytest.approx(expected,rel=1e-7)


def test_compressible_poiseuille_forward_reverse_symmetry():
    args=(350,2e-5,.127,.00033,309,287.05)
    assert compressible_poiseuille(2e5,1e5,*args)==-compressible_poiseuille(1e5,2e5,*args)
    assert compressible_poiseuille(2e5,2e5,*args)==0


def test_legacy_mean_density_is_exact_pressure_squared_identity():
    tr=DiluteGasTransport();gas=CaloricallyPerfectGas(tr.gas_constant,1005,1005-tr.gas_constant)
    legacy=TubeHalfLink(BANK,tr.viscosity(350),1,0)
    modern=replace(legacy,gas_model=MicrotubeGasModel(tr))
    for p in (200001,200100,201000):
        assert modern.directed_flow(p,2e5,350,gas).mass_flow_rate==pytest.approx(legacy.directed_flow(p,2e5,350,gas).mass_flow_rate,rel=1e-12)


@pytest.mark.parametrize('kn,regime',[(0,'continuum'),(.0009,'continuum'),(.001,'slip_onset'),(.01,'slip'),(.1,'slip'),(.101,'beyond_continuum_model')])
def test_knudsen_regime_classification(kn,regime): assert knudsen_regime(kn)==regime


def test_slip_factor_tends_to_one_as_kn_tends_to_zero():
    model=graur_silica_fit('nitrogen')
    assert model.factor(0,5)==1
    assert model.factor(1e-12,1)==pytest.approx(1,abs=2e-11)
    assert model.factor(.02,1+1e-12)==pytest.approx(model.factor(.02,1),rel=1e-12)


def check_graur(gas,a,b):
    fit=graur_silica_fit(gas)
    for kn in (.01,.03,.1,.25):
        assert fit.factor(kn,5)==pytest.approx(1+a*kn+b*kn*kn,rel=1e-14)
    assert fit.factor(.2,5)>1+a*.2


def test_second_order_slip_reproduces_graur_n2(): check_graur('nitrogen',11.668,16.626)
def test_second_order_slip_reproduces_graur_he(): check_graur('helium',10.812,9.156)
def test_second_order_slip_reproduces_graur_ar(): check_graur('argon',13.218,24.274)


def test_no_silica_or_accommodation_assumption_for_metal():
    assert GasSurfaceAccommodation().momentum_accommodation is None
    with pytest.raises(ValueError,match='conventions'):
        MicrotubeGasModel(slip=graur_silica_fit('nitrogen'))
    with pytest.raises(MicrotubeDomainError,match='unknown'):
        MicrotubeGasModel().slip_factor(1000,900,300,1e-5)
    assert MicrotubeGasModel().slip_factor(3e5,2.9e5,300,.00033)==1


def test_nusselt_fully_developed_laminar_limit(): assert laminar_entry_nusselt(0)==3.66

def test_thermal_entry_correction_tends_to_fully_developed_limit():
    assert laminar_entry_nusselt(1e-12)==pytest.approx(3.66)
    assert laminar_entry_nusselt(100)>laminar_entry_nusselt(10)>3.66


def test_no_silent_use_outside_correlation_domain():
    m=MicrotubeGasModel();tr=m.transport
    flow=2500*309*math.pi*.00033*tr.viscosity(300)/4
    d=m.diagnose(BANK,flow,3e5,2.99e5,300)
    assert d.nusselt is None and d.model_validity=='invalid'
    with pytest.raises(MicrotubeDomainError): m.require(d)
    rare=m.diagnose(BANK,1e-8,100,99,300)
    assert 'thermal_slip_not_implemented' in rare.issues
    with pytest.raises(MicrotubeDomainError): m.require(rare)
    with pytest.raises(MicrotubeDomainError): gnielinski(2000,.7)
    with pytest.raises(MicrotubeDomainError): gnielinski(1e4,.1)


@pytest.mark.parametrize('species',['air','nitrogen','argon','helium'])
def test_transport_properties_change_with_temperature(species):
    t=DiluteGasTransport(species)
    assert t.viscosity(600)>t.viscosity(300)>0
    assert t.conductivity(600)>t.conductivity(300)>0
    assert t.cp(600)>=t.cp(300)>t.gas_constant
    assert t.mean_free_path(1e5,600)>t.mean_free_path(1e5,300)
    with pytest.raises(ValueError): t.viscosity(1200)


def check_doty(fluid):
    rows=validate_gas_hydraulics(ROOT/f'examples/data/doty_1991_{fluid}_reference.csv')
    assert len(rows)==6
    for r in rows:
        assert abs(r['inverse_flow_residual_kg_s'])<1e-15
        assert r['source'] and r['predicted']>0
        assert r['model_validity']['knudsen']<.001
        assert isinstance(r['inside_experimental_uncertainty'],bool)
        assert r['thermal_validation']['predicted_UA_W_K'] is None
    return rows


def test_doty_nitrogen_hydraulics():
    rows=check_doty('nitrogen')
    for r in rows:
        assert abs(r['relative_error'])<.5


def test_doty_helium_hydraulics():
    rows=check_doty('helium')
    variable=[r for r in rows if r['model_id']=='compressible_variable_transport']
    np.testing.assert_allclose([r['predicted'] for r in variable],[3318.8793,2205.7510,5477.2102],rtol=.001)
    assert all(r['relative_error']>0 and not r['inside_experimental_uncertainty'] for r in variable)


def test_variable_film_requires_and_uses_flow_context_and_preserves_wall_balance():
    _,bank,inputs,_=load_hardware_definition((ROOT/'examples/motor_hardware.toml').read_text(),1005)
    thermal,_=build_exchanger(bank,replace(inputs,gas_model=MicrotubeGasModel()))
    energy=thermal.wall_capacity_j_k*400
    with pytest.raises(ValueError,match='context'): thermal.rates(300,energy)
    paused=dict(frequency_hz=2,passages=((0.,2e5,2e5),(0.,2e5,2e5)))
    moving=dict(frequency_hz=2,passages=((1e-4,200100,2e5),(1e-4,200100,2e5)))
    r0=thermal.rates(300,energy,context=paused);r1=thermal.rates(300,energy,context=moving)
    assert r1['gas_heat_w']>r0['gas_heat_w']>0
    assert r1['wall_energy_rate_w']==pytest.approx(r1['air_heat_w']-r1['gas_heat_w'])
    poisoned=replace(thermal,gas_wall_conductance_w_k=1e15)
    assert poisoned.rates(300,energy,context=moving)==r1


def test_hardware_toml_selects_variable_model_explicitly():
    text=(ROOT/'examples/motor_hardware.toml').read_text()
    _,_,legacy,_=load_hardware_definition(text,1005)
    _,_,modern,_=load_hardware_definition(text+'\n[gas_model]\nmode="variable_properties"\nspecies="air"\n',1005)
    assert legacy.gas_model is None
    assert isinstance(modern.gas_model,MicrotubeGasModel)


def test_turbulent_link_matches_explicit_darcy_and_rejects_high_compressibility():
    from dada_solver.exchangers.gas_correlations import darcy_smooth
    tr=DiluteGasTransport();t=300.;p=1e6;mu=tr.viscosity(t)
    bank=replace(BANK,inner_diameter_m=.003,pitch_m=.004,tube_length_m=.6)
    area=bank.dimensions()['tube_flow_area_m2'];rho=p/(tr.gas_constant*t)
    expected=1e4*area*mu/bank.inner_diameter_m
    dp=darcy_smooth(1e4)*(bank.tube_length_m/2)/bank.inner_diameter_m*expected**2/(2*rho*area**2)
    link=TubeHalfLink(bank,mu,1,0,gas_model=MicrotubeGasModel(tr))
    gas=CaloricallyPerfectGas(tr.gas_constant,1005,1005-tr.gas_constant)
    result=link.directed_flow(p+dp/2,p-dp/2,t,gas)
    assert result.mass_flow_rate==pytest.approx(expected,rel=1e-8)
    d=link.gas_model.diagnose(bank,expected,p*1.2,p*.8,t)
    assert 'large_relative_pressure_drop' in d.issues
    with pytest.raises(MicrotubeDomainError): link.gas_model.require(d)


def test_heat_fraction_diagnostics_use_time_weights_not_adaptive_sample_counts():
    from types import SimpleNamespace
    from dada_solver.exchangers.gas_diagnostics import cycle_microtube_diagnostics
    from dada_solver.geometry import CylinderVolumeLimits
    from dada_solver.state import ThermodynamicState
    _,bank,inputs,_=load_hardware_definition((ROOT/'examples/motor_hardware.toml').read_text(),1005)
    wall,_=build_exchanger(bank,replace(inputs,gas_model=MicrotubeGasModel()))
    gas=CaloricallyPerfectGas(287.,1005.,718.)
    # Deliberately nonuniform samples: the first interval occupies 10%, not 50%.
    angles=np.array([0.,.1,1.]);mass=1e-4;temp=300
    trajectory=np.tile(np.array([mass,mass*718*temp]*4+[wall.wall_capacity_j_k*400]*2)[:,None],(1,3))
    def contexts(angle,values):
        m=0. if angle<.2 else 1e-4
        return (dict(frequency_hz=2,passages=((m,2e5,2e5),(m,2e5,2e5))),)*2
    wrapper=SimpleNamespace(heat_in=wall,heat_out=wall,model=SimpleNamespace(gas=gas,angular_speed=1.),flow_contexts=contexts,
        thermal_rates=lambda a,v:(wall.rates(300,v[8],context=contexts(a,v)[0]),)*2)
    result=cycle_microtube_diagnostics(wrapper,angles,trajectory)
    for passage in result['passages'].values():
        fractions=passage['domains']['correlation_id']
        assert fractions['stagnant_radial_screening']['time_fraction']==pytest.approx(.55)
        assert sum(x['time_fraction'] for x in fractions.values())==pytest.approx(1)
        assert sum(x['absolute_heat_fraction'] for x in fractions.values())==pytest.approx(1)
        assert passage['ranges']['strouhal']['maximum']>0


def test_explicit_slip_uncertainty_scenarios_change_flow_without_metal_calibration():
    tr=DiluteGasTransport();gas=CaloricallyPerfectGas(tr.gas_constant,1005,1005-tr.gas_constant)
    bank=replace(BANK,inner_diameter_m=3e-6,pitch_m=.001,tube_length_m=.03)
    scenario=SecondOrderSlip(1.,.1,'Illustrative HS sensitivity, not measured metal accommodation','viscosity_hard_sphere')
    low=TubeHalfLink(bank,tr.viscosity(300),1,0,gas_model=MicrotubeGasModel(tr,slip=scenario))
    high=replace(low,gas_model=replace(low.gas_model,slip=replace(scenario,a1=2.,a2=.3)))
    m0=low.directed_flow(110000,100000,300,gas).mass_flow_rate
    m1=high.directed_flow(110000,100000,300,gas).mass_flow_rate
    assert 0<m0<m1
    assert low.gas_model.accommodation.momentum_accommodation is None
    # Hydraulic slip does not manufacture a thermal temperature-jump closure.
    with pytest.raises(MicrotubeDomainError,match='thermal_slip'):
        low.gas_model.require(low.gas_model.diagnose(bank,m0,110000,100000,300))
