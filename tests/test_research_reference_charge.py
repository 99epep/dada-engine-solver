"""Atmospheric filling uses the assembled model's simultaneous gas volumes."""
from dataclasses import replace
import math
import tomllib
from types import SimpleNamespace

import pytest

from dada_solver.research.charge import POLICY, maximum_total_volume
from dada_solver.research.presets import initialize_kinematics, initialize_external_stream
from dada_solver.research.schema import load_study, compile_study, candidate_for_values
from dada_solver.research.study_io import dumps


def derived(raw):
    raw['policies']['charge'] = POLICY
    raw['charge_reference'] = dict(pressure_pa=100000.,temperature_k=293.15,
                                   volume_state='maximum_total_gas_volume')
    raw['parameters'] = [r for r in raw['parameters'] if r['name'] != 'charge.total_mass_kg']
    raw['warm_start']['initial_source'] = 'uniform'


def write(path, raw):
    path.write_text(dumps(raw))
    return load_study(path)


def design(study, **updates):
    d = compile_study(study)
    p = dict(d.fixed_parameters, **{x.name:x.initial for x in d.space.parameters})
    p.update(updates)
    return d.adapter.build(p)


@pytest.fixture
def configured(tmp_path):
    path = initialize_kinematics(tmp_path/'study.toml','harmonic','harmonic')
    raw = tomllib.loads(path.read_text())
    derived(raw)
    return path, raw


@pytest.mark.parametrize('version', [2,3])
def test_explicit_inventory_unchanged(tmp_path, version):
    path = (initialize_kinematics if version==2 else initialize_external_stream)(tmp_path/'study.toml')
    study = load_study(path); before = study.study_id
    built = design(study)
    assert built.charge_diagnostics is None
    assert built.configuration.charge.total_mass == study.fixed_parameters['charge.total_mass_kg']
    assert load_study(path).study_id == before


@pytest.mark.parametrize('edit,message', [
    (lambda r:r.pop('charge_reference'),'charge_reference'),
    (lambda r:r['parameters'].append(dict(name='charge.total_mass_kg',value=.001,unit='kg')),'must not be declared'),
    (lambda r:r['parameters'].append(dict(name='charge.total_mass_kg',initial=.001,lower=.0005,upper=.002,kind='continuous',transform='linear',unit='kg')),'must not be declared'),
    (lambda r:r['warm_start'].update(initial_source='source_exact'),'source_exact'),
    (lambda r:r['charge_reference'].update(pressure_pa=0.),'positive'),
    (lambda r:r['charge_reference'].update(temperature_k=-1.),'positive'),
    (lambda r:r['charge_reference'].update(volume_state='independent_maxima'),'volume_state'),
])
def test_schema_rejections(configured, edit, message):
    path, raw = configured; edit(raw)
    with pytest.raises(ValueError,match=message): write(path,raw)


def test_nonfinite_reference(configured):
    path,raw=configured
    path.write_text(dumps(raw).replace('100000.0','inf'))
    with pytest.raises(ValueError): load_study(path)


def test_harmonic_simultaneous_maximum_and_hx(configured):
    path, raw = configured
    study = write(path,raw)
    a = design(study, **{'kinematics.small.phase_rad':.713,'kinematics.large.phase_rad':2.247})
    c = a.charge_diagnostics
    w = a.build(); m = getattr(w,'model',w); v = m.machine_volumes
    s,l = v.small_cylinder,v.large_cylinder
    amplitude = .5*math.sqrt(s.swept**2+l.swept**2+2*s.swept*l.swept*math.cos(.713-2.247))
    expected = (s.minimum+s.maximum+l.minimum+l.maximum)/2+amplitude+v.cold_heat_exchanger+v.hot_heat_exchanger
    assert c['reference_total_gas_volume_m3'] == pytest.approx(expected,rel=2e-13)
    assert expected < s.maximum+l.maximum+v.cold_heat_exchanger+v.hot_heat_exchanger
    assert m.volumes(c['reference_angle_rad']).total == c['reference_total_gas_volume_m3']
    for side,name in (('heat_in','cold_heat_exchanger'),('heat_out','hot_heat_exchanger')):
        assert getattr(v,name) == getattr(a,side).build().gas_volume_m3
    assert a.configuration.charge.temperature == 293.15
    assert a.configuration.charge.total_mass*m.gas.gas_constant*293.15/expected == pytest.approx(100000.,rel=2e-13)
    assert design(study, **{'kinematics.small.phase_rad':.713,'kinematics.large.phase_rad':2.247}).charge_diagnostics == c
    raw['screening']['samples'] *= 2
    assert design(write(path,raw), **{'kinematics.small.phase_rad':.713,'kinematics.large.phase_rad':2.247}).charge_diagnostics == c


def test_geometry_and_motion_change_inventory(configured):
    path,raw=configured; study=write(path,raw)
    a=design(study, **{'kinematics.small.phase_rad':0.,'kinematics.large.phase_rad':0.})
    b=design(study, **{'kinematics.small.phase_rad':0.,'kinematics.large.phase_rad':math.pi})
    assert b.configuration.charge.total_mass < a.configuration.charge.total_mass
    base=design(study)
    larger=design(study, **{'microtube.heat_in.tube_length_m':base.heat_in.bank.tube_length_m*2})
    delta=larger.heat_in.build().gas_volume_m3-base.heat_in.build().gas_volume_m3
    assert larger.charge_diagnostics['reference_total_gas_volume_m3']-base.charge_diagnostics['reference_total_gas_volume_m3'] == pytest.approx(delta,abs=1e-15)
    for d in (base,larger):
        c=d.charge_diagnostics
        assert d.configuration.charge.total_mass*d.configuration.gas.gas_constant*293.15/c['reference_total_gas_volume_m3'] == pytest.approx(100000.)


def test_identity_includes_reference(configured):
    path,raw=configured; a=write(path,raw)
    raw['charge_reference']['pressure_pa']=101325.
    b=write(path,raw)
    assert a.study_id != b.study_id
    assert candidate_for_values(compile_study(a),{}).candidate_id != candidate_for_values(compile_study(b),{}).candidate_id
    assert 'charge.total_mass_kg' not in a.fixed_parameters


def test_flat_and_breakpoint_maxima():
    for kink in (False,True):
        point=.123456
        def value(t): return 2.-abs(((t-point+math.pi)%(2*math.pi))-math.pi)/10 if kink else 2.
        kin=SimpleNamespace(breakpoint_angles=lambda:(point,),
            small_cylinder_volume_derivative=lambda t: .1 if ((t-point+math.pi)%(2*math.pi))-math.pi < 0 else -.1,
            large_cylinder_volume_derivative=lambda t:0.)
        m=SimpleNamespace(kinematics=kin,volumes=lambda t:SimpleNamespace(total=value(t)))
        v,a=maximum_total_volume(m)
        assert v==pytest.approx(2.,abs=1e-13)
        assert a==pytest.approx(point if kink else 0.,abs=1e-12)


@pytest.mark.parametrize('warm', [False,True])
def test_inventory_reaches_evaluator_and_report(configured, monkeypatch, warm):
    import numpy as np
    from dada_solver.campaign.evaluator import MachineEvaluator, EvaluationControl, _state_record
    path,raw=configured; study=write(path,raw); definition=compile_study(study)
    d=design(study); mass=d.configuration.charge.total_mass
    w=d.build(); caps=[w.heat_in.wall_capacity_j_k,w.heat_out.wall_capacity_j_k]
    previous=()
    if warm:
        old=np.array([mass,200.,mass,300.,mass,400.,mass,500.,caps[0]*330.,caps[1]*310.])
        state=_state_record(old,'four_gas_volumes_m_U_plus_H_i_H_o_wall_energy',
            'microtube_wall_10_state','motor',periodic=True,wall_capacities=caps)
        previous=(dict(candidate_id='old',normalized=[],status='feasible',converged=True,final_periodic_state=state),)
    def solve(wrapper,state,**kwargs):
        assert state[:8:2].sum()==pytest.approx(mass,rel=2e-15)
        if warm: np.testing.assert_allclose(state[:8],old[:8]/4,rtol=2e-15)
        return SimpleNamespace(status='interrupted',message='bounded test',history=(),
            last_complete_state=None,converged=False,backend_statistics={})
    monkeypatch.setattr('dada_solver.campaign.evaluator.solve_periodic_wall_machine',solve)
    result=MachineEvaluator(definition).evaluate_with_control(candidate_for_values(definition,{}),
        EvaluationControl(previous_records=previous))
    assert result['derived']['charge']==d.charge_diagnostics
    assert result['metrics']['total_mass_kg']==mass
    assert result['status']=='budget_exhausted'


def test_additional_hx_volume_replaces_basis_seed(configured):
    from tests.test_research_v3 import change_basis
    path,raw=configured; study=write(path,raw); a=design(study)
    def change(b):
        b['configuration']['machine_volumes']['cold_heat_exchanger']*=10
        b['heat_in']['bank']['additional_internal_volume_m3']+=.001
    change_basis(path,change)
    b=design(load_study(path))
    assert b.charge_diagnostics['reference_total_gas_volume_m3']-a.charge_diagnostics['reference_total_gas_volume_m3']==pytest.approx(.001,abs=1e-15)


@pytest.mark.parametrize('family',['slider_crank','four_bar','six_bar','free_spline','fourier_c2',
    'structured_c2_15p','ideal_piecewise','four_stage','independent_four_stage','hybrid_compact'])
def test_supported_families_use_same_reference_method(tmp_path,family):
    path=initialize_kinematics(tmp_path/'study.toml',family,family)
    raw=tomllib.loads(path.read_text()); derived(raw)
    d=design(write(path,raw)); w=d.build(); model=getattr(w,'model',w)
    c=d.charge_diagnostics
    assert model.volumes(c['reference_angle_rad']).total==c['reference_total_gas_volume_m3']
    assert c['derived_total_mass_kg']>0


def test_reservoir_warm_inventory_is_rescaled(configured,monkeypatch):
    import numpy as np
    from tests.test_research_v3 import change_basis
    from dada_solver.campaign.evaluator import MachineEvaluator,EvaluationControl,_state_record
    from dada_solver.integration import IntegrationInterrupted
    path,raw=configured; write(path,raw)
    change_basis(path,lambda b:b.update(heat_in=None,heat_out=None,warm_start=None))
    raw=tomllib.loads(path.read_text())
    raw['parameters']=[r for r in raw['parameters'] if not r['name'].startswith(('microtube.','thermal.','external_stream.'))]
    study=write(path,raw); definition=compile_study(study); d=design(study)
    mass=d.configuration.charge.total_mass
    old=np.array([mass,200.,mass,300.,mass,400.,mass,500.])
    state=_state_record(old,'four_gas_volumes_m_U','reservoir_8_state','motor',periodic=True)
    previous=(dict(candidate_id='old',normalized=[],status='feasible',converged=True,final_periodic_state=state),)
    def solve(config,initial,*args,**kwargs):
        assert initial.total_mass==pytest.approx(mass)
        np.testing.assert_allclose(initial.as_array(),old/4,rtol=2e-15)
        raise IntegrationInterrupted('bounded test')
    monkeypatch.setattr('dada_solver.campaign.evaluator.evaluate_configuration',solve)
    result=MachineEvaluator(definition).evaluate_with_control(candidate_for_values(definition,{}),
        EvaluationControl(previous_records=previous))
    assert result['metrics']['total_mass_kg']==mass


def test_derived_policy_rescale_keeps_mass_derived(configured,monkeypatch,tmp_path):
    from dada_solver.research.cli import evaluate
    from dada_solver.research.rescale import rescale
    from dada_solver.research.report import inspect,render_html
    path,raw=configured; study=write(path,raw)
    monkeypatch.setattr('dada_solver.campaign.evaluator.solve_periodic_wall_machine',lambda *a,**k:
        SimpleNamespace(status='interrupted',message='bounded test',history=(),
                        last_complete_state=None,converged=False,backend_statistics={}))
    output=tmp_path/'evaluation.json'; record=evaluate(path,output)
    destination=tmp_path/'report.html'
    render_html(inspect(output),destination)
    html=destination.read_text()
    assert 'Gas inventory [kg]' in html and 'total_mass_kg' in html
    target=tmp_path/'scaled.toml'
    rescale(output,record['candidate_id'],1.,target)
    scaled=load_study(target)
    assert 'charge.total_mass_kg' not in scaled.fixed_parameters
    assert design(scaled).charge_diagnostics==design(study).charge_diagnostics
