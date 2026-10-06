"""Schema 3 uses the current ownership and persistent campaign/evaluation engine."""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import tomllib
import pytest
from dada_solver.research.presets import initialize_external_stream
from dada_solver.research.schema import load_study,compile_study,candidate_for_values
from dada_solver.research.study_io import dumps
from dada_solver.campaign.evaluator import MachineEvaluator
from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.campaign.runner import OptimizationCampaign
from dada_solver.research import report
from dada_solver.tabulated_fluid import ideal_validation_table


def rewrite(path,transform):
    raw=tomllib.loads(path.read_text());transform(raw);path.write_text(dumps(raw))


def change_basis(path,transform):
    basis=path.with_suffix('.basis.json');raw=json.loads(basis.read_text());transform(raw)
    text=json.dumps(raw,indent=2)+'\n';basis.write_text(text)
    rewrite(path,lambda s:s['sources']['machine'].update(sha256=hashlib.sha256(text.encode()).hexdigest()))


@pytest.fixture
def study(tmp_path):
    return initialize_external_stream(tmp_path/'study.toml')


def test_schema_ownership_and_table_identity(study):
    original=compile_study(load_study(study))
    def activate(raw):
        row=next(r for r in raw['parameters'] if r['name']=='external_stream.heat_in.mass_flow_kg_s')
        value=row.pop('value');row.update(initial=value,lower=.03,upper=.08,kind='continuous',transform='linear')
    rewrite(study,activate)
    active=load_study(study)
    assert len(active.space.parameters)==1
    assert 'external_stream.heat_in.mass_flow_kg_s' not in active.fixed_parameters
    assert active.fixed_parameters['external_stream.heat_in.cp_j_kg_k']==4180.
    assert active.fixed_parameters['kinematics.small.phase_rad']<0
    table=ideal_validation_table(original.configuration.gas,rho_axis=(.01,.1,1.,10.,50.))
    change_basis(study,lambda b:b['configuration'].update(gas=table.to_data()))
    compiled=compile_study(load_study(study))
    assert compiled.identity['fluid']['table_sha256']==table.content_hash
    assert compiled.definition_id!=original.definition_id
    old=candidate_for_values(compiled,{p.name:p.initial for p in compiled.space.parameters})
    change_basis(study,lambda b:b['configuration']['gas'].update(provenance='New verified analytic provenance'))
    new=compile_study(load_study(study))
    assert candidate_for_values(new,{p.name:p.initial for p in new.space.parameters}).candidate_id!=old.candidate_id


@pytest.mark.parametrize('mutation',[
    lambda s:s['objective'].update(type='maximize_thermal_efficiency'),
    lambda s:s['policies'].update(external_loop_hydraulics='validated_liquid'),
    lambda s:s['parameters'].append(dict(name='external_stream.heat_in.pump_power_w',value=1.,unit='W')),
    lambda s:s.update(schema_version=2),
])
def test_schema_rejects_incompatible_scientific_declarations(study,mutation):
    rewrite(study,mutation)
    with pytest.raises(ValueError): load_study(study)


@pytest.fixture(scope='module')
def refrigerator(tmp_path_factory):
    path=initialize_external_stream(tmp_path_factory.mktemp('refrigerator')/'study.toml')
    definition=compile_study(load_study(path))
    result=MachineEvaluator(definition).evaluate(candidate_for_values(definition,{}))
    return path,definition,result


def test_external_boundary_refrigeration_and_margins(refrigerator):
    _,d,r=refrigerator
    assert r['status']=='feasible'
    m=r['metrics'];assert m['operating_mode']=='refrigeration'
    assert m['cooling_power_w']>0 and m['heating_power_w']>0 and m['indicated_mechanical_input_power_w']>0
    assert r['derived']['total_microtube_count'] == (r['derived']['heat_in_microtube_count'] + r['derived']['heat_out_microtube_count'])
    assert m['cooling_power_per_total_microtube_w'] == pytest.approx(m['cooling_power_w']/r['derived']['total_microtube_count'])
    assert m['cooling_cop']==pytest.approx(m['cooling_power_w']/m['indicated_mechanical_input_power_w'])
    assert m['heating_cop']==pytest.approx(m['cooling_cop']+1,abs=1e-5)
    assert abs(m['conservation']['absolute_energy_residual'])<1e-7
    assert r['periodic_convergence']['last_normalized_periodic_error']<=1
    for side in ('heat_in','heat_out'):
        stream=r['derived']['external_streams'][side]
        assert stream['capacity_rate_w_k']==stream['mass_flow_kg_s']*stream['cp_j_kg_k']
        assert 'liquid' in stream['fluid'] and stream['outlet_minimum_k']>200
    assert r['derived']['external_streams']['heat_in']['mean_heat_into_machine_w']==m['cooling_power_w']
    assert r['derived']['external_streams']['heat_out']['mean_heat_into_machine_w']==-m['heating_power_w']
    assert all('relative_margin' in row for row in r['constraints'])
    assert m['useful_mechanical_power_w'] is None


def test_cooling_power_objective_and_structured_setup(study,refrigerator,tmp_path):
    rewrite(study,lambda s:s['objective'].update(type='maximize_cooling_power',unit='W'))
    definition=compile_study(load_study(study))
    from types import SimpleNamespace
    from dada_solver.performance import OperatingMode
    value=definition.objective.evaluate(SimpleNamespace(performance=SimpleNamespace(operating_mode=OperatingMode.REFRIGERATION,cooling_power=5.)))
    assert value.value==-5. and value.available
    structured=initialize_external_stream(tmp_path/'structured.toml','structured_c2_15p','structured_c2_15p')
    assert load_study(structured).settings['small']['family']=='structured_c2_15p'


def test_table_refrigeration_cycle_parity(refrigerator,tmp_path):
    source,definition,reference=refrigerator
    path=initialize_external_stream(tmp_path/'tabulated.toml')
    table=ideal_validation_table(definition.configuration.gas,rho_axis=(.01,.1,.5,1.,2.,5.,20.,50.))
    change_basis(path,lambda b:b['configuration'].update(gas=table.to_data()))
    d=compile_study(load_study(path))
    result=MachineEvaluator(d).evaluate(candidate_for_values(d,{}))
    assert result['status']=='feasible'
    assert result['rhs_backend']['actual_backend']=='numba_tabulated'
    for key in ('cooling_power_w','heating_power_w','cooling_cop','indicated_mechanical_input_power_w'):
        assert result['metrics'][key]==pytest.approx(reference['metrics'][key],rel=2e-6,abs=1e-7)


def test_v3_sobol_resume_snapshot_and_offline_report(study,tmp_path,refrigerator):
    from tests.test_research_v2 import Clock,Evaluator
    def activate(raw):
        row=next(r for r in raw['parameters'] if r['name']=='external_stream.heat_in.wall_conductance_w_k')
        value=row.pop('value');row.update(initial=value,lower=140.,upper=160.,kind='continuous',transform='linear')
    rewrite(study,activate)
    d=compile_study(load_study(study));clock=Clock()
    run=OptimizationCampaign(d,tmp_path/'run',evaluator=Evaluator(clock),clock=clock)
    run.run(100,maximum_candidates=1)
    resumed=OptimizationCampaign.resume(tmp_path/'run',evaluator=Evaluator(clock),clock=clock)
    resumed.run(100,maximum_candidates=1)
    records=resumed.history.load()
    assert [r['sequence_index'] for r in records]==[0,1]
    assert len({r['candidate_id'] for r in records})==2
    assert CampaignDefinition.resume(tmp_path/'run').study.data['schema_version']==3
    data=report.inspect(tmp_path/'run');html=report.render_html(data,tmp_path/'report.html').read_text()
    for label in ('Current value','Relative margin','external_streams','capacity_rate_w_k','Cooling power'):
        assert label in html
    assert data['runtime_compatible']


def test_initial_fluid_domain_rejection_is_a_candidate_result(study):
    d=compile_study(load_study(study))
    table=ideal_validation_table(d.configuration.gas,rho_axis=(10.,20.,30.))
    change_basis(study,lambda b:b['configuration'].update(gas=table.to_data()))
    new=compile_study(load_study(study))
    r=MachineEvaluator(new).evaluate(candidate_for_values(new,{}))
    assert r['status']=='invalid_fluid_domain' and not r['integrated']


def test_continuous_diode_wall_events_are_explicitly_unavailable(refrigerator):
    _, definition, result = refrigerator
    assert definition.configuration.valve_model == 'continuous_ideal_diode'
    diagnostics = result['diagnostics']
    assert diagnostics['valve_events'] == []
    assert diagnostics['topology']['classification'] == 'unavailable'
    assert diagnostics['topology']['reasons'] == ['wall_integrator_does_not_record_valve_events']
    assert 'detected' in result['derived']['local_reflux']
