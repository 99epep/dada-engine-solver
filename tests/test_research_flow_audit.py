"""Optional component flow caps, retained physical guards and offline evidence."""
import copy
from dataclasses import dataclass, replace
import json
from types import SimpleNamespace as NS
import tomllib
import pytest
from dada_solver.research.cli import initialize
from dada_solver.research.presets import initialize_v2, initialize_v3
from dada_solver.research.schema import load_study, compile_study
from dada_solver.research.study_io import dumps
from dada_solver.research.margins import limiting_evidence, margin_record
from dada_solver.research.cockpit import campaign_evidence
from dada_solver.research.report import render_html
from dada_solver.sizing.constraints import MaximumAbsoluteMassFlow
from dada_solver.campaign.evaluator import MachineEvaluator
from dada_solver.exchangers.gas_correlations import MicrotubeGasModel, MicrotubeDomainError
from dada_solver.exchangers.microtube_geometry import MicrotubeBank
from dada_solver.validity import ValidityVerdict


@pytest.mark.parametrize('initializer',[initialize_v2,initialize_v3])
def test_new_presets_have_no_flow_cap_but_accept_explicit_one(tmp_path,initializer):
    path=initializer(tmp_path/'study.toml');original=load_study(path)
    assert not any(c['type']=='maximum_absolute_mass_flow' for c in original.data['constraints'])
    raw=tomllib.loads(path.read_text());raw['constraints'].append(dict(type='maximum_absolute_mass_flow',limit=.08,unit='kg/s'))
    path.write_text(dumps(raw));explicit=load_study(path)
    assert explicit.study_id!=original.study_id
    assert any(isinstance(c,MaximumAbsoluteMassFlow) and c.limit==.08 for c in compile_study(explicit).constraints)


def test_assessment_exceeds_008_without_optional_cap_and_preserves_rejections(tmp_path):
    definition=compile_study(load_study(initialize_v2(tmp_path/'study.toml')))
    @dataclass
    class Conservation: residual:float=0.
    @dataclass
    class Validity: verdict:ValidityVerdict=ValidityVerdict.VALID
    performance=NS(gas_power=50.,motor_power=50.,heat_in_power=500.,heat_out_power=-450.,thermal_efficiency=.1,
        conservation=Conservation(),cooling_power=None,cooling_cop=None,mechanical_input_power=None)
    evaluation=NS(usable=True,performance=performance,validity=Validity(),
        diagnostics=NS(pressure_extrema={'S':NS(maximum=2e5)},temperature_extrema={'S':NS(maximum=400)},
            mass_flow_extrema={'port':NS(minimum=-.12,maximum=.1)}))
    def assess(extra=()): return MachineEvaluator(definition)._assessment(evaluation,{},None,None,{'cycles_completed':1},None,extra)
    assert assess()['status']=='feasible'
    cap=MaximumAbsoluteMassFlow(.08);definition.constraints+= (cap,)
    assert assess()['status']=='converged_infeasible'
    assert cap.evaluate(evaluation).margin==pytest.approx(-.04)
    definition.constraints=tuple(c for c in definition.constraints if c is not cap)
    from dada_solver.sizing.constraints import MaximumPressure
    definition.constraints += (MaximumPressure(1.2e6),)
    evaluation.diagnostics.pressure_extrema['S'].maximum=1.3e6
    assert 'maximum_pressure' in assess()['reason']
    evaluation.diagnostics.pressure_extrema['S'].maximum=2e5;evaluation.validity.verdict=ValidityVerdict.INVALID
    assert 'valid_thermodynamic_model' in assess()['reason']
    evaluation.validity.verdict=ValidityVerdict.VALID
    assert 'microtube_model_domain' in assess([dict(name='microtube_model_domain',margin=-1,available=True,satisfied=False)])['reason']


def test_real_geometry_can_accept_high_flow_without_relaxing_domains():
    model=MicrotubeGasModel()
    bank=MicrotubeBank(50000,.8,.00033,.0001524,.00125,.001)
    valid=model.diagnose(bank,.12,2e5,1.99e5,300)
    model.require(valid)
    assert valid.model_validity=='valid'
    drop=model.diagnose(bank,.12,2e5,1e5,300)
    assert 'large_relative_pressure_drop' in drop.issues
    with pytest.raises(MicrotubeDomainError): model.require(drop)
    narrow=replace(bank,tube_count=100)
    mach=model.diagnose(narrow,.12,2e5,1.99e5,300)
    assert 'high_mach' in mach.issues
    with pytest.raises(MicrotubeDomainError):model.require(mach)
    with pytest.raises(ValueError):model.diagnose(bank,.12,2e5,1.99e5,199)


def test_suggestions_have_specific_offline_actions():
    scientific=dict(parameters=[],constraints=[],basis={})
    records=[]
    for i,(status,reason,constraint) in enumerate([
        ('converged_infeasible','mass flow exceeded','maximum_absolute_mass_flow'),
        ('invalid_exchanger','hydrodynamic_entry_unresolved',None),
        ('integration_failure','Transport temperature 199.8 K outside declared domain',None),
        ('converged_infeasible','invalid model','valid_thermodynamic_model')]):
        records.append(dict(candidate_id=str(i),status=status,reason=reason,integrated=True,converged=False,
            objective={'available':False},physical={},normalized=[],constraints=[] if constraint is None else
            [dict(name=constraint,available=True,satisfied=False)]))
    evidence=campaign_evidence(records,scientific,'my campaign',True)
    suggestions=evidence['suggestions']
    assert sum('report' in s.get('command','') for s in suggestions)==1
    actions=[s['action'] for s in suggestions if 'action' in s]
    assert len({(a['field'],a['value']) for a in actions})==5
    assert any(a['value']=='maximum_absolute_mass_flow' for a in actions)
    assert any(a['value']=='hydrodynamic_entry_domain' for a in actions)
    assert all(a['source']=='my campaign' for a in actions)


def test_limiting_evidence_uses_stored_limits_and_keeps_rejected_state_context():
    scientific={'basis':{'heat_in':{'inputs':{'gas_model':{'maximum_mach':.3,'maximum_relative_pressure_drop':.2}}}}}
    record=dict(status='invalid_exchanger',constraints=[margin_record('pressure',100,200,'maximum','Pa')],
        diagnostics={'first_microtube_failure':dict(criterion='high_mach',exchanger='heat_in',passage='inlet',
            angle_rad=.2,time_s=.1,mach=.4,relative_pressure_drop=.19)},derived={})
    before=copy.deepcopy(record)
    evidence=limiting_evidence(record,scientific)
    assert record==before
    mach=next(c for c in evidence if c['name']=='microtube.mach')
    assert mach['margin']==pytest.approx(-.1) and mach['state']=='violated'
    assert 'angle=0.2' in mach['context'] and mach['method']=='first rejected trial state'
    pressure=next(c for c in evidence if c['name']=='microtube.relative_pressure_drop')
    assert pressure['limit']==.2
    unknown=limiting_evidence(record,{'basis':{}})
    assert next(c for c in unknown if c['name']=='microtube.mach')['limit'] is None
    assert not any(c['name']=='maximum_absolute_mass_flow' for c in evidence)


def test_cycle_domain_evidence_distinguishes_sampled_extrema():
    scientific={'basis':{'heat_out':{'inputs':{'gas_model':{'maximum_mach':.3}}}}}
    record=dict(status='feasible',constraints=[],derived={'microtube_gas_domains':{'passages':{
        'Ho.inlet':{'ranges':{'mach':{'maximum':.299}}}}}})
    evidence=limiting_evidence(record,scientific)
    assert evidence[0]['near_active'] and evidence[0]['context']=='Ho.inlet'
    assert 'sampled cycle' in evidence[0]['method']


def test_report_defaults_follow_existing_ranking_and_table_retains_data(tmp_path):
    from tests.test_research_cockpit import campaign
    from dada_solver.research.report import compare
    c,_=campaign(tmp_path,4);data=compare([c.directory])
    original=list(data['selected'])
    # The last record is the objective winner; ordering in the journal is retained.
    for i,r in enumerate(data['selected']):r['objective'].update(value=-i,available=True)
    html=render_html(data,tmp_path/'report.html').read_text()
    raw=html.split('<script id="data" type="application/json">')[1].split('</script>')[0]
    embedded=json.loads(raw)
    assert embedded['comparison_default_ids']==[original[-1]['candidate_id'],original[-2]['candidate_id']]
    assert len(embedded['selected'])==4
    for token in ('max-height:640px','candidateViewport',"'Basin'",'matchesInspection','Inspect matching candidates','Signed margin','Near boundary','Scope / context'):
        assert token in html
    assert 'tr.cells[idColumn]' in html
    assert 'select.add(new Option(action.value,action.value))' in html
    assert 'No matching records in this embedded selection' in html
    data['comparison_compatible']=False
    html=render_html(data,tmp_path/'different.html').read_text()
    assert '"comparison_default_ids": []' in html


def test_obsolete_generic_mach_and_isothermality_are_not_model_boundaries():
    scientific={'basis':{'configuration':{'validity':{'maximum_mach_number':.2}}}}
    record=dict(status='converged_infeasible',constraints=[],metrics={'validity':{'maximum_mach_number':.25,
        'cold_isothermality_error':.8}})
    evidence=limiting_evidence(record,scientific)
    assert evidence==[]


@pytest.mark.parametrize('status,reason',[
    ('feasible',''),('converged_infeasible','maximum_pressure: violated'),
    ('periodic_non_convergence','periodic tolerance not reached'),
    ('budget_exhausted','deadline'),('integration_failure','other failure'),
    ('invalid_exchanger','MicrotubeDomainError: large_relative_pressure_drop'),
])
def test_historical_trial_is_not_a_final_boundary(tmp_path,status,reason):
    scientific={'basis':{'heat_out':{'inputs':{'gas_model':{'maximum_mach':.3}}}}}
    failure=dict(state_kind='rejected_trial_state',criterion='high_mach',mach=.99,
        exchanger='heat_out',passage='hot_to_small',angle_rad=0,time_s=0.)
    record=dict(status=status,reason=reason,safe_retry_used=True,constraints=[],
        diagnostics={'first_microtube_failure':failure},derived={'microtube_gas_domains':{
            'passages':{'Ho.outlet':{'ranges':{'mach':{'maximum':.007}},
                'hydraulic_upstream_ranges':{'mach':{'maximum':.006}}}}}})
    before=copy.deepcopy(record)
    evidence=limiting_evidence(record,scientific)
    assert record==before
    assert next(c for c in evidence if c['name']=='microtube.mach')['value']==.007
    assert all(c['satisfied'] for c in evidence)
    assert not any(c['value']==.99 or c['method']=='first rejected trial state' for c in evidence)


def test_matching_final_domain_rejection_keeps_first_trial_evidence():
    record=dict(status='invalid_exchanger',reason='MicrotubeDomainError: high_mach',safe_retry_used=True,
        diagnostics={'first_microtube_failure':dict(criterion='high_mach',mach=.99,exchanger='heat_out')})
    scientific={'basis':{'heat_out':{'inputs':{'gas_model':{'maximum_mach':.3}}}}}
    row=next(r for r in limiting_evidence(record,scientific) if r['name']=='microtube.mach')
    assert row['value']==.99 and row['state']=='violated'
    assert row['method']=='first rejected trial state'


def test_recovered_trial_remains_in_html_history_not_current_boundaries(tmp_path):
    from tests.test_research_cockpit import campaign
    from dada_solver.research.report import compare
    c,_=campaign(tmp_path,1);data=compare([c.directory]);record=data['selected'][0]
    record.update(status='feasible',safe_retry_used=True,diagnostics={'first_microtube_failure':dict(
        state_kind='rejected_trial_state',criterion='high_mach',mach=.99,exchanger='heat_out')})
    record['derived']={'microtube_gas_domains':{'passages':{'Ho.inlet':{'ranges':{'mach':{'maximum':.007}}}}}}
    data['scientific']['basis']={'heat_out':{'inputs':{'gas_model':{'maximum_mach':.3}}}}
    record['limiting_evidence']=limiting_evidence(record,data['scientific'])
    html=render_html(data,tmp_path/'recovered.html').read_text()
    embedded=json.loads(html.split('<script id="data" type="application/json">')[1].split('</script>')[0])['selected'][0]
    assert embedded['diagnostics']['first_microtube_failure']['mach']==.99
    mach=[c for c in embedded['limiting_evidence'] if c['name']=='microtube.mach']
    assert len(mach)==1 and mach[0]['value']==.007 and mach[0]['satisfied']
    assert 'Rejected trial history' in html and 'Recovered trial diagnostic' in html
    assert 'candidate.diagnostics?.first_microtube_failure' in html
