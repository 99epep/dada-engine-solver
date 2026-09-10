"""Fast campaign tests use a deterministic evaluator and a simulated clock."""
from dataclasses import asdict
from pathlib import Path
import json
import subprocess
import sys
import numpy as np
import pytest
from dada_solver.campaign.parameters import ContinuousParameter, ParameterSpace
from dada_solver.campaign.candidate import Candidate, content_hash
from dada_solver.campaign.strategy import SobolStrategy
from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.campaign.runner import OptimizationCampaign, parse_budget, estimated_next_seconds
from dada_solver.campaign.history import CampaignHistory, atomic_json
from dada_solver.campaign.report import make_report, elite_records
from dada_solver.campaign.evaluator import rejected, MachineEvaluator
from dada_solver.campaign.evaluator import (EvaluationControl, rescale_wall_state,
    select_warm_start)
from dada_solver.campaign.adapters import validate_ownership, microtube_design
from dada_solver.free_kinematics import FreeKinematics

ROOT = Path(__file__).resolve().parents[1]


class Clock:
    value = 0.
    def __call__(self): return self.value
    def advance(self, value): self.value += value


class Evaluator:
    def __init__(self, clock, values=None):
        self.clock, self.values, self.calls = clock, values, []
    def evaluate(self, candidate):
        i = len(self.calls); self.calls.append(candidate.candidate_id); self.clock.advance(2)
        value = self.values[i] if self.values is not None else -sum(candidate.payload['normalized'])
        result = rejected('feasible', '')
        result.update(integrated=True, converged=True, objective=dict(name='fixture', value=value, available=True),
            constraints=[dict(name='fixture_margin', margin=1., satisfied=True, available=True)],
            metrics={'indicated_power_w':-value}, periodic_cycle_count=3,
            final_periodic_state={'values':[1.,2.], 'reusable_only_as_initial_guess':True})
        return result


@pytest.fixture
def definition():
    return CampaignDefinition(ROOT/'examples/free_kinematics_campaign.toml')


@pytest.mark.parametrize('transform,bounds,expected', [('linear',(-2.,8.),3.),('log',(1.,100.),10.)])
def test_parameter_transforms(transform,bounds,expected):
    p = ContinuousParameter('x',*bounds,expected,transform)
    assert p.decode(.5) == pytest.approx(expected)
    assert p.encode(expected) == pytest.approx(.5)
    for u in (0.,.1,.9,1.): assert p.encode(p.decode(u)) == pytest.approx(u, abs=1e-14)
    assert p.decode(0) == bounds[0] and p.decode(1) == bounds[1]


@pytest.mark.parametrize('coordinates', [(-.01,), (1.01,), (float('nan'),), (0.,1.)])
def test_invalid_normalized_values(coordinates):
    with pytest.raises(ValueError): ParameterSpace((ContinuousParameter('x',1,2,1.5),)).decode(coordinates)


@pytest.mark.parametrize('args', [('x',0,2,1,'log'), ('x',2,1,1,'linear'), ('x',1,2,3,'linear')])
def test_invalid_parameter_definitions(args):
    with pytest.raises(ValueError): ContinuousParameter(*args)


def test_candidate_hash_portable_across_processes(definition):
    c = Candidate.create(definition.space, definition.space.initial_coordinates,
        families=definition.families,numerical_settings=definition.numerical_settings,definition_id=definition.definition_id)
    result = subprocess.check_output([sys.executable,'-c',
        'import hashlib,sys; print(hashlib.sha256(sys.stdin.buffer.read()).hexdigest())'],input=c.payload_json.encode())
    assert result.decode().strip() == c.candidate_id
    assert content_hash(c.payload) == c.candidate_id
    changed = dict(c.payload, numerical_settings={'different':True})
    assert content_hash(changed) != c.candidate_id
    payload = c.payload; payload['physical'].clear()
    assert c.payload['physical']


@pytest.mark.parametrize('scramble', [False,True])
def test_sobol_continuation(scramble):
    whole = SobolStrategy(4, seed=19, scramble=scramble)
    expected = [whole.next_point() for _ in range(13)]
    first = SobolStrategy(4,seed=19,scramble=scramble)
    prefix = [first.next_point() for _ in range(5)]
    state = first.state(); state.pop('type')
    resumed = SobolStrategy(**state)
    assert prefix+[resumed.next_point() for _ in range(8)] == expected
    assert resumed.index == 13


def test_resume_history_best_archive_and_phase_delta(tmp_path,definition):
    clock=Clock(); evaluator=Evaluator(clock,[-1.,-3.,-2.])
    campaign=OptimizationCampaign(definition,tmp_path,evaluator=evaluator,clock=clock)
    first=campaign.run(100,maximum_candidates=3)
    assert first['best_at_phase_end']['objective']['value'] == -3
    assert first['best_in_phase'][0]['objective']['value'] == -3
    assert first['attempted']==3 and first['integrated']==3
    assert first['best_at_phase_start'] is None
    old_bytes=(tmp_path/'history.jsonl').read_bytes()
    later=Evaluator(clock,[-4.,-2.])
    second=OptimizationCampaign.resume(tmp_path,evaluator=later,clock=clock).run(100,maximum_candidates=2)
    assert second['best_at_phase_start']['objective']['value']==-3
    assert second['best_at_phase_end']['objective']['value']==-4
    assert second['objective_improvement']==1
    assert second['next_sequence_index']==5
    assert set(later.calls).isdisjoint(evaluator.calls)
    assert (tmp_path/'history.jsonl').read_bytes().startswith(old_bytes)
    records=CampaignHistory(tmp_path).load()
    assert len(records)==5 and len({r['candidate_id'] for r in records})==5
    assert all(r['final_periodic_state'] for r in records)
    assert (tmp_path/'reports/phase_0002.txt').exists()
    delta=second['parameter_changes'][0]
    i=[p.name for p in definition.space.parameters].index(delta['name'])
    assert delta['normalized_delta']==second['best_at_phase_end']['normalized'][i]-second['best_at_phase_start']['normalized'][i]
    assert delta['physical_start']==second['best_at_phase_start']['physical'][delta['name']]


def test_duplicate_cache_uses_preexisting_exact_candidate(tmp_path,definition):
    clock=Clock(); evaluator=Evaluator(clock)
    campaign=OptimizationCampaign(definition,tmp_path,evaluator=evaluator,clock=clock)
    campaign.run(100,maximum_candidates=1)
    history=campaign.history; record=history.load()[0]
    strategy=SobolStrategy(len(definition.space.parameters),seed=definition.seed,scramble=definition.scramble,index=1)
    candidate=Candidate.create(definition.space,strategy.next_point(),families=definition.families,
        numerical_settings=definition.numerical_settings,definition_id=definition.definition_id)
    # Simulate an already completed external/initial candidate matching the next Sobol point.
    record.update(candidate.payload);record['candidate_id']=candidate.candidate_id
    # Isolated fixture journal: this is not a production history rewrite.
    (tmp_path/'history.jsonl').write_text(json.dumps(record)+'\n')
    for p in (tmp_path/'candidates').glob('*.json'): p.unlink()
    atomic_json(tmp_path/'candidates'/f'{candidate.candidate_id}.json',record)
    resumed=Evaluator(clock)
    result=OptimizationCampaign.resume(tmp_path,evaluator=resumed,clock=clock).run(100,maximum_candidates=1)
    assert not resumed.calls
    assert result['cache_hits']==1 and result['integrated']==0
    assert CampaignHistory(tmp_path).load()[-1]['objective']==record['objective']


def test_preflight_rejection_is_persisted(tmp_path,definition):
    class Invalid:
        def evaluate(self,candidate): return rejected('invalid_exchanger','tube pitch is too small')
    report=OptimizationCampaign(definition,tmp_path,evaluator=Invalid()).run(20,maximum_candidates=1)
    assert report['rejected_before_integration']==1
    assert report['integrated']==0
    record=CampaignHistory(tmp_path).load()[0]
    assert record['reason']=='tube pitch is too small' and record['objective'] is None
    assert report['best_at_phase_end'] is None


def test_real_free_preflight_and_injection(definition):
    design=definition.adapter.build(definition.space.decode(definition.space.initial_coordinates))
    assert isinstance(design.kinematics,FreeKinematics)
    model=design.build()
    assert model.kinematics.forward is design.kinematics
    settings=dict(definition.free_settings)
    settings['small_limits']={'maximum_absolute_first_derivative':0.}
    from dada_solver.campaign.adapters import FamilyDesignAdapter
    definition.adapter=FamilyDesignAdapter(definition.configuration,definition.families,settings)
    candidate=Candidate.create(definition.space,definition.space.initial_coordinates,families=definition.families,
        numerical_settings=definition.numerical_settings,definition_id=definition.definition_id)
    result=MachineEvaluator(definition).evaluate(candidate)
    assert result['status']=='invalid_kinematics' and not result['integrated']
    assert result['preflight_diagnostics'][0]['first_derivative_margin']<0


def test_preflight_rejects_invalid_parameterization(definition):
    from dada_solver.campaign.adapters import PreflightRejection
    with pytest.raises(PreflightRejection) as error:
        definition.adapter.build({'common.small_clearance_volume':-1})
    assert error.value.status=='invalid_parameterization'


@pytest.mark.parametrize('name',['common.cold_ua','common.hot_heat_exchanger_volume','common.small_to_cold_cda'])
def test_microtube_derived_quantity_ownership(name):
    with pytest.raises(ValueError,match='owns derived'):
        validate_ownership({name},{'kinematics':'free','exchanger':'microtube'}, {})


def test_wrong_family_parameters_rejected():
    with pytest.raises(ValueError,match='wrong-family'):
        validate_ownership({'four_bar.crank_ratio'},{'kinematics':'free','exchanger':'reservoir'}, {})


def test_budget_uses_robust_recent_integrations_without_sleep(tmp_path,definition):
    clock=Clock(); evaluator=Evaluator(clock)
    definition.initial_evaluation_seconds=1
    report=OptimizationCampaign(definition,tmp_path,evaluator=evaluator,clock=clock).run(5,maximum_candidates=100)
    assert report['attempted']==2 and report['actual_duration_seconds']==4
    assert report['stopping_reason']=='wall_clock_budget'
    assert estimated_next_seconds([{'integrated':True,'duration_seconds':2}],1)==pytest.approx(2.2)
    assert estimated_next_seconds([{'integrated':False,'duration_seconds':.01}],10)==10
    zero=OptimizationCampaign.resume(tmp_path,evaluator=evaluator,clock=clock).run(0)
    assert zero['attempted']==0


@pytest.mark.parametrize('text,value',[('30m',1800),('1h30m',5400),('1h 2m 3s',3723),('0.5h',1800),('10',10)])
def test_budget_parsing(text,value):
    assert parse_budget(text)==value


def test_archive_excludes_infeasible_even_with_better_score():
    records=[dict(candidate_id=str(i),status=s,objective=dict(available=True,value=v))
             for i,s,v in [(0,'feasible',-1),(1,'feasible',-2),(2,'converged_infeasible',-100),(3,'feasible',-1.5)]]
    assert [r['candidate_id'] for r in elite_records(records)]==['1','3','0']


def test_report_bound_pressure_and_suggestions():
    space=ParameterSpace((ContinuousParameter('x',1,100,10,'log'),))
    def rec(i,u,value):
        return dict(candidate_id=str(i),normalized=[u],physical=space.decode([u]),status='feasible',
            objective=dict(available=True,value=value),constraints=[dict(name='margin',margin=1,available=True,satisfied=True)],
            integrated=True,converged=True,duration_seconds=1,metrics={})
    report=make_report(space,[rec(0,.5,10)],[rec(i,.96+i*.001,10-i) for i in range(1,5)],
        requested_seconds=10,elapsed_seconds=4,elite_size=4)
    assert report['bound_pressure'][0]['side']=='upper'
    assert report['parameter_changes'][0]['normalized_delta']==pytest.approx(.464)
    assert any('bound of x' in x for x in report['suggestions'])
    assert any('diameter' in x for x in report['suggestions'])


def test_orphan_completed_record_and_torn_tail_recover_without_rerun(tmp_path,definition):
    clock=Clock(); ev=Evaluator(clock)
    campaign=OptimizationCampaign(definition,tmp_path,evaluator=ev,clock=clock)
    campaign.run(100,maximum_candidates=1)
    original=(tmp_path/'history.jsonl').read_bytes()
    # Crash after atomic candidate save but during the journal append.
    (tmp_path/'history.jsonl').write_bytes(original[:25])
    recovered=CampaignHistory(tmp_path).load()
    assert len(recovered)==1
    assert (tmp_path/'history.jsonl').read_bytes()==original
    assert list(tmp_path.glob('history_torn_tail_*.bin'))
    resumed=Evaluator(clock)
    result=OptimizationCampaign.resume(tmp_path,evaluator=resumed,clock=clock).run(100,maximum_candidates=1)
    assert result['next_sequence_index']==2
    assert set(ev.calls).isdisjoint(resumed.calls)


def test_crash_in_flight_retries_only_unfinished_candidate(tmp_path,definition):
    clock=Clock()
    class Crash:
        def evaluate(self,candidate):
            self.candidate=candidate.candidate_id
            raise KeyboardInterrupt()
    evaluator=Crash();campaign=OptimizationCampaign(definition,tmp_path,evaluator=evaluator,clock=clock)
    with pytest.raises(KeyboardInterrupt):campaign.run(100,maximum_candidates=1)
    assert campaign.history.state()['pending']['candidate_id']==evaluator.candidate
    resumed=Evaluator(clock)
    result=OptimizationCampaign.resume(tmp_path,evaluator=resumed,clock=clock).run(100,maximum_candidates=1)
    assert resumed.calls==[evaluator.candidate] and result['next_sequence_index']==1


def test_definition_changes_are_not_silently_resumed(tmp_path,definition):
    campaign=OptimizationCampaign(definition,tmp_path,evaluator=Evaluator(Clock()))
    path=tmp_path/'base.toml';path.write_text(path.read_text()+'\n# Changed snapshot\n')
    with pytest.raises(ValueError,match='changed'): OptimizationCampaign.resume(tmp_path)


def test_microtube_adapter_keeps_geometry_derived():
    from tests.test_exchanger_hardware import inputs
    from dada_solver.campaign.adapters import PreflightRejection
    bank, props = inputs()
    design = microtube_design(bank, props, 1e-4,
        {'microtube.heat_in.tube_length_m':bank.tube_length_m*2},'heat_in')
    assert design.bank.tube_length_m==bank.tube_length_m*2
    assert design.build().gas_volume_m3 > bank.dimensions()['working_gas_volume_m3']
    with pytest.raises(PreflightRejection) as error:
        microtube_design(bank,props,1e-4,{'microtube.heat_in.pitch_m':1e-8},'heat_in')
    assert error.value.status=='invalid_exchanger'


def test_four_bar_family_adapter(definition):
    from dada_solver.campaign.adapters import FamilyDesignAdapter, PreflightRejection
    from dada_solver.four_bar import FourBarKinematics
    adapter=FamilyDesignAdapter(definition.configuration,{'kinematics':'four_bar','exchanger':'reservoir'}, {})
    design=adapter.build({'four_bar.crank_ratio':.4})
    assert isinstance(design.kinematics,FourBarKinematics)
    with pytest.raises(PreflightRejection) as error:adapter.build({'four_bar.crank_ratio':-1})
    assert error.value.status=='invalid_kinematics'


@pytest.mark.parametrize('status,usable,expected',[
    ('NUMERICAL_FAILURE',False,'integration_failure'),
    ('NOT_CONVERGED',False,'periodic_non_convergence'),
    ('CONVERGED',True,'converged_infeasible')])
def test_physical_evaluator_statuses_remain_separate(monkeypatch,definition,status,usable,expected):
    from types import SimpleNamespace
    from dada_solver.sizing.evaluator import EvaluationStatus
    result=SimpleNamespace(usable=usable,status=getattr(EvaluationStatus,status),
        periodic=SimpleNamespace(history=(),message='fixture result'),performance=None,diagnostics=None,validity=None)
    monkeypatch.setattr('dada_solver.campaign.evaluator.evaluate_configuration',lambda *a,**kw:result)
    candidate=Candidate.create(definition.space,definition.space.initial_coordinates,families=definition.families,
        numerical_settings=definition.numerical_settings,definition_id=definition.definition_id)
    record=MachineEvaluator(definition).evaluate(candidate)
    assert record['status']==expected and record['integrated']
    assert record['converged']==usable
    assert all(not c['available'] for c in record['constraints'])


def test_candidate_tampering_is_rejected(definition):
    c=Candidate.create(definition.space,definition.space.initial_coordinates,families=definition.families,
        numerical_settings=definition.numerical_settings,definition_id=definition.definition_id)
    with pytest.raises(ValueError,match='identity'):Candidate(c.payload_json,'wrong')


def test_unavailable_constraints_are_not_reported_as_violations():
    space=ParameterSpace((ContinuousParameter('x',0,1,.5),))
    record=dict(candidate_id='x',normalized=[.5],physical={'x':.5},status='periodic_non_convergence',
        objective=dict(available=False,value=None),constraints=[dict(name='power',margin=None,available=False,satisfied=False)],
        integrated=True,converged=False,duration_seconds=1,metrics={},reason='maximum cycles',
        periodic_convergence=dict(improving=True))
    report=make_report(space,[],[record],requested_seconds=2,elapsed_seconds=1)
    assert report['violated_constraints']=={}
    assert report['unavailable_constraints']=={'power':1}
    assert 'power' not in ' '.join(report['suggestions'])
    assert 'converge' in ' '.join(report['suggestions'])


def test_microtube_campaign_builds_dynamic_wall_and_snapshots_hardware(tmp_path):
    definition=CampaignDefinition(ROOT/'examples/microtube_free_campaign.toml')
    physical=definition.space.decode(definition.space.initial_coordinates)
    design=definition.adapter.build(physical)
    from dada_solver.exchangers.air_wall import AirWallMotor
    assert isinstance(design.build(),AirWallMotor)
    changed=dict(physical);changed['microtube.heat_in.tube_length_m']*=1.05
    assert definition.adapter.build(changed).heat_in.build().gas_volume_m3 != design.heat_in.build().gas_volume_m3
    OptimizationCampaign(definition,tmp_path,evaluator=Evaluator(Clock()))
    assert (tmp_path/'hardware.toml').read_text()==definition.hardware_source
    assert json.loads((tmp_path/'definition.json').read_text())['hardware_configuration']==definition.hardware_source


def test_wall_energy_rescaling_preserves_specific_energy_and_temperature():
    old=np.array([1.,10.,2.,30.,3.,60.,4.,100.,600.,1200.])
    new=rescale_wall_state(old,20.,[2.,4.],[3.,2.])
    assert new[:8:2].sum()==pytest.approx(20.)
    assert new[1:8:2]/new[:8:2] == pytest.approx(old[1:8:2]/old[:8:2])
    assert new[8:10]/[3.,2.] == pytest.approx(old[8:10]/[2.,4.])


def test_warm_start_prefers_converged_then_nearest_compatible():
    def record(cid,u,status,converged=True,family='wall'):
        return dict(candidate_id=cid,normalized=[u],status=status,converged=converged,
            initial_guess_state=dict(state_layout='layout',evaluator_family=family,
                operating_direction='motor',values=[1]*10))
    records=[record('near',.49,'periodic_non_convergence',False),record('far',.2,'converged_infeasible'),record('bad',.5,'feasible',True,'other')]
    selected,distance=select_warm_start(records,[.5],'layout','wall','motor')
    assert selected['candidate_id']=='far' and distance==pytest.approx(.3)


def test_budget_exhausted_result_is_retried_not_cached(tmp_path,definition):
    class Controlled:
        def __init__(self): self.calls=0
        def evaluate_with_control(self,candidate,control):
            self.calls+=1
            result=rejected('budget_exhausted' if self.calls==1 else 'feasible','deadline')
            result.update(integrated=True,converged=self.calls>1,
                objective=dict(name='fixture',value=0.,available=True),
                constraints=[dict(name='ok',margin=1.,available=True,satisfied=True)])
            return result
    evaluator=Controlled();campaign=OptimizationCampaign(definition,tmp_path,evaluator=evaluator)
    first=campaign.run(100,maximum_candidates=1)
    second=campaign.run(100,maximum_candidates=1,retry_incomplete=True)
    records=CampaignHistory(tmp_path).load()
    assert first['status_counts']=={'budget_exhausted':1}
    assert evaluator.calls==2 and records[0]['candidate_id']==records[1]['candidate_id']
    assert second['cache_hits']==0 and second['next_sequence_index']==1


def test_deadline_control_uses_injected_clock_and_grace(tmp_path,definition):
    clock=Clock();definition.initial_evaluation_seconds=1;definition.deadline_grace_seconds=2
    class DeadlineEvaluator:
        def evaluate_with_control(self,candidate,control):
            assert control.deadline==pytest.approx(5)
            clock.advance(7)
            try: control.check()
            except Exception as error: return rejected('budget_exhausted',str(error)) | {'integrated':True}
    report=OptimizationCampaign(definition,tmp_path,evaluator=DeadlineEvaluator(),clock=clock).run(5,maximum_candidates=1)
    assert report['actual_duration_seconds']==7
    assert report['status_counts']=={'budget_exhausted':1}


def test_wall_solver_retains_only_last_complete_cycle_on_interrupt():
    from dada_solver.exchangers.wall_cycle import solve_periodic_wall_motor
    from dada_solver.integration import IntegrationInterrupted
    class Wrapper:
        heat_in=type('Heat',(),{'wall_capacity_j_k':1})()
        heat_out=type('Heat',(),{'wall_capacity_j_k':1})()
        calls=0
        def integrate_cycle(self,state,**kwargs):
            self.calls+=1
            if self.calls==2: raise IntegrationInterrupted('deadline')
            trajectory=np.zeros((15,2));trajectory[:10,0]=state
            trajectory[:10,1]=np.asarray(state)+1
            return np.array([0.,2*np.pi]),trajectory
    result=solve_periodic_wall_motor(Wrapper(),np.ones(10),maximum_cycles=5)
    assert result.status=='interrupted' and len(result.history)==1
    assert result.last_complete_state==pytest.approx(np.full(10,2.))


def test_changed_hardware_snapshot_refuses_resume(tmp_path):
    definition=CampaignDefinition(ROOT/'examples/microtube_free_campaign.toml')
    OptimizationCampaign(definition,tmp_path,evaluator=Evaluator(Clock()))
    (tmp_path/'hardware.toml').write_text((tmp_path/'hardware.toml').read_text()+'\n# changed\n')
    with pytest.raises(ValueError,match='changed'):
        OptimizationCampaign.resume(tmp_path)
