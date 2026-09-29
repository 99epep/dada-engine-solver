"""UX and persistence tests use synthetic evaluations, never a search integration."""
import io
import json
import os
from pathlib import Path

import pytest

from dada_solver.campaign.history import CampaignHistory, atomic_json
from dada_solver.campaign.runner import OptimizationCampaign
from dada_solver.research.cli import main
from dada_solver.research import report
from dada_solver.research.cockpit import campaign_evidence, reason_category, unique_prefixes
from dada_solver.research.progress import CLIProgress
from dada_solver.solver_output import capture_native_stderr, quiet_solve_ivp
from tests.test_research_cli import definition
from tests.test_campaign import Clock, Evaluator


def campaign(tmp_path, count=3):
    d=definition(tmp_path);clock=Clock();ev=Evaluator(clock)
    c=OptimizationCampaign(d,tmp_path/'campaign',evaluator=ev,clock=clock)
    c.run(100,maximum_candidates=count)
    return c,clock


def test_report_default_overwrite_compact_and_full_json(tmp_path,capsys):
    c,_=campaign(tmp_path)
    assert main(['report',str(c.directory)])==0
    text=capsys.readouterr().out;path=c.directory/'report.html'
    assert path.exists() and 'HTML:' in text
    ids=[r['candidate_id'] for r in c.history.load()]
    assert sum(candidate in text for candidate in ids)==1
    first=path.read_bytes();main(['report',str(c.directory)])
    assert path.read_bytes()==first
    capsys.readouterr()
    main(['report',str(c.directory),'--json'])
    assert len(json.loads(capsys.readouterr().out)['records'])==3
    main(['status',str(c.directory),'--list-candidates'])
    assert all(candidate in capsys.readouterr().out for candidate in [])  # Read once below.
    main(['report',str(c.directory),'--list-candidates'])
    text=capsys.readouterr().out
    assert all(candidate in text for candidate in ids)


def test_unique_and_ambiguous_prefixes():
    rows=[dict(candidate_id='abcd1234aa'),dict(candidate_id='abcd1234bb')]
    data=dict(records=rows)
    assert report.select_records(data,['abcd1234a'])==rows[:1]
    for selector in ('abcd','unknown'):
        with pytest.raises(ValueError,match='absent or ambiguous'):report.select_records(data,[selector])
    assert unique_prefixes(r['candidate_id'] for r in rows)=={'abcd1234aa':'abcd1234a','abcd1234bb':'abcd1234b'}


def test_real_fd_capture_restores_after_exception(capfd):
    error=RuntimeError('original failure')
    with pytest.raises(RuntimeError) as caught:
        with capture_native_stderr() as captured:
            os.write(2,b'native callback failed\n')
            raise error
    assert caught.value is error
    assert captured['text']=='native callback failed\n'
    os.write(2,b'restored\n')
    assert capfd.readouterr().err=='restored\n'


def test_solver_wrapper_preserves_exception_and_arguments(monkeypatch,capfd):
    error=ValueError('physical domain failure');called=[]
    def fail(*args,**kwargs):
        called.append((args,kwargs));os.write(2,b'capi_return is NULL\n');raise error
    monkeypatch.setattr('scipy.integrate.solve_ivp',fail)
    with pytest.raises(ValueError) as caught:quiet_solve_ivp('rhs',(0,1),[1],rtol=1e-8)
    assert caught.value is error and str(error)=='physical domain failure'
    assert error.native_solver_stderr=='capi_return is NULL\n'
    assert called==[(('rhs',(0,1),[1]),dict(rtol=1e-8))]
    assert capfd.readouterr().err==''


@pytest.mark.parametrize('tty',[True,False])
def test_progress_sink_and_runner_callback(tmp_path,tty):
    class Stream(io.StringIO):
        def isatty(self):return tty
    stream=Stream();sink=CLIProgress(stream);events=[]
    d=definition(tmp_path);clock=Clock()
    ev=Evaluator(clock,[-1.]*22)
    c=OptimizationCampaign(d,tmp_path/'campaign',evaluator=ev,clock=clock)
    c.run(100,maximum_candidates=22,progress_callback=lambda x:(events.append(x),sink(x)))
    assert events[0]['event']=='start' and events[0]['sobol_index']==0
    assert events[-1]['event']=='finish' and events[-1]['attempted']==22
    assert len(events)<12 and any(e['event']=='best' for e in events)
    assert events[-1]['feasible']==22
    assert ('\r' in stream.getvalue())==tty
    assert 'max 22 new candidates' in stream.getvalue()
    assert stream.getvalue().endswith('\n')
    later=[]
    OptimizationCampaign.resume(c.directory,evaluator=Evaluator(clock),clock=clock).run(100,maximum_candidates=1,progress_callback=later.append)
    assert later[0]['sobol_index']==22 and later[-1]['sobol_index']==23


def test_progress_during_long_candidate(tmp_path):
    d=definition(tmp_path);clock=Clock();events=[]
    class Controlled(Evaluator):
        def evaluate_with_control(self,candidate,control):
            for _ in range(3):clock.advance(11);control.check()
            return self.evaluate(candidate)
    c=OptimizationCampaign(d,tmp_path/'campaign',evaluator=Controlled(clock),clock=clock)
    c.run(100,maximum_candidates=1,progress_callback=events.append)
    assert len([e for e in events if e['event']=='progress' and e['attempted']==0])==3


@pytest.mark.parametrize('stage',['recovery_only','journal_and_recovery','torn','legacy'])
def test_recovery_windows_no_reintegration(tmp_path,stage):
    c,clock=campaign(tmp_path,1);record=c.history.load()[0];journal=c.history.path
    if stage=='legacy':
        (c.directory/'candidates').mkdir()
        atomic_json(c.directory/'candidates'/f"{record['candidate_id']}.json",record)
    else:atomic_json(c.directory/'recovery.json',record)
    if stage in ('recovery_only','legacy'):journal.write_text('')
    if stage=='torn':journal.write_bytes(journal.read_bytes()[:30])
    # Read-only inspection sees the completed record, never repairs source files.
    before={p:p.read_bytes() for p in c.directory.rglob('*') if p.is_file()}
    assert len(report.inspect(c.directory)['records'])==1
    assert before=={p:p.read_bytes() for p in c.directory.rglob('*') if p.is_file()}
    # Stale pending state as after a crash between journal and state publication.
    state=c.history.state();state['pending']=dict(sequence_index=0,candidate_id=record['candidate_id'],payload_json='unused')
    state['search']['index']=1;c.history.save_state(state)
    ev=Evaluator(clock)
    resumed=OptimizationCampaign.resume(c.directory,evaluator=ev,clock=clock)
    resumed.run(100,maximum_candidates=1)
    assert record['candidate_id'] not in ev.calls
    assert [r['sequence_index'] for r in resumed.history.load()]==[0,1]
    assert not (c.directory/'recovery.json').exists()
    assert len(list((c.directory/'candidates').glob('*.json')))==(1 if stage=='legacy' else 0)
    assert resumed.history.state()['pending'] is None


def test_failed_state_publication_retains_recovery(tmp_path,monkeypatch):
    d=definition(tmp_path);clock=Clock();c=OptimizationCampaign(d,tmp_path/'c',evaluator=Evaluator(clock),clock=clock)
    original=c.history.save_state
    def fail(state):
        if c.history.path.exists():raise OSError('state write failed')
        original(state)
    monkeypatch.setattr(c.history,'save_state',fail)
    with pytest.raises(OSError):c.run(100,maximum_candidates=1)
    assert (c.directory/'recovery.json').exists()
    ev=Evaluator(clock)
    OptimizationCampaign.resume(c.directory,evaluator=ev,clock=clock).run(0)
    assert not ev.calls and not (c.directory/'recovery.json').exists()


def test_funnel_bounds_categories_and_suggestions(tmp_path):
    scientific=dict(parameters=[dict(name='x',lower=0,upper=1,unit='1')],basis=dict(
        heat_in=dict(inputs=dict(gas_model=dict(transport=dict(minimum_temperature=200,maximum_temperature=1000))))))
    rows=[]
    for i in range(6):
        rows.append(dict(candidate_id=str(i),status='feasible' if i<3 else 'invalid_exchanger',
            integrated=i<3,converged=i<3,physical={'x':.99},normalized=[.99],
            objective=dict(available=True,value=-i),constraints=[],metrics={},
            reason=f'ValueError: Transport temperature {199+i/100} K outside declared domain.'))
    evidence=campaign_evidence(rows,scientific,'run with spaces',True)
    assert evidence['funnel']['attempted']==6 and evidence['funnel']['feasible']==3
    assert evidence['rejection_categories']=={'transport_temperature_below_domain':3}
    assert evidence['bounds'][1]['count']==6 and evidence['bounds'][1]['elite_count']==3
    assert all(' resume ' not in s['command'] for s in evidence['suggestions'])
    assert any('bound' in s['signal'] for s in evidence['suggestions'])
    assert "'run with spaces'" in evidence['suggestions'][0]['command']
    sparse=campaign_evidence(rows[:1],scientific,'run',True)
    assert any(s['command']=='dada-research resume run --budget 30m' for s in sparse['suggestions'])
    high=dict(rows[-1],reason='Transport temperature 1001 K outside declared domain.')
    assert reason_category(high,scientific)=='transport_temperature_above_domain'


def test_html_cockpit_commands_and_no_external_dependencies(tmp_path):
    c,_=campaign(tmp_path)
    data=report.compare([c.directory]);page=report.render_html(data,tmp_path/'cockpit.html').read_text()
    for value in ['Campaign funnel','Pressure on parameter bounds','Suggested next steps','selectedCommands(a,b)','navigator.clipboard','display_candidate_id','--plots']:
        assert value in page
    assert '<script src=' not in page
    assert data['cockpit']['funnel']['attempted']==3


def test_new_and_rescaled_execution_defaults_are_512(tmp_path):
    from dada_solver.research.presets import initialize_v3
    from dada_solver.research.schema import load_study
    from dada_solver.research.rescale import rescale
    from tests.test_research_rescale import snapshot
    path=initialize_v3(tmp_path/'study.toml')
    study=load_study(path)
    assert study.data['execution']['default_max_candidates']==512
    _,_,r=snapshot(path,tmp_path/'run')
    out=rescale(tmp_path/'run',r['candidate_id'][:8],1,tmp_path/'scaled.toml')
    assert load_study(out).data['execution']['default_max_candidates']==512
