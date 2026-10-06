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
    text=capsys.readouterr().out
    assert all(candidate in text for candidate in ids)
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
    completed=[e for e in events if e['event']=='evaluation']
    assert len(completed)==22 and sum(e['new_best'] for e in completed)==1
    assert [e['evaluation']['evaluation_number'] for e in completed]==list(range(22))
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


@pytest.mark.parametrize('stage',['recovery_only','journal_and_recovery','torn'])
def test_recovery_windows_no_reintegration(tmp_path,stage):
    c,clock=campaign(tmp_path,1);record=c.history.load()[0];journal=c.history.path
    atomic_json(c.directory/'recovery.json',record)
    if stage=='recovery_only':journal.write_bytes(b'')
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
    assert len(list((c.directory/'candidates').glob('*.json')))==0
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
    assert all(' resume ' not in s.get('command','') for s in evidence['suggestions'])
    assert any('bound' in s['signal'] for s in evidence['suggestions'])
    assert evidence['suggestions'][0]['action']['type']=='bounds'
    assert "'run with spaces'" in evidence['suggestions'][-1]['command']
    sparse=campaign_evidence(rows[:1],scientific,'run',True)
    assert any(s.get('command')=='dada-research resume run --budget 30m' for s in sparse['suggestions'])
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


def test_real_lsoda_callback_exception_is_quiet_and_preserved(capfd):
    error=ValueError('synthetic callback domain rejection')
    def rhs(t,y):raise error
    with pytest.raises(ValueError) as caught:quiet_solve_ivp(rhs,(0,1),[1.],method='LSODA')
    assert caught.value is error
    assert capfd.readouterr().err==''
    os.write(2,b'after LSODA\n')
    assert capfd.readouterr().err=='after LSODA\n'


def test_capture_restores_on_keyboard_interrupt(capfd):
    with pytest.raises(KeyboardInterrupt):
        with capture_native_stderr():
            os.write(2,b'native interrupted\n');raise KeyboardInterrupt
    os.write(2,b'restored after interrupt\n')
    assert capfd.readouterr().err=='restored after interrupt\n'


def test_execution_default_does_not_change_scientific_or_candidate_identity(tmp_path):
    import tomllib
    from dada_solver.research.study_io import dumps
    from dada_solver.research.schema import load_study,compile_study,candidate_for_values
    d=definition(tmp_path);path=tmp_path/'study.toml'
    raw=tomllib.loads(path.read_text());raw['execution']['default_max_candidates']=16
    path.write_text(dumps(raw));changed=compile_study(load_study(path))
    assert changed.study.study_id==d.study.study_id
    values={p.name:p.initial for p in d.space.parameters}
    assert candidate_for_values(changed,values).candidate_id==candidate_for_values(d,values).candidate_id




def test_recovery_mismatch_is_not_appended(tmp_path):
    c,_=campaign(tmp_path,1);record=c.history.load()[0];before=c.history.path.read_bytes()
    record['candidate_id']='bad'
    atomic_json(c.directory/'recovery.json',record)
    with pytest.raises(ValueError,match='identity'):c.history.load()
    assert c.history.path.read_bytes()==before


def test_selected_command_shell_quoting(tmp_path):
    import shlex
    import shutil
    import subprocess
    from dada_solver.research.report import HTML
    if not shutil.which('node'):pytest.skip('Optional JavaScript runtime for browser command quoting check')
    line=next(line for line in HTML.splitlines() if line.startswith('const shellQuote='))
    values=['campaign path',"campaign's path",'$(touch nope)', 'one\ntwo', 'normal/path']
    script=line+'\nconsole.log(JSON.stringify('+json.dumps(values)+'.map(shellQuote)));'
    result=subprocess.check_output(['node','-e',script],text=True)
    assert [shlex.split(value) for value in json.loads(result)]==[[v] for v in values]


def test_inspection_tolerates_live_recovery_acknowledgment(tmp_path,monkeypatch):
    from dada_solver.campaign.history import recovery_records
    path=tmp_path/'recovery.json';path.write_text('{}')
    read=Path.read_text
    def disappeared(self,*args,**kwargs):
        if self==path:raise FileNotFoundError('writer acknowledged recovery')
        return read(self,*args,**kwargs)
    monkeypatch.setattr(Path,'read_text',disappeared)
    assert list(recovery_records(tmp_path))==[]
