"""Presentation tests use fake streams and clocks, never a running real campaign."""
import copy
import io
import pytest
from dada_solver.research.progress import CLIProgress, cells
from dada_solver.campaign.runner import OptimizationCampaign
from tests.test_campaign import Clock, Evaluator
from tests.test_research_cockpit import definition


class Stream(io.StringIO):
    def __init__(self,tty=True): super().__init__();self.tty=tty
    def isatty(self): return self.tty


def event(kind='progress',status='feasible',origin=None):
    metrics=dict(cooling_cop=1.1892,cooling_power_w=286.4,indicated_mechanical_input_power_w=240.8,
                 maximum_absolute_mass_flow_kg_s=.078)
    record=dict(evaluation_number=43,candidate_id='a'*64,status=status,reason='',metrics=metrics,
                constraints=[],search_origin=origin,cache_hit=False)
    return dict(event=kind,phase_id=1,sobol_index=44,attempted=44,maximum_candidates=4096,
        elapsed_seconds=1872,budget_seconds=25200,converged=20,feasible=13,
        best=dict(candidate_id='b'*64,metrics=metrics,objective=dict(value=-1.1892)),
        evaluation=record,new_best=False)


def test_tty_refresh_and_exactly_one_permanent_line():
    stream=Stream();sink=CLIProgress(stream,width=180)
    sink(event());sink(event())
    assert '\n' not in stream.getvalue()
    before=len(stream.getvalue());sink(event('evaluation',origin=dict(region_id='basin_2',kind='sobol')))
    added=stream.getvalue()[before:]
    assert added.count('\n')==1
    assert '[0044] basin_2:sobol  feasible' in added
    assert 'COP 1.1892' in added and 'Qcold 286.4 W' in added and 'Pin 240.8 W' in added and 'mdot 0.078' in added
    assert '44/4096' in added.split('\n')[1]
    assert not added.endswith('\n')


@pytest.mark.parametrize('kind',['center','sobol'])
def test_local_origin_and_new_best(kind):
    stream=Stream(False);sink=CLIProgress(stream)
    e=event('evaluation',origin=dict(region_id='basin_1',kind=kind));e['new_best']=True
    sink(e)
    assert f'basin_1:{kind}' in stream.getvalue()
    assert stream.getvalue().count('BEST')==1 and stream.getvalue().count('\n')==1


def test_violated_constraints_with_declared_limit():
    stream=Stream(False);sink=CLIProgress(stream,scientific={'constraints':[dict(type='maximum_absolute_mass_flow',limit=.08)]})
    e=event('evaluation','converged_infeasible')
    e['evaluation']['constraints']=[dict(name='maximum_absolute_mass_flow',margin=-.004,available=True,satisfied=False),
                                   dict(name='valid_thermodynamic_model',margin=-1.,available=True,satisfied=False)]
    sink(e)
    assert 'maximum_absolute_mass_flow + valid_thermodynamic_model' in stream.getvalue()
    assert 'mdot 0.084 > 0.080' in stream.getvalue()


def test_exchanger_preserves_multiple_scientific_criteria():
    stream=Stream(False);sink=CLIProgress(stream)
    e=event('evaluation','invalid_exchanger');e['evaluation']['reason']='MicrotubeDomainError: large_relative_pressure_drop; high_mach'
    sink(e)
    assert 'large_relative_pressure_drop + high_mach' in stream.getvalue()


def test_integration_failure_domain_uses_existing_cockpit_classification():
    scientific={'basis':{'heat_in':{'inputs':{'gas_model':{'transport':{'minimum_temperature':200,'maximum_temperature':1000}}}}}}
    stream=Stream(False);sink=CLIProgress(stream,scientific=scientific)
    e=event('evaluation','integration_failure');e['evaluation']['reason']='Transport temperature 199.8 K outside declared domain'
    sink(e)
    assert 'transport_temperature_below_domain' in stream.getvalue()


def test_unknown_exception_is_short_sanitized_and_record_unchanged():
    stream=Stream(False);sink=CLIProgress(stream)
    e=event('evaluation','integration_failure');e['evaluation']['reason']='ValueError: useful cause\n'+('detail '*100)+'\x1b[31m'
    original=copy.deepcopy(e);sink(e)
    assert e==original
    assert 'useful cause' in stream.getvalue() and '…' in stream.getvalue()
    assert '\x1b' not in stream.getvalue() and stream.getvalue().count('\n')==1


def test_non_tty_has_no_periodic_lines_or_controls():
    stream=Stream(False);sink=CLIProgress(stream)
    sink(event('start'))
    for _ in range(20): sink(event())
    sink(event('evaluation'));sink(event('finish'));sink.close()
    text=stream.getvalue()
    assert text.count('\n')==3 and '\r' not in text and '\x1b' not in text
    assert '[0044] global  feasible' in text and text.splitlines()[-1].startswith('Finished ·')


@pytest.mark.parametrize('error',[KeyboardInterrupt,RuntimeError,TimeoutError])
def test_exception_cleanup_propagates(error):
    stream=Stream()
    with pytest.raises(error):
        with CLIProgress(stream) as sink:
            sink(event());raise error('fixture')
    assert stream.getvalue().endswith('\r\x1b[K\n')


def test_normal_finish_and_close_are_idempotent():
    stream=Stream()
    with CLIProgress(stream) as sink:
        sink(event());sink(event('finish'))
    assert stream.getvalue().count('\n')==1
    assert stream.getvalue().endswith('\n')
    assert 'Finished ·' in stream.getvalue()


@pytest.mark.parametrize('width',[12,40,60,80,120])
def test_terminal_width_and_priority(width):
    stream=Stream();sink=CLIProgress(stream,width=width);sink(event())
    visible=stream.getvalue().split('\x1b[K')[-1]
    assert cells(visible)<=width-1 and '\n' not in visible
    if width>=40: assert 'feasible 13' in visible
    if width==60: assert 'COP 1.1892' in visible and 'Qcold' not in visible


def test_no_champion_and_motor_metrics():
    stream=Stream(False);sink=CLIProgress(stream)
    e=event('finish');e['best']=None;sink(e)
    assert 'best —' in stream.getvalue()
    e=event('evaluation');e['evaluation']['metrics']={'indicated_thermal_efficiency':.12,'indicated_power_w':150}
    sink(e)
    assert 'eta 0.12' in stream.getvalue() and 'Pgas 150.0 W' in stream.getvalue()


def test_callback_does_not_change_global_or_local_sequence_and_resume(tmp_path):
    d=definition(tmp_path);clock=Clock()
    for local in (False,True):
        if local:
            names=[p.name for p in d.space.parameters]
            d.search_settings=dict(domain='local_regions_v1',seed=d.seed,scramble=d.scramble,radius_fraction=.2,
                regions=[dict(id=f'basin_{i+1}',center=d.space.decode([z]*len(names))) for i,z in enumerate((.3,.7))])
        streams=[]; histories=[]
        for displayed in (False,True):
            directory=tmp_path/f'{local}-{displayed}'
            stream=Stream();sink=CLIProgress(stream,width=100);streams.append(stream)
            campaign=OptimizationCampaign(d,directory,evaluator=Evaluator(clock),clock=clock)
            campaign.run(100,maximum_candidates=3,progress_callback=sink if displayed else None)
            # Reopen the same definition and durable state; scientific resume is separately covered by local-region tests.
            campaign=OptimizationCampaign(d,directory,evaluator=Evaluator(clock),clock=clock)
            campaign.run(100,maximum_candidates=3,progress_callback=sink if displayed else None)
            histories.append(campaign.history.load())
        keys=('candidate_id','normalized','physical','sequence_index','search_origin','cache_hit')
        assert [[{k:r.get(k) for k in keys} for r in h] for h in histories][0]==[{k:r.get(k) for k in keys} for r in histories[1]]
        assert streams[1].getvalue().count('] ')==6


@pytest.mark.parametrize('error,exit_code',[(KeyboardInterrupt,130),(RuntimeError,2),(TimeoutError,2)])
def test_cli_cleans_terminal_before_error(tmp_path,monkeypatch,error,exit_code):
    from dada_solver.research import progress
    from dada_solver.research.cli import main
    definition(tmp_path);stream=Stream();sink=CLIProgress(stream)
    monkeypatch.setattr(progress,'CLIProgress',lambda **kw:sink)
    def run(self,*args,progress_callback,**kwargs):
        progress_callback(event());raise error('fixture')
    monkeypatch.setattr(OptimizationCampaign,'run',run)
    with pytest.raises(SystemExit) as caught:
        main(['run',str(tmp_path/'study.toml'),'--directory',str(tmp_path/'interrupted'),'--budget','1m'])
    assert caught.value.code==exit_code
    assert stream.getvalue().endswith('\r\x1b[K\n')


def test_width_is_refreshed_from_output_descriptor(monkeypatch):
    import os
    from dada_solver.research import progress
    class Terminal(Stream):
        def fileno(self): return 17
    widths=iter([120,40]);calls=[]
    def size(fd): calls.append(fd);return os.terminal_size((next(widths),24))
    monkeypatch.setattr(progress.os,'get_terminal_size',size)
    stream=Terminal();sink=CLIProgress(stream)
    sink(event());sink(event())
    assert calls==[17,17]
    assert cells(stream.getvalue().split('\x1b[K')[-1])<=39
