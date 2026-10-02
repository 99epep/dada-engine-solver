"""Report-only curves leave campaign records intact and reuse exact replay caches."""
import copy
import json
import shutil
import subprocess
from pathlib import Path
import numpy as np
import pytest
from dada_solver.research import report
from dada_solver.research.plot_data import plot_names, candidate_plots, replay_thermal, kinematic_plots
from dada_solver.research.presets import initialize_v2, initialize_v3
from dada_solver.research.schema import load_study, compile_study, candidate_for_values
from tests.test_research_cockpit import campaign
from tests.test_research_v3 import refrigerator


def test_plot_names():
    assert plot_names('default')==('positions','heat','pressures','temperatures')
    assert plot_names(['volumes,flows','heat'])==('volumes','flows','heat')
    assert plot_names('none')==()
    assert 'mechanisms' in plot_names('all')
    with pytest.raises(ValueError,match='Unknown plot'): plot_names('nonsense')


@pytest.mark.parametrize('family',['harmonic','slider_crank','four_bar','six_bar'])
def test_motion_and_mechanism_samples(tmp_path,family):
    study=load_study(initialize_v2(tmp_path/'study.toml',family,family))
    data=kinematic_plots(study,{'physical':{}},('positions','velocity','acceleration','mechanisms'))
    for series in data['positions']['series']:
        assert min(series['values'])>=-1e-10 and max(series['values'])<=1+1e-10
    assert data['positions']['angle'][-1]==360
    for side in ('small','large'):
        motion=data['mechanisms']['sides'][side]
        assert bool(motion['joints'])==(family!='harmonic')
        if motion['joints']:
            assert all(len(points)==181 for points in motion['joints'].values())
            assert all(a in motion['joints'] and b in motion['joints'] for a,b in motion['links'])


def test_default_targets_best_two_and_explicit_escapes_html_cap(tmp_path,monkeypatch):
    c,_=campaign(tmp_path,4)
    calls=[]
    def generate(study,row,names,**kwargs):
        calls.append(row['candidate_id']);return {'pressures':{'unavailable':'test'}}
    monkeypatch.setattr('dada_solver.research.plot_data.candidate_plots',generate)
    data=report.compare([c.directory],plots='default')
    expected=[r['candidate_id'] for r in data['best'][:2]]
    assert calls==expected and len(calls)==2
    page=report.render_html(data,tmp_path/'report.html').read_text()
    embedded=json.loads(page.split('<script id="data" type="application/json">')[1].split('</script>')[0])
    assert len(embedded['selected'])==2  # ceil(4/10)=1 must not drop the second plot.
    calls.clear();selected=data['records'][-1]['candidate_id']
    data=report.compare([c.directory],[selected],plots='positions')
    assert calls==[selected]


def test_missing_state_is_explained_without_integration(tmp_path,monkeypatch):
    study=load_study(initialize_v3(tmp_path/'study.toml'))
    monkeypatch.setattr('dada_solver.exchangers.wall_cycle.solve_periodic_wall_motor',lambda *a,**k:pytest.fail('integration'))
    row={'candidate_id':'a'*64,'physical':{}}
    result=candidate_plots(study,row,('pressures',))
    assert 'No saved final periodic state' in result['pressures']['unavailable']


def test_real_wall_replay_and_cache(refrigerator,tmp_path,monkeypatch):
    path,definition,result=refrigerator
    study=load_study(path)
    candidate=candidate_for_values(definition,{})
    row=dict(candidate.payload,**result,candidate_id=candidate.candidate_id)
    before=copy.deepcopy(row);messages=[]
    plots=candidate_plots(study,row,('pressures','heat','temperatures','flows'),cache_directory=tmp_path,notify=messages.append)
    assert row==before and not plots['pressures'].get('unavailable')
    assert len(plots['pressures']['series'])==4 and len(plots['heat']['series'])==2
    assert len(plots['temperatures']['series'])==6
    assert all('gas → wall' in series['label'] for series in plots['heat']['series'])
    assert all(v>0 for series in plots['pressures']['series'] for v in series['values'])
    p=plots['heat'];cold=p['series'][0]['values']
    power=np.trapezoid(cold,np.radians(p['angle']))/(2*np.pi)
    assert -power==pytest.approx(row['metrics']['cooling_power_w'],rel=3e-4,abs=1e-5)
    assert any('one cycle' in msg for msg in messages)
    monkeypatch.setattr('dada_solver.research.plot_data.replay_thermal',lambda *a,**k:pytest.fail('cache miss'))
    again=candidate_plots(study,row,('pressures',),cache_directory=tmp_path,notify=messages.append)
    assert again['pressures']==plots['pressures'] and 'matching derived cache' in messages[-1]
    assert len(list(tmp_path.glob('*.json.gz')))==1


def test_cache_identity_includes_state_and_runtime(tmp_path,monkeypatch):
    study=load_study(initialize_v3(tmp_path/'study.toml'))
    calls=[]
    def replay(study,row,**kwargs):
        calls.append(row['candidate_id']);return {'pressures':{'series':[]}}
    monkeypatch.setattr('dada_solver.research.plot_data.replay_thermal',replay)
    monkeypatch.setattr('dada_solver.research.plot_data.runtime_identity',lambda:{'runtime':'one'})
    row={'candidate_id':'a'*64,'physical':{},'final_periodic_state':{'values':[1]}}
    candidate_plots(study,row,('pressures',),cache_directory=tmp_path)
    row['final_periodic_state']['values']=[2]
    candidate_plots(study,row,('pressures',),cache_directory=tmp_path)
    monkeypatch.setattr('dada_solver.research.plot_data.runtime_identity',lambda:{'runtime':'two'})
    candidate_plots(study,row,('pressures',),cache_directory=tmp_path)
    assert len(calls)==3


def test_cli_default_and_none(tmp_path,monkeypatch,capsys):
    from dada_solver.research.cli import main
    c,_=campaign(tmp_path,3);calls=[]
    monkeypatch.setattr('dada_solver.research.plot_data.candidate_plots',lambda s,r,n,**kw:(calls.append(n) or {}))
    assert main(['report',str(c.directory)])==0
    assert calls==[('positions','heat','pressures','temperatures')]*2
    assert (c.directory/'report.html').exists()
    calls.clear();assert main(['report',str(c.directory),'--plots','none'])==0
    assert calls==[]
    assert main(['report',str(c.directory),'--json','--plots','none'])==0
    output=capsys.readouterr().out
    assert 'HTML:' in output


def test_renderer_script_syntax(tmp_path):
    if not shutil.which('node'): pytest.skip('Optional Node syntax check')
    c,_=campaign(tmp_path,1)
    page=report.render_html(report.inspect(c.directory),tmp_path/'report.html').read_text()
    script=page.split('<script>')[1].split('</script>')[0]
    path=tmp_path/'report.js';path.write_text(script)
    subprocess.run(['node','--check',str(path)],check=True,capture_output=True)
    assert 'requestAnimationFrame' in page and 'Mechanism study angle' in page
    assert '<script src=' not in page


def test_second_selector_uses_existing_objective_ranking():
    rows=[dict(candidate_id=str(i)*64,status='feasible',objective={'available':True,'value':value}) for i,value in [(1,-2),(2,-8),(3,-4)]]
    data={'records':rows+[dict(rows[1])]}
    assert report.select_records(data,['best','second'])==[rows[1],rows[2]]
    with pytest.raises(ValueError,match='at least 2 ranked feasible'):
        report.select_records({'records':rows[:1]},['second'])


def test_round_axis_ticks_zero_and_default_plot_order(tmp_path):
    if not shutil.which('node'): pytest.skip('Optional Node runtime')
    script=Path('src/dada_solver/research/report_curves.js').read_text()
    helpers=script.split('function additionalPlots')[0]
    test=helpers+'''
const assert=require('node:assert/strict');
let a=niceAxis(-132.8,478.2,true);
assert.equal(a.step,200);assert.equal(a.low,-200);assert.equal(a.high,600);
assert.deepEqual(a.ticks,[-200,0,200,400,600]);
a=niceAxis(.0062,.0191);assert.equal(tickLabel(a.step,a.step),'0.005');
assert.equal(tickLabel(-1e-16,.1),'0');
'''
    subprocess.run(['node','-e',test],check=True,capture_output=True)
    assert script.index("positions:'")<script.index("heat:'")<script.index("pressures:'")<script.index("temperatures:'")
    assert 'class="zero-line"' in script


def test_comparison_does_not_duplicate_selected_records(tmp_path):
    c,_=campaign(tmp_path,3)
    data=report.compare([c.directory])
    assert len(data['selected'])==3
    data=report.compare([c.directory],[data['records'][0]['candidate_id']])
    assert len(data['selected'])==1
