"""Read-only search windows and navigation/export of existing report curves."""
import copy
import json
from types import SimpleNamespace
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pytest
from dada_solver.research.search_window import search_window, format_search_window
from dada_solver.research.graph_view import browse_graphs
from dada_solver.research.cli import main
from dada_solver.research.presets import initialize_kinematics


def test_local_windows_log_integer_choice_and_regions():
    rows=[dict(name='microtube.heat_in.length_m',kind='continuous',unit='m',initial=10.,lower=1.,upper=100.,transform='log'),
          dict(name='microtube.heat_in.tube_count',kind='integer',unit='1',initial=5,lower=1,upper=11,encoding='nearest_even_v1'),
          dict(name='valve.heat_in.placement',kind='choice',unit='1',initial='upstream',choices=['upstream','downstream'])]
    centers=[{r['name']:r['initial'] for r in rows},{rows[0]['name']:1.,rows[1]['name']:1,rows[2]['name']:'downstream'}]
    search=dict(domain='local_regions_v1',choice_scope='declared_choices',radius_fraction=.1,regions=[dict(id=str(i),center=c) for i,c in enumerate(centers)])
    before=copy.deepcopy((rows,search));window=search_window(rows,search)
    assert (rows,search)==before
    first=window['parameters'][0]
    assert first['windows'][0]['lower']==pytest.approx(10**.8)
    assert first['windows'][0]['upper']==pytest.approx(10**1.2)
    assert first['windows'][1]['clipped']
    assert type(window['parameters'][1]['windows'][0]['lower']) is int
    assert window['parameters'][2]['windows'][1]['choices']==['upstream','downstream']
    assert 'Exchangers' in format_search_window(window)
    assert 'no numeric distance' in format_search_window(window)


def test_global_and_empty_windows():
    row=dict(name='operation.frequency_hz',kind='continuous',unit='Hz',initial=2.,lower=1.,upper=3.,transform='linear')
    window=search_window([row],dict(domain='fixed_global_bounds'))
    assert window['parameters'][0]['windows'][0]['lower']==1
    assert 'No active parameters' in format_search_window(search_window([],dict(domain='fixed_global_bounds')))


def test_validation_search_window_no_integration(tmp_path,monkeypatch,capsys):
    path=initialize_kinematics(tmp_path/'study.toml')
    monkeypatch.setattr('dada_solver.campaign.evaluator.MachineEvaluator.evaluate_with_control',lambda *a,**k:pytest.fail('integration'))
    assert main(['validate',str(path),'--search-window','--json'])==0
    data=json.loads(capsys.readouterr().out)
    assert data['search_window']['integration_started'] is False
    assert main(['validate',str(path),'--search-window'])==0
    assert 'SEARCH WINDOW' in capsys.readouterr().out


def graph_data():
    return dict(name='Comparison',records=[dict(candidate_id=str(i)*64,plots={
        'positions':dict(angle=[0.,180.,360.],unit='fraction',series=[dict(label='small',values=[0.,1.,0.])],method='Production kinematics'),
        'volumes':dict(angle=[0.,360.],volume_unit='m^3',small=[1.,1.],large=[2.,2.]),
        'pressures':dict(unavailable='No saved periodic state')}) for i in (1,2)])


def test_graph_navigation_overlay_clean_export_no_evaluation(tmp_path,monkeypatch):
    monkeypatch.setattr('dada_solver.campaign.evaluator.MachineEvaluator.evaluate_with_control',lambda *a,**k:pytest.fail('integration'))
    data=graph_data();before=copy.deepcopy(data)
    figure,browser=browse_graphs(data,names=('positions','volumes','pressures'),save_directory=tmp_path)
    try:
        assert len(browser.ax.lines)==2
        browser.previous();assert browser.current_graph=='pressures'
        assert 'No saved periodic state' in browser.ax.texts[0].get_text()
        browser.on_key(SimpleNamespace(key='home'));assert browser.index==0
        browser.on_key(SimpleNamespace(key='right'));assert browser.current_graph=='volumes'
        assert len(browser.ax.lines)==4
        browser.save_current(tmp_path/'volumes.png')
        browser.save_current(tmp_path/'volumes-2.png')
        assert (tmp_path/'volumes.png').stat().st_size>1000
        assert (tmp_path/'volumes-2.png').exists()
        with pytest.raises(ValueError,match='already exists'):browser.save_current(tmp_path/'volumes.png')
        browser.save_current(tmp_path/'volumes.svg')
        assert 'Previous' not in (tmp_path/'volumes.svg').read_text()
        assert data==before
    finally: plt.close(figure)


def test_graphs_cli_no_show_and_selectors(monkeypatch,tmp_path):
    seen={}
    def compare(paths,selectors,**kwargs):
        seen.update(paths=paths,selectors=selectors,**kwargs);return graph_data()
    monkeypatch.setattr('dada_solver.research.report.compare',compare)
    assert main(['graphs','campaign','--candidate','best','--candidate','second','--plots','positions,volumes','--no-show'])==0
    assert seen['selectors']==['best','second']
    assert seen['plots']==('positions','volumes')
    with pytest.raises(SystemExit) as error:
        main(['graphs','campaign','--plots','none','--no-show'])
    assert error.value.code==2


def test_legacy_choice_window_matches_scheduler_without_rewriting():
    from dada_solver.campaign.parameters import parameter_from_mapping, ParameterSpace
    from dada_solver.campaign.scheduled_search import ScheduledSobol
    row=dict(name='valve.heat_in.placement',kind='choice',unit='1',initial='upstream',choices=['upstream','downstream'])
    search=dict(type='sobol',domain='local_regions_v1',seed=19,scramble=True,radius_fraction=.1,
                regions=[dict(id='local',center={row['name']:'upstream'})])
    window=search_window([row],search)
    box=window['parameters'][0]['windows'][0]
    assert box['choices']==['upstream'] and box['policy']=='legacy_encoded_window'
    space=ParameterSpace((parameter_from_mapping({k:v for k,v in row.items() if k!='unit'}),))
    scheduler=ScheduledSobol(space,search)
    for _ in range(20):
        scheduler.next_point()
        assert scheduler.last_physical[row['name']] in box['choices']
    assert 'choice_scope' not in search


def test_graphs_from_real_campaign_kinematics_only(tmp_path,monkeypatch):
    from tests.test_research_cli import definition
    from tests.test_campaign import Clock, Evaluator
    from dada_solver.campaign.runner import OptimizationCampaign
    d=definition(tmp_path);clock=Clock();directory=tmp_path/'campaign'
    OptimizationCampaign(d,directory,evaluator=Evaluator(clock),clock=clock).run(100,maximum_candidates=2)
    monkeypatch.setattr('dada_solver.research.plot_data.replay_thermal',lambda *a,**k:pytest.fail('unexpected replay'))
    monkeypatch.setattr('dada_solver.campaign.evaluator.MachineEvaluator.evaluate_with_control',lambda *a,**k:pytest.fail('integration'))
    assert main(['graphs',str(directory),'--candidate','best','--candidate','second','--plots','positions,volumes,velocity','--no-show'])==0


def test_search_window_display_rounding_preserves_data():
    row=dict(name='kinematics.small.rod_over_crank',kind='continuous',unit='1',
             initial=4.123456789123,lower=3.123456789123,upper=5.123456789123,transform='linear')
    window=search_window([row],dict(domain='fixed_global_bounds'))
    before=copy.deepcopy(window)
    text=format_search_window(window,width=80)
    assert '4.1235' in text and '4.123456789123' not in text
    assert 'kinematics.small.' not in text and 'Mechanism / motion SMALL' in text
    assert max(map(len,text.splitlines()))<=80
    assert window==before


def test_native_toolbar_graph_and_zoom_history(tmp_path):
    from matplotlib.backend_bases import NavigationToolbar2
    figure,browser=browse_graphs(graph_data(),names=('positions','volumes','pressures'),save_directory=tmp_path)
    toolbar=NavigationToolbar2(figure.canvas)
    browser.initialize_toolbar(toolbar)
    try:
        assert len(figure.axes)==1
        assert browser.current_graph=='positions'
        toolbar.forward();assert browser.current_graph=='volumes'
        toolbar.forward();assert browser.current_graph=='pressures'
        toolbar.back();assert browser.current_graph=='volumes'
        browser.ax.set_xlim(20,40);toolbar.push_current()
        toolbar.back();assert browser.ax.get_xlim()==pytest.approx((0,360),abs=20)
        toolbar.forward();assert browser.ax.get_xlim()==pytest.approx((20,40))
        toolbar.home();assert browser.current_graph=='positions'
        browser.on_key(SimpleNamespace(key='pagedown'));assert browser.current_graph=='volumes'
        assert str(tmp_path.resolve()) in figure.canvas.get_default_filename()
        assert not hasattr(browser,'buttons')
    finally: plt.close(figure)


def test_search_window_table_columns_are_preserved():
    row=dict(name='kinematics.small.minimum_primary_transmission_sine',kind='continuous',unit='1',
             initial=.3333333333,lower=.3000000001,upper=.9999999999,transform='linear')
    window=search_window([row],dict(domain='fixed_global_bounds'))
    text=format_search_window(window,width=80)
    table=[line for line in text.splitlines() if '|' in line]
    assert len(table)>=2
    assert all(len(line)<=80 for line in table)
    columns=[i for i,c in enumerate(table[0]) if c=='|']
    assert len(columns)==6
    assert all([i for i,c in enumerate(line) if c=='|']==columns for line in table)
    assert all(name in table[0] for name in ('Initial','Declared','Explored','Center','clipped'))
    narrow=format_search_window(window,width=65)
    assert max(map(len,narrow.splitlines()))<=65


def test_live_machine_uses_production_drawing_and_native_animated_svg(tmp_path,monkeypatch):
    import xml.etree.ElementTree as ET
    from dada_solver.research.machine_render import candidate_mechanisms
    from tests.test_research_machine_render import cycle
    from dada_solver.research.schema import load_study
    study=load_study(initialize_kinematics(tmp_path/'study.toml','slider_crank','slider_crank'))
    artifacts,_=candidate_mechanisms(study,{})
    payload=dict(kind='machine_animation',artifacts={s:a.data for s,a in artifacts.items()},
                 cycle=cycle(),layout='UU',motor_operation=False,method='Production drawing')
    data=dict(name='Machine',records=[dict(candidate_id='a'*64,plots={'mechanisms':payload,'positions':graph_data()['records'][0]['plots']['positions']})])
    monkeypatch.setattr('dada_solver.campaign.evaluator.MachineEvaluator.evaluate_with_control',lambda *a,**k:pytest.fail('integration'))
    figure,browser=browse_graphs(data,names=('mechanisms','positions'))
    try:
        assert browser.ax.patches and browser.ax.lines
        assert not any(text.get_visible() for text in figure.texts)
        assert not browser.ax.get_title() and not browser.ax.texts
        x0,x1,y0,y1=browser._machine()['geometry']['limits']
        assert browser.ax.bbox.width/browser.ax.bbox.height==pytest.approx((x1-x0)/(y1-y0))
        browser.update_frame(20);assert browser.frame_index==20
        browser.on_key(SimpleNamespace(key=' '));assert browser.paused
        # Intercepted native Figure.savefig must create an animation, not a PNG-in-SVG snapshot.
        path=tmp_path/'machine.svg';figure.savefig(path)
        root=ET.parse(path).getroot();ns={'s':'http://www.w3.org/2000/svg'}
        animations=root.findall('.//s:animate',ns)
        assert len(animations)==66
        assert all(a.get('dur')=='5.61s' and a.get('repeatCount')=='indefinite' for a in animations)
        assert root.findall('.//s:path',ns) and root.findall('.//s:linearGradient',ns)
        assert not root.findall('.//s:image',ns)
        ids=[node.get('id') for node in root.iter() if node.get('id')]
        assert len(ids)==len(set(ids))
        assert path.stat().st_size<15_000_000
        browser.next()
        assert browser.heading.get_visible() and browser.status.get_visible()
        assert browser.ax.get_title()=='Positions'
        browser.previous()
        assert not browser.heading.get_visible() and not browser.status.get_visible()
        browser.animation._draw_was_started=True
    finally: plt.close(figure)


def test_vector_animation_source_does_not_encode_webp(tmp_path,monkeypatch):
    from tests.test_research_machine_render import cycle
    from dada_solver.research.schema import load_study
    from dada_solver.research.plot_data import candidate_plots
    study=load_study(initialize_kinematics(tmp_path/'study.toml','slider_crank','slider_crank'))
    monkeypatch.setattr('dada_solver.research.plot_data.replay_thermal',lambda *a,**k:dict(machine_cycle=cycle()))
    monkeypatch.setattr('dada_solver.research.machine_render.render_machine_webp',lambda *a,**k:pytest.fail('WebP encoding'))
    result=candidate_plots(study,dict(candidate_id='a'*64,physical={}),('mechanisms',),mechanism_format='vector')
    assert result['mechanisms']['kind']=='machine_animation'
    assert set(result['mechanisms']['artifacts'])=={'small','large'}
    assert 'data_uri' not in result['mechanisms']


def test_search_window_numeric_pairs_stay_on_one_line():
    rows=[dict(name='kinematics.large.'+name,kind='continuous',unit='crank_radius',
               initial=initial,lower=low,upper=high,transform='linear')
          for name,initial,low,high in [('primary_ground',7.5301,7.5198,7.5438),
                                      ('primary_coupler',10.498,10.484,10.517)]]
    window=search_window(rows,dict(domain='fixed_global_bounds'))
    before=copy.deepcopy(window)
    for width in (65,80,160):
        text=format_search_window(window,width=width)
        for line in text.splitlines():
            if '|' not in line or line.startswith('Parameter'):continue
            cells=line.split('|')
            for column in (3,4):
                value=cells[column].strip()
                if value: assert '[' not in value and ']' not in value and ',' not in value and len(value.split())==2
            assert len(line)<=width
        assert window==before
    wide=format_search_window(window,width=160)
    assert '7.5198   7.5438' in wide
    assert next(line for line in wide.splitlines() if line.startswith('Parameter')).find('|')<40


def test_search_window_sign_aware_spacing():
    from dada_solver.research.search_window import _display_pair
    assert _display_pair((-2.,3.),100)=='-2   3'
    assert _display_pair((-2.,-1.),100)=='-2  -1'
    rows=[dict(name='kinematics.small.offset_over_crank',kind='continuous',unit='crank_radius',initial=-1.,lower=-2.,upper=-.5,transform='linear')]
    window=search_window(rows,dict(domain='fixed_global_bounds'))
    before=copy.deepcopy(window)
    text=format_search_window(window,width=100)
    assert '| Initial' in text and '|-1' in text and '|-2  -0.5' in text
    assert '| r' in text and '| no' in text
    assert window==before
