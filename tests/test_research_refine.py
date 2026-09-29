"""Local search preserves physical inputs and durable deterministic scheduling."""
import copy
import json
import shutil
import tomllib
import pytest
from dada_solver.campaign.candidate import Candidate
from dada_solver.campaign.history import CampaignHistory
from dada_solver.campaign.parameters import ContinuousParameter, IntegerParameter, ParameterSpace
from dada_solver.campaign.runner import OptimizationCampaign
from dada_solver.campaign.scheduled_search import ScheduledSobol
from dada_solver.campaign.strategy import SobolStrategy
from dada_solver.research.cli import main
from dada_solver.research.presets import initialize_v2
from dada_solver.research.refine import refine
from dada_solver.research.report import inspect, render_html, text_report
from dada_solver.research.schema import load_study, compile_study
from dada_solver.research.study_io import dumps
from tests.test_campaign import Clock, Evaluator
from tests.test_research_v2 import activate


def settings(centers, radius=.2):
    return dict(type='sobol',domain='local_regions_v1',seed=19,scramble=True,
        radius_fraction=radius,allocation='round_robin',evaluate_centers=True,
        regions=[dict(id=f'basin_{i+1}',center=c,source_candidate_id='a'*64,
                      source_study_id='b'*64) for i,c in enumerate(centers)])


@pytest.fixture
def source(tmp_path):
    path=initialize_v2(tmp_path/'source.toml','harmonic','harmonic')
    raw=tomllib.loads(path.read_text())
    activate(raw,'volume.swept_ratio',.2)
    path.write_text(dumps(raw))
    definition=compile_study(load_study(path))
    directory=tmp_path/'source_campaign'
    clock=Clock()
    OptimizationCampaign(definition,directory,evaluator=Evaluator(clock),clock=clock).run(100,maximum_candidates=2)
    return directory


def local(source,tmp_path,selectors=None,radius=.2):
    records=CampaignHistory(source).load()
    path=refine([source],selectors or [r['candidate_id'][:12] for r in records],radius,tmp_path/'local.toml')
    return path,compile_study(load_study(path))


def test_refine_portable_exact_center_bounds_and_provenance(source,tmp_path):
    original=inspect(source)
    path,definition=local(source,tmp_path,[original['records'][0]['candidate_id'][:12]])
    raw=tomllib.loads(path.read_text())
    assert raw['search']['regions'][0]['center']==original['records'][0]['physical']
    assert raw['execution']['default_max_candidates']==512
    for before,after in zip(original['scientific']['parameters'],raw['parameters']):
        assert {k:v for k,v in before.items() if k!='initial'}=={k:v for k,v in after.items() if k!='initial'}
    for key in ('policies','objective','constraints','kinematics','numerical'):
        assert definition.study.scientific[key]==original['scientific'][key]
    with pytest.raises(ValueError,match='already exists'):
        refine([source],['best'],.2,path)
    shutil.rmtree(source)
    assert load_study(path).study_id==definition.study.study_id


def test_radius_clipping_transforms_and_integer_center():
    space=ParameterSpace((ContinuousParameter('x',0,10,1),ContinuousParameter('log',1,100,10,'log'),IntegerParameter('n',1,11,2)))
    center={'x':1.,'log':10.,'n':2}
    strategy=ScheduledSobol(space,settings([center]))
    strategy.next_point()
    assert strategy.last_physical==center
    reference=SobolStrategy(3,seed=19,scramble=True)
    for _ in range(20):
        u=reference.next_point()
        expected=space.decode((u[0]*.3,.3+u[1]*.4,u[2]*.3))
        strategy.next_point()
        assert strategy.last_physical==pytest.approx(expected)
        assert type(strategy.last_physical['n']) is int


def test_round_robin_indexes_resume_and_duplicate_centers():
    space=ParameterSpace((ContinuousParameter('x',0,1,.5),))
    config=settings([{'x':.1},{'x':.9},{'x':.1}])
    whole=ScheduledSobol(space,config)
    expected=[]
    for _ in range(13):
        expected.append((whole.next_point(),copy.deepcopy(whole.last_origin)))
    assert [v[1]['kind'] for v in expected[:3]]==['center','center','sobol']
    assert expected[0][1]['region_ids']==['basin_1','basin_3']
    assert [v[1]['region_id'] for v in expected[2:8]]==['basin_1','basin_2','basin_3']*2
    for split in range(14):
        resumed=ScheduledSobol(space,config,index=split)
        actual=[]
        for _ in range(13-split): actual.append((resumed.next_point(),copy.deepcopy(resumed.last_origin)))
        assert actual==expected[split:]


def test_campaign_centers_cache_resume_snapshots_report(source,tmp_path):
    path,definition=local(source,tmp_path,radius=1.)
    clock=Clock(); evaluator=Evaluator(clock)
    whole=tmp_path/'whole'
    OptimizationCampaign(definition,whole,evaluator=evaluator,clock=clock).run(100,maximum_candidates=10)
    records=CampaignHistory(whole).load()
    assert [r['search_origin']['kind'] for r in records[:3]]==['center','center','sobol']
    assert records[0]['physical']==definition.search_settings['regions'][0]['center']
    assert records[2]['candidate_id']==records[3]['candidate_id']
    assert len(evaluator.calls)==len({r['candidate_id'] for r in records})==4
    split=tmp_path/'split'
    OptimizationCampaign(definition,split,evaluator=Evaluator(clock),clock=clock).run(100,maximum_candidates=3)
    OptimizationCampaign.resume(split,evaluator=Evaluator(clock),clock=clock).run(100,maximum_candidates=7)
    resumed=CampaignHistory(split).load()
    assert [(r['candidate_id'],r['search_origin']) for r in resumed]==[(r['candidate_id'],r['search_origin']) for r in records]
    data=inspect(split)
    assert [r['attempts'] for r in data['local_search']['regions']]==[5,5]
    assert [r['feasible'] for r in data['local_search']['regions']]==[5,5]
    assert 'basin_1' in text_report(data)
    render_html(data,tmp_path/'report.html')
    html=(tmp_path/'report.html').read_text()
    assert 'Local regions' in html and 'Search origin' in html and 'radius_fraction' in html


def test_inflight_center_retried_before_sobol(source,tmp_path):
    _,definition=local(source,tmp_path)
    class Crash:
        def evaluate(self,candidate): raise RuntimeError('interruption')
    directory=tmp_path/'crash';clock=Clock()
    with pytest.raises(RuntimeError,match='interruption'):
        OptimizationCampaign(definition,directory,evaluator=Crash(),clock=clock).run(100,maximum_candidates=1)
    pending=json.loads((directory/'state.json').read_text())['pending']
    OptimizationCampaign.resume(directory,evaluator=Evaluator(clock),clock=clock).run(100,maximum_candidates=3)
    records=CampaignHistory(directory).load()
    assert records[0]['candidate_id']==pending['candidate_id']
    assert [r['search_origin']['kind'] for r in records]==['center','center','sobol']


@pytest.mark.parametrize('enabled',[False,True])
def test_global_initial_is_opt_in_and_not_repeated(source,tmp_path,enabled):
    raw=tomllib.loads((source/'study.toml').read_text())
    if enabled: raw['search']['evaluate_initial']=True
    raw['sources']['machine']['path']=str(source/'basis.json')
    path=tmp_path/'global.toml';path.write_text(dumps(raw))
    definition=compile_study(load_study(path));clock=Clock();directory=tmp_path/'global'
    OptimizationCampaign(definition,directory,evaluator=Evaluator(clock),clock=clock).run(100,maximum_candidates=1)
    OptimizationCampaign.resume(directory,evaluator=Evaluator(clock),clock=clock).run(100,maximum_candidates=2)
    records=CampaignHistory(directory).load()
    if enabled:
        assert records[0]['physical']=={p.name:p.initial for p in definition.space.parameters}
        assert [r['search_origin']['kind'] for r in records]==['initial','sobol','sobol']
    else:
        strategy=SobolStrategy(len(definition.space.parameters),seed=definition.seed,scramble=definition.scramble)
        assert [r['normalized'] for r in records]==[list(strategy.next_point()) for _ in records]
        assert all('search_origin' not in r for r in records)


def test_multiple_standalones_and_incompatibility(source,tmp_path):
    data=inspect(source)
    identity=json.loads((source/'definition.json').read_text())
    artifacts=[]
    for i,record in enumerate(data['records']):
        p=tmp_path/f'evaluation_{i}.json'
        p.write_text(json.dumps(dict(artifact_type='research_evaluation_v1',name='fixture',definition=identity,record=record)))
        artifacts.append(p)
    output=tmp_path/'standalone.toml'
    assert main(['refine',*map(str,artifacts),'--radius','.2','--output',str(output)])==0
    assert len(load_study(output).data['search']['regions'])==2
    # Alter the embedded scientific identity, not physical result values.
    raw=json.loads(artifacts[1].read_text());raw['definition']['study_id']='f'*64
    artifacts[1].write_text(json.dumps(raw))
    with pytest.raises(ValueError): refine(artifacts,[],.2,tmp_path/'incompatible.toml')


@pytest.mark.parametrize('radius',[0,-.1,1.01,float('nan')])
def test_invalid_radius(source,tmp_path,radius):
    with pytest.raises(ValueError,match='radius'):
        refine([source],['best'],radius,tmp_path/'bad.toml')


def test_budget_interrupted_center_stays_pending(source,tmp_path):
    from dada_solver.campaign.evaluator import rejected
    _,definition=local(source,tmp_path)
    class Deadline:
        def evaluate(self,candidate): return rejected('budget_exhausted','fixture deadline')
    directory=tmp_path/'deadline';clock=Clock()
    OptimizationCampaign(definition,directory,evaluator=Deadline(),clock=clock).run(100,maximum_candidates=3)
    assert json.loads((directory/'state.json').read_text())['pending']['search_origin']['kind']=='center'
    OptimizationCampaign.resume(directory,evaluator=Evaluator(clock),clock=clock).run(100,maximum_candidates=3)
    rows=CampaignHistory(directory).load()
    assert rows[0]['candidate_id']==rows[1]['candidate_id']
    assert [r['search_origin']['kind'] for r in rows]==['center','center','center','sobol']
    assert [r['sequence_index'] for r in rows]==[0,0,1,2]


def test_recovery_after_journal_append_before_state_commit(source,tmp_path,monkeypatch):
    _,definition=local(source,tmp_path)
    directory=tmp_path/'recovery';clock=Clock()
    campaign=OptimizationCampaign(definition,directory,evaluator=Evaluator(clock),clock=clock)
    real_save=campaign.history.save
    def crash(record):
        real_save(record)
        raise RuntimeError('crash after append')
    monkeypatch.setattr(campaign.history,'save',crash)
    with pytest.raises(RuntimeError,match='after append'): campaign.run(100,maximum_candidates=1)
    evaluator=Evaluator(clock)
    OptimizationCampaign.resume(directory,evaluator=evaluator,clock=clock).run(100,maximum_candidates=2)
    rows=CampaignHistory(directory).load()
    assert len(rows)==3 and rows[0]['candidate_id'] not in evaluator.calls
    assert [r['search_origin']['kind'] for r in rows]==['center','center','sobol']


def test_integer_duplicates_share_identity_without_origin():
    space=ParameterSpace((IntegerParameter('n',1,3,2),))
    strategy=ScheduledSobol(space,settings([{'n':2},{'n':3}],radius=1.))
    by_value={}
    for _ in range(15):
        strategy.next_point()
        candidate=Candidate.from_physical(space,strategy.last_physical,families={},numerical_settings={},definition_id='fixture')
        value=candidate.payload['physical']['n']
        assert by_value.setdefault(value,candidate.candidate_id)==candidate.candidate_id
        assert 'search_origin' not in candidate.payload
    assert len(by_value)==3


def test_local_discrete_assembly_branches_remain_fixed(tmp_path):
    path=initialize_v2(tmp_path/'mechanism.toml','six_bar','six_bar')
    raw=tomllib.loads(path.read_text());activate(raw,'volume.swept_ratio',.1)
    center={'volume.swept_ratio':next(r['initial'] for r in raw['parameters'] if 'initial' in r)}
    raw['search']=settings([center]);path.write_text(dumps(raw))
    study=load_study(path)
    assert 'kinematics.small.primary_branch' not in center
    assert study.scientific['fixed_parameters']['kinematics.small.primary_branch'] in (-1,1)
    activate(raw,'kinematics.small.primary_branch',.1);path.write_text(dumps(raw))
    with pytest.raises(ValueError,match='must remain fixed'): load_study(path)


def test_scientifically_distinct_valid_centers_rejected(source,tmp_path):
    _,definition=local(source,tmp_path)
    other=tmp_path/'other';clock=Clock()
    OptimizationCampaign(definition,other,evaluator=Evaluator(clock),clock=clock).run(100,maximum_candidates=1)
    with pytest.raises(ValueError,match='incompatible scientific study identities'):
        refine([source,other],['best'],.2,tmp_path/'mixed.toml')


@pytest.mark.parametrize('mutation',[
    lambda s:s.update(evaluate_centers=False),
    lambda s:s.update(allocation='adaptive'),
    lambda s:s['regions'][0]['center'].update(unknown=1),
    lambda s:s['regions'].append(copy.deepcopy(s['regions'][0])),
])
def test_invalid_local_schema(source,tmp_path,mutation):
    path,_=local(source,tmp_path)
    raw=tomllib.loads(path.read_text());mutation(raw['search']);path.write_text(dumps(raw))
    with pytest.raises(ValueError): load_study(path)


def test_rescale_local_region_is_explicitly_unsupported(source,tmp_path):
    from dada_solver.research.rescale import rescale
    _,definition=local(source,tmp_path)
    directory=tmp_path/'local_campaign';clock=Clock()
    OptimizationCampaign(definition,directory,evaluator=Evaluator(clock),clock=clock).run(100,maximum_candidates=1)
    with pytest.raises(ValueError,match='rescale the global source first'):
        rescale(directory,'best',5,tmp_path/'unsafe.toml')
