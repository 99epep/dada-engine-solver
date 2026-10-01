"""Categorical valve placements share existing campaign and local-search paths."""
import itertools
import tomllib
import pytest
from dada_solver.campaign.parameters import ChoiceParameter, ParameterSpace
from dada_solver.campaign.candidate import Candidate
from dada_solver.campaign.strategy import SobolStrategy
from dada_solver.campaign.scheduled_search import ScheduledSobol
from dada_solver.campaign.runner import OptimizationCampaign
from dada_solver.research.schema import load_study, compile_study
from dada_solver.research.presets import initialize_v2
from dada_solver.research.study_io import dumps
from dada_solver.research.refine import refine
from dada_solver.research.report import inspect, render_html
from tests.test_campaign import Clock, Evaluator


def rows():
    return [dict(name=f'valve.{side}.placement', unit='1', kind='choice',
                 initial='downstream', choices=['downstream','upstream'])
            for side in ('heat_in','heat_out')]


@pytest.fixture
def study_path(tmp_path):
    path=initialize_v2(tmp_path/'study.toml','harmonic','harmonic')
    raw=tomllib.loads(path.read_text()); raw['parameters'].extend(rows())
    path.write_text(dumps(raw))
    return path


def test_choice_encoding():
    p=ChoiceParameter('placement',['downstream','upstream'],'downstream')
    assert [p.decode(u) for u in (0,.25,.49999,.5,.75,1)]==['downstream']*3+['upstream']*3
    assert p.encode('downstream')==.25 and p.encode('upstream')==.75
    for value in (-.1,1.1,float('nan')):
        with pytest.raises(ValueError): p.decode(value)
    with pytest.raises(ValueError): p.encode('unknown')
    for choices,initial in (([], 'a'),(['a','a'],'a'),(['a'],'b')):
        with pytest.raises(ValueError): ChoiceParameter('x',choices,initial)


@pytest.mark.parametrize('pair',list(itertools.product(('downstream','upstream'),repeat=2)))
def test_compilation_routes_placements(study_path,pair):
    d=compile_study(load_study(study_path))
    physical={p.name:v for p,v in zip(d.space.parameters,pair)}
    config=d.adapter.build(dict(d.fixed_parameters,**physical)).configuration
    assert (config.heat_in_valve_placement,config.heat_out_valve_placement)==pair


def test_toml_rejects_unknown_choice(study_path):
    raw=tomllib.loads(study_path.read_text());raw['parameters'][-1]['choices'].append('sideways')
    study_path.write_text(dumps(raw))
    with pytest.raises(ValueError,match='categorical'): load_study(study_path)


def test_sobol_choices_identity_and_canonical_bins(study_path):
    space=load_study(study_path).space; engine=SobolStrategy(2,seed=7,scramble=True)
    candidates=[Candidate.create(space,engine.next_point(),families={},numerical_settings={},definition_id='test') for _ in range(16)]
    assert len({tuple(c.payload['physical'].values()) for c in candidates})==4
    assert len({c.candidate_id for c in candidates})==4


def test_refine_choice_portable_report_and_resume(study_path,tmp_path):
    d=compile_study(load_study(study_path)); clock=Clock(); source=tmp_path/'source'
    OptimizationCampaign(d,source,evaluator=Evaluator(clock),clock=clock).run(100,maximum_candidates=4)
    records=inspect(source)['records']
    path=refine([source],[r['candidate_id'] for r in records],.2,tmp_path/'local.toml')
    local=compile_study(load_study(path)); search=local.study.data['search']
    assert len(search['regions'])==4
    engine=ScheduledSobol(local.space,search)
    sequence=[engine.next_point() for _ in range(20)]
    resumed=ScheduledSobol(local.space,search,index=7)
    assert [resumed.next_point() for _ in range(13)]==sequence[7:]
    for i,point in enumerate(sequence):
        assert local.space.decode(point)==search['regions'][i%4]['center']
    clock=Clock(); target=tmp_path/'local_campaign'
    OptimizationCampaign(local,target,evaluator=Evaluator(clock),clock=clock).run(100,maximum_candidates=4)
    report=inspect(target)
    html=tmp_path/'report.html'
    render_html(report,html)
    assert 'downstream' in html.read_text()
    assert report['records'][0]['physical']==records[0]['physical']


def test_mixed_coordinates_local_centers():
    from dada_solver.campaign.parameters import ContinuousParameter, IntegerParameter
    from dada_solver.research.local_search import validate_search
    from tests.test_research_refine import settings
    space=ParameterSpace((ContinuousParameter('x',0.,1.,.5),
                          IntegerParameter('count',1,10,5),
                          ChoiceParameter('placement',('downstream','upstream'),'downstream')))
    centers=[dict(x=.5,count=5,placement=value) for value in ('downstream','upstream')]
    search=settings(centers);validate_search(search,space)
    engine=ScheduledSobol(space,search)
    for i in range(18):
        point=engine.next_point();physical=space.decode(point)
        assert physical['placement']==centers[i%2]['placement']
        assert type(physical['count']) is int
        assert .3<=physical['x']<=.7
