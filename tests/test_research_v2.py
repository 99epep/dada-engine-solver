"""Declarative ownership, durable campaigns, mixed families and offline evidence."""
import copy
import hashlib
import json
import math
from pathlib import Path
import tomllib
import pytest
from dada_solver.campaign.candidate import Candidate
from dada_solver.campaign.evaluator import MachineEvaluator
from dada_solver.campaign.runner import OptimizationCampaign
from dada_solver.research.cli import main, evaluate
from dada_solver.research.presets import initialize_v2
from dada_solver.research.schema import load_study, compile_study, candidate_for_values
from dada_solver.research.study_io import dumps
from dada_solver.research.synthesis import release_coordinates
from dada_solver.research.report import inspect,compare,render_html
from dada_solver.research.margins import margin_record,enrich_constraints
from dada_solver.research.visualization import sample_motion
from tests.test_campaign import Clock,Evaluator


def rewrite(path,raw):
    path.write_text(dumps(raw));return load_study(path)


def activate(raw,name,span=.01):
    row=next(r for r in raw['parameters'] if r['name']==name); v=row.pop('value')
    row.update(initial=v,lower=v-span,upper=v+span,kind='continuous',transform='linear')
    return row


@pytest.fixture
def path(tmp_path): return initialize_v2(tmp_path/'study.toml')


@pytest.mark.parametrize('family',['harmonic','slider_crank','four_bar','six_bar','free_spline','fourier_c2','structured_c2_15p','ideal_piecewise','four_stage','independent_four_stage'])
def test_all_presets_and_fixed_evaluation_identity(tmp_path,family):
    path=initialize_v2(tmp_path/'study.toml',family,family)
    study=load_study(path); d=compile_study(study)
    assert not study.space.parameters
    c=candidate_for_values(d,{})
    assert c.payload['normalized']==[] and c.payload['physical']=={}
    assert d.families['small']==family
    assert set(d.fixed_parameters)==set(study.scientific['fixed_parameters'])


@pytest.mark.parametrize('small,large',[('four_bar','six_bar'),('structured_c2_15p','free_spline'),('slider_crank','harmonic')])
def test_independent_mixed_sides(tmp_path,small,large):
    path=initialize_v2(tmp_path/'study.toml',small,large); study=load_study(path)
    data=sample_motion(study,samples=19)
    assert len(data['sides']['small']['volume_m3'])==19
    assert data['coupling']=='independent'
    assert data['sides']['large']['physical_stroke_m'] is None
    if small=='four_bar': assert data['sides']['small']['joints']['H']


@pytest.mark.parametrize('mutation,message',[
    (lambda r:r['kinematics']['small'].update(family='unknown'),'Unknown'),
    (lambda r:r['parameters'][0].update(unit='mm'),'requires unit'),
    (lambda r:r['parameters'].append(dict(r['parameters'][0])),'duplicate'),
    (lambda r:r['parameters'][0].update(value=True),'finite'),
    (lambda r:r['parameters'][0].update(name='common.hot_ua'),'wrong-family'),
    (lambda r:r['search'].update(scramble=1),'Sobol'),
    (lambda r:r['numerical'].update(maximum_cycles=True),'positive integer'),
    (lambda r:r['screening'].update(samples=12),'360'),
    (lambda r:r['mechanical_constraints'].append(dict(side='small',metric='minimum_secondary_transmission_sine',relation='minimum',limit=.3,unit='1')),'unavailable'),
    (lambda r:r['warm_start'].update(initial_source='guess'),'uniform'),
])
def test_schema_rejects_ambiguity(path,mutation,message):
    raw=tomllib.loads(path.read_text()); mutation(raw)
    with pytest.raises(ValueError,match=message): rewrite(path,raw)


def test_fixed_and_active_share_one_candidate_vector(path):
    raw=tomllib.loads(path.read_text()); activate(raw,'kinematics.small.phase_rad')
    activate(raw,'operation.frequency_hz',.05)
    row=next(r for r in raw['parameters'] if r['name']=='microtube.heat_in.tube_count'); v=row.pop('value');row.update(kind='integer',initial=v,lower=v-10,upper=v+10,encoding='nearest_even_v1')
    study=rewrite(path,raw);d=compile_study(study)
    candidate=Candidate.create(d.space,[.1,.2,.3],families=d.families,numerical_settings=d.numerical_settings,definition_id=d.definition_id)
    assert len(candidate.payload['physical'])==3
    values=dict(d.fixed_parameters,**candidate.payload['physical'])
    design=d.adapter.build(values)
    assert design.heat_in.bank.tube_count==candidate.payload['physical']['microtube.heat_in.tube_count']
    assert abs(design.configuration.angular_speed)==pytest.approx(2*math.pi*values['operation.frequency_hz'])


@pytest.mark.parametrize('stage,sides,count',[('primary_discovery',('large',),6),('downstream_fit',('large',),9),
    ('full_local_polish',('small',),15),('paired_thermodynamic',('small','large'),30)])
def test_schema_can_release_hierarchical_coordinate_groups(tmp_path,stage,sides,count):
    path=initialize_v2(tmp_path/'study.toml','six_bar','six_bar');raw=tomllib.loads(path.read_text())
    for name in release_coordinates(stage,sides): activate(raw,name,.00001)
    study=rewrite(path,raw); d=compile_study(study)
    assert len(d.space.parameters)==count
    assert all('branch' not in p.name for p in d.space.parameters)
    assert all(f'kinematics.{side}.primary_branch' in d.fixed_parameters for side in ('small','large'))
    assert len(d.space.decode([.5]*count))==count


def test_active_branch_rejected(tmp_path):
    path=initialize_v2(tmp_path/'study.toml','six_bar','six_bar');raw=tomllib.loads(path.read_text());activate(raw,'kinematics.small.primary_branch')
    with pytest.raises(ValueError,match='must remain fixed'): rewrite(path,raw)


def test_fixed_spline_controls_or_nonredundant_active_shape(tmp_path):
    path=initialize_v2(tmp_path/'study.toml','free_spline','harmonic');raw=tomllib.loads(path.read_text())
    invalid=copy.deepcopy(raw);activate(invalid,'kinematics.small.control_0')
    with pytest.raises(ValueError,match='must remain fixed'): rewrite(path,invalid)
    raw['parameters']=[r for r in raw['parameters'] if not r['name'].startswith('kinematics.small.control_')]
    raw['kinematics']['small'].update(count=6,representation='shape_coordinates')
    for i,v in enumerate([.1,.7,-.4,.2]): raw['parameters'].append(dict(name=f'kinematics.small.shape_{i}',value=v,unit='1'))
    # A free exploratory shape need not inherit the historical two-reversal screen.
    raw['mechanical_constraints']=[]
    activate(raw,'kinematics.small.shape_0',.02)
    activate(raw,'kinematics.small.phase_rad',.01)
    study=rewrite(path,raw); assert len(study.space.parameters)==2
    assert sample_motion(study,samples=17)['sides']['small']['normalized_acceleration_per_rad2'] is not None


def test_shared_crank_has_explicit_single_phase_and_direction(tmp_path):
    path=initialize_v2(tmp_path/'study.toml','four_bar','four_bar',coupling='shared_crank');raw=tomllib.loads(path.read_text())
    activate(raw,'kinematics.shared.phase_rad')
    study=rewrite(path,raw);d=compile_study(study);p=dict(d.fixed_parameters,**{r.name:r.initial for r in d.space.parameters});kin=d.adapter.build(p).kinematics
    assert kin.small.model.crank_angle_offset==kin.large.model.crank_angle_offset
    assert kin.small.model.crank_direction==kin.large.model.crank_direction
    raw['parameters'].append(dict(name='kinematics.small.phase_rad',value=1.,unit='rad'))
    with pytest.raises(ValueError,match='wrong-family'): rewrite(path,raw)


def test_geometry_and_family_change_scientific_identity(path,tmp_path):
    a=load_study(path);raw=tomllib.loads(path.read_text());next(r for r in raw['parameters'] if r['name']=='kinematics.small.phase_rad')['value']+=.001
    b=rewrite(path,raw);assert a.study_id!=b.study_id
    c=load_study(initialize_v2(tmp_path/'other.toml','slider_crank','harmonic'));assert c.study_id!=a.study_id
    raw['study']['name']='Only a label';raw['execution']['default_budget']='5m';assert rewrite(path,raw).study_id==b.study_id


def test_mechanical_preflight_records_margin_before_integration(tmp_path,monkeypatch):
    path=initialize_v2(tmp_path/'study.toml','slider_crank','harmonic');raw=tomllib.loads(path.read_text())
    row=activate(raw,'kinematics.small.rod_over_crank',.1);row['lower']=.5
    study=rewrite(path,raw);d=compile_study(study)
    monkeypatch.setattr('dada_solver.campaign.evaluator.solve_periodic_wall_machine',lambda *a,**k:pytest.fail('integration must not start'))
    result=MachineEvaluator(d).evaluate(candidate_for_values(d,{'kinematics.small.rod_over_crank':.5}))
    assert result['status']=='invalid_kinematics' and not result['integrated']
    assert any(c['name']=='small.geometry' and c['state']=='violated' for c in result['constraints'])
    assert any(c['state']=='unavailable' for c in result['constraints'])


def test_constraint_screen_rejects_violation(path,monkeypatch):
    raw=tomllib.loads(path.read_text());raw['mechanical_constraints']=[dict(side='small',metric='maximum_absolute_first_derivative',relation='maximum',limit=.00049,unit='m^3/rad')]
    activate(raw,'volume.total_swept_m3',.0001)
    study=rewrite(path,raw);d=compile_study(study)
    row=d.space.parameters[0]
    # Find a limit that accepts the initial geometry and rejects its larger bound.
    initial=d.adapter.build(dict(d.fixed_parameters,**{row.name:row.initial})).kinematics.small.limits.swept/2
    raw['mechanical_constraints'][0]['limit']=initial*1.01
    d=compile_study(rewrite(path,raw))
    result=MachineEvaluator(d).evaluate(candidate_for_values(d,{row.name:row.upper}))
    assert not result['integrated'] and result['status']=='invalid_kinematics'
    assert any(c['margin'] is not None and c['margin']<0 for c in result['constraints'])


def test_resume_uses_portable_mechanism_snapshots(tmp_path):
    path=initialize_v2(tmp_path/'study.toml','four_bar','six_bar');raw=tomllib.loads(path.read_text());activate(raw,'kinematics.small.phase_rad')
    d=compile_study(rewrite(path,raw));clock=Clock()
    whole=OptimizationCampaign(d,tmp_path/'whole',evaluator=Evaluator(clock),clock=clock);whole.run(100,maximum_candidates=4)
    split=OptimizationCampaign(d,tmp_path/'split',evaluator=Evaluator(clock),clock=clock);result=split.run(100,maximum_candidates=2)
    assert result['suggestions']==[]
    for p in [path,path.with_suffix('.basis.json'),*tmp_path.glob('study.*.mechanism.json')]:p.unlink()
    resumed=OptimizationCampaign.resume(tmp_path/'split',evaluator=Evaluator(clock),clock=clock);resumed.run(100,maximum_candidates=2)
    assert [r['candidate_id'] for r in whole.history.load()]==[r['candidate_id'] for r in resumed.history.load()]
    assert [r['sequence_index'] for r in resumed.history.load()]==[0,1,2,3]
    assert inspect(tmp_path/'split')['runtime_compatible']


def test_cli_all_fixed_is_evaluation_not_search(path,tmp_path,monkeypatch,capsys):
    with pytest.raises(SystemExit) as error:main(['run',str(path),'--directory',str(tmp_path/'run')])
    assert error.value.code==2
    assert 'No active parameters' in capsys.readouterr().err
    assert not (tmp_path/'run').exists()
    with pytest.raises(TypeError):evaluate(path,tmp_path/'bad.json',reference=True)


def test_reports_show_constraint_evidence_without_solver(path,tmp_path,monkeypatch):
    raw=tomllib.loads(path.read_text());activate(raw,'kinematics.small.phase_rad');d=compile_study(rewrite(path,raw));clock=Clock()
    run=OptimizationCampaign(d,tmp_path/'run',evaluator=Evaluator(clock),clock=clock);run.run(100,maximum_candidates=1)
    monkeypatch.setattr(MachineEvaluator,'evaluate',lambda *a:pytest.fail('renderer must not evaluate'))
    data=inspect(tmp_path/'run');r=data['records'][0]
    assert 'kinematics.large.phase_rad' in r['resolved_parameters']
    html=render_html(data,tmp_path/'report.html').read_text()
    assert 'Current value' in html and 'Relative margin' in html and 'near active limit' in html
    assert 'consider a narrower' not in html


def test_margin_contract():
    c=margin_record('pressure',99.,100.,'maximum','Pa')
    assert c['margin']==1 and c['relative_margin']==.01 and c['near_active']
    assert margin_record('floor',2.,3.,'minimum','W')['state']=='violated'
    assert margin_record('zero_limit',0.,0.,'maximum','W')['relative_margin'] is None
    assert margin_record('missing',None,1.,'maximum','W')['state']=='unavailable'
    old=[dict(name='maximum_pressure',margin=1.,available=True,satisfied=True)]
    assert enrich_constraints(old,[dict(type='maximum_pressure',limit=100.,unit='Pa')])[0]['value']==99.


def test_refrigeration_objective_and_human_input_constraint_compile(path):
    raw=tomllib.loads(path.read_text());bpath=path.with_suffix('.basis.json');b=json.loads(bpath.read_text())
    b['configuration']['angular_speed']=abs(b['configuration']['angular_speed']);b['heat_in']=b['heat_out']=b['warm_start']=None
    bpath.write_text(json.dumps(b));raw['sources']['machine']['sha256']=hashlib.sha256(bpath.read_bytes()).hexdigest()
    raw['parameters']=[r for r in raw['parameters'] if not r['name'].startswith(('microtube.','thermal.'))]
    raw['objective']['type']='maximize_cooling_cop';raw['constraints']=[dict(type='valid_thermodynamic_model',unit='1'),dict(type='minimum_cooling_power',required_power=5.,unit='W'),dict(type='maximum_mechanical_input_power',limit=20.,unit='W')]
    d=compile_study(rewrite(path,raw))
    assert d.objective.name=='maximize_cooling_cop' and not d.configuration.motor_operation
    assert d.families['exchanger']=='reservoir'
    assert d.constraints[-1].limit==20.


def test_minimal_artifact_reference_loads_its_complete_settings(tmp_path):
    path=initialize_v2(tmp_path/'study.toml','four_bar','six_bar');raw=tomllib.loads(path.read_text())
    for side in ('small','large'):
        raw['kinematics'][side]={k:v for k,v in raw['kinematics'][side].items() if k in ('family','artifact','sha256')}
    study=rewrite(path,raw)
    assert study.settings['small']['output']=='coupler'
    data=sample_motion(study,samples=9)
    assert len(data['sides']['large']['joints']['H'])==9


def test_slider_frames_reuse_the_motion_geometry(tmp_path):
    path=initialize_v2(tmp_path/'study.toml','slider_crank','harmonic');study=load_study(path)
    data=sample_motion(study,samples=13)['sides']['small']
    rod=study.fixed_parameters['kinematics.small.rod_over_crank']
    for b,p in zip(data['joints']['B'],data['joints']['P']):assert math.dist(b,p)==pytest.approx(rod)
    assert data['physical_stroke_m'] is None
