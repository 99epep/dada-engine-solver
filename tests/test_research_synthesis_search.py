"""Direct geometry discovery, diverse retention and independent local polish."""
from dataclasses import replace
import math
import numpy as np
import pytest

from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.composed_kinematics import ComposedKinematics
from dada_solver.research.families import build_side
from dada_solver.research.motion_target import MotionTarget
from dada_solver.research.artifacts import MechanismArtifact, MechanismLibrary
from dada_solver.research.synthesis import SynthesisPlan, SynthesisRequest
from dada_solver.research.synthesis_search import SearchPolicy, GeometrySearch, geometric_distance, BOUNDS


def source(family):
    p=dict(rod_over_crank=4.5,offset_over_crank=.7,phase_rad=.4,volume_increases_with_coordinate=True) if family=='slider_crank' else dict(
        coupler=3.8,rocker=2.7,ground_x=3.2,ground_y=.2,output_along=2.6,output_normal=.4,rod_length=7.5,
        slider_origin_x=.1,slider_origin_y=-.2,axis_angle=.15,phase_rad=.5,loop_branch=1,slider_branch=1,
        crank_direction=1,volume_increases_with_coordinate=True)
    settings=dict(family=family,**({'output':'rocker'} if family=='four_bar' else {}))
    law,_=build_side(settings,p,'large',CylinderVolumeLimits(1.,2.))
    target=MotionTarget.from_kinematics(ComposedKinematics(law,law),source={'identity':'synthetic-production-'+family},samples=361)
    return target,p,settings


def policy(family):
    categories={'volume_increases_with_coordinate':[True]} if family=='slider_crank' else dict(
        output=['rocker'],loop_branch=[1],slider_branch=[1],crank_direction=[1],volume_increases_with_coordinate=[True])
    return SearchPolicy(islands=1,population=6,generations=2,samples=360,mechanical_samples=360,root_samples=360,
        retain_per_side=6,categories=categories,discovery_polish_evaluations=40)


def plan(target,family,stage='global_discovery',ids=(),constraints=()):
    return SynthesisPlan(target,SynthesisRequest(target.content_hash,family,'design_exploitation',(stage,),
        retained_family_ids=ids,mechanical_constraints=constraints))


@pytest.fixture(scope='module',params=['slider_crank','four_bar'])
def discovered(request):
    family=request.param
    target,_,_=source(family)
    p=policy(family)
    return family,target,p,plan(target,family).execute(policy=p)


def test_discovery_fit_topology_diversity_and_round_trip(discovered,tmp_path):
    family,target,p,library=discovered
    errors=[m['metadata']['evidence']['large']['fit']['position_rms'] for m in library.members]
    assert min(errors)<1e-6
    assert len(library.members)>1
    engine=GeometrySearch(plan(target,family),p,('large',),None)
    category=p.choices(family)[0]
    seed=engine.initial(np.random.default_rng(123),category,0,'large')
    baseline=engine.evaluate(seed,'large',category,{},retain=False)
    assert baseline['valid']
    assert min(errors)**2 < baseline['terms']['weighted_position_mse']*.01
    for i,member in enumerate(library.members):
        artifact=MechanismArtifact.from_data(member['mechanisms']['large'])
        evidence=member['metadata']['evidence']['large']
        assert evidence['thermodynamic'] is None
        assert member['metadata']['thermodynamic'] is None
        assert len(evidence['mechanical']['topology']['extrema'])==2
        assert {x['kind'] for x in evidence['mechanical']['topology']['extrema']}=={'minimum','maximum'}
        assert member['metadata']['search']['score']==pytest.approx(sum((member['metadata']['search']['terms']['weighted_position_mse'],member['metadata']['search']['terms']['velocity_term'])))
        for other in library.members[i+1:]:
            assert geometric_distance(family,artifact.scientific['geometry'],
                MechanismArtifact.from_data(other['mechanisms']['large']).scientific['geometry'],p.coordinates(family))>=p.cluster_distance-1e-12
    library.save(tmp_path/'library.json')
    assert MechanismLibrary.load(tmp_path/'library.json').to_data()==library.to_data()


def test_reproducible_search_and_no_thermodynamic_calls(discovered,monkeypatch):
    from dada_solver.campaign.evaluator import MachineEvaluator
    monkeypatch.setattr(MachineEvaluator,'evaluate_with_control',lambda *a,**k:pytest.fail('Thermodynamics is forbidden'))
    family,target,p,library=discovered
    assert plan(target,family).execute(policy=p).to_data()==library.to_data()
    ids=tuple(m['family_id'] for m in library.members[:2])
    polished=plan(target,family,'full_local_polish',ids).execute(policy=replace(p,polish_evaluations=5),library=library)
    assert len(polished.members)==2
    for member in polished.members:
        origin=member['metadata']['provenance']['origin']
        parent=library.member(origin['parent_family_id'])
        assert origin['parent_artifact_hash']==parent['mechanisms']['large']['content_hash']
        assert member['metadata']['provenance']['categories']==parent['metadata']['provenance']['categories']
        assert member['metadata']['search']['score']<=parent['metadata']['search']['score']+1e-12
        assert member['mechanisms']['large']['scientific']['settings']==parent['mechanisms']['large']['scientific']['settings']
    with pytest.raises(ValueError,match='Research source'):
        plan(target,'six_bar','paired_thermodynamic').execute()


def test_declared_constraints_and_budget_are_explicit():
    t,_,_=source('slider_crank');p=replace(policy('slider_crank'),max_evaluations=12,discovery_polish_evaluations=0)
    constraint=dict(metric='minimum_rod_axis_cosine',relation='minimum',limit=.5,unit='1')
    library=plan(t,'slider_crank',constraints=(constraint,)).execute(policy=p)
    for m in library.members:
        record=m['metadata']['evidence']['large']['mechanical']['constraints'][0]
        assert record['satisfied'] and record['value']>=.5
        assert m['metadata']['search']['evaluations']<=12
        assert m['metadata']['search']['stopping_reason']=='evaluation_limit'
    with pytest.raises(ValueError,match='No mechanically admissible'):
        plan(t,'slider_crank',constraints=(dict(constraint,limit=100),)).execute(policy=p)
    dimensional=dict(metric='maximum_absolute_first_derivative',relation='maximum',limit=100.,unit='m^3/rad')
    with pytest.raises(ValueError,match='volume limits'):
        plan(t,'slider_crank',constraints=(dimensional,)).execute(policy=p)
    scaled=plan(t,'slider_crank',constraints=(dimensional,)).execute(policy=p,volume_limits={'large':CylinderVolumeLimits(.1,.3)})
    assert scaled.members[0]['metadata']['evidence']['large']['mechanical']['volume_limits_m3']==[.1,.3]


def test_periodic_distance_removes_only_drawing_gauges():
    _,p,_=source('four_bar');q=dict(p)
    rotation=.8;c,s=math.cos(rotation),math.sin(rotation)
    for x,y in (('ground_x','ground_y'),('slider_origin_x','slider_origin_y')):
        q[x],q[y]=c*p[x]-s*p[y],s*p[x]+c*p[y]
    q['axis_angle']+=rotation;q['phase_rad']+=rotation
    assert geometric_distance('four_bar',p,q,BOUNDS['four_bar'])<1e-15
    q['slider_origin_x']+=2*math.cos(q['axis_angle']);q['slider_origin_y']+=2*math.sin(q['axis_angle'])
    assert geometric_distance('four_bar',p,q,BOUNDS['four_bar'])<1e-15
    q=dict(p,phase_rad=p['phase_rad']+2*math.pi)
    assert geometric_distance('four_bar',p,q,BOUNDS['four_bar'])<1e-15


@pytest.mark.parametrize('orientation',[True,False])
def test_slider_exact_stationary_points_and_batch_parity(orientation):
    from dada_solver.slider_crank import SliderMotion
    motion=SliderMotion(4.5,.7,.4,orientation)
    angles=np.linspace(0,2*math.pi,101)
    q,v=motion.normalized(angles)
    np.testing.assert_allclose(q,[motion.normalized(float(t))[0] for t in angles],atol=1e-15)
    np.testing.assert_allclose(v,[motion.normalized(float(t))[1] for t in angles],atol=1e-15)
    for point in motion.stationary_points():
        value,derivative=motion.normalized(point['angle_rad'])
        assert abs(derivative)<1e-14
        assert value==pytest.approx(1. if point['kind']=='maximum' else 0.)


@pytest.mark.parametrize('output',['rocker','coupler'])
def test_four_bar_scalar_and_batch_geometry_are_identical(output):
    _,p,settings=source('four_bar');settings['output']=output
    law,assembly=build_side(settings,p,'large',CylinderVolumeLimits(1,2))
    angles=np.linspace(0,2*math.pi,47)
    pins,derivatives=law.model._crank_many(angles)
    batch=assembly.evaluate_many(pins,derivatives)
    for i,t in enumerate(angles):
        scalar=assembly.evaluate(*law.model._crank(float(t)))
        assert batch.coordinate[i]==pytest.approx(scalar.coordinate,abs=1e-13)
        assert batch.coordinate_derivative[i]==pytest.approx(scalar.coordinate_derivative,abs=1e-13)
        for name in ('crank_pin','coupler_joint','output_point'):
            np.testing.assert_allclose(np.array(getattr(batch,name))[:,i],getattr(scalar,name),atol=1e-13)
    np.testing.assert_allclose(law.value(angles),[law.value(float(t)) for t in angles],atol=1e-13)


def test_discrete_categories_and_no_acceleration_fit():
    t,_,_=source('slider_crank');p=replace(policy('slider_crank'),generations=0,discovery_polish_evaluations=0,categories={})
    library=plan(t,'slider_crank').execute(policy=p)
    assert {m['metadata']['provenance']['categories']['volume_increases_with_coordinate'] for m in library.members}=={True,False}
    raw=t.scientific
    for side in raw['sides'].values(): side['second_derivative']=[1e5]*len(raw['angles_rad'])
    altered=MotionTarget.create(raw['angles_rad'],raw['sides'],source=raw['source'])
    other=plan(altered,'slider_crank').execute(policy=p)
    assert [m['mechanisms']['large']['scientific'] for m in library.members]==[m['mechanisms']['large']['scientific'] for m in other.members]
    with pytest.raises(ValueError):replace(p,categories={'volume_increases_with_coordinate':[1]}).choices('slider_crank')
    with pytest.raises(ValueError):replace(p,bounds={'unknown':(0,1)}).coordinates('slider_crank')


def test_cli_catalogue_animation_and_validation(discovered,tmp_path,monkeypatch):
    family,target,p,library=discovered
    matplotlib=pytest.importorskip('matplotlib');matplotlib.use('Agg')
    from dada_solver.research.cli import main
    from dada_solver.research.synthesis_catalogue import render_synthesis_catalogue
    from dada_solver.research.mechanism_view import animate_mechanism
    import matplotlib.pyplot as plt
    from dada_solver.campaign.evaluator import MachineEvaluator
    monkeypatch.setattr(MachineEvaluator,'evaluate_with_control',lambda *a,**k:pytest.fail('Thermodynamics is forbidden'))
    html=render_synthesis_catalogue(library,target,tmp_path/'catalogue.html').read_text()
    assert '<svg' in html and 'sortCatalogue' in html and 'No thermodynamic performance' in html
    for member in library.members: assert member['family_id'] in html
    a=MechanismArtifact.from_data(library.members[0]['mechanisms']['large'])
    figure,animation=animate_mechanism(a,target=target,side='large',samples=25)
    animation._func(4);figure.canvas.draw();plt.close(figure)
    target.save(tmp_path/'target.json')
    assert main(['mechanism','synthesize',str(tmp_path/'target.json'),'--family',family,'--stage','global_discovery','--validate-only'])==0
    config=tmp_path/'policy.toml'
    config.write_text('islands=1\npopulation=4\ngenerations=0\ndiscovery_polish_evaluations=0\nretain_per_side=32\n')
    assert main(['mechanism','synthesize',str(tmp_path/'target.json'),'--family',family,'--stage','global_discovery','--config',str(config),'--output',str(tmp_path/'generated.json')])==0
    assert (tmp_path/'generated.html').is_file()
    generated=MechanismLibrary.load(tmp_path/'generated.json')
    assert generated.members
    config.write_text(config.read_text()+'polish_evaluations=5\n')
    selected=generated.members[0]['family_id']
    assert main(['mechanism','synthesize',str(tmp_path/'target.json'),'--family',family,
        '--stage','full_local_polish','--config',str(config),'--library',str(tmp_path/'generated.json'),
        '--family-id',selected,'--output',str(tmp_path/'polished.json')])==0
    polished=MechanismLibrary.load(tmp_path/'polished.json')
    assert polished.members[0]['metadata']['provenance']['origin']['parent_family_id']==selected
    with pytest.raises(SystemExit):
        main(['mechanism','synthesize',str(tmp_path/'target.json'),'--family',family,'--stage','global_discovery',
              '--output',str(tmp_path/'generated.json')])


def test_position_only_target_preserves_unavailable_velocity():
    target,_,_=source('slider_crank')
    raw=target.scientific
    for side in raw['sides'].values():
        side['first_derivative']=None
        side['second_derivative']=None
    target=MotionTarget.create(raw['angles_rad'],raw['sides'],source=raw['source'])
    p=replace(policy('slider_crank'),generations=0,discovery_polish_evaluations=0)
    library=plan(target,'slider_crank').execute(policy=p)
    for member in library.members:
        assert member['metadata']['search']['terms']['normalized_velocity_mse'] is None
        assert member['metadata']['evidence']['large']['fit']['velocity_rms_per_rad'] is None
    assert len(SearchPolicy().choices('four_bar'))==32


def test_polish_rejects_target_mismatch_missing_side_and_dropped_constraints(discovered):
    family,target,p,library=discovered
    ids=(library.members[0]['family_id'],)
    with pytest.raises(ValueError,match='requested piston'):
        plan(target,family,'full_local_polish',ids).execute(policy=p,library=library,sides=('small',))
    raw=target.scientific
    altered=MotionTarget.create(raw['angles_rad'],raw['sides'],source={'identity':'another source'})
    with pytest.raises(ValueError,match='identities disagree'):
        plan(altered,family,'full_local_polish',ids).execute(policy=p,library=library)


def test_polish_preserves_settings_and_rejects_dropped_constraints():
    target,geometry,settings=source('four_bar')
    settings['envelope_frame_angle_rad']=.37
    constraint=dict(metric='minimum_rod_axis_cosine',relation='minimum',limit=.5,unit='1')
    artifact=MechanismArtifact.create('four_bar',geometry,settings=settings,constraints=(constraint,))
    p=replace(policy('four_bar'),polish_evaluations=3)
    member=plan(target,'four_bar',constraints=(constraint,)).member('chosen',{'large':artifact})
    library=MechanismLibrary((member,))
    with pytest.raises(ValueError,match='retain all source'):
        plan(target,'four_bar','full_local_polish',('chosen',)).execute(policy=p,library=library)
    result=plan(target,'four_bar','full_local_polish',('chosen',),(constraint,)).execute(policy=p,library=library)
    assert result.members[0]['mechanisms']['large']['scientific']['settings']==settings
    assert result.members[0]['metadata']['evidence']['large']['mechanical']['constraints'][0]['satisfied']
