"""Four-bar design selection is independent of piston fit and parent lineage."""
from dataclasses import replace
import copy
import json
import math

import numpy as np
import pytest

from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.research.artifacts import MechanismArtifact, MechanismLibrary
from dada_solver.research.families import build_side, side_metrics
from dada_solver.research.mechanical_screen import mechanical_screen, inherit_parent_plan, overlay_profile_constraints
from dada_solver.research.synthesis_search import GeometrySearch, CATEGORY_VALUES
from dada_solver.research.synthesis import assess_mechanism
from dada_solver.research.cli import main
from tests.test_research_synthesis_search import source, plan, policy
from tests.test_research_synthesis_architecture import SEEDS


@pytest.fixture(autouse=True)
def no_thermodynamics(monkeypatch):
    from dada_solver.campaign.evaluator import MachineEvaluator
    monkeypatch.setattr(MachineEvaluator,'evaluate_with_control',lambda *a,**k: pytest.fail('Thermodynamic integration'))


def screened(target,stage='global_discovery',ids=(),constraints=(),screen='four_bar_design'):
    p=plan(target,'four_bar',stage,ids,constraints)
    return replace(p,request=replace(p.request,mechanical_screen=screen))


def setup_engine(p=None,**kwargs):
    target,geometry,settings=source('four_bar')
    geometry=dict(geometry,rod_length=9.) if p is None else p
    engine=GeometrySearch(screened(target),replace(policy('four_bar'),**kwargs),('large',),None)
    categories={n:settings['output'] if n=='output' else geometry[n] for n in CATEGORY_VALUES['four_bar']}
    return target,geometry,settings,engine,categories


def constraint(metric,limit):return dict(metric=metric,relation='minimum',limit=limit,unit='1')


def test_packaged_profile_and_unchanged_six_bar():
    four=mechanical_screen('four_bar_design')
    assert (four['name'],four['version'],four['canonical_samples'])==('four_bar_design','four_bar_design_v1',1440)
    assert four['constraints']==[constraint('minimum_primary_transmission_sine',.30),constraint('minimum_rod_axis_cosine',.95),constraint('stroke_over_envelope',.1238)]
    assert len(mechanical_screen('six_bar_design')['constraints'])==10


@pytest.mark.parametrize('family,stage,profile',[
    ('four_bar','global_discovery','four_bar_design'),('four_bar','full_local_polish','four_bar_design'),
    ('six_bar','downstream_fit','six_bar_design'),('slider_crank','global_discovery','none')])
def test_cli_defaults(tmp_path,capsys,family,stage,profile):
    target=source('four_bar')[0];path=tmp_path/'target.json';target.save(path)
    assert main(['mechanism','synthesize',str(path),'--family',family,'--stage',stage,'--validate-only'])==0
    assert json.loads(capsys.readouterr().out)['mechanical_screen']['name']==profile


@pytest.mark.parametrize('limit',[.90,.95,.98])
def test_cli_alignment_override_replaces_default(tmp_path,capsys,limit):
    path=tmp_path/'target.json';source('four_bar')[0].save(path)
    config=tmp_path/'config.toml';config.write_text('[[mechanical_constraints]]\nmetric="minimum_rod_axis_cosine"\nrelation="minimum"\nlimit=0.97\nunit="1"\n')
    args=['mechanism','synthesize',str(path),'--family','four_bar','--stage','global_discovery','--validate-only','--config',str(config)]
    assert main(args)==0
    records=json.loads(capsys.readouterr().out)['mechanical_screen']['effective_constraints']
    assert len(records)==3 and next(r for r in records if r['metric']=='minimum_rod_axis_cosine')['limit']==.97
    assert main([*args,'--minimum-rod-axis-cosine',str(limit)])==0
    records=json.loads(capsys.readouterr().out)['mechanical_screen']['effective_constraints']
    assert len(records)==3 and next(r for r in records if r['metric']=='minimum_rod_axis_cosine')['limit']==limit
    assert main([*args,'--mechanical-screen','none','--minimum-rod-axis-cosine',str(limit)])==0
    records=json.loads(capsys.readouterr().out)['mechanical_screen']['effective_constraints']
    assert records==[constraint('minimum_rod_axis_cosine',limit)]


@pytest.mark.parametrize('value',['-0.1','1.1','nan','inf'])
def test_cli_invalid_alignment(tmp_path,value):
    path=tmp_path/'target.json';source('four_bar')[0].save(path)
    with pytest.raises(SystemExit) as error:
        main(['mechanism','synthesize',str(path),'--family','four_bar','--stage','global_discovery','--minimum-rod-axis-cosine',value,'--validate-only'])
    assert error.value.code==2


@pytest.mark.parametrize('family,screen',[('four_bar','six_bar_design'),('six_bar','four_bar_design'),('slider_crank','four_bar_design')])
def test_incompatible_profiles(tmp_path,family,screen):
    path=tmp_path/'target.json';source('four_bar')[0].save(path)
    with pytest.raises(SystemExit) as error:
        main(['mechanism','synthesize',str(path),'--family',family,'--stage','full_local_polish','--mechanical-screen',screen,'--validate-only'])
    assert error.value.code==2


def test_overlay_intervals_units_and_bool_limits():
    with pytest.raises(ValueError,match='Incompatible units'):
        overlay_profile_constraints([constraint('stroke_over_crank',1.)],[dict(constraint('stroke_over_crank',2.),unit='m')])
    with pytest.raises(ValueError,match='Contradictory'):
        screened(source('four_bar')[0],constraints=(dict(constraint('stroke_over_crank',4.),unit='crank_radius'),dict(constraint('stroke_over_crank',3.),relation='maximum',unit='crank_radius')))
    with pytest.raises(ValueError,match='finite'):
        screened(source('four_bar')[0],constraints=(constraint('minimum_rod_axis_cosine',True),))


@pytest.mark.parametrize('kind,metric',[
    ('transmission','minimum_primary_transmission_sine'),('rod','minimum_rod_axis_cosine'),('envelope','stroke_over_envelope')])
def test_real_geometry_rejection_before_dense_fit(kind,metric):
    _,p,_=source('four_bar');p=dict(p,rod_length=9.)
    if kind=='transmission':p.update(ground_x=5.48,ground_y=0.)
    elif kind=='rod':p['rod_length']=7.5
    else:p['rod_length']=16.
    _,_,_,engine,categories=setup_engine(p)
    # These fixtures genuinely close and have usable stroke; rejection is not
    # a constructor failure masquerading as a design-filter test.
    settings=dict(family='four_bar',output='rocker',envelope_frame='slider_axis')
    law,backend=build_side(settings,p,'large',CylinderVolumeLimits(1.,2.))
    assert len(law.model.stationary_points('large',samples=1440))==2
    result=engine.evaluate(engine.encode(p),'large',categories,{})
    assert not result['valid'] and result['reason']=='mechanical_constraints'
    assert engine.rejections['mechanical_constraints:'+metric]==1
    assert result['mechanical_constraints'] and not engine.archive.rows
    assert 'terms' not in result


def test_alignment_override_changes_real_admission_not_just_metadata():
    target,p,settings=source('four_bar')
    options=policy('four_bar')
    categories={n:settings['output'] if n=='output' else p[n] for n in CATEGORY_VALUES['four_bar']}
    strict=GeometrySearch(screened(target),options,('large',),None)
    assert not strict.evaluate(strict.encode(p),'large',categories,{})['valid']
    relaxed=GeometrySearch(screened(target,constraints=(constraint('minimum_rod_axis_cosine',.90),)),options,('large',),None)
    row=relaxed.evaluate(relaxed.encode(p),'large',categories,{})
    assert row['valid'] and row['terms']['weighted_position_mse']<1e-20
    result=relaxed.final_library().members[0]
    records=result['mechanisms']['large']['scientific']['constraints']
    assert len(records)==3 and next(c for c in records if c['metric']=='minimum_rod_axis_cosine')['limit']==.90


def test_frame_setting_is_explicit_and_not_a_free_coordinate():
    from dada_solver.research.families import validate_settings, parameter_specs
    settings=dict(family='four_bar',output='rocker',envelope_frame='slider_axis')
    validate_settings(settings,'large')
    assert sum(s.kind=='continuous' for s in parameter_specs(settings,'large').values())==11
    for invalid in (dict(settings,envelope_frame='other'),dict(settings,envelope_frame_angle_rad=0.)):
        with pytest.raises(ValueError,match='envelope_frame'):validate_settings(invalid,'large')


def test_reference_frame_rotation_and_axis_changes_are_reproducible():
    for side in ('small','large'):
        seed=SEEDS[side]['four_bar'];settings=seed['settings'];p=seed['parameters']
        old=MechanismArtifact.create('four_bar',p,settings=settings)
        loaded=MechanismArtifact.from_data(old.data)
        assert loaded.data==old.data
        original=assess_mechanism(old,source('four_bar')[0],side)['mechanical']['metrics']
        dynamic=dict(settings,envelope_frame='slider_axis');dynamic.pop('envelope_frame_angle_rad',None)
        law,g=build_side(dynamic,p,side,CylinderVolumeLimits(1.,2.))
        metrics=side_metrics(dynamic,law,g,1440)
        assert metrics['stroke_over_envelope']==pytest.approx(original['stroke_over_envelope'],abs=1e-13)
        rotation=.73;c,s=math.cos(rotation),math.sin(rotation);q=dict(p)
        for x,y in (('ground_x','ground_y'),('slider_origin_x','slider_origin_y')):
            q[x],q[y]=c*p[x]-s*p[y],s*p[x]+c*p[y]
        q['axis_angle']+=rotation;q['phase_rad']+=rotation
        law,g=build_side(dynamic,q,side,CylinderVolumeLimits(1.,2.))
        rotated=side_metrics(dynamic,law,g,1440)
        assert rotated['stroke_over_envelope']==pytest.approx(metrics['stroke_over_envelope'],abs=1e-12)
        assert rotated['envelope_frame_angle_rad']==-q['axis_angle']
        new=MechanismArtifact.create('four_bar',q,settings=dynamic)
        replay=assess_mechanism(MechanismArtifact.from_data(new.data),source('four_bar')[0],side)['mechanical']['metrics']
        assert replay['stroke_over_envelope']==rotated['stroke_over_envelope']


def test_final_filter_1440_and_all_margins(monkeypatch):
    from dada_solver.research import synthesis
    target,p,_,engine,categories=setup_engine()
    assert engine.evaluate(engine.encode(p),'large',categories,{})['valid']
    original=synthesis.assess_mechanism;calls=[]
    def failing(*args,**kwargs):
        calls.append(kwargs['mechanical_samples']);e=original(*args,**kwargs)
        for r in e['mechanical']['constraints']:
            if r['name']=='minimum_rod_axis_cosine':r.update(value=.94,margin=-.01,satisfied=False)
        return e
    monkeypatch.setattr(synthesis,'assess_mechanism',failing)
    with pytest.raises(ValueError,match='No mechanically admissible'):engine.final_library()
    assert calls==[1440]
    assert engine.rejections['final_mechanical_constraints:minimum_rod_axis_cosine']==1
    assert len(engine.rejection_evidence[-1]['constraints'])==3


def test_parent_contexts_strongest_requirements_and_cache_isolation():
    target,p,settings,engine,categories=setup_engine()
    parents=[]
    for name,limit in [('A',.40),('B',.60)]:
        a=MechanismArtifact.create('four_bar',p,settings=settings,constraints=(constraint('minimum_primary_transmission_sine',limit),))
        parent=dict(family_id=name,mechanisms={'large':a.data},metadata={'target_hash':target.content_hash})
        parents.append(parent)
        effective=inherit_parent_plan(engine.plan,parent,'large')
        engine.parent_plans[(name,'large')]=effective
        engine.parents[(name,'large')]=parent
    a=engine.evaluate(engine.encode(p),'large',categories,{'parent_family_id':'A','source_settings':settings})
    b=engine.evaluate(engine.encode(p),'large',categories,{'parent_family_id':'B','source_settings':settings})
    assert a['valid'] and not b['valid'] and engine.cache_hits==0
    assert len(engine.cache)==2
    assert next(c for c in engine.parent_plans[('B','large')].request.mechanical_constraints if c['metric']=='minimum_primary_transmission_sine')['limit']==.6


def test_multi_parent_polish_keeps_constraints_artifacts_and_unique_ids():
    target,p,settings,_,_=setup_engine()
    parents=[]
    for name,limit in [('A',.40),('B',.50)]:
        constraints=(constraint('minimum_primary_transmission_sine',limit),constraint('minimum_rod_axis_cosine',.97))
        # A longer rod satisfies the parent's stricter alignment requirement.
        a=MechanismArtifact.create('four_bar',dict(p,rod_length=12.),settings=settings,constraints=constraints)
        parent=plan(target,'four_bar').member(name,{'large':a})
        parents.append(parent)
    lib=MechanismLibrary(tuple(parents));before=lib.to_data()
    result=screened(target,'full_local_polish',('A','B')).execute(policy=replace(policy('four_bar'),polish_evaluations=2),library=lib)
    assert len({m['family_id'] for m in result.members})==2
    for member in result.members:
        origin=member['metadata']['provenance']['origin'];name=origin['parent_family_id']
        raw=member['mechanisms']['large']['scientific']
        assert next(c for c in raw['constraints'] if c['metric']=='minimum_primary_transmission_sine')['limit']==(.4 if name=='A' else .5)
        assert next(c for c in raw['constraints'] if c['metric']=='minimum_rod_axis_cosine')['limit']==.97
        assert raw['settings']['envelope_frame']=='slider_axis'
        assert member['metadata']['evidence']['large']['mechanical']['samples']==1440
        assert member['metadata']['provenance']['mechanical_screen']['parent_prevents_relaxation']
        assert origin['parent_artifact_hash']==lib.member(name)['mechanisms']['large']['content_hash']
        assert all(c['satisfied'] for c in member['metadata']['evidence']['large']['mechanical']['constraints'])
        MechanismArtifact.from_data(member['mechanisms']['large'])
    assert lib.to_data()==before


def test_category_coverage_has_all_32_and_records_partial_search():
    target,p,settings=source('four_bar')
    options=replace(policy('four_bar'),categories={},retain_per_side=32,discovery_polish_evaluations=0,max_evaluations=3)
    engine=GeometrySearch(screened(target),options,('large',),None)
    engine.discover()
    assert len(engine.choices)==32 and len(engine.category_counts)==32 and len(engine.tasks)==32
    assert sum(c['evaluations'] for c in engine.category_counts.values())==3
    assert engine.stop_reason=='evaluation_limit'
    assert sum(bool(t['population']) for t in engine.tasks)==3
    assert all(t['completed_generations']==0 for t in engine.tasks)
    coverage=engine.category_coverage()
    assert len(coverage)==32
    assert sum(r['islands_launched'] for r in coverage)==3
    assert all(r['stopped_before_completion'] for r in coverage)


def test_polish_screens_export_grid_before_replacing_parent(monkeypatch):
    from dada_solver.research import synthesis_search
    target,p,settings,_,categories=setup_engine()
    engine=GeometrySearch(screened(target,'full_local_polish',('parent',)),policy('four_bar'),('large',),None)
    original=synthesis_search.side_metrics;calls=[]
    def screened_metrics(settings,law,backend,samples):
        calls.append(samples)
        metrics=original(settings,law,backend,samples)
        metrics['minimum_primary_transmission_sine']=.299
        return metrics
    monkeypatch.setattr(synthesis_search,'side_metrics',screened_metrics)
    row=engine.evaluate(engine.encode(p),'large',categories,{'parent_family_id':'parent','source_settings':settings})
    assert calls==[1440]
    assert not row['valid'] and not engine.archive.rows
    assert engine.rejections['mechanical_constraints:minimum_primary_transmission_sine']==1


def test_polish_selects_best_only_after_final_admission(monkeypatch):
    from dada_solver.research import synthesis
    target,p,settings,_,categories=setup_engine()
    engine=GeometrySearch(screened(target,'full_local_polish',('parent',)),policy('four_bar'),('large',),None)
    origin={'parent_family_id':'parent','source_settings':settings}
    first=engine.evaluate(engine.encode(p),'large',categories,origin,retain=False)
    second=engine.evaluate(engine.encode(dict(p,phase_rad=p['phase_rad']+.01)),'large',categories,origin,retain=False)
    assert first['valid'] and second['valid']
    first=dict(first,score=0.);second=dict(second,score=1.)
    engine.archive.rows=[second,first]
    original=synthesis.assess_mechanism;calls=[]
    def assess(*args,**kwargs):
        evidence=original(*args,**kwargs);calls.append(args[0].scientific['geometry']['phase_rad'])
        if len(calls)==1:
            for record in evidence['mechanical']['constraints']:
                if record['name']=='minimum_primary_transmission_sine':record.update(value=.299,margin=-.001,satisfied=False)
        return evidence
    monkeypatch.setattr(synthesis,'assess_mechanism',assess)
    library=engine.final_library()
    assert len(library.members)==1 and len(calls)==2
    assert library.members[0]['mechanisms']['large']['scientific']['geometry']==second['geometry']
    assert engine.rejections['final_mechanical_constraints:minimum_primary_transmission_sine']==1
