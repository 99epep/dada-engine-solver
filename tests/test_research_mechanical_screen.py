"""Explicit design profiles, inherited per side, with canonical final screening."""
from dataclasses import replace
from pathlib import Path
import json
import tomllib
import pytest
from dada_solver.research.mechanical_screen import mechanical_screen, inherit_parent_plan
from dada_solver.research.synthesis_six_bar import SixBarStageSearch, mirror_member
from dada_solver.research.artifacts import MechanismArtifact, MechanismLibrary
from dada_solver.research.cli import main
from tests.test_research_six_bar_synthesis import geometry, target_for, plan, policy, complete_library


def screened_plan(target,stage='full_local_polish'):
    p=plan(target,stage)
    return replace(p,request=replace(p.request,mechanical_screen='six_bar_design'))


def engine_for(target=None,**kwargs):
    target=target or target_for()
    return SixBarStageSearch(screened_plan(target),policy(**kwargs),('large',),None,complete_library(target).members[0])


def test_packaged_profile_matches_documentary_reference_exactly():
    reference=tomllib.loads(Path('docs/repro/six_bar_mechanism_families/design_screen.toml').read_text())
    screen=mechanical_screen('six_bar_design')
    assert screen['constraints']==reference['constraints']
    assert screen['canonical_samples']==reference['canonical_samples']==1440
    assert len(screen['constraints'])==10
    assert 'dense_cross_check_samples' not in screen


@pytest.mark.parametrize('stage,option,expected',[
    ('downstream_fit',None,10),('downstream_fit','none',0),('primary_discovery',None,0)])
def test_cli_default_and_opt_out(tmp_path,capsys,stage,option,expected):
    source=tmp_path/'target.json';target_for().save(source)
    args=['mechanism','synthesize',str(source),'--family','six_bar','--stage',stage,'--validate-only']
    if option:args+=['--mechanical-screen',option]
    assert main(args)==0
    result=json.loads(capsys.readouterr().out)
    assert len(result['mechanical_screen']['effective_constraints'])==expected
    assert result['integration_started'] is False


def test_parent_constraints_are_inherited_per_parent_and_none_cannot_remove_them():
    target=target_for();base=complete_library(target).members[0]
    parents=[]
    for limit in (.8,.9):
        a=MechanismArtifact.create('six_bar',geometry(),constraints=[dict(metric='minimum_rod_axis_cosine',relation='minimum',limit=limit,unit='1')])
        parents.append(dict(base,mechanisms={'large':a.data}))
    for parent,limit in zip(parents,(.8,.9)):
        inherited=inherit_parent_plan(plan(target,'full_local_polish'),parent,'large')
        assert [c['limit'] for c in inherited.request.mechanical_constraints]==[limit]
        engine=SixBarStageSearch(plan(target,'full_local_polish'),policy(),('large',),None,parent)
        assert engine.plan.request.mechanical_constraints==inherited.request.mechanical_constraints


@pytest.mark.parametrize('metric,value',[
    ('minimum_secondary_transmission_sine',.1),('minimum_rod_axis_cosine',.8),
    ('crank_axis_to_EFH_clearance_over_crank',.1),('H_axis_lateral_span_over_stroke',.8)])
def test_constraints_filter_before_dense_fit_and_keep_all_margins(monkeypatch,metric,value):
    import dada_solver.research.synthesis_search as search
    original=search.side_metrics
    def failing(*args):return dict(original(*args),**{metric:value})
    monkeypatch.setattr(search,'side_metrics',failing)
    engine=engine_for();row=engine.evaluate(engine.encode(geometry()),'large',engine.choices[0],{})
    assert not row['valid'] and row['reason']=='mechanical_constraints'
    assert engine.rejections['mechanical_constraints']==1
    assert engine.rejections['mechanical_constraints:'+metric]==1
    assert len(row['mechanical_constraints'])==10
    assert any(r['name']==metric and not r['satisfied'] for r in row['mechanical_constraints'])
    assert 'terms' not in row  # Dense target fit was never reached.


def test_final_1440_rejects_a_coarse_pass_and_preserves_multiple_failures(monkeypatch):
    import dada_solver.research.synthesis as six
    engine=engine_for();row=engine.evaluate(engine.encode(geometry()),'large',engine.choices[0],{})
    assert row['valid']
    original=six.assess_mechanism
    def final(*args,**kwargs):
        assert kwargs['mechanical_samples']==1440
        evidence=original(*args,**kwargs)
        for record in evidence['mechanical']['constraints']:
            if record['name'] in ('EH_over_crank','minimum_rod_axis_cosine'):
                record.update(satisfied=False,margin=-.01,state='violated')
        return evidence
    monkeypatch.setattr(six,'assess_mechanism',final)
    with pytest.raises(ValueError,match='No admissible'):engine.final_library()
    assert engine.rejections['final_mechanical_constraints']==1
    assert engine.rejections['final_mechanical_constraints:EH_over_crank']==1
    assert engine.rejections['final_mechanical_constraints:minimum_rod_axis_cosine']==1
    assert len(engine.rejection_evidence[-1]['constraints'])==10


def test_final_profile_hashes_and_optional_fit_filter_do_not_change_score(monkeypatch,tmp_path):
    from dada_solver.campaign.evaluator import MachineEvaluator
    monkeypatch.setattr(MachineEvaluator,'evaluate_with_control',lambda *a,**k:pytest.fail('Integration forbidden'))
    engine=engine_for();row=engine.evaluate(engine.encode(geometry()),'large',engine.choices[0],{})
    library=engine.final_library();member=library.members[0]
    evidence=member['metadata']['evidence']['large']
    assert evidence['mechanical']['samples']==1440
    assert all(c['satisfied'] for c in evidence['mechanical']['constraints'])
    assert len(evidence['mechanical']['topology']['extrema'])==2
    assert evidence['mechanical']['closure']['satisfied']
    artifact=MechanismArtifact.from_data(member['mechanisms']['large'])
    assert artifact.data['provenance']['mechanical_screen']['version']=='six_bar_design_v1'
    library.save(tmp_path/'library.json');assert MechanismLibrary.load(tmp_path/'library.json').to_data()==library.to_data()
    inherited=SixBarStageSearch(plan(engine.plan.target,'full_local_polish'),policy(),('large',),None,member)
    assert len(inherited.plan.request.mechanical_constraints)==10
    same=inherited.evaluate(inherited.encode(geometry()),'large',inherited.choices[0],{})
    assert same['score']==row['score']
    assert inherited.final_library().members[0]['metadata']['search']['mechanical_screen']['activation']=='inherited'
    rejected=engine_for(maximum_position_rms=0.)
    other=dict(geometry(),primary_phase=.41)
    assert rejected.evaluate(rejected.encode(other),'large',rejected.choices[0],{})['valid']
    with pytest.raises(ValueError,match='No admissible'):rejected.final_library()
    assert rejected.rejections['final_maximum_position_rms']==1


def test_mirror_and_opposite_adaptation_inherit_constraints():
    engine=engine_for();engine.evaluate(engine.encode(geometry()),'large',engine.choices[0],{})
    member=engine.final_library().members[0];target=engine.plan.target
    mirrored=mirror_member(plan(target,'mirror_initialization'),member,'small',policy(),None)
    assert len(mirrored['mechanisms']['small']['scientific']['constraints'])==10
    assert mirrored['mechanisms']['large']==member['mechanisms']['large']
    assert all(e['mechanical']['samples']==1440 for e in mirrored['metadata']['evidence'].values())
    adapted=SixBarStageSearch(plan(target,'opposite_local_adaptation'),policy(categories={}),('small',),None,mirrored)
    assert len(adapted.plan.request.mechanical_constraints)==10


def test_catalogue_exposes_design_screen_and_fit_separately(tmp_path):
    from dada_solver.research.synthesis_catalogue import render_synthesis_catalogue
    engine=engine_for();engine.evaluate(engine.encode(geometry()),'large',engine.choices[0],{})
    text=render_synthesis_catalogue(engine.final_library(),engine.plan.target,tmp_path/'catalogue.html').read_text()
    assert 'six_bar_design_v1' in text and 'Mechanical admissibility' in text
    assert '1440' in text and 'Position RMS' in text


@pytest.mark.parametrize('metric,name,value',[
    ('minimum_secondary_transmission_sine','second_pivot_y',2.),
    ('minimum_rod_axis_cosine','piston_rod',2.),
    ('crank_axis_to_EFH_clearance_over_crank','primary_e_along',0.),
    ('H_axis_lateral_span_over_stroke','h_along_over_ef',-.5)])
def test_real_closed_two_turn_mechanisms_fail_individual_design_limits(metric,name,value):
    from dada_solver.six_bar import SixBarCylinderMechanism
    target=target_for();g=dict(geometry(),**{name:value})
    mechanism=SixBarCylinderMechanism(**g)
    assert len(mechanism.stationary_points(samples=360))==2
    constraint=next(c for c in mechanical_screen('six_bar_design')['constraints'] if c['metric']==metric)
    p=plan(target,'full_local_polish',constraints=(constraint,))
    engine=SixBarStageSearch(p,policy(),('large',),None,complete_library(target).members[0])
    row=engine.evaluate(engine.encode(g),'large',engine.choices[0],{})
    assert row['reason']=='mechanical_constraints'
    assert engine.rejections['mechanical_constraints:'+metric]==1
    assert any(r['name']==metric and r['margin']<0 for r in row['mechanical_constraints'])


def test_downstream_keeps_distinct_screened_basins():
    from tests.test_research_six_bar_synthesis import primary_library
    target=target_for();p=screened_plan(target,'downstream_fit')
    parent=primary_library(target).members[0]
    engine=SixBarStageSearch(p,policy(retain_per_side=8),('large',),None,parent)
    for rod in (8.,14.):
        g=dict(geometry(),piston_rod=rod)
        row=engine.evaluate(engine.encode(g),'large',engine.choices[0],{})
        assert row['valid']
    library=engine.final_library()
    assert len(library.members)==2
    assert all(len(m['metadata']['evidence']['large']['mechanical']['constraints'])==10 for m in library.members)


@pytest.mark.parametrize('extra_config',[False,True])
def test_cli_generation_default_screen_and_disabled_parent_inheritance(tmp_path,monkeypatch,capsys,extra_config):
    from tests.test_research_six_bar_synthesis import primary_library
    from dada_solver.campaign.evaluator import MachineEvaluator
    monkeypatch.setattr(MachineEvaluator,'evaluate_with_control',lambda *a,**k:pytest.fail('Integration forbidden'))
    target=target_for();source=tmp_path/'target.json';target.save(source)
    parents=primary_library(target);input_path=tmp_path/'parents.json';parents.save(input_path)
    expected=11 if extra_config else 10
    config_args=[]
    if extra_config:
        config=tmp_path/'search.toml'
        config.write_text('[[mechanical_constraints]]\nmetric = "stroke_over_crank"\nrelation = "maximum"\nlimit = 2.9\nunit = "crank_radius"\n')
        config_args=['--config',str(config)]
    def exact_seed(engine):
        assert len(engine.plan.request.mechanical_constraints)==expected
        categories=next(c for c in engine.choices if c=={'primary_branch':1,'second_branch':1})
        engine.evaluate(engine.encode(geometry()),'large',categories,{})
    monkeypatch.setattr(SixBarStageSearch,'discover',exact_seed)
    output=tmp_path/'output.json'
    assert main(['mechanism','synthesize',str(source),'--family','six_bar','--stage','downstream_fit',
                 '--library',str(input_path),'--family-id','primary-L-001','--output',str(output),*config_args])==0
    library=MechanismLibrary.load(output);member=library.members[0]
    assert len(member['mechanisms']['large']['scientific']['constraints'])==expected
    assert member['metadata']['evidence']['large']['mechanical']['samples']==1440
    capsys.readouterr()
    assert main(['mechanism','synthesize',str(source),'--family','six_bar','--stage','full_local_polish',
                 '--library',str(output),'--family-id',member['family_id'],'--mechanical-screen','none','--validate-only'])==0
    data=json.loads(capsys.readouterr().out)
    screen=data['parent_screens'][0]['mechanical_screen']
    assert screen['name']=='six_bar_design' and screen['activation']=='inherited'
    assert len(screen['effective_constraints'])==expected
