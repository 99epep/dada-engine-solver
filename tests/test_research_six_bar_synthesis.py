"""Hierarchical synthesis uses production geometry without thermal integration."""
from dataclasses import replace
import math
import time
import numpy as np
import pytest

from dada_solver.six_bar import SixBarCylinderMechanism, SixBarPrimaryMechanism
from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.composed_kinematics import ComposedKinematics
from dada_solver.research.artifacts import MechanismArtifact, MechanismLibrary
from dada_solver.research.families import build_side, PRIMARY_COORDINATES, DOWNSTREAM_COORDINATES, SIXBAR_CONTINUOUS
from dada_solver.research.motion_target import MotionTarget
from dada_solver.research.synthesis import SynthesisPlan, SynthesisRequest, synthesis_member
from dada_solver.research.synthesis_search import SearchPolicy, GeometrySearch, geometric_distance
from dada_solver.research.synthesis_six_bar import SixBarStageSearch, mirror_geometry, primary_evidence, relative_downstream_bounds


def geometry():
    return dict(primary_ground=3.,primary_coupler=3.5,primary_rocker=2.5,
        primary_e_along=2.5,primary_e_normal=.3,primary_phase=.4,
        second_pivot_x=4.,second_pivot_y=-3.,link_ef=4.,link_gf=4.5,
        h_along_over_ef=.6,h_normal_over_ef=.1,piston_rod=8.,
        slider_axis_offset=0.,slider_axis_angle=.1,primary_branch=1,second_branch=1)


def target_for(p=None,small=None):
    large,_=build_side(dict(family='six_bar'),p or geometry(),'large',CylinderVolumeLimits(1.,2.))
    small,_=build_side(dict(family='six_bar'),small or p or geometry(),'small',CylinderVolumeLimits(1.,2.))
    return MotionTarget.from_kinematics(ComposedKinematics(small,large),source={'identity':'synthetic-six-bar'},samples=721)


def plan(target,stage,ids=(),constraints=()):
    return SynthesisPlan(target,SynthesisRequest(target.content_hash,'six_bar','design_exploitation',(stage,),
        retained_family_ids=ids,mechanical_constraints=constraints))


def policy(**kwargs):
    values=dict(seed=42,islands=1,population=4,generations=0,samples=360,
        mechanical_samples=360,root_samples=360,retain_per_side=2,
        categories=dict(primary_branch=[1],second_branch=[1]),discovery_polish_evaluations=45,polish_evaluations=60)
    values.update(kwargs)
    return SearchPolicy(**values)


def primary_library(target,p=None,identifier='primary-L-001',constraints=()):
    p=p or geometry()
    a=MechanismArtifact.create('six_bar',{n:p[n] for n in (*PRIMARY_COORDINATES,'primary_branch')},
        settings={'component':'primary'},constraints=constraints)
    return MechanismLibrary((dict(family_id=identifier,mechanisms={'large':a.data},metadata=dict(
        target_hash=target.content_hash,evidence={'large':primary_evidence(a,target,'large',360)},
        provenance=dict(stage='primary_discovery',primary_family_id=identifier),thermodynamic=None)),))


def complete_library(target,p=None,identifier='full-L-001'):
    a=MechanismArtifact.create('six_bar',p or geometry())
    return MechanismLibrary((synthesis_member(identifier,{'large':a},target,
        provenance={'primary_family_id':'primary-L-001'}),))


@pytest.fixture(scope='module')
def downstream():
    target=target_for(small=dict(geometry(),primary_phase=.8,slider_axis_angle=.3)); parents=primary_library(target)
    started=time.monotonic()
    library=plan(target,'downstream_fit',('primary-L-001',)).execute(policy=policy(),library=parents)
    return target,parents,library,time.monotonic()-started


def test_shared_engine_hierarchy_coordinate_ownership():
    target=target_for(); parents=primary_library(target)
    stages=dict(primary_discovery=PRIMARY_COORDINATES,downstream_fit=DOWNSTREAM_COORDINATES,
                full_local_polish=SIXBAR_CONTINUOUS,opposite_local_adaptation=SIXBAR_CONTINUOUS)
    for stage,names in stages.items():
        parent=None if stage=='primary_discovery' else parents.members[0] if stage=='downstream_fit' else complete_library(target).members[0]
        engine=SixBarStageSearch(plan(target,stage),policy(),('large',),None,parent)
        assert tuple(engine.names)==names
        assert engine.discover.__func__ is GeometrySearch.discover
        assert engine.polish.__func__ is GeometrySearch.polish
        assert not any('branch' in name for name in engine.names)
    with pytest.raises(ValueError,match='Research source'):
        plan(target,'paired_thermodynamic').execute(sides=('small','large'))
    with pytest.raises(ValueError,match='paired Research source'): plan(target,'hardware_retuning').execute()


def test_primary_discovery_is_inspectable_diverse_and_not_piston_fit(tmp_path):
    target=target_for(); p=policy(population=8,generations=1,retain_per_side=8,
        categories=dict(primary_branch=[-1,1],second_branch=[-1,1]),discovery_polish_evaluations=0)
    library=plan(target,'primary_discovery').execute(policy=p)
    assert len(library.members)>1
    assert {m['mechanisms']['large']['scientific']['geometry']['primary_branch'] for m in library.members}=={-1,1}
    for member in library.members:
        a=MechanismArtifact.from_data(member['mechanisms']['large'])
        assert a.scientific['settings']==dict(family='six_bar',component='primary')
        assert set(a.scientific['geometry'])==set((*PRIMARY_COORDINATES,'primary_branch'))
        assert isinstance(a.reconstruct(),SixBarPrimaryMechanism)
        assert member['metadata']['evidence']['large']['fit'] is None
        assert 'E is not P' in member['metadata']['evidence']['large']['primary_projection']['meaning']
        assert member['metadata']['search']['released_coordinates']==list(PRIMARY_COORDINATES)
        with pytest.raises((ValueError,TypeError)):
            build_side(a.scientific['settings'],a.scientific['geometry'],'large',CylinderVolumeLimits(1.,2.))
    library.save(tmp_path/'primary.json')
    assert MechanismLibrary.load(tmp_path/'primary.json').to_data()==library.to_data()
    assert plan(target,'primary_discovery').execute(policy=p).to_data()==library.to_data()


def test_downstream_exact_fixed_primary_fit_and_constraints(downstream):
    target,parents,library,elapsed=downstream
    assert min(m['metadata']['evidence']['large']['fit']['position_rms'] for m in library.members)<1e-3
    for member in library.members:
        a=MechanismArtifact.from_data(member['mechanisms']['large']);p=a.scientific['geometry']
        assert {n:p[n] for n in PRIMARY_COORDINATES}=={n:geometry()[n] for n in PRIMARY_COORDINATES}
        assert member['metadata']['search']['released_coordinates']==list(DOWNSTREAM_COORDINATES)
        assert len(member['metadata']['evidence']['large']['mechanical']['topology']['extrema'])==2
        assert member['metadata']['provenance']['origin']['parent_artifact_hash']==parents.members[0]['mechanisms']['large']['content_hash']
        assert member['metadata']['thermodynamic'] is None
        assert member['metadata']['search']['effective_bounds']['link_ef'][1]>0
    print(f'Downstream synthetic discovery: {elapsed:.3f}s')


def test_polish_precision_and_pair_independence(downstream,monkeypatch):
    from dada_solver.campaign.evaluator import MachineEvaluator
    monkeypatch.setattr(MachineEvaluator,'evaluate_with_control',lambda *a,**k:pytest.fail('Thermodynamics forbidden'))
    target,_,library,_=downstream
    identifier=min(library.members,key=lambda m:m['metadata']['evidence']['large']['fit']['position_rms'])['family_id']
    polished=plan(target,'full_local_polish',(identifier,)).execute(policy=policy(polish_evaluations=90),library=library)
    row=polished.members[0];error=row['metadata']['evidence']['large']['fit']['position_rms']
    assert error<1e-5
    parent=library.member(identifier)['mechanisms']['large']['scientific']['geometry']
    new=row['mechanisms']['large']['scientific']['geometry']
    assert (new['primary_branch'],new['second_branch'])==(parent['primary_branch'],parent['second_branch'])
    for name,(lo,hi) in row['metadata']['search']['effective_bounds'].items():
        delta=(new[name]-parent[name]+math.pi)%(2*math.pi)-math.pi if name in ('primary_phase','slider_axis_angle') else new[name]-parent[name]
        assert abs(delta)<=.1*(hi-lo)+1e-10
    mirror=plan(target,'mirror_initialization',(row['family_id'],)).execute(policy=policy(),sides=('small',),library=polished)
    pair=mirror.members[0]
    assert pair['mechanisms']['large']==row['mechanisms']['large']
    assert pair['mechanisms']['small']['provenance']['mirror_policy']=='initialization_only'
    adapted=plan(target,'opposite_local_adaptation',(pair['family_id'],)).execute(
        policy=policy(categories={}),sides=('small',),library=mirror)
    assert adapted.members[0]['mechanisms']['large']==pair['mechanisms']['large']
    assert adapted.members[0]['mechanisms']['small']['scientific']['geometry']!=pair['mechanisms']['small']['scientific']['geometry']
    assert adapted.members[0]['metadata']['search']['released_coordinates']==list(SIXBAR_CONTINUOUS)
    assert adapted.members[0]['metadata']['evidence']['small']['fit']['position_rms']<pair['metadata']['evidence']['small']['fit']['position_rms']
    print(f'Polished synthetic position RMS: {error:.9g}')


@pytest.mark.parametrize('primary_branch',[-1,1])
@pytest.mark.parametrize('second_branch',[-1,1])
def test_production_vector_joints_and_mirror_phase(primary_branch,second_branch):
    p=geometry();p.update(primary_branch=primary_branch,second_branch=second_branch)
    # Long rod tolerates both dyad placements; no added mechanical floor.
    p['piston_rod']=20.
    original=SixBarCylinderMechanism(**p);reflected=SixBarCylinderMechanism(**mirror_geometry(p))
    angles=np.linspace(-7.,7.,361)
    np.testing.assert_allclose(reflected.normalized_motion(angles)[0],original.normalized_motion(-angles)[0],atol=3e-14)
    batch=original.joint_state(angles)
    for i in (0,90,180,360):
        scalar=original.joint_state(float(angles[i]))
        for name in ('B','C','E','F','H','P'):
            np.testing.assert_allclose([v[i] for v in batch['joints'][name]],scalar['joints'][name],atol=1e-13)
    for left,right,length in (('A','B',1.),('B','C',p['primary_coupler']),('C','D',p['primary_rocker']),
                              ('E','F',p['link_ef']),('F','G',p['link_gf']),('H','P',p['piston_rod'])):
        a=np.broadcast_arrays(*batch['joints'][left]);b=np.broadcast_arrays(*batch['joints'][right])
        np.testing.assert_allclose(np.hypot(a[0]-b[0],a[1]-b[1]),length,atol=1e-13)


def test_primary_constraints_are_applied_and_downstream_constraints_deferred():
    target=target_for();constraints=(dict(metric='minimum_primary_transmission_sine',relation='minimum',limit=.1,unit='1'),
        dict(metric='minimum_secondary_transmission_sine',relation='minimum',limit=.1,unit='1'))
    parents=primary_library(target,constraints=constraints);evidence=parents.members[0]['metadata']['evidence']['large']
    assert evidence['mechanical']['constraints'][0]['satisfied']
    assert evidence['mechanical']['deferred_constraints']==[constraints[1]]
    engine=SixBarStageSearch(plan(target,'full_local_polish',constraints=(dict(metric='minimum_rod_axis_cosine',relation='minimum',limit=1.,unit='1'),)),
        policy(),('large',),None,complete_library(target).members[0])
    row=engine.evaluate(engine.encode(geometry()),'large',engine.choices[0],{})
    assert not row['valid'] and row['reason']=='mechanical_constraints'


def test_relative_bounds_and_periodic_distance_are_reproducible():
    primary=primary_library(target_for()).members[0]['mechanisms']['large']
    a=MechanismArtifact.from_data(primary).reconstruct()
    bounds=relative_downstream_bounds(a,policy())
    assert bounds==relative_downstream_bounds(a,policy())
    p=geometry();q=dict(p,primary_phase=p['primary_phase']+2*math.pi,slider_axis_angle=p['slider_axis_angle']-2*math.pi)
    assert geometric_distance('six_bar',p,q,policy().coordinates('six_bar'))<1e-15
    changed=relative_downstream_bounds(a,policy(pivot_envelope_radius=3.))
    assert changed['second_pivot_x'][1]>bounds['second_pivot_x'][1]


def test_primary_and_complete_catalogues_animations_and_cli(downstream,tmp_path,monkeypatch):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from dada_solver.research.synthesis_catalogue import render_synthesis_catalogue
    from dada_solver.research.mechanism_view import animate_mechanism,MechanismModel
    from dada_solver.research.cli import main
    target,parents,full,_=downstream
    for index,library in enumerate((parents,full)):
        path=tmp_path/f'catalogue-{index}.html';render_synthesis_catalogue(library,target,path)
        assert library.members[0]['family_id'] in path.read_text()
        artifact=MechanismArtifact.from_data(library.members[0]['mechanisms']['large'])
        model=MechanismModel(artifact);state=model.state(.3)
        assert len(state['joints'])==(5 if index==0 else 9)
        figure,animation=animate_mechanism(artifact,target=target,side='large',samples=31)
        animation._draw_was_started=True;plt.close(figure)
    source=tmp_path/'target.json';target.save(source)
    assert main(['mechanism','synthesize',str(source),'--family','six_bar','--stage','primary_discovery','--validate-only'])==0
    parents.save(tmp_path/'parents.json')
    assert main(['mechanism','synthesize',str(source),'--family','six_bar','--stage','primary_discovery',
                 '--output',str(tmp_path/'cli.json'),'--islands','1','--population','4','--generations','0','--max-evaluations','8'])==0
    assert MechanismLibrary.load(tmp_path/'cli.json').members


def test_complete_pipeline_with_narrow_synthetic_primary_bounds(monkeypatch):
    discover=GeometrySearch.discover
    def hierarchical_only(engine):
        assert len(engine.names) in (6,9), 'A six-bar global stage cannot release fifteen coordinates.'
        return discover(engine)
    monkeypatch.setattr(GeometrySearch,'discover',hierarchical_only)
    target=target_for();p=geometry()
    bounds={n:(p[n]-1e-5,p[n]+1e-5) for n in PRIMARY_COORDINATES}
    primary=plan(target,'primary_discovery').execute(policy=policy(bounds=bounds,discovery_polish_evaluations=0))
    first=primary.members[0]
    assert set(first['mechanisms']['large']['scientific']['geometry'])==set((*PRIMARY_COORDINATES,'primary_branch'))
    downstream=plan(target,'downstream_fit',(first['family_id'],)).execute(policy=policy(),library=primary)
    best=min(downstream.members,key=lambda m:m['metadata']['evidence']['large']['fit']['position_rms'])
    polished=plan(target,'full_local_polish',(best['family_id'],)).execute(policy=policy(polish_evaluations=100),library=downstream)
    row=polished.members[0]
    assert row['metadata']['evidence']['large']['fit']['position_rms']<1e-5
    assert row['metadata']['provenance']['primary_family_id']==first['family_id']
    assert row['metadata']['provenance']['origin']['parent_family_id']==best['family_id']


def test_parent_lineages_are_never_clustered_together():
    target=target_for();first=primary_library(target)
    second=primary_library(target,dict(geometry(),primary_e_along=2.9),identifier='primary-L-002')
    parents=MechanismLibrary(first.members+second.members)
    ids=tuple(m['family_id'] for m in parents.members)
    result=plan(target,'downstream_fit',ids).execute(policy=policy(discovery_polish_evaluations=0,max_evaluations=8),library=parents)
    assert {m['metadata']['provenance']['origin']['parent_family_id'] for m in result.members}==set(ids)
    for row in result.members:
        parent=parents.member(row['metadata']['provenance']['origin']['parent_family_id'])
        p=row['mechanisms']['large']['scientific']['geometry']
        assert {n:p[n] for n in PRIMARY_COORDINATES}=={n:parent['mechanisms']['large']['scientific']['geometry'][n] for n in PRIMARY_COORDINATES}


def test_downstream_reproducibility_branches_and_no_thermal_calls(monkeypatch):
    from dada_solver.campaign.evaluator import MachineEvaluator
    monkeypatch.setattr(MachineEvaluator,'evaluate_with_control',lambda *a,**k:pytest.fail('Thermodynamics forbidden'))
    target=target_for();parents=primary_library(target)
    request=plan(target,'downstream_fit',('primary-L-001',))
    p=policy(categories=dict(primary_branch=[1],second_branch=[-1,1]),discovery_polish_evaluations=0)
    a=request.execute(policy=p,library=parents);b=request.execute(policy=p,library=parents)
    assert a.to_data()==b.to_data()
    assert {m['mechanisms']['large']['scientific']['geometry']['second_branch'] for m in a.members}=={-1,1}
    for member in a.members:
        artifact=MechanismArtifact.from_data(member['mechanisms']['large'])
        rebuilt=MechanismArtifact.from_data(artifact.data)
        assert rebuilt.content_hash==artifact.content_hash
        assert rebuilt.reconstruct().normalized_motion(.123)==artifact.reconstruct().normalized_motion(.123)


def test_an_unsuccessful_parent_does_not_destroy_other_descendants():
    target=target_for()
    bad=primary_library(target,dict(geometry(),primary_coupler=3.9,primary_rocker=1.95),identifier='near-toggle')
    good=primary_library(target)
    parents=MechanismLibrary(bad.members+good.members)
    constraint=dict(metric='minimum_primary_transmission_sine',relation='minimum',limit=.3,unit='1')
    result=plan(target,'downstream_fit',('near-toggle','primary-L-001'),constraints=(constraint,)).execute(
        policy=policy(discovery_polish_evaluations=0,max_evaluations=8),library=parents)
    assert all(m['metadata']['provenance']['origin']['parent_family_id']=='primary-L-001' for m in result.members)
    failed=result.members[0]['metadata']['search']['failed_parent_searches']
    assert failed[0]['parent_family_id']=='near-toggle'
    assert failed[0]['evaluations']==4
    assert failed[0]['rejections'].get('mechanical_constraints',0)>0
