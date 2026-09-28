"""Dense parity against data frozen before extraction, plus mechanism handoffs."""
from dataclasses import asdict, fields
import json
import math
from pathlib import Path
import numpy as np
import pytest
from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.research.families import build_side, side_metrics, SIXBAR_CONTINUOUS
from dada_solver.research.artifacts import MechanismArtifact, MechanismLibrary
from dada_solver.research.synthesis import release_coordinates, SynthesisRequest, STAGES
from dada_solver.six_bar import SixBarCylinderMechanism

ROOT=Path(__file__).resolve().parents[1]
REFERENCE=json.loads((ROOT/'tests/fixtures/research_v2/reference.json').read_text())['inputs']
SEEDS=json.loads((ROOT/'src/dada_solver/research/data/kinematics_seeds.json').read_text())['seeds']


@pytest.fixture(scope='module')
def trajectories():
    with np.load(ROOT/'tests/fixtures/research_v2/trajectories.npz') as data:
        yield data


@pytest.mark.parametrize('family,key',[('harmonic','harmonic'),('slider_crank','slider_crank'),('four_bar','four_bar'),
    ('six_bar','six_bar_1'),('free_spline','free_spline'),('fourier_c2','fourier_c2'),('structured_c2_15p','structured_c2'),('ideal_piecewise','ideal_piecewise')])
@pytest.mark.parametrize('side,index',[('small',0),('large',1)])
def test_pre_extraction_dense_cycle_parity(family,key,side,index,trajectories):
    seed=SEEDS[side][family]
    law,_=build_side(seed['settings'],seed['parameters'],side,CylinderVolumeLimits(**REFERENCE['limits'][side+'_cylinder']))
    actual=np.array([[law.value(float(t)),law.value(float(t),1)] for t in trajectories['theta']])
    np.testing.assert_allclose(actual,trajectories[key][:,[index,index+2]],rtol=2e-11,atol=1e-13)
    if key in ('structured_c2','free_spline'):
        second=[law.value(float(t),2) for t in trajectories['theta']]
        np.testing.assert_allclose(second,trajectories[key+'_second'][:,index],rtol=2e-11,atol=1e-11)


@pytest.mark.parametrize('rank',[1,4,12,50])
@pytest.mark.parametrize('side,index',[('small',0),('large',1)])
def test_known_sixbar_family_roundtrip_and_margins(rank,side,index,trajectories):
    reference=REFERENCE['sixbar_families'][str(rank)]
    artifact=MechanismArtifact.create('six_bar',reference['geometry'][side],provenance={'family':rank})
    reconstructed=MechanismArtifact.from_data(artifact.data).reconstruct()
    for t in trajectories['theta'][::8]:
        expected=SixBarCylinderMechanism(**reference['geometry'][side]).normalized_motion(float(t))
        np.testing.assert_allclose(reconstructed.normalized_motion(float(t)),expected,rtol=2e-13,atol=2e-13)
        break  # Dense trajectory parity below avoids reconstructing on every sample.
    law,g=build_side(artifact.scientific['settings'],artifact.scientific['geometry'],side,CylinderVolumeLimits(**REFERENCE['limits'][side+'_cylinder']))
    actual=np.array([[law.value(float(t)),law.value(float(t),1)] for t in trajectories['theta']])
    np.testing.assert_allclose(actual,trajectories[f'six_bar_{rank}'][:,[index,index+2]],rtol=2e-11,atol=1e-13)
    metrics=side_metrics(artifact.scientific['settings'],law,g,1440)
    for key in ('stroke_over_crank','minimum_primary_transmission_sine','minimum_secondary_transmission_sine',
                'minimum_rod_axis_cosine','EH_over_crank','H_axis_lateral_rms_over_stroke','H_axis_lateral_span_over_stroke',
                'crank_axis_to_EFH_clearance_over_crank'):
        # Historical diagnostics used sampled stroke; production refines velocity roots.
        assert metrics[key]==pytest.approx(reference['historical_metrics'][side][key],rel=2e-5,abs=2e-5)
    assert metrics['zero_crossing_count']==2
    state=reconstructed.joint_state(.7)
    assert math.dist(state['joints']['H'],state['joints']['P'])==pytest.approx(reconstructed.piston_rod)


@pytest.mark.parametrize('family',['slider_crank','four_bar','six_bar'])
def test_artifact_identity_and_roundtrip(family,tmp_path):
    seed=SEEDS['small'][family]
    artifact=MechanismArtifact.create(family,seed['parameters'],settings=seed['settings'])
    artifact.save(tmp_path/'mechanism.json')
    assert MechanismArtifact.load(tmp_path/'mechanism.json',artifact.content_hash).data==artifact.data
    raw=artifact.data; raw['provenance']['drawing_offset']=[12,34]
    assert MechanismArtifact.from_data(raw).content_hash==artifact.content_hash
    name=next(k for k in seed['parameters'] if 'phase' in k)
    raw['scientific']['geometry'][name]+=.001
    with pytest.raises(ValueError,match='hash'): MechanismArtifact.from_data(raw)


@pytest.mark.parametrize('bad',[0,2,True,1.,'1'])
def test_branches_are_scientific_categories(bad):
    seed=SEEDS['small']['six_bar']; p=dict(seed['parameters'],primary_branch=bad)
    with pytest.raises(ValueError,match='branch'): MechanismArtifact.create('six_bar',p)


def test_synthesis_releases_and_diverse_library():
    assert len(release_coordinates('primary_discovery'))==6
    assert len(release_coordinates('downstream_fit'))==9
    assert len(release_coordinates('full_local_polish'))==15
    assert len(release_coordinates('paired_thermodynamic',('small','large')))==30
    assert not any('branch' in k for k in release_coordinates('paired_thermodynamic',('small','large')))
    request=SynthesisRequest('target-id','six_bar','design_exploitation',STAGES,('compact','long_link'))
    assert request.acceleration_role=='diagnostic_only'
    with pytest.raises(ValueError,match='fresh'): SynthesisRequest('target','six_bar','fresh_island_saturation',('primary_discovery',),('known_seed',))
    with pytest.raises(ValueError,match='both'): release_coordinates('paired_thermodynamic')
    members=[]
    for rank in (1,50):
        members.append(dict(family_id=str(rank),metadata={'role':'diversity'},mechanisms={side:MechanismArtifact.create('six_bar',REFERENCE['sixbar_families'][str(rank)]['geometry'][side]).data for side in ('small','large')}))
    lib=MechanismLibrary(tuple(members));assert MechanismLibrary.from_data(lib.to_data()).to_data()==lib.to_data()


@pytest.mark.parametrize('family',['harmonic','slider_crank','fourier_c2'])
def test_analytic_acceleration(family):
    seed=SEEDS['small'][family]
    law,_=build_side(seed['settings'],seed['parameters'],'small',CylinderVolumeLimits(.1,1.1))
    for t in (.17,1.7,3.7):
        finite_difference=(law.value(t+1e-5,1)-law.value(t-1e-5,1))/(2e-5)
        assert law.value(t,2)==pytest.approx(finite_difference,rel=2e-7,abs=1e-8)


def test_fresh_island_protocol_is_explicit_and_separate():
    from dada_solver.research.synthesis import FreshIslandPolicy
    policy=FreshIslandPolicy({'primary_ground':(1.,3.)},(3,5,7),'fixed_population_and_budget_v1')
    request=SynthesisRequest('target','six_bar','fresh_island_saturation',('primary_discovery',),saturation_policy=policy)
    assert request.saturation_policy.primary_branches==(-1,1)
    with pytest.raises(ValueError,match='Only fresh'):
        SynthesisRequest('target','six_bar','design_exploitation',('primary_discovery',),saturation_policy=policy)
    with pytest.raises(ValueError,match='balance'):
        FreshIslandPolicy({'x':(1.,3.)},(3,5),'fixed',primary_branches=(1,1))
