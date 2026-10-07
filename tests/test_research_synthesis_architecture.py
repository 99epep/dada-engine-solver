"""Family-owned synthesis, periodic targets and production-geometry inspection."""
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.composed_kinematics import ComposedKinematics
from dada_solver.research.artifacts import MechanismArtifact, MechanismLibrary
from dada_solver.research.families import build_side, parameter_specs, PHYSICAL_FAMILIES
from dada_solver.research.motion_target import MotionTarget, MotionEvent, target_from_study, load_motion_target
from dada_solver.research.motion_refit import MotionRefitRequest
from dada_solver.research.synthesis import (SynthesisRequest, SynthesisPlan, synthesis_protocol,
    release_coordinates, assess_mechanism, mechanism_catalogue)
from dada_solver.research.mechanism_view import (MechanismModel, mechanism_state,
    plot_mechanism, animate_mechanism, plot_motion_comparison, plot_library_catalogue)

SEEDS = json.loads((Path(__file__).parents[1]/'src/dada_solver/research/data/kinematics_seeds.json').read_text())['seeds']


def _pyplot():
    matplotlib = pytest.importorskip('matplotlib')
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def artifact(family, side='small'):
    seed = SEEDS[side][family]
    return MechanismArtifact.create(family, seed['parameters'], settings=seed['settings'])


def target_for(family):
    laws = []
    for side in ('small', 'large'):
        seed = SEEDS[side][family]
        laws.append(build_side(seed['settings'], seed['parameters'], side, CylinderVolumeLimits(.1, .3))[0])
    return MotionTarget.from_kinematics(ComposedKinematics(*laws), source={'study_id': 'synthetic-'+family}, samples=73)


@pytest.mark.parametrize('family,count', [('slider_crank',3), ('four_bar',11), ('six_bar',15)])
def test_family_protocols_own_coordinates_and_discretes(family, count):
    protocol = synthesis_protocol(family)
    assert len(protocol.coordinates('full_local_polish')) == count
    settings = SEEDS['small'][family]['settings']
    specs = parameter_specs(settings, 'small')
    assert set(protocol.coordinates('full_local_polish')) == {k for k, s in specs.items() if s.kind == 'continuous'}
    request = SynthesisRequest('target', family, 'design_exploitation', protocol.stages)
    assert request.mechanism_family == family
    assert len(release_coordinates('paired_thermodynamic', ('small','large'), family=family)) == 2*count
    assert not release_coordinates('hardware_retuning', family=family)
    for stage in {'global_discovery','primary_discovery','downstream_fit','mirror_initialization'}-set(protocol.stages):
        with pytest.raises(ValueError):
            release_coordinates(stage, family=family)
        with pytest.raises(ValueError):
            SynthesisRequest('target', family, 'design_exploitation', (stage,))
    with pytest.raises(ValueError):
        SynthesisRequest('target', family, 'design_exploitation', tuple(reversed(protocol.stages)))


def test_default_six_bar_hierarchy_is_preserved():
    assert len(release_coordinates('primary_discovery')) == 6
    assert len(release_coordinates('downstream_fit')) == 9
    with pytest.raises(ValueError, match='both'):
        release_coordinates('paired_thermodynamic')


@pytest.mark.parametrize('family', PHYSICAL_FAMILIES)
def test_target_exact_normalization_periodicity_and_derivatives(family, tmp_path):
    target = target_for(family)
    assert target.scientific['angle_domain'] == 'study_angle_before_operation_transform'
    assert target.scientific['derivative_variable'] == 'study_angle_radians'
    assert target.scientific['angles_rad'] == pytest.approx(np.linspace(0, 2*math.pi, 73))
    for side in ('small','large'):
        raw = target.scientific['sides'][side]
        seed = SEEDS[side][family]
        law, _ = build_side(seed['settings'], seed['parameters'], side, CylinderVolumeLimits(.1,.3))
        np.testing.assert_allclose(raw['position'], [(law.value(t)-.1)/.2 for t in target.scientific['angles_rad']])
        np.testing.assert_allclose(raw['first_derivative'], [law.value(t,1)/.2 for t in target.scientific['angles_rad']])
        assert raw['position'][0] == pytest.approx(raw['position'][-1])
        assert (raw['second_derivative'] is None) == (family != 'slider_crank')
    target.save(tmp_path/'target.json')
    assert MotionTarget.load(tmp_path/'target.json').data == target.data
    assert target_for(family).content_hash == target.content_hash
    data = target.data
    data['provenance'] = {'note':'presentation only'}
    assert MotionTarget.from_data(data).content_hash == target.content_hash
    data['scientific']['sides']['small']['position'][1] += .01
    with pytest.raises(ValueError, match='hash'):
        MotionTarget.from_data(data)


def test_unavailable_derivatives_are_not_fabricated_and_events_are_extensible():
    limits = CylinderVolumeLimits(.1, .3)
    model = SimpleNamespace(small_volume_limits=limits, large_volume_limits=limits,
        small_cylinder_volume=lambda t: .2+.1*math.cos(t), large_cylinder_volume=lambda t:.2+.1*math.cos(t))
    event = MotionEvent('source_feature', .7, metadata_json='{"descriptor":"custom"}')
    target = MotionTarget.from_kinematics(model, source={'identity':'analytic-position-only'}, events={'small':[event]})
    assert target.scientific['sides']['small']['first_derivative'] is None
    assert target.scientific['sides']['large']['second_derivative'] is None
    assert target.scientific['sides']['small']['events'][-1]['kind'] == 'source_feature'
    raw = target.scientific
    raw['sides']['small']['position'][-1] -= .1
    with pytest.raises(ValueError, match='periodic'):
        MotionTarget.create(raw['angles_rad'], raw['sides'], source=raw['source'])
    raw['sides']['small']['position'][-1] += .1
    raw['sides']['small']['position'][1] = 1.1
    with pytest.raises(ValueError, match='normalized'):
        MotionTarget.create(raw['angles_rad'], raw['sides'], source=raw['source'])
    with pytest.raises(ValueError, match='angles'):
        MotionTarget.create([0, 1, 2], raw['sides'], source=raw['source'])


def test_study_extraction_features_and_refit_contract(tmp_path):
    from dada_solver.research.presets import initialize_kinematics
    from dada_solver.research.schema import load_study
    path = initialize_kinematics(tmp_path/'study.toml', 'hybrid_compact', 'hybrid_compact')
    study = load_study(path)
    target = target_from_study(study, samples=91)
    assert load_motion_target(path, samples=91).content_hash == target.content_hash
    for side in ('small','large'):
        raw = target.scientific['sides'][side]
        kinds = {event['kind'] for event in raw['events']}
        assert {'maximum','minimum','rounding','cadence_change','kink','periodic_seam'} <= kinds
        for event in raw['events']:
            if event['kind'] == 'maximum':
                from dada_solver.research.schema import compile_study
                model = compile_study(study).adapter.build(study.fixed_parameters).kinematics
                limits = getattr(model, side+'_volume_limits')
                assert (getattr(model, side+'_cylinder_volume')(event['angle_rad'])-limits.minimum)/limits.swept == pytest.approx(1.)
            if event['kind'] == 'kink':
                assert event['metadata']['neighborhood_policy'] is None
                assert 0 <= event['metadata']['normalized_position'] <= 1
                assert math.isfinite(event['metadata']['first_derivative_per_rad'])
        assert raw['second_derivative'] is None
    request = MotionRefitRequest(target)
    assert request.destination_family == 'structured_c2_15p'
    with pytest.raises(ValueError):
        MotionRefitRequest(target, destination_family='free_spline')


@pytest.mark.parametrize('family', PHYSICAL_FAMILIES)
def test_synthesis_evidence_artifact_and_diversity(family):
    target = target_for(family)
    request = SynthesisRequest(target.content_hash, family, 'design_exploitation', ('full_local_polish',))
    plan = SynthesisPlan(target, request)
    seed = SEEDS['small'][family]
    a = plan.artifact(seed['parameters'], settings=seed['settings'])
    assert MechanismArtifact.from_data(a.data).data == a.data
    evidence = assess_mechanism(a, target, 'small', mechanical_samples=360)
    assert evidence['fit']['position_rms'] < 1e-12
    assert evidence['fit']['velocity_rms_per_rad'] < 1e-12
    assert evidence['thermodynamic'] is None
    assert evidence['mechanical']['metrics']['maximum_absolute_first_derivative'] is None
    scaled = assess_mechanism(a, target, 'small', mechanical_samples=360, volume_limits=CylinderVolumeLimits(.1,.3))
    assert scaled['fit']['position_rms'] < 1e-12
    assert scaled['mechanical']['metrics']['maximum_absolute_first_derivative'] > 0
    member = plan.member('scaled', {'small':a}, volume_limits={'small':CylinderVolumeLimits(.1,.3)})
    assert member['metadata']['evidence']['small']['mechanical']['volume_limits_m3'] == [.1,.3]
    with pytest.raises(ValueError):
        plan.execute()
    with pytest.raises(ValueError, match='identity'):
        SynthesisPlan(target, SynthesisRequest('wrong',family,'design_exploitation',('full_local_polish',)))


def test_library_mixed_families_catalogue_and_storage(tmp_path):
    from dada_solver.research.synthesis import synthesis_member
    target = target_for('slider_crank')
    members = tuple(synthesis_member(family, {'small':artifact(family)}, target) for family in PHYSICAL_FAMILIES)
    library = MechanismLibrary(members)
    library.save(tmp_path/'library.json')
    loaded = MechanismLibrary.load(tmp_path/'library.json')
    assert loaded.to_data() == library.to_data()
    rows = mechanism_catalogue(loaded)
    assert [row['family_id'] for row in rows] == list(PHYSICAL_FAMILIES)
    assert all(row['fit'] and row['mechanical'] and row['thermodynamic'] is None for row in rows)
    assert loaded.member('four_bar')['mechanisms']['small']['content_hash'] == artifact('four_bar').content_hash
    with pytest.raises(ValueError):
        loaded.member('missing')
    plt = _pyplot()
    plt.close(plot_library_catalogue(loaded))


@pytest.mark.parametrize('family,joint_count', [('slider_crank',3), ('four_bar',6), ('six_bar',9)])
def test_geometry_and_views_share_production_closure(family, joint_count):
    a = artifact(family)
    model = MechanismModel(a)
    for theta in np.linspace(0., 2*math.pi, 17):
        state = mechanism_state(model, float(theta))
        joints = state['joints']
        assert len(joints) == joint_count
        assert np.isfinite(np.array(list(joints.values()))).all()
        assert math.dist(joints['A'], joints['B']) == pytest.approx(1.)
        assert state['position'] == pytest.approx(model.law.value(theta)-1.)
        if family == 'slider_crank':
            assert math.dist(joints['B'], joints['P']) == pytest.approx(model.geometry.rod_over_crank)
            q, _ = model.geometry.normalized(theta)
            assert q == pytest.approx(state['position'])
        elif family == 'six_bar':
            g = model.geometry
            assert math.dist(joints['B'], joints['C']) == pytest.approx(g.primary_coupler)
            assert math.dist(joints['C'], joints['D']) == pytest.approx(g.primary_rocker)
            assert math.dist(joints['E'], joints['F']) == pytest.approx(g.link_ef)
            assert math.dist(joints['F'], joints['G']) == pytest.approx(g.link_gf)
            assert math.dist(joints['H'], joints['P']) == pytest.approx(g.piston_rod)
            np.testing.assert_allclose(joints['P'], g.joint_state(theta)['joints']['P'])
        else:
            assembly = model.geometry
            original = assembly.evaluate(*model.law.model._crank(theta))
            np.testing.assert_allclose(joints['B'], original.crank_pin)
            assert math.dist(joints['B'], joints['C']) == pytest.approx(assembly.loop.coupler_length)
            assert math.dist(joints['C'], joints['D']) == pytest.approx(assembly.loop.rocker_length)
            assert math.dist(joints['H'], joints['P']) == pytest.approx(assembly.slider.connecting_rod_length)
            axis = np.array([math.cos(assembly.slider.axis_angle), math.sin(assembly.slider.axis_angle)])
            origin = np.array([assembly.slider.axis_origin_x, assembly.slider.axis_origin_y])
            assert np.dot(np.array(joints['P'])-origin, axis) == pytest.approx(original.coordinate)
    plt = _pyplot()
    target = target_for(family)
    figure = plot_mechanism(a, target=target, samples=25)
    figure.canvas.draw()
    plt.close(figure)
    figure, animation = animate_mechanism(a, target=target, samples=25)
    figure.canvas.draw()
    animation._func(5)
    plt.close(figure)
    figure, axes = plot_motion_comparison(a, target)
    assert all(ax.get_xlim() == (0., 2*math.pi) for ax in axes)
    figure.canvas.draw()
    plt.close(figure)


def test_cli_target_protocol_refit_and_visualization_never_integrate(tmp_path, monkeypatch, capsys):
    from dada_solver.research.cli import main
    from dada_solver.campaign.evaluator import MachineEvaluator
    from dada_solver.research.presets import initialize_kinematics
    monkeypatch.setattr(MachineEvaluator, 'evaluate_with_control', lambda *a, **k: pytest.fail('Unexpected thermodynamic evaluation'))
    path = initialize_kinematics(tmp_path/'study.toml', 'slider_crank', 'harmonic')
    target_path = tmp_path/'target.json'
    assert main(['motion-target', str(path), '--output', str(target_path), '--samples','41']) == 0
    assert MotionTarget.load(target_path).scientific['source']['candidate_id']
    for family in PHYSICAL_FAMILIES:
        stage = synthesis_protocol(family).stages[0]
        assert main(['mechanism','synthesize',str(target_path),'--family',family,'--stage',stage,'--validate-only']) == 0
        with pytest.raises(SystemExit) as error:
            main(['mechanism','synthesize',str(target_path),'--family',family,'--stage',stage])
        assert error.value.code == 2
    with pytest.raises(SystemExit):
        main(['mechanism','synthesize',str(target_path),'--family','slider_crank','--stage','downstream_fit','--validate-only'])
    assert main(['motion-refit',str(target_path),'--validate-only']) == 0
    with pytest.raises(SystemExit):
        main(['motion-refit',str(target_path)])
    with pytest.raises(SystemExit):
        main(['motion-refit',str(target_path),'--count','10','--validate-only'])
    _pyplot()
    a = artifact('four_bar')
    a.save(tmp_path/'artifact.json')
    assert main(['mechanism','visualize',str(tmp_path/'artifact.json'),'--target',str(target_path),'--no-show']) == 0
    library = MechanismLibrary(({'family_id':'retained', 'mechanisms':{'large':a.data}, 'metadata':{}},))
    library.save(tmp_path/'library.json')
    assert main(['mechanism','catalogue',str(tmp_path/'library.json')]) == 0
    assert main(['mechanism','visualize',str(tmp_path/'library.json'),'--family-id','retained','--side','large','--static','--no-show']) == 0
    with pytest.raises(SystemExit):
        main(['mechanism','visualize',str(tmp_path/'library.json'),'--no-show'])
    assert 'no thermodynamic integration started' in capsys.readouterr().out


def test_extract_stored_candidate_preserves_identity_without_replay(tmp_path, monkeypatch):
    from tests.test_research_cli import definition
    from tests.test_campaign import Clock, Evaluator
    from dada_solver.campaign.runner import OptimizationCampaign
    from dada_solver.research.cli import evaluate
    d = definition(tmp_path)
    clock = Clock()
    campaign = OptimizationCampaign(d, tmp_path/'campaign', evaluator=Evaluator(clock), clock=clock)
    campaign.run(100, maximum_candidates=2)
    records = campaign.history.load()
    with pytest.raises(ValueError, match='selector'):
        load_motion_target(tmp_path/'campaign', samples=41)
    target = load_motion_target(tmp_path/'campaign', candidate=records[1]['candidate_id'], samples=41)
    assert target.scientific['source']['candidate_id'] == records[1]['candidate_id']
    assert target.scientific['source']['source_definition_id'] == d.definition_id
    assert target.content_hash == load_motion_target(tmp_path/'campaign', candidate=records[1]['candidate_id'][:12], samples=41).content_hash
    from dada_solver.campaign.evaluator import MachineEvaluator
    monkeypatch.setattr(MachineEvaluator, 'evaluate_with_control', lambda *a, **k: {'status':'synthetic', 'metrics':{}})
    evaluate(tmp_path/'study.toml', tmp_path/'evaluation.json')
    single = load_motion_target(tmp_path/'evaluation.json', samples=41)
    assert single.scientific['source']['source_definition_id'] == d.definition_id


def test_four_bar_joint_state_respects_direction_phase_and_orientation():
    seed = SEEDS['small']['four_bar']
    params = dict(seed['parameters'], phase_rad=.71, crank_direction=-1,
                  volume_increases_with_coordinate=not seed['parameters']['volume_increases_with_coordinate'])
    a = MechanismArtifact.create('four_bar', params, settings=seed['settings'])
    model = MechanismModel(a)
    for angle in (.1, 1.3, 3.7):
        for side in ('small','large'):
            state = model.law.model.joint_state(angle, side)
            pin, _ = model.law.model._crank(angle)
            np.testing.assert_allclose(state['joints']['B'], pin)
            original = getattr(model.law.model, side+'_assembly').evaluate(*model.law.model._crank(angle))
            assert state['position'] == original.coordinate
            assert state['derivative'] == original.coordinate_derivative
    with pytest.raises(ValueError):
        model.law.model.joint_state(.1, 'invalid')
