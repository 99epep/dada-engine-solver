"""Scientific refit initialization, analytic topology and exact machine handoff."""
from dataclasses import asdict
import math
import json
from pathlib import Path
from types import SimpleNamespace
import tomllib

import numpy as np
import pytest

from dada_solver.free_kinematics import FreeMotionDefinition, _PeriodicMotion
from dada_solver.hybrid_compact_kinematics import HybridCompactKinematics
from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.research.motion_target import MotionTarget, _hybrid_events
from dada_solver.research.motion_refit import MotionRefitRequest, MotionRefitResult, RefitPolicy, STEP, plot_refit
from dada_solver.research.refit_study import refit_study, shape_bounds
from dada_solver.research.presets import initialize_kinematics
from dada_solver.research.schema import load_study, compile_study
from dada_solver.research.study_io import dumps


@pytest.mark.parametrize('count', [4,8,15,24])
def test_shape_chart_inverse_and_positive_affine_invariance(count):
    rng = np.random.default_rng(251)
    for _ in range(12):
        controls = rng.normal(size=count)
        definition = FreeMotionDefinition(tuple(controls),1.,2.)
        z = definition.to_shape_coordinates()
        rebuilt = FreeMotionDefinition.from_shape_coordinates(z,1.,2.)
        assert len(z) == count-2
        np.testing.assert_allclose(rebuilt.control_values,definition.control_values,atol=2e-14,rtol=2e-14)
        transformed = FreeMotionDefinition(tuple(7.+3.*controls),1.,2.)
        np.testing.assert_allclose(transformed.to_shape_coordinates(),z,atol=3e-14,rtol=3e-14)
    pole = np.r_[np.full(count-1,1/math.sqrt((count-1)*count)), -(count-1)/math.sqrt((count-1)*count)]
    with pytest.raises(ValueError,match='pole'):
        FreeMotionDefinition(tuple(pole),1.,2.).to_shape_coordinates()


def target_of_motion(motion, phases=(.173,.289)):
    limits = CylinderVolumeLimits(1.,2.)
    model = SimpleNamespace(small_volume_limits=limits,large_volume_limits=limits)
    for side,phase in zip(('small','large'),phases):
        for order,suffix in enumerate(('', '_derivative','_second_derivative')):
            setattr(model,side+'_cylinder_volume'+suffix,
                    lambda angle,p=phase,o=order:motion.evaluate(angle-p,o))
    return MotionTarget.from_kinematics(model,source={'test':'shifted production spline'},samples=3001)


@pytest.fixture(scope='module')
def spline_fit():
    angles = np.arange(15)*STEP
    motion = _PeriodicMotion(FreeMotionDefinition(tuple(np.cos(angles)+.08*np.cos(2*angles+.2)),1.,2.))
    target = target_of_motion(motion)
    result = MotionRefitRequest(target,policy=RefitPolicy(phase_samples=16,dense_samples=720,maximum_evaluations=180)).execute()
    return motion,target,result


def test_exact_shifted_spline_refit_and_phase_sign(spline_fit):
    source,target,result = spline_fit
    theta = np.linspace(0,2*math.pi,5001)
    for side,phase in zip(('small','large'),(.173,.289)):
        fitted = result.data['scientific']['sides'][side]
        rebuilt = _PeriodicMotion(FreeMotionDefinition.from_shape_coordinates(fitted['shape_coordinates'],1.,2.))
        assert 0 <= fitted['phase_rad'] < STEP
        # The production convention is q(theta - phase), never theta + phase.
        np.testing.assert_allclose(rebuilt.evaluate(theta-fitted['phase_rad']),source.evaluate(theta-phase),atol=2e-7,rtol=0)
        assert fitted['phase_rad'] == pytest.approx(phase,abs=2e-5)
        assert fitted['diagnostics']['position_rms'] < 4e-8
        assert len(rebuilt.stationary_points()) == fitted['diagnostics']['extrema_count'] == 2
        assert {p['kind'] for p in rebuilt.stationary_points()} == {'maximum','minimum'}
        for order in (0,1,2):
            assert math.isfinite(rebuilt.evaluate(-1e-18,order))
            np.testing.assert_allclose(rebuilt.evaluate(np.array([-1e-18,0.]),order),rebuilt.evaluate(0.,order),atol=1e-14)
    assert MotionRefitResult.from_data(result.data).data == result.data
    bad = result.data
    bad['scientific']['count'] = 10
    with pytest.raises(ValueError,match='hash'):
        MotionRefitResult.from_data(bad)


def hybrid_target(variant=0):
    model = HybridCompactKinematics(CylinderVolumeLimits(1.,2.),CylinderVolumeLimits(1.,2.),
        small_max_deg=140.+5*variant,small_down_duration_deg=175.-5*variant,large_down_duration_deg=190.+5*variant,
        large_down_rounding=.12+.02*variant,small_up_rounding=.14+.02*variant,
        small_down_kink_u=.47,small_down_kink_q=.48,large_up_kink_u=.53,large_up_kink_q=.52)
    if variant == 2:
        seeds = json.loads((Path(__file__).parents[1]/'src/dada_solver/research/data/kinematics_seeds.json').read_text())['seeds']
        parameters = dict(seeds['small']['hybrid_compact']['parameters'],**seeds['large']['hybrid_compact']['parameters'])
        model = HybridCompactKinematics(CylinderVolumeLimits(1.,2.),CylinderVolumeLimits(1.,2.),**parameters)
    events = {side:_hybrid_events(model,side) for side in ('small','large')}
    return model,MotionTarget.from_kinematics(model,source={'test':'hybrid','variant':variant},events=events,samples=1441)


@pytest.fixture(scope='module',params=[0,1,2])
def hybrid_fit(request):
    model,target = hybrid_target(request.param)
    policy = RefitPolicy(phase_samples=16,dense_samples=720,maximum_evaluations=180)
    return model,target,policy,MotionRefitRequest(target,policy=policy).execute()


def test_hybrid_features_topology_quality_and_naive_comparison(hybrid_fit):
    source,target,policy,result = hybrid_fit
    theta = np.linspace(0,2*math.pi,6001)
    for side in ('small','large'):
        raw = result.data['scientific']['sides'][side]
        d,n = raw['diagnostics'],raw['naive']['diagnostics']
        if target.scientific['source']['variant'] == 2:
            assert d['position_rms'] < n['position_rms']*.85
        assert d['terms']['weighted_position_mse'] < n['terms']['weighted_position_mse']
        assert d['position_rms'] < .006
        assert d['maximum_absolute_position_error'] < .02
        assert d['extrema_count'] == 2
        assert max(abs(x) for x in d['extrema_angular_errors_rad'].values()) < STEP/8
        assert d['feature_position_rms'] > 0 and d['velocity_rms'] > 0
        assert d['maximum_absolute_first_derivative'] > 0 and d['maximum_absolute_second_derivative'] > 0
        assert not d['target_acceleration_used']
        assert target.scientific['sides'][side]['second_derivative'] is None
        assert {'kink','rounding','cadence_change','maximum','minimum'} <= {e['kind'] for e in d['events']}
        for event in d['events']:
            if event['kind'] == 'kink':
                source_event = target.scientific['sides'][side]['events'][event['event_index']]
                m = source_event['metadata']
                width = (m['branch_end_rad']-m['branch_start_rad'])%(2*math.pi)
                assert event['radius_rad'] == pytest.approx(min(STEP/2,width/8))
        fitted = _PeriodicMotion(FreeMotionDefinition.from_shape_coordinates(raw['shape_coordinates'],1.,2.))
        exact = np.array([getattr(source,side+'_cylinder_volume')(t) for t in theta])
        error = fitted.evaluate(theta-raw['phase_rad'])-exact
        assert np.sqrt(np.mean(error**2)) == pytest.approx(d['position_rms'],rel=.01)
        print(f"{side}: naive RMS={n['position_rms']:.8g}, feature RMS={d['position_rms']:.8g}, max={d['maximum_absolute_position_error']:.8g}")


def test_refit_reproducible_and_acceleration_not_used(spline_fit):
    _,target,result = spline_fit
    raw = target.scientific
    for side in ('small','large'):
        raw['sides'][side]['second_derivative'] = [1e6]*len(raw['angles_rad'])
    altered = MotionTarget.create(raw['angles_rad'],raw['sides'],source=raw['source'])
    policy = RefitPolicy(phase_samples=16,dense_samples=720,maximum_evaluations=180)
    replay = MotionRefitRequest(target,policy=policy).execute()
    assert replay.payload_json == result.payload_json
    other = MotionRefitRequest(altered,policy=policy).execute()
    assert other.content_hash != result.content_hash
    assert other.data['scientific']['sides'] == result.data['scientific']['sides']


def test_position_only_target_keeps_velocity_unavailable(spline_fit):
    _,target,_ = spline_fit
    raw = target.scientific
    for side in ('small','large'):
        raw['sides'][side]['first_derivative'] = None
        raw['sides'][side]['second_derivative'] = None
    target = MotionTarget.create(raw['angles_rad'],raw['sides'],source=raw['source'])
    result = MotionRefitRequest(target,policy=RefitPolicy(phase_samples=8,dense_samples=360)).execute()
    for side in ('small','large'):
        d = result.data['scientific']['sides'][side]['diagnostics']
        assert d['velocity_rms'] is None
        assert d['terms']['normalized_velocity_mse'] is None
        assert all(e['velocity_error_at_event'] is None for e in d['events'])


def test_chart_search_envelope_is_a_policy_and_excludes_pole(spline_fit):
    controls = spline_fit[2].data['scientific']['sides']['small']['controls']
    low,high,metadata = shape_bounds(controls,.15)
    z = np.array(FreeMotionDefinition(tuple(controls),1.,2.).to_shape_coordinates())
    sphere = np.r_[2*z,z@z-1]/(1+z@z)
    # Check every joint corner, not just independent one-coordinate changes.
    import itertools
    signs = np.array(list(itertools.product((-1.,1.),repeat=13)))
    corners = z+signs*(high-low)/2
    squared = np.sum(corners**2,axis=1)
    points = np.column_stack((2*corners,squared-1))/(1+squared[:,None])
    angles = np.arccos(np.clip(points@sphere,-1.,1.))
    assert np.max(angles) <= .15+1e-12
    np.testing.assert_allclose((low+high)/2,z,atol=1e-15)
    assert metadata['box_inside_cap']
    assert metadata['chart_box_half_width'] > 0
    tighter_low,tighter_high,tighter = shape_bounds(controls)
    assert tighter['canonical_angular_radius_rad'] == .02
    assert np.all(tighter_low > low) and np.all(tighter_high < high)
    assert not metadata['box_is_cap']
    assert metadata['minimum_box_angle_from_excluded_pole_rad'] > 0
    with pytest.raises(ValueError): shape_bounds(controls,0.)
    near_pole = FreeMotionDefinition.from_shape_coordinates([20.]*13,1.,2.)
    with pytest.raises(ValueError,match='pole'): shape_bounds(near_pole.control_values,.15)


def no_thermodynamics(monkeypatch):
    from dada_solver.campaign.evaluator import MachineEvaluator
    monkeypatch.setattr(MachineEvaluator,'evaluate_with_control',lambda *a,**k:pytest.fail('Refit must not integrate thermodynamics'))


def test_selected_candidate_generates_portable_motion_only_study(tmp_path,monkeypatch):
    from tests.test_research_rescale import snapshot
    no_thermodynamics(monkeypatch)
    source = initialize_kinematics(tmp_path/'source.toml','hybrid_compact','hybrid_compact')
    raw = tomllib.loads(source.read_text())
    frequency = next(r for r in raw['parameters'] if r['name'] == 'operation.frequency_hz')
    value = frequency.pop('value')
    frequency.update(initial=value,lower=value/2,upper=value*2,kind='continuous',transform='linear')
    source.write_text(dumps(raw))
    study,definition,record = snapshot(source,tmp_path/'campaign',{'operation.frequency_hz':value*1.25})
    exact = dict(study.fixed_parameters,**record['physical'])
    design = definition.adapter.build(exact)
    output = tmp_path/'next/study.toml'
    report = output.with_suffix('.refit.json')
    from dada_solver.research.cli import main
    assert main(['motion-refit',str(tmp_path/'campaign'),'--candidate',record['candidate_id'][:12],
                 '--output',str(output),'--report',str(report)]) == 0
    result = MotionRefitResult.from_data(json.loads(report.read_text()))
    new = load_study(output)
    assert new.data['schema_version'] == 3
    assert len(new.space.parameters) == 28
    assert all(p.name.startswith('kinematics.') for p in new.space.parameters)
    for name,value in exact.items():
        if not name.startswith('kinematics.'):
            assert new.fixed_parameters[name] == value
    for side in ('small','large'):
        assert new.settings[side] == dict(family='free_spline',representation='shape_coordinates',count=15)
        assert {p.name for p in new.space.parameters if p.name.startswith(f'kinematics.{side}.')} == {
            f'kinematics.{side}.shape_{i}' for i in range(13)} | {f'kinematics.{side}.phase_rad'}
    rebuilt = compile_study(new).adapter.build(dict(new.fixed_parameters,**{p.name:p.initial for p in new.space.parameters}))
    assert asdict(rebuilt.configuration) == asdict(design.configuration)
    assert asdict(rebuilt.heat_in) == asdict(design.heat_in)
    assert asdict(rebuilt.heat_out) == asdict(design.heat_out)
    assert new.data['constraints'] == study.data['constraints']
    assert new.data['numerical'] == study.data['numerical']
    assert new.data['study']['parent_candidate_id'] == record['candidate_id']
    for side in ('small','large'):
        region = new.basis.data['provenance']['motion_refit']['search_regions'][side]
        assert region['canonical_angular_radius_rad'] == .02
        assert region['box_inside_cap']
        assert region['phase_radius_rad'] == pytest.approx(math.radians(.48))
        for p in new.space.parameters:
            if p.name == f'kinematics.{side}.phase_rad':
                assert p.upper-p.lower == pytest.approx(2*math.radians(.48))
    assert new.data['search']['evaluate_initial']
    assert result.data['scientific']['target']['scientific']['source']['candidate_id'] == record['candidate_id']
    from dada_solver.research.cli import main
    assert main(['validate',str(output)]) == 0
    with pytest.raises(ValueError,match='must not exist'):
        refit_study(source,output)


def test_reference_pressure_source_freezes_exact_candidate_inventory(tmp_path,monkeypatch):
    no_thermodynamics(monkeypatch)
    source = initialize_kinematics(tmp_path/'source.toml','hybrid_compact','hybrid_compact')
    raw = tomllib.loads(source.read_text())
    from tests.test_research_reference_charge import derived
    derived(raw)
    source.write_text(dumps(raw))
    original = load_study(source)
    original_design = compile_study(original).adapter.build(original.fixed_parameters)
    output = tmp_path/'next.toml'
    refit_study(source,output,policy=RefitPolicy(phase_samples=8,dense_samples=360))
    new = load_study(output)
    assert new.data['policies']['charge'] == 'explicit_inventory'
    assert 'charge_reference' not in new.data
    assert new.fixed_parameters['charge.total_mass_kg'] == original_design.configuration.charge.total_mass
    assert new.basis.data['provenance']['motion_refit']['inventory_policy_change']['source_policy'] == original.data['policies']['charge']
    values = dict(new.fixed_parameters,**{p.name:p.initial for p in new.space.parameters})
    values['kinematics.small.phase_rad'] += STEP/4
    changed = compile_study(new).adapter.build(values)
    assert changed.configuration.charge.total_mass == original_design.configuration.charge.total_mass


def test_cli_fit_report_and_plot_without_machine_invention(tmp_path,monkeypatch,spline_fit):
    no_thermodynamics(monkeypatch)
    import dada_solver.research.motion_refit as refit_module
    _,target,result = spline_fit
    # CLI plumbing reuses an already tested fit; no duplicate expensive search.
    monkeypatch.setattr(refit_module.MotionRefitRequest,'execute',lambda self:result)
    path = tmp_path/'target.json'; target.save(path)
    from dada_solver.research.cli import main
    assert main(['motion-refit',str(path),'--validate-only']) == 0
    report = tmp_path/'report.json'
    assert main(['motion-refit',str(path),'--report',str(report)]) == 0
    assert MotionRefitResult.load(report).payload_json == result.payload_json
    with pytest.raises(SystemExit): main(['motion-refit',str(path),'--output',str(tmp_path/'study.toml')])
    assert not (tmp_path/'study.toml').exists()
    with pytest.raises(SystemExit): main(['motion-refit',str(path),'--report',str(report)])
    matplotlib = pytest.importorskip('matplotlib'); matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure = plot_refit(result)
    assert len(figure.axes) == 6
    for row in (0,3):
        nodes = next(c for c in figure.axes[row].collections if c.get_label() == '15 spline nodes')
        assert len(nodes.get_offsets()) == 15
    figure.canvas.draw(); plt.close(figure)
    assert main(['motion-refit',str(path),'--report',str(tmp_path/'second.json'),'--plot','--no-show']) == 0


def test_analytic_topology_rejects_extra_extrema():
    motion = _PeriodicMotion(FreeMotionDefinition(tuple(np.cos(3*np.arange(15)*STEP)),1.,2.))
    points = motion.stationary_points()
    assert len(points) == 6
    assert sum(p['kind'] == 'maximum' for p in points) == 3
    target = target_of_motion(motion)
    with pytest.raises(ValueError,match='topology-preserving'):
        MotionRefitRequest(target,policy=RefitPolicy(phase_samples=8,dense_samples=360)).execute()
