"""Free laws, analytic feasibility and generic thermodynamic integration."""
from dataclasses import asdict, replace
import json
import math
from pathlib import Path
import numpy as np
import pytest
from dada_solver.free_kinematics import FreeMotionDefinition, FreeKinematicsConfiguration, FreeKinematics
from dada_solver.kinematics import KinematicsModel, ReversedVolumeKinematics
from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_model, build_initial_state, build_periodic_solver, initial_valve_topology
from dada_solver.results import extract_cycle_diagnostics
from dada_solver.sizing.design import apply_design_point, DesignPoint, DesignParameter

ROOT = Path(__file__).resolve().parents[1]


def motion(**kwargs):
    return FreeMotionDefinition((0, .2, 1, .4, -.3, .1), 1e-5, 1e-3, **kwargs)


def free():
    return FreeKinematics(FreeKinematicsConfiguration(motion(),
        FreeMotionDefinition((1, .3, -.2, .2), 2e-5, 2e-3)))


def test_periodic_wrap_and_analytic_derivative_continuity():
    model = free()
    assert isinstance(model, KinematicsModel)
    for side in ('small', 'large'):
        for suffix in ('', '_derivative', '_second_derivative'):
            function = getattr(model, side+'_cylinder_volume'+suffix)
            assert function(0) == function(2*math.pi) == function(-2*math.pi)
            assert function(-1e-8) == pytest.approx(function(1e-8), abs=1e-10)
            angles = np.linspace(-8, 8, 101)
            assert function(angles) == pytest.approx([function(float(x)) for x in angles], abs=1e-16)
        x = .713
        h = 1e-5
        v = getattr(model, side+'_cylinder_volume')
        d = getattr(model, side+'_cylinder_volume_derivative')
        dd = getattr(model, side+'_cylinder_volume_second_derivative')
        assert d(x) == pytest.approx((v(x+h)-v(x-h))/(2*h), rel=1e-8)
        assert dd(x) == pytest.approx((d(x+h)-d(x-h))/(2*h), rel=1e-8)


def test_true_extrema_scaling_independence_and_no_clipping():
    model = free()
    for internal in (model._small, model._large):
        from scipy.interpolate import CubicSpline
        values = internal.definition.control_values
        spline = CubicSpline(np.linspace(0, 2*math.pi, len(values)+1), (*values, values[0]), bc_type='periodic')
        roots = spline.derivative().roots(extrapolate=False)
        roots = roots[np.isfinite(roots)]
        values = internal.evaluate(np.r_[0, roots])
        assert min(values) == pytest.approx(internal.definition.minimum_volume, abs=1e-15)
        assert max(values) == pytest.approx(internal.definition.maximum_volume, abs=1e-15)
        grid = np.linspace(0, 2*math.pi, 10001)
        assert np.min(internal.evaluate(grid)) >= internal.definition.minimum_volume-1e-15
        assert np.max(internal.evaluate(grid)) <= internal.definition.maximum_volume+1e-15
    changed = FreeKinematics(replace(model.configuration, small=replace(motion(), control_values=(1,0,0,.5))))
    grid = np.linspace(0, 6, 101)
    np.testing.assert_array_equal(changed.large_cylinder_volume(grid), model.large_cylinder_volume(grid))
    assert not np.allclose(changed.small_cylinder_volume(grid), model.small_cylinder_volume(grid))


def test_normalization_and_serialization_remove_affine_redundancy():
    a = motion()
    b = replace(a, control_values=tuple(3*x+7 for x in a.control_values))
    np.testing.assert_allclose(a.control_values, b.control_values, atol=2e-15)
    restored = FreeMotionDefinition(**json.loads(json.dumps(asdict(a))))
    np.testing.assert_allclose(restored.control_values, a.control_values, atol=2e-15)
    chart = FreeMotionDefinition.from_shape_coordinates([.1, -.2, .4, .7], 1e-5, 1e-3)
    assert len(chart.control_values) == 6
    assert sum(chart.control_values) == pytest.approx(0, abs=1e-15)
    assert np.linalg.norm(chart.control_values) == pytest.approx(1)
    m = FreeKinematics(FreeKinematicsConfiguration(a, a))
    n = FreeKinematics(FreeKinematicsConfiguration(restored, restored))
    np.testing.assert_allclose(m.small_cylinder_volume(np.arange(7)), n.small_cylinder_volume(np.arange(7)), atol=1e-17)


@pytest.mark.parametrize('values', [(1,1,1,1), (0,1,2), (0,1,float('nan'),2), ((0,1),(1,0))])
def test_invalid_controls(values):
    with pytest.raises(ValueError):
        replace(motion(), control_values=values)


@pytest.mark.parametrize('kwargs', [dict(minimum_volume=0), dict(maximum_volume=1e-6),
    dict(maximum_absolute_first_derivative=-1), dict(maximum_absolute_second_derivative=float('inf'))])
def test_invalid_limits(kwargs):
    with pytest.raises(ValueError):
        replace(motion(), **kwargs)


def test_derivative_constraint_extrema_and_explicit_infeasibility():
    m = free()
    d = m.diagnostics[0]
    grid = np.unique(np.r_[np.linspace(0, 2*math.pi, 100001), np.linspace(0, 2*math.pi, 7)])
    assert d.maximum_absolute_first_derivative == pytest.approx(np.max(np.abs(m.small_cylinder_volume_derivative(grid))), rel=1e-7)
    assert d.maximum_absolute_second_derivative == pytest.approx(np.max(np.abs(m.small_cylinder_volume_second_derivative(grid))), rel=1e-7)
    invalid = FreeKinematics(replace(m.configuration, small=replace(m.configuration.small,
        maximum_absolute_first_derivative=d.maximum_absolute_first_derivative/2)))
    assert invalid.diagnostics[0].first_derivative_margin < 0
    with pytest.raises(ValueError, match='derivative limits'):
        invalid.require_feasible()
    np.testing.assert_array_equal(m.small_cylinder_volume(grid), invalid.small_cylinder_volume(grid))


def test_motor_reversal_and_sizing_volume_update():
    config = load_simulation_configuration(ROOT/'examples/motor_free_kinematics.toml')
    model = build_model(config)
    assert isinstance(model.kinematics, ReversedVolumeKinematics)
    forward = FreeKinematics(config.free_kinematics)
    assert model.kinematics.small_cylinder_volume(.7) == forward.small_cylinder_volume(-.7)
    assert model.kinematics.small_cylinder_volume_second_derivative(.7) == forward.small_cylinder_volume_second_derivative(-.7)
    changed = apply_design_point(config, DesignPoint({DesignParameter.SMALL_SWEPT_VOLUME: .002}))
    assert changed.free_kinematics.small.limits.swept == pytest.approx(.002)
    assert build_model(changed).kinematics.small_volume_limits == changed.machine_volumes.small_cylinder


@pytest.mark.parametrize('filename', ['motor_free_kinematics.toml', 'motor_demonstrator_original_325c.toml'])
def test_thermodynamic_cycle_and_periodic_diagnostics_through_interface(filename):
    config = load_simulation_configuration(ROOT/'examples'/filename)
    model = build_model(config)
    assert isinstance(model.kinematics, KinematicsModel)
    # One full cycle through the actual periodic evaluator, without claiming convergence.
    solver = replace(build_periodic_solver(config, model), maximum_cycles=1)
    initial = build_initial_state(config, model)
    result = solver.solve(initial, initial_valve_topology())
    assert result.final_cycle is not None and result.final_cycle.completed
    assert len(result.history) == 1
    assert np.all(np.isfinite(result.final_cycle.states))
    np.testing.assert_allclose(np.sum(result.final_cycle.states[::2], axis=0), initial.total_mass, atol=1e-12)
    assert extract_cycle_diagnostics(result.final_cycle, model) is not None


def test_sizing_preserves_explicit_kinematic_rejection_without_running_cycle():
    from dada_solver.sizing.evaluator import ThermodynamicSizingEvaluator, EvaluationStatus, evaluate_configuration
    c = load_simulation_configuration(ROOT/'examples/motor_free_kinematics.toml')
    c = replace(c, free_kinematics=replace(c.free_kinematics,
        small=replace(c.free_kinematics.small, maximum_absolute_first_derivative=0)))
    evaluator = ThermodynamicSizingEvaluator(c)
    result = evaluator.evaluate(DesignPoint({}))
    assert result.status is EvaluationStatus.INVALID_KINEMATICS
    assert not result.usable
    assert result.kinematic_diagnostics[0].first_derivative_margin < 0
    assert result.periodic.final_cycle is None
    assert evaluate_configuration(c).status is result.status
    assert evaluator.evaluate(DesignPoint({})) is result


def test_free_shape_survives_existing_similarity_transforms():
    from dada_solver.similarity import scale_volume_at_constant_inventory
    c = load_simulation_configuration(ROOT/'examples/motor_free_kinematics.toml')
    changed = scale_volume_at_constant_inventory(c, 2)
    assert changed.free_kinematics.small.control_values == c.free_kinematics.small.control_values
    assert changed.free_kinematics.small.limits.swept == pytest.approx(2*c.free_kinematics.small.limits.swept)
    assert build_model(changed).kinematics.small_cylinder_volume(.5) == pytest.approx(2*build_model(c).kinematics.small_cylinder_volume(.5))


def test_free_family_allows_multiple_extrema_without_phase_prescription():
    definition = FreeMotionDefinition((1,-1,1,-1,1,-1), 1e-5, 1e-3)
    model = FreeKinematics(FreeKinematicsConfiguration(definition, motion()))
    values = model.small_cylinder_volume(np.arange(6)*math.pi/3)
    np.testing.assert_allclose(values[::2], 1e-3, atol=1e-15)
    np.testing.assert_allclose(values[1::2], 1e-5, atol=1e-15)
    model.require_feasible()
