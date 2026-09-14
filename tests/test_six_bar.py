"""Independent six-bar kinematics and synthesis/factory regression checks."""
from dataclasses import replace
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pytest

from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_model
from dada_solver.kinematics import KinematicsModel, ReversedVolumeKinematics
from dada_solver.six_bar import IndependentSixBarVolumeKinematics, load_six_bar_mechanism

ROOT = Path(__file__).resolve().parents[1]
S = ROOT/'outputs/small_sixbar_stage2f_r4_freeH_tightaxis.json'
L = ROOT/'outputs/large_sixbar_stageL1.json'


@pytest.fixture(scope='module')
def mechanisms():
    return (load_six_bar_mechanism(S, restart=0, second_branch=1),
            load_six_bar_mechanism(L))


@pytest.fixture(scope='module')
def pair(mechanisms):
    config = load_simulation_configuration(ROOT/'examples/motor_demonstrator_original_325c.toml')
    limits = config.machine_volumes
    return config, IndependentSixBarVolumeKinematics(*mechanisms, limits.small_cylinder, limits.large_cylinder)


@pytest.mark.parametrize('side', [0, 1])
def test_periodic_and_smooth_analytic_derivative(mechanisms, side):
    m = mechanisms[side]
    for theta in np.linspace(-7, 7, 81):
        x, dx = m.slider_position_and_derivative(theta)
        np.testing.assert_allclose(m.slider_position_and_derivative(theta+2*math.pi), (x, dx), atol=2e-13)
        h = 1e-5
        numeric = (m.slider_position_and_derivative(theta+h)[0]-m.slider_position_and_derivative(theta-h)[0])/(2*h)
        assert dx == pytest.approx(numeric, abs=2e-8)
    assert m.slider_position_and_derivative(0) == m.slider_position_and_derivative(2*math.pi)
    np.testing.assert_allclose(m.slider_position_and_derivative(-1e-9),
                               m.slider_position_and_derivative(1e-9), atol=2e-8, rtol=0)


@pytest.mark.parametrize('side', [0, 1])
def test_matches_original_synthesis_position_and_derivative(mechanisms, side):
    path = [S, L][side]
    data = json.loads(path.read_text())
    candidate = (next(r['best_dense'] for r in data['runs']
                      if r['restart'] == 0 and r['second_branch'] == 1)
                 if side == 0 else data['best'])
    script = ['search_small_sixbar_stage2f', 'search_large_sixbar_stageL1'][side]
    spec = importlib.util.spec_from_file_location(script, ROOT/'examples'/f'{script}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    m = mechanisms[side]
    theta = np.linspace(0, 2*math.pi, 1440, endpoint=False)
    slider, derivative = np.array([m.slider_position_and_derivative(t) for t in theta]).T
    # The historical script normalizes against this grid's sampled extrema.
    sampled_q = 1-(slider-slider.min())/np.ptp(slider)
    sampled_dq = -derivative/np.ptp(slider)
    args = [np.array([candidate['parameters'][n] for n in module.NAMES])]
    if side == 1:
        args.append(m.primary_branch)
    result = module.evaluate(*args, m.second_branch, theta, sampled_q, sampled_dq, maximum_eh=100)
    assert result['position_rms'] < 1e-13
    assert result['derivative_rms'] < 1e-13
    assert result['stroke_over_crank'] == pytest.approx(candidate['stroke_over_crank'], abs=1e-12)
    continuous_q = np.array([m.normalized_motion(t)[0] for t in theta])
    assert np.max(np.abs(continuous_q-sampled_q)) < 1.3e-6


@pytest.mark.parametrize('side', ['small', 'large'])
def test_volume_limits_and_combined_evaluation(pair, side):
    _, kin = pair
    m = getattr(kin, side)
    limits = getattr(kin, f'{side}_volume_limits')
    volume = getattr(kin, f'{side}_cylinder_volume')
    assert volume(m.minimum_angle) == limits.maximum
    assert volume(m.maximum_angle) == limits.minimum
    values = np.array([volume(t) for t in np.linspace(0, 2*math.pi, 10001)])
    assert np.all(values >= limits.minimum)
    assert np.all(values <= limits.maximum)
    assert getattr(kin, f'{side}_physical_stroke') is None
    assert isinstance(kin, KinematicsModel)
    assert kin.breakpoint_angles() == ()
    t = .32
    assert kin.cylinder_volumes_and_derivatives(t) == (
        kin.small_cylinder_volume(t), kin.large_cylinder_volume(t),
        kin.small_cylinder_volume_derivative(t), kin.large_cylinder_volume_derivative(t))


def test_motor_factory_reverses_injected_sixbar_exactly_once(pair):
    config, kin = pair
    model = build_model(config, kinematics=kin)
    assert isinstance(model.kinematics, ReversedVolumeKinematics)
    assert model.kinematics.forward is kin
    for side in ('small', 'large'):
        for t in np.linspace(0, 2*math.pi, 30):
            assert getattr(model.kinematics, f'{side}_cylinder_volume')(t) == getattr(kin, f'{side}_cylinder_volume')(-t)
            assert getattr(model.kinematics, f'{side}_cylinder_volume_derivative')(t) == -getattr(kin, f'{side}_cylinder_volume_derivative')(-t)


def test_explicit_selection_and_invalid_geometry(mechanisms):
    assert mechanisms[0].primary_branch == -1
    assert mechanisms[1].primary_branch == 1
    assert mechanisms[0].second_branch == 1
    assert mechanisms[1].second_branch == -1
    with pytest.raises(ValueError, match='both'):
        load_six_bar_mechanism(S, restart=0)
    with pytest.raises(ValueError, match='absent'):
        load_six_bar_mechanism(S, restart=999, second_branch=1)
    with pytest.raises(ValueError, match='branches'):
        replace(mechanisms[0], second_branch=0)
    with pytest.raises(ValueError, match='piston rod'):
        replace(mechanisms[0], piston_rod=.01)
    with pytest.raises(ValueError, match='finite'):
        replace(mechanisms[0], primary_phase=float('nan'))
