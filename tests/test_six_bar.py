"""Independent six-bar kinematics and synthesis/factory regression checks."""
from dataclasses import replace
import json
import math
from pathlib import Path

import numpy as np
import pytest

from dada_solver.factory import build_model
from dada_solver.kinematics import KinematicsModel, ReversedVolumeKinematics
from dada_solver.six_bar import IndependentSixBarVolumeKinematics, SixBarCylinderMechanism, load_six_bar_mechanism
from tests.synthetic_machine import configuration as synthetic_configuration

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def mechanisms():
    seeds=json.loads((ROOT/'src/dada_solver/research/data/kinematics_seeds.json').read_text())['seeds']
    return tuple(SixBarCylinderMechanism(**seeds[side]['six_bar']['parameters'])
                 for side in ('small','large'))


@pytest.fixture(scope='module')
def pair(mechanisms):
    config = synthetic_configuration(motor=True)
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


def test_explicit_selection_and_invalid_geometry(mechanisms,tmp_path):
    # The currently supported loader selects an explicit run without fallback.
    from dataclasses import asdict
    parameters=asdict(mechanisms[0])
    primary=parameters.pop('primary_branch'); secondary=parameters.pop('second_branch')
    candidate=dict(feasible=True,parameters=parameters,primary=dict(assembly_branch=primary),second_branch=secondary)
    path=tmp_path/'mechanism.json'
    path.write_text(json.dumps(dict(length_unit='crank_radius',best=candidate,
        runs=[dict(restart=3,second_branch=secondary,best_dense=candidate)])))
    restored=load_six_bar_mechanism(path,restart=3,second_branch=secondary)
    assert restored == mechanisms[0]
    with pytest.raises(ValueError, match='both'):
        load_six_bar_mechanism(path, restart=3)
    with pytest.raises(ValueError, match='absent'):
        load_six_bar_mechanism(path, restart=999, second_branch=secondary)
    with pytest.raises(ValueError, match='branches'):
        replace(mechanisms[0], second_branch=0)
    with pytest.raises(ValueError, match='piston rod'):
        replace(mechanisms[0], piston_rod=.01)
    with pytest.raises(ValueError, match='finite'):
        replace(mechanisms[0], primary_phase=float('nan'))
