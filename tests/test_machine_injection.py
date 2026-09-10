from dataclasses import replace
from pathlib import Path
import numpy as np
import pytest
from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_model
from dada_solver.free_kinematics import FreeKinematics
from dada_solver.machine import MachineDesign
from dada_solver.kinematics import ReversedVolumeKinematics

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('family', ['free', 'four_bar'])
@pytest.mark.parametrize('motor', [False, True])
def test_direct_injection_skips_selector_and_reverses_once(family, motor):
    name = 'motor_free_kinematics.toml' if family == 'free' else 'motor_demonstrator_original_325c.toml'
    source = load_simulation_configuration(ROOT/'examples'/name)
    source = replace(source, angular_speed=abs(source.angular_speed))
    injected = FreeKinematics(source.free_kinematics) if family == 'free' else build_model(source).kinematics
    # A different, valid selector makes accidental rebuilding observable.
    target = replace(source, kinematics_type='harmonic_example', small_phase_offset_degrees=1,
                     angular_speed=-source.angular_speed if motor else source.angular_speed)
    result = MachineDesign(configuration=target, kinematics=injected).build()
    for theta in np.linspace(0, 2*np.pi, 31):
        sign = -1 if motor else 1
        assert result.kinematics.small_cylinder_volume(theta) == pytest.approx(injected.small_cylinder_volume(sign*theta), abs=1e-15)
        assert result.kinematics.large_cylinder_volume_derivative(theta) == pytest.approx(sign*injected.large_cylinder_volume_derivative(sign*theta), abs=1e-15)
    if not motor:
        assert result.kinematics is injected


def test_already_reversed_and_mismatched_volumes_are_rejected():
    config = load_simulation_configuration(ROOT/'examples/motor_free_kinematics.toml')
    free = FreeKinematics(config.free_kinematics)
    with pytest.raises(ValueError, match='already reversed'):
        MachineDesign(config, kinematics=ReversedVolumeKinematics(free)).build()
    wrong = FreeKinematics(replace(config.free_kinematics, small=replace(config.free_kinematics.small, maximum_volume=.1)))
    with pytest.raises(ValueError, match='volume limits'):
        MachineDesign(config, kinematics=wrong).build()
