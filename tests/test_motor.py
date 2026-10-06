"""Motor signs, geometry reversal, and a conservative periodic motor cycle."""
from tests.synthetic_machine import configuration as synthetic_configuration, configuration_data
from dada_solver.research.study_io import dumps

from dataclasses import astuple, replace
import math
from types import SimpleNamespace

import numpy as np
import pytest

from dada_solver.factory import build_initial_state, build_model
from dada_solver.performance import OperatingMode, calculate_cycle_performance
from dada_solver.sizing.constraints import MinimumCoolingPower, MinimumMotorPower
from dada_solver.sizing.design import DesignParameter, DesignPoint, DesignVariable, apply_design_point
from dada_solver.sizing.objectives import MaximizeMotorPower, MaximizeThermalEfficiency
from dada_solver.sizing.configuration import load_sizing_problem
from tests.test_periodic import create_static_cycle


@pytest.mark.parametrize("family", ["harmonic_example", "ideal_piecewise_linear", "published_e0_opposed", "published_f65_opposed", "shared_crank_rocker"])
def test_reversal_preserves_origin_hardware_and_reverses_volume_rates(family):
    config = synthetic_configuration(family)
    # Unequal UA values catch an accidental exchange of hardware as well.
    config = replace(config, cold_thermal_conductance=11.0, hot_thermal_conductance=23.0)
    forward = build_model(config)
    motor_config = replace(config, angular_speed=-config.angular_speed)
    motor = build_model(motor_config)
    assert astuple(motor.volumes(0.0)) == pytest.approx(astuple(forward.volumes(0.0)), rel=1e-12)
    assert build_initial_state(motor_config, motor).as_array() == pytest.approx(build_initial_state(config, forward).as_array(), rel=1e-12)
    assert motor.hot_small_valve == forward.hot_small_valve
    assert motor.cold_large_valve == forward.cold_large_valve
    assert motor.machine_volumes == forward.machine_volumes
    assert motor.cold_heat_transfer.conductance == 11.0
    assert motor.hot_heat_transfer.conductance == 23.0
    assert motor.cold_heat_transfer.reservoir_temperature == config.hot_reservoir_temperature
    assert motor.hot_heat_transfer.reservoir_temperature == config.cold_reservoir_temperature
    for angle in (0.12, 1.2, 3.3, 5.9):
        assert astuple(motor.volumes(angle)) == pytest.approx(astuple(forward.volumes(-angle)), rel=1e-12)
        assert motor.cylinder_volume_rates(angle) == pytest.approx(
            -np.array(forward.cylinder_volume_rates(-angle))
        )
        # Independent finite difference also checks the derivative chain rule.
        h = 1e-6
        numerical = np.array([
            (motor.kinematics.small_cylinder_volume(angle+h) - motor.kinematics.small_cylinder_volume(angle-h))/(2*h),
            (motor.kinematics.large_cylinder_volume(angle+h) - motor.kinematics.large_cylinder_volume(angle-h))/(2*h),
        ]) * abs(config.angular_speed)
        assert motor.cylinder_volume_rates(angle) == pytest.approx(numerical, rel=1e-6, abs=1e-10)


def test_reversed_piecewise_discontinuities_and_origin():
    config = synthetic_configuration('ideal_piecewise_linear', motor=True)
    model = build_model(config)
    assert model.volumes(0.0).large_cylinder == config.machine_volumes.large_cylinder.maximum
    assert model.kinematics.breakpoint_angles() == pytest.approx(
        np.radians([108.0, 180.0, 288.0])
    )


@pytest.mark.parametrize("speed", [0.0, float("nan"), float("inf"), -float("inf")])
def test_nonfinite_or_zero_speed_is_rejected(speed):
    config = synthetic_configuration('ideal_piecewise_linear', motor=True)
    with pytest.raises(ValueError, match="non-zero"):
        replace(config, angular_speed=speed)


@pytest.mark.parametrize("qi,qo,work,mode", [
    (100.0, -60.0, 40.0, OperatingMode.MOTOR),
    (100.0, -140.0, -40.0, OperatingMode.REFRIGERATION),
    (-100.0, 140.0, 40.0, OperatingMode.NON_REFRIGERATION),
    (0.0, 0.0, 0.0, OperatingMode.NON_REFRIGERATION),
])
def test_performance_uses_received_heat_and_positive_elapsed_time(ideal_gas, qi, qo, work, mode):
    _, state, cycle = create_static_cycle(ideal_gas)
    cycle = replace(cycle, cold_heat=np.array([7.0, 7.0+qi]),
                    hot_heat=np.array([9.0, 9.0+qo]), gas_work=np.array([3.0, 3.0+work]))
    speed = 2.0*math.pi if mode is OperatingMode.REFRIGERATION else -2.0*math.pi
    performance = calculate_cycle_performance(cycle, state, speed)
    assert performance.operating_mode is mode
    assert performance.gas_power == pytest.approx(work)
    assert performance.heat_in_power == qi
    assert performance.heat_out_power == qo
    assert performance.conservation.absolute_energy_residual == pytest.approx(0.0)
    if mode is OperatingMode.MOTOR:
        assert performance.thermal_efficiency == pytest.approx(0.4)
        assert performance.motor_power == 40.0
        assert performance.cooling_cop is None
        assert performance.heating_cop is None
        evaluation = SimpleNamespace(performance=performance)
        assert not MinimumCoolingPower(1.0).evaluate(evaluation).available
        assert MinimumMotorPower(30.0).evaluate(evaluation).satisfied
        assert MaximizeThermalEfficiency().evaluate(evaluation).value == pytest.approx(-0.4)
        assert MaximizeMotorPower().evaluate(evaluation).value == -40.0
    else:
        assert performance.thermal_efficiency is None
        assert performance.motor_power is None


def test_motor_speed_design_bounds_keep_rotation_direction():
    config = synthetic_configuration('ideal_piecewise_linear', motor=True)
    variable = DesignVariable(DesignParameter.ANGULAR_SPEED, -2.0, -0.5, -1.0)
    changed = apply_design_point(config, DesignPoint({variable.parameter: -1.5}))
    assert changed.angular_speed == -1.5
    for lower, upper in [(-1.0, 1.0), (-1.0, 0.0), (0.0, 1.0)]:
        with pytest.raises(ValueError, match="zero"):
            DesignVariable(DesignParameter.ANGULAR_SPEED, lower, upper, 0.0)


def test_failed_motor_is_not_reported_as_a_refrigerator(ideal_gas):
    _, state, cycle = create_static_cycle(ideal_gas)
    cycle = replace(cycle, cold_heat=np.array([0.0, 100.0]),
                    hot_heat=np.array([0.0, -140.0]), gas_work=np.array([0.0, -40.0]))
    performance = calculate_cycle_performance(cycle, state, -2.0*math.pi)
    assert performance.operating_mode is OperatingMode.NON_REFRIGERATION
    assert performance.cooling_cop is None
    assert performance.heating_cop is None
    assert performance.thermal_efficiency is None
    assert not MinimumCoolingPower(1.0).evaluate(SimpleNamespace(performance=performance)).available


@pytest.mark.parametrize("objective", ["maximize_thermal_efficiency", "maximize_motor_power"])
def test_motor_sizing_configuration_loads(objective, tmp_path):
    base = configuration_data(); base['operation']['angular_speed'] = -10.
    (tmp_path / 'base.toml').write_text(dumps(base))
    path = tmp_path / 'motor_sizing.toml'
    path.write_text(f'''
[problem]
base_configuration = "base.toml"
[[variables]]
parameter = "angular_speed"
lower_bound = -2.0
upper_bound = -0.5
initial_value = -1.0
[objective]
type = "{objective}"
[[constraints]]
type = "minimum_motor_power"
required_power = 100.0
[optimizer]
objective_scale = 1.0
unavailable_objective_penalty = 1000000.0
unavailable_constraint_margin = -1.0
maximum_iterations = 10
function_tolerance = 1e-6
[optimizer.constraint_scales]
minimum_motor_power = 100.0
''')
    loaded = load_sizing_problem(path)
    assert loaded.problem.objective.name == objective
    assert isinstance(loaded.problem.constraints[0], MinimumMotorPower)
    assert loaded.problem.evaluator.base_configuration.motor_operation
