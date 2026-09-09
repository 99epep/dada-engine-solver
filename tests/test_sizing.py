from dataclasses import replace
import math
from pathlib import Path
from types import SimpleNamespace

import pytest
from dada_solver.performance import OperatingMode

from dada_solver.configuration import load_simulation_configuration
from dada_solver.sizing.constraints import (
    MaximumMachNumber,
    MaximumPressure,
    MinimumCoolingPower,
    MaximumPistonGasForce,
    MaximumIsothermalityError,
    RequireNominalCycleTopology,
    RequireValidThermodynamicModel,
    CompleteCoolingTaskWithinTime,
    MaximumMechanicalInputPower,
)
from dada_solver.sizing.design import (
    DesignParameter,
    DesignPoint,
    DesignVariable,
    apply_design_point,
)
from dada_solver.sizing.objectives import (
    MaximizeCoolingCop,
    MinimizeChargingPressure,
    MinimizeTotalSweptVolume,
    MinimizeTotalUa,
)
from dada_solver.mechanics import PistonFaceAreas
from dada_solver.topology import CycleTopologyClassification
from dada_solver.validity import ValidityVerdict
from dada_solver.thermal_load_configuration import load_cooling_cell_scenario


EXAMPLE_PATH = (
    Path(__file__).parents[1] / "examples" / "harmonic_controlled_example.toml"
)
LOAD_PATH = Path(__file__).parents[1] / "examples" / "cooling_cell_reference_load.toml"


def test_design_point_changes_supported_parameters_without_mutating_base() -> None:
    base = load_simulation_configuration(EXAMPLE_PATH)
    point = DesignPoint(
        {
            DesignParameter.SMALL_CLEARANCE_VOLUME: 1.5e-4,
            DesignParameter.SMALL_SWEPT_VOLUME: 3.5e-4,
            DesignParameter.COLD_UA: 12.0,
            DesignParameter.COLD_RESERVOIR_TEMPERATURE: 265.0,
            DesignParameter.HOT_RESERVOIR_TEMPERATURE: 305.0,
            DesignParameter.LARGE_TO_HOT_CDA: 8.0e-6,
            DesignParameter.CHARGE_PRESSURE: 3.0e5,
            DesignParameter.ANGULAR_SPEED: 15.0,
            DesignParameter.SMALL_PHASE_OFFSET: 0.5 * math.pi,
        }
    )

    changed = apply_design_point(base, point)

    assert changed.machine_volumes.small_cylinder.minimum == pytest.approx(1.5e-4)
    assert changed.machine_volumes.small_cylinder.swept == pytest.approx(3.5e-4)
    assert changed.cold_thermal_conductance == pytest.approx(12.0)
    assert changed.cold_reservoir_temperature == pytest.approx(265.0)
    assert changed.hot_reservoir_temperature == pytest.approx(305.0)
    assert changed.hydraulics.large_to_hot_cda == pytest.approx(8.0e-6)
    assert changed.charge.pressure == pytest.approx(3.0e5)
    assert changed.angular_speed == pytest.approx(15.0)
    assert changed.small_phase_offset_degrees == pytest.approx(90.0)
    assert base.machine_volumes.small_cylinder.minimum == pytest.approx(1.0e-4)
    assert base.cold_thermal_conductance == pytest.approx(5.0)


def test_pressure_and_mass_cannot_both_be_design_variables() -> None:
    base = load_simulation_configuration(EXAMPLE_PATH)
    point = DesignPoint(
        {
            DesignParameter.CHARGE_PRESSURE: 2.0e5,
            DesignParameter.TOTAL_GAS_MASS: 0.002,
        }
    )

    with pytest.raises(ValueError, match="cannot both"):
        apply_design_point(base, point)


def test_clearance_ratio_scales_with_designed_swept_volume_and_excludes_hx() -> None:
    base = load_simulation_configuration(EXAMPLE_PATH)
    point = DesignPoint(
        {
            DesignParameter.SMALL_SWEPT_VOLUME: 4.0e-4,
            DesignParameter.SMALL_CLEARANCE_RATIO: 0.025,
            DesignParameter.LARGE_SWEPT_VOLUME: 8.0e-4,
            DesignParameter.LARGE_CLEARANCE_RATIO: 0.10,
        }
    )

    changed = apply_design_point(base, point)

    assert changed.machine_volumes.small_cylinder.minimum == pytest.approx(1.0e-5)
    assert changed.machine_volumes.large_cylinder.minimum == pytest.approx(8.0e-5)
    assert changed.machine_volumes.cold_heat_exchanger == pytest.approx(
        base.machine_volumes.cold_heat_exchanger
    )
    assert changed.machine_volumes.hot_heat_exchanger == pytest.approx(
        base.machine_volumes.hot_heat_exchanger
    )


def test_common_lambda_design_keeps_piecewise_targets_equal() -> None:
    path = Path(__file__).parents[1] / "examples" / "cooling_cell_machine_exploratory.toml"
    base = load_simulation_configuration(path)
    changed = apply_design_point(
        base,
        DesignPoint(
            {
                DesignParameter.COMMON_LAMBDA_TARGET: 0.7,
                DesignParameter.ADIABATIC_SECTOR_FRACTION: 0.15,
            }
        ),
    )

    assert changed.small_lambda_target == pytest.approx(0.7)
    assert changed.large_lambda_target == pytest.approx(0.7)
    assert changed.adiabatic_sector_fraction == pytest.approx(0.15)


def test_individual_lambda_design_decouples_piecewise_targets() -> None:
    base = load_simulation_configuration(
        Path(__file__).parents[1] / "examples" / "cooling_cell_machine_exploratory.toml"
    )

    changed = apply_design_point(
        base,
        DesignPoint(
            {
                DesignParameter.SMALL_LAMBDA_TARGET: 0.5,
                DesignParameter.LARGE_LAMBDA_TARGET: 0.8,
            }
        ),
    )

    assert changed.small_lambda_target == pytest.approx(0.5)
    assert changed.large_lambda_target == pytest.approx(0.8)


def test_common_and_individual_lambda_targets_are_mutually_exclusive() -> None:
    base = load_simulation_configuration(
        Path(__file__).parents[1] / "examples" / "cooling_cell_machine_exploratory.toml"
    )

    with pytest.raises(ValueError, match="cannot be combined"):
        apply_design_point(
            base,
            DesignPoint(
                {
                    DesignParameter.COMMON_LAMBDA_TARGET: 0.6,
                    DesignParameter.LARGE_LAMBDA_TARGET: 0.7,
                }
            ),
        )


def test_kinematic_design_parameters_reject_incompatible_law() -> None:
    base = load_simulation_configuration(EXAMPLE_PATH)
    with pytest.raises(ValueError, match="ideal_piecewise_linear"):
        apply_design_point(
            base,
            DesignPoint({DesignParameter.ADIABATIC_SECTOR_FRACTION: 0.2}),
        )


def test_absolute_clearance_and_ratio_cannot_both_be_designed() -> None:
    base = load_simulation_configuration(EXAMPLE_PATH)
    point = DesignPoint(
        {
            DesignParameter.SMALL_CLEARANCE_VOLUME: 1.0e-5,
            DesignParameter.SMALL_CLEARANCE_RATIO: 0.05,
        }
    )

    with pytest.raises(ValueError, match="cannot both"):
        apply_design_point(base, point)


def test_design_variable_validates_bounds() -> None:
    with pytest.raises(ValueError, match="within"):
        DesignVariable(
            DesignParameter.COLD_UA,
            lower_bound=1.0,
            upper_bound=10.0,
            initial_value=20.0,
        )


def test_zero_ua_is_an_explicitly_supported_design_limit() -> None:
    variable = DesignVariable(
        DesignParameter.COLD_UA,
        lower_bound=0.0,
        upper_bound=10.0,
        initial_value=0.0,
    )
    point = DesignPoint({DesignParameter.COLD_UA: 0.0})

    assert variable.lower_bound == 0.0
    assert point.values[DesignParameter.COLD_UA] == 0.0


def test_physical_constraint_limits_must_be_positive() -> None:
    with pytest.raises(ValueError, match="positive"):
        MaximumPressure(0.0)


def test_objectives_remain_interchangeable() -> None:
    configuration = load_simulation_configuration(EXAMPLE_PATH)
    evaluation = SimpleNamespace(
        configuration=configuration,
        performance=SimpleNamespace(cooling_cop=2.5),
    )

    assert MaximizeCoolingCop().evaluate(evaluation).value == pytest.approx(-2.5)
    assert MinimizeTotalUa().evaluate(evaluation).value == pytest.approx(10.0)
    swept = MinimizeTotalSweptVolume().evaluate(evaluation).value
    assert swept == pytest.approx(6.0e-4)


def test_charging_pressure_objective_uses_actual_theta_zero_filling_volume() -> None:
    configuration = load_simulation_configuration(EXAMPLE_PATH)
    mass_configuration = replace(
        configuration,
        charge=replace(configuration.charge, pressure=None, total_mass=0.002),
    )
    evaluation = SimpleNamespace(configuration=mass_configuration)
    filling_volume = 1.0e-4 + 6.0e-4 + 1.0e-4 + 1.0e-4
    expected = (
        0.002
        * configuration.gas.gas_constant
        * configuration.charge.temperature
        / filling_volume
    )

    value = MinimizeChargingPressure().evaluate(evaluation)

    assert value.value == pytest.approx(expected)


def test_constraints_use_positive_feasibility_margin_and_preserve_unavailable() -> None:
    evaluation = SimpleNamespace(
        performance=SimpleNamespace(cooling_power=1200.0, operating_mode=OperatingMode.REFRIGERATION),
        diagnostics=SimpleNamespace(
            pressure_extrema={
                "S": SimpleNamespace(maximum=4.0e5),
                "L": SimpleNamespace(maximum=2.0e5),
            },
        ),
        validity=SimpleNamespace(maximum_mach_number=None),
    )

    cooling = MinimumCoolingPower(1000.0).evaluate(evaluation)
    pressure = MaximumPressure(5.0e5).evaluate(evaluation)
    mach = MaximumMachNumber(0.2).evaluate(evaluation)

    assert cooling.satisfied and cooling.margin == pytest.approx(200.0)
    assert pressure.satisfied and pressure.margin == pytest.approx(1.0e5)
    assert not mach.available
    assert not mach.satisfied
    assert mach.margin is None

    force = MaximumPistonGasForce(
        limit=1000.0,
        piston_areas=PistonFaceAreas(small=2.0e-3, large=3.0e-3),
    ).evaluate(evaluation)
    assert force.satisfied
    assert force.margin == pytest.approx(200.0)


def test_validity_and_topology_constraints_are_explicit() -> None:
    evaluation = SimpleNamespace(
        diagnostics=SimpleNamespace(
            topology=SimpleNamespace(
                classification=CycleTopologyClassification.NON_NOMINAL
            )
        ),
        validity=SimpleNamespace(
            verdict=ValidityVerdict.INDETERMINATE,
            cold_isothermality_error=0.02,
            hot_isothermality_error=0.03,
        ),
    )

    isothermality = MaximumIsothermalityError(0.04).evaluate(evaluation)
    topology = RequireNominalCycleTopology().evaluate(evaluation)
    validity = RequireValidThermodynamicModel().evaluate(evaluation)

    assert isothermality.satisfied
    assert not topology.satisfied
    assert not validity.satisfied


def test_cooling_task_and_human_power_constraints_use_machine_results() -> None:
    scenario = load_cooling_cell_scenario(LOAD_PATH)
    evaluation = SimpleNamespace(
        performance=SimpleNamespace(
            cooling_power=350.0,
            operating_mode=OperatingMode.REFRIGERATION,
            mechanical_input_power=145.0,
        )
    )

    task = CompleteCoolingTaskWithinTime(scenario.reference_task).evaluate(evaluation)
    human = MaximumMechanicalInputPower(150.0).evaluate(evaluation)

    assert task.satisfied
    assert task.margin == pytest.approx(
        350.0 - scenario.reference_task.required_average_cooling_power
    )
    assert human.satisfied
    assert human.margin == pytest.approx(5.0)
