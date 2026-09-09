from pathlib import Path

from dada_solver.sizing.configuration import load_sizing_problem
from dada_solver.sizing.design import DesignParameter
from dada_solver.sizing.objectives import MinimizeTotalSweptVolume


EXAMPLE_PATH = (
    Path(__file__).parents[1] / "examples" / "sizing_controlled_example.toml"
)
COOLING_CELL_SPACE = (
    Path(__file__).parents[1] / "examples" / "cooling_cell_sensitivity_space.toml"
)


def test_complete_sizing_configuration_loads_without_running_simulation() -> None:
    loaded = load_sizing_problem(EXAMPLE_PATH)

    assert isinstance(loaded.problem.objective, MinimizeTotalSweptVolume)
    assert tuple(variable.parameter for variable in loaded.problem.variables) == (
        DesignParameter.SMALL_SWEPT_VOLUME,
        DesignParameter.LARGE_SWEPT_VOLUME,
        DesignParameter.CHARGE_PRESSURE,
    )
    assert len(loaded.problem.constraints) == 7
    assert set(loaded.optimization_settings.constraint_scales) == {
        constraint.name for constraint in loaded.problem.constraints
    }
    assert loaded.base_configuration_path.name == "harmonic_controlled_example.toml"
    assert loaded.cooling_cell_scenario is not None
    assert (
        loaded.cooling_cell_scenario.reference_task.required_average_cooling_power
        > 347.0
    )


def test_cooling_cell_space_preserves_user_selected_si_bounds() -> None:
    loaded = load_sizing_problem(COOLING_CELL_SPACE)
    variables = {item.parameter: item for item in loaded.problem.variables}

    assert variables[DesignParameter.ANGULAR_SPEED].lower_bound == 3.141592653589793
    assert variables[DesignParameter.ANGULAR_SPEED].upper_bound == 157.07963267948966
    assert variables[DesignParameter.ANGULAR_SPEED].initial_value == 5.0
    assert variables[DesignParameter.CHARGE_PRESSURE].lower_bound == 1.0e5
    assert variables[DesignParameter.CHARGE_PRESSURE].upper_bound == 2.0e6
    assert variables[DesignParameter.SMALL_SWEPT_VOLUME].lower_bound == 1.0e-5
    assert variables[DesignParameter.SMALL_SWEPT_VOLUME].upper_bound == 5.0e-3
    assert variables[DesignParameter.SMALL_SWEPT_VOLUME].initial_value == 4.0e-4
    assert variables[DesignParameter.LARGE_SWEPT_VOLUME].initial_value == 8.0e-4
    assert variables[DesignParameter.SMALL_CLEARANCE_RATIO].lower_bound == 0.005
    assert variables[DesignParameter.SMALL_CLEARANCE_RATIO].upper_bound == 0.10
    assert variables[DesignParameter.SMALL_CLEARANCE_RATIO].initial_value == 0.01
    assert variables[DesignParameter.LARGE_CLEARANCE_RATIO].initial_value == 0.01
    assert variables[DesignParameter.COLD_HEAT_EXCHANGER_VOLUME].lower_bound == 1.0e-6
    assert variables[DesignParameter.COLD_HEAT_EXCHANGER_VOLUME].upper_bound == 2.0e-3
    assert variables[DesignParameter.COLD_HEAT_EXCHANGER_VOLUME].initial_value == 2.0e-5
    assert variables[DesignParameter.HOT_HEAT_EXCHANGER_VOLUME].initial_value == 2.0e-5
    assert variables[DesignParameter.COLD_UA].lower_bound == 0.1
    assert variables[DesignParameter.COLD_UA].upper_bound == 500.0
    assert variables[DesignParameter.LARGE_TO_HOT_CDA].lower_bound == 1.0e-8
    assert variables[DesignParameter.LARGE_TO_HOT_CDA].upper_bound == 1.0e-3
    assert variables[DesignParameter.LARGE_TO_HOT_CDA].initial_value == 1.0e-5
    assert variables[DesignParameter.SMALL_TO_COLD_CDA].initial_value == 1.0e-5
    assert variables[DesignParameter.HOT_TO_SMALL_VALVE_CDA].initial_value == 4.0e-6
    assert variables[DesignParameter.COLD_TO_LARGE_VALVE_CDA].initial_value == 4.0e-6
    assert loaded.problem.evaluator.base_configuration.cold_reservoir_temperature == 268.15
    assert loaded.problem.evaluator.base_configuration.hot_reservoir_temperature == 298.15
