from pathlib import Path

import pytest

from dada_solver.performance import (
    ConservationReport,
    CyclePerformance,
    OperatingMode,
)
from dada_solver.thermal_load_configuration import load_cooling_cell_scenario
from dada_solver.thermal_loads import assess_cooling_task


SCENARIO_PATH = (
    Path(__file__).parents[1] / "examples" / "cooling_cell_reference_load.toml"
)


def refrigeration_performance(cooling_power: float, mechanical_power: float):
    frequency = 1.0
    cooling_energy = cooling_power / frequency
    mechanical_energy = mechanical_power / frequency
    hot_energy = -(cooling_energy + mechanical_energy)
    return CyclePerformance(
        cold_heat_per_cycle=cooling_energy,
        hot_heat_per_cycle=hot_energy,
        gas_work_per_cycle=-mechanical_energy,
        mechanical_input_per_cycle=mechanical_energy,
        cooling_cop=cooling_power / mechanical_power,
        heating_cop=-hot_energy / mechanical_energy,
        cooling_power=cooling_power,
        heating_power=-hot_energy,
        mechanical_input_power=mechanical_power,
        operating_mode=OperatingMode.REFRIGERATION,
        conservation=ConservationReport(0.0, 0.0, 0.0, 0.0),
    )


def test_reference_water_load_matches_expected_order_of_magnitude() -> None:
    scenario = load_cooling_cell_scenario(SCENARIO_PATH)

    assert scenario.load.sensible_energy == pytest.approx(83_600.0)
    assert scenario.load.phase_change_energy == pytest.approx(333_500.0)
    assert scenario.load.ideal_cooling_energy == pytest.approx(417_100.0)
    assert scenario.reference_task.required_average_cooling_power == pytest.approx(
        347.5833333333333
    )
    assert scenario.ambitious_task.required_average_cooling_power == pytest.approx(
        463.44444444444446
    )
    assert scenario.reference_task.required_cop(150.0) == pytest.approx(
        2.317222222222222
    )
    assert scenario.ambitious_task.required_cop(150.0) == pytest.approx(
        3.0896296296296297
    )


def test_fixed_operating_point_reports_time_and_required_target_power() -> None:
    scenario = load_cooling_cell_scenario(SCENARIO_PATH)
    performance = refrigeration_performance(350.0, 150.0)

    assessment = assess_cooling_task(
        performance, scenario.reference_task, scenario.human_power
    )

    assert assessment.ideal_estimated_completion_time == pytest.approx(
        417_100.0 / 350.0
    )
    assert assessment.required_pedaling_power_for_target_time == pytest.approx(
        scenario.reference_task.required_average_cooling_power
        / performance.cooling_cop
    )
    assert assessment.operable_at_nominal_human_power
    assert assessment.target_met_at_operating_point


def test_machine_requiring_excess_power_is_not_scaled_without_resimulation() -> None:
    scenario = load_cooling_cell_scenario(SCENARIO_PATH)
    performance = refrigeration_performance(500.0, 200.0)

    assessment = assess_cooling_task(
        performance, scenario.reference_task, scenario.human_power
    )

    assert assessment.ideal_estimated_completion_time == pytest.approx(834.2)
    assert not assessment.operable_at_nominal_human_power
    assert not assessment.target_met_at_operating_point
