import pytest

from dada_solver.performance import (
    ConservationReport,
    CyclePerformance,
    OperatingMode,
)
from dada_solver.thermal_loads import assess_cooling_task


def test_load_configuration_parses_declared_si_values(tmp_path):
    from dada_solver.research.study_io import dumps
    from dada_solver.thermal_load_configuration import load_cooling_cell_scenario
    path = tmp_path / 'load.toml'
    path.write_text(dumps(dict(
        thermal_load=dict(material_name='synthetic material', mass=2.,
            initial_temperature=300., target_temperature=290.,
            initial_phase_specific_heat=100., phase_change_specific_energy=1000.,
            phase_change_fraction=.5),
        targets=dict(reference_time=10., ambitious_time=5.),
        human_power=dict(minimum_power=50., nominal_power=150., maximum_power=200.),
        environment=dict(ambient_temperature=305.))))
    loaded = load_cooling_cell_scenario(path)
    assert loaded.load.sensible_energy == 2000.
    assert loaded.load.phase_change_energy == 1000.
    assert loaded.reference_task.allowed_time == 10.
    assert loaded.ambitious_task.allowed_time == 5.
    assert loaded.human_power.nominal_power == 150.
    assert loaded.ambient_temperature == 305.


def scenario():
    from dada_solver.thermal_loads import LumpedPhaseChangeLoad, CoolingTask, HumanPowerEnvelope
    from dada_solver.thermal_load_configuration import CoolingCellScenario
    load=LumpedPhaseChangeLoad('synthetic material',1.,300.,290.,100.,1000.)
    return CoolingCellScenario(load,CoolingTask(load,10.),CoolingTask(load,5.),
                               HumanPowerEnvelope(50.,150.,200.),300.)


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


def test_fixed_operating_point_reports_time_and_required_target_power() -> None:
    case = scenario()
    performance = refrigeration_performance(350.0, 150.0)

    assessment = assess_cooling_task(
        performance, case.reference_task, case.human_power
    )

    assert assessment.ideal_estimated_completion_time == pytest.approx(
        2000.0 / 350.0
    )
    assert assessment.required_pedaling_power_for_target_time == pytest.approx(
        case.reference_task.required_average_cooling_power
        / performance.cooling_cop
    )
    assert assessment.operable_at_nominal_human_power
    assert assessment.target_met_at_operating_point


def test_machine_requiring_excess_power_is_not_scaled_without_resimulation() -> None:
    case = scenario()
    performance = refrigeration_performance(500.0, 200.0)

    assessment = assess_cooling_task(
        performance, case.reference_task, case.human_power
    )

    assert assessment.ideal_estimated_completion_time == pytest.approx(4.0)
    assert not assessment.operable_at_nominal_human_power
    assert not assessment.target_met_at_operating_point
