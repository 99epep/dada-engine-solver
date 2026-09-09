import pytest

from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_initial_state, build_model
from dada_solver.similarity import (
    required_power_similarity_factor,
    scale_capacity_and_speed,
    scale_extensive_machine,
    scale_volume_at_constant_inventory,
)
from dada_solver.performance import ConservationReport, CyclePerformance, OperatingMode


def test_extensive_similarity_preserves_intensive_initial_state() -> None:
    configuration = load_simulation_configuration(
        "examples/cooling_cell_mechanical_candidate_001.toml"
    )
    factor = 3.5
    scaled = scale_extensive_machine(configuration, factor)
    model = build_model(configuration)
    scaled_model = build_model(scaled)
    state = build_initial_state(configuration, model)
    scaled_state = build_initial_state(scaled, scaled_model)

    assert scaled_state.total_mass == pytest.approx(state.total_mass * factor)
    assert scaled_state.temperatures(scaled.gas) == pytest.approx(
        state.temperatures(configuration.gas)
    )
    assert scaled_state.pressures(scaled.gas, scaled_model.volumes(0.0)) == pytest.approx(
        state.pressures(configuration.gas, model.volumes(0.0))
    )
    assert scaled.four_bar_ground_distance == pytest.approx(
        configuration.four_bar_ground_distance * factor ** (1.0 / 3.0)
    )


def test_similarity_rejects_non_positive_factor() -> None:
    configuration = load_simulation_configuration(
        "examples/cooling_cell_mechanical_candidate_001.toml"
    )
    with pytest.raises(ValueError, match="Capacity factor"):
        scale_extensive_machine(configuration, 0.0)


def test_speed_similarity_scales_rates_without_scaling_inventory() -> None:
    configuration = load_simulation_configuration(
        "examples/cooling_cell_mechanical_candidate_001.toml"
    )
    scaled = scale_capacity_and_speed(configuration, 1.0, 3.0)

    assert scaled.machine_volumes == configuration.machine_volumes
    assert scaled.charge == configuration.charge
    assert scaled.angular_speed == pytest.approx(configuration.angular_speed * 3.0)
    assert scaled.cold_thermal_conductance == pytest.approx(
        configuration.cold_thermal_conductance * 3.0
    )
    assert scaled.hydraulics.cold_to_large_valve_cda == pytest.approx(
        configuration.hydraulics.cold_to_large_valve_cda * 3.0
    )


def test_required_power_similarity_factor_uses_refrigerating_reference() -> None:
    performance = CyclePerformance(
        cold_heat_per_cycle=20.0,
        hot_heat_per_cycle=-30.0,
        gas_work_per_cycle=-10.0,
        mechanical_input_per_cycle=10.0,
        cooling_cop=2.0,
        heating_cop=3.0,
        cooling_power=20.0,
        heating_power=30.0,
        mechanical_input_power=10.0,
        operating_mode=OperatingMode.REFRIGERATION,
        conservation=ConservationReport(0.0, 0.0, 0.0, 0.0),
    )

    assert required_power_similarity_factor(performance, 100.0) == pytest.approx(5.0)


def test_pressure_volume_similarity_preserves_inventory_and_scales_pressure() -> None:
    configuration = load_simulation_configuration(
        "examples/domestic_freezer_helium_seed.toml"
    )
    scaled = scale_volume_at_constant_inventory(configuration, 0.2)
    model = build_model(configuration)
    scaled_model = build_model(scaled)
    state = build_initial_state(configuration, model)
    scaled_state = build_initial_state(scaled, scaled_model)

    assert scaled_state.total_mass == pytest.approx(state.total_mass)
    assert scaled.charge.pressure == pytest.approx(configuration.charge.pressure / 0.2)
    assert scaled.machine_volumes.large_cylinder.swept == pytest.approx(
        configuration.machine_volumes.large_cylinder.swept * 0.2
    )
    assert scaled.hydraulics.large_to_hot_cda == pytest.approx(
        configuration.hydraulics.large_to_hot_cda * 0.2
    )
    reference_pressures = state.pressures(configuration.gas, model.volumes(0.0))
    scaled_pressures = scaled_state.pressures(scaled.gas, scaled_model.volumes(0.0))
    assert scaled_pressures == pytest.approx(reference_pressures / 0.2)


def test_pressure_volume_similarity_requires_pressure_charge() -> None:
    configuration = load_simulation_configuration(
        "examples/cooling_cell_mechanical_candidate_001.toml"
    )

    with pytest.raises(ValueError, match="pressure-based charge"):
        scale_volume_at_constant_inventory(configuration, 0.2)
