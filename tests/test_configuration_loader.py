from pathlib import Path
import math

import numpy as np
import pytest

from dada_solver.configuration import ChargeConfiguration, load_simulation_configuration
from dada_solver.factory import build_initial_state, build_model, build_periodic_solver


EXAMPLE_PATH = (
    Path(__file__).parents[1] / "examples" / "harmonic_controlled_example.toml"
)


def test_complete_example_configuration_builds_model_and_uniform_charge() -> None:
    configuration = load_simulation_configuration(EXAMPLE_PATH)
    model = build_model(configuration)
    state = build_initial_state(configuration, model)

    assert configuration.example_data
    assert configuration.angular_speed == pytest.approx(10.0)
    assert model.volumes(0.0).large_cylinder == pytest.approx(
        configuration.machine_volumes.large_cylinder.maximum
    )
    assert state.total_mass == pytest.approx(
        configuration.charge.total_mass
        or configuration.charge.pressure
        * model.volumes(0.0).total
        / (configuration.gas.gas_constant * configuration.charge.temperature)
    )
    assert build_periodic_solver(configuration, model).maximum_cycles == 100


def test_total_mass_charge_resolves_equivalent_uniform_pressure() -> None:
    charge = ChargeConfiguration(temperature=300.0, total_mass=0.01)
    configuration = load_simulation_configuration(EXAMPLE_PATH)

    pressure = charge.resolved_pressure(configuration.gas, 8.0e-4)

    assert pressure == pytest.approx(
        0.01 * configuration.gas.gas_constant * 300.0 / 8.0e-4
    )


def test_charge_requires_exactly_one_inventory_definition() -> None:
    with pytest.raises(ValueError, match="exactly one"):
        ChargeConfiguration(temperature=300.0, pressure=2.0e5, total_mass=0.01)


def test_harmonic_kinematics_requires_explicit_example_marker(tmp_path: Path) -> None:
    text = EXAMPLE_PATH.read_text(encoding="utf-8").replace(
        "example_data = true", "example_data = false"
    )
    path = tmp_path / "not_marked.toml"
    path.write_text(text, encoding="utf-8")

    with pytest.raises(ValueError, match="example_data"):
        load_simulation_configuration(path)


def test_piecewise_linear_cooling_cell_configuration_builds() -> None:
    path = Path(__file__).parents[1] / "examples" / "cooling_cell_machine_exploratory.toml"
    configuration = load_simulation_configuration(path)
    model = build_model(configuration)

    assert configuration.kinematics_type == "ideal_piecewise_linear"
    assert configuration.small_lambda_target == pytest.approx(0.7)
    assert configuration.large_lambda_target == pytest.approx(0.7)
    assert configuration.cold_reservoir_temperature == pytest.approx(268.15)
    assert configuration.adiabatic_sector_fraction == pytest.approx(0.15)
    assert configuration.numerical.integration_method == "Radau"
    assert model.kinematics.breakpoint_angles() == pytest.approx(
        (0.3 * math.pi, math.pi, 1.3 * math.pi)
    )
    state = build_initial_state(configuration, model)
    theta_zero_volume = 0.000004 + 0.000808 + 0.00002 + 0.00002
    expected_mass = (
        200000.0 * theta_zero_volume
        / (configuration.gas.gas_constant * configuration.charge.temperature)
    )
    assert state.total_mass == pytest.approx(expected_mass)


def test_33_liter_exploratory_configuration_preserves_selected_geometry() -> None:
    path = Path(__file__).parents[1] / "examples" / "cooling_cell_33_liter_exploratory.toml"
    configuration = load_simulation_configuration(path)

    assert configuration.angular_speed == pytest.approx(2.0 * math.pi)
    assert configuration.machine_volumes.large_cylinder.swept == pytest.approx(0.033)
    assert configuration.machine_volumes.small_cylinder.swept == pytest.approx(0.0165)
    assert configuration.machine_volumes.cold_heat_exchanger == pytest.approx(0.00066)
    assert configuration.machine_volumes.hot_heat_exchanger == pytest.approx(0.00066)
    assert configuration.cold_thermal_conductance == pytest.approx(125.0)
    assert configuration.hydraulics.small_to_cold_cda == pytest.approx(0.00075)
    assert configuration.hydraulics.cold_to_large_valve_cda == pytest.approx(0.00150)


def test_published_f65_configuration_builds_shared_crank_kinematics(
    tmp_path: Path,
) -> None:
    source = (
        Path(__file__).parents[1]
        / "examples"
        / "cooling_cell_high_cop_candidate.toml"
    ).read_text(encoding="utf-8")
    start = source.index("[kinematics]")
    end = source.index("\n[validity]", start)
    source = source[:start] + """[kinematics]
type = "published_f65_opposed"
ground_distance = 0.1
connecting_rod_to_projected_stroke_ratio = 5.0
crank_angle_offset_degrees = 241.5
crank_direction = 1
""" + source[end:]
    path = tmp_path / "f65.toml"
    path.write_text(source, encoding="utf-8")

    configuration = load_simulation_configuration(path)
    model = build_model(configuration)

    assert configuration.kinematics_type == "published_f65_opposed"
    assert model.kinematics.crank_radius == pytest.approx(0.04496)
    small, large = model.kinematics.slider_states(0.63)
    np.testing.assert_allclose(small.crank_pin, large.crank_pin)


def test_published_e0_configuration_builds_shared_crank_kinematics() -> None:
    path = Path(__file__).parents[1] / "examples" / "cooling_cell_e0_reference.toml"
    configuration = load_simulation_configuration(path)
    model = build_model(configuration)

    assert configuration.kinematics_type == "published_e0_opposed"
    assert configuration.valve_model == "continuous_ideal_diode"
    assert configuration.heat_in_valve_placement == "downstream"
    assert configuration.heat_out_valve_placement == "downstream"
    assert configuration.numerical.integration_method == "LSODA"
    assert model.kinematics.crank_radius == pytest.approx(0.073)
    small, large = model.kinematics.slider_states(1.27)
    np.testing.assert_allclose(small.crank_pin, large.crank_pin)


def test_valve_placement_parsing_and_validation(tmp_path: Path) -> None:
    source = EXAMPLE_PATH.read_text(encoding="utf-8").replace(
        '[valves.hot_to_small]\n',
        '[valves]\nmodel = "continuous_ideal_diode"\nheat_in_placement = "upstream"\nheat_out_placement = "downstream"\n\n[valves.hot_to_small]\n'
    )
    path = tmp_path / 'placements.toml'
    path.write_text(source, encoding='utf-8')
    configuration = load_simulation_configuration(path)
    assert configuration.heat_in_valve_placement == 'upstream'
    assert configuration.heat_out_valve_placement == 'downstream'
    path.write_text(source.replace('heat_in_placement = "upstream"',
                                   'heat_in_placement = "sideways"'), encoding='utf-8')
    with pytest.raises(ValueError, match='placement'):
        load_simulation_configuration(path)
    path.write_text(source.replace('model = "continuous_ideal_diode"',
                                   'model = "discrete_hysteretic"'), encoding='utf-8')
    with pytest.raises(ValueError, match='unsupported for discrete_hysteretic'):
        load_simulation_configuration(path)


def test_general_shared_crank_configuration_builds_independent_loops(tmp_path) -> None:
    source_path = Path(__file__).parents[1] / "examples" / "cooling_cell_e0_reference.toml"
    source = source_path.read_text(encoding="utf-8")
    start = source.index("[kinematics]\n")
    end = source.index("\n[validity]", start)
    replacement = """[kinematics]
type = "shared_crank_rocker"
ground_distance = 0.1
crank_ratio = 0.4
crank_angle_offset_degrees = 304.375
crank_direction = -1

[kinematics.small_four_bar]
coupler_ratio = 1.0
rocker_ratio = 1.17
output_along_ratio = 0.653805
output_normal_ratio = -0.970072
assembly_branch = 1
slider_rod_ratio = 5.0

[kinematics.large_four_bar]
coupler_ratio = 1.02
rocker_ratio = 1.15
output_along_ratio = 0.653805
output_normal_ratio = 0.970072
assembly_branch = -1
slider_rod_ratio = 5.0
"""
    path = tmp_path / "general_four_bar.toml"
    path.write_text(source[:start] + replacement + source[end:], encoding="utf-8")

    configuration = load_simulation_configuration(path)
    model = build_model(configuration)

    assert configuration.shared_four_bar_design is not None
    assert model.kinematics.crank_radius == pytest.approx(0.04)
    assert model.kinematics.small_assembly.loop.coupler_length == pytest.approx(0.1)
    assert model.kinematics.large_assembly.loop.coupler_length == pytest.approx(0.102)
