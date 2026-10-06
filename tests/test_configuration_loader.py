from tests.synthetic_machine import CONFIGURATION, configuration_data
from pathlib import Path

import numpy as np
import pytest

from dada_solver.configuration import ChargeConfiguration, load_simulation_configuration
from dada_solver.factory import build_initial_state, build_model, build_periodic_solver


def test_complete_example_configuration_builds_model_and_uniform_charge() -> None:
    configuration = load_simulation_configuration(CONFIGURATION)
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
    assert build_periodic_solver(configuration, model).maximum_cycles == 10


def test_total_mass_charge_resolves_equivalent_uniform_pressure() -> None:
    charge = ChargeConfiguration(temperature=300.0, total_mass=0.01)
    configuration = load_simulation_configuration(CONFIGURATION)

    pressure = charge.resolved_pressure(configuration.gas, 8.0e-4)

    assert pressure == pytest.approx(
        0.01 * configuration.gas.gas_constant * 300.0 / 8.0e-4
    )


def test_charge_requires_exactly_one_inventory_definition() -> None:
    with pytest.raises(ValueError, match="exactly one"):
        ChargeConfiguration(temperature=300.0, pressure=2.0e5, total_mass=0.01)


def test_harmonic_kinematics_requires_explicit_example_marker(tmp_path: Path) -> None:
    text = CONFIGURATION.read_text(encoding="utf-8").replace(
        "example_data = true", "example_data = false"
    )
    path = tmp_path / "not_marked.toml"
    path.write_text(text, encoding="utf-8")

    with pytest.raises(ValueError, match="example_data"):
        load_simulation_configuration(path)


def test_published_f65_configuration_builds_shared_crank_kinematics(tmp_path):
    from dada_solver.research.study_io import dumps
    data = configuration_data()
    data['kinematics'] = dict(type='published_f65_opposed', ground_distance=.1,
        connecting_rod_to_projected_stroke_ratio=5., crank_angle_offset_degrees=0., crank_direction=1)
    path = tmp_path / 'f65.toml'; path.write_text(dumps(data))
    configuration = load_simulation_configuration(path)
    model = build_model(configuration)
    assert configuration.kinematics_type == 'published_f65_opposed'
    assert model.kinematics.crank_radius == pytest.approx(.04496)
    small, large = model.kinematics.slider_states(.63)
    np.testing.assert_allclose(small.crank_pin, large.crank_pin)


def test_published_e0_configuration_builds_shared_crank_kinematics(tmp_path):
    from dada_solver.research.study_io import dumps
    data = configuration_data()
    data['valves']['model'] = 'continuous_ideal_diode'
    data['kinematics'] = dict(type='published_e0_opposed', ground_distance=.1,
        connecting_rod_to_projected_stroke_ratio=5., crank_angle_offset_degrees=0., crank_direction=1)
    path = tmp_path / 'e0.toml'; path.write_text(dumps(data))
    configuration = load_simulation_configuration(path)
    model = build_model(configuration)
    assert configuration.kinematics_type == 'published_e0_opposed'
    assert configuration.valve_model == 'continuous_ideal_diode'
    assert model.kinematics.crank_radius == pytest.approx(.073)
    small, large = model.kinematics.slider_states(1.27)
    np.testing.assert_allclose(small.crank_pin, large.crank_pin)


def test_valve_placement_parsing_and_validation(tmp_path: Path) -> None:
    source = CONFIGURATION.read_text(encoding="utf-8").replace(
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
    source = CONFIGURATION.read_text()
    start = source.index('[kinematics]\n')
    end = len(source)
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


def test_validity_loader_uses_only_declared_keys(tmp_path):
    import tomllib
    from dada_solver.research.study_io import dumps

    data = tomllib.loads(CONFIGURATION.read_text())
    data["validity"] = {
        "maximum_compressibility_deviation": 0.01,
        "maximum_cp_variation": 0.02,
    }
    path = tmp_path / "configuration.toml"
    path.write_text(dumps(data))
    configuration = load_simulation_configuration(path)
    assert configuration.validity.maximum_cp_variation == 0.02
    data["validity"]["unexpected_validity_field"] = 0.03
    path.write_text(dumps(data))
    with pytest.raises(ValueError, match="Unknown validity keys"):
        load_simulation_configuration(path)
