from pathlib import Path

from dada_solver.exchangers.configuration import load_exchanger_screening_configuration
from dada_solver.exchangers.coupled_configuration import load_coupled_exchanger_problem
from dada_solver.exchangers.optimization import ExchangerObjective


SCREENING_TOML = '''[metadata]
example_data = false
[working_gas]
gas_constant = 287.0
heat_capacity_cp = 1004.5
heat_capacity_cv = 717.5
[transport]
dynamic_viscosity = 0.000018
thermal_conductivity = 0.026
[operating_point]
absolute_mass_flow_rate = 0.0001
pressure = 100000.0
temperature = 290.0
[thermal_resistances]
wall_thickness = 0.001
wall_thermal_conductivity = 200.0
external_conductance = 100.0
[requirements]
minimum_overall_conductance = 5.0
maximum_pressure_drop = 10000.0
maximum_mach_number = 0.3
maximum_gas_volume = 0.001
maximum_pressure_drop_fraction = 0.1
[parallel_channel_grid]
channel_counts = [20, 40]
channel_widths = [0.01, 0.02]
channel_heights = [0.001, 0.002]
channel_lengths = [0.1, 0.2]
'''


def test_load_declared_exchanger_grid_and_properties(tmp_path):
    path = tmp_path / "screening.toml"
    path.write_text(SCREENING_TOML)
    loaded = load_exchanger_screening_configuration(path)

    assert not loaded.example_data
    assert loaded.channel_counts == (20, 40)
    assert loaded.channel_lengths == (0.1, 0.2)
    assert loaded.transport.dynamic_viscosity == 0.000018
    assert loaded.requirements.maximum_mach_number == 0.3
    assert loaded.requirements.maximum_pressure_drop_fraction == 0.1


def test_coupled_loader_resolves_relative_paths_and_independent_searches(tmp_path):
    (tmp_path / "machine.toml").write_bytes(
        (Path(__file__).parent / "data" / "sizing_machine.toml").read_bytes()
    )
    (tmp_path / "cold.toml").write_text(SCREENING_TOML)
    (tmp_path / "hot.toml").write_text(
        SCREENING_TOML.replace("external_conductance = 100.0", "external_conductance = 200.0")
    )
    path = tmp_path / "coupled.toml"
    path.write_text('''[metadata]
example_data = false
[problem]
simulation_configuration = "machine.toml"
cold_exchanger_configuration = "cold.toml"
hot_exchanger_configuration = "hot.toml"
[coupling]
maximum_iterations = 3
relative_tolerance = 0.001
flow_under_relaxation = 0.5
exchanger_maximum_iterations = 10
exchanger_function_tolerance = 0.000001
hydraulic_continuation_steps = 2
[objectives]
cold = "minimum_gas_volume"
hot = "minimum_pressure_drop"
[hydraulic_partition]
cold_inlet_core_fraction = 0.4
cold_inlet_collector_loss_coefficient = 0.0
cold_outlet_collector_loss_coefficient = 0.0
hot_inlet_core_fraction = 0.6
hot_inlet_collector_loss_coefficient = 0.0
hot_outlet_collector_loss_coefficient = 0.0
''')
    loaded = load_coupled_exchanger_problem(path)

    assert not loaded.example_data
    assert loaded.cold_search.objective is ExchangerObjective.MINIMUM_GAS_VOLUME
    assert loaded.hot_search.objective is ExchangerObjective.MINIMUM_PRESSURE_DROP
    assert loaded.cold_search.thermal_resistances.external_conductance == 100.0
    assert loaded.hot_search.thermal_resistances.external_conductance == 200.0
    assert loaded.cold_search.inlet_core_resistance_fraction == 0.4
    assert loaded.hot_search.inlet_core_resistance_fraction == 0.6
    assert loaded.settings.hydraulic_continuation_steps == 2
