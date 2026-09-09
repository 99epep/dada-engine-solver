"""Interchangeability and regression checks at construction boundaries."""
from dataclasses import dataclass, replace
import ast
import math
from pathlib import Path
import tomllib
import numpy as np
import pytest
from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_model, build_initial_state, initial_valve_topology
from dada_solver.four_bar import FourBarKinematics, shared_crank_rocker_kinematics
from dada_solver.kinematics import KinematicsModel
from dada_solver.exchangers.base import ExchangerComponents, ExchangerModel, connect_exchangers
from dada_solver.exchangers.microtube import MicrotubeExchanger
from dada_solver.exchangers.hardware import HardwareInputs, build_exchanger, TubeHalfLink, connect_hardware
from dada_solver.exchangers.microtube_geometry import MicrotubeBank
from dada_solver.exchangers.air_wall import AirWallMotor
from dada_solver.heat_transfer import ReservoirHeatTransfer
from dada_solver.hydraulics import CompressibleOrifice
from dada_solver.state import ThermodynamicState

ROOT = Path(__file__).resolve().parents[1]


def config():
    return load_simulation_configuration(ROOT/'examples/motor_demonstrator_original_325c.toml')


def test_four_bar_protocol_matches_legacy_trajectory_and_branches():
    c = config()
    legacy = shared_crank_rocker_kinematics(c.shared_four_bar_design,
        c.machine_volumes.small_cylinder, c.machine_volumes.large_cylinder,
        ground_distance=c.four_bar_ground_distance,
        crank_angle_offset=math.radians(c.four_bar_crank_angle_offset_degrees),
        crank_direction=-c.four_bar_crank_direction)
    current = build_model(c).kinematics
    assert isinstance(current, FourBarKinematics)
    assert isinstance(current, KinematicsModel)
    for theta in np.linspace(-2*math.pi, 4*math.pi, 361):
        np.testing.assert_array_equal(current.cylinder_volumes_and_derivatives(theta),
                                      legacy.cylinder_volumes_and_derivatives(theta))
        assert current.slider_states(theta) == legacy.slider_states(theta)
    assert current.small_physical_stroke == legacy.small_physical_stroke


@dataclass(frozen=True)
class StaticTestExchanger:
    """Test-only non-geometric family: no tube, material, wall or air parameters."""
    volume: float
    temperature: float

    def build(self):
        return ExchangerComponents(self.volume, CompressibleOrifice(1e-5, .03),
            CompressibleOrifice(2e-5, .03), ReservoirHeatTransfer(3, self.temperature),
            validity_domain=('test_fixture_only',))


def test_non_microtube_static_family_runs_through_generic_connector():
    c = config();base = build_model(c)
    incoming = StaticTestExchanger(3e-5, 598.15)
    assert isinstance(incoming, ExchangerModel)
    model = connect_exchangers(base, incoming, StaticTestExchanger(5e-5, 298.15))
    assert model.kinematics is base.kinematics
    assert model.machine_volumes.cold_heat_exchanger == 3e-5
    assert model.machine_volumes.hot_heat_exchanger == 5e-5
    initial = build_initial_state(c, model)
    derivative = model.evaluate(.3, initial, initial_valve_topology())
    assert np.all(np.isfinite(derivative.state_derivative))


@dataclass(frozen=True)
class NonlinearTestWall:
    """Deliberately lacks linear conductance or air-stream attributes."""
    wall_capacity_j_k: float = 100

    def rates(self, gas_temperature_k, wall_energy_j):
        difference = wall_energy_j/self.wall_capacity_j_k-gas_temperature_k
        gas_heat = .001*difference**3
        return dict(gas_heat_w=gas_heat, air_heat_w=0., wall_energy_rate_w=-gas_heat)


@dataclass(frozen=True)
class WallTestExchanger:
    def build(self):
        return ExchangerComponents(3e-5, CompressibleOrifice(1e-5, .03),
            CompressibleOrifice(2e-5, .03), wall_thermal=NonlinearTestWall())


def test_wall_wrapper_consumes_rates_without_air_or_linear_film_assumptions():
    c = config()
    wrapper = connect_exchangers(build_model(c), WallTestExchanger(), WallTestExchanger())
    state = np.r_[build_initial_state(c, wrapper.model).as_array(), 40000, 30000, np.zeros(5)]
    rates = wrapper.derivative(.3, state)*wrapper.model.angular_speed
    assert rates[1:8:2].sum()+rates[8:10].sum() == pytest.approx(-rates[14], abs=1e-9)
    assert rates[:8:2].sum() == pytest.approx(0, abs=1e-15)


def test_microtube_connector_and_dynamic_rhs_preserve_previous_equations():
    c = config();base = build_model(c)
    data = tomllib.loads((ROOT/'examples/motor_hardware_parallel_325c.toml').read_text())
    bank = MicrotubeBank(**data['geometry'])
    hi = HardwareInputs(**data['properties'], **data['heat_in'])
    ho = HardwareInputs(**data['properties'], **data['heat_out'])
    wrapper, report = connect_hardware(base, bank, bank, hi, ho,
        heat_in_valve_cda_m2=c.hydraulics.cold_to_large_valve_cda,
        heat_out_valve_cda_m2=c.hydraulics.hot_to_small_valve_cda)
    thermal, original_report = build_exchanger(bank, hi)
    assert report['H_i'] == original_report
    assert isinstance(MicrotubeExchanger(bank, hi, 1e-4), ExchangerModel)
    state = np.r_[build_initial_state(c, wrapper.model).as_array(), 450*thermal.wall_capacity_j_k,
        330*thermal.wall_capacity_j_k, np.zeros(5)]
    gas = ThermodynamicState.from_array(state[:8])
    for angle in np.linspace(0, 2*math.pi, 51):
        # Original linear-film formulation, before the generic rates adapter.
        old = replace(wrapper.model,
            cold_heat_transfer=ReservoirHeatTransfer(wrapper.heat_in.gas_wall_conductance_w_k,
                state[8]/wrapper.heat_in.wall_capacity_j_k),
            hot_heat_transfer=ReservoirHeatTransfer(wrapper.heat_out.gas_wall_conductance_w_k,
                state[9]/wrapper.heat_out.wall_capacity_j_k))
        expected = old.evaluate(angle, gas, initial_valve_topology()).state_derivative
        np.testing.assert_allclose(wrapper.derivative(angle, state)[:8]*old.angular_speed,
            expected, rtol=3e-15, atol=1e-12)


def test_generic_physics_has_no_concrete_kinematics_imports():
    for filename in ('dynamics.py', 'integration.py', 'periodic.py', 'exchangers/air_wall.py'):
        tree = ast.parse((ROOT/'src/dada_solver'/filename).read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                assert 'four_bar' not in (node.module or '')
                assert 'free_kinematics' not in (node.module or '')


def test_mixed_storage_is_explicitly_unsupported():
    with pytest.raises(ValueError, match='state-layout'):
        connect_exchangers(build_model(config()), WallTestExchanger(), StaticTestExchanger(1e-5, 300))


def test_machine_composes_independent_families_without_a_registry():
    from dada_solver.machine import MachineDesign
    c = load_simulation_configuration(ROOT/'examples/motor_free_kinematics.toml')
    model = MachineDesign(c, StaticTestExchanger(4e-5, 598.15),
                          StaticTestExchanger(6e-5, 298.15)).build()
    assert model.machine_volumes.cold_heat_exchanger == 4e-5
    assert isinstance(model.kinematics, KinematicsModel)
    assert model.kinematics.small_physical_stroke is None
    with pytest.raises(ValueError, match='both'):
        MachineDesign(c, StaticTestExchanger(4e-5, 598.15))


def test_generic_closure_without_reservoir_has_unavailable_optional_diagnostic(ideal_gas):
    from tests.test_periodic import create_static_cycle
    from dada_solver.heat_transfer import PrescribedHeatRate
    from dada_solver.validity import assess_cycle_validity
    from dada_solver.sizing.constraints import MaximumIsothermalityError
    from types import SimpleNamespace
    model, _, cycle = create_static_cycle(ideal_gas)
    model = replace(model, cold_heat_transfer=PrescribedHeatRate(0), hot_heat_transfer=PrescribedHeatRate(0))
    report = assess_cycle_validity(cycle, model, config().validity)
    assert report.cold_isothermality_error is None
    assert report.hot_isothermality_error is None
    assert not MaximumIsothermalityError(.1).evaluate(SimpleNamespace(validity=report)).available
