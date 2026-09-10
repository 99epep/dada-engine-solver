from dataclasses import replace
from pathlib import Path
import numpy as np
import pytest
from scipy.integrate import solve_ivp

from dada_solver.exchangers.air_wall import AirWallExchanger, AirWallMotor
from dada_solver.exchangers.wall_cycle import wall_cycle_performance
from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_model, build_initial_state


def exchanger(flow=.05):
    return AirWallExchanger(10, 20, 100, flow, 1005, 448.15)


def test_finite_air_capacity_and_energy_conservation():
    hx = exchanger()
    rates = hx.rates(300, 350*100)
    assert 350 < rates['air_outlet_temperature_k'] < 448.15
    assert rates['air_heat_w'] == pytest.approx(rates['wall_energy_rate_w']+rates['gas_heat_w'])
    assert rates['air_heat_w'] < .05*1005*(448.15-350)


def test_paused_isolated_gas_and_wall_relax_without_losing_energy():
    hx = exchanger(0)
    gas_capacity = 20
    def rhs(t, state):
        rates = hx.rates(state[0]/gas_capacity, state[1])
        return [rates['gas_heat_w'], rates['wall_energy_rate_w']]
    initial = [300*gas_capacity, 400*100]
    solution = solve_ivp(rhs, (0, 100), initial, rtol=1e-10, atol=1e-10)
    assert np.sum(solution.y[:, -1]) == pytest.approx(sum(initial), rel=1e-12)
    equilibrium = sum(initial)/(gas_capacity+100)
    assert solution.y[0, -1]/gas_capacity == pytest.approx(equilibrium, rel=1e-6)
    assert hx.rates(300, 40000)['air_outlet_temperature_k'] is None


def test_coupled_motor_instantaneous_total_energy_balance():
    config = load_simulation_configuration(Path(__file__).resolve().parents[1]/'examples/motor_demonstrator_piecewise.toml')
    model = build_model(config)
    wrapper = AirWallMotor(model, exchanger(), replace(exchanger(), air_inlet_temperature_k=298.15))
    state = np.r_[build_initial_state(config, model).as_array(), 40000, 30000, np.zeros(5)]
    derivative = wrapper.derivative(.4, state)*model.angular_speed
    stored_energy_rate = derivative[1:8:2].sum()+derivative[8:10].sum()
    assert stored_energy_rate == pytest.approx(derivative[10]+derivative[11]-derivative[14], abs=1e-10)
    assert derivative[:8:2].sum() == pytest.approx(0, abs=1e-15)


@pytest.mark.parametrize('field,value', [('wall_capacity_j_k',0),('air_mass_flow_kg_s',-1),('air_cp_j_kg_k',float('nan'))])
def test_invalid_inputs(field,value):
    with pytest.raises(ValueError):
        replace(exchanger(), **{field:value})


def test_steady_wall_matches_series_resistance_limit():
    hx = exchanger()
    capacity_rate = hx.air_mass_flow_kg_s*hx.air_cp_j_kg_k
    external = capacity_rate*(-np.expm1(-hx.air_wall_conductance_w_k/capacity_rate))
    gas_temperature = 310
    wall_temperature = (external*hx.air_inlet_temperature_k+hx.gas_wall_conductance_w_k*gas_temperature)/(external+hx.gas_wall_conductance_w_k)
    rates = hx.rates(gas_temperature, wall_temperature*hx.wall_capacity_j_k)
    expected = (hx.air_inlet_temperature_k-gas_temperature)/(1/external+1/hx.gas_wall_conductance_w_k)
    assert rates['wall_energy_rate_w'] == pytest.approx(0, abs=1e-10)
    assert rates['gas_heat_w'] == pytest.approx(expected)


def test_wall_efficiency_uses_external_air_heat_boundary():
    from types import SimpleNamespace
    trajectory=np.zeros((15,2))
    trajectory[:8,0]=trajectory[:8,1]=[1,100,1,100,1,100,1,100]
    trajectory[8:10,0]=[1000,1000];trajectory[8:10,1]=[1005,1005]
    trajectory[10:12,1]=[100,-70]
    trajectory[12:14,1]=[999,-999]  # Gas-wall quadratures must not set efficiency.
    trajectory[14,1]=20
    performance=wall_cycle_performance(SimpleNamespace(model=SimpleNamespace(signed_angular_speed=-2*np.pi)),trajectory)
    assert performance.thermal_efficiency==pytest.approx(.2)
    assert performance.gas_power==pytest.approx(20)
