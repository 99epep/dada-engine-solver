from dataclasses import replace
from pathlib import Path
import tomllib
import pytest

from dada_solver.exchangers.hardware import HardwareInputs, TubeHalfLink, build_exchanger, connect_hardware
from dada_solver.exchangers.microtube import MicrotubeExchanger
from dada_solver.exchangers.microtube_geometry import MicrotubeBank
from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_model

ROOT = Path(__file__).resolve().parents[1]


def inputs():
    data = tomllib.loads((ROOT/'examples/motor_hardware.toml').read_text())
    return MicrotubeBank(**data['geometry']), HardwareInputs(**data['properties'], **data['heat_in'])


def test_half_tube_hydraulics_matches_full_poiseuille(ideal_gas):
    bank, props = inputs()
    link = TubeHalfLink(bank, props.gas_viscosity_pa_s, 1, 0)
    result = link.directed_flow(200001, 199999, 350, ideal_gas)
    rho = 200000/(ideal_gas.gas_constant*350)
    full = bank.laminar_tube_loss(result.mass_flow_rate, density_kg_m3=rho, viscosity_pa_s=props.gas_viscosity_pa_s)
    assert full['signed_tube_pressure_drop_pa'] == pytest.approx(4)
    assert link.bidirectional_flow(199999,200001,350,350,ideal_gas).mass_flow_rate == -result.mass_flow_rate
    assert link.directed_flow(199999,200001,350,ideal_gas).mass_flow_rate == 0
    assert replace(link, header_loss_coefficient=100).directed_flow(200001,199999,350,ideal_gas).mass_flow_rate < result.mass_flow_rate


def test_resistance_capacity_and_fan_energy():
    bank, props = inputs()
    hx, report = build_exchanger(bank, props)
    assert 1/hx.gas_wall_conductance_w_k+1/hx.air_wall_conductance_w_k == pytest.approx(1/report['overall_static_conductance_w_k'])
    assert hx.wall_capacity_j_k == pytest.approx(report['wall_material_volume_m3']*8000*500+2)
    _, changed = build_exchanger(bank, replace(props, fan_total_efficiency=.25))
    assert changed['fan_electrical_power_w'] == pytest.approx(2*report['fan_electrical_power_w'])
    _, conductive = build_exchanger(bank, replace(props, metal_conductivity_w_m_k=150))
    assert conductive['metal_resistance_k_w'] == pytest.approx(report['metal_resistance_k_w']/10)


def test_hardware_replaces_hold_up_and_all_four_port_losses():
    bank, props = inputs()
    model = build_model(load_simulation_configuration(ROOT/'examples/motor_demonstrator_four_bar.toml'))
    coupled, report = connect_hardware(model, bank, bank, props, replace(props, air_inlet_temperature_k=298.15), heat_in_valve_cda_m2=.0001, heat_out_valve_cda_m2=.0002)
    assert coupled.model.machine_volumes.cold_heat_exchanger == bank.dimensions()['working_gas_volume_m3']
    assert coupled.model.kinematics is model.kinematics
    assert isinstance(coupled.model.small_cold_link, TubeHalfLink)
    assert coupled.model.cold_large_valve.flow_model.valve_cda_m2 == .0001
    assert coupled.model.hot_small_valve.flow_model.valve_cda_m2 == .0002


@pytest.mark.parametrize('placement', ['downstream', 'upstream'])
def test_microtube_valve_cda_moves_with_placement(placement):
    bank, props = inputs()
    components = MicrotubeExchanger(bank, props, .0001, placement).build()
    expected = (None, .0001) if placement == 'downstream' else (.0001, None)
    assert (components.inlet.valve_cda_m2, components.outlet.valve_cda_m2) == expected


@pytest.mark.parametrize('field,value',[('fan_total_efficiency',2),('gas_nusselt',0),('extra_wall_capacity_j_k',-1),('air_density_kg_m3',float('nan'))])
def test_invalid_inputs(field,value):
    _, props = inputs()
    with pytest.raises(ValueError):
        replace(props, **{field:value})
