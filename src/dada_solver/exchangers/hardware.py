"""Geometry-to-motor connection with explicit transport and material scenarios.

Laminar equivalent-passage laws are screening assumptions, not a Doty fit.
Half the tube/header loss is assigned to each side of the gas storage node.
"""
from dataclasses import dataclass, replace
import math

from dada_solver.exchangers.air_wall import AirWallExchanger, AirWallMotor
from dada_solver.exchangers.microtube_geometry import MicrotubeBank
from dada_solver.hydraulics import CompressibleOrifice, FlowResult


@dataclass(frozen=True)
class TubeHalfLink:
    bank: MicrotubeBank
    viscosity_pa_s: float
    core_loss_multiplier: float
    header_loss_coefficient: float
    valve_cda_m2: float | None = None

    def __post_init__(self):
        for value in (self.viscosity_pa_s, self.core_loss_multiplier):
            if not math.isfinite(value) or value <= 0:
                raise ValueError('Viscosity and loss multiplier must be positive.')
        if not math.isfinite(self.header_loss_coefficient) or self.header_loss_coefficient < 0:
            raise ValueError('Header loss coefficient must be nonnegative.')
        if self.valve_cda_m2 is not None and (not math.isfinite(self.valve_cda_m2) or self.valve_cda_m2 <= 0):
            raise ValueError('Valve CdA must be positive when supplied.')

    def directed_flow(self, upstream_pressure, downstream_pressure, upstream_temperature, gas):
        if any(not math.isfinite(v) or v <= 0 for v in (upstream_pressure, downstream_pressure, upstream_temperature)):
            raise ValueError('Pressures and temperature must be positive.')
        if downstream_pressure >= upstream_pressure:
            return FlowResult(0, False)
        rho = (upstream_pressure+downstream_pressure)/(2*gas.gas_constant*upstream_temperature)
        area = self.bank.dimensions()['tube_flow_area_m2']
        linear = self.core_loss_multiplier*64*self.viscosity_pa_s*self.bank.tube_length_m/(rho*self.bank.tube_count*math.pi*self.bank.inner_diameter_m**4)
        # Header K refers explicitly to total tube-passage velocity.
        quadratic = self.header_loss_coefficient/(4*rho*area**2)
        if self.valve_cda_m2 is not None:
            quadratic += 1/(2*rho*self.valve_cda_m2**2)
        dp = upstream_pressure-downstream_pressure
        flow = 2*dp/(linear+math.sqrt(linear*linear+4*quadratic*dp))
        cap = CompressibleOrifice(min(area, self.valve_cda_m2 or area)).directed_flow(
            upstream_pressure, downstream_pressure, upstream_temperature, gas)
        return FlowResult(min(flow, cap.mass_flow_rate), flow >= cap.mass_flow_rate and cap.is_choked)

    def bidirectional_flow(self, first_pressure, second_pressure, first_temperature, second_temperature, gas):
        if first_pressure >= second_pressure:
            return self.directed_flow(first_pressure, second_pressure, first_temperature, gas)
        result = self.directed_flow(second_pressure, first_pressure, second_temperature, gas)
        return FlowResult(-result.mass_flow_rate, result.is_choked)


@dataclass(frozen=True)
class HardwareInputs:
    metal_conductivity_w_m_k: float
    metal_density_kg_m3: float
    metal_cp_j_kg_k: float
    extra_wall_capacity_j_k: float
    gas_conductivity_w_m_k: float
    gas_viscosity_pa_s: float
    gas_nusselt: float
    air_conductivity_w_m_k: float
    air_viscosity_pa_s: float
    air_nusselt: float
    air_density_kg_m3: float
    air_cp_j_kg_k: float
    air_mass_flow_kg_s: float
    air_inlet_temperature_k: float
    air_poiseuille_number: float
    air_minor_loss_coefficient: float
    fan_total_efficiency: float
    core_loss_multiplier: float
    header_loss_coefficient: float

    def __post_init__(self):
        nonnegative = {'extra_wall_capacity_j_k', 'air_minor_loss_coefficient', 'header_loss_coefficient'}
        for key, value in vars(self).items():
            if not math.isfinite(value) or (value < 0 if key in nonnegative else value <= 0):
                raise ValueError(f'Invalid hardware input: {key}')
        if self.fan_total_efficiency > 1:
            raise ValueError('Total fan efficiency must not exceed one.')


def build_exchanger(bank: MicrotubeBank, inputs: HardwareInputs):
    dimensions = bank.dimensions()
    inner = bank.inner_diameter_m
    outer = inner+2*bank.wall_thickness_m
    area_out = bank.tube_count*math.pi*outer*bank.tube_length_m
    face = dimensions['core_width_m']*dimensions['core_height_m']
    free_area = face-bank.tube_count*math.pi*outer**2/4
    wetted = bank.tube_count*math.pi*outer+2*(dimensions['core_width_m']+dimensions['core_height_m'])
    air_diameter = 4*free_area/wetted
    gas_resistance = inner/(inputs.gas_nusselt*inputs.gas_conductivity_w_m_k*dimensions['tube_internal_area_m2'])
    air_resistance = air_diameter/(inputs.air_nusselt*inputs.air_conductivity_w_m_k*area_out)
    metal_resistance = math.log(outer/inner)/(2*math.pi*inputs.metal_conductivity_w_m_k*bank.tube_count*bank.tube_length_m)
    capacity = dimensions['wall_material_volume_m3']*inputs.metal_density_kg_m3*inputs.metal_cp_j_kg_k+inputs.extra_wall_capacity_j_k
    exchanger = AirWallExchanger(1/(gas_resistance+metal_resistance/2),
        1/(air_resistance+metal_resistance/2), capacity, inputs.air_mass_flow_kg_s,
        inputs.air_cp_j_kg_k, inputs.air_inlet_temperature_k)
    velocity = inputs.air_mass_flow_kg_s/(inputs.air_density_kg_m3*free_area)
    pressure_drop = (inputs.air_poiseuille_number*inputs.air_viscosity_pa_s*bank.tube_length_m*velocity/(2*air_diameter**2)
                     +inputs.air_minor_loss_coefficient*inputs.air_density_kg_m3*velocity**2/2)
    report = dict(**dimensions, gas_film_resistance_k_w=gas_resistance,
        metal_resistance_k_w=metal_resistance, air_film_resistance_k_w=air_resistance,
        wall_capacity_j_k=capacity, overall_static_conductance_w_k=1/(gas_resistance+metal_resistance+air_resistance),
        air_hydraulic_diameter_m=air_diameter, air_velocity_m_s=velocity,
        air_reynolds=inputs.air_density_kg_m3*velocity*air_diameter/inputs.air_viscosity_pa_s,
        air_pressure_drop_pa=pressure_drop,
        fan_electrical_power_w=pressure_drop*inputs.air_mass_flow_kg_s/inputs.air_density_kg_m3/inputs.fan_total_efficiency,
        assumptions='constant_properties_and_Nusselt; equivalent_laminar_air_passage; not_Doty_calibrated')
    return exchanger, report


def connect_hardware(model, heat_in_bank, heat_out_bank, heat_in_inputs, heat_out_inputs, *, heat_in_valve_cda_m2, heat_out_valve_cda_m2):
    # Compatibility constructor; geometry-specific assembly stays at this boundary.
    from dada_solver.exchangers.microtube import MicrotubeExchanger
    from dada_solver.exchangers.base import connect_exchangers
    incoming = MicrotubeExchanger(heat_in_bank, heat_in_inputs, heat_in_valve_cda_m2)
    outgoing = MicrotubeExchanger(heat_out_bank, heat_out_inputs, heat_out_valve_cda_m2)
    hi_report = dict(incoming.build().metadata)
    ho_report = dict(outgoing.build().metadata)
    return connect_exchangers(model, incoming, outgoing), dict(H_i=hi_report, H_o=ho_report,
        total_fan_electrical_power_w=hi_report['fan_electrical_power_w']+ho_report['fan_electrical_power_w'])
