"""Geometry-to-motor connection with explicit transport and material scenarios.

Laminar equivalent-passage laws are screening assumptions, not a Doty fit.
Half the tube/header loss is assigned to each side of the gas storage node.
"""
from dataclasses import dataclass, replace
import math
import tomllib

from dada_solver.exchangers.air_wall import AirWallExchanger, AirWallMotor
from dada_solver.exchangers.microtube_geometry import MicrotubeBank
from dada_solver.hydraulics import CompressibleOrifice, FlowResult
from dada_solver.exchangers.gas_correlations import MicrotubeGasModel, MicrotubeDomainError, compressible_poiseuille, darcy_smooth


def load_hardware_definition(source: str, gas_heat_capacity_cp: float):
    """Parse a frozen microtube hardware TOML snapshot."""
    data = tomllib.loads(source)
    if 'air_sizing' in data:
        sizing = data['air_sizing']
        factors = [sizing['reference_peak_internal_mass_flow_kg_s'],
            sizing['expected_peak_margin'], sizing['capacity_rate_ratio_to_expected_peak'],
            sizing['working_gas_cp_j_kg_k']]
        if not all(math.isfinite(value) and value > 0 for value in factors):
            raise ValueError('Air sizing factors must be finite and positive.')
        if not math.isclose(sizing['working_gas_cp_j_kg_k'], gas_heat_capacity_cp):
            raise ValueError('Air sizing heat capacity does not match the working gas.')
        data['properties']['air_mass_flow_kg_s'] = math.prod(factors)/data['properties']['air_cp_j_kg_k']
    bank = MicrotubeBank(**data['geometry'])
    heat_in = HardwareInputs(**data['properties'], **data['heat_in'])
    heat_out = HardwareInputs(**data['properties'], **data['heat_out'])
    if 'gas_model' in data:
        from dada_solver.exchangers.gas_transport import DiluteGasTransport
        from dada_solver.exchangers.gas_correlations import GasSurfaceAccommodation, SecondOrderSlip
        settings = dict(data['gas_model'])
        mode = settings.pop('mode', 'variable_properties')
        if mode != 'variable_properties': raise ValueError('gas_model.mode must be variable_properties; omit the section for legacy.')
        transport = DiluteGasTransport(settings.pop('species', 'air'))
        accommodation = GasSurfaceAccommodation(**settings.pop('accommodation', {}))
        slip = settings.pop('slip', None)
        model = MicrotubeGasModel(transport, accommodation,
            SecondOrderSlip(**slip) if slip is not None else None, **settings)
        heat_in = replace(heat_in, gas_model=model)
        heat_out = replace(heat_out, gas_model=model)
    return data, bank, heat_in, heat_out


@dataclass(frozen=True)
class TubeHalfLink:
    bank: MicrotubeBank
    viscosity_pa_s: float
    core_loss_multiplier: float
    header_loss_coefficient: float
    valve_cda_m2: float | None = None
    gas_model: MicrotubeGasModel | None = None

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
        if self.gas_model is not None:
            return self._gas_flow(upstream_pressure, downstream_pressure, upstream_temperature, gas)
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

    def _gas_flow(self, pin, pout, temperature, gas):
        """Isothermal pressure-squared Poiseuille; unchanged header/valve network.

        Each link owns L/2 and half the total header K. The minor-loss terms
        remain mean-density closures, explicitly separate from tube friction.
        """
        from scipy.optimize import brentq
        tr = self.gas_model.transport
        if not math.isclose(tr.gas_constant, gas.gas_constant, rel_tol=.005):
            raise ValueError('Transport species does not match thermodynamic gas.')
        mu = tr.viscosity(temperature)
        length = self.bank.tube_length_m/2
        diameter = self.bank.inner_diameter_m
        area = self.bank.tube_count*math.pi*diameter**2/4
        rho = (pin+pout)/(2*gas.gas_constant*temperature)
        factor = self.gas_model.slip_factor(pin,pout,temperature,diameter)
        nominal = compressible_poiseuille(pin,pout,temperature,mu,length,diameter,
                                          self.bank.tube_count,gas.gas_constant)*factor/self.core_loss_multiplier
        linear = (pin-pout)/nominal
        quadratic = self.header_loss_coefficient/(4*rho*area**2)
        if self.valve_cda_m2 is not None: quadratic += 1/(2*rho*self.valve_cda_m2**2)
        dp = pin-pout
        flow = 2*dp/(linear+math.sqrt(linear**2+4*quadratic*dp))
        re = flow*diameter/(area*mu)
        if re >= 2300:
            # The transition bridge is only a root-bracketing device. A root in
            # transition is rejected below, never certified as a correlation.
            def residual(m):
                reynolds = m*diameter/(area*mu)
                if reynolds < 2300: loss = linear*m
                else:
                    f0 = 64/2300
                    f = (f0+(darcy_smooth(4000)-f0)*(reynolds-2300)/1700
                         if reynolds<4000 else darcy_smooth(reynolds))
                    loss = self.core_loss_multiplier*f*length/diameter*m*m/(2*rho*area*area)
                return loss+quadratic*m*m-dp
            upper = min(flow,5e6*area*mu/diameter)
            if residual(upper)<0: raise MicrotubeDomainError('Required flow exceeds turbulent domain.')
            flow = brentq(residual,0,upper,xtol=1e-15)
            re = flow*diameter/(area*mu)
            if 2300<=re<4000: raise MicrotubeDomainError('Transition flow has no validated hydraulic closure.')
        diagnostics = self.gas_model.diagnose(self.bank,flow,pin,pout,temperature)
        failures = [x for x in diagnostics.issues if x in ('high_mach','beyond_continuum_model')]
        if re>=4000:
            failures += [x for x in diagnostics.issues if x in ('large_relative_pressure_drop','thermal_slip_not_implemented')]
        if failures and self.gas_model.domain_policy=='reject': raise MicrotubeDomainError('; '.join(failures))
        cap = CompressibleOrifice(min(area,self.valve_cda_m2 or area)).directed_flow(pin,pout,temperature,gas)
        return FlowResult(min(flow,cap.mass_flow_rate),flow>=cap.mass_flow_rate and cap.is_choked)

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
    gas_model: MicrotubeGasModel | None = None

    def __post_init__(self):
        nonnegative = {'extra_wall_capacity_j_k', 'air_minor_loss_coefficient', 'header_loss_coefficient'}
        if self.gas_model is not None and not isinstance(self.gas_model, MicrotubeGasModel):
            raise ValueError('gas_model must be MicrotubeGasModel or None (legacy).')
        for key, value in vars(self).items():
            if key == 'gas_model': continue
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
    if inputs.gas_model is not None:
        from dada_solver.exchangers.gas_film import MicrotubeGasFilm
        exchanger = replace(exchanger, gas_film=MicrotubeGasFilm(
            bank, inputs.gas_model, metal_resistance/2))
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
    if inputs.gas_model is not None:
        report['assumptions'] = 'variable_internal_transport; quasi_steady_gas_film; lumped_wall; legacy_external_air_film; not_Doty_calibrated'
        report['static_conductance_role'] = 'Legacy reference only; dynamic internal conductance comes from instantaneous flow'
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
