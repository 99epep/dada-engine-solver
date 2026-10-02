"""Microtube internal physics with a declared external thermal boundary."""
from dataclasses import dataclass, asdict
import math
from .external_stream import ExternalFluidStream, ExternalStreamWallExchanger
from .gas_film import MicrotubeGasFilm
from .gas_correlations import MicrotubeGasModel
from .hardware import TubeHalfLink
from .microtube_geometry import MicrotubeBank
from .base import ExchangerComponents


@dataclass(frozen=True)
class ExternalStreamHardwareInputs:
    metal_conductivity_w_m_k: float
    metal_density_kg_m3: float
    metal_cp_j_kg_k: float
    extra_wall_capacity_j_k: float
    gas_conductivity_w_m_k: float
    gas_viscosity_pa_s: float
    gas_nusselt: float
    core_loss_multiplier: float
    header_loss_coefficient: float
    external_stream: ExternalFluidStream
    gas_model: MicrotubeGasModel | None = None

    def __post_init__(self):
        if not isinstance(self.external_stream, ExternalFluidStream): raise ValueError('Expected external stream.')
        if self.gas_model is not None and not isinstance(self.gas_model, MicrotubeGasModel): raise ValueError('Unsupported gas film.')
        for key, value in vars(self).items():
            if key in ('external_stream','gas_model'): continue
            if not math.isfinite(value) or (value < 0 if key in ('extra_wall_capacity_j_k','header_loss_coefficient') else value <= 0):
                raise ValueError('Invalid internal hardware input: '+key)


@dataclass(frozen=True)
class ExternalStreamMicrotubeExchanger:
    bank: MicrotubeBank
    inputs: ExternalStreamHardwareInputs
    outlet_valve_cda_m2: float | None = None
    valve_placement: str = 'downstream'

    def __post_init__(self):
        if self.bank.circular_collectors:
            object.__setattr__(self, "outlet_valve_cda_m2", self.bank.dimensions()["conduit_area_m2"])
        if self.valve_placement not in ('upstream','downstream'): raise ValueError('Invalid valve placement.')

    def build(self):
        bank, inputs = self.bank, self.inputs
        dimensions = bank.dimensions()
        inner = bank.inner_diameter_m
        outer = inner + 2*bank.wall_thickness_m
        gas_resistance = inner/(inputs.gas_nusselt*inputs.gas_conductivity_w_m_k*dimensions['tube_internal_area_m2'])
        metal_resistance = math.log(outer/inner)/(2*math.pi*inputs.metal_conductivity_w_m_k*bank.tube_count*bank.tube_length_m)
        capacity = dimensions['wall_material_volume_m3']*inputs.metal_density_kg_m3*inputs.metal_cp_j_kg_k+inputs.extra_wall_capacity_j_k
        film = MicrotubeGasFilm(bank, inputs.gas_model, metal_resistance/2) if inputs.gas_model else None
        thermal = ExternalStreamWallExchanger(1/(gas_resistance+metal_resistance/2), capacity, inputs.external_stream, film)
        def passage(valve=None):
            if self.bank.circular_collectors: valve = None  # Ideal diode: direction is handled by the network.
            return TubeHalfLink(bank, inputs.gas_viscosity_pa_s, inputs.core_loss_multiplier,
                inputs.header_loss_coefficient, valve, inputs.gas_model)
        report = dict(dimensions, gas_film_resistance_k_w=gas_resistance,
            metal_resistance_k_w=metal_resistance, wall_capacity_j_k=capacity,
            external_stream=asdict(inputs.external_stream), external_loop_losses='excluded; no hydraulic or pump/fan model')
        if bank.circular_collectors:
            report.update(valve_model='ideal_diode_no_hydraulic_loss', valve_cda_m2=self.outlet_valve_cda_m2,
                header_loss_coefficient=inputs.header_loss_coefficient,
                header_loss_model='lumped_coefficient_based_on_total_tube_velocity; independent_of_frustum_geometry')
        return ExchangerComponents(dimensions['working_gas_volume_m3'],
            passage(self.outlet_valve_cda_m2 if self.valve_placement=='upstream' else None),
            passage(self.outlet_valve_cda_m2 if self.valve_placement=='downstream' else None),
            wall_thermal=thermal, metadata=tuple(report.items()),
            validity_domain=('declared_external_conductance', 'finite_capacity_stream',
                'variable_internal_transport' if film else 'constant_internal_properties',
                'external_hydraulics_unmodelled', 'quasi_steady_not_pulse_calibrated'))
