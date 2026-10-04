"""Current microtube family: explicit design in, derived components out."""
from dataclasses import dataclass
from dada_solver.exchangers.base import ExchangerComponents
from dada_solver.exchangers.hardware import HardwareInputs, TubeHalfLink, build_exchanger
from dada_solver.exchangers.microtube_geometry import MicrotubeBank


@dataclass(frozen=True, slots=True)
class MicrotubeExchanger:
    bank: MicrotubeBank
    inputs: HardwareInputs
    outlet_valve_cda_m2: float | None = None
    valve_placement: str = 'downstream'

    def __post_init__(self):
        if self.bank.circular_collectors:
            object.__setattr__(self, "outlet_valve_cda_m2", self.bank.dimensions()["conduit_area_m2"])
        if self.valve_placement not in {'upstream', 'downstream'}:
            raise ValueError("Valve placement must be 'upstream' or 'downstream'.")

    def build(self) -> ExchangerComponents:
        thermal, report = build_exchanger(self.bank, self.inputs)
        if self.bank.circular_collectors:
            report.update(valve_model='ideal_diode_no_hydraulic_loss', valve_cda_m2=self.outlet_valve_cda_m2,
                header_loss_coefficient=self.inputs.header_loss_coefficient,
                header_loss_model='lumped_coefficient_based_on_total_tube_velocity; independent_of_frustum_geometry')
        def passage(half,valve=None):
            if self.bank.circular_collectors: valve = None  # Ideal diode: direction is handled by the network.
            return TubeHalfLink(self.bank, self.inputs.gas_viscosity_pa_s,
                self.inputs.core_loss_multiplier, self.inputs.header_loss_coefficient, valve,
                self.inputs.gas_model, half)
        if self.inputs.gas_model is not None:
            report.update(laminar_thermal_correlation=('bennett_2020_combined_entry_constant_wall'
                if self.inputs.gas_model.thermal_entry else 'fully_developed_3_66_screening'),
                laminar_hydraulic_correlation='compressible_poiseuille_plus_shah_london_1978_eq192',
                transition_thermal_correlation='bennett_gnielinski_transition_interpolation',
                transition_hydraulic_correlation='linear_darcy_complete_segment_endpoints_2300_4000',
                entrance_density='arithmetic_mean_pressure_at_upstream_temperature',
                entrance_segmentation='physical_halves; mirrored_on_reverse; no_midpoint_restart',
                header_loss_scope='manifold_contraction_exit; excludes_tube_profile_development')
        inlet_valve = self.outlet_valve_cda_m2 if self.valve_placement == 'upstream' else None
        outlet_valve = self.outlet_valve_cda_m2 if self.valve_placement == 'downstream' else None
        return ExchangerComponents(report['working_gas_volume_m3'], passage(0,inlet_valve),
            passage(1,outlet_valve), wall_thermal=thermal,
            metadata=tuple(report.items()),
            validity_domain=(('variable_internal_transport', 'instantaneous_correlation_domain',
                'quasi_steady_not_pulse_calibrated') if self.inputs.gas_model is not None else
                ('constant_properties_and_Nusselt', 'laminar_internal_tubes',
                'equivalent_laminar_external_passage', 'pulse_and_distribution_unvalidated',
                'not_empirically_Doty_calibrated')))
