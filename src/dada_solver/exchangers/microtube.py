"""Current microtube family: explicit design in, derived components out."""
from dataclasses import dataclass
from dada_solver.exchangers.base import ExchangerComponents
from dada_solver.exchangers.hardware import HardwareInputs, TubeHalfLink, build_exchanger
from dada_solver.exchangers.microtube_geometry import MicrotubeBank


@dataclass(frozen=True, slots=True)
class MicrotubeExchanger:
    bank: MicrotubeBank
    inputs: HardwareInputs
    outlet_valve_cda_m2: float
    valve_placement: str = 'downstream'

    def __post_init__(self):
        if self.valve_placement not in {'upstream', 'downstream'}:
            raise ValueError("Valve placement must be 'upstream' or 'downstream'.")

    def build(self) -> ExchangerComponents:
        thermal, report = build_exchanger(self.bank, self.inputs)
        def passage(valve=None):
            return TubeHalfLink(self.bank, self.inputs.gas_viscosity_pa_s,
                self.inputs.core_loss_multiplier, self.inputs.header_loss_coefficient, valve,
                self.inputs.gas_model)
        inlet_valve = self.outlet_valve_cda_m2 if self.valve_placement == 'upstream' else None
        outlet_valve = self.outlet_valve_cda_m2 if self.valve_placement == 'downstream' else None
        return ExchangerComponents(report['working_gas_volume_m3'], passage(inlet_valve),
            passage(outlet_valve), wall_thermal=thermal,
            metadata=tuple(report.items()),
            validity_domain=(('variable_internal_transport', 'instantaneous_correlation_domain',
                'quasi_steady_not_pulse_calibrated') if self.inputs.gas_model is not None else
                ('constant_properties_and_Nusselt', 'laminar_internal_tubes',
                'equivalent_laminar_external_passage', 'pulse_and_distribution_unvalidated',
                'not_empirically_Doty_calibrated')))
