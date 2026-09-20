"""Internal variable-property film attached to the existing lumped wall."""
from dataclasses import dataclass
import math
from .gas_correlations import MicrotubeGasModel
from .microtube_geometry import MicrotubeBank


@dataclass(frozen=True)
class MicrotubeGasFilm:
    bank: MicrotubeBank
    model: MicrotubeGasModel
    half_wall_resistance_k_w: float

    def evaluate(self,temperature,wall_temperature,context):
        """Area-average two local port flow estimates, without net-flow cancellation.

        Each owns half the heat-transfer area. The thermal development length
        remains the physical full tube length, not a fictitious new entrance at
        the gas-storage node. This is a lumped approximation, not an axial solver.
        """
        conductance=0.; diagnostics=[]
        k=self.model.transport.conductivity(temperature)
        area=self.bank.tube_internal_area_m2
        for flow,p1,p2 in context['passages']:
            d=self.model.diagnose(self.bank,flow,p1,p2,temperature,frequency=context['frequency_hz'])
            self.model.require(d)
            conductance += .5*area*d.nusselt*k/self.bank.inner_diameter_m
            diagnostics.append(d)
        overall=1/(1/conductance+self.half_wall_resistance_k_w)
        return overall, tuple(diagnostics)
