"""Parametric tube-bank geometry and fully developed laminar screening.

These calculations are not a Doty calibration or a thermal cycle closure.
Header dimensions and losses are explicit engineering inputs.
"""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class MicrotubeBank:
    tube_count: int
    tube_length_m: float
    inner_diameter_m: float
    wall_thickness_m: float
    pitch_m: float
    header_depth_m: float
    additional_internal_volume_m3: float = 0.0

    def __post_init__(self):
        positive = (self.tube_length_m, self.inner_diameter_m,
                    self.wall_thickness_m, self.pitch_m, self.header_depth_m)
        if (type(self.tube_count) is not int or self.tube_count < 1
                or any(not math.isfinite(v) or v <= 0 for v in positive)
                or not math.isfinite(self.additional_internal_volume_m3)
                or self.additional_internal_volume_m3 < 0
                or self.pitch_m <= self.inner_diameter_m + 2*self.wall_thickness_m):
            raise ValueError('Expected positive geometry, nonoverlapping tubes and nonnegative additional volume.')

    def dimensions(self) -> dict:
        """Square-pitch rectangular packing with two full-face gas plenums.

        Dimensions describe internal fluid spaces, not pressure-vessel walls.
        Unoccupied positions in the final row remain part of the envelope.
        """
        columns = math.ceil(math.sqrt(self.tube_count))
        rows = math.ceil(self.tube_count / columns)
        width, height = columns*self.pitch_m, rows*self.pitch_m
        flow_area = self.tube_count * math.pi*self.inner_diameter_m**2/4
        tubes = flow_area*self.tube_length_m
        headers = 2*width*height*self.header_depth_m
        return dict(
            tube_flow_area_m2=flow_area,
            tube_internal_area_m2=self.tube_count*math.pi*self.inner_diameter_m*self.tube_length_m,
            tube_gas_volume_m3=tubes, header_gas_volume_m3=headers,
            working_gas_volume_m3=tubes+headers+self.additional_internal_volume_m3,
            core_width_m=width, core_height_m=height,
            fluid_envelope_length_m=self.tube_length_m+2*self.header_depth_m,
            wall_material_volume_m3=self.tube_count*math.pi/4*((self.inner_diameter_m+2*self.wall_thickness_m)**2-self.inner_diameter_m**2)*self.tube_length_m,
        )

    def laminar_tube_loss(self, mass_flow_kg_s: float, *, density_kg_m3: float,
                          viscosity_pa_s: float) -> dict:
        """Signed Poiseuille loss; excludes headers, entry and pulsation effects.

        Returns Reynolds number for an independent applicability check. The
        caller must not use this laminar estimate as a transition correlation.
        """
        if (not math.isfinite(mass_flow_kg_s)
                or any(not math.isfinite(v) or v <= 0 for v in (density_kg_m3, viscosity_pa_s))):
            raise ValueError('Expected finite flow and positive density and viscosity.')
        reynolds = 4*abs(mass_flow_kg_s)/(self.tube_count*math.pi*self.inner_diameter_m*viscosity_pa_s)
        pressure_drop = 128*viscosity_pa_s*self.tube_length_m*mass_flow_kg_s/(density_kg_m3*self.tube_count*math.pi*self.inner_diameter_m**4)
        return dict(signed_tube_pressure_drop_pa=pressure_drop, reynolds_number=reynolds,
                    status='laminar_assumption_requires_verification')
