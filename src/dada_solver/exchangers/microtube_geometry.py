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
    pitch_m: float | None = None
    header_depth_m: float | None = None
    additional_internal_volume_m3: float = 0.0
    pitch_ratio: float | None = None
    collector_half_angle_deg: float | None = None
    conduit_area_ratio: float = 1.0

    def __post_init__(self):
        positive = (self.tube_length_m, self.inner_diameter_m,
                    self.wall_thickness_m)
        if self.circular_collectors:
            if self.pitch_ratio is None or self.pitch_m is not None or self.header_depth_m is not None:
                raise ValueError("Circular collectors require pitch_ratio and replace pitch_m/header_depth_m.")
            if (not math.isfinite(self.collector_half_angle_deg) or not 0 < self.collector_half_angle_deg < 90):
                raise ValueError("collector_half_angle_deg must lie strictly between 0 and 90 degrees.")
        else:
            positive += (self.header_depth_m,)
            if self.conduit_area_ratio != 1.0:
                raise ValueError("conduit_area_ratio requires circular collectors.")
        if not math.isfinite(self.conduit_area_ratio) or self.conduit_area_ratio < 1:
            raise ValueError("conduit_area_ratio must be finite and at least 1.")
        if self.pitch_ratio is not None:
            if self.pitch_m is not None or not math.isfinite(self.pitch_ratio) or self.pitch_ratio <= 1:
                raise ValueError("pitch_ratio must exceed 1 and cannot be supplied with pitch_m.")
        elif self.pitch_m is None or not math.isfinite(self.pitch_m) or self.pitch_m <= self.outer_diameter_m:
            raise ValueError("Expected positive geometry, nonoverlapping tubes and nonnegative additional volume.")
        if (type(self.tube_count) is not int or self.tube_count < 1
                or any(v is None or not math.isfinite(v) or v <= 0 for v in positive)
                or not math.isfinite(self.additional_internal_volume_m3)
                or self.additional_internal_volume_m3 < 0):
            raise ValueError('Expected positive geometry, nonoverlapping tubes and nonnegative additional volume.')
        # Candidate-owned derived geometry; dataclass serialization stays input-only.
        dimensions = self._calculate_dimensions()
        object.__setattr__(self, '_prepared_dimensions', tuple(dimensions.items()))
        object.__setattr__(self, '_tube_flow_area_m2', dimensions['tube_flow_area_m2'])
        object.__setattr__(self, '_tube_internal_area_m2', dimensions['tube_internal_area_m2'])

    @property
    def circular_collectors(self):
        return self.collector_half_angle_deg is not None

    @property
    def outer_diameter_m(self):
        return self.inner_diameter_m + 2*self.wall_thickness_m

    @property
    def effective_pitch_m(self):
        return self.pitch_m if self.pitch_ratio is None else self.pitch_ratio*self.outer_diameter_m

    @property
    def tube_flow_area_m2(self):
        return self._tube_flow_area_m2

    @property
    def tube_internal_area_m2(self):
        return self._tube_internal_area_m2

    def dimensions(self) -> dict:
        """Return an independent reporting dictionary of prepared geometry."""
        return dict(self._prepared_dimensions)

    def _calculate_dimensions(self) -> dict:
        """Full-face internal gas plenums; no external interstitial volume.

        Ratio input selects staggered triangular packing bounded by tube edges.
        Absolute pitch selects the square envelope.
        """
        if self.circular_collectors:
            return self._circular_dimensions()
        columns = math.ceil(math.sqrt(self.tube_count))
        rows = math.ceil(self.tube_count / columns)
        if self.pitch_ratio is None:
            width, height = columns*self.pitch_m, rows*self.pitch_m
        else:
            pitch = self.effective_pitch_m
            # Only the last row can be incomplete. Earlier odd rows extend
            # half a pitch beyond full even rows; include actual occupied sites.
            last_count = self.tube_count - (rows-1)*columns
            right = last_count-1 + .5*((rows-1) % 2)
            if rows >= 2:
                right = max(right, columns-1)
            if rows >= 3:
                right = max(right, columns-.5)
            width = right*pitch + self.outer_diameter_m
            height = (rows-1)*math.sqrt(3)/2*pitch + self.outer_diameter_m
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

    def _circular_dimensions(self):
        """Continuous triangular-cell area, not a discrete packing guarantee.

        Both collectors are internal gas frusta. The external interstitial
        fluid volume is not part of the thermodynamic inventory.
        """
        pitch = self.effective_pitch_m
        face = self.tube_count*math.sqrt(3)/2*pitch**2
        bundle = math.sqrt(4*face/math.pi)
        flow_area = self.tube_count*math.pi*self.inner_diameter_m**2/4
        conduit_area = self.conduit_area_ratio*flow_area
        conduit = math.sqrt(4*conduit_area/math.pi)
        if conduit >= bundle:
            raise ValueError('Circular collector requires conduit_diameter_m < bundle_diameter_m; reduce conduit_area_ratio.')
        height = (bundle-conduit)/(2*math.tan(math.radians(self.collector_half_angle_deg)))
        one = math.pi*height/12*(bundle**2+bundle*conduit+conduit**2)
        tubes = flow_area*self.tube_length_m
        return dict(
            geometry_model='circular_triangular_frustum_v1',
            packing_method='continuous_triangular_cell_area; no_discrete_edge_correction',
            pitch_m=pitch, pitch_ratio=self.pitch_ratio, outer_diameter_m=self.outer_diameter_m,
            area_per_tube_m2=math.sqrt(3)/2*pitch**2,
            bundle_face_area_m2=face, bundle_diameter_m=bundle,
            conduit_area_m2=conduit_area, conduit_diameter_m=conduit,
            conduit_area_ratio=self.conduit_area_ratio,
            collector_half_angle_deg=self.collector_half_angle_deg,
            collector_height_m=height, single_collector_gas_volume_m3=one,
            tube_flow_area_m2=flow_area,
            tube_internal_area_m2=self.tube_count*math.pi*self.inner_diameter_m*self.tube_length_m,
            tube_gas_volume_m3=tubes, header_gas_volume_m3=2*one,
            additional_internal_volume_m3=self.additional_internal_volume_m3,
            working_gas_volume_m3=tubes+2*one+self.additional_internal_volume_m3,
            # Bounding-box extents only; never multiply these to obtain face area.
            core_width_m=bundle, core_height_m=bundle,
            fluid_envelope_length_m=self.tube_length_m+2*height,
            wall_material_volume_m3=self.tube_count*math.pi/4*(self.outer_diameter_m**2-self.inner_diameter_m**2)*self.tube_length_m)

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
