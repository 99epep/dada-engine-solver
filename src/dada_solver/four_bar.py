"""Planar shared-crank four-bar kinematics with an explicit slider output."""

from __future__ import annotations

from dataclasses import dataclass, field
import math

import numpy as np
from scipy.optimize import brentq

from dada_solver.geometry import CylinderVolumeLimits


def _positive(name: str, value: float) -> None:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive.")


def _unit_from_angle(angle: float) -> np.ndarray:
    return np.array((math.cos(angle), math.sin(angle)), dtype=float)


def _left_normal(vector: np.ndarray) -> np.ndarray:
    return np.array((-vector[1], vector[0]), dtype=float)


@dataclass(frozen=True, slots=True)
class FourBarLoop:
    """One coupler-rocker loop driven by the shared crank pin B."""

    coupler_length: float
    rocker_length: float
    rocker_pivot_x: float
    rocker_pivot_y: float
    assembly_branch: int = 1
    singularity_tolerance: float = 1.0e-10

    def __post_init__(self) -> None:
        _positive("Coupler length", self.coupler_length)
        _positive("Rocker length", self.rocker_length)
        for name, value in (
            ("Rocker pivot X", self.rocker_pivot_x),
            ("Rocker pivot Y", self.rocker_pivot_y),
        ):
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite.")
        if self.assembly_branch not in (-1, 1):
            raise ValueError("Assembly branch must be -1 or 1.")
        _positive("Singularity tolerance", self.singularity_tolerance)


@dataclass(frozen=True, slots=True)
class RockerOutputPoint:
    """Output point in the local frame from fixed pivot D toward joint C."""

    along_rocker: float
    normal_to_rocker: float = 0.0

    def __post_init__(self) -> None:
        if not math.isfinite(self.along_rocker) or not math.isfinite(
            self.normal_to_rocker
        ):
            raise ValueError("Rocker output coordinates must be finite.")


@dataclass(frozen=True, slots=True)
class CouplerOutputPoint:
    """Output point in the local frame from crank pin B toward joint C."""

    along_coupler: float
    normal_to_coupler: float = 0.0

    def __post_init__(self) -> None:
        if not math.isfinite(self.along_coupler) or not math.isfinite(
            self.normal_to_coupler
        ):
            raise ValueError("Coupler output coordinates must be finite.")


@dataclass(frozen=True, slots=True)
class SliderConstraint:
    """Finite rod from the rocker output point to a prismatic slider."""

    axis_origin_x: float
    axis_origin_y: float
    axis_angle: float
    connecting_rod_length: float
    assembly_branch: int = 1

    def __post_init__(self) -> None:
        for name, value in (
            ("Slider origin X", self.axis_origin_x),
            ("Slider origin Y", self.axis_origin_y),
            ("Slider axis angle", self.axis_angle),
        ):
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite.")
        _positive("Connecting-rod length", self.connecting_rod_length)
        if self.assembly_branch not in (-1, 1):
            raise ValueError("Slider assembly branch must be -1 or 1.")


@dataclass(frozen=True, slots=True)
class SliderState:
    """Position and angular derivative of one physical piston slider."""

    coordinate: float
    coordinate_derivative: float
    crank_pin: tuple[float, float]
    coupler_joint: tuple[float, float]
    output_point: tuple[float, float]
    four_bar_cross_product: float
    connecting_rod_transverse_margin: float


@dataclass(frozen=True, slots=True)
class SliderRodSensitivity:
    """Finite-rod effect relative to the output-point axial projection."""

    connecting_rod_to_projected_stroke_ratio: float
    projected_stroke: float
    slider_stroke: float
    maximum_normalized_shape_error: float


@dataclass(frozen=True, slots=True)
class SharedCrankRockerDesign:
    """Dimensionless two-loop rocker-output design around one shared crank."""

    crank_ratio: float
    small_coupler_ratio: float
    small_rocker_ratio: float
    small_output_along_ratio: float
    small_output_normal_ratio: float
    large_coupler_ratio: float
    large_rocker_ratio: float
    large_output_along_ratio: float
    large_output_normal_ratio: float
    small_pivot_x_ratio: float = 1.0
    small_pivot_y_ratio: float = 0.0
    large_pivot_x_ratio: float = 1.0
    large_pivot_y_ratio: float = 0.0
    small_four_bar_branch: int = 1
    large_four_bar_branch: int = -1
    small_slider_rod_ratio: float = 5.0
    large_slider_rod_ratio: float = 5.0

    def __post_init__(self) -> None:
        for name, value in (
            ("Crank ratio", self.crank_ratio),
            ("Small coupler ratio", self.small_coupler_ratio),
            ("Small rocker ratio", self.small_rocker_ratio),
            ("Large coupler ratio", self.large_coupler_ratio),
            ("Large rocker ratio", self.large_rocker_ratio),
            ("Small slider-rod ratio", self.small_slider_rod_ratio),
            ("Large slider-rod ratio", self.large_slider_rod_ratio),
        ):
            _positive(name, value)
        for value in (
            self.small_output_along_ratio,
            self.small_output_normal_ratio,
            self.large_output_along_ratio,
            self.large_output_normal_ratio,
            self.small_pivot_x_ratio,
            self.small_pivot_y_ratio,
            self.large_pivot_x_ratio,
            self.large_pivot_y_ratio,
        ):
            if not math.isfinite(value):
                raise ValueError("Four-bar design coordinates must be finite.")
        if self.small_four_bar_branch not in (-1, 1) or self.large_four_bar_branch not in (-1, 1):
            raise ValueError("Four-bar design branches must be -1 or 1.")


@dataclass(frozen=True, slots=True)
class FourBarSliderAssembly:
    """One four-bar and finite connecting rod driven by a supplied crank pin."""

    loop: FourBarLoop
    output: RockerOutputPoint | CouplerOutputPoint
    slider: SliderConstraint

    def evaluate(
        self,
        crank_pin: np.ndarray,
        crank_pin_derivative: np.ndarray,
    ) -> SliderState:
        pin_x, pin_y = float(crank_pin[0]), float(crank_pin[1])
        pin_dx, pin_dy = (
            float(crank_pin_derivative[0]),
            float(crank_pin_derivative[1]),
        )
        pivot_x = self.loop.rocker_pivot_x
        pivot_y = self.loop.rocker_pivot_y
        delta_x = pivot_x - pin_x
        delta_y = pivot_y - pin_y
        distance = math.hypot(delta_x, delta_y)
        scale = max(self.loop.coupler_length, self.loop.rocker_length, distance)
        tolerance = self.loop.singularity_tolerance * scale
        if distance <= tolerance:
            raise ValueError("Four-bar circle centers coincide or nearly coincide.")
        if distance > self.loop.coupler_length + self.loop.rocker_length + tolerance:
            raise ValueError("Four-bar loop cannot close: links are too short.")
        if distance < abs(self.loop.coupler_length - self.loop.rocker_length) - tolerance:
            raise ValueError("Four-bar loop cannot close: one link contains the other.")

        direction_x = delta_x / distance
        direction_y = delta_y / distance
        along = (
            self.loop.coupler_length**2
            - self.loop.rocker_length**2
            + distance**2
        ) / (2.0 * distance)
        height_squared = self.loop.coupler_length**2 - along**2
        if height_squared <= tolerance**2:
            raise ValueError("Four-bar loop is at or too close to a toggle singularity.")
        height = math.sqrt(height_squared)
        branch = self.loop.assembly_branch
        joint_x = pin_x + along * direction_x - branch * height * direction_y
        joint_y = pin_y + along * direction_y + branch * height * direction_x

        coupler_x = joint_x - pin_x
        coupler_y = joint_y - pin_y
        rocker_x = joint_x - pivot_x
        rocker_y = joint_y - pivot_y
        determinant = coupler_x * rocker_y - coupler_y * rocker_x
        if abs(determinant) <= tolerance * scale:
            raise ValueError("Four-bar velocity closure is singular.")
        velocity_rhs = coupler_x * pin_dx + coupler_y * pin_dy
        joint_dx = velocity_rhs * rocker_y / determinant
        joint_dy = -velocity_rhs * rocker_x / determinant

        if isinstance(self.output, RockerOutputPoint):
            output_unit_x = rocker_x / self.loop.rocker_length
            output_unit_y = rocker_y / self.loop.rocker_length
            output_unit_dx = joint_dx / self.loop.rocker_length
            output_unit_dy = joint_dy / self.loop.rocker_length
            output_origin_x, output_origin_y = pivot_x, pivot_y
            along = self.output.along_rocker
            normal_offset = self.output.normal_to_rocker
            output_origin_dx = output_origin_dy = 0.0
        else:
            output_unit_x = coupler_x / self.loop.coupler_length
            output_unit_y = coupler_y / self.loop.coupler_length
            output_unit_dx = (joint_dx - pin_dx) / self.loop.coupler_length
            output_unit_dy = (joint_dy - pin_dy) / self.loop.coupler_length
            output_origin_x, output_origin_y = pin_x, pin_y
            along = self.output.along_coupler
            normal_offset = self.output.normal_to_coupler
            output_origin_dx, output_origin_dy = pin_dx, pin_dy
        output_x = (
            output_origin_x + along * output_unit_x - normal_offset * output_unit_y
        )
        output_y = (
            output_origin_y + along * output_unit_y + normal_offset * output_unit_x
        )
        output_dx = (
            output_origin_dx
            + along * output_unit_dx
            - normal_offset * output_unit_dy
        )
        output_dy = (
            output_origin_dy
            + along * output_unit_dy
            + normal_offset * output_unit_dx
        )

        axis_x = math.cos(self.slider.axis_angle)
        axis_y = math.sin(self.slider.axis_angle)
        relative_x = output_x - self.slider.axis_origin_x
        relative_y = output_y - self.slider.axis_origin_y
        longitudinal = relative_x * axis_x + relative_y * axis_y
        transverse = -relative_x * axis_y + relative_y * axis_x
        transverse_derivative = -output_dx * axis_y + output_dy * axis_x
        longitudinal_derivative = output_dx * axis_x + output_dy * axis_y
        margin_squared = self.slider.connecting_rod_length**2 - transverse**2
        if margin_squared <= tolerance**2:
            raise ValueError("Slider connecting rod is at or beyond a toggle singularity.")
        margin = math.sqrt(margin_squared)
        slider_branch = self.slider.assembly_branch
        coordinate = longitudinal + slider_branch * margin
        coordinate_derivative = (
            longitudinal_derivative
            - slider_branch * transverse * transverse_derivative / margin
        )
        return SliderState(
            coordinate=coordinate,
            coordinate_derivative=coordinate_derivative,
            crank_pin=(pin_x, pin_y),
            coupler_joint=(joint_x, joint_y),
            output_point=(output_x, output_y),
            four_bar_cross_product=determinant,
            connecting_rod_transverse_margin=margin,
        )


@dataclass(frozen=True, slots=True)
class _NormalizedAssembly:
    assembly: FourBarSliderAssembly
    minimum_coordinate: float
    maximum_coordinate: float
    volume_increases_with_coordinate: bool

    @property
    def stroke(self) -> float:
        return self.maximum_coordinate - self.minimum_coordinate


@dataclass(frozen=True, slots=True)
class SharedCrankFourBarVolumeKinematics:
    """Two independent four-bar outputs driven by one physical crank pin."""

    crank_radius: float
    small_assembly: FourBarSliderAssembly
    large_assembly: FourBarSliderAssembly
    small_volume_limits: CylinderVolumeLimits
    large_volume_limits: CylinderVolumeLimits
    small_volume_increases_with_coordinate: bool = True
    large_volume_increases_with_coordinate: bool = True
    crank_angle_offset: float = 0.0
    crank_direction: int = 1
    extrema_scan_count: int = 721
    _small: _NormalizedAssembly = field(init=False, repr=False)
    _large: _NormalizedAssembly = field(init=False, repr=False)

    def __post_init__(self) -> None:
        _positive("Crank radius", self.crank_radius)
        if not math.isfinite(self.crank_angle_offset):
            raise ValueError("Crank angle offset must be finite.")
        if self.crank_direction not in (-1, 1):
            raise ValueError("Crank direction must be -1 or 1.")
        if self.extrema_scan_count < 16:
            raise ValueError("Extrema scan count must be at least 16.")
        object.__setattr__(
            self,
            "_small",
            self._normalize(
                self.small_assembly, self.small_volume_increases_with_coordinate
            ),
        )
        object.__setattr__(
            self,
            "_large",
            self._normalize(
                self.large_assembly, self.large_volume_increases_with_coordinate
            ),
        )

    def _crank(self, theta: float) -> tuple[np.ndarray, np.ndarray]:
        angle = self.crank_direction * theta + self.crank_angle_offset
        pin = self.crank_radius * np.array((math.cos(angle), math.sin(angle)))
        derivative = (
            self.crank_direction
            * self.crank_radius
            * np.array((-math.sin(angle), math.cos(angle)))
        )
        return pin, derivative

    def slider_states(self, theta: float) -> tuple[SliderState, SliderState]:
        pin, derivative = self._crank(theta)
        return (
            self.small_assembly.evaluate(pin, derivative),
            self.large_assembly.evaluate(pin, derivative),
        )

    def _normalize(
        self, assembly: FourBarSliderAssembly, increases: bool
    ) -> _NormalizedAssembly:
        angles = np.linspace(0.0, 2.0 * math.pi, self.extrema_scan_count, endpoint=False)
        values = np.empty(angles.size)
        derivatives = np.empty(angles.size)
        for index, angle in enumerate(angles):
            pin, pin_derivative = self._crank(float(angle))
            state = assembly.evaluate(pin, pin_derivative)
            values[index] = state.coordinate
            derivatives[index] = state.coordinate_derivative
        roots: list[float] = []
        for index, angle in enumerate(angles):
            next_index = (index + 1) % angles.size
            end = float(angles[next_index])
            if next_index == 0:
                end = 2.0 * math.pi
            left = derivatives[index]
            right = derivatives[next_index]
            if left == 0.0:
                roots.append(float(angle))
            elif left * right < 0.0:
                roots.append(
                    brentq(
                        lambda candidate: assembly.evaluate(
                            *self._crank(candidate)
                        ).coordinate_derivative,
                        float(angle),
                        end,
                    )
                    % (2.0 * math.pi)
                )
        extrema = [assembly.evaluate(*self._crank(angle)).coordinate for angle in roots]
        extrema.extend((float(np.min(values)), float(np.max(values))))
        minimum = min(extrema)
        maximum = max(extrema)
        if maximum - minimum <= np.finfo(float).eps * max(1.0, abs(minimum), abs(maximum)):
            raise ValueError("Four-bar slider has zero usable stroke.")
        return _NormalizedAssembly(assembly, minimum, maximum, increases)

    def _volume_and_derivative(
        self,
        theta: float,
        normalized: _NormalizedAssembly,
        limits: CylinderVolumeLimits,
    ) -> tuple[float, float]:
        state = normalized.assembly.evaluate(*self._crank(theta))
        fraction = (
            state.coordinate - normalized.minimum_coordinate
        ) / normalized.stroke
        derivative = state.coordinate_derivative / normalized.stroke
        if not normalized.volume_increases_with_coordinate:
            fraction = 1.0 - fraction
            derivative = -derivative
        swept = limits.swept
        return limits.minimum + swept * fraction, swept * derivative

    def cylinder_volumes_and_derivatives(
        self, theta: float
    ) -> tuple[float, float, float, float]:
        small_volume, small_derivative = self._volume_and_derivative(
            theta, self._small, self.small_volume_limits
        )
        large_volume, large_derivative = self._volume_and_derivative(
            theta, self._large, self.large_volume_limits
        )
        return small_volume, large_volume, small_derivative, large_derivative

    def small_cylinder_volume(self, theta: float) -> float:
        return self._volume_and_derivative(
            theta, self._small, self.small_volume_limits
        )[0]

    def large_cylinder_volume(self, theta: float) -> float:
        return self._volume_and_derivative(
            theta, self._large, self.large_volume_limits
        )[0]

    def small_cylinder_volume_derivative(self, theta: float) -> float:
        return self._volume_and_derivative(
            theta, self._small, self.small_volume_limits
        )[1]

    def large_cylinder_volume_derivative(self, theta: float) -> float:
        return self._volume_and_derivative(
            theta, self._large, self.large_volume_limits
        )[1]

    def breakpoint_angles(self) -> tuple[float, ...]:
        return ()

    @property
    def small_physical_stroke(self) -> float:
        return self._small.stroke

    @property
    def large_physical_stroke(self) -> float:
        return self._large.stroke


def shared_crank_rocker_kinematics(
    design: SharedCrankRockerDesign,
    small_volume_limits: CylinderVolumeLimits,
    large_volume_limits: CylinderVolumeLimits,
    *,
    ground_distance: float,
    crank_angle_offset: float = 0.0,
    crank_direction: int = 1,
    extrema_scan_count: int = 721,
) -> SharedCrankFourBarVolumeKinematics:
    """Build two independently parameterized rocker-output loops on one crank.

    Every geometric coordinate in ``design`` is normalized by the common
    ground-distance scale. The two horizontal slider axes are centered on the
    transverse excursions of their respective output points. This centering is
    a declared construction rule, not an optimized thermodynamic parameter.
    """

    _positive("Ground distance", ground_distance)
    crank_radius = design.crank_ratio * ground_distance

    def make_assembly(
        coupler_ratio: float,
        rocker_ratio: float,
        pivot_x_ratio: float,
        pivot_y_ratio: float,
        output_along_ratio: float,
        output_normal_ratio: float,
        four_bar_branch: int,
        slider_rod_ratio: float,
    ) -> FourBarSliderAssembly:
        loop = FourBarLoop(
            coupler_length=coupler_ratio * ground_distance,
            rocker_length=rocker_ratio * ground_distance,
            rocker_pivot_x=pivot_x_ratio * ground_distance,
            rocker_pivot_y=pivot_y_ratio * ground_distance,
            assembly_branch=four_bar_branch,
        )
        output = RockerOutputPoint(
            output_along_ratio * ground_distance,
            output_normal_ratio * ground_distance,
        )
        probe = FourBarSliderAssembly(
            loop,
            output,
            SliderConstraint(0.0, 0.0, 0.0, 100.0 * ground_distance),
        )
        minimum_x = minimum_y = math.inf
        maximum_x = maximum_y = -math.inf
        for angle in np.linspace(
            0.0, 2.0 * math.pi, extrema_scan_count, endpoint=False
        ):
            pin = crank_radius * np.array((math.cos(angle), math.sin(angle)))
            derivative = crank_radius * np.array((-math.sin(angle), math.cos(angle)))
            point_x, point_y = probe.evaluate(pin, derivative).output_point
            minimum_x = min(minimum_x, point_x)
            maximum_x = max(maximum_x, point_x)
            minimum_y = min(minimum_y, point_y)
            maximum_y = max(maximum_y, point_y)
        projected_stroke = maximum_x - minimum_x
        _positive("Projected output stroke", projected_stroke)
        return FourBarSliderAssembly(
            loop,
            output,
            SliderConstraint(
                axis_origin_x=0.0,
                axis_origin_y=0.5 * (minimum_y + maximum_y),
                axis_angle=0.0,
                connecting_rod_length=slider_rod_ratio * projected_stroke,
            ),
        )

    small = make_assembly(
        design.small_coupler_ratio,
        design.small_rocker_ratio,
        design.small_pivot_x_ratio,
        design.small_pivot_y_ratio,
        design.small_output_along_ratio,
        design.small_output_normal_ratio,
        design.small_four_bar_branch,
        design.small_slider_rod_ratio,
    )
    large = make_assembly(
        design.large_coupler_ratio,
        design.large_rocker_ratio,
        design.large_pivot_x_ratio,
        design.large_pivot_y_ratio,
        design.large_output_along_ratio,
        design.large_output_normal_ratio,
        design.large_four_bar_branch,
        design.large_slider_rod_ratio,
    )
    return SharedCrankFourBarVolumeKinematics(
        crank_radius=crank_radius,
        small_assembly=small,
        large_assembly=large,
        small_volume_limits=small_volume_limits,
        large_volume_limits=large_volume_limits,
        crank_angle_offset=crank_angle_offset,
        crank_direction=crank_direction,
        extrema_scan_count=extrema_scan_count,
    )


def published_f65_opposed_kinematics(
    small_volume_limits: CylinderVolumeLimits,
    large_volume_limits: CylinderVolumeLimits,
    *,
    ground_distance: float,
    connecting_rod_to_projected_stroke_ratio: float = 5.0,
    crank_angle_offset: float = 0.0,
    crank_direction: int = 1,
    extrema_scan_count: int = 721,
) -> SharedCrankFourBarVolumeKinematics:
    """Build the historically opposed F65 pair with parallel slider axes.

    The four-bar ratios and coupler point are published. The ground distance,
    finite-rod ratio, and centering of each slider axis are explicit modelling
    choices. Slider axes are parallel to fixed link AD and pass through the
    midpoint of the output point's transverse excursion.
    """

    _positive("Ground distance", ground_distance)
    _positive(
        "Connecting-rod to projected-stroke ratio",
        connecting_rod_to_projected_stroke_ratio,
    )
    scale = ground_distance / 100.0
    crank_radius = 44.96 * scale

    def make_assembly(four_bar_branch: int) -> FourBarSliderAssembly:
        normal_coordinate = -88.65 * scale * four_bar_branch
        loop = FourBarLoop(
            coupler_length=117.91 * scale,
            rocker_length=141.73 * scale,
            rocker_pivot_x=ground_distance,
            rocker_pivot_y=0.0,
            assembly_branch=four_bar_branch,
        )
        output = CouplerOutputPoint(94.4 * scale, normal_coordinate)
        probe = FourBarSliderAssembly(
            loop,
            output,
            SliderConstraint(0.0, 0.0, 0.0, 100.0 * ground_distance),
        )
        angles = np.linspace(0.0, 2.0 * math.pi, extrema_scan_count, endpoint=False)
        points = []
        for angle in angles:
            pin = crank_radius * np.array((math.cos(angle), math.sin(angle)))
            derivative = crank_radius * np.array((-math.sin(angle), math.cos(angle)))
            points.append(probe.evaluate(pin, derivative).output_point)
        point_array = np.asarray(points)
        projected_stroke = float(np.ptp(point_array[:, 0]))
        axis_y = 0.5 * (
            float(np.min(point_array[:, 1])) + float(np.max(point_array[:, 1]))
        )
        return FourBarSliderAssembly(
            loop,
            output,
            SliderConstraint(
                axis_origin_x=0.0,
                axis_origin_y=axis_y,
                axis_angle=0.0,
                connecting_rod_length=(
                    connecting_rod_to_projected_stroke_ratio * projected_stroke
                ),
            ),
        )

    return SharedCrankFourBarVolumeKinematics(
        crank_radius=crank_radius,
        small_assembly=make_assembly(1),
        large_assembly=make_assembly(-1),
        small_volume_limits=small_volume_limits,
        large_volume_limits=large_volume_limits,
        small_volume_increases_with_coordinate=True,
        # Reflection already opposes the two slider motions. Reversing this
        # volume convention again would make both cylinder volumes move nearly
        # together and would not represent the published opposed construction.
        large_volume_increases_with_coordinate=True,
        crank_angle_offset=crank_angle_offset,
        crank_direction=crank_direction,
        extrema_scan_count=extrema_scan_count,
    )


def published_e0_opposed_kinematics(
    small_volume_limits: CylinderVolumeLimits,
    large_volume_limits: CylinderVolumeLimits,
    *,
    ground_distance: float,
    connecting_rod_to_projected_stroke_ratio: float = 5.0,
    crank_angle_offset: float = 0.0,
    crank_direction: int = 1,
    extrema_scan_count: int = 721,
) -> SharedCrankFourBarVolumeKinematics:
    """Build the published E10, later optimizer-ranked E0, opposed pair."""

    return shared_crank_rocker_kinematics(
        SharedCrankRockerDesign(
            crank_ratio=0.73,
            small_coupler_ratio=1.0,
            small_rocker_ratio=1.103728,
            small_output_along_ratio=0.9753102035400876,
            small_output_normal_ratio=-0.5167064010195661,
            large_coupler_ratio=1.0,
            large_rocker_ratio=1.103728,
            large_output_along_ratio=0.9753102035400876,
            large_output_normal_ratio=0.5167064010195661,
            small_slider_rod_ratio=connecting_rod_to_projected_stroke_ratio,
            large_slider_rod_ratio=connecting_rod_to_projected_stroke_ratio,
        ),
        small_volume_limits,
        large_volume_limits,
        ground_distance=ground_distance,
        crank_angle_offset=crank_angle_offset,
        crank_direction=crank_direction,
        extrema_scan_count=extrema_scan_count,
    )


def slider_rod_sensitivity(
    assembly: FourBarSliderAssembly,
    crank_radius: float,
    *,
    sample_count: int = 1441,
) -> SliderRodSensitivity:
    """Measure finite-rod waveform distortion for one configured assembly."""

    _positive("Crank radius", crank_radius)
    if sample_count < 16:
        raise ValueError("Sensitivity sample count must be at least 16.")
    projected = np.empty(sample_count)
    slider = np.empty(sample_count)
    for index, angle in enumerate(
        np.linspace(0.0, 2.0 * math.pi, sample_count, endpoint=False)
    ):
        pin = crank_radius * np.array((math.cos(angle), math.sin(angle)))
        derivative = crank_radius * np.array((-math.sin(angle), math.cos(angle)))
        state = assembly.evaluate(pin, derivative)
        axis = _unit_from_angle(assembly.slider.axis_angle)
        projected[index] = float(np.dot(state.output_point, axis))
        slider[index] = state.coordinate
    projected_stroke = float(np.ptp(projected))
    slider_stroke = float(np.ptp(slider))
    projected_normalized = (projected - np.min(projected)) / projected_stroke
    slider_normalized = (slider - np.min(slider)) / slider_stroke
    return SliderRodSensitivity(
        connecting_rod_to_projected_stroke_ratio=(
            assembly.slider.connecting_rod_length / projected_stroke
        ),
        projected_stroke=projected_stroke,
        slider_stroke=slider_stroke,
        maximum_normalized_shape_error=float(
            np.max(np.abs(slider_normalized - projected_normalized))
        ),
    )


# Family-facing name; retain the original concrete type for animation/sizing.
FourBarKinematics = SharedCrankFourBarVolumeKinematics
