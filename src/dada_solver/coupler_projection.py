"""Fast projection-only coupler-point four-bar kinematics."""

from __future__ import annotations

from dataclasses import dataclass, field
import math

import numpy as np
from scipy.optimize import brentq

from dada_solver.geometry import CylinderVolumeLimits


def _positive(name: str, value: float) -> None:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive.")


@dataclass(frozen=True, slots=True)
class CouplerProjectionSideDesign:
    coupler_ratio: float
    rocker_ratio: float
    pivot_x_ratio: float
    pivot_y_ratio: float
    output_along_ratio: float
    output_normal_ratio: float
    axis_angle_degrees: float
    assembly_branch: int

    def __post_init__(self) -> None:
        _positive("coupler ratio", self.coupler_ratio)
        _positive("rocker ratio", self.rocker_ratio)
        for name, value in (
            ("pivot x ratio", self.pivot_x_ratio),
            ("pivot y ratio", self.pivot_y_ratio),
            ("output along ratio", self.output_along_ratio),
            ("output normal ratio", self.output_normal_ratio),
            ("axis angle", self.axis_angle_degrees),
        ):
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite.")
        if self.assembly_branch not in (-1, 1):
            raise ValueError("assembly_branch must be -1 or 1.")


@dataclass(frozen=True, slots=True)
class SharedCrankCouplerProjectionDesign:
    small: CouplerProjectionSideDesign
    large: CouplerProjectionSideDesign


@dataclass(frozen=True, slots=True)
class _ProjectionSide:
    design: CouplerProjectionSideDesign
    scan_count: int = 721
    _minimum: float = field(init=False, repr=False)
    _maximum: float = field(init=False, repr=False)

    def __post_init__(self) -> None:
        angles = np.linspace(0.0, 2.0 * math.pi, self.scan_count, endpoint=False)
        values = np.empty(self.scan_count)
        derivatives = np.empty(self.scan_count)
        for index, angle in enumerate(angles):
            values[index], derivatives[index] = self.evaluate(float(angle))

        roots: list[float] = []
        for index, angle in enumerate(angles):
            next_index = (index + 1) % self.scan_count
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
                        lambda candidate: self.evaluate(candidate)[1],
                        float(angle),
                        end,
                    )
                    % (2.0 * math.pi)
                )

        extrema = [self.evaluate(angle)[0] for angle in roots]
        extrema.extend((float(np.min(values)), float(np.max(values))))
        minimum = min(extrema)
        maximum = max(extrema)
        if maximum - minimum <= 1.0e-12:
            raise ValueError("Projected four-bar has zero usable stroke.")
        object.__setattr__(self, "_minimum", minimum)
        object.__setattr__(self, "_maximum", maximum)

    @property
    def stroke(self) -> float:
        return self._maximum - self._minimum

    def evaluate(self, theta: float) -> tuple[float, float]:
        d = self.design

        pin_x = math.cos(theta)
        pin_y = math.sin(theta)
        pin_dx = -math.sin(theta)
        pin_dy = math.cos(theta)

        delta_x = d.pivot_x_ratio - pin_x
        delta_y = d.pivot_y_ratio - pin_y
        distance = math.hypot(delta_x, delta_y)
        if distance <= 1.0e-12:
            raise ValueError("Four-bar circle centers coincide.")

        direction_x = delta_x / distance
        direction_y = delta_y / distance
        along = (
            d.coupler_ratio**2
            - d.rocker_ratio**2
            + distance**2
        ) / (2.0 * distance)
        height_squared = d.coupler_ratio**2 - along**2
        if height_squared <= 1.0e-20:
            raise ValueError("Four-bar cannot close or is at toggle.")
        height = math.sqrt(height_squared)

        branch = d.assembly_branch
        joint_x = pin_x + along * direction_x - branch * height * direction_y
        joint_y = pin_y + along * direction_y + branch * height * direction_x

        coupler_x = joint_x - pin_x
        coupler_y = joint_y - pin_y
        rocker_x = joint_x - d.pivot_x_ratio
        rocker_y = joint_y - d.pivot_y_ratio

        determinant = coupler_x * rocker_y - coupler_y * rocker_x
        if abs(determinant) <= 1.0e-12:
            raise ValueError("Four-bar velocity closure is singular.")

        velocity_rhs = coupler_x * pin_dx + coupler_y * pin_dy
        joint_dx = velocity_rhs * rocker_y / determinant
        joint_dy = -velocity_rhs * rocker_x / determinant

        unit_x = coupler_x / d.coupler_ratio
        unit_y = coupler_y / d.coupler_ratio
        unit_dx = (joint_dx - pin_dx) / d.coupler_ratio
        unit_dy = (joint_dy - pin_dy) / d.coupler_ratio

        output_x = (
            pin_x
            + d.output_along_ratio * unit_x
            - d.output_normal_ratio * unit_y
        )
        output_y = (
            pin_y
            + d.output_along_ratio * unit_y
            + d.output_normal_ratio * unit_x
        )
        output_dx = (
            pin_dx
            + d.output_along_ratio * unit_dx
            - d.output_normal_ratio * unit_dy
        )
        output_dy = (
            pin_dy
            + d.output_along_ratio * unit_dy
            + d.output_normal_ratio * unit_dx
        )

        axis = math.radians(d.axis_angle_degrees)
        axis_x = math.cos(axis)
        axis_y = math.sin(axis)
        return (
            output_x * axis_x + output_y * axis_y,
            output_dx * axis_x + output_dy * axis_y,
        )

    def normalized(self, theta: float) -> tuple[float, float]:
        coordinate, derivative = self.evaluate(theta)
        return (
            (coordinate - self._minimum) / self.stroke,
            derivative / self.stroke,
        )


@dataclass(frozen=True, slots=True)
class SharedCrankCouplerProjectionKinematics:
    design: SharedCrankCouplerProjectionDesign
    small_volume_limits: CylinderVolumeLimits
    large_volume_limits: CylinderVolumeLimits
    _small: _ProjectionSide = field(init=False, repr=False)
    _large: _ProjectionSide = field(init=False, repr=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "_small", _ProjectionSide(self.design.small))
        object.__setattr__(self, "_large", _ProjectionSide(self.design.large))

    @property
    def small_physical_stroke(self):
        return None

    @property
    def large_physical_stroke(self):
        return None

    def breakpoint_angles(self) -> tuple[float, ...]:
        return ()

    @staticmethod
    def _volume(side, limits, theta):
        fraction, derivative = side.normalized(theta)
        return (
            limits.minimum + limits.swept * fraction,
            limits.swept * derivative,
        )

    def cylinder_volumes_and_derivatives(self, theta: float):
        sv, sd = self._volume(self._small, self.small_volume_limits, theta)
        lv, ld = self._volume(self._large, self.large_volume_limits, theta)
        return sv, lv, sd, ld

    def small_cylinder_volume(self, theta: float) -> float:
        return self._volume(self._small, self.small_volume_limits, theta)[0]

    def large_cylinder_volume(self, theta: float) -> float:
        return self._volume(self._large, self.large_volume_limits, theta)[0]

    def small_cylinder_volume_derivative(self, theta: float) -> float:
        return self._volume(self._small, self.small_volume_limits, theta)[1]

    def large_cylinder_volume_derivative(self, theta: float) -> float:
        return self._volume(self._large, self.large_volume_limits, theta)[1]


def shared_crank_coupler_projection_kinematics(
    design: SharedCrankCouplerProjectionDesign,
    small_volume_limits: CylinderVolumeLimits,
    large_volume_limits: CylinderVolumeLimits,
) -> SharedCrankCouplerProjectionKinematics:
    return SharedCrankCouplerProjectionKinematics(
        design=design,
        small_volume_limits=small_volume_limits,
        large_volume_limits=large_volume_limits,
    )
