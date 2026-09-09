"""Interfaces for imposed, periodic cylinder-volume laws."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Protocol, runtime_checkable

from dada_solver.geometry import CylinderVolumeLimits


@runtime_checkable
class KinematicsModel(Protocol):
    """Imposed volume boundary, independent of forces or mechanism family.

    Angles are radians, one cycle is 2*pi. Only first derivatives are required
    by thermodynamics; smooth backends may expose analytic second derivatives.
    """

    @property
    def small_volume_limits(self) -> CylinderVolumeLimits: ...

    @property
    def large_volume_limits(self) -> CylinderVolumeLimits: ...

    @property
    def small_physical_stroke(self) -> float | None: ...

    @property
    def large_physical_stroke(self) -> float | None: ...

    def small_cylinder_volume(self, theta: float) -> float:
        """Return small-cylinder volume in m^3."""

    def small_cylinder_volume_derivative(self, theta: float) -> float:
        """Return dV_S/dtheta in m^3/rad."""

    def large_cylinder_volume(self, theta: float) -> float:
        """Return large-cylinder volume in m^3."""

    def large_cylinder_volume_derivative(self, theta: float) -> float:
        """Return dV_L/dtheta in m^3/rad."""

    def breakpoint_angles(self) -> tuple[float, ...]:
        """Return derivative-discontinuity angles strictly inside one cycle."""


class KinematicConstraintViolation(ValueError):
    """An explicitly infeasible motion with individually inspectable diagnostics."""
    def __init__(self, message, diagnostics):
        super().__init__(message)
        self.diagnostics = tuple(diagnostics)


class SmoothKinematicsModel(KinematicsModel, Protocol):
    """Optional analytic acceleration capability; not required by thermodynamics."""
    def small_cylinder_volume_second_derivative(self, theta: float) -> float: ...
    def large_cylinder_volume_second_derivative(self, theta: float) -> float: ...


# Backward-compatible public name.
VolumeKinematics = KinematicsModel


@dataclass(frozen=True, slots=True)
class ReversedVolumeKinematics:
    """Evaluate V(-phi), with increasing cycle progress phi = |omega| t.

    The initial geometric origin is preserved. Derivatives and discontinuity
    positions are transformed as well as volumes; time never runs backwards.
    """

    forward: VolumeKinematics

    def small_cylinder_volume(self, theta: float) -> float:
        return self.forward.small_cylinder_volume(-theta)

    def large_cylinder_volume(self, theta: float) -> float:
        return self.forward.large_cylinder_volume(-theta)

    def small_cylinder_volume_derivative(self, theta: float) -> float:
        return -self.forward.small_cylinder_volume_derivative(-theta)

    def large_cylinder_volume_derivative(self, theta: float) -> float:
        return -self.forward.large_cylinder_volume_derivative(-theta)

    @property
    def small_volume_limits(self):
        return self.forward.small_volume_limits

    @property
    def large_volume_limits(self):
        return self.forward.large_volume_limits

    @property
    def small_physical_stroke(self):
        return getattr(self.forward, 'small_physical_stroke', None)

    @property
    def large_physical_stroke(self):
        return getattr(self.forward, 'large_physical_stroke', None)

    @property
    def diagnostics(self):
        return self.forward.diagnostics

    def small_cylinder_volume_second_derivative(self, theta):
        return self.forward.small_cylinder_volume_second_derivative(-theta)

    def large_cylinder_volume_second_derivative(self, theta):
        return self.forward.large_cylinder_volume_second_derivative(-theta)

    def breakpoint_angles(self) -> tuple[float, ...]:
        provider = getattr(self.forward, "breakpoint_angles", None)
        if provider is None:
            return ()
        return tuple(sorted(2.0 * math.pi - angle for angle in provider()))

    def cylinder_volumes_and_derivatives(
        self, theta: float
    ) -> tuple[float, float, float, float]:
        provider = getattr(self.forward, "cylinder_volumes_and_derivatives", None)
        if provider is not None:
            small, large, ds, dl = provider(-theta)
            return small, large, -ds, -dl
        return (
            self.small_cylinder_volume(theta),
            self.large_cylinder_volume(theta),
            self.small_cylinder_volume_derivative(theta),
            self.large_cylinder_volume_derivative(theta),
        )


@dataclass(frozen=True, slots=True)
class HarmonicVolumeKinematics:
    """Generic controlled example; not a validated DADA linkage model."""

    small_limits: CylinderVolumeLimits
    large_limits: CylinderVolumeLimits
    small_phase_offset: float

    @property
    def small_volume_limits(self):
        return self.small_limits

    @property
    def large_volume_limits(self):
        return self.large_limits

    @property
    def small_physical_stroke(self):
        return None

    @property
    def large_physical_stroke(self):
        return None

    def small_cylinder_volume(self, theta: float) -> float:
        midpoint = 0.5 * (self.small_limits.minimum + self.small_limits.maximum)
        amplitude = 0.5 * self.small_limits.swept
        return midpoint + amplitude * math.cos(theta - self.small_phase_offset)

    def small_cylinder_volume_derivative(self, theta: float) -> float:
        amplitude = 0.5 * self.small_limits.swept
        return -amplitude * math.sin(theta - self.small_phase_offset)

    def large_cylinder_volume(self, theta: float) -> float:
        midpoint = 0.5 * (self.large_limits.minimum + self.large_limits.maximum)
        amplitude = 0.5 * self.large_limits.swept
        return midpoint + amplitude * math.cos(theta)

    def large_cylinder_volume_derivative(self, theta: float) -> float:
        amplitude = 0.5 * self.large_limits.swept
        return -amplitude * math.sin(theta)

    def breakpoint_angles(self) -> tuple[float, ...]:
        return ()

    def cylinder_volumes_and_derivatives(
        self, theta: float
    ) -> tuple[float, float, float, float]:
        """Return V_S, V_L, dV_S/dtheta and dV_L/dtheta together."""

        small_midpoint = 0.5 * (self.small_limits.minimum + self.small_limits.maximum)
        small_amplitude = 0.5 * self.small_limits.swept
        large_midpoint = 0.5 * (self.large_limits.minimum + self.large_limits.maximum)
        large_amplitude = 0.5 * self.large_limits.swept
        small_angle = theta - self.small_phase_offset
        return (
            small_midpoint + small_amplitude * math.cos(small_angle),
            large_midpoint + large_amplitude * math.cos(theta),
            -small_amplitude * math.sin(small_angle),
            -large_amplitude * math.sin(theta),
        )


@dataclass(frozen=True, slots=True)
class IdealPiecewiseLinearVolumeKinematics:
    """Ideal four-sector receiver-cycle kinematics, not a linkage model."""

    small_limits: CylinderVolumeLimits
    large_limits: CylinderVolumeLimits
    small_lambda_target: float
    large_lambda_target: float
    adiabatic_sector_fraction: float = 0.25

    @property
    def small_volume_limits(self):
        return self.small_limits

    @property
    def large_volume_limits(self):
        return self.large_limits

    @property
    def small_physical_stroke(self):
        return None

    @property
    def large_physical_stroke(self):
        return None

    def __post_init__(self) -> None:
        for name, value in (
            ("small Lambda target", self.small_lambda_target),
            ("large Lambda target", self.large_lambda_target),
        ):
            if not math.isfinite(value) or not 0.0 < value < 1.0:
                raise ValueError(f"{name} must lie strictly between zero and one.")
        if (
            not math.isfinite(self.adiabatic_sector_fraction)
            or not 0.0 < self.adiabatic_sector_fraction < 0.5
        ):
            raise ValueError("Adiabatic sector fraction must lie strictly between 0 and 0.5.")

    @property
    def transfer_sector_fraction(self) -> float:
        return 0.5 - self.adiabatic_sector_fraction

    def breakpoint_angles(self) -> tuple[float, ...]:
        fractions = self._breakpoint_fractions()
        return tuple(2.0 * math.pi * value for value in fractions[1:-1])

    def small_cylinder_volume(self, theta: float) -> float:
        value, _ = self._normalized_volume_and_derivative(theta, small=True)
        return self.small_limits.minimum + value * self.small_limits.swept

    def small_cylinder_volume_derivative(self, theta: float) -> float:
        _, derivative = self._normalized_volume_and_derivative(theta, small=True)
        return derivative * self.small_limits.swept

    def large_cylinder_volume(self, theta: float) -> float:
        value, _ = self._normalized_volume_and_derivative(theta, small=False)
        return self.large_limits.minimum + value * self.large_limits.swept

    def large_cylinder_volume_derivative(self, theta: float) -> float:
        _, derivative = self._normalized_volume_and_derivative(theta, small=False)
        return derivative * self.large_limits.swept

    def cylinder_volumes_and_derivatives(
        self, theta: float
    ) -> tuple[float, float, float, float]:
        """Return both imposed volume laws from one sector lookup."""

        index, local, width = self._sector(theta)
        small_values = (0.0, 0.0, 1.0 - self.small_lambda_target, 1.0, 0.0)
        large_values = (1.0, 1.0 - self.large_lambda_target, 0.0, 0.0, 1.0)

        def interpolate(values):
            difference = values[index + 1] - values[index]
            return values[index] + local * difference, difference / (
                width * 2.0 * math.pi
            )

        small, small_derivative = interpolate(small_values)
        large, large_derivative = interpolate(large_values)
        return (
            self.small_limits.minimum + small * self.small_limits.swept,
            self.large_limits.minimum + large * self.large_limits.swept,
            small_derivative * self.small_limits.swept,
            large_derivative * self.large_limits.swept,
        )

    def _breakpoint_fractions(self) -> tuple[float, ...]:
        adiabatic = self.adiabatic_sector_fraction
        transfer = self.transfer_sector_fraction
        return (0.0, adiabatic, adiabatic + transfer, 2.0 * adiabatic + transfer, 1.0)

    def _normalized_volume_and_derivative(
        self, theta: float, *, small: bool
    ) -> tuple[float, float]:
        index, local, width = self._sector(theta)
        small_target_volume = 1.0 - self.small_lambda_target
        large_target_volume = 1.0 - self.large_lambda_target
        # Study origin and receiver-cycle order:
        # I: L compression; II: L -> H -> S; III: S expansion; IV: S -> C -> L.
        values = (
            (0.0, 0.0, small_target_volume, 1.0, 0.0)
            if small
            else (1.0, large_target_volume, 0.0, 0.0, 1.0)
        )
        value = values[index] + local * (values[index + 1] - values[index])
        derivative = (values[index + 1] - values[index]) / (width * 2.0 * math.pi)
        return value, derivative

    def _sector(self, theta: float) -> tuple[int, float, float]:
        cycle_fraction = (theta % (2.0 * math.pi)) / (2.0 * math.pi)
        breaks = self._breakpoint_fractions()
        for index in range(4):
            if cycle_fraction < breaks[index + 1] or index == 3:
                width = breaks[index + 1] - breaks[index]
                return index, (cycle_fraction - breaks[index]) / width, width
        raise RuntimeError("Unreachable piecewise-kinematics sector.")
