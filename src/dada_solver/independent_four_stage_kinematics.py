"""Independent four-stage linear volume laws for the two cylinders."""

from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass
import math

from dada_solver.geometry import CylinderVolumeLimits


TAU = 2.0 * math.pi


@dataclass(frozen=True, slots=True)
class IndependentFourStageVolumeKinematics:
    """Four linear stages per cylinder with independent timing.

    Parameters are expressed in forward motor-cycle time.

    The large cylinder defines the global origin:
        L = (1, b_l, 0, a_l, 1)
        at local times (0, t1_l, t2_l, t3_l, 1).

    The small cylinder owns an independent cyclic origin t0_s:
        S = (b_s, 1, a_s, 0, b_s)
        at local times (0, t1_s, t2_s, t3_s, 1),

    where its local time is (global_motor_time - t0_s) mod 1.
    """

    small_volume_limits: CylinderVolumeLimits
    large_volume_limits: CylinderVolumeLimits
    t0_s: float
    t1_l: float
    t2_l: float
    t3_l: float
    t1_s: float
    t2_s: float
    t3_s: float
    a_l: float
    b_l: float
    a_s: float
    b_s: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.t0_s) or not -0.5 <= self.t0_s < 0.5:
            raise ValueError("Require -0.5 <= t0_s < 0.5.")
        self._validate_knots("large", self.t1_l, self.t2_l, self.t3_l)
        self._validate_knots("small", self.t1_s, self.t2_s, self.t3_s)
        for name, value in (
            ("a_l", self.a_l), ("b_l", self.b_l),
            ("a_s", self.a_s), ("b_s", self.b_s),
        ):
            if not math.isfinite(value) or not 0.0 <= value <= 1.0:
                raise ValueError(
                    f"{name} must be a finite swept-volume fraction in [0, 1]."
                )

    @staticmethod
    def _validate_knots(name: str, t1: float, t2: float, t3: float) -> None:
        if not all(math.isfinite(value) for value in (t1, t2, t3)):
            raise ValueError(f"{name} stage times must be finite.")
        if not 0.0 < t1 < t2 < t3 < 1.0:
            short = name[0]
            raise ValueError(
                f"Require 0 < t1_{short} < t2_{short} < t3_{short} < 1."
            )

    @property
    def small_physical_stroke(self):
        return None

    @property
    def large_physical_stroke(self):
        return None

    @staticmethod
    def _interpolate(local_time, knots, values, limits):
        index = min(bisect_right(knots, local_time) - 1, 3)
        width = knots[index + 1] - knots[index]
        slope = (values[index + 1] - values[index]) / width
        fraction = values[index] + slope * (local_time - knots[index])
        volume = limits.minimum + fraction * limits.swept
        derivative = -slope * limits.swept / TAU
        return volume, derivative

    def cylinder_volumes_and_derivatives(self, theta):
        if not math.isfinite(theta):
            raise ValueError("Angle must be finite.")

        global_time = (-theta / TAU) % 1.0
        small_time = (global_time - self.t0_s) % 1.0

        small, dsmall = self._interpolate(
            small_time,
            (0.0, self.t1_s, self.t2_s, self.t3_s, 1.0),
            (self.b_s, 1.0, self.a_s, 0.0, self.b_s),
            self.small_volume_limits,
        )
        large, dlarge = self._interpolate(
            global_time,
            (0.0, self.t1_l, self.t2_l, self.t3_l, 1.0),
            (1.0, self.b_l, 0.0, self.a_l, 1.0),
            self.large_volume_limits,
        )
        return small, large, dsmall, dlarge

    def small_cylinder_volume(self, theta):
        return self.cylinder_volumes_and_derivatives(theta)[0]

    def large_cylinder_volume(self, theta):
        return self.cylinder_volumes_and_derivatives(theta)[1]

    def small_cylinder_volume_derivative(self, theta):
        return self.cylinder_volumes_and_derivatives(theta)[2]

    def large_cylinder_volume_derivative(self, theta):
        return self.cylinder_volumes_and_derivatives(theta)[3]

    @staticmethod
    def _study_angle_from_motor_time(motor_time):
        motor_time %= 1.0
        return (TAU * (1.0 - motor_time)) % TAU

    def breakpoint_angles(self):
        motor_times = {
            self.t1_l,
            self.t2_l,
            self.t3_l,
            self.t0_s % 1.0,
            (self.t0_s + self.t1_s) % 1.0,
            (self.t0_s + self.t2_s) % 1.0,
            (self.t0_s + self.t3_s) % 1.0,
        }
        angles = {
            self._study_angle_from_motor_time(value)
            for value in motor_times
        }
        return tuple(sorted(angle for angle in angles if 0.0 < angle < TAU))
