"""Independent-phase wrapper around production canonical periodic splines."""
from dataclasses import dataclass
import math
import numpy as np
from dada_solver.free_kinematics import FreeKinematics

@dataclass(frozen=True)
class PhaseShiftedFreeKinematics:
    """Production FreeKinematics with independent piston phase shifts.

    The generic FreeKinematics API remains authoritative.  The combined scalar
    evaluator below merely evaluates its already-built cubic coefficients
    directly, avoiding NumPy allocation on every LSODA RHS call.
    """

    base: FreeKinematics
    small_phase_rad: float = 0.0
    large_phase_rad: float = 0.0

    @property
    def small_volume_limits(self):
        return self.base.small_volume_limits

    @property
    def large_volume_limits(self):
        return self.base.large_volume_limits

    @property
    def small_physical_stroke(self):
        return None

    @property
    def large_physical_stroke(self):
        return None

    @property
    def diagnostics(self):
        return self.base.diagnostics

    def require_feasible(self):
        self.base.require_feasible()

    @staticmethod
    def _scalar_pair(motion, theta):
        """Return volume and dV/dtheta from production spline coefficients."""
        wrapped = float(theta) % (2.0 * math.pi)
        knots = motion._knots
        n = len(knots) - 1
        # Equally spaced knots by FreeKinematics construction.
        index = int(wrapped * n / (2.0 * math.pi))
        if index >= n:
            index = n - 1
        dx = wrapped - knots[index]
        coeff = motion._coefficients
        a = coeff[0][index]
        b = coeff[1][index]
        c = coeff[2][index]
        d = coeff[3][index]
        raw = ((a * dx + b) * dx + c) * dx + d
        draw = (3.0 * a * dx + 2.0 * b) * dx + c
        volume = (
            (raw - motion._minimum) * motion._scale
            + motion.definition.minimum_volume
        )
        derivative = draw * motion._scale
        return volume, derivative

    def small_cylinder_volume(self, theta):
        if np.ndim(theta):
            return self.base.small_cylinder_volume(
                np.asarray(theta) - self.small_phase_rad
            )
        return self._scalar_pair(
            self.base._small, float(theta) - self.small_phase_rad
        )[0]

    def large_cylinder_volume(self, theta):
        if np.ndim(theta):
            return self.base.large_cylinder_volume(
                np.asarray(theta) - self.large_phase_rad
            )
        return self._scalar_pair(
            self.base._large, float(theta) - self.large_phase_rad
        )[0]

    def small_cylinder_volume_derivative(self, theta):
        if np.ndim(theta):
            return self.base.small_cylinder_volume_derivative(
                np.asarray(theta) - self.small_phase_rad
            )
        return self._scalar_pair(
            self.base._small, float(theta) - self.small_phase_rad
        )[1]

    def large_cylinder_volume_derivative(self, theta):
        if np.ndim(theta):
            return self.base.large_cylinder_volume_derivative(
                np.asarray(theta) - self.large_phase_rad
            )
        return self._scalar_pair(
            self.base._large, float(theta) - self.large_phase_rad
        )[1]

    def small_cylinder_volume_second_derivative(self, theta):
        return self.base.small_cylinder_volume_second_derivative(
            np.asarray(theta) - self.small_phase_rad
        )

    def large_cylinder_volume_second_derivative(self, theta):
        return self.base.large_cylinder_volume_second_derivative(
            np.asarray(theta) - self.large_phase_rad
        )

    def cylinder_volumes_and_derivatives(self, theta):
        # LSODA / CompiledWallRHS uses scalar angles here.
        if np.ndim(theta):
            return (
                self.small_cylinder_volume(theta),
                self.large_cylinder_volume(theta),
                self.small_cylinder_volume_derivative(theta),
                self.large_cylinder_volume_derivative(theta),
            )
        sv, ds = self._scalar_pair(
            self.base._small, float(theta) - self.small_phase_rad
        )
        lv, dl = self._scalar_pair(
            self.base._large, float(theta) - self.large_phase_rad
        )
        return sv, lv, ds, dl

    def breakpoint_angles(self):
        return ()

