"""Normalized Fourier/C2 law, independent of search policy."""
from dataclasses import dataclass, field
import math
import numpy as np
from scipy.optimize import minimize_scalar
from dada_solver.geometry import CylinderVolumeLimits
HARMONICS = 4

@dataclass(frozen=True)
class FourierVolumeKinematics:
    """Periodic full-stroke Fourier kinematics, smooth to all orders."""

    small_volume_limits: CylinderVolumeLimits
    large_volume_limits: CylinderVolumeLimits
    small_coefficients: tuple[float, ...]
    large_coefficients: tuple[float, ...]
    harmonics: int = HARMONICS

    _small_offset: float = field(init=False, repr=False, compare=False)
    _small_scale: float = field(init=False, repr=False, compare=False)
    _large_offset: float = field(init=False, repr=False, compare=False)
    _large_scale: float = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        if type(self.harmonics) is not int or self.harmonics < 1:
            raise ValueError('Harmonics must be a positive integer.')
        expected = 2 * self.harmonics
        if len(self.small_coefficients) != expected or len(self.large_coefficients) != expected:
            raise ValueError(f"Expected {expected} coefficients per piston.")
        for coeffs in (self.small_coefficients, self.large_coefficients):
            if not np.all(np.isfinite(np.asarray(coeffs, float))):
                raise ValueError("Fourier coefficients must be finite.")
            if np.linalg.norm(np.asarray(coeffs, float)) < 1e-12:
                raise ValueError("Fourier coefficient vector must be non-zero.")

        so, ss = self._normalization(self.small_coefficients)
        lo, ls = self._normalization(self.large_coefficients)
        object.__setattr__(self, "_small_offset", so)
        object.__setattr__(self, "_small_scale", ss)
        object.__setattr__(self, "_large_offset", lo)
        object.__setattr__(self, "_large_scale", ls)

    @property
    def small_physical_stroke(self):
        return None

    @property
    def large_physical_stroke(self):
        return None

    def breakpoint_angles(self):
        # No derivative discontinuity exists.
        return ()

    def _raw(self, t, coeffs):
        t = np.asarray(t, dtype=float)
        out = np.zeros_like(t)
        for k in range(1, self.harmonics + 1):
            a = coeffs[2 * (k - 1)]
            b = coeffs[2 * (k - 1) + 1]
            w = 2.0 * math.pi * k
            out = out + a * np.cos(w * t) + b * np.sin(w * t)
        return out

    def _raw_dt(self, t, coeffs):
        t = np.asarray(t, dtype=float)
        out = np.zeros_like(t)
        for k in range(1, self.harmonics + 1):
            a = coeffs[2 * (k - 1)]
            b = coeffs[2 * (k - 1) + 1]
            w = 2.0 * math.pi * k
            out = out + (-a * w * np.sin(w * t) + b * w * np.cos(w * t))
        return out

    def _normalization(self, coeffs):
        coeffs = tuple(float(x) for x in coeffs)
        grid = np.linspace(0.0, 1.0, 2049, endpoint=False)
        raw = self._raw(grid, coeffs)

        i_min = int(np.argmin(raw))
        i_max = int(np.argmax(raw))
        dx = 1.0 / len(grid)

        def wrapped(x):
            return float(self._raw(np.asarray([x % 1.0]), coeffs)[0])

        def local_extreme(index, sign):
            center = grid[index]
            result = minimize_scalar(
                lambda x: sign * wrapped(x),
                bounds=(center - 2.0 * dx, center + 2.0 * dx),
                method="bounded",
                options={"xatol": 1e-13},
            )
            return wrapped(result.x)

        minimum = local_extreme(i_min, +1.0)
        maximum = local_extreme(i_max, -1.0)
        span = maximum - minimum
        if not math.isfinite(span) or span <= 1e-9:
            raise ValueError("Degenerate Fourier motion.")
        return minimum, 1.0 / span

    def normalized_fractions(self, t):
        sr = self._raw(t, self.small_coefficients)
        lr = self._raw(t, self.large_coefficients)
        s = (sr - self._small_offset) * self._small_scale
        l = (lr - self._large_offset) * self._large_scale
        return s, l

    def normalized_fraction_derivatives(self, t):
        sd = self._raw_dt(t, self.small_coefficients) * self._small_scale
        ld = self._raw_dt(t, self.large_coefficients) * self._large_scale
        return sd, ld

    def cylinder_volumes_and_derivatives(self, theta):
        if not math.isfinite(theta):
            raise ValueError("Angle must be finite.")
        t = (-theta / (2.0 * math.pi)) % 1.0
        s, l = self.normalized_fractions(t)
        sd, ld = self.normalized_fraction_derivatives(t)

        s = float(s)
        l = float(l)
        sd = float(sd)
        ld = float(ld)

        sv = self.small_volume_limits.minimum + s * self.small_volume_limits.swept
        lv = self.large_volume_limits.minimum + l * self.large_volume_limits.swept
        # dt/dtheta = -1/(2*pi)
        svd = -sd * self.small_volume_limits.swept / (2.0 * math.pi)
        lvd = -ld * self.large_volume_limits.swept / (2.0 * math.pi)
        return sv, lv, svd, lvd

    def small_cylinder_volume(self, theta):
        return self.cylinder_volumes_and_derivatives(theta)[0]

    def large_cylinder_volume(self, theta):
        return self.cylinder_volumes_and_derivatives(theta)[1]

    def small_cylinder_volume_derivative(self, theta):
        return self.cylinder_volumes_and_derivatives(theta)[2]

    def large_cylinder_volume_derivative(self, theta):
        return self.cylinder_volumes_and_derivatives(theta)[3]

    def _second(self, theta, coefficients, scale, limits):
        t=(-theta/(2*math.pi)) % 1.0
        value=sum(-k*k*(coefficients[2*(k-1)]*math.cos(2*math.pi*k*t)
                       +coefficients[2*(k-1)+1]*math.sin(2*math.pi*k*t))
                  for k in range(1,self.harmonics+1))
        return value*scale*limits.swept

    def small_cylinder_volume_second_derivative(self, theta):
        return self._second(theta,self.small_coefficients,self._small_scale,self.small_volume_limits)

    def large_cylinder_volume_second_derivative(self, theta):
        return self._second(theta,self.large_coefficients,self._large_scale,self.large_volume_limits)
