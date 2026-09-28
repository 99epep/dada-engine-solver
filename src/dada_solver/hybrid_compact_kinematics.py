"""Nine-coordinate compact hybrid C2 law extracted from the historical motor study.

Rounded linear BP branches and two-quintic HP kinks are distinct from the
structured C2 15p representation. The large maximum fixes the motor-angle gauge
at zero; study angle is minus motor angle. Physical stroke remains unspecified.
The historical scalar path, vector path and sampled monotonicity guard are kept.
"""
from dataclasses import dataclass, field
import math
import numpy as np
from dada_solver.geometry import CylinderVolumeLimits

MIN_BRANCH_DEG = 35.0
MAX_BRANCH_DEG = 325.0
ROUND_MIN = 0.008
ROUND_MAX = 0.48
KINK_MIN = 0.04
KINK_MAX = 0.96

def _smoothstep3(x):
    """C1 velocity ramp: s(0)=0, s(1)=1, s'(0)=s'(1)=0."""
    return x*x*(3.0 - 2.0*x)


def _rounded_linear_scalar(u, rounding):
    """Monotone 0->1 ramp, C2 in position, with constant-speed middle."""
    if u <= 0.0:
        return 0.0, 0.0
    if u >= 1.0:
        return 1.0, 0.0

    r = float(rounding)
    vmax = 1.0 / (1.0 - r)

    if u < r:
        x = u / r
        q = vmax * r * (x**3 - 0.5*x**4)
        dq = vmax * _smoothstep3(x)
        return q, dq

    if u <= 1.0 - r:
        return vmax * (u - 0.5*r), vmax

    w = 1.0 - u
    x = w / r
    q_from_end = vmax * r * (x**3 - 0.5*x**4)
    dq = vmax * _smoothstep3(x)
    return 1.0 - q_from_end, dq


def _rounded_linear(u, rounding):
    a = np.asarray(u, dtype=float)
    if a.ndim == 0:
        return _rounded_linear_scalar(float(a), rounding)

    q = np.empty_like(a)
    dq = np.empty_like(a)
    r = float(rounding)
    vmax = 1.0 / (1.0 - r)

    left = a < r
    middle = (a >= r) & (a <= 1.0-r)
    right = a > 1.0-r

    if np.any(left):
        x = a[left] / r
        q[left] = vmax*r*(x**3 - 0.5*x**4)
        dq[left] = vmax*(x*x*(3.0-2.0*x))

    if np.any(middle):
        q[middle] = vmax*(a[middle] - 0.5*r)
        dq[middle] = vmax

    if np.any(right):
        x = (1.0-a[right]) / r
        q[right] = 1.0 - vmax*r*(x**3 - 0.5*x**4)
        dq[right] = vmax*(x*x*(3.0-2.0*x))

    q[a <= 0.0] = 0.0
    dq[a <= 0.0] = 0.0
    q[a >= 1.0] = 1.0
    dq[a >= 1.0] = 0.0
    return q, dq


def _quintic_segment_coefficients(y0, y1, width, slope0, slope1):
    """Local x in [0,1], zero global second derivative at both ends."""
    h = float(width)
    c0 = float(y0)
    c1 = h*float(slope0)
    c2 = 0.0

    rhs0 = float(y1) - c0 - c1
    rhs1 = h*float(slope1) - c1
    rhs2 = 0.0

    # Solve:
    # c3+c4+c5 = rhs0
    # 3c3+4c4+5c5 = rhs1
    # 6c3+12c4+20c5 = rhs2
    c3 = 10.0*rhs0 - 4.0*rhs1 + 0.5*rhs2
    c4 = -15.0*rhs0 + 7.0*rhs1 - rhs2
    c5 = 6.0*rhs0 - 3.0*rhs1 + 0.5*rhs2
    return (c0, c1, c2, c3, c4, c5)


def _poly_pair(coeffs, x, width):
    c0, c1, c2, c3, c4, c5 = coeffs
    y = ((((c5*x + c4)*x + c3)*x + c2)*x + c1)*x + c0
    dy_dx = (((5.0*c5*x + 4.0*c4)*x + 3.0*c3)*x + 2.0*c2)*x + c1
    return y, dy_dx / width


def _kink_coefficients(kink_u, kink_q):
    u = float(kink_u)
    q = float(kink_q)
    d1 = q/u
    d2 = (1.0-q)/(1.0-u)
    # Harmonic mean keeps the shared knot slope between the two secants and
    # proved monotone over the admitted (u,q) box.
    slope = 2.0*d1*d2/(d1+d2)
    left = _quintic_segment_coefficients(0.0, q, u, 0.0, slope)
    right = _quintic_segment_coefficients(q, 1.0, 1.0-u, slope, 0.0)
    return left, right, slope


def _kink_transition_scalar(u, kink_u, left, right):
    if u <= 0.0:
        return 0.0, 0.0
    if u >= 1.0:
        return 1.0, 0.0
    ku = float(kink_u)
    if u <= ku:
        return _poly_pair(left, u/ku, ku)
    width = 1.0-ku
    return _poly_pair(right, (u-ku)/width, width)


def _kink_transition(u, kink_u, left, right):
    a = np.asarray(u, dtype=float)
    if a.ndim == 0:
        return _kink_transition_scalar(float(a), kink_u, left, right)

    q = np.empty_like(a)
    dq = np.empty_like(a)
    ku = float(kink_u)
    mask = a <= ku

    if np.any(mask):
        x = a[mask]/ku
        c0,c1,c2,c3,c4,c5 = left
        q[mask] = ((((c5*x+c4)*x+c3)*x+c2)*x+c1)*x+c0
        dq[mask] = (
            (((5*c5*x+4*c4)*x+3*c3)*x+2*c2)*x+c1
        )/ku

    if np.any(~mask):
        width = 1.0-ku
        x = (a[~mask]-ku)/width
        c0,c1,c2,c3,c4,c5 = right
        q[~mask] = ((((c5*x+c4)*x+c3)*x+c2)*x+c1)*x+c0
        dq[~mask] = (
            (((5*c5*x+4*c4)*x+3*c3)*x+2*c2)*x+c1
        )/width

    q[a <= 0.0] = 0.0
    dq[a <= 0.0] = 0.0
    q[a >= 1.0] = 1.0
    dq[a >= 1.0] = 0.0
    return q, dq


@dataclass(frozen=True, slots=True)
class HybridCompactKinematics:
    small_limits: CylinderVolumeLimits
    large_limits: CylinderVolumeLimits

    small_max_deg: float
    small_down_duration_deg: float
    large_down_duration_deg: float

    large_down_rounding: float
    small_up_rounding: float

    small_down_kink_u: float
    small_down_kink_q: float
    large_up_kink_u: float
    large_up_kink_q: float

    _small_down_left: tuple = field(init=False, repr=False)
    _small_down_right: tuple = field(init=False, repr=False)
    _large_up_left: tuple = field(init=False, repr=False)
    _large_up_right: tuple = field(init=False, repr=False)

    def __post_init__(self):
        for x in (self.small_down_duration_deg, self.large_down_duration_deg):
            if not MIN_BRANCH_DEG <= x <= MAX_BRANCH_DEG:
                raise ValueError("Branch duration outside admissible range.")

        for x in (self.large_down_rounding, self.small_up_rounding):
            if not ROUND_MIN <= x <= ROUND_MAX:
                raise ValueError("Rounding fraction outside admissible range.")

        for x in (
            self.small_down_kink_u, self.small_down_kink_q,
            self.large_up_kink_u, self.large_up_kink_q,
        ):
            if not KINK_MIN <= x <= KINK_MAX:
                raise ValueError("Kink coordinate outside admissible range.")

        sl, sr, _ = _kink_coefficients(
            self.small_down_kink_u, self.small_down_kink_q
        )
        ll, lr, _ = _kink_coefficients(
            self.large_up_kink_u, self.large_up_kink_q
        )
        object.__setattr__(self, "_small_down_left", sl)
        object.__setattr__(self, "_small_down_right", sr)
        object.__setattr__(self, "_large_up_left", ll)
        object.__setattr__(self, "_large_up_right", lr)

        # Cheap monotonicity guard for the derived C2 quintics.
        grid = np.linspace(0.0, 1.0, 129)
        _, sd = _kink_transition(
            grid, self.small_down_kink_u, sl, sr
        )
        _, lu = _kink_transition(
            grid, self.large_up_kink_u, ll, lr
        )
        if float(np.min(sd)) < -1e-10 or float(np.min(lu)) < -1e-10:
            raise ValueError("Derived kink transition is not monotone.")

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

    @property
    def diagnostics(self):
        return ()

    def require_feasible(self):
        return None

    def breakpoint_angles(self):
        # Position, velocity and acceleration are continuous at all internal joins.
        return ()

    def _small_normalized(self, phi_deg):
        phi = np.asarray(phi_deg, dtype=float) % 360.0
        max_deg = float(self.small_max_deg) % 360.0
        from_max = np.mod(phi-max_deg, 360.0)
        down = from_max <= self.small_down_duration_deg

        q = np.empty_like(phi)
        dq = np.empty_like(phi)

        if np.any(down):
            u = from_max[down]/self.small_down_duration_deg
            y, dy = _kink_transition(
                u, self.small_down_kink_u,
                self._small_down_left, self._small_down_right
            )
            q[down] = 1.0-y
            dq[down] = -dy*(180.0/math.pi)/self.small_down_duration_deg

        if np.any(~down):
            width = 360.0-self.small_down_duration_deg
            u = (from_max[~down]-self.small_down_duration_deg)/width
            y, dy = _rounded_linear(u, self.small_up_rounding)
            q[~down] = y
            dq[~down] = dy*(180.0/math.pi)/width

        if q.ndim == 0:
            return float(q), float(dq)
        return q, dq

    def _large_normalized(self, phi_deg):
        phi = np.asarray(phi_deg, dtype=float) % 360.0
        down = phi <= self.large_down_duration_deg

        q = np.empty_like(phi)
        dq = np.empty_like(phi)

        if np.any(down):
            u = phi[down]/self.large_down_duration_deg
            y, dy = _rounded_linear(u, self.large_down_rounding)
            q[down] = 1.0-y
            dq[down] = -dy*(180.0/math.pi)/self.large_down_duration_deg

        if np.any(~down):
            width = 360.0-self.large_down_duration_deg
            u = (phi[~down]-self.large_down_duration_deg)/width
            y, dy = _kink_transition(
                u, self.large_up_kink_u,
                self._large_up_left, self._large_up_right
            )
            q[~down] = y
            dq[~down] = dy*(180.0/math.pi)/width

        if q.ndim == 0:
            return float(q), float(dq)
        return q, dq

    def _small_normalized_scalar(self, phi):
        phi = float(phi) % 360.0
        from_max = (phi-float(self.small_max_deg)) % 360.0

        if from_max <= self.small_down_duration_deg:
            u = from_max/self.small_down_duration_deg
            y, dy = _kink_transition_scalar(
                u, self.small_down_kink_u,
                self._small_down_left, self._small_down_right
            )
            return (
                1.0-y,
                -dy*(180.0/math.pi)/self.small_down_duration_deg,
            )

        width = 360.0-self.small_down_duration_deg
        u = (from_max-self.small_down_duration_deg)/width
        y, dy = _rounded_linear_scalar(u, self.small_up_rounding)
        return y, dy*(180.0/math.pi)/width

    def _large_normalized_scalar(self, phi):
        phi = float(phi) % 360.0
        if phi <= self.large_down_duration_deg:
            u = phi/self.large_down_duration_deg
            y, dy = _rounded_linear_scalar(u, self.large_down_rounding)
            return (
                1.0-y,
                -dy*(180.0/math.pi)/self.large_down_duration_deg,
            )

        width = 360.0-self.large_down_duration_deg
        u = (phi-self.large_down_duration_deg)/width
        y, dy = _kink_transition_scalar(
            u, self.large_up_kink_u,
            self._large_up_left, self._large_up_right
        )
        return y, dy*(180.0/math.pi)/width

    def _pair(self, theta, small):
        if np.ndim(theta) == 0:
            phi = (-float(theta)*180.0/math.pi) % 360.0
            if small:
                q, dq = self._small_normalized_scalar(phi)
                lim = self.small_limits
            else:
                q, dq = self._large_normalized_scalar(phi)
                lim = self.large_limits
            return lim.minimum + q*lim.swept, -dq*lim.swept

        phi = np.mod(-np.degrees(np.asarray(theta, dtype=float)), 360.0)
        if small:
            q, dq = self._small_normalized(phi)
            lim = self.small_limits
        else:
            q, dq = self._large_normalized(phi)
            lim = self.large_limits
        return lim.minimum + q*lim.swept, -dq*lim.swept

    def small_cylinder_volume(self, theta):
        return self._pair(theta, True)[0]

    def large_cylinder_volume(self, theta):
        return self._pair(theta, False)[0]

    def small_cylinder_volume_derivative(self, theta):
        return self._pair(theta, True)[1]

    def large_cylinder_volume_derivative(self, theta):
        return self._pair(theta, False)[1]

    def cylinder_volumes_and_derivatives(self, theta):
        sv, ds = self._pair(theta, True)
        lv, dl = self._pair(theta, False)
        return sv, lv, ds, dl


def params_to_kinematics(limits, p):
    return HybridCompactKinematics(
        limits.small_cylinder,
        limits.large_cylinder,
        float(p["small_max_deg"]) % 360.0,
        float(p["small_down_duration_deg"]),
        float(p["large_down_duration_deg"]),
        float(p["large_down_rounding"]),
        float(p["small_up_rounding"]),
        float(p["small_down_kink_u"]),
        float(p["small_down_kink_q"]),
        float(p["large_up_kink_u"]),
        float(p["large_up_kink_q"]),
    )


