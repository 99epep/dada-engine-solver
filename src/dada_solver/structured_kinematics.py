"""Structured C2 15-coordinate law and its validated scalar polynomial path.

Extracted from the source16 fitting and thermodynamic 260 K lineage.
The historical motor-time definition maps to study angle t = -theta/(2*pi).
"""
import math
import numpy as np
from scipy.interpolate import BPoly, PPoly

def _harmonic(a, b):
    a = max(float(a), 1e-12)
    b = max(float(b), 1e-12)
    return 2.0*a*b/(a+b)


def _bpoly_from_nodes(nodes):
    """nodes = [(u, q, dq/du, d2q/du2), ...]."""
    x = [float(n[0]) for n in nodes]
    yi = [
        [float(n[1]), float(n[2]), float(n[3])]
        for n in nodes
    ]
    return BPoly.from_derivatives(x, yi)


def _bp_branch(q_mid, acc_start, acc_end):
    """Increasing 0->1 branch with one global bowing parameter."""
    q_mid = float(q_mid)

    sec_left = q_mid / 0.5
    sec_right = (1.0-q_mid) / 0.5
    slope_mid = _harmonic(sec_left, sec_right)

    return _bpoly_from_nodes([
        (0.0, 0.0, 0.0, float(acc_start)),
        (0.5, q_mid, slope_mid, 0.0),
        (1.0, 1.0, 0.0, float(acc_end)),
    ])


def _kink_branch(kink_u, kink_q, width_rel, acc_start, acc_end):
    """Increasing 0->1 branch with a tunably hard internal C2 turn."""
    ku = float(kink_u)
    kq = float(kink_q)
    wr = float(width_rel)

    # Full transition width scales with the available distance to the nearer
    # endpoint.  Hence wr can approach zero without crossing an endpoint.
    full_w = wr * min(ku, 1.0-ku)
    half = 0.5*full_w
    ul = ku-half
    ur = ku+half

    # "Outer" slopes implied by the two sides of the kink.
    m1 = kq/ku
    m2 = (1.0-kq)/(1.0-ku)
    mk = _harmonic(m1, m2)

    # Side values lie on the two asymptotic secants, while the exact kink point
    # is retained at (ku,kq).  As width -> 0, the C2 turn becomes arbitrarily
    # concentrated around the kink.
    ql = kq - m1*half
    qr = kq + m2*half

    return _bpoly_from_nodes([
        (0.0, 0.0, 0.0, float(acc_start)),
        (ul, ql, m1, 0.0),
        (ku, kq, mk, 0.0),
        (ur, qr, m2, 0.0),
        (1.0, 1.0, 0.0, float(acc_end)),
    ])


class StructuredMotion15:
    def __init__(self, p):
        self.p = dict(p)
        for name,value in self.p.items():
            if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
                raise ValueError(f'{name} must be finite.')
            if name.endswith('duration_deg') and not 0 < value < 360:
                raise ValueError(f'{name} must lie strictly between 0 and 360 degrees.')
            if name.endswith('curvature') and value <= 0:
                raise ValueError(f'{name} must be positive.')
            if name.endswith(('_u','_q')) and not 0 < value < 1:
                raise ValueError(f'{name} must lie strictly between zero and one.')
            if name.endswith('width_rel') and not 0 < value < 2:
                raise ValueError(f'{name} must lie strictly between zero and two.')

        sd = self.p["small_down_duration_deg"]/360.0
        su = 1.0-sd
        ld = self.p["large_down_duration_deg"]/360.0
        lu = 1.0-ld

        # Global normalized motor-time extrema curvatures:
        # q''(max) = -Amax ; q''(min) = +Amin.
        smax = self.p["small_max_curvature"]
        smin = self.p["small_min_curvature"]
        lmax = self.p["large_max_curvature"]
        lmin = self.p["large_min_curvature"]

        # Small down (non-BP), physical q: 1 -> 0.
        # Helper y=1-q increases 0->1.
        self.small_down = _kink_branch(
            self.p["small_down_kink_u"],
            self.p["small_down_kink_q"],
            self.p["small_down_kink_width_rel"],
            +smax*sd*sd,
            -smin*sd*sd,
        )

        # Small up (BP), physical q increases 0->1.
        self.small_up = _bp_branch(
            self.p["small_up_bp_mid_q"],
            +smin*su*su,
            -smax*su*su,
        )

        # Large down (BP), physical q: 1 -> 0, helper y=1-q.
        self.large_down = _bp_branch(
            self.p["large_down_bp_mid_q"],
            +lmax*ld*ld,
            -lmin*ld*ld,
        )

        # Large up (non-BP), physical q increases 0->1.
        self.large_up = _kink_branch(
            self.p["large_up_kink_u"],
            self.p["large_up_kink_q"],
            self.p["large_up_kink_width_rel"],
            +lmin*lu*lu,
            -lmax*lu*lu,
        )

    @staticmethod
    def _eval(poly, u):
        u = np.asarray(u, dtype=float)
        return (
            np.asarray(poly(u), dtype=float),
            np.asarray(poly.derivative(1)(u), dtype=float),
            np.asarray(poly.derivative(2)(u), dtype=float),
        )

    def normalized(self, t):
        t = np.asarray(t, dtype=float) % 1.0
        p = self.p

        # Small piston.
        smax_t = (p["small_max_deg"] % 360.0)/360.0
        sd = p["small_down_duration_deg"]/360.0
        from_smax = np.mod(t-smax_t, 1.0)
        sdown = from_smax <= sd

        sq = np.empty_like(t)
        sv = np.empty_like(t)
        sa = np.empty_like(t)

        if np.any(sdown):
            u = from_smax[sdown]/sd
            y, dy, ddy = self._eval(self.small_down, u)
            sq[sdown] = 1.0-y
            sv[sdown] = -dy/sd
            sa[sdown] = -ddy/(sd*sd)

        if np.any(~sdown):
            su = 1.0-sd
            u = (from_smax[~sdown]-sd)/su
            y, dy, ddy = self._eval(self.small_up, u)
            sq[~sdown] = y
            sv[~sdown] = dy/su
            sa[~sdown] = ddy/(su*su)

        # Large piston. Gauge: large maximum at t=0.
        ld = p["large_down_duration_deg"]/360.0
        ldown = t <= ld

        lq = np.empty_like(t)
        lv = np.empty_like(t)
        la = np.empty_like(t)

        if np.any(ldown):
            u = t[ldown]/ld
            y, dy, ddy = self._eval(self.large_down, u)
            lq[ldown] = 1.0-y
            lv[ldown] = -dy/ld
            la[ldown] = -ddy/(ld*ld)

        if np.any(~ldown):
            lu = 1.0-ld
            u = (t[~ldown]-ld)/lu
            y, dy, ddy = self._eval(self.large_up, u)
            lq[~ldown] = y
            lv[~ldown] = dy/lu
            la[~ldown] = ddy/(lu*lu)

        return (
            np.vstack((sq, lq)),
            np.vstack((sv, lv)),
            np.vstack((sa, la)),
        )

    def monotonicity_penalty(self, n=160):
        u = np.linspace(0.0, 1.0, n)
        derivatives = [
            self.small_down.derivative(1)(u),
            self.small_up.derivative(1)(u),
            self.large_down.derivative(1)(u),
            self.large_up.derivative(1)(u),
        ]
        return np.concatenate([
            np.minimum(np.asarray(v, dtype=float), 0.0)
            for v in derivatives
        ])


class _FastPPoly:
    """Allocation-free scalar PPoly evaluator for the ODE hot path."""

    __slots__ = ("x", "c", "dc")

    def __init__(self, bpoly):
        pp = PPoly.from_bernstein_basis(bpoly)
        dp = pp.derivative(1)
        self.x = np.asarray(pp.x, dtype=float)
        self.c = np.asarray(pp.c, dtype=float)
        self.dc = np.asarray(dp.c, dtype=float)

    @staticmethod
    def _horner(coeffs, z):
        y = 0.0
        for a in coeffs:
            y = y*z + float(a)
        return y

    def pair(self, u):
        u = float(u)
        if u <= self.x[0]:
            i = 0
            z = 0.0
        elif u >= self.x[-1]:
            i = len(self.x)-2
            z = self.x[-1]-self.x[-2]
        else:
            i = int(np.searchsorted(self.x, u, side="right")-1)
            i = max(0, min(i, len(self.x)-2))
            z = u-self.x[i]

        return (
            self._horner(self.c[:, i], z),
            self._horner(self.dc[:, i], z),
        )


class StructuredKinematics15:
    """Solver-compatible adapter around StructuredMotion15.

    Vector paths use the validated BPoly representation.
    Scalar RHS paths use preconverted power-basis polynomials and manual Horner
    evaluation to avoid repeated scipy interpolation overhead.
    """

    __slots__ = (
        "small_limits", "large_limits", "params", "motion",
        "_small_max_t", "_sd", "_ld",
        "_small_down", "_small_up", "_large_down", "_large_up",
    )

    def __init__(self, limits, params):
        self.small_limits = limits.small_cylinder
        self.large_limits = limits.large_cylinder
        self.params = dict(params)
        self.motion = StructuredMotion15(self.params)

        self._small_max_t = (self.params["small_max_deg"] % 360.0)/360.0
        self._sd = self.params["small_down_duration_deg"]/360.0
        self._ld = self.params["large_down_duration_deg"]/360.0

        self._small_down = _FastPPoly(self.motion.small_down)
        self._small_up = _FastPPoly(self.motion.small_up)
        self._large_down = _FastPPoly(self.motion.large_down)
        self._large_up = _FastPPoly(self.motion.large_up)

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
        # All joins are C2. No discontinuous solver breakpoint is needed.
        return ()

    def _small_scalar_q_dqdt(self, t):
        from_max = (float(t)-self._small_max_t) % 1.0

        if from_max <= self._sd:
            u = from_max/self._sd
            y, dy = self._small_down.pair(u)
            return 1.0-y, -dy/self._sd

        su = 1.0-self._sd
        u = (from_max-self._sd)/su
        y, dy = self._small_up.pair(u)
        return y, dy/su

    def _large_scalar_q_dqdt(self, t):
        t = float(t) % 1.0

        if t <= self._ld:
            u = t/self._ld
            y, dy = self._large_down.pair(u)
            return 1.0-y, -dy/self._ld

        lu = 1.0-self._ld
        u = (t-self._ld)/lu
        y, dy = self._large_up.pair(u)
        return y, dy/lu

    def _pair(self, theta, small):
        if np.ndim(theta) == 0:
            t = (-float(theta)/(2.0*math.pi)) % 1.0
            if small:
                q, dqdt = self._small_scalar_q_dqdt(t)
                lim = self.small_limits
            else:
                q, dqdt = self._large_scalar_q_dqdt(t)
                lim = self.large_limits

            v = lim.minimum + q*lim.swept
            dvdtheta = -(dqdt/(2.0*math.pi))*lim.swept
            return v, dvdtheta

        a = np.asarray(theta, dtype=float)
        t = np.mod(-a/(2.0*math.pi), 1.0)
        q, dqdt, _ = self.motion.normalized(t)
        row = 0 if small else 1
        lim = self.small_limits if small else self.large_limits
        return (
            lim.minimum + q[row]*lim.swept,
            -(dqdt[row]/(2.0*math.pi))*lim.swept,
        )

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

    def _second(self, theta, side, limits):
        t=np.mod(-np.asarray(theta)/(2*math.pi),1.)
        value=self.motion.normalized(np.atleast_1d(t))[2][side]*limits.swept/(2*math.pi)**2
        return float(value[0]) if t.ndim==0 else value

    def small_cylinder_volume_second_derivative(self, theta):
        return self._second(theta,0,self.small_limits)

    def large_cylinder_volume_second_derivative(self, theta):
        return self._second(theta,1,self.large_limits)
