#!/usr/bin/env python3
"""Nine-parameter hybrid compact motor-motion search at DeltaT = 260 K.

The motion family deliberately treats the two thermodynamic sides differently.

BP-like exchange branches:
    large piston: maximum -> minimum
    small piston: minimum -> maximum

These use a nearly linear ramp with C2 rounded ends.  One rounding fraction per
branch controls how much of the branch is spent accelerating/decelerating.

Opposite-side branches:
    small piston: maximum -> minimum
    large piston: minimum -> maximum

These use a monotone C2 two-segment quintic law passing through one mobile
interior point (u_k, q_k).  The point creates a clear but smooth change of slope
without allowing oscillations or extra extrema.

Nine optimized parameters for BOTH pistons:
    small_max_deg
    small_down_duration_deg
    large_down_duration_deg
    large_down_rounding
    small_up_rounding
    small_down_kink_u
    small_down_kink_q
    large_up_kink_u
    large_up_kink_q

A common cycle rotation is redundant, so the large-piston maximum is fixed at
motor angle 0 deg.

Storage remains compact:
  * history.jsonl: 9 parameters, compact metrics, 10-state warm start
  * report.json: compact champion and campaign summary
  * best_diagnostics.json: full solver result only for the current champion
  * best_motion.csv: phase-aligned source16 and hybrid motion
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
from types import SimpleNamespace
import argparse
import hashlib
import json
import math
import time

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import qmc

from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.integration import IntegrationInterrupted
from dada_solver.wall_backend import WallBackendSettings

from compare_motor_motion_laws_stage7A5 import ROOT, _evaluate
from optimize_motor_temperature_point import base_geometry, feasibility
from refine_motor_four_stage_hx9d_variable_gas import _load_basis
from refine_motor_fourier_c2_260k import build_design as fourier_build_design
from optimize_motor_free_spline_260k_v3 import (
    _fixed_inventory_design,
    _load_history,
    _make_kinematics,
    _normalized_motion,
    _safe_uniform_initial,
    _shape_diagnostics,
    _source_best,
    _source_harmonics,
)

POWER_FLOOR_W = 25.0
DEFAULT_SOURCE = ROOT / "outputs/motor_spline_from_linear_260k/report.json"
DEFAULT_FOURIER = ROOT / "outputs/motor_fourier_c2_8h_refine/report.json"
DEFAULT_OUTPUT = ROOT / "outputs/motor_hybrid_compact_260k"

MIN_BRANCH_DEG = 35.0
MAX_BRANCH_DEG = 325.0

ROUND_MIN = 0.008
ROUND_MAX = 0.48
KINK_MIN = 0.04
KINK_MAX = 0.96

ANGLE_RADII_DEG = (16.0, 11.0, 7.0, 4.5, 2.5, 1.3, 0.65, 0.30, 0.12)
ROUND_RADII = (0.080, 0.060, 0.045, 0.032, 0.022, 0.014, 0.009, 0.005, 0.0025)
KINK_RADII = (0.120, 0.090, 0.065, 0.045, 0.030, 0.020, 0.012, 0.007, 0.0035)

EVALS_PER_RADIUS = 56
WIDE_EVERY = 32
WIDE_ANGLE_DEG = 25.0
WIDE_ROUND = 0.12
WIDE_KINK = 0.18
DIRECT_WARM_DISTANCE = 0.55

VELOCITY_FIT_WEIGHT = 0.025


def resolve(path):
    return path if path.is_absolute() else ROOT / path


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


def _source_kinematics(best, limits):
    return _make_kinematics(
        limits,
        np.asarray(best["small_controls"], dtype=float),
        np.asarray(best["large_controls"], dtype=float),
        float(best["small_phase_deg"]),
        float(best["large_phase_deg"]),
    )


def _source_aligned_motion(source, large_max_deg, t):
    tref = np.mod(np.asarray(t, dtype=float) + large_max_deg/360.0, 1.0)
    return _normalized_motion(source, tref)


def _aligned_shape_distance(source, large_max_deg, candidate, n=4096):
    t = np.linspace(0.0, 1.0, n, endpoint=False)
    ref, dref = _source_aligned_motion(source, large_max_deg, t)
    cur, dcur = _normalized_motion(candidate, t)
    dv = cur-ref
    dd = dcur-dref
    return {
        "small_position_rms": float(np.sqrt(np.mean(dv[0]**2))),
        "large_position_rms": float(np.sqrt(np.mean(dv[1]**2))),
        "combined_position_rms": float(np.sqrt(np.mean(dv**2))),
        "small_velocity_rms": float(np.sqrt(np.mean(dd[0]**2))),
        "large_velocity_rms": float(np.sqrt(np.mean(dd[1]**2))),
        "combined_velocity_rms": float(np.sqrt(np.mean(dd**2))),
        "maximum_position_difference": float(np.max(np.abs(dv))),
    }


def _fit_seed(source, best, limits):
    d = best["shape_diagnostics"]
    lmax = float(d["large_max_deg"])

    small_max = (float(d["small_max_deg"])-lmax) % 360.0
    small_down = (
        float(d["small_min_deg"])-float(d["small_max_deg"])
    ) % 360.0
    large_down = (
        float(d["large_min_deg"])-lmax
    ) % 360.0

    # x: 3 angles, 2 rounding fractions, 4 kink coordinates.
    x0 = np.array([
        small_max,
        small_down,
        large_down,
        0.15,
        0.15,
        0.50,
        0.50,
        0.50,
        0.50,
    ], dtype=float)

    lo = np.array([
        0.0,
        MIN_BRANCH_DEG,
        MIN_BRANCH_DEG,
        ROUND_MIN,
        ROUND_MIN,
        KINK_MIN,
        KINK_MIN,
        KINK_MIN,
        KINK_MIN,
    ])
    hi = np.array([
        360.0,
        MAX_BRANCH_DEG,
        MAX_BRANCH_DEG,
        ROUND_MAX,
        ROUND_MAX,
        KINK_MAX,
        KINK_MAX,
        KINK_MAX,
        KINK_MAX,
    ])

    t = np.linspace(0.0, 1.0, 1440, endpoint=False)
    ref, dref = _source_aligned_motion(source, lmax, t)

    def unpack(x):
        return {
            "small_max_deg": float(x[0]) % 360.0,
            "small_down_duration_deg": float(x[1]),
            "large_down_duration_deg": float(x[2]),
            "large_down_rounding": float(x[3]),
            "small_up_rounding": float(x[4]),
            "small_down_kink_u": float(x[5]),
            "small_down_kink_q": float(x[6]),
            "large_up_kink_u": float(x[7]),
            "large_up_kink_q": float(x[8]),
        }

    def residual(x):
        kin = params_to_kinematics(limits, unpack(x))
        cur, dcur = _normalized_motion(kin, t)
        pos = (cur-ref).ravel()
        vel = VELOCITY_FIT_WEIGHT*(dcur-dref).ravel()
        return np.concatenate((pos, vel))

    fit = least_squares(
        residual,
        x0,
        bounds=(lo, hi),
        xtol=1e-11,
        ftol=1e-11,
        gtol=1e-11,
        max_nfev=3000,
    )
    p = unpack(fit.x)
    kin = params_to_kinematics(limits, p)
    return p, kin, {
        "success": bool(fit.success),
        "cost": float(fit.cost),
        "nfev": int(fit.nfev),
        "velocity_fit_weight": VELOCITY_FIT_WEIGHT,
        "shape_distance_vs_source16_aligned":
            _aligned_shape_distance(source, lmax, kin),
    }


def _reflect(x, lo, hi):
    w = hi-lo
    y = (x-lo) % (2.0*w)
    return lo + (2.0*w-y if y > w else y)


def _propose(center, u, angle_radius, round_radius, kink_radius):
    p = dict(center)
    u = np.asarray(u, dtype=float)

    p["small_max_deg"] = (
        float(center["small_max_deg"]) + angle_radius*(2*u[0]-1)
    ) % 360.0

    p["small_down_duration_deg"] = _reflect(
        float(center["small_down_duration_deg"])
        + angle_radius*(2*u[1]-1),
        MIN_BRANCH_DEG, MAX_BRANCH_DEG,
    )
    p["large_down_duration_deg"] = _reflect(
        float(center["large_down_duration_deg"])
        + angle_radius*(2*u[2]-1),
        MIN_BRANCH_DEG, MAX_BRANCH_DEG,
    )

    p["large_down_rounding"] = _reflect(
        float(center["large_down_rounding"])
        + round_radius*(2*u[3]-1),
        ROUND_MIN, ROUND_MAX,
    )
    p["small_up_rounding"] = _reflect(
        float(center["small_up_rounding"])
        + round_radius*(2*u[4]-1),
        ROUND_MIN, ROUND_MAX,
    )

    names = (
        "small_down_kink_u",
        "small_down_kink_q",
        "large_up_kink_u",
        "large_up_kink_q",
    )
    for j, name in enumerate(names):
        p[name] = _reflect(
            float(center[name]) + kink_radius*(2*u[5+j]-1),
            KINK_MIN, KINK_MAX,
        )

    return p


def _candidate_id(p):
    return hashlib.sha256(
        json.dumps(p, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _phase_delta(a, b):
    return ((float(a)-float(b)+180.0) % 360.0)-180.0


def _param_distance(p, r):
    q = r["parameters"]
    terms = [
        _phase_delta(p["small_max_deg"], q["small_max_deg"])/30.0,
        (p["small_down_duration_deg"]-q["small_down_duration_deg"])/30.0,
        (p["large_down_duration_deg"]-q["large_down_duration_deg"])/30.0,
        (p["large_down_rounding"]-q["large_down_rounding"])/0.10,
        (p["small_up_rounding"]-q["small_up_rounding"])/0.10,
    ]
    for name in (
        "small_down_kink_u",
        "small_down_kink_q",
        "large_up_kink_u",
        "large_up_kink_q",
    ):
        terms.append((p[name]-q[name])/0.15)
    return float(np.linalg.norm(terms))


def _compact_result(r):
    c = r.get("convergence") or {}
    v = r.get("validity") or {}
    return {
        "status": r.get("status"),
        "indicated_thermal_efficiency": r.get("indicated_thermal_efficiency"),
        "indicated_power_w": r.get("indicated_power_w"),
        "heat_input_w": r.get("heat_input_w"),
        "heat_out_w": r.get("heat_out_w"),
        "maximum_tube_reynolds": r.get("maximum_tube_reynolds"),
        "maximum_tube_mach_number": r.get("maximum_tube_mach_number"),
        "maximum_absolute_mass_flow_kg_s":
            r.get("maximum_absolute_mass_flow_kg_s"),
        "maximum_pressure_pa": r.get("maximum_pressure_pa"),
        "maximum_temperature_k": r.get("maximum_temperature_k"),
        "total_mass_kg": r.get("total_mass_kg"),
        "validity_verdict": v.get("verdict"),
        "failed_criteria": v.get("failed_criteria"),
        "cycles_completed": c.get("cycles_completed"),
        "last_normalized_periodic_error":
            c.get("last_normalized_periodic_error"),
        "normalized_periodic_error_ratio":
            c.get("normalized_periodic_error_ratio"),
        "message": r.get("message") if r.get("status") != "converged" else None,
    }


def _eligible(r):
    return bool(
        r.get("feasible")
        and r.get("result", {}).get("status") == "converged"
        and r["result"].get("indicated_thermal_efficiency") is not None
    )


def _eta(r):
    return float(r["result"]["indicated_thermal_efficiency"])


def _backend_compact(s):
    if not s:
        return None
    cache = s.get("exact_kinematics_cache") or {}
    return {
        "actual_backend": s.get("actual_backend"),
        "fallback_calls": s.get("fallback_calls"),
        "exact_cache_hit_fraction": cache.get("hit_fraction"),
    }


def _save_motion(path, source, lmax, candidate):
    t = np.linspace(0.0, 1.0, 2881)
    ref, dref = _source_aligned_motion(source, lmax, t)
    cur, dcur = _normalized_motion(candidate, t)
    np.savetxt(
        path,
        np.column_stack((
            t, 360.0*t,
            ref[0], ref[1], cur[0], cur[1],
            dref[0], dref[1], dcur[0], dcur[1],
        )),
        delimiter=",",
        header=(
            "motor_time_fraction,motor_angle_deg,"
            "source16_aligned_small,source16_aligned_large,"
            "hybrid_small,hybrid_large,"
            "source16_aligned_dsmall_dt,source16_aligned_dlarge_dt,"
            "hybrid_dsmall_dt,hybrid_dlarge_dt"
        ),
        comments="",
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-report", type=Path, default=DEFAULT_SOURCE)
    ap.add_argument("--fourier-report", type=Path, default=DEFAULT_FOURIER)
    ap.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--evaluations", type=int, default=512)
    ap.add_argument("--budget-seconds", type=float, default=21600)
    ap.add_argument("--candidate-seconds", type=float, default=150)
    ap.add_argument(
        "--evaluations-per-radius", type=int, default=EVALS_PER_RADIUS
    )
    ap.add_argument("--seed", type=int, default=26092393)
    args = ap.parse_args()

    sp = resolve(args.source_report)
    fp = resolve(args.fourier_report)
    out = resolve(args.output_directory)
    out.mkdir(parents=True, exist_ok=True)

    sr = json.loads(sp.read_text())
    sb = sr.get("best_feasible")
    if sb is None:
        raise RuntimeError("Source report has no best_feasible.")

    fr = json.loads(fp.read_text())
    fb = _source_best(fr)
    h = _source_harmonics(fb)

    cd, base = _load_basis()
    geometry = base_geometry(base)
    machine = fourier_build_design(
        base,
        geometry,
        fb["thermo"],
        np.asarray(fb["small_coefficients"]),
        np.asarray(fb["large_coefficients"]),
        h,
    )
    total_mass = float(sb["result"]["total_mass_kg"])
    machine = _fixed_inventory_design(machine, total_mass)
    limits = machine.configuration.machine_volumes

    source = _source_kinematics(sb, limits)
    source_lmax = float(sb["shape_diagnostics"]["large_max_deg"])
    fit_params, fit_kin, fit_info = _fit_seed(source, sb, limits)

    print(json.dumps({
        "hybrid_shape_fit": {
            "parameters": fit_params,
            "combined_position_rms_vs_source16":
                fit_info["shape_distance_vs_source16_aligned"][
                    "combined_position_rms"
                ],
            "combined_velocity_rms_vs_source16":
                fit_info["shape_distance_vs_source16_aligned"][
                    "combined_velocity_rms"
                ],
            "maximum_position_difference_vs_source16":
                fit_info["shape_distance_vs_source16_aligned"][
                    "maximum_position_difference"
                ],
        }
    }, separators=(",", ":")), flush=True)

    definition = SimpleNamespace(
        wall_numerical_settings=cd.wall_numerical_settings,
        wall_backend=WallBackendSettings("numba"),
        exact_kinematics_cache=True,
        shared_replay=True,
    )

    identity = {
        "study": "260 K 9-parameter hybrid compact motion v1",
        "degrees_of_freedom_total": 9,
        "global_phase_gauge": "large maximum fixed at motor angle 0 deg",
        "bp_branches": {
            "large_down": "rounded linear C2",
            "small_up": "rounded linear C2",
        },
        "opposite_branches": {
            "small_down": "two C2 quintics through mobile (u_k,q_k)",
            "large_up": "two C2 quintics through mobile (u_k,q_k)",
        },
        "source_report": str(sp),
        "source_sha256":
            hashlib.sha256(sp.read_bytes()).hexdigest(),
        "fourier_report": str(fp),
        "fourier_sha256":
            hashlib.sha256(fp.read_bytes()).hexdigest(),
        "source_candidate_id": sb["candidate_id"],
        "source_large_max_deg_for_alignment": source_lmax,
        "fixed_total_mass_kg": total_mass,
        "fixed_thermo": fb["thermo"],
        "parameter_bounds": {
            "branch_deg": [MIN_BRANCH_DEG, MAX_BRANCH_DEG],
            "rounding_fraction": [ROUND_MIN, ROUND_MAX],
            "kink_u_q": [KINK_MIN, KINK_MAX],
        },
        "angle_radii_deg": list(ANGLE_RADII_DEG),
        "round_radii": list(ROUND_RADII),
        "kink_radii": list(KINK_RADII),
        "evaluations_per_radius": args.evaluations_per_radius,
        "wide_probe_every": WIDE_EVERY,
        "wide_angle_deg": WIDE_ANGLE_DEG,
        "wide_round": WIDE_ROUND,
        "wide_kink": WIDE_KINK,
        "velocity_fit_weight": VELOCITY_FIT_WEIGHT,
        "sobol_seed": args.seed,
        "storage_policy":
            "compact history/report; full solver result only for champion",
    }

    dp = out/"definition.json"
    hp = out/"history.jsonl"
    rp = out/"report.json"
    bp = out/"best_diagnostics.json"

    if dp.exists() and json.loads(dp.read_text()) != identity:
        raise RuntimeError(
            "Campaign definition changed; use a new output directory."
        )
    if not dp.exists():
        dp.write_text(json.dumps(identity, indent=2)+"\n")

    history = _load_history(hp)
    valid = [r for r in history if _eligible(r)]
    best = max(valid, key=_eta) if valid else None
    finished = {
        r["candidate_id"]
        for r in history
        if r.get("result", {}).get("status") != "interrupted"
    }

    points = qmc.Sobol(
        d=9, scramble=True, seed=args.seed
    ).random_base2(16)

    start = time.monotonic()
    deadline = start+args.budget_seconds
    completed = 0
    pi = max(
        (r.get("proposal_index", -1) for r in history),
        default=-1,
    )+1

    while completed < args.evaluations and time.monotonic() < deadline:
        idx = pi
        pi += 1

        if not history and idx == 0:
            params = dict(fit_params)
            kind = "source16_shape_fit_seed"
            ri = -1
            angle_radius = 0.0
            round_radius = 0.0
            kink_radius = 0.0
        else:
            searched = sum(
                1
                for r in history
                if r.get("kind") != "source16_shape_fit_seed"
            )
            ri = min(
                searched//args.evaluations_per_radius,
                len(ANGLE_RADII_DEG)-1,
            )

            if searched > 0 and searched % WIDE_EVERY == WIDE_EVERY-1:
                center = fit_params
                angle_radius = WIDE_ANGLE_DEG
                round_radius = WIDE_ROUND
                kink_radius = WIDE_KINK
                kind = "fit_seed_wide_probe"
            else:
                center = (
                    best["parameters"]
                    if best is not None
                    else fit_params
                )
                angle_radius = ANGLE_RADII_DEG[ri]
                round_radius = ROUND_RADII[ri]
                kink_radius = KINK_RADII[ri]
                kind = "hybrid_refinement"

            params = _propose(
                center,
                points[idx % len(points)],
                angle_radius,
                round_radius,
                kink_radius,
            )

        cid = _candidate_id(params)
        if cid in finished:
            continue

        try:
            kin = params_to_kinematics(limits, params)
        except ValueError:
            continue

        design = replace(machine, kinematics=kin)
        before = time.monotonic()
        cand_deadline = min(
            deadline, before+args.candidate_seconds
        )

        reusable = [
            r for r in history
            if r.get("result", {}).get("status") == "converged"
            and r.get("last_complete_state") is not None
        ]

        warm_distance = None
        if reusable:
            warm = min(
                reusable,
                key=lambda r: _param_distance(params, r),
            )
            warm_distance = _param_distance(params, warm)
            warm_source = warm["candidate_id"]

            if warm_distance <= DIRECT_WARM_DISTANCE:
                initial = np.asarray(
                    warm["last_complete_state"], dtype=float
                )
                warm_mode = "nearest_periodic_state"
            else:
                initial = _safe_uniform_initial(
                    design,
                    wall_source=np.asarray(
                        warm["last_complete_state"], dtype=float
                    ),
                )
                warm_mode = "safe_uniform_for_distant_shape"
        else:
            initial = _safe_uniform_initial(
                design,
                wall_source=np.asarray(
                    sb["last_complete_state"], dtype=float
                ),
            )
            warm_source = sb["candidate_id"]
            warm_mode = "safe_uniform_with_source16_wall_energy"

        def progress(_):
            now = time.monotonic()
            if now >= deadline:
                raise IntegrationInterrupted(
                    "Campaign budget exhausted."
                )
            if now >= cand_deadline:
                raise IntegrationInterrupted(
                    "Candidate budget exhausted."
                )

        box = {}

        def observe(periodic):
            box["statistics"] = periodic.backend_statistics

        safe_retry = False
        try:
            try:
                full, state = _evaluate(
                    f"hybrid_compact_260k_{idx}",
                    design,
                    definition,
                    initial_state=initial,
                    progress_callback=progress,
                    periodic_observer=observe,
                )
            except MicrotubeDomainError:
                safe_retry = True
                safe = _safe_uniform_initial(
                    design, wall_source=initial
                )
                full, state = _evaluate(
                    f"hybrid_compact_260k_{idx}_safe",
                    design,
                    definition,
                    initial_state=safe,
                    progress_callback=progress,
                    periodic_observer=observe,
                )
                warm_mode = (
                    "safe_uniform_retry_after_domain_error"
                )
        except MicrotubeDomainError as exc:
            full = {
                "status": "invalid_exchanger",
                "message": str(exc),
            }
            state = None
        except IntegrationInterrupted as exc:
            full = {
                "status": "interrupted",
                "message": str(exc),
            }
            state = None
        except (ValueError, RuntimeError, ArithmeticError) as exc:
            full = {
                "status": "integration_failure",
                "message": str(exc),
            }
            state = None

        feasible, reasons = (
            feasibility(full, POWER_FLOOR_W)
            if full.get("status") == "converged"
            else (False, [])
        )

        compact = _compact_result(full)
        shape = _shape_diagnostics(kin)
        distance = _aligned_shape_distance(
            source, source_lmax, kin, 2048
        )

        rec = {
            "index": len(history),
            "proposal_index": idx,
            "candidate_id": cid,
            "kind": kind,
            "radius_index": ri,
            "angle_radius_deg": angle_radius,
            "round_radius": round_radius,
            "kink_radius": kink_radius,
            "parameters": {
                k: float(v) for k, v in params.items()
            },
            "shape_diagnostics": shape,
            "shape_distance_vs_source16_aligned": distance,
            "result": compact,
            "feasible": feasible,
            "physical_constraint_failures": reasons,
            "elapsed_seconds": time.monotonic()-before,
            "warm_start_source": warm_source,
            "warm_start_mode": warm_mode,
            "warm_start_distance": warm_distance,
            "backend": _backend_compact(
                box.get("statistics")
            ),
            "safe_retry_used": safe_retry,
            "last_complete_state":
                state.tolist() if state is not None else None,
        }

        with hp.open("a") as f:
            f.write(
                json.dumps(rec, separators=(",", ":"))+"\n"
            )
            f.flush()

        history.append(rec)
        completed += 1

        if compact.get("status") != "interrupted":
            finished.add(cid)

        moved = False
        if feasible and (
            best is None or _eta(rec) > _eta(best)
        ):
            best = rec
            moved = True
            _save_motion(
                out/"best_motion.csv",
                source,
                source_lmax,
                kin,
            )
            bp.write_text(
                json.dumps({
                    "candidate_id": cid,
                    "parameters": rec["parameters"],
                    "shape_diagnostics": shape,
                    "shape_distance_vs_source16_aligned":
                        distance,
                    "full_result": full,
                }, indent=2)+"\n"
            )

        report = {
            "best_feasible": best,
            "source16": {
                "candidate_id": sb["candidate_id"],
                "efficiency":
                    sb["result"][
                        "indicated_thermal_efficiency"
                    ],
                "power_W":
                    sb["result"]["indicated_power_w"],
                "shape_diagnostics":
                    sb["shape_diagnostics"],
                "phase_alignment_large_max_deg":
                    source_lmax,
            },
            "fitted_seed": {
                "parameters": fit_params,
                "fit": fit_info,
                "shape_diagnostics":
                    _shape_diagnostics(fit_kin),
            },
            "attempted_total": len(history),
            "completed_this_run": completed,
            "requested_duration_seconds":
                args.budget_seconds,
            "actual_duration_seconds":
                time.monotonic()-start,
            "definition_file": "definition.json",
            "best_diagnostics_file":
                "best_diagnostics.json",
        }
        rp.write_text(
            json.dumps(report, indent=2)+"\n"
        )

        print(json.dumps({
            "index": rec["index"],
            "proposal_index": idx,
            "kind": kind,
            "radius_index": ri,
            "status": compact.get("status"),
            "feasible": feasible,
            "efficiency":
                compact.get(
                    "indicated_thermal_efficiency"
                ),
            "power_W":
                compact.get("indicated_power_w"),
            "small_max_deg":
                params["small_max_deg"],
            "small_min_deg": (
                params["small_max_deg"]
                + params["small_down_duration_deg"]
            ) % 360.0,
            "large_max_deg": 0.0,
            "large_min_deg":
                params["large_down_duration_deg"],
            "large_down_rounding":
                params["large_down_rounding"],
            "small_up_rounding":
                params["small_up_rounding"],
            "small_down_kink": [
                params["small_down_kink_u"],
                params["small_down_kink_q"],
            ],
            "large_up_kink": [
                params["large_up_kink_u"],
                params["large_up_kink_q"],
            ],
            "rms_vs_source16":
                distance["combined_position_rms"],
            "center_moved": moved,
            "best_efficiency":
                _eta(best) if best else None,
            "actual_backend":
                (rec["backend"] or {}).get(
                    "actual_backend"
                ),
            "elapsed_seconds":
                rec["elapsed_seconds"],
        }, separators=(",", ":")), flush=True)

    print(f"Saved {rp}", flush=True)


if __name__ == "__main__":
    main()
