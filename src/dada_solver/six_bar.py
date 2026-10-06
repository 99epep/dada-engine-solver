"""Independent six-bar cylinder motions in the common study-angle convention.

All lengths are in crank-radius units. No physical crank scale, spatial mirror,
extra cylinder phase, or motor-direction reversal is inferred here.
"""
from __future__ import annotations

from dataclasses import dataclass, fields, field
import math

import numpy as np
from scipy.optimize import brentq

from dada_solver.geometry import CylinderVolumeLimits


def _rr(origin, velocity, pivot, moving_length, fixed_length, branch):
    """Circle closure and analytic velocity for one driven RR dyad."""
    x, y = origin
    dx, dy = velocity
    px, py = pivot
    sx, sy = px-x, py-y
    distance = math.hypot(sx, sy)
    if distance <= 1e-12:
        raise ValueError('Six-bar dyad centers coincide.')
    ux, uy = sx/distance, sy/distance
    along = (moving_length**2-fixed_length**2+distance**2)/(2*distance)
    height2 = moving_length**2-along**2
    if height2 <= 1e-10:
        raise ValueError('Six-bar dyad cannot close or is at toggle.')
    height = math.sqrt(height2)
    jx = x+along*ux-branch*height*uy
    jy = y+along*uy+branch*height*ux
    ax, ay = jx-x, jy-y
    bx, by = jx-px, jy-py
    determinant = ax*by-ay*bx
    if abs(determinant) <= 1e-10:
        raise ValueError('Six-bar dyad velocity closure is singular.')
    rhs = ax*dx+ay*dy
    return (jx, jy), (rhs*by/determinant, -rhs*bx/determinant)


def _rigid_point(origin, velocity, joint, joint_velocity, along, normal):
    """Point coordinates expressed as fractions of the origin-joint vector."""
    x, y = origin
    dx, dy = velocity
    ux, uy = joint[0]-x, joint[1]-y
    vx, vy = joint_velocity[0]-dx, joint_velocity[1]-dy
    return ((x+along*ux-normal*uy, y+along*uy+normal*ux),
            (dx+along*vx-normal*vy, dy+along*vy+normal*vx))


@dataclass(frozen=True, slots=True)
class SixBarCylinderMechanism:
    """A-B-C-D, E on BC, E-F-G, H on EF, and the positive H-P slider branch.

    Primary E coordinates are lengths; secondary H coordinates are fractions
    of EF, matching Stage 2F/L1. Extrema are found once from analytic velocity
    roots bracketed on a 1440-interval cycle scan. Geometry is checked eagerly
    on that scan and on every subsequent evaluation; no motion is clipped.
    """
    primary_ground: float
    primary_coupler: float
    primary_rocker: float
    primary_e_along: float
    primary_e_normal: float
    primary_phase: float
    second_pivot_x: float
    second_pivot_y: float
    link_ef: float
    link_gf: float
    h_along_over_ef: float
    h_normal_over_ef: float
    piston_rod: float
    slider_axis_offset: float
    slider_axis_angle: float
    primary_branch: int
    second_branch: int
    slider_min: float = field(init=False)
    slider_max: float = field(init=False)
    minimum_angle: float = field(init=False)
    maximum_angle: float = field(init=False)

    def __post_init__(self):
        for item in fields(self):
            if item.init and not math.isfinite(getattr(self, item.name)):
                raise ValueError(f'{item.name} must be finite.')
        for name in ('primary_ground', 'primary_coupler', 'primary_rocker',
                     'link_ef', 'link_gf', 'piston_rod'):
            if getattr(self, name) <= 0:
                raise ValueError(f'{name} must be positive.')
        if self.primary_branch not in (-1, 1) or self.second_branch not in (-1, 1):
            raise ValueError('Six-bar assembly branches must be -1 or +1.')
        g, c, r = self.primary_ground, self.primary_coupler, self.primary_rocker
        if abs(g-1) <= abs(c-r)+1e-10 or c+r <= g+1+1e-10:
            raise ValueError('Primary loop cannot close over a full revolution.')
        angles = np.linspace(0., 2*math.pi, 1441)
        derivatives = [self.slider_position_and_derivative(float(a))[1] for a in angles]
        roots = []
        for a, b, da, db in zip(angles[:-1], angles[1:], derivatives[:-1], derivatives[1:]):
            if da == 0:
                roots.append(float(a))
            if da*db < 0:
                roots.append(brentq(lambda t: self.slider_position_and_derivative(t)[1],
                                    float(a), float(b), xtol=1e-14))
        if len(roots) < 2:
            raise ValueError('Could not bracket six-bar slider travel extrema.')
        minimum = min(roots, key=lambda t: self.slider_position_and_derivative(t)[0])
        maximum = max(roots, key=lambda t: self.slider_position_and_derivative(t)[0])
        lo = self.slider_position_and_derivative(minimum)[0]
        hi = self.slider_position_and_derivative(maximum)[0]
        if hi-lo <= 1e-10:
            raise ValueError('Six-bar slider has zero usable stroke.')
        for name, value in (('slider_min', lo), ('slider_max', hi),
                            ('minimum_angle', minimum), ('maximum_angle', maximum)):
            object.__setattr__(self, name, value)

    @property
    def stroke_over_crank(self) -> float:
        return self.slider_max-self.slider_min

    def slider_position_and_derivative(self, theta: float) -> tuple[float, float]:
        return self._evaluate(theta)

    def joint_state(self, theta: float):
        """Return joints and mechanical indicators from the same closure as the RHS.

        Coordinates are in crank-radius units in the stored local frame.
        These instantaneous indicators do not certify full-cycle feasibility.
        """
        return self._evaluate(theta, frames=True)

    def _evaluate(self, theta, *, frames=False):
        if not math.isfinite(theta):
            raise ValueError('Study angle must be finite.')
        angle = theta % (2*math.pi) + self.primary_phase
        b = (math.cos(angle), math.sin(angle))
        bd = (-b[1], b[0])
        c, cd = _rr(b, bd, (self.primary_ground, 0.), self.primary_coupler,
                    self.primary_rocker, self.primary_branch)
        e, ed = _rigid_point(b, bd, c, cd, self.primary_e_along/self.primary_coupler,
                            self.primary_e_normal/self.primary_coupler)
        f, fd = _rr(e, ed, (self.second_pivot_x, self.second_pivot_y), self.link_ef,
                    self.link_gf, self.second_branch)
        h, hd = _rigid_point(e, ed, f, fd, self.h_along_over_ef, self.h_normal_over_ef)
        ax, ay = math.cos(self.slider_axis_angle), math.sin(self.slider_axis_angle)
        longitudinal = h[0]*ax+h[1]*ay
        transverse = -h[0]*ay+h[1]*ax-self.slider_axis_offset
        longitudinal_d = hd[0]*ax+hd[1]*ay
        transverse_d = -hd[0]*ay+hd[1]*ax
        margin2 = self.piston_rod**2-transverse**2
        if margin2 <= 1e-10:
            raise ValueError('Six-bar piston rod cannot close or is at toggle.')
        margin = math.sqrt(margin2)
        position = longitudinal+margin
        derivative = longitudinal_d-transverse*transverse_d/margin
        if not frames:
            return position, derivative
        def sine(u, v, length_product):
            return abs(u[0]*v[1]-u[1]*v[0])/length_product
        p = (position*ax-self.slider_axis_offset*ay,
             position*ay+self.slider_axis_offset*ax)
        return dict(joints=dict(A=(0.,0.), B=b, C=c, D=(self.primary_ground,0.),
                               E=e, F=f, G=(self.second_pivot_x,self.second_pivot_y), H=h, P=p),
                    position=position, derivative=derivative, transverse=transverse,
                    primary_transmission_sine=sine((c[0]-b[0],c[1]-b[1]),
                        (c[0]-self.primary_ground,c[1]),self.primary_coupler*self.primary_rocker),
                    secondary_transmission_sine=sine((f[0]-e[0],f[1]-e[1]),
                        (f[0]-self.second_pivot_x,f[1]-self.second_pivot_y),self.link_ef*self.link_gf),
                    rod_axis_cosine=margin/self.piston_rod)

    def normalized_motion(self, theta: float) -> tuple[float, float]:
        slider, derivative = self.slider_position_and_derivative(theta)
        return (1-(slider-self.slider_min)/self.stroke_over_crank,
                -derivative/self.stroke_over_crank)


@dataclass(frozen=True, slots=True)
class IndependentSixBarVolumeKinematics:
    """Independent S/L mechanisms, with the phases already in study angle."""
    small: SixBarCylinderMechanism
    large: SixBarCylinderMechanism
    small_volume_limits: CylinderVolumeLimits
    large_volume_limits: CylinderVolumeLimits

    @property
    def small_physical_stroke(self) -> None:
        return None

    @property
    def large_physical_stroke(self) -> None:
        return None

    def breakpoint_angles(self) -> tuple[float, ...]:
        return ()

    @staticmethod
    def _volume(mechanism, limits, theta):
        q, dq = mechanism.normalized_motion(theta)
        return limits.minimum+q*limits.swept, dq*limits.swept

    def small_cylinder_volume(self, theta: float) -> float:
        return self._volume(self.small, self.small_volume_limits, theta)[0]

    def small_cylinder_volume_derivative(self, theta: float) -> float:
        return self._volume(self.small, self.small_volume_limits, theta)[1]

    def large_cylinder_volume(self, theta: float) -> float:
        return self._volume(self.large, self.large_volume_limits, theta)[0]

    def large_cylinder_volume_derivative(self, theta: float) -> float:
        return self._volume(self.large, self.large_volume_limits, theta)[1]

    def cylinder_volumes_and_derivatives(self, theta: float):
        s, ds = self._volume(self.small, self.small_volume_limits, theta)
        l, dl = self._volume(self.large, self.large_volume_limits, theta)
        return s, l, ds, dl
