"""Four shared linear stages specified in forward motor-cycle time."""
from dataclasses import dataclass
import math
from bisect import bisect_right
from dada_solver.geometry import CylinderVolumeLimits


@dataclass(frozen=True)
class FourStageVolumeKinematics:
    """Seven independent parameters; zero level means minimum enclosed volume.

    Motor-time knots: L=(1, b_l, 0, a_l), S=(b_s, 1, a_s, 0).
    This object adapts motor time to study angle (-theta); build_model then
    applies its usual motor reversal once. It must be used for motor operation.
    """
    small_volume_limits: CylinderVolumeLimits
    large_volume_limits: CylinderVolumeLimits
    t1: float
    t2: float
    t3: float
    a_l: float
    b_l: float
    a_s: float
    b_s: float

    def __post_init__(self):
        if not 0 < self.t1 < self.t2 < self.t3 < 1:
            raise ValueError('Require 0 < t1 < t2 < t3 < 1.')
        if any(not math.isfinite(v) or not 0 <= v <= 1 for v in (self.a_l,self.b_l,self.a_s,self.b_s)):
            raise ValueError('Intermediate swept-volume fractions must lie in [0, 1].')

    @property
    def small_physical_stroke(self): return None
    @property
    def large_physical_stroke(self): return None

    def breakpoint_angles(self):
        return tuple(2*math.pi*(1-t) for t in (self.t3,self.t2,self.t1))

    def cylinder_volumes_and_derivatives(self, theta):
        if not math.isfinite(theta): raise ValueError('Angle must be finite.')
        t = (-theta/(2*math.pi)) % 1
        knots = (0,self.t1,self.t2,self.t3,1)
        i = min(bisect_right(knots,t)-1,3)
        width = knots[i+1]-knots[i]
        values = []
        derivatives = []
        for q, limits in (((self.b_s,1,self.a_s,0,self.b_s),self.small_volume_limits),
                          ((1,self.b_l,0,self.a_l,1),self.large_volume_limits)):
            slope = (q[i+1]-q[i])/width
            values.append(limits.minimum+(q[i]+slope*(t-knots[i]))*limits.swept)
            derivatives.append(-slope*limits.swept/(2*math.pi))
        return *values,*derivatives

    def small_cylinder_volume(self, theta): return self.cylinder_volumes_and_derivatives(theta)[0]
    def large_cylinder_volume(self, theta): return self.cylinder_volumes_and_derivatives(theta)[1]
    def small_cylinder_volume_derivative(self, theta): return self.cylinder_volumes_and_derivatives(theta)[2]
    def large_cylinder_volume_derivative(self, theta): return self.cylinder_volumes_and_derivatives(theta)[3]
