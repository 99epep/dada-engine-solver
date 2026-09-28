"""Finite centered/offset slider-crank laws extracted from the K2 evaluation.

Lengths use crank-radius units; phase is theta + phase_rad.
"""
from dataclasses import dataclass
import math
from dada_solver.geometry import CylinderVolumeLimits

@dataclass(frozen=True)
class SliderMotion:
    rod_over_crank: float
    offset_over_crank: float
    phase_rad: float
    volume_increases_with_coordinate: bool

    def __post_init__(self):
        if type(self.volume_increases_with_coordinate) is not bool:
            raise ValueError('Volume direction must be an explicit boolean.')
        if not all(math.isfinite(v) for v in (self.rod_over_crank,self.offset_over_crank,self.phase_rad)) or self.rod_over_crank<=1+abs(self.offset_over_crank):
            raise ValueError('Invalid slider-crank closure')

    def normalized(self, theta):
        l,e=self.rod_over_crank,self.offset_over_crank
        a=theta+self.phase_rad
        z=e-math.sin(a); root=math.sqrt(l*l-z*z)
        x=math.cos(a)+root; dx=-math.sin(a)+z*math.cos(a)/root
        low=math.sqrt((l-1)**2-e*e); high=math.sqrt((l+1)**2-e*e)
        q=(x-low)/(high-low); dq=dx/(high-low)
        return (q,dq) if self.volume_increases_with_coordinate else (1-q,-dq)

    @property
    def stroke_over_crank(self):
        l,e=self.rod_over_crank,self.offset_over_crank
        return math.sqrt((l+1)**2-e*e)-math.sqrt((l-1)**2-e*e)

    def normalized_second_derivative(self, theta):
        l,e=self.rod_over_crank,self.offset_over_crank
        a=theta+self.phase_rad
        z=e-math.sin(a); c=math.cos(a); s=math.sin(a)
        root=math.sqrt(l*l-z*z)
        ddx=-c+(-c*c-z*s)/root-z*z*c*c/root**3
        value=ddx/self.stroke_over_crank
        return value if self.volume_increases_with_coordinate else -value

    def joint_state(self, theta):
        """Original coordinate convention, reconstructed from the normalized law."""
        q,_=self.normalized(theta)
        if not self.volume_increases_with_coordinate: q=1-q
        low=math.sqrt((self.rod_over_crank-1)**2-self.offset_over_crank**2)
        angle=theta+self.phase_rad
        return dict(A=(0.,0.),B=(math.cos(angle),math.sin(angle)),
                    P=(low+q*self.stroke_over_crank,self.offset_over_crank))


@dataclass(frozen=True)
class SliderCrankKinematics:
    small: SliderMotion
    large: SliderMotion
    small_volume_limits: CylinderVolumeLimits
    large_volume_limits: CylinderVolumeLimits
    small_physical_stroke = None
    large_physical_stroke = None

    def breakpoint_angles(self): return ()
    def _side(self, side, theta):
        q,dq=getattr(self,side).normalized(theta)
        limits=getattr(self,side+'_volume_limits')
        return limits.minimum+limits.swept*q,limits.swept*dq
    def small_cylinder_volume(self,theta): return self._side('small',theta)[0]
    def large_cylinder_volume(self,theta): return self._side('large',theta)[0]
    def small_cylinder_volume_derivative(self,theta): return self._side('small',theta)[1]
    def large_cylinder_volume_derivative(self,theta): return self._side('large',theta)[1]
    def small_cylinder_volume_second_derivative(self,theta):
        return self.small.normalized_second_derivative(theta)*self.small_volume_limits.swept
    def large_cylinder_volume_second_derivative(self,theta):
        return self.large.normalized_second_derivative(theta)*self.large_volume_limits.swept
    def cylinder_volumes_and_derivatives(self,theta):
        s,ds=self._side('small',theta);l,dl=self._side('large',theta)
        return s,l,ds,dl
