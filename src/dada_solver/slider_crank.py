"""Finite centered/offset slider-crank laws extracted from the K2 evaluation.

Lengths use crank-radius units; phase is theta + phase_rad.
"""
from dataclasses import dataclass
import math
import numpy as np
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
        operations = np if np.ndim(theta) else math
        z=e-operations.sin(a); root=operations.sqrt(l*l-z*z)
        x=operations.cos(a)+root; dx=-operations.sin(a)+z*operations.cos(a)/root
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
        operations = np if np.ndim(theta) else math
        z=e-operations.sin(a); c=operations.cos(a); s=operations.sin(a)
        root=operations.sqrt(l*l-z*z)
        ddx=-c+(-c*c-z*s)/root-z*z*c*c/root**3
        value=ddx/self.stroke_over_crank
        return value if self.volume_increases_with_coordinate else -value

    def stationary_points(self):
        """The two exact collinear crank/rod dead centers in study angle."""
        maximum = (math.asin(self.offset_over_crank/(self.rod_over_crank+1))-self.phase_rad) % (2*math.pi)
        minimum = (math.pi+math.asin(self.offset_over_crank/(self.rod_over_crank-1))-self.phase_rad) % (2*math.pi)
        points = [dict(angle_rad=maximum,kind='maximum' if self.volume_increases_with_coordinate else 'minimum'),
                  dict(angle_rad=minimum,kind='minimum' if self.volume_increases_with_coordinate else 'maximum')]
        return tuple(sorted(points,key=lambda p:p['angle_rad']))

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
