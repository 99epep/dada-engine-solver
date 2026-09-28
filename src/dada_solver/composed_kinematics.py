"""Compose independently selected cylinder laws without imposing motor direction."""
from dataclasses import dataclass


@dataclass(frozen=True)
class CylinderLaw:
    """One side of an existing backend, with explicit optional physical scaling."""
    model: object
    side: str
    physical_stroke: float | None = None

    @property
    def limits(self):
        return getattr(self.model,self.side+'_volume_limits')

    def value(self, theta, order=0):
        suffix=('', '_derivative', '_second_derivative')[order]
        method=getattr(self.model,self.side+'_cylinder_volume'+suffix,None)
        if method is None:
            raise NotImplementedError('Analytic acceleration is unavailable for this family.')
        return method(theta)


@dataclass(frozen=True)
class ComposedKinematics:
    small: CylinderLaw
    large: CylinderLaw
    mechanical_diagnostics: tuple = ()
    mechanical_metrics: tuple = ()

    @property
    def small_volume_limits(self): return self.small.limits
    @property
    def large_volume_limits(self): return self.large.limits
    @property
    def small_physical_stroke(self): return self.small.physical_stroke
    @property
    def large_physical_stroke(self): return self.large.physical_stroke
    @property
    def diagnostics(self): return ()

    def require_feasible(self):
        for law in (self.small,self.large):
            check=getattr(law.model,'require_feasible',None)
            if check is not None: check()

    def small_cylinder_volume(self,theta): return self.small.value(theta)
    def large_cylinder_volume(self,theta): return self.large.value(theta)
    def small_cylinder_volume_derivative(self,theta): return self.small.value(theta,1)
    def large_cylinder_volume_derivative(self,theta): return self.large.value(theta,1)
    def small_cylinder_volume_second_derivative(self,theta): return self.small.value(theta,2)
    def large_cylinder_volume_second_derivative(self,theta): return self.large.value(theta,2)

    def cylinder_volumes_and_derivatives(self,theta):
        return (self.small.value(theta),self.large.value(theta),
                self.small.value(theta,1),self.large.value(theta,1))

    def breakpoint_angles(self):
        return tuple(sorted(set(self.small.model.breakpoint_angles()) |
                            set(self.large.model.breakpoint_angles())))
