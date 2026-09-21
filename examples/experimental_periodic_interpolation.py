"""Experimental smooth-motion interpolation; never selected by production code."""
import bisect
import math
import numpy as np
from scipy.interpolate import CubicSpline

TAU=2*math.pi


class PeriodicVolumeInterpolant:
    """Derive volume and angular velocity from the same periodic cubic.

    The caller must establish a smooth fixed assembly branch. Absence of declared
    breakpoints alone is not proof of smoothness. The exact provider remains the
    source of all coefficients. No approximate path is registered in campaigns.
    """
    def __init__(self, forward, intervals):
        if intervals < 8 or intervals != int(intervals):
            raise ValueError('At least eight integer intervals are required.')
        if tuple(getattr(forward,'breakpoint_angles',lambda:())()):
            raise ValueError('Cannot interpolate across declared kinematic breakpoints.')
        self.forward=forward
        self.intervals=int(intervals)
        grid=np.linspace(0.,TAU,self.intervals+1)
        exact=np.array([forward.cylinder_volumes_and_derivatives(float(a)) for a in grid])
        if not np.all(np.isfinite(exact)) or np.any(exact[:,:2]<=0):
            raise ValueError('Exact mechanism must produce finite positive volumes.')
        scale=np.maximum(np.max(np.abs(exact),axis=0),np.finfo(float).tiny)
        if np.any(np.abs(exact[-1]-exact[0])>1e-10*scale):
            raise ValueError('Exact mechanism does not close periodically.')
        values=exact[:,:2].copy();values[-1]=values[0]
        spline=CubicSpline(grid,values,bc_type='periodic',axis=0)
        self.knots=tuple(float(x) for x in grid)
        self.coefficients=tuple(tuple(tuple(float(v) for v in spline.c[:,j,k]) for k in (0,1))
                                for j in range(self.intervals))
        # Check polynomial extrema, not just positive grid values.
        for k in (0,1):
            scalar=CubicSpline(grid,values[:,k],bc_type='periodic')
            roots=scalar.derivative().roots(extrapolate=False)
            roots=roots[np.isfinite(roots)&(roots>=0)&(roots<=TAU)]
            if np.min(scalar(np.r_[grid,roots]))<=0:
                raise ValueError('Interpolated volume is not strictly positive.')

    def __getattr__(self,name):return getattr(self.forward,name)

    def cylinder_volumes_and_derivatives(self,theta):
        if not math.isfinite(theta):raise ValueError('Cycle angle must be finite.')
        angle=theta%TAU
        j=bisect.bisect_right(self.knots,angle)-1
        dx=angle-self.knots[j]
        a,b,c,d=self.coefficients[j][0]
        v_s=((a*dx+b)*dx+c)*dx+d;d_s=(3*a*dx+2*b)*dx+c
        a,b,c,d=self.coefficients[j][1]
        v_l=((a*dx+b)*dx+c)*dx+d;d_l=(3*a*dx+2*b)*dx+c
        return v_s,v_l,d_s,d_l

    def small_cylinder_volume(self,a):return self.cylinder_volumes_and_derivatives(a)[0]
    def large_cylinder_volume(self,a):return self.cylinder_volumes_and_derivatives(a)[1]
    def small_cylinder_volume_derivative(self,a):return self.cylinder_volumes_and_derivatives(a)[2]
    def large_cylinder_volume_derivative(self,a):return self.cylinder_volumes_and_derivatives(a)[3]
    def breakpoint_angles(self):return ()
