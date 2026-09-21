"""Experimental interpolants derive consistent periodic velocities from volume."""
import sys
from pathlib import Path
import math
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'examples'))
from experimental_periodic_interpolation import PeriodicVolumeInterpolant


class Smooth:
    def cylinder_volumes_and_derivatives(self,a):return 2+math.cos(a),3+math.sin(a),-math.sin(a),math.cos(a)
    def breakpoint_angles(self):return ()


def test_periodicity_derivative_consistency_and_resolution_convergence():
    coarse=PeriodicVolumeInterpolant(Smooth(),32);fine=PeriodicVolumeInterpolant(Smooth(),64)
    assert fine.cylinder_volumes_and_derivatives(0)==fine.cylinder_volumes_and_derivatives(2*math.pi)
    angles=(np.arange(4099)+.5)*2*math.pi/4099
    errors=[]
    for p in (coarse,fine):
        errors.append(np.max(np.abs(np.array([p.cylinder_volumes_and_derivatives(a) for a in angles])-
                                   np.array([Smooth().cylinder_volumes_and_derivatives(a) for a in angles])),axis=0))
    assert np.all(errors[1][:2]<errors[0][:2]/10)
    assert np.all(errors[1][2:]<errors[0][2:]/6)
    a=.54321;h=1e-5
    before=np.array(fine.cylinder_volumes_and_derivatives(a-h))
    after=np.array(fine.cylinder_volumes_and_derivatives(a+h))
    np.testing.assert_allclose((after[:2]-before[:2])/(2*h),fine.cylinder_volumes_and_derivatives(a)[2:],atol=1e-10,rtol=0)


def test_declared_discontinuities_are_rejected():
    class Piecewise(Smooth):
        def breakpoint_angles(self):return (math.pi,)
    with pytest.raises(ValueError,match='breakpoints'):PeriodicVolumeInterpolant(Piecewise(),32)
