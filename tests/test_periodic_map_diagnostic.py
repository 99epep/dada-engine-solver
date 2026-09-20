"""Validate the diagnostic coordinate reduction independently of the engine."""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'examples'))
from diagnose_periodic_map import MASS, tangent_coordinates, central_jacobian, ExactAngleProvider


def test_mass_constraint_and_known_linear_spectrum():
    anchor = np.array([.002, 500., .03, 8000., .0002, 50., .0001, 30., 10000., 20000.])
    scales, basis = tangent_coordinates(anchor)
    np.testing.assert_allclose(basis.T@basis, np.eye(9), atol=1e-15)
    np.testing.assert_allclose((scales[:, None]*basis)[MASS].sum(axis=0), 0., atol=1e-17)
    expected = np.diag(np.linspace(.05, .95, 9))
    def physical_map(x):
        return anchor+scales*(basis@(expected@(basis.T@((x-anchor)/scales))))
    def reduced_map(z):
        return basis.T@((physical_map(anchor+scales*(basis@z))-anchor)/scales)
    measured = central_jacobian(reduced_map, 9, 1e-3)
    np.testing.assert_allclose(measured, expected, atol=3e-13)
    np.testing.assert_allclose(np.sort(np.linalg.eigvals(measured)), np.diag(expected), atol=3e-13)


def test_cache_reuses_only_identical_consecutive_angles():
    class Provider:
        def __init__(self): self.calls = 0
        def cylinder_volumes_and_derivatives(self, angle):
            self.calls += 1
            return angle, angle*2, 1., 2.
    forward = Provider()
    cache = ExactAngleProvider(forward, True)
    a = .5
    b = np.nextafter(a, 1.)
    for angle in [a, a, b, b, a]:
        assert cache.cylinder_volumes_and_derivatives(angle) == (angle, angle*2, 1., 2.)
    assert forward.calls == 3
    assert cache.repeats == 2
