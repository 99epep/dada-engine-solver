import pytest

from dada_solver.geometry import CylinderVolumeLimits


def test_closure_fraction_uses_published_definition() -> None:
    limits = CylinderVolumeLimits(minimum=1.0e-4, maximum=5.0e-4)

    assert limits.swept == pytest.approx(4.0e-4)
    assert limits.closure_fraction(limits.maximum) == pytest.approx(0.0)
    assert limits.closure_fraction(limits.minimum) == pytest.approx(1.0)


def test_invalid_cylinder_limits_are_rejected() -> None:
    with pytest.raises(ValueError, match="must exceed"):
        CylinderVolumeLimits(minimum=2.0e-4, maximum=2.0e-4)

