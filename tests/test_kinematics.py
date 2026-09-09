import math

import pytest

from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.kinematics import HarmonicVolumeKinematics, IdealPiecewiseLinearVolumeKinematics


def test_harmonic_example_respects_limits_and_large_origin() -> None:
    small = CylinderVolumeLimits(1.0e-4, 3.0e-4)
    large = CylinderVolumeLimits(2.0e-4, 6.0e-4)
    kinematics = HarmonicVolumeKinematics(small, large, math.pi)

    assert kinematics.large_cylinder_volume(0.0) == pytest.approx(large.maximum)
    assert kinematics.large_cylinder_volume(math.pi) == pytest.approx(large.minimum)
    assert kinematics.small_cylinder_volume(0.0) == pytest.approx(small.minimum)
    assert kinematics.large_cylinder_volume_derivative(0.0) == pytest.approx(0.0)


def test_ideal_piecewise_linear_kinematics_respects_lambda_and_four_sectors() -> None:
    small = CylinderVolumeLimits(1.0, 5.0)
    large = CylinderVolumeLimits(2.0, 12.0)
    kinematics = IdealPiecewiseLinearVolumeKinematics(small, large, 0.7, 0.7, 0.25)
    quarter = 0.5 * math.pi

    assert kinematics.small_cylinder_volume(0.0) == pytest.approx(small.minimum)
    assert kinematics.large_cylinder_volume(0.0) == pytest.approx(large.maximum)
    assert kinematics.small_cylinder_volume(quarter) == pytest.approx(small.minimum)
    assert large.closure_fraction(kinematics.large_cylinder_volume(quarter)) == pytest.approx(0.7)
    assert small.closure_fraction(kinematics.small_cylinder_volume(2.0 * quarter)) == pytest.approx(0.7)
    assert kinematics.large_cylinder_volume(2.0 * quarter) == pytest.approx(large.minimum)
    assert kinematics.small_cylinder_volume(3.0 * quarter) == pytest.approx(small.maximum)
    assert kinematics.large_cylinder_volume(3.0 * quarter) == pytest.approx(large.minimum)
    assert kinematics.breakpoint_angles() == pytest.approx((quarter, 2.0 * quarter, 3.0 * quarter))


def test_longer_transfer_sectors_preserve_periodicity_and_targets() -> None:
    limits = CylinderVolumeLimits(1.0, 2.0)
    kinematics = IdealPiecewiseLinearVolumeKinematics(limits, limits, 0.7, 0.7, 0.15)

    assert kinematics.transfer_sector_fraction == pytest.approx(0.35)
    assert kinematics.small_cylinder_volume(0.0) == pytest.approx(
        kinematics.small_cylinder_volume(2.0 * math.pi)
    )
