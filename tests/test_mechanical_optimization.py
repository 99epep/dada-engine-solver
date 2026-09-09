import pytest

from dada_solver.four_bar import shared_crank_rocker_kinematics
from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.mechanical_optimization import (
    ReflectedRockerParameters,
    scale_for_large_bore_to_stroke_ratio,
)


def test_mechanism_scale_can_be_selected_from_large_piston_proportions() -> None:
    small = CylinderVolumeLimits(0.0002772, 0.0279972)
    large = CylinderVolumeLimits(0.00033, 0.03333)
    mechanism = shared_crank_rocker_kinematics(
        ReflectedRockerParameters(0.4, 1.0, 1.17, -56.0, 304.375).design(),
        small,
        large,
        ground_distance=0.1,
    )

    scaled = scale_for_large_bore_to_stroke_ratio(
        mechanism,
        current_ground_distance=0.1,
        target_large_bore_to_stroke_ratio=1.0,
    )

    assert scaled.large_bore == pytest.approx(scaled.large_stroke)
    assert scaled.ground_distance == pytest.approx(0.384, rel=0.01)
    assert scaled.small_bore_to_stroke_ratio < 1.0
