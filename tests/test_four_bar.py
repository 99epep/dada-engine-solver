import math

import numpy as np
import pytest

from dada_solver.four_bar import (
    CouplerOutputPoint,
    FourBarLoop,
    FourBarSliderAssembly,
    RockerOutputPoint,
    SharedCrankRockerDesign,
    SharedCrankFourBarVolumeKinematics,
    SliderConstraint,
    published_f65_opposed_kinematics,
    published_e0_opposed_kinematics,
    shared_crank_rocker_kinematics,
    slider_rod_sensitivity,
)
from dada_solver.geometry import CylinderVolumeLimits


def assembly(pivot_x: float, pivot_y: float, branch: int = 1):
    return FourBarSliderAssembly(
        loop=FourBarLoop(
            coupler_length=0.12,
            rocker_length=0.11,
            rocker_pivot_x=pivot_x,
            rocker_pivot_y=pivot_y,
            assembly_branch=branch,
        ),
        output=RockerOutputPoint(along_rocker=0.14),
        slider=SliderConstraint(
            axis_origin_x=0.05,
            axis_origin_y=0.16 if branch > 0 else -0.16,
            axis_angle=0.0,
            connecting_rod_length=0.50,
            assembly_branch=1,
        ),
    )


def kinematics():
    return SharedCrankFourBarVolumeKinematics(
        crank_radius=0.03,
        small_assembly=assembly(0.10, 0.0, 1),
        large_assembly=assembly(-0.11, 0.01, -1),
        small_volume_limits=CylinderVolumeLimits(1.0e-4, 5.0e-4),
        large_volume_limits=CylinderVolumeLimits(2.0e-4, 1.0e-3),
        large_volume_increases_with_coordinate=False,
    )


def test_four_bar_closure_and_shared_crank_pin() -> None:
    model = kinematics()
    small, large = model.slider_states(0.73)

    np.testing.assert_allclose(small.crank_pin, large.crank_pin)
    for state, item in (
        (small, model.small_assembly),
        (large, model.large_assembly),
    ):
        pin = np.asarray(state.crank_pin)
        joint = np.asarray(state.coupler_joint)
        pivot = np.array((item.loop.rocker_pivot_x, item.loop.rocker_pivot_y))
        assert np.linalg.norm(joint - pin) == pytest.approx(item.loop.coupler_length)
        assert np.linalg.norm(joint - pivot) == pytest.approx(item.loop.rocker_length)


def test_general_design_supports_independent_loops_on_shared_crank() -> None:
    limits = CylinderVolumeLimits(1.0e-4, 5.0e-4)
    design = SharedCrankRockerDesign(
        crank_ratio=0.25,
        small_coupler_ratio=0.90,
        small_rocker_ratio=0.85,
        small_output_along_ratio=0.72,
        small_output_normal_ratio=-0.18,
        large_coupler_ratio=1.05,
        large_rocker_ratio=0.95,
        large_output_along_ratio=0.81,
        large_output_normal_ratio=0.22,
    )
    model = shared_crank_rocker_kinematics(
        design, limits, limits, ground_distance=0.1
    )
    small, large = model.slider_states(0.91)

    np.testing.assert_allclose(small.crank_pin, large.crank_pin)
    assert model.small_assembly.loop.coupler_length == pytest.approx(0.09)
    assert model.large_assembly.loop.coupler_length == pytest.approx(0.105)
    assert model.small_physical_stroke > 0.0
    assert model.large_physical_stroke > 0.0


def test_slider_derivative_matches_centered_difference() -> None:
    model = kinematics()
    angle = 1.37
    step = 1.0e-6
    state = model.slider_states(angle)[0]
    difference = (
        model.slider_states(angle + step)[0].coordinate
        - model.slider_states(angle - step)[0].coordinate
    ) / (2.0 * step)

    assert state.coordinate_derivative == pytest.approx(difference, rel=2.0e-8)


def test_published_f65_coupler_point_and_derivative() -> None:
    """Exercise the published F65 geometry without inventing its slider layout."""
    scale = 1.0e-3
    item = FourBarSliderAssembly(
        loop=FourBarLoop(
            coupler_length=117.91 * scale,
            rocker_length=141.73 * scale,
            rocker_pivot_x=100.0 * scale,
            rocker_pivot_y=0.0,
        ),
        output=CouplerOutputPoint(
            along_coupler=94.4 * scale,
            normal_to_coupler=-88.65 * scale,
        ),
        # This long rod is only a test fixture, not part of the published geometry.
        slider=SliderConstraint(0.0, -0.30, 0.0, 0.50),
    )
    angle = 1.13
    step = 1.0e-6
    crank = lambda value: (
        44.96 * scale * np.array((math.cos(value), math.sin(value))),
        44.96 * scale * np.array((-math.sin(value), math.cos(value))),
    )
    state = item.evaluate(*crank(angle))
    before = item.evaluate(*crank(angle - step)).coordinate
    after = item.evaluate(*crank(angle + step)).coordinate

    pin = np.asarray(state.crank_pin)
    joint = np.asarray(state.coupler_joint)
    unit = (joint - pin) / (117.91 * scale)
    expected_point = pin + 94.4 * scale * unit - 88.65 * scale * np.array(
        (-unit[1], unit[0])
    )
    np.testing.assert_allclose(state.output_point, expected_point, atol=1.0e-14)
    assert state.coordinate_derivative == pytest.approx(
        (after - before) / (2.0 * step), rel=2.0e-8
    )


def test_volume_laws_are_periodic_and_respect_configured_extrema() -> None:
    model = kinematics()
    angles = np.linspace(0.0, 2.0 * math.pi, 4001)
    small = np.array([model.small_cylinder_volume(angle) for angle in angles])
    large = np.array([model.large_cylinder_volume(angle) for angle in angles])

    assert small[0] == pytest.approx(small[-1], abs=1.0e-14)
    assert large[0] == pytest.approx(large[-1], abs=1.0e-14)
    assert np.min(small) == pytest.approx(1.0e-4, abs=2.0e-10)
    assert np.max(small) == pytest.approx(5.0e-4, abs=2.0e-10)
    assert np.min(large) == pytest.approx(2.0e-4, abs=2.0e-10)
    assert np.max(large) == pytest.approx(1.0e-3, abs=2.0e-10)


def test_volume_derivative_matches_centered_difference() -> None:
    model = kinematics()
    angle = 2.11
    step = 1.0e-6
    expected = (
        model.small_cylinder_volume(angle + step)
        - model.small_cylinder_volume(angle - step)
    ) / (2.0 * step)

    assert model.small_cylinder_volume_derivative(angle) == pytest.approx(
        expected, rel=2.0e-8
    )


def test_toggle_configuration_is_rejected() -> None:
    item = FourBarSliderAssembly(
        loop=FourBarLoop(0.10, 0.10, 0.20, 0.0),
        output=RockerOutputPoint(0.10),
        slider=SliderConstraint(0.0, 0.0, 0.0, 0.20),
    )

    with pytest.raises(ValueError, match="toggle"):
        item.evaluate(np.array((0.0, 0.0)), np.array((0.0, 0.1)))


def test_published_f65_opposed_pair_shares_pin_and_has_parallel_sliders() -> None:
    limits = CylinderVolumeLimits(1.0e-5, 1.0e-4)
    model = published_f65_opposed_kinematics(
        limits,
        limits,
        ground_distance=0.10,
        connecting_rod_to_projected_stroke_ratio=5.0,
    )
    small, large = model.slider_states(0.81)

    np.testing.assert_allclose(small.crank_pin, large.crank_pin)
    assert model.small_assembly.slider.axis_angle == pytest.approx(0.0)
    assert model.large_assembly.slider.axis_angle == pytest.approx(0.0)
    assert model.small_assembly.slider.axis_origin_y == pytest.approx(
        -model.large_assembly.slider.axis_origin_y
    )
    angles = np.linspace(0.0, 2.0 * math.pi, 1441)
    small = np.array([model.small_cylinder_volume(angle) for angle in angles])
    large = np.array([model.large_cylinder_volume(angle) for angle in angles])
    small_minimum_angle = angles[int(np.argmin(small))]
    large_maximum_angle = angles[int(np.argmax(large))]
    separation = abs(
        (small_minimum_angle - large_maximum_angle + math.pi) % (2.0 * math.pi)
        - math.pi
    )
    assert math.degrees(separation) < 20.0


def test_longer_f65_piston_rod_reduces_projection_error() -> None:
    limits = CylinderVolumeLimits(1.0e-5, 1.0e-4)
    errors = []
    for ratio in (3.0, 5.0, 8.0, 12.0):
        model = published_f65_opposed_kinematics(
            limits,
            limits,
            ground_distance=0.10,
            connecting_rod_to_projected_stroke_ratio=ratio,
        )
        errors.append(
            slider_rod_sensitivity(
                model.small_assembly, model.crank_radius
            ).maximum_normalized_shape_error
        )

    assert all(left > right for left, right in zip(errors, errors[1:]))


def test_published_e0_opposed_pair_has_expected_geometry_and_motion() -> None:
    limits = CylinderVolumeLimits(1.0e-5, 1.0e-4)
    model = published_e0_opposed_kinematics(
        limits,
        limits,
        ground_distance=0.10,
        connecting_rod_to_projected_stroke_ratio=5.0,
    )

    assert model.crank_radius == pytest.approx(0.073)
    assert model.small_assembly.loop.coupler_length == pytest.approx(0.10)
    assert model.small_assembly.loop.rocker_length == pytest.approx(0.1103728)
    small_state, large_state = model.slider_states(0.83)
    np.testing.assert_allclose(small_state.crank_pin, large_state.crank_pin)
    angles = np.linspace(0.0, 2.0 * math.pi, 1441)
    small = np.array([model.small_cylinder_volume(angle) for angle in angles])
    large = np.array([model.large_cylinder_volume(angle) for angle in angles])
    separation = abs(
        (angles[int(np.argmin(small))] - angles[int(np.argmax(large))] + math.pi)
        % (2.0 * math.pi)
        - math.pi
    )
    assert math.degrees(separation) == pytest.approx(37.75, abs=0.3)
