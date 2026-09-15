"""Independent four-stage chronology, phase release and breakpoints."""

import math
from pathlib import Path

import numpy as np
import pytest

from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_model
from dada_solver.four_stage_kinematics import FourStageVolumeKinematics
from dada_solver.independent_four_stage_kinematics import (
    IndependentFourStageVolumeKinematics,
)


TAU = 2.0 * math.pi


@pytest.fixture
def limits():
    config = load_simulation_configuration(
        Path(__file__).resolve().parents[1]
        / "examples"
        / "motor_demonstrator_original_325c.toml"
    )
    return config, config.machine_volumes


def test_shared_seed_reproduces_existing_four_stage_law(limits):
    config, volumes = limits
    shared = FourStageVolumeKinematics(
        volumes.small_cylinder,
        volumes.large_cylinder,
        0.30,
        0.55,
        0.80,
        0.40,
        0.70,
        0.60,
        0.20,
    )
    independent = IndependentFourStageVolumeKinematics(
        volumes.small_cylinder,
        volumes.large_cylinder,
        t0_s=0.0,
        t1_l=0.30,
        t2_l=0.55,
        t3_l=0.80,
        t1_s=0.30,
        t2_s=0.55,
        t3_s=0.80,
        a_l=0.40,
        b_l=0.70,
        a_s=0.60,
        b_s=0.20,
    )

    shared_model = build_model(config, kinematics=shared)
    independent_model = build_model(config, kinematics=independent)

    for theta in np.linspace(0.0, TAU, 101):
        np.testing.assert_allclose(
            shared_model.kinematics.cylinder_volumes_and_derivatives(theta),
            independent_model.kinematics.cylinder_volumes_and_derivatives(theta),
            rtol=0.0,
            atol=1e-14,
        )


def test_small_origin_and_knots_are_independent(limits):
    config, volumes = limits
    independent = IndependentFourStageVolumeKinematics(
        volumes.small_cylinder,
        volumes.large_cylinder,
        t0_s=0.10,
        t1_l=0.30,
        t2_l=0.55,
        t3_l=0.80,
        t1_s=0.25,
        t2_s=0.50,
        t3_s=0.75,
        a_l=0.40,
        b_l=0.70,
        a_s=0.60,
        b_s=0.20,
    )
    k = build_model(config, kinematics=independent).kinematics

    s = volumes.small_cylinder
    l = volumes.large_cylinder

    def value(limits_, fraction):
        return limits_.minimum + fraction * limits_.swept

    for t, expected in (
        (0.10, 0.20),
        (0.35, 1.00),
        (0.60, 0.60),
        (0.85, 0.00),
        (1.10, 0.20),
    ):
        assert k.small_cylinder_volume(t * TAU) == pytest.approx(
            value(s, expected)
        )

    for t, expected in (
        (0.00, 1.00),
        (0.30, 0.70),
        (0.55, 0.00),
        (0.80, 0.40),
        (1.00, 1.00),
    ):
        assert k.large_cylinder_volume(t * TAU) == pytest.approx(
            value(l, expected)
        )

    expected_breaks = np.array(
        [0.10, 0.30, 0.35, 0.55, 0.60, 0.80, 0.85]
    ) * TAU
    np.testing.assert_allclose(k.breakpoint_angles(), expected_breaks)


def test_derivatives_match_finite_difference_away_from_breakpoints(limits):
    config, volumes = limits
    independent = IndependentFourStageVolumeKinematics(
        volumes.small_cylinder,
        volumes.large_cylinder,
        t0_s=-0.08,
        t1_l=0.31,
        t2_l=0.57,
        t3_l=0.82,
        t1_s=0.28,
        t2_s=0.53,
        t3_s=0.77,
        a_l=0.42,
        b_l=0.66,
        a_s=0.58,
        b_s=0.24,
    )
    k = build_model(config, kinematics=independent).kinematics

    for side in ("small", "large"):
        volume = getattr(k, f"{side}_cylinder_volume")
        derivative = getattr(k, f"{side}_cylinder_volume_derivative")
        for t in (0.03, 0.18, 0.44, 0.72, 0.94):
            theta = t * TAU
            h = 1e-6
            numerical = (volume(theta + h) - volume(theta - h)) / (2.0 * h)
            assert derivative(theta) == pytest.approx(numerical, rel=1e-8)


@pytest.mark.parametrize(
    "updates",
    [
        {"t0_s": 0.5},
        {"t0_s": -0.5001},
        {"t1_l": 0.0},
        {"t2_l": 0.1},
        {"t3_s": 1.0},
        {"a_s": -0.01},
        {"b_l": float("nan")},
    ],
)
def test_rejects_invalid_parameterization(limits, updates):
    _, volumes = limits
    values = dict(
        t0_s=0.0,
        t1_l=0.30,
        t2_l=0.55,
        t3_l=0.80,
        t1_s=0.30,
        t2_s=0.55,
        t3_s=0.80,
        a_l=0.40,
        b_l=0.70,
        a_s=0.60,
        b_s=0.20,
    )
    values.update(updates)
    with pytest.raises(ValueError):
        IndependentFourStageVolumeKinematics(
            volumes.small_cylinder,
            volumes.large_cylinder,
            **values,
        )
