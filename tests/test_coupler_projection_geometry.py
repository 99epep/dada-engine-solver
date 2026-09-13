import math

import numpy as np
import pytest

from dada_solver.coupler_projection import (
    load_shared_crank_coupler_projection_geometry,
)
from dada_solver.geometry import CylinderVolumeLimits


GEOMETRY = """
schema_version = 1
angle_convention = "study_angle_radians"
method = "projection_only"

[small]
coupler_ratio = 1.6563469344036357
rocker_ratio = 5.1468443338063725
rocker_pivot_x_ratio = -4.554542857952098
rocker_pivot_y_ratio = -2.117480304002033
output_along_ratio = 1.5159800513453636
output_normal_ratio = -0.4124220620237704
axis_angle_shared_crank_rad = -2.3048052319993504
assembly_branch = -1

[large]
coupler_ratio = 1.8978981818881382
rocker_ratio = 8.562558314779313
rocker_pivot_x_ratio = 6.5305460224395135
rocker_pivot_y_ratio = 5.217865580266911
output_along_ratio = 1.8306390741406238
output_normal_ratio = 0.33169924587809335
axis_angle_shared_crank_rad = 0.40824128311504243
assembly_branch = 1
"""


def test_synthesis_geometry_loads_as_periodic_shared_crank_kinematics(tmp_path):
    path = tmp_path / "geometry.toml"
    path.write_text(GEOMETRY)
    small = CylinderVolumeLimits(8.316831683168316e-6, 8.4e-4)
    large = CylinderVolumeLimits(9.900990099009901e-6, 1.0e-3)
    model = load_shared_crank_coupler_projection_geometry(path, small, large)

    at_zero = model.cylinder_volumes_and_derivatives(0.0)
    at_wrap = model.cylinder_volumes_and_derivatives(2.0 * math.pi)
    assert at_wrap == pytest.approx(at_zero, rel=0.0, abs=2.0e-18)

    angles = np.linspace(0.0, 2.0 * math.pi, 4001)
    small_values = np.array([model.small_cylinder_volume(x) for x in angles])
    large_values = np.array([model.large_cylinder_volume(x) for x in angles])
    assert small_values.min() == pytest.approx(small.minimum, abs=2.0e-9)
    assert small_values.max() == pytest.approx(small.maximum, abs=2.0e-9)
    assert large_values.min() == pytest.approx(large.minimum, abs=2.0e-9)
    assert large_values.max() == pytest.approx(large.maximum, abs=2.0e-9)


def test_synthesis_geometry_rejects_a_competing_angle_convention(tmp_path):
    path = tmp_path / "geometry.toml"
    path.write_text(GEOMETRY.replace("study_angle_radians", "motor_angle_radians"))
    limits = CylinderVolumeLimits(1.0e-6, 2.0e-6)
    with pytest.raises(ValueError, match="study_angle_radians"):
        load_shared_crank_coupler_projection_geometry(path, limits, limits)
