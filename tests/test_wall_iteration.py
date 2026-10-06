import numpy as np
import pytest
from dada_solver.exchangers.wall_iteration import extrapolate_wall_energies


def test_linear_contraction_predicts_fixed_point():
    result=extrapolate_wall_energies([100,200],[150,250],[175,275],[10,10])
    assert result == pytest.approx([200,300])


def test_diverging_or_flat_sequence_is_not_accelerated():
    assert extrapolate_wall_energies([100,200],[110,200],[130,200],[10,10]) is None


def test_change_is_bounded_and_capacity_checked():
    result=extrapolate_wall_energies([100,100],[110,110],[119.9,119.9],[1,1])
    assert np.max(result-119.9)<=30
    with pytest.raises(ValueError):
        extrapolate_wall_energies([1,1],[2,2],[3,3],[0,1])
