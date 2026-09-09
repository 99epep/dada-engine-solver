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


@pytest.mark.parametrize("configuration", ["motor_demonstrator_lower_lambda_trial.toml", "motor_demonstrator_opposed_lambda_trial.toml"])
def test_independent_four_bar_trial_preserves_volume_limits(configuration):
    import math
    from pathlib import Path
    from dada_solver.configuration import load_simulation_configuration
    from dada_solver.factory import build_model
    root=Path(__file__).resolve().parents[1]
    config=load_simulation_configuration(root/'examples'/configuration)
    design=config.shared_four_bar_design
    assert abs(design.small_output_normal_ratio) != abs(design.large_output_normal_ratio)
    model=build_model(config)
    for angle in np.linspace(0,2*math.pi,721):
        v=model.volumes(float(angle))
        for side in ['small','large']:
            limits=getattr(config.machine_volumes,side+'_cylinder')
            volume=getattr(v,side+'_cylinder')
            assert limits.minimum*(1-1e-8)<=volume<=limits.maximum*(1+1e-8)
