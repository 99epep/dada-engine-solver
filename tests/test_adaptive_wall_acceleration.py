"""Trend-sensitive proposals and rollback of unsuccessful wall guesses."""
from dataclasses import replace
from types import SimpleNamespace
import numpy as np
import pytest

from dada_solver.exchangers.wall_iteration import (
    AdaptiveWallAccelerationSettings, adaptive_wall_proposal)
from dada_solver.exchangers.wall_cycle import solve_periodic_wall_motor, WallCycleNumericalSettings
from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
from dada_solver.integration import IntegrationInterrupted


def test_regular_contraction_gets_partial_signed_jump():
    capacities=np.array([10.,20.])
    temperatures=np.array([[300,350],[304,346],[306,344],[307,343]])
    proposal=adaptive_wall_proposal(temperatures*capacities,capacities)
    assert proposal.correction_kelvin==pytest.approx((.99,-.99))
    assert proposal.wall_energies==pytest.approx(np.array([307.99,342.01])*capacities)
    assert proposal.confidence==pytest.approx((1,1))


@pytest.mark.parametrize('increments',[(4,3,2.9),(4,1,-1),(1,2,4),(1,.999,.998001),(1,1,1),(1e-9,5e-10,2.5e-10)])
def test_inflection_oscillation_growth_flatness_and_noise_are_not_boosted(increments):
    temperatures=300+np.r_[0,np.cumsum(increments)]
    assert adaptive_wall_proposal(np.column_stack((temperatures,temperatures)),[1,1]) is None


def test_ratio_drift_reduces_jump_confidence():
    regular=np.array([300,304,306,307])
    changing=np.array([300,304,306,307.1])
    stable=adaptive_wall_proposal(np.column_stack((regular,regular)),[1,1])
    curved=adaptive_wall_proposal(np.column_stack((changing,changing)),[1,1],
        settings=AdaptiveWallAccelerationSettings(minimum_confidence=.1))
    assert 0<curved.confidence[0]<stable.confidence[0]
    # Compare to the undamped latest-ratio prediction for the curved sequence.
    raw=1.1*.55/(1-.55)
    assert curved.correction_kelvin[0]<.99*raw


def test_each_wall_is_gated_independently():
    temperatures=np.array([[300,350],[304,354],[306,355],[307,354]])
    proposal=adaptive_wall_proposal(temperatures,[1,1])
    assert proposal.eligible==(True,False)
    assert proposal.wall_energies[1]==354


def test_jump_is_bounded_by_recent_speed_and_absolute_temperature():
    for increment,ratio in ((1,.97),(100,.9)):
        temperatures=300+np.r_[0,np.cumsum([increment,increment*ratio,increment*ratio**2])]
        proposal=adaptive_wall_proposal(np.column_stack((temperatures,temperatures)),[1,1])
        assert 0<proposal.correction_kelvin[0]<=min(30,20*increment*ratio**2)+1e-12


@pytest.mark.parametrize('kwargs',[dict(damping=0),dict(damping=1.1),dict(maximum_contraction_ratio=1),
    dict(relative_ratio_drift_limit=float('nan')),dict(maximum_residual_growth=.5),dict(recovery_cycles=3)])
def test_invalid_acceleration_controls(kwargs):
    with pytest.raises(ValueError): AdaptiveWallAccelerationSettings(**kwargs)


class ContractingWrapper:
    heat_in=SimpleNamespace(wall_capacity_j_k=1.)
    heat_out=SimpleNamespace(wall_capacity_j_k=1.)

    def __init__(self, failure=None):
        self.inputs=[]; self.last=None;self.failure=failure;self.failed=False

    def integrate_cycle(self,state,**kwargs):
        state=np.asarray(state).copy();self.inputs.append(state)
        boosted=self.last is not None and not np.array_equal(state,self.last)
        end=state.copy();end[8:10]=(state[8:10]+200)/2
        if boosted and self.failure and not self.failed:
            self.failed=True
            if self.failure=='domain': raise MicrotubeDomainError('Unsupported accelerated guess')
            if self.failure=='interrupt': raise IntegrationInterrupted('Deadline during boosted cycle')
            if self.failure=='growth': end[8:10]=state[8:10]+100
        self.last=end.copy()
        trajectory=np.zeros((15,2));trajectory[:10,0]=state;trajectory[:10,1]=end
        return np.array([0.,2*np.pi]),trajectory


def solve(wrapper,maximum_cycles=30,**kwargs):
    return solve_periodic_wall_motor(wrapper,np.r_[np.ones(8),100,100],maximum_cycles=maximum_cycles,
        settings=WallCycleNumericalSettings(accelerate_walls=True),
        adaptive_acceleration=AdaptiveWallAccelerationSettings(),**kwargs)


def test_acceleration_preserves_all_gas_values_and_requires_physical_convergence():
    wrapper=ContractingWrapper();result=solve(wrapper)
    assert result.converged and len(wrapper.inputs)>4
    assert result.history[3]['wall_initial_guess_extrapolated']
    assert result.history[4]['wall_acceleration_accepted']
    assert result.history[-1]['normalized_state_error']<=1
    for state in wrapper.inputs: np.testing.assert_array_equal(state[:8],np.ones(8))
    np.testing.assert_array_equal(result.last_complete_state,result.trajectory[:10,-1])


@pytest.mark.parametrize('failure',['growth','domain'])
def test_bad_jump_returns_to_anchor_and_counts_the_cost(failure):
    wrapper=ContractingWrapper(failure);callbacks=[]
    result=solve(wrapper,maximum_cycles=6,cycle_callback=lambda cycle,*_:callbacks.append(cycle))
    anchor=(wrapper.inputs[3].copy());anchor[8:10]=(anchor[8:10]+200)/2
    assert not np.array_equal(wrapper.inputs[4],anchor)
    np.testing.assert_array_equal(wrapper.inputs[5],anchor)
    assert callbacks==[1,2,3,4,6]
    assert result.status=='maximum_cycles'
    assert len(wrapper.inputs)==6
    assert any(r.get('phase')=='wall_acceleration_rollback' for r in result.solver_statistics)
    np.testing.assert_array_equal(result.last_complete_state,result.trajectory[:10,-1])


def test_interruption_during_boost_retains_preboost_physical_endpoint():
    wrapper=ContractingWrapper('interrupt');result=solve(wrapper)
    assert result.status=='interrupted' and len(result.history)==4
    np.testing.assert_array_equal(result.last_complete_state,result.trajectory[:10,-1])
    assert not np.array_equal(result.last_complete_state[8:10],wrapper.inputs[-1][8:10])


def test_last_budgeted_cycle_never_returns_an_extrapolated_guess():
    wrapper=ContractingWrapper();result=solve(wrapper,maximum_cycles=4)
    assert not result.history[-1].get('wall_initial_guess_extrapolated',False)
    np.testing.assert_array_equal(result.last_complete_state,result.trajectory[:10,-1])


class TransientWrapper(ContractingWrapper):
    def __init__(self, outcome):
        super().__init__();self.outcome=outcome

    def integrate_cycle(self,state,**kwargs):
        if len(self.inputs)==5 and self.outcome=='interrupt':
            raise IntegrationInterrupted('Deadline during transient validation')
        if len(self.inputs)==5 and self.outcome=='domain':
            raise MicrotubeDomainError('Invalid transient validation')
        angles,trajectory=super().integrate_cycle(state,**kwargs)
        if len(self.inputs)==5: trajectory[1,-1]=1.05
        if len(self.inputs)==6:
            trajectory[1,-1]=1.10 if self.outcome=='growth' else 1.049
        self.last=trajectory[:10,-1].copy()
        return angles,trajectory


@pytest.mark.parametrize('outcome',['contracts','growth','interrupt','domain','budget'])
def test_plausible_transient_needs_validation_before_retaining_it(outcome):
    wrapper=TransientWrapper(outcome);callbacks=[]
    result=solve(wrapper,maximum_cycles=5 if outcome=='budget' else 6,
        cycle_callback=lambda cycle,*_:callbacks.append(cycle))
    assert result.history[4]['wall_acceleration_pending_validation']
    if outcome=='contracts':
        assert result.history[5]['wall_acceleration_accepted']
        assert callbacks==[1,2,3,4,6]
        assert result.last_complete_state[1]==1.049
    else:
        assert callbacks==[1,2,3,4]
        assert result.last_complete_state[1]==1
        assert result.last_complete_state[8]==193.75
        if outcome in ('growth','domain'):
            assert any(r['phase']=='wall_acceleration_rollback' for r in result.solver_statistics)
    np.testing.assert_array_equal(result.last_complete_state,result.trajectory[:10,-1])


def test_default_waits_when_two_contracting_walls_have_unequal_reliability():
    temperatures=np.array([[300,350],[304,354],[306,356],[307,357.1]])
    assert adaptive_wall_proposal(temperatures,[1,1]) is None


def test_adaptive_campaign_option_is_explicit_and_identity_owned(tmp_path):
    from pathlib import Path
    from dada_solver.campaign.definition import CampaignDefinition
    root=Path(__file__).resolve().parents[1]
    source=root/'examples/motor_mechanics_stage7A5_edge.toml'
    original=CampaignDefinition(source)
    assert original.adaptive_wall_acceleration is None
    assert 'adaptive_wall_acceleration' not in original.numerical_settings
    configured=tmp_path/'adaptive.toml'
    configured.write_text(source.read_text()+'\n[adaptive_wall_acceleration]\ndamping = 0.97\n')
    adapted=CampaignDefinition(configured,base_path=original.base_path,
        hardware_path=root/'examples/motor_hardware_parallel_325c.toml')
    assert adapted.adaptive_wall_acceleration.damping==.97
    assert adapted.numerical_settings['adaptive_wall_acceleration']['damping']==.97
    assert adapted.numerical_settings!=original.numerical_settings
    assert adapted.definition_id!=original.definition_id
    assert adapted.wall_numerical_settings==original.wall_numerical_settings


def test_adaptive_campaign_rejects_reservoir_family(tmp_path):
    from pathlib import Path
    from dada_solver.campaign.definition import CampaignDefinition
    root=Path(__file__).resolve().parents[1]
    source=root/'examples/free_kinematics_campaign.toml'
    original=CampaignDefinition(source)
    assert original.wall_numerical_settings is None
    configured=tmp_path/'adaptive.toml'
    configured.write_text(source.read_text()+'\n[adaptive_wall_acceleration]\n')
    with pytest.raises(ValueError,match='wall evaluator'):
        CampaignDefinition(configured,base_path=original.base_path)
