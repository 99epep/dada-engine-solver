"""Anderson remains an initial-guess experiment with physical-cycle certification."""
from dataclasses import replace
import numpy as np
import pytest
from dada_solver.exchangers.wall_iteration import (
    MassConservingCoordinates, PeriodicMapPair, AndersonAccelerationSettings,
    anderson_proposal, MASS_INDICES, AdaptiveWallAccelerationSettings)
from dada_solver.exchangers.wall_cycle import solve_periodic_wall_motor,WallCycleNumericalSettings
from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
from dada_solver.integration import IntegrationInterrupted

ANCHOR=np.array([.002,500.,.03,8000.,.0002,50.,.0001,30.,10000.,20000.])


def test_coordinates_round_trip_mass_and_owned_history():
    c=MassConservingCoordinates.from_anchor(ANCHOR);z=np.linspace(-.01,.01,9)
    np.testing.assert_allclose(c.basis.T@c.basis,np.eye(9),atol=1e-15)
    np.testing.assert_allclose((c.scales[:,None]*c.basis)[MASS_INDICES].sum(axis=0),0,atol=1e-17)
    state=c.reconstruct(z);np.testing.assert_allclose(c.project(state),z,atol=2e-16)
    assert state[MASS_INDICES].sum()==pytest.approx(c.total_mass,rel=1e-15)
    p=PeriodicMapPair(state,ANCHOR,1.);state[:]=99
    assert not p.initial_state.flags.writeable and p.initial_state[0]!=99


@pytest.mark.parametrize('kw',[dict(memory=1),dict(memory=True),dict(memory=11),dict(memory=3.5),
    dict(minimum_history=1),dict(minimum_history=4),dict(damping=0),dict(damping=1.01),
    dict(damping=True),dict(maximum_coefficient_l1=.9),dict(maximum_residual_growth=.5),
    dict(maximum_relative_state_correction=float('nan')),dict(maximum_relative_state_correction=0)])
def test_invalid_controls(kw):
    with pytest.raises(ValueError):AndersonAccelerationSettings(**kw)


class LinearWrapper:
    def __init__(self,failure=None):
        self.c=MassConservingCoordinates.from_anchor(ANCHOR)
        self.diag=np.array([.95,.93,.025,.01,.008,.005,.002,.001,0.])
        self.inputs=[];self.ends=[];self.failure=failure;self.failed=False
    def integrate_cycle(self,state,**kwargs):
        state=state.copy();self.inputs.append(state)
        end=self.c.reconstruct(self.diag*self.c.project(state))
        if len(self.inputs)==3 and self.failure and not self.failed:
            self.failed=True
            if self.failure=='domain':raise MicrotubeDomainError('Test accelerated domain failure.')
            if self.failure=='interrupt':raise IntegrationInterrupted('Test accelerated deadline.')
            if self.failure=='growth':end=state.copy();end[8]*=2
            if self.failure=='nonfinite':end[8]=np.nan
        self.ends.append(end.copy())
        trajectory=np.zeros((15,2));trajectory[:10]=np.array([state,end]).T
        return np.array([0.,2*np.pi]),trajectory


def solve(w,acceleration=None,budget=300,callback=None):
    return solve_periodic_wall_motor(w,w.c.reconstruct(np.linspace(.1,-.1,9)),maximum_cycles=budget,
        settings=WallCycleNumericalSettings(accelerate_walls=False),anderson_acceleration=acceleration,
        cycle_callback=callback)


def test_linear_map_reduces_calls_and_certifies_with_ordinary_cycle():
    plain=LinearWrapper();a=solve(plain)
    fast=LinearWrapper();b=solve(fast,AndersonAccelerationSettings(memory=4,damping=1))
    assert a.converged and b.converged
    assert len(fast.inputs)<len(plain.inputs)
    assert b.history[-1]['normalized_state_error']<=1
    np.testing.assert_array_equal(b.last_complete_state,b.trajectory[:10,-1])
    for state in fast.inputs:assert state[MASS_INDICES].sum()==pytest.approx(ANCHOR[MASS_INDICES].sum(),rel=1e-14)


def test_duplicate_and_collinear_history_reduces_rank_or_skips():
    c=MassConservingCoordinates.from_anchor(ANCHOR);d=np.zeros(9);d[-1]=.01
    x=c.reconstruct(d);y=c.reconstruct(.5*d)
    duplicate=[PeriodicMapPair(x,y,10.)]*4
    assert anderson_proposal(duplicate,settings=AndersonAccelerationSettings(memory=4)).state is None
    pairs=[PeriodicMapPair(c.reconstruct(a*d),c.reconstruct(a*.5*d),10.) for a in (8,4,2,1)]
    p=anderson_proposal(pairs,settings=AndersonAccelerationSettings(memory=4))
    assert p.state is not None and p.statistics['effective_memory']==2
    assert np.all(np.isfinite(p.state))


def aggressive_pairs():
    # Two positive map pairs extrapolate one wall energy below zero.
    x=ANCHOR.copy();x[9]=20000;y=x.copy();y[9]=14000
    u=x.copy();u[9]=14000;v=x.copy();v[9]=9500
    return [PeriodicMapPair(x,y,1.),PeriodicMapPair(u,v,1.)]


def test_uniform_correction_bound_and_positivity_rejection():
    pairs=aggressive_pairs();settings=AndersonAccelerationSettings(damping=1,maximum_relative_state_correction=.1)
    p=anderson_proposal(pairs,settings=settings)
    assert p.state is not None
    assert p.statistics['raw_relative_correction']>.1
    assert p.statistics['applied_relative_correction']==pytest.approx(.1)
    c=MassConservingCoordinates.from_anchor(pairs[-1].final_state)
    unlimited=anderson_proposal(pairs,settings=replace(settings,maximum_relative_state_correction=10.))
    assert unlimited.state is None and unlimited.statistics['status']=='rejected_preflight'
    assert p.state[MASS_INDICES].sum()==pytest.approx(c.total_mass,rel=1e-15)
    alpha=np.array(p.statistics['coefficients'])
    raw=np.column_stack([c.project(pair.final_state) for pair in pairs])@alpha
    np.testing.assert_allclose(c.project(p.state),raw*p.statistics['correction_scale'],atol=1e-15)


@pytest.mark.parametrize('failure',['growth','domain','nonfinite'])
def test_rejected_trial_restores_endpoint_and_counts_attempt(failure):
    w=LinearWrapper(failure);callbacks=[]
    r=solve(w,AndersonAccelerationSettings(memory=3),budget=4,callback=lambda cycle,*_:callbacks.append(cycle))
    assert r.anderson_statistics['proposal_count']>=1
    assert r.anderson_statistics['rollback_count']==1
    np.testing.assert_array_equal(w.inputs[3],w.ends[1])
    assert callbacks==[1,2,4]
    assert r.anderson_statistics['attempted_cycles']==4
    assert r.anderson_statistics['retained_cycles']==3
    assert r.anderson_statistics['domain_failure_count']==(failure=='domain')
    decisions=[s for s in r.solver_statistics if s['phase']=='anderson_decision']
    assert [s['cycle'] for s in decisions if s['status']=='proposed']==[2]


def test_interrupt_returns_physical_endpoint_not_proposal():
    w=LinearWrapper('interrupt');r=solve(w,AndersonAccelerationSettings(),budget=4)
    assert r.status=='interrupted'
    np.testing.assert_array_equal(r.last_complete_state,w.ends[-1])
    assert not np.array_equal(r.last_complete_state,w.inputs[-1])
    assert r.anderson_statistics['attempted_cycles']==3


def test_last_budget_and_mutually_exclusive_modes():
    w=LinearWrapper();r=solve(w,AndersonAccelerationSettings(),budget=2)
    assert r.anderson_statistics['proposal_count']==0
    np.testing.assert_array_equal(r.last_complete_state,w.ends[-1])
    with pytest.raises(ValueError,match='mutually exclusive'):
        solve_periodic_wall_motor(w,ANCHOR,maximum_cycles=2,
            adaptive_acceleration=AdaptiveWallAccelerationSettings(),anderson_acceleration=AndersonAccelerationSettings())


def test_ordinary_domain_failure_propagates():
    class Broken(LinearWrapper):
        def integrate_cycle(self,*args,**kwargs):raise MicrotubeDomainError('Ordinary failure.')
    with pytest.raises(MicrotubeDomainError,match='Ordinary failure'):solve(Broken(),AndersonAccelerationSettings())


def test_coefficient_cap_and_incompatible_historical_mass_skip():
    pairs=aggressive_pairs()
    p=anderson_proposal(pairs,settings=AndersonAccelerationSettings(maximum_coefficient_l1=2))
    assert p.state is None and p.statistics['status']=='skipped_coefficients'
    altered=pairs[0].initial_state.copy();altered[0]*=1.01
    p=anderson_proposal([PeriodicMapPair(altered,pairs[0].final_state,1.),pairs[1]])
    assert p.state is None and p.statistics['status']=='rejected_preflight'


def test_slow_scalar_mode_requires_coefficients_above_experimental_cap():
    x=ANCHOR.copy();x[9]=40000
    y=x.copy();y[9]=39000
    end=y.copy();end[9]=38050
    p=anderson_proposal([PeriodicMapPair(x,y,10.),PeriodicMapPair(y,end,9.5)])
    assert p.state is None and p.statistics['status']=='skipped_coefficients'
    assert p.statistics['coefficient_l1']==pytest.approx(39.)


def test_fixed_wall_extrapolation_is_suppressed():
    a=LinearWrapper();b=LinearWrapper();initial=a.c.reconstruct(np.linspace(.1,-.1,9))
    settings=AndersonAccelerationSettings(memory=4,damping=1)
    x=solve(a,settings)
    y=solve_periodic_wall_motor(b,initial,maximum_cycles=300,
        settings=WallCycleNumericalSettings(accelerate_walls=True),anderson_acceleration=settings)
    assert x.history==y.history
    np.testing.assert_array_equal(x.last_complete_state,y.last_complete_state)


def test_moderate_non_normal_growth_is_retained():
    class Moderate(LinearWrapper):
        def integrate_cycle(self,state,**kwargs):
            angles,trajectory=super().integrate_cycle(state,**kwargs)
            if len(self.inputs)==3:
                previous=np.max(abs(self.ends[1]-self.inputs[1])/(1e-12+1e-6*np.maximum(abs(self.ends[1]),abs(self.inputs[1]))))
                target=1.5*previous
                trajectory[8,-1]=(state[8]+target*1e-12)/(1-target*1e-6)
                self.ends[-1]=trajectory[:10,-1].copy()
            return angles,trajectory
    w=Moderate();r=solve(w,AndersonAccelerationSettings(),budget=3)
    assert r.anderson_statistics['accepted_trial_count']==1
    assert r.anderson_statistics['rollback_count']==0
    assert r.history[-1]['normalized_state_error']>r.history[-2]['normalized_state_error']
