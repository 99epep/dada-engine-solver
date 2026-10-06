"""Exact report contracts and one reconstruction per final sample."""
from dataclasses import asdict, replace
import numpy as np
import pytest
from dada_solver.diagnostic_replay import replay_wall_trajectory
from dada_solver.exchangers.wall_cycle import WallDiagnosticCycle
from dada_solver.exchangers.gas_diagnostics import cycle_microtube_diagnostics
from dada_solver.results import extract_cycle_diagnostics
from dada_solver.validity import assess_cycle_validity
from dada_solver.dynamics import ThermodynamicModel

from tests.test_solver_acceleration import variable_wrapper
from dada_solver.configuration import ValidityThresholds
from dada_solver.state import UniformCharge


@pytest.mark.parametrize('variable_film', [False, True])
def test_shared_facts_preserve_every_report_field(variable_film,monkeypatch):
    w=variable_wrapper()
    if not variable_film:
        w=replace(w,heat_in=replace(w.heat_in,gas_film=None),
                  heat_out=replace(w.heat_out,gas_film=None))
    angles=np.linspace(0.,2*np.pi,17)
    columns=[]
    for angle in angles:
        gas=UniformCharge(2e5,350.).create_state(w.model.gas,w.model.volumes(angle))
        columns.append(np.r_[gas.as_array(),w.heat_in.wall_capacity_j_k*400.,
                              w.heat_out.wall_capacity_j_k*330.,np.zeros(5)])
    trajectory=np.array(columns).T
    thresholds=ValidityThresholds(maximum_cp_variation=.02,maximum_compressibility_deviation=.01)
    cycle=WallDiagnosticCycle(angles,trajectory)
    def heats(i,a,state):
        rates=w.thermal_rates(a,trajectory[:,i]);return tuple(r['gas_heat_w'] for r in rates)
    expected=extract_cycle_diagnostics(cycle,w.model,heats)
    validity=assess_cycle_validity(cycle,w.model,thresholds)
    domains=cycle_microtube_diagnostics(w,angles,trajectory)
    original=ThermodynamicModel.instantaneous_point;calls=[]
    def counted(self,*args,**kwargs):
        calls.append(1);return original(self,*args,**kwargs)
    monkeypatch.setattr(ThermodynamicModel,'instantaneous_point',counted)
    replay=replay_wall_trajectory(w,cycle,trajectory)
    assert asdict(extract_cycle_diagnostics(cycle,w.model,replay=replay))==asdict(expected)
    assert asdict(assess_cycle_validity(cycle,w.model,thresholds,replay=replay))==asdict(validity)
    assert cycle_microtube_diagnostics(w,angles,trajectory,replay=replay)==domains
    override=extract_cycle_diagnostics(cycle,w.model,lambda *_:(123.,-45.),replay=replay)
    assert override.cold_heat_rate_extrema.minimum==override.cold_heat_rate_extrema.maximum==123.
    assert override.hot_heat_rate_extrema.minimum==override.hot_heat_rate_extrema.maximum==-45.
    assert len(calls)==len(angles)
    assert not replay.pressures.flags.writeable
    changed=trajectory.copy();changed[0,0]*=1.001
    with pytest.raises(ValueError,match='Cycle states do not match'):
        replay_wall_trajectory(w,cycle,changed)
    with pytest.raises(ValueError,match='different model or trajectory'):
        replay.require(wrapper=replace(w))
