"""Exact report contracts and one reconstruction per final sample."""
from dataclasses import asdict, replace
import json
from pathlib import Path
import sys
import numpy as np
import pytest
from dada_solver.diagnostic_replay import replay_wall_trajectory
from dada_solver.exchangers.wall_cycle import WallDiagnosticCycle
from dada_solver.exchangers.gas_diagnostics import cycle_microtube_diagnostics
from dada_solver.results import extract_cycle_diagnostics
from dada_solver.validity import assess_cycle_validity
from dada_solver.dynamics import ThermodynamicModel

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'examples'))
from benchmark_solver_acceleration import DEFAULT_MANIFEST, legacy_design
from refine_motor_four_stage_hx9d_variable_gas import _load_basis,_build_design


@pytest.mark.parametrize('name',['production_warm','nearby_warm','production_cold','smooth_four_bar','legacy'])
def test_shared_facts_preserve_every_report_field(name,monkeypatch):
    manifest=json.loads(DEFAULT_MANIFEST.read_text())
    case=next(c for c in manifest['cases'] if c['name']==name)
    _,base=_load_basis()
    design=legacy_design() if name=='legacy' else _build_design(base,case['parameters'])
    if case.get('motion')=='base':design=replace(design,kinematics=base.kinematics)
    w=design.build()
    with np.load(ROOT/'outputs/solver_acceleration_stage4/numba_final'/(name+'_0.npz')) as f:
        ids=np.unique(np.r_[0,np.linspace(0,len(f['angles'])-1,51,dtype=int),len(f['angles'])-1])
        angles=f['angles'][ids];trajectory=f['trajectory'][:,ids]
    cycle=WallDiagnosticCycle(angles,trajectory)
    def heats(i,a,state):
        rates=w.thermal_rates(a,trajectory[:,i]);return tuple(r['gas_heat_w'] for r in rates)
    expected=extract_cycle_diagnostics(cycle,w.model,heats)
    validity=assess_cycle_validity(cycle,w.model,design.configuration.validity)
    domains=cycle_microtube_diagnostics(w,angles,trajectory)
    original=ThermodynamicModel.instantaneous_point;calls=[]
    def counted(self,*args,**kwargs):
        calls.append(1);return original(self,*args,**kwargs)
    monkeypatch.setattr(ThermodynamicModel,'instantaneous_point',counted)
    replay=replay_wall_trajectory(w,cycle,trajectory)
    assert asdict(extract_cycle_diagnostics(cycle,w.model,replay=replay))==asdict(expected)
    assert asdict(assess_cycle_validity(cycle,w.model,design.configuration.validity,replay=replay))==asdict(validity)
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
