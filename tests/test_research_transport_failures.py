"""A transport-domain rejection never escapes the evaluator into the runner."""
from types import SimpleNamespace
import numpy as np
import pytest
from dada_solver.exchangers.gas_transport import DiluteGasTransport
from dada_solver.exchangers.wall_cycle import WallCycleResult
from dada_solver.campaign.runner import OptimizationCampaign
from dada_solver.campaign.history import CampaignHistory
from dada_solver.research.presets import initialize_v3
from dada_solver.research.schema import load_study,compile_study
from dada_solver.research.study_io import dumps
from tests.test_research_v2 import activate


@pytest.mark.parametrize('temperature,category',[(99.9,'below'),(1000.1,'above')])
@pytest.mark.parametrize('stage',['integration','replay','diagnostics'])
def test_transport_rejection_then_next_candidate(tmp_path,monkeypatch,stage,temperature,category):
    path=initialize_v3(tmp_path/'study.toml')
    raw=load_study(path).data;activate(raw,'volume.swept_ratio',.01)
    path.write_text(dumps(raw));definition=compile_study(load_study(path))
    calls=[]
    def solve(wrapper,state,**kwargs):
        calls.append(state.copy())
        if len(calls)>1:
            return WallCycleResult('maximum_cycles','test stop',(),None,None,None)
        if stage=='integration': DiluteGasTransport().viscosity(temperature)
        values=np.r_[state,np.zeros(5)]
        return WallCycleResult('converged','test convergence',({'normalized_state_error':.5},),
            np.array([0.,2*np.pi]),np.column_stack([values,values]),state,{},{})
    monkeypatch.setattr('dada_solver.campaign.evaluator.solve_periodic_wall_motor',solve)
    def fail(*args,**kwargs): DiluteGasTransport().conductivity(temperature)
    if stage=='replay':
        monkeypatch.setattr('dada_solver.diagnostic_replay.replay_wall_trajectory',fail)
    if stage=='diagnostics':
        monkeypatch.setattr('dada_solver.diagnostic_replay.replay_wall_trajectory',lambda *a,**kw:SimpleNamespace())
        monkeypatch.setattr('dada_solver.campaign.evaluator.extract_cycle_diagnostics',fail)
    destination=tmp_path/'campaign'
    OptimizationCampaign(definition,destination).run(60,maximum_candidates=2)
    records=CampaignHistory(destination).load()
    assert len(calls)==2 and len(records)==2
    first=records[0]
    assert first['status']=='invalid_fluid_domain' and first['objective'] is None
    assert first['integrated'] and first['converged']==(stage!='integration')
    failure=first['diagnostics']['transport_failure']
    assert failure['category']==f'transport_temperature_{category}_domain'
    assert failure['temperature_k']==temperature and failure['minimum_temperature_k']==100
    assert failure['phase']==('integration' if stage=='integration' else 'postprocessing')
    assert records[1]['status']=='periodic_non_convergence'
