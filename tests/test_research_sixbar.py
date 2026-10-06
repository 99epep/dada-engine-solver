"""Physical parity is evaluated at fixed inputs, independently of search order."""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from dada_solver.campaign.evaluator import MachineEvaluator, EvaluationControl
from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
from dada_solver.research.schema import load_study, compile_study, candidate_for_values
from dada_solver.research.presets import initialize_kinematics

ROOT = Path(__file__).resolve().parents[1]
def template(directory):
    path = initialize_kinematics(directory/'study.toml', 'six_bar', 'six_bar')
    raw = load_study(path).data
    for row in raw['parameters']:
        if row['name'] == 'operation.frequency_hz':
            value = row.pop('value')
            row.update(initial=value, lower=value*.99, upper=value*1.01, kind='continuous', transform='linear')
    from dada_solver.research.study_io import dumps
    path.write_text(dumps(raw))
    return path


@pytest.fixture(scope='module')
def definition(tmp_path_factory):
    return compile_study(load_study(template(tmp_path_factory.mktemp('sixbar'))))


def values(definition):
    return {p.name:p.initial for p in definition.space.parameters}




def test_safe_domain_retry_preserves_wall_state_and_deadline(definition,monkeypatch):
    calls=[]
    def solve(wrapper,state,**kwargs):
        calls.append((state.copy(),kwargs['progress_callback'].__self__.deadline))
        if len(calls)==1: raise MicrotubeDomainError('synthetic initial-state failure')
        return SimpleNamespace(status='interrupted',message='synthetic deadline',history=(),
                               last_complete_state=None,converged=False,backend_statistics={})
    monkeypatch.setattr('dada_solver.campaign.evaluator.solve_periodic_wall_machine',solve)
    result=MachineEvaluator(definition).evaluate_with_control(candidate_for_values(definition,values(definition)),
        EvaluationControl(deadline=50,clock=lambda:0))
    assert result['status']=='budget_exhausted' and result['safe_retry_used']
    assert len(calls)==2 and calls[0][1]==calls[1][1]==50
    np.testing.assert_array_equal(calls[0][0][8:],calls[1][0][8:])
    assert not np.array_equal(calls[0][0][:8],calls[1][0][:8])


def test_nonconvergence_has_no_objective_or_available_physical_constraints(definition,monkeypatch):
    monkeypatch.setattr('dada_solver.campaign.evaluator.solve_periodic_wall_machine',lambda *a,**k:
        SimpleNamespace(status='maximum_cycles',message='not periodic',history=(),
                        last_complete_state=None,converged=False,backend_statistics={}))
    result=MachineEvaluator(definition).evaluate(candidate_for_values(definition,values(definition)))
    assert result['status']=='periodic_non_convergence' and result['objective'] is None
    assert all(not c['available'] for c in result['constraints'] if c['name'] == 'valid_thermodynamic_model')
    assert all(c['available'] for c in result['constraints'] if c['name'].startswith(('small.', 'large.')))


def test_unexpected_programming_error_is_not_physical_rejection(definition,monkeypatch):
    def broken(*a,**k): raise TypeError('programming error')
    monkeypatch.setattr('dada_solver.campaign.evaluator.solve_periodic_wall_machine',broken)
    with pytest.raises(TypeError,match='programming error'):
        MachineEvaluator(definition).evaluate(candidate_for_values(definition,values(definition)))


def test_real_sobol_candidate_is_durably_serializable(definition,tmp_path):
    from dada_solver.campaign.runner import OptimizationCampaign
    campaign = OptimizationCampaign(definition,tmp_path/'real_sobol')
    summary = campaign.run(180,maximum_candidates=1)
    assert summary['attempted'] == 1
    records = campaign.history.load()
    assert len(records) == 1
    assert records[0]['sequence_index'] == 0
    assert records[0]['status'] in ('feasible','converged_infeasible')
    assert campaign.history.state()['pending'] is None
    json.dumps(records,allow_nan=False)
