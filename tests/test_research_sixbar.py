"""Physical parity is evaluated at fixed inputs, independently of search order."""
from dataclasses import asdict
import json
from pathlib import Path
from types import SimpleNamespace
import time

import numpy as np
import pytest

from dada_solver.campaign.evaluator import MachineEvaluator, EvaluationControl
from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
from dada_solver.research.schema import load_study, compile_study, candidate_for_values
from dada_solver.research.sixbar import PARAMETERS

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT/'src/dada_solver/research/data/sixbar-thermo5d.toml'


@pytest.fixture(scope='module')
def definition():
    return compile_study(load_study(TEMPLATE))


def values(definition):
    return {name:definition.study.basis.data['reference_parameters'][spec[2]] for name,spec in PARAMETERS.items()}


def test_stored_reference_physical_metrics(definition):
    if not definition.numerical_settings['wall_backend'].get('numba_available'):
        pytest.skip('Stored reference uses Numba; the independent legacy parity checks still exercise the available backend.')
    result=MachineEvaluator(definition).evaluate(candidate_for_values(definition,values(definition)))
    json.dumps(result, allow_nan=False)
    reference=definition.study.basis.data['reference_result']
    assert result['status']=='feasible'
    # Same-backend wrapper parity: the existing RHS equivalence tolerance is
    # 2e-11 relative / 1e-11 absolute. This is not a physical accuracy claim.
    for key in ('indicated_power_w','indicated_thermal_efficiency','heat_input_w',
                'maximum_pressure_pa','maximum_temperature_k','maximum_absolute_mass_flow_kg_s'):
        assert result['metrics'][key]==pytest.approx(reference[key],rel=2e-11,abs=1e-11)
    assert result['periodic_cycle_count']==reference['cycles_completed']
    assert result['periodic_convergence']['last_normalized_periodic_error']<=1
    assert result['metrics']['useful_mechanical_power_w'] is None
    assert result['warm_start_source_status']=='external_reference_initial_guess'
    assert result['derived']['external_air_capacity_diagnostics']['heat_in']['external_to_peak_internal_capacity_rate_ratio'] > 1
    assert len(result['derived']['local_reflux']['minimum_signed_flows_kg_s']) == 4


@pytest.mark.parametrize('variant', [1,2])
def test_varied_inputs_match_original_evaluator(definition,variant,monkeypatch,tmp_path):
    if not (ROOT/'outputs/motor_spline_from_linear_260k/report.json').exists():
        pytest.skip('Optional historical parity oracle requires the source-checkout output artifacts.')
    pytest.importorskip('matplotlib')
    monkeypatch.chdir(ROOT)
    monkeypatch.setenv('MPLCONFIGDIR',str(tmp_path/'mpl'))
    monkeypatch.syspath_prepend(str(ROOT/'examples'))
    import optimize_sixbar_pairs_thermo5d_3952_hlat25 as legacy
    seed,old_definition,mass,geometry,_,_=legacy.make_3952_basis()
    pair=legacy.load_pair(1)
    p=values(definition)
    p['volume.swept_ratio'] *= .98 if variant==1 else 1.02
    p['microtube.heat_in.tube_count'] += 17*variant
    p['microtube.heat_out.tube_count'] -= 23*variant
    p['microtube.heat_in.tube_length_m'] *= 1.01
    p['microtube.heat_out.tube_length_m'] *= .99
    old_p={PARAMETERS[k][2]:v for k,v in p.items()}
    old_design=legacy.build_design(seed,geometry,mass,old_p,pair['small'],pair['large'])
    new_design=definition.adapter.build(p)
    assert asdict(new_design.configuration)==asdict(old_design.configuration)
    assert asdict(new_design.heat_in)==asdict(old_design.heat_in)
    assert asdict(new_design.heat_out)==asdict(old_design.heat_out)
    for side in ('small','large'):
        for suffix in ('cylinder_volume','cylinder_volume_derivative'):
            method=f'{side}_{suffix}'
            # Compare the injected study-angle laws before the single motor reversal.
            a=[getattr(new_design.kinematics,method)(t) for t in np.linspace(0,2*np.pi,361)]
            b=[getattr(old_design.kinematics,method)(t) for t in np.linspace(0,2*np.pi,361)]
            np.testing.assert_allclose(a,b,rtol=2e-13,atol=2e-13)
    initial=np.array(definition.initial_wall_state['values'])
    hw=legacy.hardware(old_design)
    caps=np.array([hw[side]['wall_capacity_j_k'] for side in ('H_i','H_o')])
    initial[8:10] *= caps/np.array(definition.initial_wall_state['wall_capacities_j_k'])
    expected=legacy.evaluate(label='research_parity',design=old_design,definition=old_definition,
                             initial=initial,deadline=time.monotonic()+180,candidate_seconds=180)
    actual=MachineEvaluator(definition).evaluate(candidate_for_values(definition,p))
    json.dumps(actual, allow_nan=False)
    assert actual['converged']==(expected['compact']['status']=='converged')
    assert (actual['status']=='feasible')==expected['feasible']
    for key in ('indicated_power_w','indicated_thermal_efficiency','heat_input_w',
                'maximum_pressure_pa','maximum_temperature_k','maximum_absolute_mass_flow_kg_s'):
        assert actual['metrics'][key]==pytest.approx(expected['compact'][key],rel=2e-11,abs=1e-11)
    failures=[c['name'] for c in actual['constraints'] if not c['available'] or not c['satisfied']]
    assert failures==expected['reasons']


def test_safe_domain_retry_preserves_wall_state_and_deadline(definition,monkeypatch):
    calls=[]
    def solve(wrapper,state,**kwargs):
        calls.append((state.copy(),kwargs['progress_callback'].__self__.deadline))
        if len(calls)==1: raise MicrotubeDomainError('synthetic initial-state failure')
        return SimpleNamespace(status='interrupted',message='synthetic deadline',history=(),
                               last_complete_state=None,converged=False,backend_statistics={})
    monkeypatch.setattr('dada_solver.campaign.evaluator.solve_periodic_wall_motor',solve)
    result=MachineEvaluator(definition).evaluate_with_control(candidate_for_values(definition,values(definition)),
        EvaluationControl(deadline=50,clock=lambda:0))
    assert result['status']=='budget_exhausted' and result['safe_retry_used']
    assert len(calls)==2 and calls[0][1]==calls[1][1]==50
    np.testing.assert_array_equal(calls[0][0][8:],calls[1][0][8:])
    assert not np.array_equal(calls[0][0][:8],calls[1][0][:8])


def test_nonconvergence_has_no_available_objective_or_constraints(definition,monkeypatch):
    monkeypatch.setattr('dada_solver.campaign.evaluator.solve_periodic_wall_motor',lambda *a,**k:
        SimpleNamespace(status='maximum_cycles',message='not periodic',history=(),
                        last_complete_state=None,converged=False,backend_statistics={}))
    result=MachineEvaluator(definition).evaluate(candidate_for_values(definition,values(definition)))
    assert result['status']=='periodic_non_convergence' and result['objective'] is None
    assert all(not c['available'] for c in result['constraints'])


def test_unexpected_programming_error_is_not_physical_rejection(definition,monkeypatch):
    def broken(*a,**k): raise TypeError('programming error')
    monkeypatch.setattr('dada_solver.campaign.evaluator.solve_periodic_wall_motor',broken)
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
