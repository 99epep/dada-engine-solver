"""Stored historical physical results, independently of historical search order."""
import hashlib
import json
from pathlib import Path
import tomllib
import pytest
from dada_solver.research.presets import initialize_v2
from dada_solver.research.schema import load_study,compile_study,candidate_for_values
from dada_solver.research.study_io import dumps
from dada_solver.campaign.evaluator import MachineEvaluator

FIXTURES=Path(__file__).parent/'fixtures/research_v2'
METRICS=('indicated_power_w','indicated_thermal_efficiency','heat_input_w','maximum_pressure_pa','maximum_temperature_k','maximum_absolute_mass_flow_kg_s')


@pytest.mark.parametrize('family',['slider_crank','four_bar','fourier_c2','free_spline'])
def test_historical_thermal_machine_and_initial_state(tmp_path,family):
    case=json.loads((FIXTURES/'thermodynamic.json').read_text())[family]
    path=initialize_v2(tmp_path/'study.toml',family,family)
    raw=tomllib.loads(path.read_text());basis=path.with_suffix('.basis.json')
    basis.write_text(json.dumps(case['basis']))
    raw['sources']['machine']['sha256']=hashlib.sha256(basis.read_bytes()).hexdigest()
    raw['parameters']=[r for r in raw['parameters'] if r['name'].startswith('kinematics.')]
    raw['policies']['outlet_valve_cda']='fixed_source_cda'
    raw['warm_start']['initial_source']='source_exact'
    raw['numerical']['backend']=case['backend']
    path.write_text(dumps(raw));definition=compile_study(load_study(path))
    if case['backend']=='numba' and not definition.numerical_settings['wall_backend'].get('numba_available'):
        pytest.skip('Stored historical case uses the optional Numba backend.')
    result=MachineEvaluator(definition).evaluate(candidate_for_values(definition,{}))
    expected_status='feasible' if case['metrics']['indicated_power_w']>=25 else 'converged_infeasible'
    assert result['status']==expected_status,result['reason']
    if expected_status=='converged_infeasible':assert result['reason']=='minimum_motor_power: violated'
    tolerance=1e-5 if family in ('slider_crank','four_bar') else 5e-8
    for key in METRICS:
        # Roundoff in volume partition/common-frame normalization perturbs LSODA
        # steps and sampled extrema. K2 measured maxima are 2.45e-6 / 5.01e-6;
        # use 1e-5 for these historical integrated results, not the 2e-11 RHS
        # tolerance. Modern fixed-inventory references are substantially closer.
        assert result['metrics'][key]==pytest.approx(case['metrics'][key],rel=tolerance,abs=1e-10)
    assert result['periodic_cycle_count']==case['cycles']
    assert result['periodic_convergence']['last_normalized_periodic_error']<=1
    assert result['metrics']['useful_mechanical_power_w'] is None
    assert result['metrics']['cooling_power_w'] is None
    assert len(result['derived']['kinematic_metrics'])==2
    assert all('value' in c and 'limit' in c and 'relative_margin' in c for c in result['constraints'])


@pytest.mark.parametrize('family',['six_bar','structured_c2_15p'])
def test_packaged_historical_reference(tmp_path,family):
    path=initialize_v2(tmp_path/'study.toml',family,family,champion=family=='structured_c2_15p')
    definition=compile_study(load_study(path))
    if not definition.numerical_settings['wall_backend'].get('numba_available'):pytest.skip('Stored reference uses Numba.')
    result=MachineEvaluator(definition).evaluate(candidate_for_values(definition,{}))
    if family=='six_bar':
        ref=json.loads((Path(__file__).parents[1]/'src/dada_solver/research/data/rank01_basis.json').read_text())['reference_result']
    else:ref=definition.study.basis.data['provenance']['historical_result']
    assert result['status']=='feasible',result['reason']
    for key in METRICS:assert result['metrics'][key]==pytest.approx(ref[key],rel=2e-11,abs=1e-11)
    assert result['periodic_cycle_count']==ref['cycles_completed']
