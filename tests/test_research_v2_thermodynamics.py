"""Versioned current-physics thermal regression at fixed inputs."""
import json
from pathlib import Path
import pytest
from dada_solver.research.presets import initialize_kinematics
from dada_solver.research.schema import load_study,compile_study,candidate_for_values
from dada_solver.campaign.evaluator import MachineEvaluator

METRICS=('indicated_power_w','indicated_thermal_efficiency','heat_input_w','maximum_pressure_pa','maximum_temperature_k','maximum_absolute_mass_flow_kg_s')


@pytest.mark.parametrize('family',['six_bar','structured_c2_15p'])
def test_packaged_current_physics_reference(tmp_path,family):
    path=initialize_kinematics(tmp_path/'study.toml',family,family,champion=family=='structured_c2_15p')
    definition=compile_study(load_study(path))
    if not definition.numerical_settings['wall_backend'].get('numba_available'):pytest.skip('Stored reference uses Numba.')
    result=MachineEvaluator(definition).evaluate(candidate_for_values(definition,{}))
    case=json.loads((Path(__file__).parent/'data/developing_entry_thermal_reference.json').read_text())['cases'][family]
    ref=dict(case['metrics'],cycles_completed=case['cycles'])
    assert result['status']=='feasible',result['reason']
    for key in METRICS:assert result['metrics'][key]==pytest.approx(ref[key],rel=2e-11,abs=1e-11)
    assert result['periodic_cycle_count']==ref['cycles_completed']
