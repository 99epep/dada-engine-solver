"""Freeze bounded historical-machine replays under the new dilute transport law.

No search. Historical reference files remain untouched. This tool is deliberately
separate from the CoolProp property oracle: thermodynamic results use DADA's EOS.
"""
import hashlib
import json
from pathlib import Path
import tempfile
from dada_solver.research.presets import initialize_v2
from dada_solver.research.schema import load_study,compile_study,candidate_for_values
from dada_solver.research.study_io import dumps
from dada_solver.campaign.evaluator import MachineEvaluator

ROOT=Path(__file__).resolve().parents[1]
METRICS=('indicated_power_w','indicated_thermal_efficiency','heat_input_w',
         'maximum_pressure_pa','maximum_temperature_k','maximum_absolute_mass_flow_kg_s')


def evaluate_case(family,directory):
    path=initialize_v2(directory/'study.toml',family,family,champion=family=='structured_c2_15p')
    raw=load_study(path).data
    if family in ('fourier_c2','free_spline'):
        case=json.loads((ROOT/'tests/fixtures/research_v2/thermodynamic.json').read_text())[family]
        basis=path.with_suffix('.basis.json');basis.write_text(json.dumps(case['basis']))
        raw['sources']['machine']['sha256']=hashlib.sha256(basis.read_bytes()).hexdigest()
        raw['parameters']=[r for r in raw['parameters'] if r['name'].startswith('kinematics.')]
        raw['constraints'].append(dict(type='minimum_motor_power',required_power=25.,unit='W'))
        raw['policies']['outlet_valve_cda']='fixed_source_cda'
        raw['warm_start']['initial_source']='source_exact';raw['numerical']['backend']=case['backend']
        old=case['metrics']
    else:
        basis=load_study(path).basis.data
        old=(json.loads((ROOT/'src/dada_solver/research/data/rank01_basis.json').read_text())['reference_result']
             if family=='six_bar' else basis['provenance']['historical_result'])
    path.write_text(dumps(raw));definition=compile_study(load_study(path))
    result=MachineEvaluator(definition).evaluate(candidate_for_values(definition,{}))
    if not result['converged']: raise RuntimeError(result['reason'])
    metrics={k:result['metrics'][k] for k in METRICS}
    return dict(metrics=metrics,cycles=result['periodic_cycle_count'],status=result['status'],
                relative_change_from_historical={k:metrics[k]/old[k]-1 for k in METRICS})


if __name__=='__main__':
    cases={}
    for family in ('fourier_c2','free_spline','six_bar','structured_c2_15p'):
        with tempfile.TemporaryDirectory() as directory:
            cases[family]=evaluate_case(family,Path(directory))
        print(family,cases[family]['status'],flush=True)
    output=ROOT/'tests/data/transport_v2_thermal_reference.json'
    from dada_solver import numerical_primitives
    output.write_text(json.dumps(dict(correlation_version='dilute_species_v2',
        scope='Same historical geometry and numerical settings, new transport; no optimization.',
        numerical_primitives_sha256=hashlib.sha256(Path(numerical_primitives.__file__).read_bytes()).hexdigest(),
        cases=cases),indent=2)+'\n')
