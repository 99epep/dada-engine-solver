"""Compare saved fixed/adaptive replays without claiming exact equivalence.

Both input directories must use the same frozen manifest and tighter periodic
settings. The accuracy manifest supplies a separately stated comparison scale;
this offline summary does not alter convergence or rerun the solver.
"""
import argparse
import json
from pathlib import Path
import numpy as np


def audit(reference, candidate, accuracy_manifest):
    metadata = [json.loads((p/'environment.json').read_text()) for p in (reference, candidate)]
    if metadata[0]['manifest_sha256'] != metadata[1]['manifest_sha256']:
        raise ValueError('The two replays must use the same frozen manifest.')
    manifest = json.loads(accuracy_manifest.read_text())
    settings = manifest['numerical_settings']
    rtol, atol = settings['periodic_relative_tolerance'], settings['periodic_absolute_tolerance']
    cases = {}
    for path in sorted(candidate.glob('*_0.json')):
        b = json.loads(path.read_text())
        a = json.loads((reference/path.name).read_text())
        x, y = (np.asarray(row['final_state']) for row in (a, b))
        if any(row['result']['status'] != 'converged' for row in (a, b)):
            raise ValueError('Accuracy comparison requires two converged trajectories.')
        distance = float(np.max(np.abs(x-y)/(atol+rtol*np.maximum(np.abs(x), np.abs(y)))))
        relative = {key: abs(a['result'][key]-b['result'][key])/abs(a['result'][key])
                    for key in ('indicated_power_w', 'indicated_thermal_efficiency', 'total_mass_kg')}
        validity_equal = all(a['result'][section][key] == b['result'][section][key]
            for section, key in (('validity', 'verdict'), ('validity', 'failed_criteria'),
                                ('microtube_gas_domains', 'model_validity'),
                                ('microtube_gas_domains', 'failed_criteria')))
        runs = {}
        for label, row in (('fixed', a), ('adaptive', b)):
            runs[label] = dict(elapsed_seconds=row['elapsed_seconds'],
                cycles_completed=row['result']['convergence']['cycles_completed'],
                final_normalized_residual=row['result']['convergence']['last_normalized_periodic_error'],
                indicated_power_w=row['result']['indicated_power_w'],
                indicated_thermal_efficiency=row['result']['indicated_thermal_efficiency'],
                proposal_computation_seconds=sum(s.get('elapsed_seconds', 0) for s in row['solver_statistics']
                    if s['phase'] == 'wall_acceleration_decision'))
        cases[b['case']] = dict(runs=runs, comparison_scaled_state_distance=distance,
            relative_output_differences=relative, validity_classifications_equal=validity_equal,
            state_within_comparison_scale=distance <= 1,
            power_and_efficiency_within_relative_comparison_scale=all(
                relative[k] <= rtol for k in ('indicated_power_w', 'indicated_thermal_efficiency')))
    return dict(accuracy_manifest=str(accuracy_manifest), comparison_rtol=rtol,
        comparison_atol=atol, cases=cases,
        limitation='One timing per case; state periodic tolerances are not output-error guarantees.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=Path)
    parser.add_argument('candidate', type=Path)
    parser.add_argument('accuracy_manifest', type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.reference, args.candidate, args.accuracy_manifest), indent=2))
