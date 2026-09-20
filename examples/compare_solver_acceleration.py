"""Compare frozen benchmark physics exactly; report timing separately."""
import argparse
import json
from pathlib import Path
import numpy as np


def compare(reference, candidate):
    checks=[]; interruptions=[]
    for path in sorted(candidate.glob('*_*.json')):
        if path.name in ('environment.json','summary.json'): continue
        row=json.loads(path.read_text())
        if 'case' not in row: continue
        original=reference/f"{row['case']}_0.json"
        baseline=json.loads(original.read_text())
        if row['case']=='interrupted':
            # Wall-clock callback locations change with throughput. Check the
            # status and retained physical state, not identical cycle counts.
            assert row['result']['status']==baseline['result']['status']=='interrupted'
            assert row['exception']==baseline['exception']
            if row['final_state'] is None:
                assert row['result']['convergence']['cycles_completed']==0
            else:
                with np.load(path.with_suffix('.npz')) as saved:
                    np.testing.assert_array_equal(row['final_state'],saved['trajectory'][:10,-1])
            interruptions.append(dict(file=path.name,
                retained_cycles=row['result']['convergence']['cycles_completed']))
            continue
        for field in ('result','final_state','exception'):
            if row[field]!=baseline[field]:
                raise AssertionError(f'{path.name}: changed {field}')
        trajectory=path.with_suffix('.npz')
        if trajectory.exists():
            with np.load(trajectory) as actual, np.load(original.with_suffix('.npz')) as expected:
                for field in ('angles','trajectory'):
                    np.testing.assert_array_equal(actual[field],expected[field],err_msg=str(path))
        checks.append(path.name)
    if not checks: raise ValueError('No benchmark cases found.')
    old=json.loads((reference/'summary.json').read_text())
    new=json.loads((candidate/'summary.json').read_text())
    ratios={key:old[key]['median_seconds']/value['median_seconds'] for key,value in new.items() if key not in ('interrupted','invalid_domain')}
    result=dict(exact_comparisons=checks,interruption_semantics_checks=interruptions,median_speedup=ratios,
                note='Observed maximum is reported separately; three repeats do not estimate population tail latency.')
    print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference',type=Path);parser.add_argument('candidate',type=Path)
    args=parser.parse_args();compare(args.reference,args.candidate)
