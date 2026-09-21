"""Require identical public reports, final states and complete trajectories."""
import argparse
import json
from pathlib import Path
import numpy as np


def compare(reference, candidate):
    rows=[]
    for path in sorted(candidate.glob('*_*.json')):
        actual=json.loads(path.read_text())
        if 'case' not in actual: continue
        expected=json.loads((reference/(actual['case']+'_0.json')).read_text())
        row=dict(case=actual['case'],repeat=actual['repeat'],
            report_exact=expected['result']==actual['result'],
            final_state_exact=expected['final_state']==actual['final_state'],
            exception_exact=expected['exception']==actual['exception'])
        a=reference/(actual['case']+'_0.npz');b=path.with_suffix('.npz')
        row['trajectory_presence_exact']=a.exists()==b.exists()
        if a.exists() and b.exists():
            with np.load(a) as x, np.load(b) as y:
                row['angles_exact']=np.array_equal(x['angles'],y['angles'])
                row['all_15_rows_exact']=np.array_equal(x['trajectory'],y['trajectory'])
        rows.append(row)
    passed=bool(rows) and all(all(v for k,v in r.items() if k not in ('case','repeat')) for r in rows)
    return dict(passed=passed,runs=rows)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('reference',type=Path);p.add_argument('candidate',type=Path)
    a=p.parse_args();result=compare(a.reference,a.candidate)
    print(json.dumps(result,indent=2))
    if not result['passed']:raise SystemExit(1)
