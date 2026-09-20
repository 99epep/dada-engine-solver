"""Report physical differences and timings; never label compiled replays exact."""
import argparse
import json
from pathlib import Path
import numpy as np


def compare(reference,candidate,manifest):
    identities=[json.loads((p/'environment.json').read_text()) for p in (reference,candidate)]
    if identities[0]['manifest_sha256']!=identities[1]['manifest_sha256']:
        raise ValueError('Both replays must use the same frozen numerical experiment.')
    settings=json.loads(manifest.read_text())['numerical_settings']
    rtol=settings['periodic_relative_tolerance'];atol=settings['periodic_absolute_tolerance']
    results=[]
    for path in sorted(candidate.glob('*_*.json')):
        b=json.loads(path.read_text())
        if 'case' not in b: continue
        a=json.loads((reference/(b['case']+'_0.json')).read_text())
        x,y=a['result'],b['result']
        row=dict(case=b['case'],repeat=b['repeat'],status_equal=x['status']==y['status'],
            exception_equal=a['exception']==b['exception'],elapsed_seconds=b['elapsed_seconds'])
        if b['case']=='invalid_domain': row['message_equal']=x['message']==y['message']
        if b['final_state'] is not None:
            u,v=np.asarray(a['final_state']),np.asarray(b['final_state'])
            row['comparison_scaled_state_distance']=float(np.max(abs(u-v)/(atol+rtol*np.maximum(abs(u),abs(v)))))
            row['relative_output_differences']={k:abs(x[k]-y[k])/abs(x[k])
                for k in ('indicated_power_w','indicated_thermal_efficiency','total_mass_kg')}
            row['cycles']=[x['convergence']['cycles_completed'],y['convergence']['cycles_completed']]
            row['validity_equal']=all((x[s]==y[s] if x[s] is None or y[s] is None else x[s][k]==y[s][k]) for s,k in (
                ('validity','verdict'),('validity','failed_criteria'),
                ('microtube_gas_domains','model_validity'),('microtube_gas_domains','failed_criteria')))
            row['state_within_comparison_scale']=row['comparison_scaled_state_distance']<=1
            row['outputs_within_relative_comparison_scale']=max(row['relative_output_differences'].values())<=rtol
        results.append(row)
    if not results: raise ValueError('No comparable benchmark cases found.')
    old=json.loads((reference/'summary.json').read_text());new=json.loads((candidate/'summary.json').read_text())
    speedup={k:old[k]['median_seconds']/v['median_seconds'] for k,v in new.items()
             if k not in ('interrupted','invalid_domain')}
    passed=all(all(row.get(key,True) for key in (
        'status_equal','exception_equal','message_equal','validity_equal',
        'state_within_comparison_scale','outputs_within_relative_comparison_scale')) for row in results)
    return dict(passed=passed,comparison_rtol=rtol,comparison_atol=atol,runs=results,median_speedup=speedup,
        note='Rounding can change adaptive ODE steps. Periodic state tolerance is not an output-error bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('reference',type=Path);p.add_argument('candidate',type=Path);p.add_argument('manifest',type=Path)
    a=p.parse_args();result=compare(a.reference,a.candidate,a.manifest)
    print(json.dumps(result,indent=2))
    if not result['passed']: raise SystemExit(1)
