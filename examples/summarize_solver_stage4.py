"""Summarize frozen Stage-4 timings, cache behavior and residual attribution."""
import argparse
import json
from pathlib import Path


def summarize(root):
    micro=json.loads((root/'microprofile.json').read_text())['cases']
    costs={}
    for case in ('production_warm','production_cold','smooth_four_bar'):
        row=json.loads((root/'profile'/(case+'_1.json')).read_text())
        backend=row['backend_statistics'];records=row['solver_statistics']
        sections={k:sum(v) for k,v in row['sections_seconds'].items()}
        segments=[s for s in records if s['phase']=='solver_segment']
        solver=sum(s['integration_elapsed_seconds'] for s in segments)
        rhs=sum(s['rhs_elapsed_seconds'] for s in segments)
        preflight=sum(s.get('elapsed_seconds',0) for s in records if s['phase']=='integration_preflight')
        calls=backend['calls'];native_estimate=calls*micro[case]['native_seconds_per_call']
        costs[case]=dict(elapsed_seconds=row['elapsed_seconds'],
            candidate_and_model_build_seconds=row['candidate_seconds']+sections.get('build',0),
            periodic_integration_seconds=sections['periodic_integration'],
            integration_segments_seconds=solver,rhs_callback_seconds=rhs,
            kinematics_seconds=backend['kinematics_seconds'],
            kinematics_microseconds_per_call=1e6*backend['kinematics_seconds']/calls,
            kinematics_fraction_of_integration=backend['kinematics_seconds']/sections['periodic_integration'],
            compiled_dispatch_seconds=backend['dispatch_seconds'],
            native_kernel_seconds_estimate=native_estimate,
            python_compiled_crossing_seconds_estimate=backend['dispatch_seconds']-native_estimate,
            python_adapter_seconds=backend['total_rhs_seconds']-backend['kinematics_seconds']-backend['dispatch_seconds'],
            rhs_outer_wrapper_and_preflight_bookkeeping_seconds=rhs+preflight-backend['total_rhs_seconds'],
            integrator_callback_and_segment_reporting_residual_seconds=solver-rhs,
            periodic_control_residual_seconds=sections['periodic_integration']-solver-preflight-backend['preparation_seconds'],
            integration_preflight_seconds=preflight,backend_preparation_seconds=backend['preparation_seconds'],
            generic_diagnostics_seconds=sections['generic_diagnostics'],
            microtube_diagnostics_seconds=sections['microtube_diagnostics'],
            validity_seconds=sections['generic_validity'],
            report_materialization_and_other_seconds=row['elapsed_seconds']-row['candidate_seconds']-sum(sections.values()),
            json_serialization_seconds=row['report_serialization_seconds'],
            calls=calls,**{k:sum(s[k] for s in segments) for k in ('nfev','njev','nlu','accepted_steps')})
    cache={}
    for name in ('numba_final','cache_empty','cache_populated'):
        cache[name]={}
        for key in ('production_warm_0','production_warm_1','nearby_warm_0'):
            row=json.loads((root/name/(key+'.json')).read_text());b=row['backend_statistics']
            cache[name][key]=dict(elapsed_seconds=row['elapsed_seconds'],first_rhs_seconds=b['first_call_seconds'],
                backend_preparation_seconds=b['preparation_seconds'],hits=b['dispatcher_cache_hits'],
                misses=b['dispatcher_cache_misses'],cache_namespace=b['cache_namespace'])
    return dict(costs=costs,cache=cache,
        limitations=[
            'Profile timings are instrumented, using the second replay to exclude compilation.',
            'Native estimates use final-cycle sample microbenchmarks; cold-transient instruction mixes can differ.',
            'Native batch cost includes native call/allocation/loop overhead; crossing estimates are approximate.',
            'Integrator residual includes Python/SciPy callback bookkeeping and segment statistics, not LSODA alone.',
            'JSON serialization is measured separately after elapsed_seconds; report materialization is included.',
            'Dispatcher cache counters are cumulative within each process, not cache loads per candidate.'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path)
    a=p.parse_args();print(json.dumps(summarize(a.directory),indent=2))
