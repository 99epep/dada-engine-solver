"""Bounded post-Stage-5 Anderson experiment; no campaign parsing or persistence."""
import argparse
from dataclasses import asdict,replace
import json
from pathlib import Path
import time
from types import SimpleNamespace
import numpy as np
from benchmark_solver_acceleration import ROOT,DEFAULT_MANIFEST,digest,environment
from benchmark_solver_stage4 import build_frozen_design
from refine_motor_four_stage_hx9d_variable_gas import _load_basis
from compare_motor_motion_laws_stage7A5 import _evaluate
from dada_solver.exchangers.wall_cycle import WallCycleNumericalSettings,solve_periodic_wall_motor,wall_cycle_performance
from dada_solver.exchangers.wall_iteration import AndersonAccelerationSettings
from dada_solver.wall_backend import WallBackendSettings,WallRHS,backend_identity
from dada_solver.campaign.evaluator import json_values


def run(output,screen=False,memory=3,damping=.8,repeats=3,tight=False):
    if output.exists():raise FileExistsError('Use a new Anderson experiment directory.')
    output.mkdir(parents=True)
    manifest_path=ROOT/'outputs/solver_acceleration_stage2/manifest_tight_periodic.json' if tight else DEFAULT_MANIFEST
    manifest=json.loads(manifest_path.read_text())
    for name,expected in manifest['input_hashes'].items():
        if digest(ROOT/name)!=expected:raise ValueError('Changed frozen input: '+name)
    _,base=_load_basis();backend=WallBackendSettings('numba')
    warm=next(c for c in json.loads(DEFAULT_MANIFEST.read_text())['cases'] if c['name']=='production_warm')
    rhs=WallRHS(build_frozen_design(base,warm).build(),backend)
    started=time.perf_counter();rhs(0.,np.r_[warm['initial_state'],np.zeros(5)])
    metadata=dict(environment=environment(),backend=backend_identity(backend),excluded_compilation_seconds=time.perf_counter()-started,
        source_hashes={str(p.relative_to(ROOT)):digest(p) for p in sorted((ROOT/'src').rglob('*.py'))},
        manifest=str(manifest_path.relative_to(ROOT)),manifest_sha256=digest(manifest_path),
        script_sha256=digest(Path(__file__)),screen=screen,tight=tight,
        note='Exact cache and shared replay enabled; ordinary frozen budget retained, including for plain iteration.')
    (output/'environment.json').write_text(json.dumps(metadata,indent=2)+'\n')
    cases=[c for c in manifest['cases'] if c['name'] in (('production_cold','smooth_four_bar') if screen or tight else
           ('production_warm','nearby_warm','production_cold','smooth_four_bar'))]
    if not screen and not tight:
        extra=ROOT/'outputs/solver_acceleration_stage5/manifest_sixbar.json'
        cases+=json.loads(extra.read_text())['cases'];metadata['six_bar_manifest_sha256']=digest(extra)
        (output/'environment.json').write_text(json.dumps(metadata,indent=2)+'\n')
    policies=[('plain',None),('production',None)]
    if screen:policies += [(f'anderson_m{m}_b{b}',AndersonAccelerationSettings(memory=m,damping=b)) for m in (2,3,4) for b in (.5,.8,1.)]
    else:policies=[('production',None),(f'anderson_m{memory}_b{damping}',AndersonAccelerationSettings(memory=memory,damping=damping))]
    rows=[]
    for case in cases:
        for repeat in range(1 if screen else repeats):
            # Alternate complete-validation order without overlapping timers.
            ordered=policies if repeat%2==0 else list(reversed(policies))
            for policy,acceleration in ordered:
                design=build_frozen_design(base,case)
                raw=dict(manifest['numerical_settings']);raw['integration_absolute_tolerances']=tuple(raw['integration_absolute_tolerances'])
                settings=replace(WallCycleNumericalSettings(**raw),accelerate_walls=policy!='plain')
                initial=np.array(case['initial_state']);records=[];capture=[]
                started=time.perf_counter();failure=None;result=None;performance=None
                try:
                    if screen:
                        wrapper=design.build()
                        periodic=solve_periodic_wall_motor(wrapper,initial,maximum_cycles=case['maximum_cycles'],settings=settings,
                            backend=backend,anderson_acceleration=acceleration,statistics_callback=records.append)
                        capture.append(periodic);state=periodic.last_complete_state
                    else:
                        definition=SimpleNamespace(wall_numerical_settings=settings,wall_backend=backend,anderson_acceleration=acceleration)
                        result,state=_evaluate(case['name'],design,definition,initial_state=initial,
                            statistics_callback=records.append,periodic_observer=capture.append)
                    elapsed=time.perf_counter()-started
                    periodic=capture[0]
                    if periodic.converged:performance=json_values(asdict(wall_cycle_performance(design.build(),periodic.trajectory)))
                except (ValueError,RuntimeError,ArithmeticError) as error:
                    elapsed=time.perf_counter()-started;failure=dict(type=type(error).__name__,message=str(error));state=None
                    periodic=capture[0] if capture else None
                counts=periodic.anderson_statistics if periodic and acceleration else {}
                attempts=counts.get('attempted_cycles',max([r.get('cycle',0) for r in records]+[len(periodic.history) if periodic else 0]))
                retained=counts.get('retained_cycles',len(periodic.history) if periodic else 0)
                decisions=[r for r in records if r['phase']=='anderson_decision']
                periodic_time=elapsed if screen else sum(r['elapsed_seconds'] for r in records if r['phase']=='periodic_integration')
                row=dict(case=case['name'],policy=policy,repeat=repeat,configuration=json_values(asdict(design.configuration)),
                    parameters=case['parameters'],initial_state=initial.tolist(),numerical_settings=asdict(settings),
                    anderson_settings=asdict(acceleration) if acceleration else None,
                    maximum_cycles=case['maximum_cycles'],status=periodic.status if periodic else 'failure',exception=failure,
                    attempted_cycles=attempts,retained_cycles=retained,anderson_statistics=counts,
                    elapsed_seconds=elapsed,periodic_seconds=periodic_time,decision_seconds=sum(r['elapsed_seconds'] for r in decisions),
                    history=list(periodic.history) if periodic else [],solver_statistics=records,
                    final_state=state.tolist() if state is not None else None,
                    final_normalized_residual=next((h['normalized_state_error'] for h in reversed(periodic.history) if h.get('retained',True)),None) if periodic else None,
                    result=result,performance=performance,backend_statistics=periodic.backend_statistics if periodic else {})
                key=f'{case["name"]}_{policy}_{repeat}'
                (output/(key+'.json')).write_text(json.dumps(row,indent=2)+'\n');rows.append(row)
                if periodic and periodic.trajectory is not None:
                    np.savez_compressed(output/(key+'.npz'),angles=periodic.angles,trajectory=periodic.trajectory)
                print(case['name'],policy,repeat,row['status'],attempts,retained,f'{periodic_time:.3f}s',
                      'rollbacks',counts.get('rollback_count',0),flush=True)
    (output/'summary.json').write_text(json.dumps([dict(case=r['case'],policy=r['policy'],repeat=r['repeat'],status=r['status'],
        attempts=r['attempted_cycles'],retained=r['retained_cycles'],periodic_seconds=r['periodic_seconds'],
        decision_seconds=r['decision_seconds'],counts=r['anderson_statistics']) for r in rows],indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--screen',action='store_true');p.add_argument('--tight',action='store_true')
    p.add_argument('--memory',type=int,default=3);p.add_argument('--damping',type=float,default=.8)
    p.add_argument('--repeats',type=int,default=3)
    a=p.parse_args();run(a.output,a.screen,a.memory,a.damping,a.repeats,a.tight)
