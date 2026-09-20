"""Frozen Stage-4 production-backend replay, without process-global patching."""
import argparse
from dataclasses import replace,asdict
import json
from pathlib import Path
import time
from types import SimpleNamespace
import numpy as np
from benchmark_solver_acceleration import DEFAULT_MANIFEST,ROOT,digest,environment,legacy_design
from refine_motor_four_stage_hx9d_variable_gas import _load_basis,_build_design
from dada_solver.exchangers.wall_cycle import WallCycleNumericalSettings
from dada_solver.wall_backend import WallBackendSettings,backend_identity
from dada_solver.integration import IntegrationInterrupted
from compare_motor_motion_laws_stage7A5 import _evaluate


def run(output,manifest_path=DEFAULT_MANIFEST,backend='python',repeats=3,cases=None,profile=False,**backend_options):
    if output.exists(): raise FileExistsError('Use a new measurement directory.')
    manifest=json.loads(manifest_path.read_text())
    for name,expected in manifest['input_hashes'].items():
        if digest(ROOT/name)!=expected: raise ValueError('Changed frozen input: '+name)
    output.mkdir(parents=True)
    preparation=time.perf_counter();_,base=_load_basis()
    mode=WallBackendSettings(name=backend,profile=profile,**backend_options)
    identity=backend_identity(mode)
    metadata=dict(environment=environment(),manifest_sha256=digest(manifest_path),backend=identity,
        preparation_seconds=time.perf_counter()-preparation,
        source_hashes={str(p.relative_to(ROOT)):digest(p) for p in sorted((ROOT/'src').rglob('*.py'))})
    (output/'environment.json').write_text(json.dumps(metadata,indent=2)+'\n')
    rows=[]
    for case in manifest['cases']:
        if cases and case['name'] not in cases: continue
        for repeat in range(repeats):
            start=time.perf_counter();records=[];capture=[]
            design=legacy_design() if case.get('legacy') else _build_design(base,case['parameters'])
            if case.get('motion')=='base': design=replace(design,kinematics=base.kinematics)
            design=replace(design,configuration=replace(design.configuration,
                numerical=replace(design.configuration.numerical,maximum_cycles=case['maximum_cycles'])))
            candidate_seconds=time.perf_counter()-start
            raw=dict(manifest['numerical_settings']);raw['integration_absolute_tolerances']=tuple(raw['integration_absolute_tolerances'])
            settings=WallCycleNumericalSettings() if case.get('legacy') else WallCycleNumericalSettings(**raw)
            definition=SimpleNamespace(wall_numerical_settings=settings,wall_backend=mode)
            def progress(_):
                if case.get('interrupt_first_progress'): raise IntegrationInterrupted('Frozen benchmark deadline at first progress callback.')
            try:
                result,state=_evaluate(case['name'],design,definition,initial_state=np.asarray(case['initial_state']),
                    progress_callback=progress,statistics_callback=records.append,periodic_observer=capture.append,
                    measure_rhs_time=profile)
                error=None
            except (ValueError,RuntimeError,ArithmeticError) as exc:
                result=dict(status=type(exc).__name__,message=str(exc));state=None;error=type(exc).__name__
            elapsed=time.perf_counter()-start
            key=f"{case['name']}_{repeat}"
            if capture and capture[0].trajectory is not None:
                np.savez_compressed(output/(key+'.npz'),angles=capture[0].angles,trajectory=capture[0].trajectory)
            sections={}
            for record in records:
                if record['phase'] in ('build','periodic_integration','generic_diagnostics','generic_validity','microtube_diagnostics'):
                    sections.setdefault(record['phase'],[]).append(record['elapsed_seconds'])
            row=dict(case=case['name'],repeat=repeat,elapsed_seconds=elapsed,candidate_seconds=candidate_seconds,
                sections_seconds=sections,solver_statistics=records,result=result,
                final_state=None if state is None else state.tolist(),exception=error,
                backend_statistics=capture[0].backend_statistics if capture else next(
                    (r for r in reversed(records) if r['phase']=='rhs_backend'),{}))
            before=time.perf_counter();json.dumps(row);row['report_serialization_seconds']=time.perf_counter()-before
            (output/(key+'.json')).write_text(json.dumps(row,indent=2)+'\n');rows.append(row)
            print(key,result['status'],f'{elapsed:.3f}s',result.get('indicated_power_w'),flush=True)
    summary={name:dict(median_seconds=float(np.median([r['elapsed_seconds'] for r in rows if r['case']==name])),
        maximum_seconds=max(r['elapsed_seconds'] for r in rows if r['case']==name),
        repeats=sum(r['case']==name for r in rows)) for name in {r['case'] for r in rows}}
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--manifest',type=Path,default=DEFAULT_MANIFEST)
    p.add_argument('--backend',choices=('python','numba'),default='python');p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--cases',nargs='+');p.add_argument('--profile',action='store_true')
    p.add_argument('--disk-cache',action='store_true');p.add_argument('--cache-directory')
    a=p.parse_args();run(a.output,a.manifest,a.backend,a.repeats,a.cases,a.profile,
        disk_cache=a.disk_cache,cache_directory=a.cache_directory)
