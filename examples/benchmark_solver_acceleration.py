"""Frozen, serial Stage-2 replays; never runs or writes an optimization campaign.

Run with PYTHONPATH=src. Freeze once, then use distinct output directories for
baseline and each patch. Timers wrap coarse boundaries only, in this process.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import hashlib
import inspect
import json
from pathlib import Path
import platform
import subprocess
import time
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import scipy

from dada_solver.configuration import load_simulation_configuration
from dada_solver.exchangers.hardware import load_hardware_definition
from dada_solver.exchangers.microtube import MicrotubeExchanger
from dada_solver.exchangers.wall_cycle import WallCycleNumericalSettings
from dada_solver.integration import IntegrationInterrupted
from dada_solver.machine import MachineDesign
import dada_solver.exchangers.gas_diagnostics as gas_diagnostics
import compare_motor_motion_laws_stage7A5 as comparison
from refine_motor_four_stage_hx9d_variable_gas import _load_basis, _build_design

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT/'outputs/solver_acceleration_stage2/manifest.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def environment():
    return dict(platform=platform.platform(), python=platform.python_version(),
                numpy=np.__version__, scipy=scipy.__version__,
                debian_version=Path('/etc/debian_version').read_text().strip(),
                commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())


def legacy_design():
    config = load_simulation_configuration(ROOT/'examples/motor_demonstrator_original_325c.toml')
    _, bank, hi, ho = load_hardware_definition(
        (ROOT/'examples/motor_hardware_parallel_325c.toml').read_text(), config.gas.heat_capacity_cp)
    return MachineDesign(config,
        MicrotubeExchanger(bank,hi,config.hydraulics.cold_to_large_valve_cda),
        MicrotubeExchanger(bank,ho,config.hydraulics.hot_to_small_valve_cda))


def freeze(path):
    if path.exists():
        raise FileExistsError('A frozen manifest must not be overwritten.')
    artifact = json.loads((ROOT/'outputs/solver_performance_stage1_20260917.json').read_text())
    definition, base = _load_basis()
    params = artifact['parameters']
    design = _build_design(base, params)
    initial = np.asarray(artifact['initial_state'])
    nearby = dict(params, a_l=params['a_l']+.0002)
    nearby_design = _build_design(base,nearby)
    nearby_initial = initial.copy()
    target = comparison._uniform_wall_initial(nearby_design.configuration,nearby_design.build())
    nearby_initial[:8] *= target[:8:2].sum()/initial[:8:2].sum()
    # Use the study-angle family from the base design, avoiding double reversal.
    smooth = replace(design,kinematics=base.kinematics)
    cold = comparison._uniform_wall_initial(design.configuration,design.build())
    invalid = initial.copy(); invalid[5] = invalid[4]*design.configuration.gas.heat_capacity_cv*2000
    legacy = json.loads((ROOT/'outputs/motor_parallel_325c_screening.json').read_text())
    legacy_state = np.asarray(legacy['final_state'])
    wrapper = legacy_design().build()
    for i,side,wall in ((8,'H_i',wrapper.heat_in),(9,'H_o',wrapper.heat_out)):
        legacy_state[i] *= wall.wall_capacity_j_k/legacy['hardware'][side]['wall_capacity_j_k']
    cases = []
    def add(name,state,**extra):
        cases.append(dict(name=name,parameters=params,initial_state=state.tolist(),
                          maximum_cycles=design.configuration.numerical.maximum_cycles,**extra))
    add('production_warm',initial)
    add('nearby_warm',nearby_initial); cases[-1]['parameters']=nearby
    add('production_cold',cold)
    add('invalid_domain',invalid)
    add('interrupted',cold,interrupt_first_progress=True)
    add('smooth_four_bar',comparison._uniform_wall_initial(smooth.configuration,smooth.build()),motion='base')
    add('legacy',legacy_state,legacy=True); cases[-1]['maximum_cycles']=1
    paths = list((ROOT/'examples').glob('*.toml'))
    paths += [ROOT/'outputs'/x for x in ('solver_performance_stage1_20260917.json',
        'motor_four_stage_thermo3d/report.json','motor_mechanics_stage7A5_edge/report.json',
        'motor_parallel_325c_screening.json')]
    report=json.loads((ROOT/'outputs/motor_mechanics_stage7A5_edge/report.json').read_text())
    paths.append(ROOT/'outputs/motor_mechanics_stage7A5_edge/candidates'/
                 (report['best_at_phase_end']['candidate_id']+'.json'))
    manifest=dict(schema_version=1,environment=environment(),
        source_candidate=artifact['source_candidate'],
        input_hashes={str(p.relative_to(ROOT)):digest(p) for p in sorted(paths)},
        numerical_settings=artifact['numerical_settings'],cases=cases)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(manifest,indent=2)+'\n')
    print('Frozen',path,flush=True)


def run(manifest_path, output, names, repeats, measure_rhs=False, adaptive_wall=False, adaptive_controls=None):
    manifest=json.loads(manifest_path.read_text())
    for name, expected in manifest['input_hashes'].items():
        if digest(ROOT/name)!=expected: raise ValueError('Changed frozen input: '+name)
    if output.exists(): raise FileExistsError('Use a new measurement directory.')
    output.mkdir(parents=True)
    definition,base=_load_basis()
    settings=dict(manifest['numerical_settings'])
    settings['integration_absolute_tolerances']=tuple(settings['integration_absolute_tolerances'])
    definition.wall_numerical_settings=WallCycleNumericalSettings(**settings)
    adaptive = None
    if adaptive_wall:
        from dada_solver.exchangers.wall_iteration import AdaptiveWallAccelerationSettings
        controls = json.loads(adaptive_controls.read_text()) if adaptive_controls else {}
        adaptive = AdaptiveWallAccelerationSettings(**controls)
    rows=[]
    supports_statistics = "statistics_callback" in inspect.signature(comparison.solve_periodic_wall_motor).parameters
    metadata=dict(environment=environment(),manifest_sha256=digest(manifest_path),
        adaptive_wall_acceleration=asdict(adaptive) if adaptive is not None else None,
        source_hashes={str(p.relative_to(ROOT)):digest(p) for p in sorted((ROOT/'src').rglob('*.py'))})
    (output/'environment.json').write_text(json.dumps(metadata,indent=2)+'\n')
    for case in manifest['cases']:
        if names and case['name'] not in names: continue
        for repeat in range(repeats):
            start=time.perf_counter()
            design=legacy_design() if case.get('legacy') else _build_design(base,case['parameters'])
            if case.get('motion')=='base': design=replace(design,kinematics=base.kinematics)
            design=replace(design,configuration=replace(design.configuration,
                numerical=replace(design.configuration.numerical,maximum_cycles=case['maximum_cycles'])))
            control=SimpleNamespace(wall_numerical_settings=(WallCycleNumericalSettings()
                if case.get('legacy') else definition.wall_numerical_settings))
            sections={}; captured={}; solver_statistics=[]; cycle_endpoints=[]
            def timed(name,function):
                def wrapped(*args,**kwargs):
                    before=time.perf_counter()
                    if name=='integration' and supports_statistics:
                        kwargs['statistics_callback']=solver_statistics.append
                        kwargs['measure_rhs_time']=measure_rhs
                        if adaptive is not None:
                            kwargs['adaptive_acceleration']=adaptive
                        kwargs['cycle_callback']=lambda cycle,end,error,item: cycle_endpoints.append(
                            dict(cycle=cycle,state=end.tolist(),normalized_state_error=error))
                    try:
                        result=function(*args,**kwargs)
                        if name=='integration': captured['periodic']=result
                        return result
                    finally: sections.setdefault(name,[]).append(time.perf_counter()-before)
                return wrapped
            def progress(_):
                if case.get('interrupt_first_progress'):
                    raise IntegrationInterrupted('Frozen benchmark deadline at first progress callback.')
            try:
                with patch.object(MachineDesign,'build',timed('build',MachineDesign.build)), \
                     patch.object(comparison,'solve_periodic_wall_motor',timed('integration',comparison.solve_periodic_wall_motor)), \
                     patch.object(comparison,'extract_cycle_diagnostics',timed('generic_diagnostics',comparison.extract_cycle_diagnostics)), \
                     patch.object(comparison,'assess_cycle_validity',timed('generic_validity',comparison.assess_cycle_validity)), \
                     patch.object(gas_diagnostics,'cycle_microtube_diagnostics',timed('microtube_diagnostics',gas_diagnostics.cycle_microtube_diagnostics)):
                    result,state=comparison._evaluate(case['name'],design,control,
                        initial_state=np.asarray(case['initial_state']),progress_callback=progress)
                error=None
            except (ValueError,RuntimeError,ArithmeticError) as exc:
                result=dict(status=type(exc).__name__,message=str(exc));state=None;error=type(exc).__name__
            elapsed=time.perf_counter()-start
            key=f"{case['name']}_{repeat}"
            periodic=captured.get('periodic')
            if periodic is not None and periodic.trajectory is not None:
                np.savez_compressed(output/(key+'.npz'),angles=periodic.angles,trajectory=periodic.trajectory)
            row=dict(case=case['name'],repeat=repeat,elapsed_seconds=elapsed,sections_seconds=sections,
                solver_statistics=solver_statistics, rhs_timing_enabled=measure_rhs,
                cycle_endpoints=cycle_endpoints,
                result=result,final_state=None if state is None else state.tolist(),exception=error)
            serialization_started=time.perf_counter()
            json.dumps(row)
            row['report_serialization_seconds']=time.perf_counter()-serialization_started
            rows.append(row)
            (output/(key+'.json')).write_text(json.dumps(row,indent=2)+'\n')
            print(key, result['status'],f'{elapsed:.3f}s',result.get('indicated_power_w'),flush=True)
    summary={name:dict(median_seconds=float(np.median([r['elapsed_seconds'] for r in rows if r['case']==name])),
                      maximum_seconds=max(r['elapsed_seconds'] for r in rows if r['case']==name),
                      repeats=sum(r['case']==name for r in rows)) for name in {r['case'] for r in rows}}
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,default=DEFAULT_MANIFEST)
    parser.add_argument('--freeze',action='store_true')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--cases',nargs='+')
    parser.add_argument('--repeats',type=int,default=3)
    parser.add_argument('--adaptive-wall',action='store_true',help='Opt into safeguarded four-endpoint wall acceleration.')
    parser.add_argument('--adaptive-controls',type=Path,help='JSON controls for a separately identified adaptive experiment.')
    parser.add_argument('--measure-rhs',action='store_true',help='Opt-in per-RHS timer; report separately from uninstrumented latency.')
    args=parser.parse_args()
    if args.freeze: freeze(args.manifest)
    else:
        if args.output is None: parser.error('--output is required')
        run(args.manifest,args.output,args.cases,args.repeats,args.measure_rhs,
            args.adaptive_wall or args.adaptive_controls is not None,args.adaptive_controls)
