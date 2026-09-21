"""Read-only feasibility-contract check on saved exact benchmark evaluations."""
from dataclasses import replace
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from benchmark_solver_acceleration import ROOT,DEFAULT_MANIFEST
from refine_motor_four_stage_hx9d_variable_gas import _load_basis,_build_design
from dada_solver.campaign.evaluator import MachineEvaluator,json_values
from dada_solver.exchangers.wall_cycle import WallDiagnosticCycle,wall_cycle_performance
from dada_solver.sizing.evaluator import EvaluationStatus
from dada_solver.validity import ValidityReport,ValidityVerdict
from dada_solver.topology import CycleTopologyClassification


def check():
    definition,base=_load_basis();helper=MachineEvaluator(definition)
    manifest=json.loads(DEFAULT_MANIFEST.read_text());results=[]
    for case in manifest['cases']:
        if case['name'] not in ('production_warm','nearby_warm','production_cold','smooth_four_bar'):continue
        design=_build_design(base,case['parameters'])
        if case.get('motion')=='base':design=replace(design,kinematics=base.kinematics)
        wrapper=design.build();reports=[]
        for directory in ('p1_baseline','p2_shared'):
            path=ROOT/'outputs/solver_acceleration_stage5'/directory/(case['name']+'_0')
            row=json.loads(path.with_suffix('.json').read_text());r=row['result']
            with np.load(path.with_suffix('.npz')) as f:
                cycle=WallDiagnosticCycle(f['angles'],f['trajectory'])
                performance=wall_cycle_performance(wrapper,f['trajectory'])
            d=dict(r['diagnostics'])
            for key in ('pressure_extrema','temperature_extrema','mass_flow_extrema'):
                d[key]={k:SimpleNamespace(**v) for k,v in d[key].items()}
            topology=dict(d['topology']);topology['classification']=CycleTopologyClassification(topology['classification'])
            d['topology']=SimpleNamespace(**topology)
            v=dict(r['validity']);v['verdict']=ValidityVerdict(v['verdict'])
            evaluation=SimpleNamespace(configuration=design.configuration,usable=True,status=EvaluationStatus.CONVERGED,
                periodic=SimpleNamespace(message=r['message']),performance=performance,diagnostics=SimpleNamespace(**d),
                validity=ValidityReport(**v),model=wrapper.model,cycle=cycle)
            domains=r['microtube_gas_domains']
            ok=domains['model_validity']=='valid' and r['maximum_tube_mach_number']<=design.configuration.validity.maximum_mach_number
            extra=[dict(name='microtube_model_domain',margin=1. if ok else -1.,satisfied=ok,available=True)]
            reports.append(json_values(helper._assessment(evaluation,{},None,None,r['convergence'],row['final_state'],extra)))
        results.append(dict(case=case['name'],exact=reports[0]==reports[1],status=reports[1]['status'],
            constraints=reports[1]['constraints'],objective=reports[1]['objective'],conservation=reports[1]['metrics']['conservation']))
    return dict(passed=all(r['exact'] for r in results),runs=results,
        note='Configured assessment only: no campaign, history, cache or warm-start record is written.')

if __name__=='__main__':
    result=check();print(json.dumps(result,indent=2))
    if not result['passed']:raise SystemExit(1)
