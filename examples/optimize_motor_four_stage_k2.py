"""Persistent seven-variable four-stage motor search on fixed K2 hardware.

Run with --budget-seconds 480; rerunning resumes the append-only history.
No mechanical efficiency or altered thermodynamic tolerances are introduced.
"""
from dataclasses import replace
import argparse
import json
import math
from pathlib import Path
import time
import hashlib
import numpy as np
from scipy.stats import qmc
from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.four_stage_kinematics import FourStageVolumeKinematics
from dada_solver.integration import IntegrationInterrupted
from compare_motor_motion_laws_stage7A5 import A5_CAMPAIGN, ROOT, _candidate_design_and_mass, _same_inventory_design, _evaluate
from optimize_motor_exchanger_asymmetry_stageK1 import _scaled_exchanger
from optimize_motor_piecewise_stageP3 import _feasibility

NAMES = ('t1','t2','t3','a_l','b_l','a_s','b_s')


def _symmetry_diagnostics(parameters):
    """Return passive diagnostics for S/L intermediate-volume symmetry."""
    delta_a = float(parameters['a_l'] - parameters['a_s'])
    delta_b = float(parameters['b_l'] - parameters['b_s'])
    return {
        'delta_a': delta_a,
        'delta_b': delta_b,
        'rms': math.sqrt((delta_a * delta_a + delta_b * delta_b) / 2.0),
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--budget-seconds',type=float,default=480)
    parser.add_argument('--evaluations',type=int,default=32)
    parser.add_argument('--candidate-seconds',type=float,default=90,
                        help='Cooperative per-candidate limit, independent of maximum cycles.')
    args=parser.parse_args()
    directory=ROOT/'outputs/motor_four_stage_k2'
    directory.mkdir(exist_ok=True)
    history_path=directory/'history.jsonl'
    reference=json.loads((ROOT/'outputs/motor_champion_sixbar_k2.json').read_text())
    definition=CampaignDefinition(A5_CAMPAIGN)
    base,mass,_,_,_= _candidate_design_and_mass(definition)
    basis=reference['thermodynamic_basis']
    base=replace(base,heat_in=_scaled_exchanger(base.heat_in,basis['stage_k2_k_i']),
                 heat_out=_scaled_exchanger(base.heat_out,basis['stage_k2_k_o']))
    limits=base.configuration.machine_volumes
    target=np.genfromtxt(ROOT/'outputs/motor_champion_motion_target.csv',delimiter=',',names=True)
    angle=target['theta_rad']; s=target['small_fraction_0_1']; l=target['large_fraction_0_1']
    origin=angle[l.argmax()]
    t1=((origin-angle[s.argmax()])%(2*math.pi))/(2*math.pi)
    t2=((origin-angle[l.argmin()])%(2*math.pi))/(2*math.pi)
    t3=((origin-angle[s.argmin()])%(2*math.pi))/(2*math.pi)
    at=lambda values,t: float(np.interp((origin-2*math.pi*t)%(2*math.pi),angle,values))
    seed=np.array([t1,t2,t3,at(l,t3),at(l,t1),at(s,t2),at(s,0)])
    identity=dict(names=NAMES,seed=seed.tolist(),minimum_stage_fraction=.02,
                  basis=basis, source_reference_sha256=hashlib.sha256((ROOT/'outputs/motor_champion_sixbar_k2.json').read_bytes()).hexdigest(),
                  strategy='Seed, alternating global Sobol and incumbent-centered Sobol; seed 2026; all seven coordinates vary.',
                  minimum_indicated_power_W=40, angle_convention='Forward motor time',
                  volume_zero='K2 minimum enclosed volume', hydraulic_loss_multiplier=1)
    identity=json.loads(json.dumps(identity))
    definition_path=directory/'definition.json'
    if definition_path.exists() and json.loads(definition_path.read_text()) != identity:
        raise ValueError('Search definition changed; use a new output directory.')
    definition_path.write_text(json.dumps(identity,indent=2)+'\n')
    history=[json.loads(line) for line in history_path.read_text().splitlines()] if history_path.exists() else []
    eligible=[r for r in history if r['feasible']]
    best=max(eligible,key=lambda r:r['result']['indicated_thermal_efficiency']) if eligible else None
    initial_best=best
    points=qmc.Sobol(d=7,scramble=True,seed=2026).random_base2(12)
    start=time.monotonic(); deadline=start+args.budget_seconds
    completed=0
    def progress(_):
        if time.monotonic()>=deadline: raise IntegrationInterrupted('Phase wall-clock budget exhausted.')
        if time.monotonic()>=candidate_deadline: raise IntegrationInterrupted('Candidate wall-clock budget exhausted; result remains unknown.')
    index=max((r['index'] for r in history),default=-1)+1
    while completed<args.evaluations and time.monotonic()<deadline:
        u=points[index]
        if index==0: params=seed.copy(); kind='target_extrema_seed'
        elif index%4==0:
            # Uniform ordered knots with a 2% minimum duration for all stages.
            remaining=.92; gaps=[]
            for j in range(3):
                portion=remaining*(1-(1-u[j])**(1/(3-j)))
                gaps.append(.02+portion); remaining-=portion
            params=np.r_[np.cumsum(gaps),u[3:]]; kind='global_sobol'
        else:
            center=np.array(list(best['parameters'].values())) if best else seed
            radius=.12 if index<16 else .06
            params=center+(2*u-1)*radius; kind='incumbent_neighborhood_sobol'
        params=np.asarray(params)
        if (np.any(params[3:]<0) or np.any(params[3:]>1)
            or np.min(np.diff(np.r_[0,params[:3],1]))<.02):
            index+=1; continue
        mapping=dict(zip(NAMES,map(float,params)))
        candidate_id=hashlib.sha256(json.dumps(mapping,sort_keys=True).encode()).hexdigest()
        if any(r['candidate_id']==candidate_id and r['result']['status']!='interrupted' for r in history):
            index+=1; continue
        kin=FourStageVolumeKinematics(limits.small_cylinder,limits.large_cylinder,**mapping)
        design=_same_inventory_design(base,mass,kin)
        reusable=[r for r in history if r['result']['status']=='converged' and r['last_complete_state']]
        source=min(reusable,key=lambda r:np.linalg.norm(params-np.array(list(r['parameters'].values())))) if reusable else None
        initial=source['last_complete_state'] if source else reference['last_complete_state']
        before=time.monotonic()
        candidate_deadline=min(deadline,before+args.candidate_seconds)
        try:
            result,state=_evaluate(f'four_stage_{index}',design,definition,initial_state=np.array(initial),progress_callback=progress)
        except (ValueError,RuntimeError) as exc:
            result={'status':'integration_failure','message':str(exc)};state=None
        feasible,reasons=_feasibility(result) if result['status']=='converged' else (False,[])
        record=dict(index=index,candidate_id=candidate_id,kind=kind,parameters=mapping,
                    symmetry=_symmetry_diagnostics(mapping),result=result,
                    feasible=feasible,physical_constraint_failures=reasons,
                    elapsed_seconds=time.monotonic()-before,
                    warm_start_source=source['candidate_id'] if source else 'sixbar_reference',
                    last_complete_state=state.tolist() if state is not None else None)
        with history_path.open('a') as stream:
            stream.write(json.dumps(record)+'\n');stream.flush()
        history.append(record);completed+=1
        if feasible and (best is None or result['indicated_thermal_efficiency']>best['result']['indicated_thermal_efficiency']): best=record
        print(json.dumps(dict(index=index,status=result['status'],feasible=feasible,
                              efficiency=result.get('indicated_thermal_efficiency'),power_W=result.get('indicated_power_w'),
                              elapsed_seconds=record['elapsed_seconds'])),flush=True)
        report=dict(best_feasible=best,
                    best_symmetry=_symmetry_diagnostics(best['parameters']) if best else None,
                    phase_start_best=initial_best,attempted_total=len(history),
                    requested_duration_seconds=args.budget_seconds,actual_duration_seconds=time.monotonic()-start,
                    reference_efficiency=basis['reference_efficiency'],reference_sixbar_efficiency=reference['result']['indicated_thermal_efficiency'])
        (directory/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        if time.monotonic()>=deadline: break
        index+=1
    print('Saved '+str(directory/'report.json'),flush=True)


if __name__=='__main__': main()
