"""Isolate centered conventional rod shortening at fixed crank phase in K2."""
from dataclasses import asdict, replace
import json
import time
import numpy as np
from evaluate_slider_crank_k2 import (SliderMotion, SliderCrankKinematics, CampaignDefinition,
    A5_CAMPAIGN, ROOT, _candidate_design_and_mass, _scaled_exchanger,
    _same_inventory_design, _evaluate)


def main():
    definition=CampaignDefinition(A5_CAMPAIGN)
    base,mass,_,_,_=_candidate_design_and_mass(definition)
    k2=json.loads((ROOT/'outputs/motor_exchanger_asymmetry_stageK2.json').read_text())['best_feasible']
    fitted=json.loads((ROOT/'outputs/slider_crank_target_comparison.json').read_text())
    prior=json.loads((ROOT/'outputs/slider_crank_k2.json').read_text())
    baseline=prior['results']['conventional']
    hardware=replace(base,heat_in=_scaled_exchanger(base.heat_in,k2['k_i']),heat_out=_scaled_exchanger(base.heat_out,k2['k_o']))
    report=dict(experiment='Conventional centered rods shortened at fixed phase, K2',
        phase_degrees=baseline['small_minus_large_crank_angle_degrees'],
        fixed_total_mass_kg=mass,baseline=baseline,baseline_inputs=fitted,
        hardware_basis={k:prior[k] for k in ('heat_in','heat_out','wall_numerical_settings')},
        numerical=asdict(base.configuration.numerical),results=[],
        limitations=['Fixed phase, not a phase optimization','Indicated work; mechanical losses unknown','No cylinder layout changes modeled'])
    state=np.array(baseline['last_complete_state']); source='conventional_10'
    started=time.monotonic()
    for ratio in (6.,4.,3.,2.):
        motions=[]
        for side in ('small','large'):
            raw=fitted['sides'][side]['conventional']
            motions.append(SliderMotion(ratio,0.,raw['phase_rad'],False))
        kin=SliderCrankKinematics(*motions,base.configuration.machine_volumes.small_cylinder,base.configuration.machine_volumes.large_cylinder)
        design=_same_inventory_design(hardware,mass,kin)
        tick=time.monotonic()
        result,new_state=_evaluate(f'conventional_rod_{ratio}',design,definition,initial_state=state)
        result.update(rod_over_crank=ratio,warm_start_source=source,elapsed_seconds=time.monotonic()-tick,
                      last_complete_state=new_state.tolist() if new_state is not None else None)
        report['results'].append(result)
        if result['status']=='converged':state=new_state;source=f'conventional_{ratio}'
        report['elapsed_seconds']=time.monotonic()-started
        path=ROOT/'outputs/short_rods_k2.json';temp=path.with_suffix('.tmp')
        temp.write_text(json.dumps(report,indent=2)+'\n');temp.replace(path)
        print(ratio,result['status'],result.get('indicated_thermal_efficiency'),result.get('indicated_power_w'),result['convergence']['cycles_completed'],flush=True)


if __name__=='__main__':main()
