"""Fine 245–255 degree harmonic phase scan with unchanged K2 wall physics."""
from dataclasses import replace, asdict
import json
import math
import time

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.kinematics import HarmonicVolumeKinematics
from compare_motor_motion_laws_stage7A5 import A5_CAMPAIGN, ROOT, _candidate_design_and_mass, _same_inventory_design, _evaluate
from optimize_motor_exchanger_asymmetry_stageK1 import _scaled_exchanger, _hardware_metrics


def main():
    definition=CampaignDefinition(A5_CAMPAIGN)
    base,mass,_,_,_= _candidate_design_and_mass(definition)
    k2=json.loads((ROOT/'outputs/motor_exchanger_asymmetry_stageK2.json').read_text())['best_feasible']
    base=replace(base,heat_in=_scaled_exchanger(base.heat_in,k2['k_i']),heat_out=_scaled_exchanger(base.heat_out,k2['k_o']))
    report=dict(experiment='Fine harmonic phase scan with K2 dynamic-wall exchangers',
        phase_min_degrees=245.,phase_max_degrees=255.,phase_step_degrees=.5,
        total_mass_kg=mass,angular_speed=base.configuration.angular_speed,
        heat_in=_hardware_metrics(base.heat_in),heat_out=_hardware_metrics(base.heat_out),
        wall_numerical_settings=asdict(definition.wall_numerical_settings),
        maximum_cycles=base.configuration.numerical.maximum_cycles,
        efficiency_boundary='External-source heat; indicated gas work; mechanical losses unknown',results=[])
    output=ROOT/'outputs/motor_harmonic_k2_fine_scan.json'
    state=None; source=None; started=time.monotonic()
    for phase in [245.+.5*i for i in range(21)]:
        kin=HarmonicVolumeKinematics(base.configuration.machine_volumes.small_cylinder,base.configuration.machine_volumes.large_cylinder,math.radians(phase))
        design=_same_inventory_design(base,mass,kin)
        tick=time.monotonic()
        result,new_state=_evaluate(f'harmonic_k2_{phase:.1f}',design,definition,initial_state=state)
        result.update(phase_degrees=phase,elapsed_seconds=time.monotonic()-tick,warm_start_phase_degrees=source,
                      last_complete_state=new_state.tolist() if new_state is not None else None)
        report['results'].append(result)
        if result['status']=='converged':
            state=new_state; source=phase
        valid=[r for r in report['results'] if r['status']=='converged' and r.get('indicated_thermal_efficiency') is not None and r.get('validity',{}).get('verdict')=='valid']
        report['best']=max(valid,key=lambda r:r['indicated_thermal_efficiency']) if valid else None
        report['elapsed_seconds']=time.monotonic()-started
        temporary=output.with_suffix('.tmp'); temporary.write_text(json.dumps(report,indent=2)+'\n'); temporary.replace(output)
        print(phase,result['status'],result.get('indicated_thermal_efficiency'),result.get('indicated_power_w'),result['convergence']['cycles_completed'],flush=True)
    print('BEST',report['best']['phase_degrees'],report['best']['indicated_thermal_efficiency'],flush=True)


if __name__=='__main__': main()
