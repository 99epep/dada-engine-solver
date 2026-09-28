"""Evaluate the centered and offset fitted slider-crank laws with K2 physics."""
from dataclasses import dataclass, replace, asdict
import json
import argparse
import math

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.geometry import CylinderVolumeLimits
from compare_motor_motion_laws_stage7A5 import A5_CAMPAIGN, ROOT, _candidate_design_and_mass, _same_inventory_design, _evaluate
from evaluate_motor_champion_four_bar_k2 import _rescaled_initial_state
from optimize_motor_exchanger_asymmetry_stageK1 import _scaled_exchanger, _hardware_metrics


from dada_solver.slider_crank import SliderMotion, SliderCrankKinematics




def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--family', choices=['inverted','inverted_offset','conventional'], action='append')
    args=parser.parse_args()
    definition=CampaignDefinition(A5_CAMPAIGN)
    base,mass,saved,_,_= _candidate_design_and_mass(definition)
    k2=json.loads((ROOT/'outputs/motor_exchanger_asymmetry_stageK2.json').read_text())['best_feasible']
    inputs=json.loads((ROOT/'outputs/slider_crank_target_comparison.json').read_text())
    hardware=replace(base,heat_in=_scaled_exchanger(base.heat_in,k2['k_i']),heat_out=_scaled_exchanger(base.heat_out,k2['k_o']))
    report=dict(experiment='Fitted finite slider-crank motions in K2',input=inputs,
                total_mass_kg=mass,heat_in=_hardware_metrics(hardware.heat_in),heat_out=_hardware_metrics(hardware.heat_out),
                wall_numerical_settings=asdict(definition.wall_numerical_settings),results={})
    output=ROOT/'outputs/slider_crank_k2.json'
    if output.exists():
        previous=json.loads(output.read_text())
        if previous.get('input') == inputs:
            report['results']=previous['results']
    for family in (args.family or ['inverted','inverted_offset']):
        motions=[]
        for side in ['small','large']:
            raw=inputs['sides'][side][family]
            motions.append(SliderMotion(**{k:raw[k] for k in SliderMotion.__dataclass_fields__}))
        kin=SliderCrankKinematics(*motions,base.configuration.machine_volumes.small_cylinder,base.configuration.machine_volumes.large_cylinder)
        design=_same_inventory_design(hardware,mass,kin)
        initial=_rescaled_initial_state(saved,base,design)
        result,state=_evaluate('slider_crank_'+family,design,definition,initial_state=initial)
        result['last_complete_state']=state.tolist() if state is not None else None
        result['small_minus_large_crank_angle_degrees']=math.degrees(motions[0].phase_rad-motions[1].phase_rad)%360
        report['results'][family]=result
        (ROOT/'outputs/slider_crank_k2.json').write_text(json.dumps(report,indent=2)+'\n')
        print(family,result.get('indicated_thermal_efficiency'),result.get('indicated_power_w'),result['convergence'],flush=True)


if __name__=='__main__':main()
