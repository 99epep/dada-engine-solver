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


@dataclass(frozen=True)
class SliderMotion:
    rod_over_crank: float
    offset_over_crank: float
    phase_rad: float
    volume_increases_with_coordinate: bool

    def __post_init__(self):
        if not all(math.isfinite(v) for v in (self.rod_over_crank,self.offset_over_crank,self.phase_rad)) or self.rod_over_crank<=1+abs(self.offset_over_crank):
            raise ValueError('Invalid slider-crank closure')

    def normalized(self, theta):
        l,e=self.rod_over_crank,self.offset_over_crank
        a=theta+self.phase_rad
        z=e-math.sin(a); root=math.sqrt(l*l-z*z)
        x=math.cos(a)+root; dx=-math.sin(a)+z*math.cos(a)/root
        low=math.sqrt((l-1)**2-e*e); high=math.sqrt((l+1)**2-e*e)
        q=(x-low)/(high-low); dq=dx/(high-low)
        return (q,dq) if self.volume_increases_with_coordinate else (1-q,-dq)


@dataclass(frozen=True)
class SliderCrankKinematics:
    small: SliderMotion
    large: SliderMotion
    small_volume_limits: CylinderVolumeLimits
    large_volume_limits: CylinderVolumeLimits
    small_physical_stroke = None
    large_physical_stroke = None

    def breakpoint_angles(self): return ()
    def _side(self, side, theta):
        q,dq=getattr(self,side).normalized(theta)
        limits=getattr(self,side+'_volume_limits')
        return limits.minimum+limits.swept*q,limits.swept*dq
    def small_cylinder_volume(self,theta): return self._side('small',theta)[0]
    def large_cylinder_volume(self,theta): return self._side('large',theta)[0]
    def small_cylinder_volume_derivative(self,theta): return self._side('small',theta)[1]
    def large_cylinder_volume_derivative(self,theta): return self._side('large',theta)[1]
    def cylinder_volumes_and_derivatives(self,theta):
        s,ds=self._side('small',theta);l,dl=self._side('large',theta)
        return s,l,ds,dl


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
