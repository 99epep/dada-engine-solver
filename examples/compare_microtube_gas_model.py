"""Frozen six-bar K2 comparison: legacy, variable transport, thermal entry."""
from dataclasses import replace, asdict
import json
import time
import numpy as np
from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.campaign.evaluator import json_values
from dada_solver.exchangers.gas_correlations import MicrotubeGasModel
from dada_solver.six_bar import IndependentSixBarVolumeKinematics,load_six_bar_mechanism
from compare_motor_motion_laws_stage7A5 import A5_CAMPAIGN,ROOT,_candidate_design_and_mass,_same_inventory_design,_evaluate
from optimize_motor_exchanger_asymmetry_stageK1 import _scaled_exchanger


def main():
    reference=json.loads((ROOT/'outputs/motor_champion_sixbar_k2.json').read_text())
    definition=CampaignDefinition(A5_CAMPAIGN)
    base,mass,*_=_candidate_design_and_mass(definition)
    v=base.configuration.machine_volumes;b=reference['thermodynamic_basis']
    motion=IndependentSixBarVolumeKinematics(
        load_six_bar_mechanism(ROOT/'outputs/small_sixbar_stage2f_r4_freeH_tightaxis.json',restart=0,second_branch=1),
        load_six_bar_mechanism(ROOT/'outputs/large_sixbar_stageL1.json'),v.small_cylinder,v.large_cylinder)
    base=_same_inventory_design(base,mass,motion)
    base=replace(base,heat_in=_scaled_exchanger(base.heat_in,b['stage_k2_k_i']),
                 heat_out=_scaled_exchanger(base.heat_out,b['stage_k2_k_o']))
    records=[]
    for label,settings in [('legacy',None),('variable_transport_developed',MicrotubeGasModel(thermal_entry=False)),
                            ('variable_transport_entry',MicrotubeGasModel())]:
        design=replace(base,heat_in=replace(base.heat_in,inputs=replace(base.heat_in.inputs,gas_model=settings)),
                       heat_out=replace(base.heat_out,inputs=replace(base.heat_out.inputs,gas_model=settings)))
        started=time.monotonic()
        try:
            result,state=_evaluate(label,design,definition,initial_state=np.array(reference['last_complete_state']))
        except (ValueError,RuntimeError) as error:
            result={'status':'model_domain_or_integration_failure','message':str(error)};state=None
        record=dict(label=label,gas_model=asdict(settings) if settings else None,
                    elapsed_seconds=time.monotonic()-started,result=result,
                    last_complete_state=state.tolist() if state is not None else None)
        records.append(record)
        report=dict(reference='motor_champion_sixbar_k2.json',configuration=json_values(asdict(base.configuration)),
                    wall_numerical_settings=asdict(definition.wall_numerical_settings),
                    heat_in=json_values(asdict(base.heat_in)),heat_out=json_values(asdict(base.heat_out)),
                    physical_changes='Internal transport and thermal entry only; external-air model unchanged',
                    mechanical_losses='unknown',records=records)
        (ROOT/'outputs/microtube_gas_model_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
        print(label,json.dumps({k:val for k,val in result.items() if k in ['status','message','indicated_thermal_efficiency','indicated_power_w','convergence']}),flush=True)


if __name__=='__main__': main()
