import json,math,tomllib
from pathlib import Path
import numpy as np
from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_model
from dada_solver.exchangers.hardware import HardwareInputs,connect_hardware
from dada_solver.exchangers.microtube_geometry import MicrotubeBank
for name,config_path in [('motor_doubled','examples/motor_demonstrator_four_bar.toml'),('motor_lower_lambda','examples/motor_demonstrator_lower_lambda_trial.toml')]:
 r=json.loads(Path('outputs/'+name+'_screening.json').read_text());d=r['hardware_inputs'];b=MicrotubeBank(**d['geometry']);c=load_simulation_configuration(config_path)
 w,_=connect_hardware(build_model(c),b,b,HardwareInputs(**d['properties'],**d['heat_in']),HardwareInputs(**d['properties'],**d['heat_out']),heat_in_valve_cda_m2=c.hydraulics.cold_to_large_valve_cda,heat_out_valve_cda_m2=c.hydraulics.hot_to_small_valve_cda)
 state=np.asarray(r['final_state'])
 angles,y=w.integrate_cycle(state,rtol=1e-9,atol=np.array([1e-14,1e-9]*4+[1e-8]*2+[1e-9]*5),maximum_step_angle=math.pi/720)
 error=float(np.max(np.abs(y[:10,-1]-state)/(1e-12+1e-6*np.maximum(np.abs(y[:10,-1]),np.abs(state)))))
 check=dict(method='one_unaccelerated_refined_cycle_from_reference_state',scaled_periodic_error=error,indicated_power_w=float(y[14,-1]*2),power_difference_w=float(y[14,-1]*2-r['last_cycle_indicated_power_w']))
 Path('outputs/'+name+'_refinement.json').write_text(json.dumps(check,indent=2));print(name,check,flush=True)
