"""Compare saved thermal trials and pressure differences during forward transfer."""
import json
from pathlib import Path
import tempfile
import numpy as np
from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_model
from dada_solver.exchangers.hardware import HardwareInputs, connect_hardware
from dada_solver.exchangers.microtube_geometry import MicrotubeBank
from dada_solver.state import ThermodynamicState

results={}
for name in ('motor_325c','motor_original_325c','motor_parallel_325c'):
    report=json.loads(Path('outputs/'+name+'_screening.json').read_text())
    with tempfile.NamedTemporaryFile(mode='w',suffix='.toml') as stream:
        stream.write(report['configuration_source']);stream.flush()
        config=load_simulation_configuration(Path(stream.name))
    data=report['hardware_inputs'];bank=MicrotubeBank(**data['geometry'])
    wrapper,_=connect_hardware(build_model(config),bank,bank,
        HardwareInputs(**data['properties'],**data['heat_in']),
        HardwareInputs(**data['properties'],**data['heat_out']),
        heat_in_valve_cda_m2=config.hydraulics.cold_to_large_valve_cda,
        heat_out_valve_cda_m2=config.hydraulics.hot_to_small_valve_cda)
    trajectory=np.load('outputs/'+name+'_trajectory.npz')
    pressures=[]
    for angle,state in zip(trajectory['angles'],trajectory['states'].T):
        _,p=ThermodynamicState.from_array(state[:8]).temperatures_and_pressures(config.gas,wrapper.model.volumes(float(angle)))
        pressures.append(p)
    p=np.asarray(pressures);flows=trajectory['ports'];drops={}
    for label,up,mid,down,cols in [('H_i',0,2,1,(0,1)),('H_o',1,3,0,(2,3))]:
        active=(flows[:,cols[0]]>1e-8)&(flows[:,cols[1]]>1e-8)
        drops[label]=dict(maximum_forward_path_pressure_drop_pa=float(np.max((p[:,up]-p[:,down])[active])),
            maximum_inlet_half_pressure_drop_pa=float(np.max((p[:,up]-p[:,mid])[active])),
            maximum_outlet_half_including_valve_pressure_drop_pa=float(np.max((p[:,mid]-p[:,down])[active])),
            scope='both adjacent ports flowing forward above 1e-8 kg/s; includes tube, header and outlet valve; individual maxima need not coincide')
    results[name]=dict(indicated_power_w=report['last_cycle_indicated_power_w'],heat_input_w=report['last_cycle_external_heat_input_w'],
        indicated_efficiency=report['indicated_thermal_efficiency'],pressure_diagnostics=drops,
        geometry=bank.dimensions(),overall_static_conductance_w_k=report['hardware']['H_i']['overall_static_conductance_w_k'],
        external_air_capacity_diagnostics=report['external_air_capacity_diagnostics'])
Path('outputs/motor_parallel_comparison.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2))
