"""Reproducible four-bar geometry/thermal/hydraulic/fan screening.

Run: PYTHONPATH=src python3 examples/motor_hardware_screening.py
"""
from pathlib import Path
import argparse
import csv
import json
import math
import time as timer
import numpy as np

from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_model, build_initial_state, initial_valve_topology
from dada_solver.state import ThermodynamicState
from dada_solver.exchangers.hardware import connect_hardware, load_hardware_definition
from dada_solver.exchangers.gas_diagnostics import cycle_microtube_diagnostics
from dada_solver.exchangers.duty import summarize_port
from dada_solver.exchangers.wall_cycle import solve_periodic_wall_motor, WallCycleNumericalSettings

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--hardware', type=Path, default=Path('examples/motor_hardware.toml'))
parser.add_argument('--configuration', type=Path, default=Path('examples/motor_demonstrator_four_bar.toml'))
parser.add_argument('--output-prefix', type=Path, default=Path('outputs/motor_hardware'))
parser.add_argument('--checkpoint-every', type=int, default=10)
parser.add_argument('--exclude-external-losses', action='store_true')
parser.add_argument('--warm-start', type=Path)
parser.add_argument('--accelerate-walls', action='store_true')
parser.add_argument('--plots', action='store_true')
parser.add_argument('--restart', type=Path, help='Use a previous report as a thermal initial guess')
parser.add_argument('--maximum-cycles', type=int, default=250)
args = parser.parse_args()
if args.restart and args.warm_start:
    parser.error('Choose restart or warm start, not both')
if args.maximum_cycles < 1 or args.checkpoint_every < 1:
    parser.error('Maximum cycles must be positive')
ROOT = Path(__file__).resolve().parents[1]
args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
def output_path(suffix):
    return Path(str(args.output_prefix)+suffix)
run_start=timer.perf_counter()
checkpoint_seconds=0.0
config = load_simulation_configuration(args.configuration)
data, bank, hi, ho = load_hardware_definition(args.hardware.read_text(), config.gas.heat_capacity_cp)
base = build_model(config)
wrapper, hardware = connect_hardware(base, bank, bank, hi, ho,
    heat_in_valve_cda_m2=config.hydraulics.cold_to_large_valve_cda,
    heat_out_valve_cda_m2=config.hydraulics.hot_to_small_valve_cda)
# Fixed gas inventory from the original seed; distribute it over the new volumes.
initial = build_initial_state(config, base)
volumes = wrapper.model.volumes(0)
from dada_solver.state import UniformCharge
pressure = initial.total_mass*config.gas.gas_constant*config.charge.temperature/volumes.total
state = np.r_[UniformCharge(pressure, config.charge.temperature).create_state(config.gas,volumes).as_array(),
              wrapper.heat_in.wall_capacity_j_k*hi.air_inlet_temperature_k,
              wrapper.heat_out.wall_capacity_j_k*ho.air_inlet_temperature_k]
if args.restart:
    previous = json.loads(args.restart.read_text())
    if previous['hardware_inputs'] != data:
        raise ValueError('Restart hardware inputs do not match')
    state = np.asarray(previous['final_state'], dtype=float)
reference_mass=initial.total_mass
if args.warm_start:
    previous=json.loads(args.warm_start.read_text())
    state=np.asarray(previous['final_state'],dtype=float)
    state[8] *= wrapper.heat_in.wall_capacity_j_k/previous['hardware']['H_i']['wall_capacity_j_k']
    state[9] *= wrapper.heat_out.wall_capacity_j_k/previous['hardware']['H_o']['wall_capacity_j_k']
    reference_mass=float(state[:8:2].sum())
def completed_cycle(cycle, end, error, history_item):
    global checkpoint_seconds
    history_item['scaled_periodic_error'] = history_item.pop('normalized_state_error')
    if cycle == 1 or cycle % 10 == 0 or error <= 1:
        print(history_item,flush=True)
    if cycle % args.checkpoint_every == 0 or error <= 1:
        checkpoint_start=timer.perf_counter()
        output_path('_checkpoint.json').write_text(json.dumps(dict(
            hardware_inputs=data, final_state=end.tolist(), cycle=cycle,
            scaled_periodic_error=error, status='initial_guess_only')))
        checkpoint_seconds += timer.perf_counter()-checkpoint_start
periodic = solve_periodic_wall_motor(wrapper, state, maximum_cycles=args.maximum_cycles,
    settings=WallCycleNumericalSettings(accelerate_walls=args.accelerate_walls),
    cycle_callback=completed_cycle)
if periodic.trajectory is None:
    raise RuntimeError(periodic.message)
angles, trajectory = periodic.angles, periodic.trajectory
state = periodic.last_complete_state
history_log = list(periodic.history)
converged = periodic.converged
stored = lambda x: float(np.sum(x[1:8:2])+np.sum(x[8:10]))
start = trajectory[:10,0]; end = trajectory[:10,-1]
residual = stored(end)-stored(start)-trajectory[10:12,-1].sum()+trajectory[14,-1]
frequency = config.angular_speed/(-2*math.pi)
power = float(trajectory[14,-1]*frequency)
qi = float(trajectory[10,-1]*frequency)
qo = float(trajectory[11,-1]*frequency)
fans = 0.0 if args.exclude_external_losses else hardware['total_fan_electrical_power_w']
max_reynolds=0.0
max_mach=0.0
port_history=[]
for angle, values in zip(angles, trajectory.T):
    gas_state=ThermodynamicState.from_array(values[:8])
    temperatures, pressures = gas_state.temperatures_and_pressures(config.gas,wrapper.model.volumes(float(angle)))
    flows=wrapper.model.evaluate(float(angle),gas_state,initial_valve_topology()).flows
    port_history.append([flows.small_to_cold, flows.cold_to_large, flows.large_to_hot, flows.hot_to_small])
    for flow, index in ((flows.small_to_cold,2),(flows.cold_to_large,2),(flows.large_to_hot,3),(flows.hot_to_small,3)):
        area=bank.dimensions()['tube_flow_area_m2']
        reynolds=abs(flow)*bank.inner_diameter_m/(area*hi.gas_viscosity_pa_s)
        density=pressures[index]/(config.gas.gas_constant*temperatures[index])
        mach=abs(flow)/(density*area*math.sqrt(config.gas.heat_capacity_cp/config.gas.heat_capacity_cv*config.gas.gas_constant*temperatures[index]))
        max_reynolds=max(max_reynolds,reynolds)
        max_mach=max(max_mach,mach)
gas_domains = cycle_microtube_diagnostics(wrapper, angles, trajectory)
if gas_domains is not None:
    passages = gas_domains['passages'].values()
    max_reynolds = max(p['hydraulic_upstream_ranges']['reynolds']['maximum'] for p in passages)
    max_mach = max(p['hydraulic_upstream_ranges']['mach']['maximum'] for p in passages)
time=angles/wrapper.model.angular_speed
port_history=np.asarray(port_history)
port_names=['S_to_H_i','H_i_to_L','L_to_H_o','H_o_to_S']
with output_path('_ports.csv').open('w',newline='') as stream:
    writer=csv.writer(stream)
    writer.writerow(['time_s']+[name+'_mass_flow_kg_s' for name in port_names])
    writer.writerows(zip(time,*port_history.T,strict=True))
output_path('_trajectory.npz').parent.mkdir(parents=True,exist_ok=True)
np.savez_compressed(output_path('_trajectory.npz'),angles=angles,states=trajectory,ports=port_history)
working_temperatures=trajectory[1:8:2]/(trajectory[:8:2]*config.gas.heat_capacity_cv)
mean_temperature_differences=dict(
    heat_in=float(np.trapz(hi.air_inlet_temperature_k-working_temperatures[2],time)/(time[-1]-time[0])),
    heat_out=float(np.trapz(working_temperatures[3]-ho.air_inlet_temperature_k,time)/(time[-1]-time[0])))
air_capacity_diagnostics={}
sampled_thermal = [wrapper.thermal_rates(float(angle), values) for angle, values in zip(angles, trajectory.T)]
for name, exchanger, columns, gas_index, wall_index in (
    ('heat_in',wrapper.heat_in,(0,1),2,8),
    ('heat_out',wrapper.heat_out,(2,3),3,9)):
    peak=float(np.max(np.abs(port_history[:,columns])))
    capacity=exchanger.air_mass_flow_kg_s*exchanger.air_cp_j_kg_k
    outlet_changes=[abs(rates[gas_index-2]['air_heat_w'])/capacity for rates in sampled_thermal]
    air_capacity_diagnostics[name]=dict(
        peak_internal_mass_flow_kg_s=peak,
        external_to_peak_internal_capacity_rate_ratio=capacity/(peak*config.gas.heat_capacity_cp) if peak else None,
        maximum_external_air_temperature_change_k=max(outlet_changes),
        scope='sampled_cycle; fixed_external_flow; finite_film_resistance_retained')
report=dict(mean_inlet_air_gas_temperature_differences_k=mean_temperature_differences,
    external_air_capacity_diagnostics=air_capacity_diagnostics,
    wall_time_seconds=timer.perf_counter()-run_start,checkpoint_write_seconds=checkpoint_seconds,
    external_air_losses='excluded_by_user' if args.exclude_external_losses else 'estimated',
    accounted_fan_power_w=fans,
    warm_start_source=str(args.warm_start) if args.warm_start else None,
    port_duty={name:summarize_port(time,port_history[:,i]) for i,name in enumerate(port_names)},
    periodic_status='converged' if converged else 'maximum_cycles',
    hardware_inputs=data, hardware=hardware, convergence_history=history_log,
    restart_source=str(args.restart) if args.restart else None,
    frequency_hz=frequency, final_state=state.tolist(), total_energy_residual_j=float(residual),
    total_mass_drift_kg=float(state[:8:2].sum()-reference_mass),
    configuration_source=args.configuration.read_text(),
    last_cycle_indicated_power_w=power,last_cycle_external_heat_input_w=qi,
    last_cycle_external_heat_out_w=qo,
    fan_electrical_power_w=None if args.exclude_external_losses else fans,
    indicated_power_minus_fans_w=power-fans if converged and not args.exclude_external_losses else None,
    power_after_accounted_external_losses_w=power-fans if converged else None,
    indicated_thermal_efficiency=power/qi if converged and power>0 and qi>0 and qo<0 else None,
    useful_shaft_power_w=None,
    maximum_tube_reynolds=max_reynolds, maximum_tube_mach_at_exchanger_state=max_mach,
    microtube_gas_domains=gas_domains,
    applicability=('variable_transport_and_entry; pulse_and_distribution_unvalidated' if gas_domains is not None
                   else 'screening_only; constant_Nusselt; entrance_pulse_and_distribution_unvalidated'),
    laminar_reynolds_screen_passed=bool(max_reynolds<2300),
    mach_screen_passed=bool(max_mach <= config.validity.maximum_mach_number),
    mechanical_losses='unavailable', calibration='Doty_geometry_reference; no_empirical_UA_fit')
output=output_path('_screening.json')
output.write_text(json.dumps(report,indent=2))
print(json.dumps({k:v for k,v in report.items() if k not in ('hardware','hardware_inputs','convergence_history','final_state','configuration_source','port_duty')},indent=2))

if args.plots:
    from dada_solver.exchangers.trial_plots import save_trial_plots
    sampled_events=save_trial_plots(wrapper,angles,trajectory,port_history,args.output_prefix,
        external_losses_excluded=args.exclude_external_losses)
    report['sampled_valve_events']=sampled_events
    output.write_text(json.dumps(report,indent=2))
