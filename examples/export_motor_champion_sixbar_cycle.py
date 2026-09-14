"""Export the reference six-bar K2 cycle for illustration (English labels, SI units)."""
from dataclasses import replace
import csv
import hashlib
import json
import math
import numpy as np

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.dynamics import ValveTopology
from dada_solver.exchangers.wall_cycle import solve_periodic_wall_motor, wall_cycle_performance
from dada_solver.heat_transfer import PrescribedHeatRate
from dada_solver.six_bar import IndependentSixBarVolumeKinematics, load_six_bar_mechanism
from dada_solver.state import ThermodynamicState
from dada_solver.valves import ValveState
from compare_motor_motion_laws_stage7A5 import A5_CAMPAIGN, ROOT, _candidate_design_and_mass, _same_inventory_design
from optimize_motor_exchanger_asymmetry_stageK1 import _scaled_exchanger
from evaluate_motor_champion_sixbar import SMALL, LARGE


def main():
    source = ROOT/'outputs/motor_champion_sixbar_k2.json'
    reference = json.loads(source.read_text())
    for side, path in [('small', SMALL), ('large', LARGE)]:
        if hashlib.sha256(path.read_bytes()).hexdigest() != reference['sources'][side]['sha256']:
            raise ValueError('Synthesis geometry changed since the reference evaluation.')
    definition = CampaignDefinition(A5_CAMPAIGN)
    base, mass, _, _, _ = _candidate_design_and_mass(definition)
    limits = base.configuration.machine_volumes
    kin = IndependentSixBarVolumeKinematics(
        load_six_bar_mechanism(SMALL, restart=0, second_branch=1),
        load_six_bar_mechanism(LARGE), limits.small_cylinder, limits.large_cylinder)
    design = _same_inventory_design(base, mass, kin)
    basis = reference['thermodynamic_basis']
    design = replace(design,
        heat_in=_scaled_exchanger(design.heat_in, basis['stage_k2_k_i']),
        heat_out=_scaled_exchanger(design.heat_out, basis['stage_k2_k_o']))
    wrapper = design.build()
    periodic = solve_periodic_wall_motor(wrapper, reference['last_complete_state'],
        maximum_cycles=design.configuration.numerical.maximum_cycles,
        settings=definition.wall_numerical_settings)
    if not periodic.converged:
        raise RuntimeError('Export requires a converged physical cycle.')
    model = wrapper.model
    stored = ValveTopology(ValveState.CLOSED, ValveState.CLOSED)

    def row(angle, values):
        gas = ThermodynamicState.from_array(values[:8])
        volumes = model.volumes(angle)
        temperatures, pressures = gas.temperatures_and_pressures(model.gas, volumes)
        hi = wrapper.heat_in.rates(temperatures[2], values[8])
        ho = wrapper.heat_out.rates(temperatures[3], values[9])
        instantaneous = replace(model, cold_heat_transfer=PrescribedHeatRate(hi['gas_heat_w']),
                                hot_heat_transfer=PrescribedHeatRate(ho['gas_heat_w']))
        rates = instantaneous.evaluate(angle, gas, stored)
        topology = model.effective_topology(angle, gas, stored)
        out = dict(time_s=angle/model.angular_speed, cycle_fraction=angle/(2*math.pi),
                   motor_angle_deg=math.degrees(angle), study_angle_deg=-math.degrees(angle))
        for i, region in enumerate(('S', 'L', 'Hi', 'Ho')):
            out.update({f'{region}_volume_m3': volumes.as_tuple()[i],
                        f'{region}_pressure_Pa': pressures[i],
                        f'{region}_temperature_K': temperatures[i],
                        f'{region}_mass_kg': values[2*i],
                        f'{region}_internal_energy_J': values[2*i+1]})
        for region, index, exchanger, thermal in [('Hi',8,wrapper.heat_in,hi),('Ho',9,wrapper.heat_out,ho)]:
            out.update({f'{region}_wall_temperature_K': values[index]/exchanger.wall_capacity_j_k,
                        f'{region}_wall_energy_J': values[index],
                        f'{region}_wall_to_gas_heat_W': thermal['gas_heat_w'],
                        f'{region}_external_to_wall_heat_W': thermal['air_heat_w'],
                        f'{region}_air_inlet_temperature_K': exchanger.air_inlet_temperature_k,
                        f'{region}_air_outlet_temperature_K': thermal['air_outlet_temperature_k']})
        for name, field in [('S_to_Hi','small_to_cold'),('Hi_to_L','cold_to_large'),
                            ('L_to_Ho','large_to_hot'),('Ho_to_S','hot_to_small')]:
            out[f'{name}_mass_flow_kg_s'] = getattr(rates.flows,field)
        out['Hi_to_L_valve_open'] = int(topology.cold_to_large == ValveState.OPEN)
        out['Ho_to_S_valve_open'] = int(topology.hot_to_small == ValveState.OPEN)
        vs, vl = model.cylinder_volume_rates(angle)
        out.update(S_volume_rate_m3_s=vs, L_volume_rate_m3_s=vl,
                   S_indicated_power_W=pressures[0]*vs, L_indicated_power_W=pressures[1]*vl,
                   total_indicated_power_W=rates.gas_work_rate,
                   cumulative_external_Hi_heat_J=values[10], cumulative_external_Ho_heat_J=values[11],
                   cumulative_Hi_gas_heat_J=values[12], cumulative_Ho_gas_heat_J=values[13],
                   cumulative_indicated_work_J=values[14])
        return out

    files = []
    for suffix, angles, trajectory in [
        ('raw',periodic.angles,periodic.trajectory),
        ('',np.linspace(0,2*math.pi,1441),np.array([
            np.interp(np.linspace(0,2*math.pi,1441),periodic.angles,v) for v in periodic.trajectory]))]:
        path = ROOT/f'outputs/champion-6bar_thermodynamic_cycle{"_"+suffix if suffix else ""}.csv'
        rows = [row(float(a),v) for a,v in zip(angles,trajectory.T)]
        with path.open('w',newline='') as stream:
            writer = csv.DictWriter(stream,fieldnames=rows[0])
            writer.writeheader(); writer.writerows(rows)
        files.append(dict(file=path.name, rows=len(rows), columns=len(rows[0])))
    performance = wall_cycle_performance(wrapper,periodic.trajectory)
    metadata = dict(reference_file=source.name, reference_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                    hydraulic_loss_multiplier=1, files=files,
                    indicated_thermal_efficiency=performance.thermal_efficiency,
                    indicated_power_W=performance.gas_power, convergence_history=periodic.history,
                    period_s=2*math.pi/model.angular_speed,
                    sampling='Main CSV: 0.25 degree grid, linear interpolation of conservative states and quadratures; rates recomputed. Raw CSV: original adaptive solver samples.',
                    signs='Heat positive into gas/wall as named; indicated work positive on expansion; signed flows positive in column-name direction.',
                    angle_convention='Motor angle increases with time; study angle is its negative. No additional motion transformation.',
                    units='SI units encoded in column names. Valve open: 1, closed: 0.',
                    endpoint='Both cycle endpoints retained, without forcing the final state equal to the initial state.',
                    mechanical_losses='unknown')
    (ROOT/'outputs/champion-6bar_thermodynamic_cycle_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(metadata,indent=2))


if __name__ == '__main__':
    main()
