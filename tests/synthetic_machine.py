"""Small current-contract inputs, without application or campaign results."""
from dataclasses import asdict, replace
from pathlib import Path
import tomllib

from dada_solver.configuration import load_simulation_configuration
from dada_solver.free_kinematics import FreeKinematicsConfiguration, FreeMotionDefinition
from dada_solver.four_bar import SharedCrankRockerDesign
from dada_solver.exchangers.hardware import HardwareInputs
from dada_solver.exchangers.microtube_geometry import MicrotubeBank
from dada_solver.research.study_io import dumps

CONFIGURATION = Path(__file__).parent / 'data' / 'sizing_machine.toml'


def configuration(family='harmonic_example', *, motor=False):
    base = load_simulation_configuration(CONFIGURATION)
    base = replace(base, angular_speed=-10. if motor else 10.,
                   valve_model='continuous_ideal_diode')
    if family == 'free':
        motions = [FreeMotionDefinition((1., 0., -1., 0.), v.minimum, v.maximum)
                   for v in (base.machine_volumes.small_cylinder,
                             base.machine_volumes.large_cylinder)]
        return replace(base, kinematics_type='free',
                       free_kinematics=FreeKinematicsConfiguration(*motions))
    if family == 'shared_crank_rocker':
        return replace(base, kinematics_type=family, four_bar_ground_distance=.1,
                       four_bar_crank_angle_offset_degrees=0., four_bar_crank_direction=1,
                       shared_four_bar_design=SharedCrankRockerDesign(
                           .25, .90, .85, .72, -.18, 1.05, .95, .81, .22))
    if family in ('published_e0_opposed', 'published_f65_opposed'):
        return replace(base, kinematics_type=family, four_bar_ground_distance=.1,
                       four_bar_connecting_rod_ratio=5.,
                       four_bar_crank_angle_offset_degrees=0., four_bar_crank_direction=1)
    if family == 'ideal_piecewise_linear':
        return replace(base, kinematics_type=family, small_lambda_target=.6,
                       large_lambda_target=.6, adiabatic_sector_fraction=.2)
    return base


def configuration_data():
    return tomllib.loads(CONFIGURATION.read_text())


def hardware():
    bank = MicrotubeBank(256, .2, .001, .0001, .0015, .002)
    inputs = HardwareInputs(
        metal_conductivity_w_m_k=15., metal_density_kg_m3=8000.,
        metal_cp_j_kg_k=500., extra_wall_capacity_j_k=2.,
        gas_conductivity_w_m_k=.026, gas_viscosity_pa_s=1.8e-5, gas_nusselt=3.66,
        air_conductivity_w_m_k=.026, air_viscosity_pa_s=1.8e-5, air_nusselt=3.66,
        air_density_kg_m3=1.2, air_cp_j_kg_k=1005., air_mass_flow_kg_s=.02,
        air_inlet_temperature_k=400., air_poiseuille_number=64.,
        air_minor_loss_coefficient=1., fan_total_efficiency=.5,
        core_loss_multiplier=1., header_loss_coefficient=0.)
    return bank, inputs


def hardware_text():
    bank, inputs = hardware()
    props = asdict(inputs)
    props.pop('gas_model')
    inlet = props.pop('air_inlet_temperature_k')
    geometry = {k: v for k, v in asdict(bank).items() if v is not None}
    return dumps(dict(geometry=geometry, properties=props,
                      heat_in=dict(air_inlet_temperature_k=inlet),
                      heat_out=dict(air_inlet_temperature_k=300.)))


def campaign_file(directory, *, microtube=False):
    """Write only the inputs required by the campaign loader/identity tests."""
    directory.mkdir(parents=True, exist_ok=True)
    data = configuration_data()
    data['operation']['angular_speed'] = -10.
    data['valves']['model'] = 'continuous_ideal_diode'
    data['kinematics'] = dict(
        type='shared_crank_rocker', ground_distance=.1, crank_ratio=.25,
        crank_angle_offset_degrees=0., crank_direction=1,
        small_four_bar=dict(coupler_ratio=.90, rocker_ratio=.85,
                            output_along_ratio=.72, output_normal_ratio=-.18),
        large_four_bar=dict(coupler_ratio=1.05, rocker_ratio=.95,
                            output_along_ratio=.81, output_normal_ratio=.22))
    (directory / 'base.toml').write_text(dumps(data))
    raw = dict(campaign=dict(base_configuration='base.toml'),
               families=dict(kinematics='free', exchanger='microtube' if microtube else 'reservoir'),
               free=dict(small_coordinates=[.2, -.3], large_coordinates=[-.1, .4]),
               parameters=[dict(name='free.small.0', lower=-.5, upper=.5, initial=.2),
                           dict(name='free.large.0', lower=-.5, upper=.5, initial=-.1)],
               objective=dict(type='maximize_thermal_efficiency'),
               constraints=[dict(type='valid_thermodynamic_model')],
               search=dict(type='sobol', seed=7, scramble=True))
    if microtube:
        (directory / 'hardware.toml').write_text(hardware_text())
        raw['campaign']['hardware_configuration'] = 'hardware.toml'
        raw['parameters'].append(dict(name='microtube.heat_in.tube_length_m',
                                      lower=.1, upper=.3, initial=.2))
        raw['wall_numerical'] = dict(integration_method='LSODA',
                                    integration_absolute_tolerances=[1e-9] * 15,
                                    accelerate_walls=True)
    path = directory / 'campaign.toml'
    path.write_text(dumps(raw))
    return path


def connect_microtube_machine(model, heat_in_bank, heat_out_bank, heat_in_inputs,
                              heat_out_inputs, *, heat_in_valve_cda_m2, heat_out_valve_cda_m2):
    from dada_solver.exchangers.microtube import MicrotubeExchanger
    from dada_solver.exchangers.base import connect_exchangers
    incoming = MicrotubeExchanger(heat_in_bank, heat_in_inputs, heat_in_valve_cda_m2,
                                 model.heat_in_valve_placement)
    outgoing = MicrotubeExchanger(heat_out_bank, heat_out_inputs, heat_out_valve_cda_m2,
                                 model.heat_out_valve_placement)
    reports = {name: dict(exchanger.build().metadata)
               for name, exchanger in (('H_i', incoming), ('H_o', outgoing))}
    reports['total_fan_electrical_power_w'] = sum(
        reports[name]['fan_electrical_power_w'] for name in ('H_i', 'H_o'))
    return connect_exchangers(model, incoming, outgoing), reports
