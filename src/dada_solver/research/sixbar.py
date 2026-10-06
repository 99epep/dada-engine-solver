"""Reconstruct explicit current machine configuration and exchanger inputs."""
from dada_solver import configuration as cfg
from dada_solver.exchangers.gas_correlations import MicrotubeGasModel, GasSurfaceAccommodation, SecondOrderSlip
from dada_solver.exchangers.gas_transport import DiluteGasTransport
from dada_solver.exchangers.hardware import HardwareInputs
from dada_solver.exchangers.microtube import MicrotubeExchanger
from dada_solver.exchangers.microtube_geometry import MicrotubeBank
from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.four_bar import SharedCrankRockerDesign
from dada_solver.geometry import CylinderVolumeLimits, MachineVolumes


def configuration_from_data(data):
    raw = dict(data)
    for name in ('hydraulic_flow_models', 'humidity_screening', 'free_kinematics', 'shared_coupler_projection_design'):
        if raw.get(name) is not None:
            raise ValueError(f'Machine basis does not support {name}.')
    volumes = dict(raw['machine_volumes'])
    for side in ('small_cylinder', 'large_cylinder'):
        volumes[side] = CylinderVolumeLimits(**volumes[side])
    raw['machine_volumes'] = MachineVolumes(**volumes)
    constructors = dict(gas=CaloricallyPerfectGas, charge=cfg.ChargeConfiguration,
                        hydraulics=cfg.HydraulicNetworkConfiguration,
                        hot_to_small_valve=cfg.ValveThresholdConfiguration,
                        cold_to_large_valve=cfg.ValveThresholdConfiguration,
                        validity=cfg.ValidityThresholds, numerical=cfg.NumericalConfiguration,
                        shared_four_bar_design=SharedCrankRockerDesign)
    for name, constructor in constructors.items():
        if name=='gas' and raw[name].get('model')=='tabulated_single_phase_rho_u':
            from dada_solver.tabulated_fluid import TabulatedFluid
            raw[name]=TabulatedFluid.from_data(raw[name])
        elif raw.get(name) is not None: raw[name] = constructor(**raw[name])
    return cfg.SimulationConfiguration(**raw)


def exchanger_from_data(data):
    raw = dict(data)
    family=raw.pop('family')
    if family not in ('microtube_air_wall','external_stream_wall'): raise ValueError('Unknown exchanger family.')
    inputs = dict(raw['inputs'])
    if inputs.get('gas_model') is not None:
        model = dict(inputs['gas_model'])
        model['transport'] = DiluteGasTransport(**model['transport'])
        model['accommodation'] = GasSurfaceAccommodation(**model['accommodation'])
        if model['slip'] is not None: model['slip'] = SecondOrderSlip(**model['slip'])
        inputs['gas_model'] = MicrotubeGasModel(**model)
    if family=='external_stream_wall':
        from dada_solver.exchangers.external_stream import ExternalFluidStream
        from dada_solver.exchangers.external_microtube import ExternalStreamHardwareInputs, ExternalStreamMicrotubeExchanger
        inputs['external_stream']=ExternalFluidStream(**inputs['external_stream'])
        raw['inputs']=ExternalStreamHardwareInputs(**inputs)
    else: raw['inputs'] = HardwareInputs(**inputs)
    raw['bank'] = MicrotubeBank(**raw['bank'])
    result = (ExternalStreamMicrotubeExchanger if family=='external_stream_wall' else MicrotubeExchanger)(**raw)
    result.build()
    return result
