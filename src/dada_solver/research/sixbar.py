"""Typed fixed-pair input bundle and five-coordinate study policy.

The physical mechanism, exchanger correlations and cycle solver remain in their
existing production modules. This adapter never imports historical examples.
"""
from dataclasses import dataclass, replace
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from dada_solver import configuration as cfg
from dada_solver.campaign.adapters import PreflightRejection
from dada_solver.campaign.candidate import canonical_json
from dada_solver.exchangers.gas_correlations import MicrotubeGasModel, GasSurfaceAccommodation, SecondOrderSlip
from dada_solver.exchangers.gas_transport import DiluteGasTransport
from dada_solver.exchangers.hardware import HardwareInputs
from dada_solver.exchangers.microtube import MicrotubeExchanger
from dada_solver.exchangers.microtube_geometry import MicrotubeBank
from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.four_bar import SharedCrankRockerDesign
from dada_solver.geometry import CylinderVolumeLimits, MachineVolumes
from dada_solver.machine import MachineDesign
from dada_solver.six_bar import SixBarCylinderMechanism, IndependentSixBarVolumeKinematics

PARAMETERS = {
    'volume.swept_ratio': ('continuous', '1', 'swept_ratio'),
    'microtube.heat_in.tube_count': ('integer', '1', 'n_i'),
    'microtube.heat_in.tube_length_m': ('continuous', 'm', 'length_i_m'),
    'microtube.heat_out.tube_count': ('integer', '1', 'n_o'),
    'microtube.heat_out.tube_length_m': ('continuous', 'm', 'length_o_m'),
}


def configuration_from_data(data):
    raw = dict(data)
    for name in ('hydraulic_flow_models', 'humidity_screening', 'free_kinematics', 'shared_coupler_projection_design'):
        if raw.get(name) is not None:
            raise ValueError(f'Fixed six-bar basis does not support {name}.')
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
    family=raw.pop('family','microtube_air_wall')
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


@dataclass(frozen=True)
class SixBarStudyBasis:
    source: str
    sha256: str
    configuration: cfg.SimulationConfiguration
    small: SixBarCylinderMechanism
    large: SixBarCylinderMechanism
    heat_in: MicrotubeExchanger
    heat_out: MicrotubeExchanger

    @property
    def data(self):
        return json.loads(self.source)


def load_basis(path, expected_sha256):
    source = Path(path).read_bytes()
    digest = hashlib.sha256(source).hexdigest()
    if digest != expected_sha256:
        raise ValueError(f'Basis SHA-256 mismatch: {path}. Restore the source or start a new study.')
    data = json.loads(source)
    canonical_json(data)  # Reject NaN/Infinity even though json.loads accepts them.
    required = {'schema_version', 'configuration', 'small', 'large', 'heat_in', 'heat_out',
                'geometry', 'warm_start', 'wall_settings', 'reference_parameters', 'reference_result', 'provenance'}
    if set(data) != required or type(data['schema_version']) is not int or data['schema_version'] != 1:
        raise ValueError('Unsupported or incomplete fixed-six-bar basis schema.')
    geometry_keys = {'total_swept_m3','small_clearance_ratio','large_clearance_ratio',
                     'seed_ni','seed_no','seed_li','seed_lo','seed_hi_cda','seed_ho_cda'}
    if set(data['geometry']) != geometry_keys or any(
            isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) or v <= 0
            for v in data['geometry'].values()):
        raise ValueError('Basis geometry and valve-scaling references must be explicit, positive and finite.')
    if any(type(data['geometry'][k]) is not int for k in ('seed_ni','seed_no')):
        raise ValueError('Basis reference tube counts must be integers.')
    for side in ('small','large'):
        if any(type(data[side][k]) is not int for k in ('primary_branch','second_branch')):
            raise ValueError('Six-bar assembly branches must be integer -1 or +1.')
    if data['provenance'].get('length_unit') != 'crank_radius' or data['provenance'].get('angle_unit') != 'rad':
        raise ValueError('Basis mechanism units must be crank_radius and rad.')
    warm = data['warm_start']
    state = np.asarray(warm['values'], dtype=float)
    caps = np.asarray(warm['wall_capacities_j_k'], dtype=float)
    if state.shape != (10,) or caps.shape != (2,) or not np.all(np.isfinite(state)) or np.any(state <= 0) or not np.all(np.isfinite(caps)) or np.any(caps <= 0):
        raise ValueError('Basis warm state requires ten positive finite states and two wall capacities.')
    configuration = configuration_from_data(data['configuration'])
    if not configuration.motor_operation or configuration.charge.total_mass is None:
        raise ValueError('The fixed-pair protocol requires motor direction and fixed gas inventory.')
    if not np.isclose(state[:8:2].sum(), configuration.charge.total_mass, rtol=1e-10, atol=0):
        raise ValueError('Warm-state gas inventory differs from the fixed basis inventory.')
    try:
        small = SixBarCylinderMechanism(**data['small'])
        large = SixBarCylinderMechanism(**data['large'])
    except (ValueError, TypeError, ArithmeticError) as error:
        raise ValueError(f'Invalid fixed six-bar geometry: {error}') from error
    return SixBarStudyBasis(source.decode(), digest, configuration, small, large,
                            exchanger_from_data(data['heat_in']), exchanger_from_data(data['heat_out']))


class SixBarStudyAdapter:
    def __init__(self, basis, configuration=None):
        self.basis = basis
        self.configuration = configuration or basis.configuration

    def build(self, physical):
        if set(physical) != set(PARAMETERS):
            raise PreflightRejection('invalid_parameterization', 'The fixed-pair protocol requires exactly its five owned coordinates.')
        p = {PARAMETERS[k][2]: v for k, v in physical.items()}
        if any(not np.isfinite(v) or v <= 0 for v in p.values()) or any(type(p[k]) is not int for k in ('n_i', 'n_o')):
            raise PreflightRejection('invalid_parameterization', 'Positive finite coordinates and integer tube counts are required.')
        b = self.basis
        g = b.data['geometry']
        # Preserve the reference arithmetic, including subtraction for S.
        large_swept = g['total_swept_m3'] / (1.0 + p['swept_ratio'])
        small_swept = g['total_swept_m3'] - large_swept
        smin = g['small_clearance_ratio'] * small_swept
        lmin = g['large_clearance_ratio'] * large_swept
        small = CylinderVolumeLimits(smin, smin + small_swept)
        large = CylinderVolumeLimits(lmin, lmin + large_swept)
        if large.maximum > 0.066:
            raise PreflightRejection('invalid_parameterization', 'Large-cylinder maximum enclosed volume exceeds 0.066 m^3.')
        config = replace(self.configuration, machine_volumes=replace(self.configuration.machine_volumes,
                         small_cylinder=small, large_cylinder=large))
        try:
            hi = replace(b.heat_in, bank=replace(b.heat_in.bank, tube_count=p['n_i'], tube_length_m=p['length_i_m']),
                         outlet_valve_cda_m2=g['seed_hi_cda'] * p['n_i'] / g['seed_ni'])
            ho = replace(b.heat_out, bank=replace(b.heat_out.bank, tube_count=p['n_o'], tube_length_m=p['length_o_m']),
                         outlet_valve_cda_m2=g['seed_ho_cda'] * p['n_o'] / g['seed_no'])
            hi.build(); ho.build()
        except (ValueError, ArithmeticError) as error:
            raise PreflightRejection('invalid_exchanger', str(error)) from error
        return MachineDesign(config, heat_in=hi, heat_out=ho,
                             kinematics=IndependentSixBarVolumeKinematics(b.small, b.large, small, large))
