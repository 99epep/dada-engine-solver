"""Portable machine inputs independent of the selected kinematic family."""
from dataclasses import dataclass, replace
import hashlib
import json
from pathlib import Path
import math
import numpy as np
from dada_solver.campaign.candidate import canonical_json
from dada_solver.campaign.adapters import PreflightRejection
from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.machine import MachineDesign
from .sixbar import configuration_from_data, exchanger_from_data
from .families import ParameterSpec


@dataclass(frozen=True)
class MachineBasis:
    source: str
    sha256: str
    configuration: object
    heat_in: object
    heat_out: object

    @property
    def data(self): return json.loads(self.source)


def load_machine_basis(path, expected_sha256):
    source=Path(path).read_bytes()
    digest=hashlib.sha256(source).hexdigest()
    if digest!=expected_sha256: raise ValueError('Machine basis SHA-256 mismatch.')
    data=json.loads(source); canonical_json(data)
    required={'schema_version','configuration','heat_in','heat_out','geometry','wall_settings','warm_start','provenance'}
    if set(data)!=required or type(data['schema_version']) is not int or data['schema_version'] not in (2,3):
        raise ValueError('Unsupported machine basis schema.')
    if data['schema_version']==2 and (data['configuration']['gas'].get('model') or any((data[side] or {}).get('family')=='external_stream_wall' for side in ('heat_in','heat_out'))):
        raise ValueError('New fluid/exchanger definitions require machine schema 3.')
    configuration=configuration_from_data(data['configuration'])
    if (data['heat_in'] is None)!=(data['heat_out'] is None): raise ValueError('Provide both exchangers or neither.')
    hi=exchanger_from_data(data['heat_in']) if data['heat_in'] is not None else None
    ho=exchanger_from_data(data['heat_out']) if data['heat_out'] is not None else None
    if configuration.charge.total_mass is None: raise ValueError('V2 machine studies currently require explicit gas inventory.')
    if data['warm_start'] is not None:
        warm=data['warm_start']; values=np.asarray(warm['values'],float)
        if hi is None or values.shape!=(10,) or np.any(values<=0) or not np.all(np.isfinite(values)):
            raise ValueError('Source warm state requires ten positive finite wall-model states.')
        caps=np.asarray(warm['wall_capacities_j_k'],float)
        if caps.shape!=(2,) or np.any(caps<=0) or not np.all(np.isfinite(caps)): raise ValueError('Invalid source wall capacities.')
        if not np.isclose(values[:8:2].sum(),configuration.charge.total_mass,rtol=1e-10,atol=0): raise ValueError('Source warm state inventory mismatch.')
    return MachineBasis(source.decode(),digest,configuration,hi,ho)


MACHINE_SPECS={
    'volume.swept_ratio':ParameterSpec(positive=True),
    'volume.total_swept_m3':ParameterSpec('m^3',positive=True),
    'volume.small_clearance_ratio':ParameterSpec(positive=True),
    'volume.large_clearance_ratio':ParameterSpec(positive=True),
    'operation.frequency_hz':ParameterSpec('Hz',positive=True),
    'charge.total_mass_kg':ParameterSpec('kg',positive=True)}
BANK_SPECS={'tube_count':ParameterSpec(kind='integer',positive=True),
    **{k:ParameterSpec('m',positive=True) for k in ('tube_length_m','inner_diameter_m','wall_thickness_m','pitch_m','header_depth_m')}}
INPUT_SPECS={'air_inlet_temperature_k':ParameterSpec('K',positive=True),
    'air_mass_flow_kg_s':ParameterSpec('kg/s',positive=True),
    'metal_conductivity_w_m_k':ParameterSpec('W/(m*K)',positive=True),
    'metal_density_kg_m3':ParameterSpec('kg/m^3',positive=True),
    'metal_cp_j_kg_k':ParameterSpec('J/(kg*K)',positive=True)}


def machine_parameters(basis):
    c=basis.configuration; v=c.machine_volumes
    defaults={'volume.swept_ratio':v.small_cylinder.swept/v.large_cylinder.swept,
        'volume.total_swept_m3':v.small_cylinder.swept+v.large_cylinder.swept,
        'volume.small_clearance_ratio':v.small_cylinder.minimum/v.small_cylinder.swept,
        'volume.large_clearance_ratio':v.large_cylinder.minimum/v.large_cylinder.swept,
        'operation.frequency_hz':abs(c.angular_speed)/(2*math.pi),'charge.total_mass_kg':c.charge.total_mass}
    specs=dict(MACHINE_SPECS)
    for side in ('heat_in','heat_out'):
        exchanger=getattr(basis,side)
        if exchanger is None: continue
        for name,spec in BANK_SPECS.items():
            key=f'microtube.{side}.{name}'; specs[key]=spec; defaults[key]=getattr(exchanger.bank,name)
        input_specs={k:v for k,v in INPUT_SPECS.items() if not k.startswith('air_')} if hasattr(exchanger.inputs,'external_stream') else INPUT_SPECS
        for name,spec in input_specs.items():
            key=f'thermal.{side}.{name}'; specs[key]=spec; defaults[key]=getattr(exchanger.inputs,name)
        if hasattr(exchanger.inputs,'external_stream'):
            for name,spec in STREAM_SPECS.items():
                key=f'external_stream.{side}.{name}'; specs[key]=spec
                defaults[key]=getattr(exchanger.inputs.external_stream,name)
    return specs,defaults


def build_machine(basis, configuration, physical, policies):
    p=physical
    try:
        large=p['volume.total_swept_m3']/(1+p['volume.swept_ratio']); small=p['volume.total_swept_m3']-large
        smin=small*p['volume.small_clearance_ratio']; lmin=large*p['volume.large_clearance_ratio']
        volumes=replace(configuration.machine_volumes,small_cylinder=CylinderVolumeLimits(smin,smin+small),large_cylinder=CylinderVolumeLimits(lmin,lmin+large))
        config=replace(configuration,machine_volumes=volumes,angular_speed=math.copysign(2*math.pi*p['operation.frequency_hz'],configuration.angular_speed),
            charge=replace(configuration.charge,total_mass=p['charge.total_mass_kg'])
                if policies['charge']=='explicit_inventory' else configuration.charge)
    except (ValueError,ArithmeticError) as error: raise PreflightRejection('invalid_parameterization',str(error)) from error
    exchangers=[]
    for side in ('heat_in','heat_out'):
        source=getattr(basis,side)
        if source is None: exchangers.append(None); continue
        bank={k:p[f'microtube.{side}.{k}'] for k in BANK_SPECS}
        inputs={k:p[f'thermal.{side}.{k}'] for k in INPUT_SPECS if f'thermal.{side}.{k}' in p}
        if hasattr(source.inputs,'external_stream'):
            inputs['external_stream']=replace(source.inputs.external_stream,
                **{k:p[f'external_stream.{side}.{k}'] for k in STREAM_SPECS})
        g=basis.data['geometry']; tag='hi' if side=='heat_in' else 'ho'; count='ni' if side=='heat_in' else 'no'
        cda=source.outlet_valve_cda_m2
        if policies['outlet_valve_cda']=='source_cda_times_count_ratio_v1':
            cda=g[f'seed_{tag}_cda']*bank['tube_count']/g[f'seed_{count}']
        try:
            exchanger=replace(source,bank=replace(source.bank,**bank),inputs=replace(source.inputs,**inputs),outlet_valve_cda_m2=cda)
            exchanger.build(); exchangers.append(exchanger)
        except (ValueError,ArithmeticError) as error: raise PreflightRejection('invalid_exchanger',str(error)) from error
    if basis.data['schema_version']==3 and all(x is not None and hasattr(x.inputs,'external_stream') for x in exchangers):
        inlet,outlet=(x.inputs.external_stream.inlet_temperature_k for x in exchangers)
        # The reservoir closures remain diagnostic references when walls supply heat.
        # Keep those references consistent with the actually declared boundaries.
        try:
            config=replace(config,cold_reservoir_temperature=outlet if config.motor_operation else inlet,
                hot_reservoir_temperature=inlet if config.motor_operation else outlet)
        except ValueError as error: raise PreflightRejection('invalid_exchanger',str(error)) from error
    return MachineDesign(config,heat_in=exchangers[0],heat_out=exchangers[1])


STREAM_SPECS={
    'inlet_temperature_k':ParameterSpec('K',positive=True),
    'mass_flow_kg_s':ParameterSpec('kg/s',positive=True),
    'cp_j_kg_k':ParameterSpec('J/(kg*K)',positive=True),
    'wall_conductance_w_k':ParameterSpec('W/K',positive=True)}
