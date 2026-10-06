"""Research schema 2: independently owned fixed/active coordinates and motion families."""
from dataclasses import dataclass, asdict, replace
import copy
import json
import math
from pathlib import Path
import tomllib

from dada_solver.campaign.adapters import PreflightRejection
from dada_solver.campaign.candidate import canonical_json, content_hash
from dada_solver.campaign.definition import CampaignDefinition, runtime_identity
from dada_solver.campaign.history import atomic_json, atomic_text
from dada_solver.campaign.parameters import ParameterSpace, parameter_from_mapping
from dada_solver.campaign.runner import parse_budget
from dada_solver.composed_kinematics import ComposedKinematics
from dada_solver.exchangers.wall_cycle import WallCycleNumericalSettings
from dada_solver.sizing.configuration import _load_constraint, _load_objective
from dada_solver.sizing.constraints import MaximumMechanicalInputPower
from dada_solver.wall_backend import WallBackendSettings, backend_identity
from .artifacts import MechanismArtifact
from .families import validate_settings, build_side, side_metrics, ParameterSpec
from .machine_basis import load_machine_basis, machine_parameters, build_machine, validate_microtube_coordinate
from .margins import PHYSICAL_CONSTRAINTS, validate_mechanical_constraint, margin_record
from .schema import keys, positive

OBJECTIVE_UNITS = {
    'maximize_thermal_efficiency': '1',
    'maximize_cooling_cop': '1',
    'maximize_motor_power': 'W',
    'maximize_cooling_power': 'W',
    'maximize_cooling_power_per_total_microtube': 'W/microtube',
    'maximize_cooling_cop_times_power_per_total_microtube': 'W/microtube',
}
COOLING_OBJECTIVES = frozenset((
    'maximize_cooling_cop', 'maximize_cooling_power',
    'maximize_cooling_power_per_total_microtube',
    'maximize_cooling_cop_times_power_per_total_microtube',
))

POLICIES=dict(volume_partition='total_swept_and_clearance_ratios',charge='explicit_inventory',
    outlet_valve_cda='source_cda_times_count_ratio_v1',mechanical_losses='unknown',
    useful_power='unavailable',external_air_aerodynamic_losses='excluded_from_balance',
    fan_consumption='excluded_from_balance',local_reflux='retain_signed_flows')


POLICIES_V3={k:v for k,v in POLICIES.items() if k not in ('external_air_aerodynamic_losses','fan_consumption')}
POLICIES_V3.update(external_loop_hydraulics='unmodelled',external_pump_fan_consumption='excluded_from_balance')

@dataclass(frozen=True)
class StudyV2:
    path: Path
    source: str
    basis: object
    space: ParameterSpace
    study_id: str
    scientific: dict
    fixed_parameters: dict
    artifacts: dict
    settings: dict
    mechanical_constraints: tuple

    @property
    def data(self): return tomllib.loads(self.source)


def load_study_v2(path, *, basis_path=None, artifact_directory=None):
    path=Path(path); source=path.read_text(); raw=tomllib.loads(source); canonical_json(raw)
    keys(raw,('schema_version','study','sources','kinematics','parameters','policies','objective','constraints',
              'mechanical_constraints','screening','search','numerical','warm_start','execution'),'study schema_version = 2',('charge_reference',))
    if type(raw['schema_version']) is not int or raw['schema_version'] not in (2,3): raise ValueError('Expected schema_version = 2.')
    keys(raw['study'],('name','protocol','purpose'),'study',('parent_candidate_id',))
    if raw['study']['protocol']!='machine_design': raise ValueError('V2 protocol must be machine_design.')
    if any(not isinstance(x,str) or not x for x in raw['study'].values()): raise ValueError('Study labels must be nonempty strings.')
    keys(raw['sources'],('machine',),'sources'); ref=raw['sources']['machine']
    keys(ref,('path','sha256'),'sources.machine')
    basis=load_machine_basis(basis_path or path.parent/ref['path'],ref['sha256'])
    if basis.data['schema_version']!=raw['schema_version']: raise ValueError('Study and machine schema versions must match.')
    specs,defaults=machine_parameters(basis, [p.get('name') for p in raw['parameters']])
    expected_policy=POLICIES_V3 if raw['schema_version']==3 else POLICIES
    policy=raw['policies']
    from .charge import POLICY as REFERENCE_CHARGE
    if dict(policy,outlet_valve_cda=expected_policy['outlet_valve_cda'],charge='explicit_inventory')!=expected_policy or policy.get('charge') not in ('explicit_inventory',REFERENCE_CHARGE) or policy.get('outlet_valve_cda') not in ('source_cda_times_count_ratio_v1','fixed_source_cda','geometry_conduit_area_v1'):
        raise ValueError('Unsupported V2 physical policy.')
    if policy['outlet_valve_cda']=='geometry_conduit_area_v1' and any(f'microtube.{side}.collector_half_angle_deg' not in specs for side in ('heat_in','heat_out')):
        raise ValueError('geometry_conduit_area_v1 requires circular collectors on both exchangers.')
    if policy['charge']==REFERENCE_CHARGE:
        from dada_solver.fluids import CaloricallyPerfectGas
        if type(basis.configuration.gas) is not CaloricallyPerfectGas:
            raise ValueError('Reference-pressure inventory currently requires CaloricallyPerfectGas; no real-gas ideal filling approximation is implicit.')
        if 'charge_reference' not in raw: raise ValueError('Derived charge requires charge_reference.')
        reference=raw['charge_reference']
        keys(reference,('pressure_pa','temperature_k','volume_state'),'charge_reference')
        for key in ('pressure_pa','temperature_k'): positive(reference[key],'charge_reference.'+key)
        if reference['volume_state']!='maximum_total_gas_volume':
            raise ValueError('charge_reference.volume_state must be maximum_total_gas_volume.')
        if any(row.get('name')=='charge.total_mass_kg' for row in raw['parameters']):
            raise ValueError('charge.total_mass_kg is derived and must not be declared under the reference-pressure policy.')
        if raw['warm_start'].get('initial_source')=='source_exact':
            raise ValueError('source_exact is incompatible with geometry-derived reference-pressure charge; use uniform.')
        specs.pop('charge.total_mass_kg'); defaults.pop('charge.total_mass_kg')
    elif 'charge_reference' in raw:
        raise ValueError('charge_reference requires the reference-pressure charge policy.')
    kin=raw['kinematics']; keys(kin,('coupling','small','large'),'kinematics')
    if kin['coupling'] not in ('independent','shared_crank'): raise ValueError('Choose independent or shared_crank coupling.')
    settings={}; artifacts={}; mechanical=list(raw['mechanical_constraints'])
    for side in ('small','large'):
        cfg=copy.deepcopy(kin[side])
        if ('artifact' in cfg)!=('sha256' in cfg): raise ValueError('An artifact requires its content hash.')
        if 'artifact' in cfg:
            artifact=MechanismArtifact.load(Path(artifact_directory)/(side+'.mechanism.json') if artifact_directory else path.parent/cfg['artifact'],cfg['sha256'])
            if artifact.scientific['settings']['family']!=cfg['family']: raise ValueError('Mechanism artifact has a different family.')
            explicit={k:v for k,v in cfg.items() if k not in ('artifact','sha256')}
            for k,v in explicit.items():
                if artifact.scientific['settings'].get(k)!=v:
                    raise ValueError('Artifact settings cannot be relabelled; create a new artifact for changed scale/output.')
            cfg=dict(artifact.scientific['settings'])
            defaults.update({f'kinematics.{side}.{k}':v for k,v in artifact.scientific['geometry'].items()})
            mechanical.extend(dict(row,side=side) for row in artifact.scientific['constraints'])
            artifacts[side]=artifact
        settings[side]=cfg
        specs.update({f'kinematics.{side}.{k}':v for k,v in validate_settings(cfg,side).items()})
    shared=kin['coupling']=='shared_crank'
    if shared:
        if any(s['family']!='four_bar' for s in settings.values()): raise ValueError('Shared-crank coupling currently requires two four-bar assemblies in a common frame.')
        for key,unit,kind in [('phase_rad','rad','continuous'),('crank_direction','1','branch')]:
            values=[defaults.pop(f'kinematics.{side}.{key}',None) for side in ('small','large')]
            for side in ('small','large'): specs.pop(f'kinematics.{side}.{key}')
            if all(v is not None for v in values):
                if values[0]!=values[1]: raise ValueError('Shared-crank artifacts must already use one crank phase/direction and common frame.')
                defaults[f'kinematics.shared.{key}']=values[0]
            specs[f'kinematics.shared.{key}']=ParameterSpec(unit,kind)
        if settings['small'].get('crank_radius_m')!=settings['large'].get('crank_radius_m'):
            raise ValueError('A shared crank has one physical radius.')
    fixed=dict(defaults); active=[]; seen=set()
    for row in raw['parameters']:
        name=row.get('name')
        if name not in specs or name in seen: raise ValueError(f'Unknown, duplicate or wrong-family parameter: {name}')
        seen.add(name); spec=specs[name]
        if 'value' in row:
            keys(row,('name','value','unit'),name)
            spec.validate(row['value']); fixed[name]=row['value']
            validate_microtube_coordinate(name,row['value'])
        elif spec.kind=='choice':
            keys(row,('name','initial','choices','unit','kind'),name)
            if row['kind']!='choice': raise ValueError(f'Incorrect parameter kind for {name}.')
            parameter=parameter_from_mapping({k:v for k,v in row.items() if k!='unit'})
            for value in parameter.choices: spec.validate(value)
            active.append(parameter); fixed.pop(name,None)
        else:
            if spec.kind not in ('continuous','integer'): raise ValueError(f'{name} must remain fixed; branches/categories are not continuous search coordinates.')
            mode='encoding' if spec.kind=='integer' else 'transform'
            keys(row,('name','initial','lower','upper','unit','kind',mode),name)
            if row['kind']!=spec.kind: raise ValueError(f'Incorrect parameter kind for {name}.')
            for k in ('initial','lower','upper'):
                spec.validate(row[k]); validate_microtube_coordinate(name,row[k])
            active.append(parameter_from_mapping({k:v for k,v in row.items() if k!='unit'})); fixed.pop(name,None)
        if row['unit']!=spec.unit: raise ValueError(f'{name} requires unit {spec.unit}.')
    if set(fixed)|{p.name for p in active}!=set(specs):
        raise ValueError(f'Missing family coordinates: {sorted(set(specs)-set(fixed)-{p.name for p in active})}')
    for name,value in fixed.items():
        specs[name].validate(value); validate_microtube_coordinate(name,value)
    space=ParameterSpace(tuple(active),allow_empty=True)
    constraints=raw['constraints']; names=set()
    for row in constraints:
        kind=row.get('type')
        if kind not in PHYSICAL_CONSTRAINTS or kind in names: raise ValueError(f'Unknown/duplicate physical constraint {kind}.')
        names.add(kind); field,_,unit=PHYSICAL_CONSTRAINTS[kind]
        keys(row,{'type','unit'}|({field} if field else set()),kind)
        if row['unit']!=unit: raise ValueError(f'{kind} requires unit {unit}.')
        if field: positive(row[field],kind)
    if 'valid_thermodynamic_model' not in names: raise ValueError('An explicit model-validity constraint is required.')
    objective=raw['objective']; keys(objective,('type','unit'),'objective')
    if objective['type'] not in OBJECTIVE_UNITS or objective['unit'] != OBJECTIVE_UNITS[objective['type']]:
        raise ValueError('Unsupported objective or unit.')
    if basis.configuration.motor_operation != (objective['type'] not in COOLING_OBJECTIVES):
        raise ValueError('The objective must match the basis operating direction.')
    signatures=set()
    for row in mechanical:
        if row.get('side') not in settings: raise ValueError('Mechanical constraint needs small or large side.')
        validate_mechanical_constraint(row,settings[row['side']]['family'])
        signature=(row['side'],row['metric'],row['relation'])
        if signature in signatures: raise ValueError('Duplicate mechanical constraint, including artifact constraints.')
        signatures.add(signature)
    screen=raw['screening']; keys(screen,('samples',),'screening',('maximum_large_enclosed_volume_m3',))
    if type(screen['samples']) is not int or screen['samples']<360: raise ValueError('Mechanical screen requires at least 360 samples per cycle.')
    if 'maximum_large_enclosed_volume_m3' in screen:
        positive(screen['maximum_large_enclosed_volume_m3'],'maximum_large_enclosed_volume_m3')
    from .local_search import validate_search
    validate_search(raw['search'],space)
    num=raw['numerical']; keys(num,('maximum_cycles','backend','candidate_budget_seconds','domain_error_retry'),'numerical')
    if type(num['maximum_cycles']) is not int or num['maximum_cycles']<1: raise ValueError('maximum_cycles must be a positive integer.')
    WallBackendSettings(num['backend']); positive(num['candidate_budget_seconds'],'candidate_budget_seconds')
    if num['domain_error_retry'] not in ('none','once_safe_uniform_state_with_source_wall_temperatures'): raise ValueError('Unknown safe-retry policy.')
    keys(raw['warm_start'],('initial_source',),'warm_start')
    if raw['warm_start']['initial_source'] not in ('uniform','source_exact'): raise ValueError('Choose uniform or source_exact initial state.')
    if raw['warm_start']['initial_source']=='source_exact' and (basis.data['warm_start'] is None or 'charge.total_mass_kg' not in fixed or fixed['charge.total_mass_kg']!=basis.configuration.charge.total_mass):
        raise ValueError('Exact source warm state requires its fixed original inventory.')
    execution=raw['execution']; keys(execution,('default_budget','default_max_candidates','initial_evaluation_seconds','deadline_grace_seconds'),'execution')
    positive(parse_budget(execution['default_budget']),'default_budget'); positive(execution['initial_evaluation_seconds'],'initial_evaluation_seconds')
    if type(execution['default_max_candidates']) is not int or execution['default_max_candidates']<1: raise ValueError('default_max_candidates must be a positive integer.')
    if isinstance(execution['deadline_grace_seconds'],bool) or not isinstance(execution['deadline_grace_seconds'],(float,int)) or not math.isfinite(execution['deadline_grace_seconds']) or execution['deadline_grace_seconds']<0: raise ValueError('Invalid deadline grace.')
    scientific=copy.deepcopy(raw); scientific.pop('execution'); scientific['study'].pop('name'); scientific['sources']['machine'].pop('path')
    scientific['basis']=basis.data; scientific['fixed_parameters']=fixed
    scientific['mechanical_constraints']=mechanical
    for side,artifact in artifacts.items():
        scientific['kinematics'][side].pop('artifact')
        scientific['kinematics'][side]['mechanism']=artifact.scientific
    study=StudyV2(path,source,basis,space,content_hash(scientific),scientific,fixed,artifacts,settings,tuple(mechanical))
    # Configuration validation constructs and screens once, never integrates.
    V2Adapter(study).build(dict(fixed,**{p.name:p.initial for p in active})).build()
    return study


class V2Adapter:
    def __init__(self,study,configuration=None):
        self.study=study; self.configuration=configuration or study.basis.configuration

    def build(self,physical):
        study=self.study; raw=study.data
        expected=set(study.fixed_parameters)|{p.name for p in study.space.parameters}
        if set(physical)!=expected: raise PreflightRejection('invalid_parameterization','Candidate parameters do not match the declared ownership.')
        design=build_machine(study.basis,self.configuration,physical,raw['policies'])
        rows=[]; laws=[]; measured=[]
        ceiling=raw['screening'].get('maximum_large_enclosed_volume_m3')
        if ceiling is not None:
            rows.append(margin_record('maximum_large_enclosed_volume',design.configuration.machine_volumes.large_cylinder.maximum,ceiling,'maximum','m^3',method='exact volume limits'))
        for side in ('small','large'):
            prefix=f'kinematics.{side}.'; p={k[len(prefix):]:v for k,v in physical.items() if k.startswith(prefix)}
            if raw['kinematics']['coupling']=='shared_crank':
                p.update({key:physical[f'kinematics.shared.{key}'] for key in ('phase_rad','crank_direction')})
            try:
                law,geometry=build_side(study.settings[side],p,side,getattr(design.configuration.machine_volumes,side+'_cylinder'))
                metrics=side_metrics(study.settings[side],law,geometry,raw['screening']['samples'])
                measured.append(dict(side=side,family=study.settings[side]['family'],values=metrics,
                    samples=raw['screening']['samples'],method='production geometry with sampled motion screens; not a continuous proof'))
                laws.append(law)
                for declaration in study.mechanical_constraints:
                    if declaration['side']!=side: continue
                    metric=declaration['metric']
                    method='analytic spline extrema' if study.settings[side]['family']=='free_spline' and metric.startswith('maximum_absolute_') else f'sampled full cycle ({raw["screening"]["samples"]} angles); not a continuous guarantee'
                    if study.settings[side]['family']=='slider_crank' and metric in ('closure_margin','minimum_rod_axis_cosine','stroke_over_crank'):
                        method='analytic slider-crank full-cycle bound'
                    rows.append(margin_record(f'{side}.{metric}.{declaration["relation"]}',metrics.get(metric),declaration['limit'],declaration['relation'],declaration['unit'],method=method))
            except (ValueError,ArithmeticError) as error:
                rows.append(margin_record(f'{side}.geometry',False,True,'equal','1',method='production closure / motion domain'))
                raise PreflightRejection('invalid_kinematics',str(error),rows) from error
        if any(not r['satisfied'] for r in rows):
            raise PreflightRejection('invalid_kinematics','Mechanical or volume constraint violated before integration.',rows)
        design=replace(design,kinematics=ComposedKinematics(*laws,tuple(rows),tuple(measured)))
        from .charge import POLICY, apply_reference_charge
        if raw['policies']['charge']==POLICY:
            try: design=apply_reference_charge(design,raw['charge_reference'])
            except (ValueError,ArithmeticError) as error:
                raise PreflightRejection('invalid_parameterization',str(error)) from error
        return design


class ResearchDefinitionV2(CampaignDefinition):
    def __init__(self,study):
        self.study=study; raw=study.data
        self.path,self.source=study.path,study.source; self.space=study.space
        self.fixed_parameters=study.fixed_parameters
        self.configuration=replace(study.basis.configuration,numerical=replace(study.basis.configuration.numerical,maximum_cycles=raw['numerical']['maximum_cycles']))
        self.adapter=V2Adapter(study,self.configuration)
        self.families={side:study.settings[side]['family'] for side in ('small','large')}
        self.families.update(kinematics='composed',exchanger='microtube' if study.basis.heat_in is not None else 'reservoir')
        from dada_solver.campaign.objectives import (
            MaximizeCoolingPowerPerTotalMicrotube,
            MaximizeCoolingCopTimesPowerPerTotalMicrotube,
        )
        microtube_objectives = {
            'maximize_cooling_power_per_total_microtube': MaximizeCoolingPowerPerTotalMicrotube(),
            'maximize_cooling_cop_times_power_per_total_microtube': MaximizeCoolingCopTimesPowerPerTotalMicrotube(),
        }
        self.objective=microtube_objectives.get(raw['objective']['type'])
        if self.objective is None:
            self.objective=_load_objective({'type':raw['objective']['type']})
        self.constraints=tuple(MaximumMechanicalInputPower(row['limit']) if row['type']=='maximum_mechanical_input_power' else _load_constraint({k:v for k,v in row.items() if k!='unit'},None) for row in raw['constraints'])
        self.wall_numerical_settings=WallCycleNumericalSettings(**study.basis.data['wall_settings']) if study.basis.heat_in is not None else None
        self.wall_backend=WallBackendSettings(raw['numerical']['backend']); self.adaptive_wall_acceleration=None
        self.initial_wall_state=study.basis.data['warm_start'] if raw['warm_start']['initial_source']=='source_exact' else None
        self.safe_domain_retry=raw['numerical']['domain_error_retry']!='none'
        self.candidate_budget_seconds=raw['numerical']['candidate_budget_seconds']
        self.numerical_settings=dict(asdict(self.configuration.numerical),wall_cycle=asdict(self.wall_numerical_settings) if self.wall_numerical_settings else None,
            wall_backend=backend_identity(self.wall_backend),research=raw['numerical'],screening=raw['screening'])
        self.search_settings=raw['search']
        self.seed,self.scramble=raw['search']['seed'],raw['search']['scramble']; self.elite_size=5
        for name in ('initial_evaluation_seconds','deadline_grace_seconds'): setattr(self,name,raw['execution'][name])
        self.maximum_candidates=raw['execution']['default_max_candidates']
        self.identity=dict(schema_version=raw['schema_version'],definition_kind='research_v'+str(raw['schema_version']),study_id=study.study_id,scientific=study.scientific,runtime=runtime_identity(),wall_backend=backend_identity(self.wall_backend))
        if raw['schema_version']==3:
            self.families['exchanger']='external_stream_wall' if study.basis.heat_in is not None and hasattr(study.basis.heat_in.inputs,'external_stream') else self.families['exchanger']
            fluid=self.configuration.gas
            if hasattr(fluid,'identity'):
                self.identity['fluid']=fluid.identity
                self.numerical_settings['fluid']=fluid.identity
                self.identity['compiled_kernel']='tabulated_rho_u_wall_v1'
            else: self.identity['compiled_kernel']='calorically_perfect_wall_v1'
        self.definition_id=content_hash(self.identity)

    def write_snapshots(self,directory):
        atomic_text(directory/'basis.json',self.study.basis.source)
        atomic_text(directory/'study.toml',self.study.source)
        atomic_json(directory/'study.json',dict(schema_version=self.study.data['schema_version'],study_id=self.study.study_id,definition_id=self.definition_id,name=self.study.data['study']['name'],scientific=self.study.scientific))
        for side,artifact in self.study.artifacts.items(): atomic_json(directory/(side+'.mechanism.json'),artifact.data)
