"""Strict v1 study schema compiled to the existing campaign definition contract."""
from dataclasses import asdict, dataclass
import copy
import math
from pathlib import Path
import tomllib

from dada_solver.campaign.candidate import Candidate, canonical_json, content_hash
from dada_solver.campaign.definition import CampaignDefinition, runtime_identity
from dada_solver.campaign.history import atomic_json, atomic_text
from dada_solver.campaign.parameters import ParameterSpace, parameter_from_mapping
from dada_solver.campaign.runner import parse_budget
from dada_solver.exchangers.wall_cycle import WallCycleNumericalSettings
from dada_solver.sizing.configuration import _load_constraint, _load_objective
from dada_solver.wall_backend import WallBackendSettings, backend_identity
from .sixbar import PARAMETERS, SixBarStudyAdapter, load_basis

POLICIES = dict(volume_partition='fixed_total_swept_and_clearance_ratios_v1',
    charge='fixed_inventory_v1', outlet_valve_cda='source_cda_times_count_ratio_v1',
    valve_scaling_reference='reference.geometry', external_air='source_finite_mass_flow_and_film_v1',
    external_air_aerodynamic_losses='excluded_from_balance', fan_consumption='excluded_from_balance',
    mechanical_losses='unknown', useful_power='unavailable', isothermality='diagnostic_only', local_reflux='retain_signed_flows')
CONSTRAINTS = {'minimum_motor_power': ('required_power', 'W'), 'maximum_pressure': ('limit', 'Pa'),
               'maximum_temperature': ('limit', 'K'), 'maximum_absolute_mass_flow': ('limit', 'kg/s'),
               'valid_thermodynamic_model': (None, '1')}


def keys(data, required, label, optional=()):
    if not isinstance(data, dict): raise ValueError(f'{label} must be a table.')
    missing, extra = set(required)-set(data), set(data)-set(required)-set(optional)
    if missing or extra:
        raise ValueError(f'{label}: missing keys {sorted(missing)}; unknown keys {sorted(extra)}.')


def positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f'{name} must be finite and positive.')


@dataclass(frozen=True)
class StudyDefinition:
    path: Path
    source: str
    basis: object
    space: ParameterSpace
    study_id: str
    scientific: dict

    @property
    def data(self):
        return tomllib.loads(self.source)


def load_study(path, *, basis_path=None, artifact_directory=None):
    path = Path(path)
    source = path.read_text()
    raw = tomllib.loads(source)
    if type(raw.get('schema_version')) is int and raw['schema_version'] in (2,3):
        from .schema_v2 import load_study_v2
        return load_study_v2(path,basis_path=basis_path,artifact_directory=artifact_directory)
    canonical_json(raw)
    keys(raw, ('schema_version','study','sources','families','fixed','policies','objective',
               'constraints','parameters','search','numerical','warm_start','execution'), 'study file')
    if type(raw['schema_version']) is not int or raw['schema_version'] != 1:
        raise ValueError('Only study schema_version = 1 is supported.')
    keys(raw['study'], ('name','protocol','purpose','parent_candidate_id'), 'study')
    if raw['study']['protocol'] != 'fixed_pair_thermo5d':
        raise ValueError('This first validation release supports fixed_pair_thermo5d.')
    if any(not isinstance(v,str) or not v for v in raw['study'].values()):
        raise ValueError('Study labels and parent identity must be nonempty strings.')
    keys(raw['sources'], ('reference',), 'sources')
    reference = raw['sources']['reference']
    keys(reference, ('path','sha256','schema'), 'sources.reference')
    if reference['schema'] != 'fixed_sixbar_basis_v1': raise ValueError('Unsupported reference schema.')
    if raw['families'] != dict(kinematics='independent_six_bar', exchanger='microtube_air_wall'):
        raise ValueError('This protocol requires independent six-bar motion and microtube air/wall exchangers.')
    if raw['policies'] != POLICIES:
        raise ValueError('Unsupported or ambiguous derived study policies; see the generated template.')
    basis = load_basis(basis_path or path.parent/reference['path'], reference['sha256'])
    b, config = basis.data, basis.configuration
    fixed = raw['fixed']
    expected = dict(gas_inventory_kg=config.charge.total_mass,
        total_swept_volume_m3=b['geometry']['total_swept_m3'],
        small_clearance_ratio=b['geometry']['small_clearance_ratio'], large_clearance_ratio=b['geometry']['large_clearance_ratio'],
        frequency_hz=abs(config.angular_speed)/(2*math.pi), hot_air_inlet_K=basis.heat_in.inputs.air_inlet_temperature_k,
        cold_air_inlet_K=basis.heat_out.inputs.air_inlet_temperature_k,
        mechanisms='source_pair_exact_geometry_and_branches', mechanism_length_unit='crank_radius', mechanism_angle_unit='rad')
    keys(fixed, expected, 'fixed')
    for name, value in expected.items():
        if isinstance(value, str):
            if fixed[name] != value: raise ValueError(f'fixed.{name} must be {value}.')
        else:
            positive(fixed[name], f'fixed.{name}')
            if not math.isclose(fixed[name], value, rel_tol=1e-14, abs_tol=0):
                raise ValueError(f'fixed.{name} differs from the immutable source basis; create a new basis for changed physics.')
    if not 2 <= fixed['frequency_hz'] <= 10:
        raise ValueError('The demonstrator frequency must remain within 2–10 Hz.')
    objective = dict(type='maximize_thermal_efficiency', unit='1', heat_boundary='external_air_input', work_boundary='indicated_gas_work')
    if raw['objective'] != objective: raise ValueError('The objective must use indicated work / external-air heat input.')
    seen = set()
    for row in raw['constraints']:
        kind = row.get('type')
        if kind not in CONSTRAINTS or kind in seen: raise ValueError(f'Unsupported or duplicate constraint: {kind}')
        seen.add(kind)
        field, unit = CONSTRAINTS[kind]
        required = {'type','unit'} | ({field} if field else set()) | ({'boundary'} if kind=='minimum_motor_power' else set())
        keys(row, required, f'constraint {kind}')
        if row['unit'] != unit: raise ValueError(f'{kind} requires unit {unit}.')
        if field: positive(row[field], kind)
        if kind == 'minimum_motor_power' and row['boundary'] != 'indicated_gas_power':
            raise ValueError('Useful power is unavailable; the power floor is indicated gas power.')
    if seen != set(CONSTRAINTS): raise ValueError('All five declared physical constraints are required in this protocol.')
    parameters = []
    for row in raw['parameters']:
        name = row.get('name')
        if name not in PARAMETERS: raise ValueError(f'Unknown or derived parameter: {name}')
        kind, unit, _ = PARAMETERS[name]
        keys(row, ('name','kind','unit','lower','upper','initial','encoding' if kind=='integer' else 'transform'), name)
        if row['kind'] != kind or row['unit'] != unit: raise ValueError(f'{name} requires {kind} coordinates in {unit}.')
        for field in ('lower','upper','initial'): positive(row[field], f'{name}.{field}')
        parameters.append(parameter_from_mapping({k:v for k,v in row.items() if k != 'unit'}))
    space = ParameterSpace(tuple(parameters))
    if {p.name for p in parameters} != set(PARAMETERS): raise ValueError('Exactly five owned active parameters are required.')
    keys(raw['search'], ('type','seed','scramble','domain'), 'search')
    search = raw['search']
    if search['type'] != 'sobol' or search['domain'] != 'fixed_global_bounds' or type(search['scramble']) is not bool or type(search['seed']) is not int or search['seed'] < 0:
        raise ValueError('Search requires bounded Sobol, a nonnegative integer seed and boolean scramble.')
    num = raw['numerical']
    keys(num, ('wall_settings','maximum_cycles','backend','backend_fallback','exact_kinematics_cache','shared_replay','candidate_budget_seconds','domain_error_retry'), 'numerical')
    if num['wall_settings'] != 'source_exact' or num['backend_fallback'] != 'record_actual_backend' or num['exact_kinematics_cache'] is not True or num['shared_replay'] is not True:
        raise ValueError('Use exact source wall settings, exact caching, shared replay and recorded backend fallback.')
    if type(num['maximum_cycles']) is not int or num['maximum_cycles'] < 1: raise ValueError('maximum_cycles must be a positive integer.')
    positive(num['candidate_budget_seconds'], 'candidate_budget_seconds')
    if num['domain_error_retry'] not in ('none', 'once_safe_uniform_state_with_source_wall_temperatures'):
        raise ValueError('Unknown domain-error retry policy.')
    WallBackendSettings(num['backend'])
    if raw['warm_start'] != dict(initial_source='reference.initial_state', source_wall_capacities='reference.selected_candidate.hardware', subsequent_policy='campaign_compatible_nearest_state'):
        raise ValueError('Unsupported warm-start policy.')
    execution = raw['execution']
    keys(execution, ('default_budget','default_max_candidates','initial_evaluation_seconds','deadline_grace_seconds'), 'execution')
    positive(parse_budget(execution['default_budget']), 'default_budget')
    if type(execution['default_max_candidates']) is not int or execution['default_max_candidates'] < 1: raise ValueError('default_max_candidates must be a positive integer.')
    positive(execution['initial_evaluation_seconds'], 'initial_evaluation_seconds')
    if isinstance(execution['deadline_grace_seconds'], bool) or not isinstance(execution['deadline_grace_seconds'], (int,float)) or not math.isfinite(execution['deadline_grace_seconds']) or execution['deadline_grace_seconds'] < 0:
        raise ValueError('deadline_grace_seconds must be finite and nonnegative.')
    scientific = copy.deepcopy(raw)
    scientific.pop('execution')
    scientific['study'].pop('name')
    scientific['sources']['reference'].pop('path')
    scientific['basis'] = b
    study_id = content_hash(scientific)
    # Cheap construction only: never run integration during validation.
    SixBarStudyAdapter(basis).build({p.name:p.initial for p in parameters}).build()
    return StudyDefinition(path, source, basis, space, study_id, scientific)


class ResearchCampaignDefinition(CampaignDefinition):
    """Compile study inputs to the runner/evaluator contract without a second engine."""
    def __init__(self, study):
        from dataclasses import replace
        self.study = study
        raw = study.data
        self.path, self.source = study.path, study.source
        self.space = study.space
        self.families = dict(kinematics='independent_six_bar', exchanger='microtube')
        self.adapter = SixBarStudyAdapter(study.basis)
        self.configuration = study.basis.configuration
        self.fixed_parameters = {}
        self.objective = _load_objective({'type': raw['objective']['type']})
        self.constraints = tuple(_load_constraint({k:v for k,v in row.items() if k not in ('unit','boundary')}, None) for row in raw['constraints'])
        self.wall_numerical_settings = WallCycleNumericalSettings(**study.basis.data['wall_settings'])
        self.wall_backend = WallBackendSettings(raw['numerical']['backend'])
        self.adaptive_wall_acceleration = None
        self.configuration = replace(self.configuration, numerical=replace(self.configuration.numerical, maximum_cycles=raw['numerical']['maximum_cycles']))
        self.adapter = SixBarStudyAdapter(study.basis, self.configuration)
        self.initial_wall_state = study.basis.data['warm_start']
        self.safe_domain_retry = raw['numerical']['domain_error_retry'] != 'none'
        self.candidate_budget_seconds = raw['numerical']['candidate_budget_seconds']
        self.numerical_settings = dict(asdict(self.configuration.numerical), wall_cycle=asdict(self.wall_numerical_settings),
            wall_backend=backend_identity(self.wall_backend), research=raw['numerical'])
        self.seed, self.scramble = raw['search']['seed'], raw['search']['scramble']
        self.elite_size = 5
        self.initial_evaluation_seconds = raw['execution']['initial_evaluation_seconds']
        self.deadline_grace_seconds = raw['execution']['deadline_grace_seconds']
        self.maximum_candidates = raw['execution']['default_max_candidates']
        self.identity = dict(schema_version=1, definition_kind='research_v1', study_id=study.study_id,
            scientific=study.scientific, runtime=runtime_identity(), wall_backend=backend_identity(self.wall_backend))
        self.definition_id = content_hash(self.identity)

    def write_snapshots(self, directory):
        atomic_text(directory/'basis.json', self.study.basis.source)
        atomic_text(directory/'study.toml', self.study.source)
        atomic_json(directory/'study.json', dict(schema_version=1, study_id=self.study.study_id,
            definition_id=self.definition_id, name=self.study.data['study']['name'], scientific=self.study.scientific))


def compile_study(study):
    if study.data['schema_version'] in (2,3):
        from .schema_v2 import ResearchDefinitionV2
        return ResearchDefinitionV2(study)
    return ResearchCampaignDefinition(study)


def candidate_for_values(definition, values):
    """Keep requested physical values exact; encoding is for identity/distance only."""
    coordinates = definition.space.encode(values)
    payload = dict(schema_version=1, definition_id=definition.definition_id,
                   normalized=coordinates, physical=dict(values), families=definition.families,
                   numerical_settings=definition.numerical_settings)
    return Candidate(canonical_json(payload), content_hash(payload))
