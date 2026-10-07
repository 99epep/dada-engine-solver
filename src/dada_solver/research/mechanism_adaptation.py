"""Generate a portable paired thermodynamic study; never integrate a cycle."""
import copy
import hashlib
import json
import math
from pathlib import Path
import tempfile

import numpy as np

from .artifacts import MechanismArtifact, MechanismLibrary
from .families import parameter_specs
from .machine_basis import machine_parameters
from .motion_target import research_motion_source
from .schema import compile_study, load_study
from .study_io import dumps
from .synthesis_search import BOUNDS, BOUNDS_VERSION
from dada_solver.campaign.candidate import content_hash

BOUNDS_POLICY = 'centered_mechanical_local_regions_v1'
DEFAULT_RADIUS = .1


def _source_warm_start(source, identity, basis, design, raw):
    """Use the normal wall-state compatibility test, otherwise retain the source policy."""
    provenance = dict(policy=raw['warm_start']['initial_source'], selected_candidate_state=False)
    if (raw['warm_start']['initial_source']=='source_exact'
            and basis['configuration']['charge']['total_mass']!=design.configuration.charge.total_mass):
        raw['warm_start']['initial_source']='uniform'
        provenance=dict(policy='uniform',selected_candidate_state=False,
            reason='Inherited source state has a different inventory; use uniform unless a compatible current candidate state is available.')
    if Path(source).suffix == '.toml' or design.heat_in is None:
        return provenance
    from .report import inspect, select_records
    from dada_solver.campaign.evaluator import select_warm_start
    record, = select_records(inspect(source), [identity['candidate_id']])
    direction = 'motor' if design.configuration.motor_operation else 'receiver'
    selected, _ = select_warm_start([record], record['normalized'],
        'four_gas_volumes_m_U_plus_H_i_H_o_wall_energy', 'microtube_wall_10_state', direction)
    if selected is None:
        return dict(provenance, reason=provenance.get('reason', '') +
                    ' No compatible selected-candidate wall state; source policy retained.')
    state = selected.get('initial_guess_state') or selected['final_periodic_state']
    values = np.asarray(state['values'], float)
    caps = np.asarray(state['wall_capacities_j_k'], float)
    if (values.shape != (10,) or caps.shape != (2,) or not np.all(np.isfinite(values))
            or not np.all(values > 0) or not np.all(np.isfinite(caps)) or not np.all(caps > 0)
            or not np.isclose(values[:8:2].sum(), design.configuration.charge.total_mass, rtol=1e-10, atol=0)):
        return dict(provenance, reason=provenance.get('reason', '') +
                    ' Selected state fails fixed-inventory/state checks; source policy retained.')
    basis['warm_start'] = dict(values=values.tolist(), wall_capacities_j_k=caps.tolist(),
                              source_candidate_id=identity['candidate_id'])
    # load_machine_basis validates the stored inventory against this configuration.
    basis['configuration']['charge']['total_mass'] = design.configuration.charge.total_mass
    raw['warm_start']['initial_source'] = 'source_exact'
    return dict(policy='source_exact', selected_candidate_state=True,
                source_candidate_id=identity['candidate_id'], reusable_as_initial_guess_only=True,
                domain_error_retry=raw['numerical']['domain_error_retry'])


def paired_thermodynamic(source, library, family_id, output, *, candidate=None, radius=DEFAULT_RADIUS,
                         mechanical_constraints=()):
    """Freeze an exact Research source and release only a selected mechanical pair.

    Bounds are a numerical search policy, not a physical domain. Each coordinate
    gets a centered reference box with the width of the synthesis default bounds;
    positive lengths cap its half-width at half the current length. The existing
    refine Sobol policy explores +/- radius in normalized box coordinates. Angles
    remain unwrapped around the selected value, including across the cycle seam.
    """
    if isinstance(radius, bool) or not math.isfinite(radius) or not 0 < radius <= .5:
        raise ValueError('Paired thermodynamic radius must lie in (0, 0.5].')
    output = Path(output)
    if output.suffix != '.toml':
        raise ValueError('Paired thermodynamic output must end in .toml.')
    library = library if isinstance(library, MechanismLibrary) else MechanismLibrary.load(library)
    member = library.member(family_id)
    if set(member['mechanisms']) != {'small', 'large'}:
        raise ValueError('Paired thermodynamics requires a complete SMALL/LARGE library member.')
    artifacts = {side: MechanismArtifact.from_data(data) for side, data in member['mechanisms'].items()}
    for artifact in artifacts.values():
        if artifact.scientific['settings'].get('component') == 'primary':
            raise ValueError('Paired thermodynamics requires complete piston mechanisms, not primaries.')
    basis_path = output.with_suffix('.basis.json')
    artifact_paths = {side: output.with_name(output.stem+'.'+side+'.mechanism.json') for side in artifacts}
    paths = [output, basis_path, *artifact_paths.values()]
    if len({p.resolve() for p in paths}) != len(paths) or any(p.exists() for p in paths):
        raise ValueError('Study and associated input paths must be distinct and must not exist.')
    with research_motion_source(source, candidate=candidate) as (study, values, identity):
        raw = copy.deepcopy(study.data)
        domain = dict(source=identity, parameters=[copy.deepcopy(row) for row in raw['parameters']
                     if 'initial' in row and not row['name'].startswith('kinematics.')])
        domain['content_hash'] = content_hash(domain)
        physical = dict(study.fixed_parameters, **values)
        source_design = compile_study(study).adapter.build(physical)
        nonkinematic = {name: value for name, value in physical.items() if not name.startswith('kinematics.')}
        specs, _ = machine_parameters(study.basis, nonkinematic)
        charge_change = None
        if raw['policies']['charge'] != 'explicit_inventory':
            charge_change = dict(source_policy=raw['policies']['charge'], destination_policy='explicit_inventory',
                                 reason='Freeze selected candidate inventory during kinematics-only adaptation.')
            raw['policies']['charge'] = 'explicit_inventory'
            raw.pop('charge_reference', None)
            nonkinematic['charge.total_mass_kg'] = source_design.configuration.charge.total_mass
        rows = [dict(name=name, value=value, unit=specs[name].unit) for name, value in nonkinematic.items()]
        raw['kinematics'] = dict(coupling='independent')
        centers = {}
        reference_boxes = {}
        # Effective source constraints include those embedded in source artifacts.
        raw['mechanical_constraints'] = copy.deepcopy(list(study.mechanical_constraints))
        for side, artifact in artifacts.items():
            settings = artifact.scientific['settings']
            family = settings['family']
            from .margins import validate_mechanical_constraint
            for requirement in mechanical_constraints:
                validate_mechanical_constraint(requirement,family,scoped=False)
                scoped=dict(requirement,side=side)
                existing=next((r for r in raw['mechanical_constraints']
                    if (r['side'],r['metric'],r['relation'])==(side,scoped['metric'],scoped['relation'])),None)
                if existing is None:
                    raw['mechanical_constraints'].append(scoped)
                elif scoped['relation']=='equal' and scoped['limit']!=existing['limit']:
                    raise ValueError('Conflicting synthesis/source mechanical equality.')
                else:
                    combine=min if scoped['relation']=='maximum' else max
                    existing['limit']=combine(existing['limit'],scoped['limit'])
            raw['kinematics'][side] = dict(family=family, artifact=artifact_paths[side].name,
                                           sha256=artifact.content_hash)
            owned = parameter_specs(settings, side)
            reference_boxes[side] = {}
            for name, value in artifact.scientific['geometry'].items():
                spec = owned[name]
                if spec.kind != 'continuous':
                    continue  # Artifact supplies fixed branches and orientations verbatim.
                low, high = BOUNDS[family][name]
                half = (high-low)/2
                if spec.positive:
                    half = min(half, value/2)
                lower, upper = value-half, value+half
                key = f'kinematics.{side}.{name}'
                rows.append(dict(name=key, initial=value, lower=lower, upper=upper,
                                 kind='continuous', transform='linear', unit=spec.unit))
                centers[key] = value
                reference_boxes[side][name] = [lower, upper]
            topology = dict(side=side, metric='zero_crossing_count', relation='equal', limit=2, unit='1')
            matching = [r for r in raw['mechanical_constraints'] if r['side'] == side and r['metric'] == 'zero_crossing_count']
            if matching and any(r != topology for r in matching):
                raise ValueError(f'{side}: paired thermodynamics requires exactly two piston reversals.')
            if not matching:
                raw['mechanical_constraints'].append(topology)
        raw['parameters'] = rows
        raw['study']['name'] += ' — paired thermodynamic adaptation'
        raw['study']['parent_candidate_id'] = identity['candidate_id']
        # Source IDs trace the machine; the transformed center is identified by
        # its artifact hashes and exact mechanical coordinates in the basis.
        raw['search'] = dict(type='sobol', domain='local_regions_v1', seed=raw['search']['seed'],
            scramble=raw['search']['scramble'], radius_fraction=radius, allocation='round_robin',
            evaluate_centers=True, regions=[dict(id='selected_mechanical_pair',
                source_candidate_id=identity['candidate_id'], source_study_id=identity['study_id'], center=centers)])
        basis = study.basis.data
        warm = _source_warm_start(source, identity, basis, source_design, raw)
        provenance = dict(stage='paired_thermodynamic', source=identity, family_id=family_id,
            hardware_source_domain=domain,
            library_hash=library.to_data()['content_hash'],
            artifacts={s:a.content_hash for s,a in artifacts.items()},
            selected_member_metadata=member['metadata'], inventory_policy_change=charge_change,
            synthesis_mechanical_constraints=list(mechanical_constraints),
            warm_start=warm, bounds_policy=BOUNDS_POLICY, synthesis_bounds_policy=BOUNDS_VERSION,
            reference_boxes=reference_boxes, radius_fraction=radius,
            thermodynamic_objective=raw['objective'], initial_mechanical_coordinates=centers)
        if Path(source).suffix != '.toml':
            from .report import inspect, select_records
            record, = select_records(inspect(source), [identity['candidate_id']])
            provenance['source_assessment'] = {key:record.get(key) for key in ('objective','metrics','status')}
        basis['provenance']['paired_thermodynamic'] = provenance
        basis_text = json.dumps(basis, indent=2, allow_nan=False)+'\n'
        raw['sources']['machine'] = dict(path=basis_path.name, sha256=hashlib.sha256(basis_text.encode()).hexdigest())
        contents = {basis_path:basis_text, **{artifact_paths[s]:json.dumps(a.data, indent=2)+'\n' for s,a in artifacts.items()}}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for path, text in contents.items():
                (root/path.name).write_text(text)
            temporary_study = root/output.name
            temporary_study.write_text(dumps(raw))
            generated = load_study(temporary_study)
            compile_study(generated).adapter.build(dict(generated.fixed_parameters, **centers))
        contents[output] = dumps(raw)
        for path, text in contents.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open('x') as stream:
                stream.write(text)
    return output
