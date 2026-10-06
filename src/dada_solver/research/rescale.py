"""Explicit capacity scaling of a stored candidate; no similarity is assumed."""
import copy
from dataclasses import asdict, replace
import hashlib
import json
import math
from pathlib import Path
import tempfile

from .report import inspect, select_records
from .snapshots import stored_study
from .schema import load_study, compile_study
from .machine_basis import machine_parameters
from .study_io import dumps


def extensive(name):
    return name in ('volume.total_swept_m3', 'charge.total_mass_kg') or (
        name.startswith('microtube.') and name.endswith('.additional_internal_volume_m3')) or (
        name.startswith('external_stream.') and name.rsplit('.', 1)[-1] in
        ('mass_flow_kg_s', 'wall_conductance_w_k')) or (
        name.startswith('thermal.') and name.endswith('.air_mass_flow_kg_s'))


def quantized_microtube_count(count, factor):
    """Nearest integer to N sqrt(s), ties to even; never clamp to one."""
    ideal = count * math.sqrt(factor)
    if not math.isfinite(ideal):
        raise ValueError('Microtube count scaling overflows.')
    result = round(ideal)
    if result < 1:
        raise ValueError('Scaled microtube count must be at least one.')
    return result


def scaled_microtube_dimensions(count, length, factor):
    """Preserve total physical microtube length after integer quantization."""
    if not math.isfinite(factor) or factor <= 0:
        raise ValueError('Capacity factor must be finite and positive.')
    new_count = quantized_microtube_count(count, factor)
    new_length = length if factor == 1 else length * (factor * count / new_count)
    if not math.isfinite(new_length) or new_length <= 0:
        raise ValueError('Scaled microtube length must be finite and positive.')
    return new_count, new_length


def rescale(source, candidate, factor, output, *, mode='capacity'):
    """Create a new portable schema-3 study centered on a selected candidate.

    Active bounds follow the capacity rule; intensive bounds and all constraints
    stay unchanged. Microtube count and length share total-length growth.
    Header packing is rebuilt by MicrotubeBank, not forced into exact similarity.
    """
    if mode != 'capacity': raise ValueError('Only capacity scaling is supported.')
    if isinstance(factor, bool) or not math.isfinite(factor) or factor <= 0:
        raise ValueError('Capacity factor must be finite and positive.')
    output = Path(output)
    if output.suffix != '.toml': raise ValueError('Study output must end in .toml.')
    basis_path = output.with_suffix('.basis.json')
    data = inspect(source)
    version = data['scientific']['schema_version']
    if type(version) is not int or version != 3:
        raise ValueError('Capacity rescale requires a schema-3 machine study.')
    record, = select_records(data, [candidate])
    changes = []
    def scale(value, name):
        new = value * factor
        if not math.isfinite(new): raise ValueError(f'Scaling overflows {name}.')
        changes.append(dict(quantity=name, before=value, after=new))
        return new
    with stored_study(data) as study:
        raw = copy.deepcopy(study.data)
        source_search_domain = raw['search']['domain']
        if source_search_domain == 'local_regions_v1':
            raw['search'] = dict(type='sobol', domain='fixed_global_bounds',
                seed=raw['search']['seed'], scramble=raw['search']['scramble'])
        raw['execution']['default_max_candidates'] = 512
        physical = dict(study.fixed_parameters, **record['physical'])
        design = compile_study(study).adapter.build(physical)
        b = copy.deepcopy(study.basis.data)
        if b['configuration'].get('hydraulic_flow_models') is not None:
            raise ValueError('Custom hydraulic flow models cannot be safely capacity-scaled.')
        # Freeze candidate inputs into the new basis as well as its declarations.
        config = asdict(design.configuration)
        config['gas'] = b['configuration']['gas']
        b['configuration'] = config
        for side in ('small_cylinder', 'large_cylinder'):
            for key in ('minimum', 'maximum'):
                config['machine_volumes'][side][key] = scale(config['machine_volumes'][side][key], f'basis.{side}.{key}')
        for key in ('cold_heat_exchanger', 'hot_heat_exchanger'):
            config['machine_volumes'][key] = scale(config['machine_volumes'][key], 'basis.'+key)
        config['charge']['total_mass'] = scale(config['charge']['total_mass'], 'basis.charge.total_mass')
        for key in ('cold_thermal_conductance', 'hot_thermal_conductance'):
            config[key] = scale(config[key], 'basis.'+key)
        for key, value in config['hydraulics'].items():
            if 'cda' in key: config['hydraulics'][key] = scale(value, 'basis.hydraulics.'+key)
            elif key != 'orifice_pressure_regularization':
                raise ValueError(f'Unknown hydraulic quantity cannot be safely scaled: {key}')
        microtube_scaling = {}
        for side in ('heat_in', 'heat_out'):
            exchanger = getattr(design, side)
            if exchanger is None: continue
            original = b[side]
            b[side] = asdict(exchanger)
            if 'family' in original: b[side]['family'] = original['family']
            row = b[side]
            old_count, old_length = exchanger.bank.tube_count, exchanger.bank.tube_length_m
            new_count, new_length = scaled_microtube_dimensions(old_count, old_length, factor)
            row['bank'].update(tube_count=new_count, tube_length_m=new_length)
            microtube_scaling[side] = dict(ideal_count=old_count*math.sqrt(factor),
                count=new_count, length_factor=1.0 if factor == 1 else factor*old_count/new_count)
            for key, before, after in (('tube_count', old_count, new_count),
                                       ('tube_length_m', old_length, new_length)):
                changes.append(dict(quantity=f'microtube.{side}.{key}', before=before, after=after))
            row['bank']['additional_internal_volume_m3'] = scale(row['bank']['additional_internal_volume_m3'], side+'.additional_internal_volume_m3')
            row['inputs']['extra_wall_capacity_j_k'] = scale(row['inputs']['extra_wall_capacity_j_k'], side+'.extra_wall_capacity_j_k')
            old_cda = exchanger.outlet_valve_cda_m2
            if exchanger.bank.circular_collectors:
                row['outlet_valve_cda_m2'] = replace(exchanger.bank,
                    tube_count=new_count, tube_length_m=new_length).dimensions()['conduit_area_m2']
            elif raw['policies']['outlet_valve_cda'] == 'source_cda_times_count_ratio_v1':
                row['outlet_valve_cda_m2'] = old_cda * new_count / old_count
            else:
                row['outlet_valve_cda_m2'] = old_cda * factor
            if not math.isfinite(row['outlet_valve_cda_m2']):
                raise ValueError(f'Scaling overflows {side}.outlet_valve_cda_m2.')
            changes.append(dict(quantity=side+'.outlet_valve_cda_m2',
                before=old_cda, after=row['outlet_valve_cda_m2']))
            if 'external_stream' in row['inputs']:
                for key in ('mass_flow_kg_s', 'wall_conductance_w_k'):
                    stream = row['inputs']['external_stream']
                    stream[key] = scale(stream[key], side+'.external_stream.'+key)
            else:
                row['inputs']['air_mass_flow_kg_s'] = scale(row['inputs']['air_mass_flow_kg_s'], side+'.air_mass_flow_kg_s')
            # Count-ratio policy already scales CdA through the candidate count.
            # Preserve its original reference slope, avoiding a second factor.
        b['geometry']['total_swept_m3'] = physical['volume.total_swept_m3'] * factor
        for side in ('small', 'large'):
            b['geometry'][side+'_clearance_ratio'] = physical['volume.'+side+'_clearance_ratio']
        # A source warm state is an initial guess only; scaling every conservative
        # component preserves its intensive state. Convergence is still required.
        if b['warm_start'] is not None:
            b['warm_start']['values'] = [v * factor for v in b['warm_start']['values']]
            b['warm_start']['wall_capacities_j_k'] = [v * factor for v in b['warm_start']['wall_capacities_j_k']]
            old_mass = study.basis.configuration.charge.total_mass
            if old_mass != design.configuration.charge.total_mass or raw['policies']['charge']!='explicit_inventory':
                b['warm_start'] = None
                raw['warm_start']['initial_source'] = 'uniform'
        specs, _ = machine_parameters(study.basis)
        declared = {r['name'] for r in raw['parameters']}
        for name, value in physical.items():
            if name not in declared and name in specs:
                raw['parameters'].append(dict(name=name, value=value, unit=specs[name].unit))
        for row in raw['parameters']:
            name = row['name']; initial_key = 'value' if 'value' in row else 'initial'
            row[initial_key] = physical[name]
            parts = name.split('.')
            coupled = (len(parts) == 3 and parts[0] == 'microtube'
                and parts[1] in microtube_scaling and parts[2] in ('tube_count', 'tube_length_m'))
            if coupled:
                side, coordinate = parts[1:]
                row[initial_key] = b[side]['bank'][coordinate]
                if initial_key == 'initial':
                    for key in ('lower', 'upper'):
                        row[key] = (quantized_microtube_count(row[key], factor)
                            if coordinate == 'tube_count' else
                            row[key] * microtube_scaling[side]['length_factor'])
            elif extensive(name):
                row[initial_key] = scale(row[initial_key], name)
                if initial_key == 'initial':
                    for key in ('lower', 'upper'):
                        row[key] *= factor
            if initial_key == 'initial' and row.get('kind') != 'choice':
                if (not all(math.isfinite(row[k]) for k in ('initial', 'lower', 'upper'))
                        or not row['lower'] < row['upper']
                        or not row['lower'] <= row['initial'] <= row['upper']):
                    raise ValueError(f'Scaled search bounds are invalid or collapse for {name}.')
        raw['study']['name'] += f' — capacity ×{factor:g}'
        raw['study']['parent_candidate_id'] = record['candidate_id']
        b['provenance']['capacity_scaling'] = dict(version=2, mode=mode, factor=factor,
            source_candidate_id=record['candidate_id'], source_study_id=data['study_id'],
            source_definition_id=data['definition_id'], source=str(source),
            source_basis_sha256=study.basis.sha256, changes=changes,
            ordinary_capacity='extensive quantities multiplied by factor',
            microtubes=dict(rule='N_ideal=N*sqrt(s); N_next=round_ties_to_even(N_ideal); L_next=s*N*L/N_next',
                sides=microtube_scaling, diameters='inner and outer unchanged'),
            search=dict(source_domain=source_search_domain, output_domain=raw['search']['domain']),
            active_bounds='ordinary extensive bounds multiplied by s; count endpoints quantized; length bounds multiplied by selected length factor; intensive bounds unchanged',
            constraints='unchanged; physical limits are not automatically relaxed',
            header_geometry='rebuilt by production bank geometry from quantized count; not forced capacity similarity',
            warm_start='original source guess scaled when inventory-compatible; not a cached solution')
        artifacts = {}
        for side, artifact in study.artifacts.items():
            path = output.with_name(output.stem+'.'+side+'.mechanism.json')
            artifacts[path] = json.dumps(artifact.data, indent=2)+'\n'
            raw['kinematics'][side]['artifact'] = path.name
        basis_text = json.dumps(b, indent=2, allow_nan=False)+'\n'
        raw['sources']['machine'] = dict(path=basis_path.name, sha256=hashlib.sha256(basis_text.encode()).hexdigest())
        contents = {output: dumps(raw), basis_path: basis_text, **artifacts}
        if any(p.exists() for p in contents): raise ValueError('Study or associated input already exists; choose a new output name.')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for path, text in contents.items(): (root/path.name).write_text(text)
            load_study(root/output.name)
        output.parent.mkdir(parents=True, exist_ok=True)
        for path, text in contents.items():
            with path.open('x') as stream: stream.write(text)
    return output
