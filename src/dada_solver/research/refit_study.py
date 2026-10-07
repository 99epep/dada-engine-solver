"""Portable motion-only Research study from an exact selected machine candidate."""
import copy
import hashlib
import json
import math
from pathlib import Path
import tempfile

from .motion_target import research_motion_source, target_from_study, MotionTarget
from .motion_refit import MotionRefitRequest, NAMES
from .schema import compile_study, load_study
from .machine_basis import machine_parameters
from .study_io import dumps


DEFAULT_RADIUS = .1


def local_bounds(name, value, radius):
    """Local numeric search policy; open structural domains are never crossed."""
    if name == 'small_max_deg': low,high=value-180.,value+180.
    elif name.endswith('duration_deg'): low,high=.1,359.9
    elif name.endswith('curvature'): low,high=value/2,value*1.5
    elif name.endswith('width_rel'): low,high=.001,1.999
    else: low,high=.0001,.9999
    half=radius*(high-low)
    return max(low,value-half),min(high,value+half)


def refit_study(source, output, *, candidate=None, report=None, policy=None,
                radius=DEFAULT_RADIUS):
    """Fit geometry, freeze the exact source machine, validate and write a study.

    A reference-pressure source inventory is frozen to its computed candidate
    mass: a motion-only search must not refill the machine when extrema change.
    No periodic integration or thermodynamic optimization is performed.
    """
    output = Path(output)
    if output.suffix != '.toml': raise ValueError('Refit study output must end in .toml.')
    basis_path = output.with_suffix('.basis.json')
    paths = [output,basis_path]+([] if report is None else [Path(report)])
    if len({p.resolve() for p in paths}) != len(paths) or any(p.exists() for p in paths):
        raise ValueError('Refit study/report paths must be distinct and must not exist.')
    if isinstance(radius,bool) or not math.isfinite(radius) or not 0 < radius <= .5:
        raise ValueError('Structured study search radius must lie in (0, 0.5].')
    with research_motion_source(source,candidate=candidate) as (study,values,identity):
        if study.data['kinematics']['coupling'] != 'independent':
            raise ValueError('Motion refit study requires independent source cylinder laws.')
        target = target_from_study(study,values,samples=1441)
        raw_target = target.scientific
        raw_target['source'].update(identity)
        target = MotionTarget.create(raw_target['angles_rad'],raw_target['sides'],source=raw_target['source'],provenance=target.data['provenance'])
        raw = copy.deepcopy(study.data)
        physical = dict(study.fixed_parameters,**values)
        design = compile_study(study).adapter.build(physical)
        nonkinematic = {name:value for name,value in physical.items() if not name.startswith('kinematics.')}
        specs, _ = machine_parameters(study.basis,nonkinematic)
        charge_change = None
        if raw['policies']['charge'] != 'explicit_inventory':
            charge_change = dict(source_policy=raw['policies']['charge'],destination_policy='explicit_inventory',
                                 reason='freeze selected candidate inventory during motion-only optimization')
            raw['policies']['charge'] = 'explicit_inventory'
            raw.pop('charge_reference',None)
            nonkinematic['charge.total_mass_kg'] = design.configuration.charge.total_mass
        rows = [dict(name=name,value=value,unit=specs[name].unit) for name,value in nonkinematic.items()]
        # Constraints stay intact; unsupported family-specific constraints are
        # rejected by the normal study loader, never dropped or reinterpreted.
        seeds = {name:physical[key] for name in NAMES for side in ('small','large')
                 if (key := f'kinematics.{side}.{name}') in physical}
        result = MotionRefitRequest(target,initial_parameters=seeds,**({} if policy is None else {'policy':policy})).execute()
        search_regions = {}
        fitted=result.data['scientific']['parameters']
        for side in ('small','large'):
            raw['kinematics'][side] = dict(family='structured_c2_15p')
        for name,value in fitted.items():
            side='small' if name.startswith('small_') else 'large'
            low,high=local_bounds(name,value,radius)
            key=f'kinematics.{side}.{name}'
            rows.append(dict(name=key,kind='continuous',unit='deg' if name.endswith('_deg') else '1',
                             transform='linear',initial=value,lower=low,upper=high))
            search_regions[key]=dict(lower=low,upper=high,initial=value)
        raw['parameters'] = rows
        raw['study']['name'] += ' — motion refit'
        raw['study']['parent_candidate_id'] = identity['candidate_id']
        raw['search'] = dict(type='sobol',domain='fixed_global_bounds',seed=raw['search']['seed'],
                             scramble=raw['search']['scramble'],evaluate_initial=True)
        basis = study.basis.data
        basis['provenance']['motion_refit'] = dict(source=identity,target_hash=target.content_hash,
            result_hash=result.content_hash,search_regions=search_regions,radius_fraction=radius,
            search_policy='structured_local_coordinate_box_v1',inventory_policy_change=charge_change)
        basis_text = json.dumps(basis,indent=2,allow_nan=False)+'\n'
        raw['sources']['machine'] = dict(path=basis_path.name,sha256=hashlib.sha256(basis_text.encode()).hexdigest())
        study_text = dumps(raw)
        # Validate the complete portable inputs before creating any user file.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/basis_path.name).write_text(basis_text)
            (root/output.name).write_text(study_text)
            new_study = load_study(root/output.name)
            definition = compile_study(new_study)
            definition.adapter.build(dict(new_study.fixed_parameters,**{p.name:p.initial for p in new_study.space.parameters}))
        contents = {basis_path:basis_text,output:study_text}
        if report is not None: contents[Path(report)] = json.dumps(result.data,indent=2,allow_nan=False)+'\n'
        for path,text in contents.items():
            path.parent.mkdir(parents=True,exist_ok=True)
            with path.open('x') as stream: stream.write(text)
    return output, result
