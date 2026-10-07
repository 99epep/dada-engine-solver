"""Portable motion-only Research study from an exact selected machine candidate."""
import copy
import hashlib
import json
import math
from pathlib import Path
import tempfile

import numpy as np
from dada_solver.free_kinematics import FreeMotionDefinition
from .motion_target import research_motion_source, target_from_study, MotionTarget
from .motion_refit import MotionRefitRequest, COUNT, STEP
from .schema import compile_study, load_study
from .machine_basis import machine_parameters
from .study_io import dumps


def shape_bounds(controls, angular_radius=.15):
    """Chart bounding box of a spherical cap in canonical control space.

    The cap radius is a configurable search-region policy, not a physical domain.
    Its stereographic image is a ball. This returns the ball's coordinate box;
    joint box corners can exceed the cap radius, which is reported explicitly.
    """
    if isinstance(angular_radius,bool) or not math.isfinite(angular_radius) or not 0 < angular_radius < math.pi/2:
        raise ValueError('Canonical shape search radius must lie in (0, pi/2) radians.')
    definition = FreeMotionDefinition(tuple(controls),1.,2.)
    z = np.asarray(definition.to_shape_coordinates())
    square = float(z@z)
    sphere = np.r_[2*z,square-1]/(1+square)
    denominator = math.cos(angular_radius)-sphere[-1]
    if denominator <= 128*np.finfo(float).eps:
        raise ValueError('Shape search cap touches the excluded chart pole; choose a smaller radius.')
    center = sphere[:-1]/denominator
    radius = math.sin(angular_radius)/denominator
    lower, upper = center-radius, center+radius
    if not (np.all(np.isfinite(lower)) and np.all(np.isfinite(upper)) and np.all(lower < z) and np.all(z < upper)):
        raise ValueError('Shape search box is numerically unresolved.')
    corner_norm = float(np.linalg.norm(np.maximum(abs(lower),abs(upper))))
    return lower, upper, dict(method='stereographic image of canonical spherical cap, coordinate envelope',
        canonical_angular_radius_rad=angular_radius,chart_ball_center=center.tolist(),chart_ball_radius=radius,
        box_is_cap=False,minimum_box_angle_from_excluded_pole_rad=2*math.atan2(1.,corner_norm))


def refit_study(source, output, *, candidate=None, report=None, policy=None,
                shape_radius=.15, phase_radius_fraction=.5):
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
    if isinstance(shape_radius,bool) or not math.isfinite(shape_radius) or not 0 < shape_radius < math.pi/2:
        raise ValueError('Canonical shape search radius must lie in (0, pi/2) radians.')
    if isinstance(phase_radius_fraction,bool) or not math.isfinite(phase_radius_fraction) or not 0 < phase_radius_fraction <= 1:
        raise ValueError('Phase search radius fraction must lie in (0, 1].')
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
        result = MotionRefitRequest(target,**({} if policy is None else {'policy':policy})).execute()
        search_regions = {}
        for side in ('small','large'):
            fitted = result.data['scientific']['sides'][side]
            low,high,metadata = shape_bounds(fitted['controls'],shape_radius)
            search_regions[side] = dict(metadata,phase_radius_rad=STEP*phase_radius_fraction)
            for i,value in enumerate(fitted['shape_coordinates']):
                rows.append(dict(name=f'kinematics.{side}.shape_{i}',kind='continuous',unit='1',transform='linear',
                                 initial=value,lower=float(low[i]),upper=float(high[i])))
            phase = fitted['phase_rad']
            rows.append(dict(name=f'kinematics.{side}.phase_rad',kind='continuous',unit='rad',transform='linear',
                             initial=phase,lower=phase-STEP*phase_radius_fraction,upper=phase+STEP*phase_radius_fraction))
            raw['kinematics'][side] = dict(family='free_spline',representation='shape_coordinates',count=COUNT)
        raw['parameters'] = rows
        raw['study']['name'] += ' — motion refit'
        raw['study']['parent_candidate_id'] = identity['candidate_id']
        raw['search'] = dict(type='sobol',domain='fixed_global_bounds',seed=raw['search']['seed'],
                             scramble=raw['search']['scramble'],evaluate_initial=True)
        basis = study.basis.data
        basis['provenance']['motion_refit'] = dict(source=identity,target_hash=target.content_hash,
            result_hash=result.content_hash,search_regions=search_regions,inventory_policy_change=charge_change)
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
