"""Small portable starting studies; these are references, not family optima."""
import copy
import hashlib
import json
from importlib.resources import files
from pathlib import Path
from .artifacts import MechanismArtifact
from .families import parameter_specs, PHYSICAL_FAMILIES
from .machine_basis import load_machine_basis, machine_parameters
from .schema_v2 import POLICIES
from .study_io import dumps


def initialize_v2(output, small='harmonic', large='harmonic', *, coupling='independent', champion=False):
    output=Path(output)
    if output.suffix!='.toml': raise ValueError('Study output must end in .toml.')
    resources=files('dada_solver.research').joinpath('data')
    seed_data=json.loads(resources.joinpath('kinematics_seeds.json').read_text())
    basis_file=output.with_suffix('.basis.json')
    artifact_paths={side:output.with_name(output.stem+'.'+side+'.mechanism.json') for side in ('small','large')}
    if any(p.exists() for p in (output,basis_file,*artifact_paths.values())): raise ValueError('Study or associated input already exists; choose a new output name.')
    if champion: small=large='structured_c2_15p'
    basis_text=resources.joinpath('structured3952_machine_basis.json' if champion else 'machine_basis_v2.json').read_text()
    digest=hashlib.sha256(basis_text.encode()).hexdigest()
    # Build and validate all content before publishing a complete template.
    import tempfile
    with tempfile.TemporaryDirectory() as directory:
        path=Path(directory)/'basis.json';path.write_text(basis_text)
        basis=load_machine_basis(path,digest)
    specs,defaults=machine_parameters(basis)
    declarations=[dict(name=k,value=v,unit=specs[k].unit) for k,v in defaults.items()]
    kinematics=dict(coupling=coupling); artifacts={}; mechanical=[]
    for side,family in (('small',small),('large',large)):
        if family not in seed_data['seeds'][side]: raise ValueError(f'Unknown family: {family}')
        seed=copy.deepcopy(seed_data['seeds'][side][family]); settings=seed['settings']; parameters=seed['parameters']
        owned=parameter_specs(settings,side)
        if family in PHYSICAL_FAMILIES:
            constraints=list(seed.get('constraints',[]))
            if family=='six_bar':
                constraints=[dict(metric=m,relation=r,limit=v,unit=u) for m,r,v,u in (
                    ('stroke_over_crank','minimum',1.,'crank_radius'),('stroke_over_crank','maximum',3.,'crank_radius'),
                    ('minimum_primary_transmission_sine','minimum',.30,'1'),('minimum_secondary_transmission_sine','minimum',.30,'1'),
                    ('minimum_rod_axis_cosine','minimum',.95,'1'),('EH_over_crank','maximum',7.,'crank_radius'),
                    ('crank_axis_to_EFH_clearance_over_crank','minimum',.5,'crank_radius'),
                    ('H_axis_lateral_rms_over_stroke','maximum',.25,'1'),('H_axis_lateral_span_over_stroke','maximum',.65,'1'),
                    ('zero_crossing_count','equal',2,'1'))]
            artifact=MechanismArtifact.create(family,parameters,settings=settings,constraints=constraints,
                provenance=dict(description='Historical geometry seed; no new optimization',sources=seed_data['provenance']))
            artifacts[side]=artifact
            settings=dict(settings,artifact=artifact_paths[side].name,sha256=artifact.content_hash)
        else:
            if family in ('free_spline','fourier_c2','structured_c2_15p'):
                mechanical.append(dict(side=side,metric='zero_crossing_count',relation='equal',limit=2,unit='1'))
        kinematics[side]=settings
        for name,value in parameters.items():
            if name not in owned: continue
            if coupling=='shared_crank' and name in ('phase_rad','crank_direction'):
                if side=='large': continue
                key=f'kinematics.shared.{name}'
            else: key=f'kinematics.{side}.{name}'
            declarations.append(dict(name=key,value=value,unit=owned[name].unit))
    raw=dict(schema_version=2,study=dict(name='Dada-Engine Research — '+('candidate 3952' if champion else small+' / '+large),
        protocol='machine_design',purpose='bounded_family_evaluation'),sources=dict(machine=dict(path=basis_file.name,sha256=digest)),
        kinematics=kinematics,parameters=declarations,policies=POLICIES,
        objective=dict(type='maximize_thermal_efficiency',unit='1'),
        constraints=[dict(type='minimum_motor_power',required_power=25.,unit='W'),dict(type='maximum_pressure',limit=1200000.,unit='Pa'),
            dict(type='maximum_temperature',limit=850.,unit='K'),dict(type='maximum_absolute_mass_flow',limit=.08,unit='kg/s'),dict(type='valid_thermodynamic_model',unit='1')],
        mechanical_constraints=mechanical,screening=dict(samples=1440,maximum_large_enclosed_volume_m3=.066),
        search=dict(type='sobol',seed=3965000,scramble=True,domain='fixed_global_bounds'),
        numerical=dict(maximum_cycles=100,backend='numba',candidate_budget_seconds=180.,domain_error_retry='once_safe_uniform_state_with_source_wall_temperatures'),
        warm_start=dict(initial_source='source_exact' if champion or (small==large=='six_bar') else 'uniform'),
        execution=dict(default_budget='2m',default_max_candidates=2,initial_evaluation_seconds=30.,deadline_grace_seconds=5.))
    source=dumps(raw)
    # Validate in isolation so bad family combinations cannot leave partial inputs.
    from .schema import load_study
    with tempfile.TemporaryDirectory() as directory:
        directory=Path(directory)
        (directory/output.name).write_text(source); (directory/basis_file.name).write_text(basis_text)
        for side,artifact in artifacts.items(): artifact.save(directory/artifact_paths[side].name)
        load_study(directory/output.name)
    output.parent.mkdir(parents=True,exist_ok=True)
    basis_file.write_text(basis_text)
    for side,artifact in artifacts.items(): artifact.save(artifact_paths[side])
    output.write_text(source)
    return output
