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
    basis_text=resources.joinpath('structured3952_machine_basis.json' if champion else 'hybrid_compact_machine_basis.json' if small==large=='hybrid_compact' else 'machine_basis_v2.json').read_text()
    basis_data=json.loads(basis_text)
    for key in ('maximum_pressure_equalization_error','maximum_mach_number'):
        basis_data['configuration']['validity'].pop(key,None)
    basis_text=json.dumps(basis_data,indent=2)+'\n'
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
            mechanical.extend(dict(row,side=side) for row in constraints)
            artifact=MechanismArtifact.create(family,parameters,settings=settings,constraints=[],
                provenance=dict(description='Historical geometry seed; no new optimization',sources=seed_data['provenance']))
            artifacts[side]=artifact
            settings=dict(settings,artifact=artifact_paths[side].name,sha256=artifact.content_hash)
        else:
            if family in ('free_spline','fourier_c2','structured_c2_15p','hybrid_compact'):
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
        constraints=[dict(type='valid_thermodynamic_model',unit='1')],
        mechanical_constraints=mechanical,screening=dict(samples=1440),
        search=dict(type='sobol',seed=3965000,scramble=True,domain='fixed_global_bounds'),
        numerical=dict(maximum_cycles=100,backend='numba',candidate_budget_seconds=180.,domain_error_retry='once_safe_uniform_state_with_source_wall_temperatures'),
        warm_start=dict(initial_source='source_exact' if champion or (small==large and small in ('six_bar','hybrid_compact')) else 'uniform'),
        execution=dict(default_budget='2m',default_max_candidates=512,initial_evaluation_seconds=30.,deadline_grace_seconds=5.))
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


def initialize_v3(output, small='harmonic', large='harmonic', *, mode='refrigeration'):
    """Bounded thermal-boundary fixture, not a human-cell optimized design.

    Liquid Cp and conductance are explicit scenario inputs, not correlations.
    The selected V2 motion coordinates remain separately editable.
    """
    import tempfile
    import tomllib
    from .schema_v2 import POLICIES_V3
    from .sixbar import exchanger_from_data
    from .schema import load_study
    output=Path(output)
    basis_file=output.with_suffix('.basis.json')
    if output.suffix!='.toml' or output.exists() or basis_file.exists():
        raise ValueError('Choose an unused study.toml and associated basis filename.')
    if mode not in ('motor','refrigeration'): raise ValueError('Unknown operating mode.')
    with tempfile.TemporaryDirectory() as temporary:
        directory=Path(temporary); path=directory/output.name
        initialize_v2(path,small,large)
        raw=tomllib.loads(path.read_text()); basis=json.loads(path.with_suffix('.basis.json').read_text())
        raw['schema_version']=basis['schema_version']=3
        raw['study'].update(name='External-stream '+mode+' validation',purpose='bounded_thermal_boundary_validation')
        raw['policies']=POLICIES_V3.copy();raw['policies']['outlet_valve_cda']='fixed_source_cda'
        raw['warm_start']['initial_source']='uniform';basis['warm_start']=None
        c=basis['configuration']
        if mode=='refrigeration':
            c['angular_speed']=2*3.141592653589793*.2
            # A benign validation compression ratio keeps startup inside 200–1000 K.
            for side in ('small_cylinder','large_cylinder'):
                limits=c['machine_volumes'][side];swept=limits['maximum']-limits['minimum']
                limits.update(minimum=.3*swept,maximum=1.3*swept)
            c['cold_reservoir_temperature']=278.15;c['hot_reservoir_temperature']=298.15
            c['heat_in_valve_placement']=c['heat_out_valve_placement']='downstream'
            c['charge']['temperature']=288.15
            raw['objective']=dict(type='maximize_cooling_cop',unit='1')
            raw['constraints']=[row for row in raw['constraints'] if row['type']!='minimum_motor_power']
            raw['constraints'].append(dict(type='maximum_mechanical_input_power',limit=100.,unit='W'))
        for side in ('heat_in','heat_out'):
            old=basis[side];thermal=exchanger_from_data(old).build().wall_thermal
            inputs={k:v for k,v in old['inputs'].items() if not k.startswith('air_') and k!='fan_total_efficiency'}
            inputs['external_stream']=dict(fluid='declared water-like liquid scenario',
                inlet_temperature_k=(278.15 if side=='heat_in' else 298.15) if mode=='refrigeration' else thermal.air_inlet_temperature_k,
                mass_flow_kg_s=.05,cp_j_kg_k=4180.,wall_conductance_w_k=150.)
            old['family']='external_stream_wall';old['inputs']=inputs
        # Rebuild owned defaults, retaining all selected V2 motion coordinates.
        text=json.dumps(basis,indent=2)+'\n';digest=hashlib.sha256(text.encode()).hexdigest()
        path.with_suffix('.basis.json').write_text(text)
        specs,defaults=machine_parameters(load_machine_basis(path.with_suffix('.basis.json'),digest))
        raw['parameters']=[dict(name=k,value=v,unit=specs[k].unit) for k,v in defaults.items()]+[
            row for row in raw['parameters'] if row['name'].startswith('kinematics.')]
        if mode=='refrigeration' and small==large=='harmonic':
            for row in raw['parameters']:
                if row['name']=='kinematics.small.phase_rad': row['value']=-3.141592653589793/2
        raw['sources']['machine']['sha256']=digest
        path.write_text(dumps(raw));load_study(path)
        output.parent.mkdir(parents=True,exist_ok=True)
        # Mechanism artifacts, when selected, are portable siblings of the study.
        for child in directory.iterdir():
            target=output.parent/child.name
            if target.exists(): raise ValueError('Associated input already exists: '+str(target))
        for child in directory.iterdir(): (output.parent/child.name).write_bytes(child.read_bytes())
    return output
