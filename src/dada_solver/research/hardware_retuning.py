"""Local hardware studies around a frozen, thermodynamically adapted mechanism."""
import copy
import hashlib
import json
import math
from pathlib import Path
import tempfile

from dada_solver.campaign.candidate import content_hash
from dada_solver.campaign.parameters import parameter_from_mapping, ChoiceParameter
from .artifacts import MechanismArtifact
from .families import PHYSICAL_FAMILIES, parameter_specs
from .machine_basis import machine_parameters
from .mechanism_adaptation import _source_warm_start
from .motion_target import research_motion_source
from .schema import load_study, compile_study
from .study_io import dumps

GROUPS = {
    'exchangers': ('microtube.', 'thermal.'),
    'volumes': ('volume.',),
    'frequency': ('operation.frequency_hz',),
    'charge': ('charge.',),
    'valves': ('valve.',),
    'external-stream': ('external_stream.',),
}
POLICY = 'source_active_hardware_local_regions_v1'


def _domain(provenance, source_study):
    paired = provenance['paired_thermodynamic']
    origins=[paired['source']]
    if 'source' in provenance.get('motion_refit',{}):
        origins.append(provenance['motion_refit']['source'])
    saved = paired.get('hardware_source_domain')
    if saved is not None:
        if set(saved)!={'source','parameters','content_hash'} or content_hash({k:v for k,v in saved.items() if k!='content_hash'})!=saved['content_hash']:
            raise ValueError('Hardware source-domain content hash mismatch.')
        if saved['source'] not in origins:
            raise ValueError('Hardware domain and paired source identities disagree.')
    if source_study is not None:
        selector=None if Path(source_study).suffix=='.toml' else (saved or paired)['source']['candidate_id']
        with research_motion_source(source_study,candidate=selector) as (study,_,identity):
            origin=next((o for o in origins if o['study_id']==identity['study_id']),None)
            if origin is None:
                raise ValueError('Explicit hardware source must match the original paired source study identity.')
            rows=[copy.deepcopy(r) for r in study.data['parameters'] if 'initial' in r and not r['name'].startswith('kinematics.')]
        if saved is not None and saved['parameters'] and (origin!=saved['source'] or rows!=saved['parameters']):
            raise ValueError('Explicit source differs from the recorded hardware domain.')
        saved=dict(source=origin,parameters=rows)
        saved['content_hash']=content_hash(saved)
    if saved is None:
        raise ValueError('The paired study has no portable hardware domains. Provide --source-study with its original Research source; no bounds will be invented.')
    rows=saved['parameters']
    if (not isinstance(rows,list) or len({r['name'] for r in rows})!=len(rows)
            or any('initial' not in r or r['name'].startswith('kinematics.') for r in rows)):
        raise ValueError('Hardware domains must be unique active non-kinematic Research declarations.')
    for row in rows:
        parameter_from_mapping({k:v for k,v in row.items() if k!='unit'})
    return copy.deepcopy(saved)


def _effective_bounds(row, initial, radius):
    p=parameter_from_mapping({k:v for k,v in row.items() if k!='unit'})
    center=p.encode(initial)
    if isinstance(p,ChoiceParameter):
        # Category order is an encoding, not a geometric distance.
        return dict(kind='choice',initial=initial,choices=list(p.choices),policy='declared_choices')
    low,high=max(0.,center-radius),min(1.,center+radius)
    return dict(kind=row['kind'],initial=initial,normalized_center=center,
                normalized_lower=low,normalized_upper=high,lower=p.decode(low),upper=p.decode(high))


def hardware_retuning(source, output, *, candidate=None, scope=None, groups=(), parameters=(),
                      radius=.1, source_study=None, expected_family=None):
    """Generate a portable study using original domains and exact current values.

    Original parameter bounds/encodings remain declarations; local Sobol clips
    their normalized intervals around the current candidate. Declared choices
    stay categorical and all are eligible. Only originally active hardware can
    be reopened. No scaling, capacity transformation or integration takes place.
    """
    if isinstance(radius,bool) or not isinstance(radius,(int,float)) or not math.isfinite(radius) or not 0<radius<=1:
        raise ValueError('Hardware retuning radius must lie in (0, 1].')
    if scope not in (None,'source-active') or scope is not None and (groups or parameters):
        raise ValueError('Choose --scope source-active or targeted --group/--parameter selections, not both.')
    if set(groups)-set(GROUPS):
        raise ValueError(f'Unknown hardware groups: {sorted(set(groups)-set(GROUPS))}.')
    output=Path(output)
    if output.suffix!='.toml': raise ValueError('Hardware retuning output must end in .toml.')
    basis_path=output.with_suffix('.basis.json')
    artifact_paths={s:output.with_name(output.stem+'.'+s+'.mechanism.json') for s in ('small','large')}
    paths=[output,basis_path,*artifact_paths.values()]
    if len({p.resolve() for p in paths})!=len(paths) or any(p.exists() for p in paths):
        raise ValueError('Retuning study and associated inputs must be distinct and must not exist.')
    with research_motion_source(source,candidate=candidate) as (study,values,identity):
        basis=study.basis.data
        provenance=basis['provenance']
        if provenance.get('paired_thermodynamic',{}).get('stage')!='paired_thermodynamic':
            raise ValueError('Hardware retuning requires a paired_thermodynamic study/candidate or its retuning descendant.')
        if study.data['kinematics']['coupling']!='independent' or any(
                settings['family'] not in PHYSICAL_FAMILIES for settings in study.settings.values()):
            raise ValueError('Hardware retuning requires two independent complete physical mechanisms.')
        if expected_family is not None and any(s['family']!=expected_family for s in study.settings.values()):
            raise ValueError('Source mechanisms do not match the requested synthesis family.')
        domain=_domain(provenance,source_study)
        raw=copy.deepcopy(study.data)
        physical=dict(study.fixed_parameters,**values)
        design=compile_study(study).adapter.build(physical)
        machine_specs,_=machine_parameters(study.basis,physical)
        available={r['name']:r for r in domain['parameters'] if r['name'] in physical and r['name'] in machine_specs}
        excluded=[r['name'] for r in domain['parameters'] if r['name'] not in available]
        targeted=bool(groups or parameters)
        selected=set(available) if not targeted else set(parameters)|{
            n for n in available if any(n.startswith(prefix) for g in groups for prefix in GROUPS[g])}
        if set(parameters)-set(available):
            raise ValueError(f'Parameters must be currently available and originally active hardware: {sorted(set(parameters)-set(available))}.')
        if not selected:
            raise ValueError('No originally active hardware parameters match this selection; fixed source parameters cannot be reopened.')
        rows=[]; centers={}; effective={}; artifacts={}; links={}
        for name,value in physical.items():
            if name.startswith('kinematics.'): continue
            if name not in selected:
                rows.append(dict(name=name,value=value,unit=machine_specs[name].unit))
        # Preserve source active-coordinate ordering, including Sobol ownership.
        for original in domain['parameters']:
            name=original['name']
            if name not in selected: continue
            value=physical[name]
            row=copy.deepcopy(original); row['initial']=value
            effective[name]=_effective_bounds(row,value,radius)
            rows.append(row); centers[name]=value
        for side in ('small','large'):
            cfg=study.settings[side]
            geometry={name:physical[f'kinematics.{side}.{name}'] for name in parameter_specs(cfg,side)}
            previous=study.artifacts.get(side)
            if previous is None:
                # Standalone evaluation snapshots contain scientific mechanics,
                # but need not retain a filesystem artifact reference.
                constraints=[{k:v for k,v in r.items() if k!='side'} for r in study.mechanical_constraints if r['side']==side]
                parent_hash=None
            else:
                constraints=previous.scientific['constraints']; parent_hash=previous.content_hash
            if previous is not None and previous.scientific['geometry']==geometry:
                artifact=previous
            else:
                artifact=MechanismArtifact.create(cfg['family'],geometry,settings=cfg,constraints=constraints,
                    provenance=dict(previous.data['provenance'] if previous is not None else {},
                        frozen_from_candidate=identity,stage='hardware_retuning',parent_artifact_hash=parent_hash))
            artifacts[side]=artifact
            links[side]=dict(hash=artifact.content_hash,parent_artifact_hash=parent_hash)
            raw['kinematics'][side]=dict(family=cfg['family'],artifact=artifact_paths[side].name,sha256=artifact.content_hash)
        raw['parameters']=rows
        raw['mechanical_constraints']=copy.deepcopy(list(study.mechanical_constraints))
        raw['study']['name']+=' — hardware retuning'
        raw['study']['parent_candidate_id']=identity['candidate_id']
        raw['search']=dict(type='sobol',domain='local_regions_v1',seed=raw['search']['seed'],scramble=raw['search']['scramble'],
            radius_fraction=radius,allocation='round_robin',evaluate_centers=True,choice_scope='declared_choices',
            regions=[dict(id='hardware_retuning',source_candidate_id=identity['candidate_id'],
                          source_study_id=identity['study_id'],center=centers)])
        warm=_source_warm_start(source,identity,basis,design,raw)
        previous_retuning=provenance.get('hardware_retuning')
        ancestors=[] if previous_retuning is None else previous_retuning.get('ancestors',[])+[
            {k:previous_retuning[k] for k in ('source','selected_parameters','radius_fraction','effective_bounds','mechanisms')}]
        entry=dict(stage='hardware_retuning',policy=POLICY,source=identity,
            abstract_source=provenance['paired_thermodynamic']['source'],
            paired_source=identity if previous_retuning is None else previous_retuning['paired_source'],
            mechanisms=links,hardware_source_domain=domain,selected_parameters=list(centers),
            selection=dict(scope='source-active' if not targeted else 'targeted',groups=list(groups),parameters=list(parameters)),
            excluded_unavailable_parameters=excluded,radius_fraction=radius,effective_bounds=effective,
            warm_start=warm,ancestors=ancestors)
        if Path(source).suffix!='.toml':
            from .report import inspect,select_records
            record,=select_records(inspect(source),[identity['candidate_id']])
            entry['source_assessment']={k:record.get(k) for k in ('objective','metrics','status')}
        # Retuning descendants retain the complete original parameter domain.
        provenance['paired_thermodynamic']['hardware_source_domain']=domain
        provenance['hardware_retuning']=entry
        basis_text=json.dumps(basis,indent=2,allow_nan=False)+'\n'
        raw['sources']['machine']=dict(path=basis_path.name,sha256=hashlib.sha256(basis_text.encode()).hexdigest())
        contents={basis_path:basis_text,**{artifact_paths[s]:json.dumps(a.data,indent=2)+'\n' for s,a in artifacts.items()},output:dumps(raw)}
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for path,text in contents.items(): (root/path.name).write_text(text)
            generated=load_study(root/output.name)
            compile_study(generated).adapter.build(dict(generated.fixed_parameters,**centers))
        for path,text in contents.items():
            path.parent.mkdir(parents=True,exist_ok=True)
            with path.open('x') as stream: stream.write(text)
    return output
