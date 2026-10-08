"""Portable Research declaration editing, independent of terminal UI and evaluation.

Domains are numerical search policies. The local scheduler radius is separate
from those domains and applies to every active coordinate, including hardware.
"""
from dataclasses import dataclass
import copy
import hashlib
import json
import math
from pathlib import Path
import tempfile

from dada_solver.campaign.parameters import parameter_from_mapping, ParameterSpace
from .artifacts import MechanismArtifact
from .families import parameter_specs, PHYSICAL_FAMILIES, ParameterSpec
from .hardware_retuning import GROUPS as HARDWARE_GROUPS, _effective_bounds
from .machine_basis import machine_parameters
from .schema import load_study, compile_study, candidate_for_values
from .study_io import dumps
from .synthesis_search import BOUNDS, BOUNDS_VERSION

POLICY = 'research_declaration_editor_v1'
GROUPS = ('mechanisms','kinematics-small','kinematics-large',*HARDWARE_GROUPS)


def freeze_declaration(row, *, value=None):
    return dict(name=row['name'],unit=row['unit'],
                value=(row['value'] if 'value' in row else row['initial']) if value is None else value)


def continuous_declaration(name, unit, initial, lower, upper, *, transform='linear'):
    return dict(name=name,unit=unit,kind='continuous',initial=initial,
                lower=lower,upper=upper,transform=transform)


def integer_declaration(name, unit, initial, lower, upper):
    return dict(name=name,unit=unit,kind='integer',initial=initial,
                lower=lower,upper=upper,encoding='nearest_even_v1')


def initial_center(rows):
    return {r['name']:r['initial'] for r in rows if 'initial' in r}


@dataclass(frozen=True)
class EditorParameter:
    name: str
    unit: str
    kind: str
    value: object
    active: bool
    source_of_value: str
    editable: bool
    releasable: bool
    reason: str = ''


class StudyEditor:
    """An in-memory edit session; only save writes, after full staged preflight."""
    def __init__(self, source):
        self.source = Path(source)
        if self.source.suffix != '.toml':
            raise ValueError('Study editing requires a study.toml, not a running campaign.')
        # A stored campaign snapshot must never be overwritten or edited in place.
        self.study = load_study(self.source)
        self.definition = compile_study(self.study)
        self.definition.adapter.build(dict(self.study.fixed_parameters,
                                          **initial_center(self.study.data['parameters'])))
        self.raw = copy.deepcopy(self.study.data)
        self.basis = self.study.basis.data
        self.specs, defaults = machine_parameters(self.study.basis,self.study.fixed_parameters.keys() | initial_center(self.raw['parameters']).keys())
        for side in ('small','large'):
            self.specs.update({f'kinematics.{side}.{n}':s for n,s in parameter_specs(self.study.settings[side],side).items()})
        if self.raw['kinematics']['coupling']=='shared_crank':
            for n,kind,unit in (('phase_rad','continuous','rad'),('crank_direction','branch','1')):
                for side in ('small','large'): self.specs.pop(f'kinematics.{side}.{n}')
                self.specs[f'kinematics.shared.{n}']=ParameterSpec(unit,kind)
        self.base_values = dict(defaults,**self.study.fixed_parameters,**initial_center(self.raw['parameters']))
        self.domains = {r['name']:copy.deepcopy(r) for r in self.raw['parameters'] if 'initial' in r}
        self.domain_sources = {n:'current declared domain' for n in self.domains}
        self.recenter = False
        self.search_request = None
        self.fixed_by_contract = {}

    def _rows(self):
        return {r['name']:r for r in self.raw['parameters']}

    def values(self):
        result = dict(self.base_values)
        for row in self.raw['parameters']:
            result[row['name']]=row.get('value',row.get('initial'))
        return result

    def inventory(self):
        rows=self._rows(); values=self.values(); result=[]
        for name,spec in self.specs.items():
            if name not in values: continue
            row=rows.get(name)
            side=name.split('.')[1] if name.startswith('kinematics.') else None
            mechanical=side in self.study.settings and self.study.settings[side]['family'] in PHYSICAL_FAMILIES
            category=spec.kind in ('branch','boolean') or mechanical and spec.kind=='choice'
            inert_charge=name=='charge.total_mass_kg' and self.raw['policies']['charge']!='explicit_inventory'
            releasable=spec.kind in ('continuous','integer','choice') and not category and not inert_charge
            source=('active declaration' if row and 'initial' in row else
                    'explicit fixed declaration' if row else
                    'mechanism artifact' if side in self.study.artifacts else 'machine basis')
            result.append(EditorParameter(name,spec.unit,spec.kind,values[name],bool(row and 'initial' in row),
                source,not category and not inert_charge,releasable,
                'Fixed scientific category' if category else
                'Charge follows the reference-inventory policy' if inert_charge else
                'Raw spline controls are fixed by the Research contract' if spec.kind=='fixed_continuous' else ''))
        return tuple(result)

    def select(self, *, groups=(), parameters=(), releasable_only=False):
        unknown=set(groups)-set(GROUPS)
        if unknown: raise ValueError(f'Unknown study groups: {sorted(unknown)}.')
        inventory={p.name:p for p in self.inventory()}
        if set(parameters)-set(inventory):
            raise ValueError(f'Unknown study parameters: {sorted(set(parameters)-set(inventory))}.')
        selected=set(parameters)
        for group in groups:
            if group in HARDWARE_GROUPS:
                selected.update(n for n in inventory if any(n.startswith(prefix) for prefix in HARDWARE_GROUPS[group]))
            else:
                sides=('small','large') if group=='mechanisms' else (group.split('-')[1],)
                for side in sides:
                    if group=='mechanisms' and self.study.settings[side]['family'] not in PHYSICAL_FAMILIES: continue
                    selected.update(n for n in inventory if n.startswith(f'kinematics.{side}.') or n.startswith('kinematics.shared.'))
        # Explicit selections never silently skip a scientific category.
        if releasable_only:
            locked=[n for n in parameters if not inventory[n].releasable]
            if locked: raise ValueError(f'Parameters must remain fixed under the current Research contract: {locked}.')
            selected={n for n in selected if inventory[n].releasable}
        return [p.name for p in self.inventory() if p.name in selected]

    def _replace(self, row):
        for i,current in enumerate(self.raw['parameters']):
            if current['name']==row['name']:
                self.raw['parameters'][i]=row; return
        self.raw['parameters'].append(row)

    def resolve_domain(self, name, *, explicit=None):
        """Current → recorded original → centered synthesis box; explicit edits override."""
        entry=next((p for p in self.inventory() if p.name==name),None)
        if entry is None: raise ValueError(f'Unknown study parameter {name}.')
        if not entry.releasable: raise ValueError(f'{name}: {entry.reason or "not releasable"}.')
        value=self.values()[name]; spec=self.specs[name]
        row=copy.deepcopy(explicit or self.domains.get(name))
        origin='explicit user domain' if explicit else self.domain_sources.get(name)
        provenance=self.basis['provenance']
        if row is None:
            # Editor descendants retain domains of coordinates subsequently frozen.
            recorded=provenance.get('study_edit',{}).get('declared_domains',{})
            if name in recorded:
                row=copy.deepcopy(recorded[name]); origin='recorded study editor domain'
        if row is None:
            saved=(provenance.get('paired_thermodynamic',{}).get('hardware_source_domain') or
                   provenance.get('hardware_retuning',{}).get('hardware_source_domain'))
            if saved:
                from dada_solver.campaign.candidate import content_hash
                if saved.get('content_hash') != content_hash({k:v for k,v in saved.items() if k!='content_hash'}):
                    raise ValueError('Recorded hardware source-domain hash mismatch.')
                row=next((copy.deepcopy(r) for r in saved['parameters'] if r['name']==name),None)
                if row: origin='recorded original hardware domain'
        if row is None and name in provenance.get('motion_refit',{}).get('search_regions',{}):
            box=provenance['motion_refit']['search_regions'][name]
            row=continuous_declaration(name,spec.unit,value,box['lower'],box['upper'])
            origin='recorded motion refit search region'
        if row is None and name.startswith('kinematics.'):
            _,side,coordinate=name.split('.')
            box=provenance.get('paired_thermodynamic',{}).get('reference_boxes',{}).get(side,{}).get(coordinate)
            if box is not None and box[0]<=value<=box[1]:
                row=continuous_declaration(name,spec.unit,value,*box); origin='paired thermodynamic reference box'
            else:
                family=self.study.settings['small']['family'] if side=='shared' else self.study.settings.get(side,{}).get('family')
                if family in BOUNDS and coordinate in BOUNDS[family]:
                    # Same centered reference-box policy as mechanism adapt. Never
                    # wrap angles or clamp a current initial into another interval.
                    low,high=BOUNDS[family][coordinate]; half=(high-low)/2
                    if spec.positive: half=min(half,value/2)
                    row=continuous_declaration(name,spec.unit,value,value-half,value+half)
                    origin=f'centered synthesis reference box ({BOUNDS_VERSION})'
        if row is None:
            raise ValueError(f'{name}: no reliable search domain. Provide explicit lower/upper bounds (or declared choices).')
        row.update(name=name,unit=spec.unit,initial=value)
        if row.get('kind')!=spec.kind:
            raise ValueError(f'{name}: recorded domain kind differs from current {spec.kind}.')
        try:
            parameter_from_mapping({k:v for k,v in row.items() if k!='unit'})
            for v in row.get('choices', [row.get('lower'),row.get('upper')]): spec.validate(v)
        except (ValueError,TypeError) as error:
            raise ValueError(f'{name}: {origin} cannot contain the exact current value {value!r}: {error}. Provide an explicit domain; the initial will not be clamped.') from error
        return row,origin

    def release(self, *, groups=(), parameters=(), domains=None):
        names=self.select(groups=groups,parameters=parameters,releasable_only=True)
        selected=self.select(groups=groups)
        self.fixed_by_contract.update({p.name:p.reason for p in self.inventory() if p.name in selected and not p.releasable})
        if not names: raise ValueError('No releasable parameters match the selection.')
        domains=domains or {}
        if set(domains)-set(names): raise ValueError('Explicit domains must belong to selected parameters.')
        # Resolve every domain before changing any declaration.
        planned=[(n,*self.resolve_domain(n,explicit=domains.get(n))) for n in names]
        for name,row,origin in planned:
            self._replace(row); self.domains[name]=copy.deepcopy(row); self.domain_sources[name]=origin
        return names

    def freeze(self, *, groups=(), parameters=()):
        names=self.select(groups=groups,parameters=parameters)
        for name in names:
            row=self._rows().get(name)
            if row and 'initial' in row:
                self.domains[name]=copy.deepcopy(row)
                self._replace(freeze_declaration(row))
        return names

    def edit_parameter(self, name, **changes):
        entry=next((p for p in self.inventory() if p.name==name),None)
        if entry is None or not entry.editable:
            raise ValueError(f'{name}: not editable under the current Research contract.')
        row=copy.deepcopy(self._rows().get(name,dict(name=name,unit=entry.unit,value=entry.value)))
        allowed={'initial','lower','upper','transform','choices'} if 'initial' in row else {'value'}
        if set(changes)-allowed: raise ValueError(f'{name}: editable fields are {sorted(allowed)}.')
        row.update(changes)
        self.specs[name].validate(row.get('value',row.get('initial')))
        if 'initial' in row:
            parameter_from_mapping({k:v for k,v in row.items() if k!='unit'})
            for v in row.get('choices',[row.get('lower'),row.get('upper')]):
                self.specs[name].validate(v)
            self.domains[name]=copy.deepcopy(row)
        self._replace(row)

    def configure_search(self, mode, *, radius=.05, recenter=False):
        if mode not in ('local','global'): raise ValueError('Choose local or global search.')
        if isinstance(radius,bool) or not isinstance(radius,(int,float)) or not math.isfinite(radius) or not 0<radius<=1:
            raise ValueError('Local scheduler radius must lie in (0, 1].')
        self.search_request=(mode,radius)
        self.recenter=bool(recenter)

    def edit_execution(self, *, budget=None, max_candidates=None, seed=None, scramble=None):
        if budget is not None:
            from dada_solver.campaign.runner import parse_budget
            if parse_budget(budget)<=0: raise ValueError('Execution budget must be positive.')
            self.raw['execution']['default_budget']=budget
        if max_candidates is not None:
            if type(max_candidates) is not int or max_candidates<1: raise ValueError('Maximum candidates must be a positive integer.')
            self.raw['execution']['default_max_candidates']=max_candidates
        if seed is not None:
            if type(seed) is not int or seed<0: raise ValueError('Sobol seed must be a nonnegative integer.')
            self.raw['search']['seed']=seed
        if scramble is not None:
            if type(scramble) is not bool: raise ValueError('Sobol scramble must be boolean.')
            self.raw['search']['scramble']=scramble

    def _search(self):
        from .local_search import validate_search
        from dada_solver.campaign.candidate import content_hash
        search=copy.deepcopy(self.raw['search']); center=initial_center(self.raw['parameters'])
        before=initial_center(self.study.data['parameters'])
        mode,radius=self.search_request or ('local' if search['domain']=='local_regions_v1' else 'global',search.get('radius_fraction',.05))
        if not center or mode=='global':
            if search['domain']=='fixed_global_bounds': return search,False
            return dict(type='sobol',domain='fixed_global_bounds',seed=search['seed'],scramble=search['scramble'],evaluate_initial=True),bool(search.get('regions'))
        rows=[r for r in self.raw['parameters'] if 'initial' in r]
        space=ParameterSpace(tuple(parameter_from_mapping({k:v for k,v in r.items() if k!='unit'}) for r in rows))
        preserve=search['domain']=='local_regions_v1' and center==before and not self.recenter
        if preserve:
            try: validate_search(dict(search,radius_fraction=radius,choice_scope='declared_choices'),space)
            except ValueError: preserve=False
        if preserve: return dict(search,radius_fraction=radius,choice_scope='declared_choices'),False
        if len(search.get('regions',[]))>1 and not self.recenter:
            raise ValueError('Editing invalidates multiple local regions. Explicitly confirm recentering (--recenter) to create one region at current initials.')
        source_initial=initial_center(self.study.data['parameters'])
        source_candidate=candidate_for_values(self.definition,source_initial)
        region=dict(id='study_initial_'+content_hash(center)[:12],source_study_id=self.study.study_id,
                    source_candidate_id=source_candidate.candidate_id,center=center)
        return dict(type='sobol',domain='local_regions_v1',seed=search['seed'],scramble=search['scramble'],
                    radius_fraction=radius,allocation='round_robin',evaluate_centers=True,
                    choice_scope='declared_choices',regions=[region]),bool(search.get('regions'))

    def review(self):
        before=initial_center(self.study.data['parameters']); after=initial_center(self.raw['parameters'])
        values=self.values(); search,replaced=self._search()
        changes={n:dict(before=self.base_values[n],after=v) for n,v in values.items() if v!=self.base_values.get(n)}
        changed_bounds={n:r for n,r in self._rows().items() if 'initial' in r and
                        {k:v for k,v in r.items() if k!='initial'}!={k:v for k,v in {x['name']:x for x in self.study.data['parameters']}.get(n,{}).items() if k!='initial'}}
        warnings=[]
        if not after: warnings.append('No active parameters: use evaluate rather than run.')
        if replaced: warnings.append('Local regions will be replaced by one unevaluated initial center.' if search['domain']=='local_regions_v1' else 'Local regions will be removed for full declared-bound search.')
        if self.raw['kinematics']['coupling']=='shared_crank': warnings.append('Shared-crank phase belongs to both cylinders; coupling is unchanged.')
        effective={}
        for r in self.raw['parameters']:
            if 'initial' not in r: continue
            if search['domain']=='local_regions_v1':
                effective[r['name']]=[_effective_bounds(r,region['center'][r['name']],search['radius_fraction']) for region in search['regions']]
            else: effective[r['name']]=[copy.deepcopy(r)]
            p=parameter_from_mapping({k:v for k,v in r.items() if k!='unit'})
            if r['kind']!='choice' and min(p.encode(r['initial']),1-p.encode(r['initial']))<.05:
                warnings.append(f"{r['name']}: initial within 5% of a declared boundary.")
        artifacts=self._artifacts()
        artifact_changes={s:dict(before=self.study.artifacts[s].content_hash,after=a.content_hash) for s,a in artifacts.items() if a.content_hash!=self.study.artifacts[s].content_hash}
        if artifact_changes: warnings.append('Fixed effective geometry changes require new mechanical artifact hashes.')
        return dict(source=str(self.source),source_study_id=self.study.study_id,policy=POLICY,
                    active_before=len(before),active_after=len(after),activated_parameters=[n for n in after if n not in before],
                    frozen_parameters=[n for n in before if n not in after],changed_initial_values=changes,
                    changed_bounds=changed_bounds,search=search,effective_bounds=effective,
                    replaced_regions=replaced,artifact_changes=artifact_changes,warnings=warnings,
                    fixed_by_contract=dict(self.fixed_by_contract),
                    execution_configuration_change={k:dict(before=self.study.data['execution'].get(k),after=v)
                        for k,v in self.raw['execution'].items() if v!=self.study.data['execution'].get(k)},
                    center_evidence='unevaluated study initials; source_candidate_id is the computed source initial configuration identity, not an evaluated result')

    def _artifacts(self):
        values=self.values(); active=initial_center(self.raw['parameters']); result={}
        shared=self.raw['kinematics']['coupling']=='shared_crank'
        for side,old in self.study.artifacts.items():
            geometry={n:values[f'kinematics.{"shared" if shared and n in ("phase_rad","crank_direction") else side}.{n}']
                      for n in parameter_specs(self.study.settings[side],side)}
            fixed_changed=any(v!=old.scientific['geometry'][n] and
                f'kinematics.{"shared" if shared and n in ("phase_rad","crank_direction") else side}.{n}' not in active for n,v in geometry.items())
            result[side]=old if not fixed_changed else MechanismArtifact.create(old.scientific['settings']['family'],geometry,
                settings=old.scientific['settings'],constraints=old.scientific['constraints'],
                provenance=dict(old.data['provenance'],study_edit=dict(policy=POLICY,source_study_id=self.study.study_id,parent_artifact_hash=old.content_hash)))
        return result

    def save(self, output):
        output=Path(output)
        if output.suffix!='.toml': raise ValueError('Edited study output must end in .toml.')
        if (output.parent/'definition.json').exists(): raise ValueError('Do not write edited studies into a running campaign directory.')
        review=self.review(); raw=copy.deepcopy(self.raw); basis=copy.deepcopy(self.basis)
        basis_path=output.with_suffix('.basis.json'); contents={}
        artifacts=self._artifacts()
        for side,a in artifacts.items():
            path=output.with_name(output.stem+'.'+side+'.mechanism.json')
            raw['kinematics'][side]['artifact']=path.name; raw['kinematics'][side]['sha256']=a.content_hash
            contents[path]=json.dumps(a.data,indent=2,allow_nan=False)+'\n'
        raw['search']=review['search']
        previous=basis['provenance'].get('study_edit')
        entry={k:copy.deepcopy(review[k]) for k in ('source_study_id','policy','activated_parameters','frozen_parameters',
            'changed_initial_values','changed_bounds','artifact_changes','center_evidence','execution_configuration_change')}
        entry['search_configuration_change']=dict(before=self.study.data['search'],after=raw['search'])
        entry['declared_domains']=copy.deepcopy(self.domains)
        entry['domain_sources']=copy.deepcopy(self.domain_sources)
        entry['previous_edits']=[] if previous is None else previous.get('previous_edits',[])+[{k:v for k,v in previous.items() if k not in ('previous_edits','declared_domains')}]
        basis['provenance']['study_edit']=entry
        basis_text=json.dumps(basis,indent=2,allow_nan=False)+'\n'
        raw['sources']['machine']=dict(path=basis_path.name,sha256=hashlib.sha256(basis_text.encode()).hexdigest())
        contents[basis_path]=basis_text; contents[output]=dumps(raw)
        if len({p.resolve() for p in contents})!=len(contents) or any(p.exists() or p.is_symlink() for p in contents):
            raise ValueError('Edited study and all associated inputs must be new, distinct paths; nothing will be overwritten.')
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for p,text in contents.items(): (root/p.name).write_text(text)
            generated=load_study(root/output.name)
            compile_study(generated).adapter.build(dict(generated.fixed_parameters,**initial_center(generated.data['parameters'])))
            # Validate before touching the destination. Publish TOML last; roll
            # back only files created by this operation on any write failure.
            output.parent.mkdir(parents=True,exist_ok=True)
            written=[]
            try:
                for p,text in contents.items():
                    with p.open('x') as stream:
                        written.append(p); stream.write(text)
            except BaseException:
                for p in reversed(written): p.unlink()
                raise
        return output
