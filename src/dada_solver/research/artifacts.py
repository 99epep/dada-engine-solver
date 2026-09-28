"""Portable scientific mechanism artifacts, separate from drawings and search policy."""
from dataclasses import dataclass
import json
from pathlib import Path
from dada_solver.campaign.candidate import content_hash, canonical_json
from dada_solver.campaign.history import atomic_json
from dada_solver.geometry import CylinderVolumeLimits
from .families import PHYSICAL_FAMILIES, validate_settings, build_side


@dataclass(frozen=True)
class MechanismArtifact:
    payload_json: str

    @property
    def data(self): return json.loads(self.payload_json)
    @property
    def scientific(self): return self.data['scientific']
    @property
    def content_hash(self): return self.data['content_hash']

    @classmethod
    def create(cls, family, geometry, *, settings=None, constraints=(), provenance=None):
        settings=dict(settings or {},family=family)
        if family not in PHYSICAL_FAMILIES: raise ValueError('Artifact family must describe a physical mechanism.')
        if set(settings)&{'artifact','sha256'}: raise ValueError('Nested artifact references are not supported.')
        specs=validate_settings(settings,'small')
        if set(geometry)!=set(specs): raise ValueError('Artifact geometry must contain every family coordinate exactly once.')
        for name,value in geometry.items(): specs[name].validate(value)
        # Construct production closure eagerly. No thermodynamics or optimizer.
        build_side(settings,geometry,'small',CylinderVolumeLimits(1.,2.))
        from .margins import validate_mechanical_constraint
        for row in constraints: validate_mechanical_constraint(row,family,scoped=False)
        scientific=dict(schema_version=1,settings=settings,geometry=geometry,
            length_convention='crank_radius',angle_convention='study_radians',constraints=list(constraints))
        return cls(canonical_json(dict(schema_version=1,artifact_type='mechanism',scientific=scientific,
            content_hash=content_hash(scientific),provenance=provenance or {})))

    @classmethod
    def from_data(cls, data):
        if set(data)!={'schema_version','artifact_type','scientific','content_hash','provenance'} or type(data['schema_version']) is not int or data['schema_version']!=1 or data['artifact_type']!='mechanism':
            raise ValueError('Unsupported mechanism artifact schema.')
        raw=data['scientific']
        if set(raw)!={'schema_version','settings','geometry','length_convention','angle_convention','constraints'} or type(raw['schema_version']) is not int or raw['schema_version']!=1:
            raise ValueError('Unsupported scientific mechanism schema.')
        if raw['length_convention']!='crank_radius' or raw['angle_convention']!='study_radians':
            raise ValueError('Mechanism units or coordinate convention are unsupported.')
        result=cls.create(raw['settings']['family'],raw['geometry'],settings=raw['settings'],constraints=raw['constraints'],provenance=data['provenance'])
        if result.content_hash!=data['content_hash']: raise ValueError('Mechanism content hash mismatch.')
        return result

    @classmethod
    def load(cls, path, expected_hash=None):
        result=cls.from_data(json.loads(Path(path).read_text()))
        if expected_hash is not None and result.content_hash!=expected_hash: raise ValueError('Referenced mechanism content hash mismatch.')
        return result

    def save(self,path):
        path=Path(path)
        if path.exists(): raise ValueError('Mechanism artifact already exists; choose a new path.')
        path.parent.mkdir(parents=True,exist_ok=True)
        atomic_json(path,self.data)

    def reconstruct(self):
        raw=self.scientific
        return build_side(raw['settings'],raw['geometry'],'small',CylinderVolumeLimits(1.,2.))[1]


@dataclass(frozen=True)
class MechanismLibrary:
    """Retain several named families; no implicit winner or score-based pruning."""
    members: tuple

    def __post_init__(self):
        ids=[m['family_id'] for m in self.members]
        if not ids or len(ids)!=len(set(ids)): raise ValueError('Library requires unique nonempty family identifiers.')
        for member in self.members:
            if set(member)!={'family_id','mechanisms','metadata'} or not isinstance(member['family_id'],str) or not member['family_id']:
                raise ValueError('Invalid library member.')
            if not member['mechanisms']: raise ValueError('Each family needs at least one mechanism.')
            for raw in member['mechanisms'].values(): MechanismArtifact.from_data(raw)

    def to_data(self):
        scientific=dict(schema_version=1,artifact_type='mechanism_library',members=list(self.members))
        return dict(scientific,content_hash=content_hash(scientific))

    @classmethod
    def from_data(cls,data):
        if data.get('schema_version')!=1 or data.get('artifact_type')!='mechanism_library' or set(data)!={'schema_version','artifact_type','members','content_hash'}:
            raise ValueError('Unsupported mechanism library.')
        if content_hash({k:v for k,v in data.items() if k!='content_hash'})!=data['content_hash']: raise ValueError('Library hash mismatch.')
        return cls(tuple(data['members']))
