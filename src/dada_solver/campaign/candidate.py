"""Content-addressed candidates with deterministic portable JSON payloads."""
from dataclasses import dataclass
import hashlib
import json


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def content_hash(value):
    return hashlib.sha256(canonical_json(value).encode('utf-8')).hexdigest()


@dataclass(frozen=True, slots=True)
class Candidate:
    """Immutable canonical payload; no randomized Python hash or mutable mappings."""
    payload_json: str
    candidate_id: str

    def __post_init__(self):
        if content_hash(json.loads(self.payload_json)) != self.candidate_id:
            raise ValueError('Candidate payload does not match its SHA-256 identity.')

    @classmethod
    def create(cls, space, coordinates, *, families, numerical_settings, definition_id):
        coordinates = tuple(0.0 if float(x) == 0 else float(x) for x in coordinates)
        payload = dict(schema_version=1, definition_id=definition_id,
            normalized=coordinates, physical=space.decode(coordinates),
            families=families, numerical_settings=numerical_settings)
        return cls(canonical_json(payload), content_hash(payload))

    @property
    def payload(self):
        return json.loads(self.payload_json)
