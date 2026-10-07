"""Strict current study schema and exact candidate identities."""
from pathlib import Path
import tomllib
import math
from dada_solver.campaign.candidate import Candidate, canonical_json, content_hash


def keys(data, required, label, optional=()):
    if not isinstance(data, dict): raise ValueError(f'{label} must be a table.')
    missing, extra = set(required)-set(data), set(data)-set(required)-set(optional)
    if missing or extra:
        raise ValueError(f'{label}: missing keys {sorted(missing)}; unknown keys {sorted(extra)}.')


def positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f'{name} must be finite and positive.')


def load_study(path, *, basis_path=None, artifact_directory=None):
    path = Path(path)
    from .validation_errors import validation_location
    with validation_location(path):
        raw = tomllib.loads(path.read_text())
        with validation_location(path,('schema_version',)):
            if type(raw.get('schema_version')) is not int or raw['schema_version'] != 3:
                raise ValueError('Only study schema_version = 3 is supported.')
        from .study_schema import load_current_study
        return load_current_study(path, basis_path=basis_path, artifact_directory=artifact_directory)


def compile_study(study):
    from .study_schema import ResearchDefinition
    return ResearchDefinition(study)


def candidate_for_values(definition, values):
    """Keep requested physical values exact; encoding is for identity/distance only."""
    coordinates = definition.space.encode(values)
    payload = dict(schema_version=1, definition_id=definition.definition_id,
                   normalized=coordinates, physical=dict(values), families=definition.families,
                   numerical_settings=definition.numerical_settings)
    return Candidate(canonical_json(payload), content_hash(payload))
