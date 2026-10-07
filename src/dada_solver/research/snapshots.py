"""Reconstruct stored scientific inputs without requiring execution resume."""
from contextlib import contextmanager
import copy
import hashlib
import json
from pathlib import Path
import tempfile

from .schema import load_study
from .study_io import dumps


@contextmanager
def stored_study(data):
    """Use verified campaign snapshots or the embedded evaluation definition.

    Source/runtime mismatch remains visible in inspection. This context only
    reconstructs inputs; callers may separately request a bounded report replay.
    """
    path = Path(data['source'])
    if path.is_dir():
        study = load_study(path/'study.toml', basis_path=path/'basis.json', artifact_directory=path)
        if study.study_id != data['study_id']:
            raise ValueError('Stored study snapshot differs from the inspected scientific identity.')
        yield study
        return
    raw = copy.deepcopy(data['scientific'])
    if raw['schema_version'] != 3:
        raise ValueError('Standalone reconstruction requires a schema-3 evaluation or a campaign directory.')
    basis = raw.pop('basis')
    fixed = raw.pop('fixed_parameters')
    raw['study']['name'] = data['name']
    raw['execution'] = dict(default_budget='3m', default_max_candidates=1,
        initial_evaluation_seconds=30., deadline_grace_seconds=5.)
    from .families import validate_settings
    from .artifacts import MechanismArtifact
    artifacts={}
    for side in ('small', 'large'):
        cfg = raw['kinematics'][side]
        if 'mechanism' in cfg:
            mechanism=cfg['mechanism']
            artifact=MechanismArtifact.create(mechanism['settings']['family'],mechanism['geometry'],
                settings=mechanism['settings'],constraints=mechanism['constraints'],
                provenance=dict(reconstructed_from_evaluation_study_id=data['study_id']))
            if 'sha256' in cfg and cfg['sha256']!=artifact.content_hash:
                raise ValueError('Embedded evaluation mechanism hash mismatch.')
            artifacts[side]=artifact
            raw['kinematics'][side]=dict(family=mechanism['settings']['family'],
                artifact=side+'.mechanism.json',sha256=artifact.content_hash)
        specs = validate_settings(artifacts[side].scientific['settings'] if side in artifacts else raw['kinematics'][side], side)
        declared = {r['name'] for r in raw['parameters']}
        for key, spec in specs.items():
            name = f'kinematics.{side}.{key}'
            if name in fixed and name not in declared:
                raw['parameters'].append(dict(name=name, value=fixed[name], unit=spec.unit))
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        text = json.dumps(basis, indent=2, allow_nan=False)+'\n'
        (root/'basis.json').write_text(text)
        for side,artifact in artifacts.items():
            (root/(side+'.mechanism.json')).write_text(json.dumps(artifact.data,indent=2)+'\n')
        raw['sources']['machine'] = dict(path='basis.json', sha256=hashlib.sha256(text.encode()).hexdigest())
        (root/'study.toml').write_text(dumps(raw))
        yield load_study(root/'study.toml')
