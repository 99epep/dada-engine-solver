"""Frozen campaign inputs, source snapshots and explicit evaluator selection."""
from dataclasses import asdict, replace
from pathlib import Path
import hashlib
import json
import platform
import tomllib
import numpy as np
import scipy
from dada_solver.configuration import load_simulation_configuration
from dada_solver.campaign.parameters import ParameterSpace, ContinuousParameter
from dada_solver.campaign.candidate import content_hash
from dada_solver.campaign.adapters import validate_ownership, FamilyDesignAdapter
from dada_solver.sizing.configuration import _load_objective, _load_constraint


def runtime_identity():
    root = Path(__file__).resolve().parents[1]
    digest = hashlib.sha256()
    for path in sorted(root.rglob('*.py')):
        digest.update(str(path.relative_to(root)).encode()); digest.update(path.read_bytes())
    return dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                source_sha256=digest.hexdigest())


class CampaignDefinition:
    def __init__(self, path, *, base_path=None):
        self.path = Path(path)
        self.source = self.path.read_text()
        self.data = tomllib.loads(self.source)
        raw = self.data
        self.base_path = Path(base_path) if base_path is not None else self.path.parent/raw['campaign']['base_configuration']
        self.base_source = self.base_path.read_text()
        self.configuration = load_simulation_configuration(self.base_path)
        if raw.get('numerical'):
            self.configuration = replace(self.configuration,
                numerical=replace(self.configuration.numerical, **raw['numerical']))
        self.numerical_settings = asdict(self.configuration.numerical)
        self.space = ParameterSpace(tuple(ContinuousParameter(**p) for p in raw['parameters']))
        self.families = raw['families']
        self.free_settings = raw.get('free', {})
        validate_ownership({p.name for p in self.space.parameters}, self.families, self.free_settings)
        # Keep the first production campaign on the existing eight-state evaluator.
        # One-wall hardware campaigns need a compatible performance adapter, not a
        # silently substituted static UA model or new convergence equations here.
        if self.families['exchanger'] != 'reservoir':
            raise ValueError('This campaign evaluator supports reservoir closures; dynamic-wall campaigns require an evaluator adapter.')
        if self.families['kinematics'] == 'free':
            from dada_solver.free_kinematics import FreeMotionDefinition
            for side in ('small','large'):
                limits = getattr(self.configuration.machine_volumes, side+'_cylinder')
                FreeMotionDefinition.from_shape_coordinates(self.free_settings[side+'_coordinates'],
                    limits.minimum, limits.maximum, **self.free_settings.get(side+'_limits', {}))
        self.objective = _load_objective(raw['objective'])
        self.constraints = tuple(_load_constraint(x, None) for x in raw.get('constraints', []))
        if len({x.name for x in self.constraints}) != len(self.constraints):
            raise ValueError('Constraint names must be unique.')
        self.adapter = FamilyDesignAdapter(self.configuration, self.families, self.free_settings)
        search = raw.get('search', {})
        if search.get('type', 'sobol') != 'sobol': raise ValueError('Only Sobol exploration is implemented.')
        self.seed, self.scramble = int(search.get('seed', 0)), bool(search.get('scramble', True))
        self.elite_size = int(raw['campaign'].get('elite_size', 5))
        self.initial_evaluation_seconds = float(raw['campaign'].get('initial_evaluation_seconds', 1))
        self.maximum_candidates = raw['campaign'].get('maximum_candidates')
        if self.elite_size < 1 or not np.isfinite(self.initial_evaluation_seconds) or self.initial_evaluation_seconds <= 0:
            raise ValueError('Archive size and initial time estimate must be positive.')
        if self.maximum_candidates is not None and (not isinstance(self.maximum_candidates, int) or self.maximum_candidates < 1):
            raise ValueError('Maximum candidates must be a positive integer.')
        self.identity = dict(schema_version=1, campaign=raw, base_configuration=self.base_source,
                             runtime=runtime_identity())
        self.definition_id = content_hash(self.identity)

    @classmethod
    def resume(cls, directory):
        directory = Path(directory)
        definition = cls(directory/'campaign.toml', base_path=directory/'base.toml')
        recorded = json.loads((directory/'definition.json').read_text())
        if recorded['definition_id'] != definition.definition_id:
            raise ValueError('Campaign definition, source code or runtime changed; create a new campaign directory.')
        return definition
