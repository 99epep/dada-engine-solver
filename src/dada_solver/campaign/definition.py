"""Frozen campaign inputs, source snapshots and explicit evaluator selection."""
from dataclasses import asdict, replace
from pathlib import Path
import hashlib
import json
import platform
import tomllib
import numpy as np
import scipy
import math
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
    def __init__(self, path, *, base_path=None, hardware_path=None):
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
        self.hardware_source = None
        self.hardware_data = self.hardware_bank = self.heat_in_inputs = self.heat_out_inputs = None
        if self.families['exchanger'] == 'microtube':
            configured = raw['campaign'].get('hardware_configuration')
            if configured is None:
                raise ValueError('Microtube campaigns require campaign.hardware_configuration.')
            selected = Path(hardware_path) if hardware_path is not None else self.path.parent/configured
            self.hardware_source = selected.read_text()
            from dada_solver.exchangers.hardware import load_hardware_definition
            (self.hardware_data, self.hardware_bank, self.heat_in_inputs,
             self.heat_out_inputs) = load_hardware_definition(
                self.hardware_source, self.configuration.gas.heat_capacity_cp)
            from dada_solver.exchangers.wall_cycle import WallCycleNumericalSettings
            wall = dict(raw.get('wall_numerical', {}))
            wall.setdefault('integration_method', self.configuration.numerical.integration_method)
            wall.setdefault('maximum_step_angle_radians',
                math.radians(self.configuration.numerical.maximum_step_angle_degrees))
            if 'integration_absolute_tolerances' in wall:
                wall['integration_absolute_tolerances'] = tuple(wall['integration_absolute_tolerances'])
            self.wall_numerical_settings = WallCycleNumericalSettings(**wall)
            self.numerical_settings = dict(self.numerical_settings,
                wall_cycle=asdict(self.wall_numerical_settings))
        else:
            self.wall_numerical_settings = None
        self.free_settings = raw.get('free', {})
        validate_ownership({p.name for p in self.space.parameters}, self.families, self.free_settings)
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
        self.adapter = FamilyDesignAdapter(self.configuration, self.families, self.free_settings,
            hardware_bank=self.hardware_bank, heat_in_inputs=self.heat_in_inputs,
            heat_out_inputs=self.heat_out_inputs)
        search = raw.get('search', {})
        if search.get('type', 'sobol') != 'sobol': raise ValueError('Only Sobol exploration is implemented.')
        self.seed, self.scramble = int(search.get('seed', 0)), bool(search.get('scramble', True))
        self.elite_size = int(raw['campaign'].get('elite_size', 5))
        self.initial_evaluation_seconds = float(raw['campaign'].get('initial_evaluation_seconds', 1))
        self.deadline_grace_seconds = float(raw['campaign'].get('deadline_grace_seconds', 12.0))
        self.maximum_candidates = raw['campaign'].get('maximum_candidates')
        if self.elite_size < 1 or not np.isfinite(self.initial_evaluation_seconds) or self.initial_evaluation_seconds <= 0:
            raise ValueError('Archive size and initial time estimate must be positive.')
        if self.maximum_candidates is not None and (not isinstance(self.maximum_candidates, int) or self.maximum_candidates < 1):
            raise ValueError('Maximum candidates must be a positive integer.')
        if not np.isfinite(self.deadline_grace_seconds) or self.deadline_grace_seconds < 0:
            raise ValueError('Deadline grace must be finite and nonnegative.')
        self.identity = dict(schema_version=2, campaign=raw, base_configuration=self.base_source,
                             hardware_configuration=self.hardware_source,
                             runtime=runtime_identity())
        self.definition_id = content_hash(self.identity)

    @classmethod
    def resume(cls, directory):
        directory = Path(directory)
        hardware = directory/'hardware.toml'
        definition = cls(directory/'campaign.toml', base_path=directory/'base.toml',
                         hardware_path=hardware if hardware.exists() else None)
        recorded = json.loads((directory/'definition.json').read_text())
        if recorded['definition_id'] != definition.definition_id:
            raise ValueError('Campaign definition, source code or runtime changed; create a new campaign directory.')
        return definition
