"""Current parsers reject superseded syntax without migration or ignored keys."""
import importlib
import json

import pytest

from dada_solver.configuration import ValidityThresholds, load_simulation_configuration
from dada_solver.campaign.history import CampaignHistory, read_journal
from dada_solver.research.cli import main
from dada_solver.research.presets import initialize_v2
from dada_solver.research.schema import load_study
from dada_solver.research.study_io import dumps
from dada_solver.sizing.configuration import load_sizing_problem
from tests.synthetic_machine import CONFIGURATION


@pytest.mark.parametrize('key', ['maximum_pressure_equalization_error', 'maximum_mach_number'])
def test_configuration_rejects_unowned_validity_keys(tmp_path, key):
    with pytest.raises(TypeError):
        ValidityThresholds(.01, .01, **{key: .2})
    path = tmp_path / 'machine.toml'
    path.write_text(CONFIGURATION.read_text().replace('[validity]', f'[validity]\n{key} = 0.2'))
    with pytest.raises(ValueError, match='Unknown validity keys'):
        load_simulation_configuration(path)


@pytest.mark.parametrize('location', ['constraints', 'scales'])
def test_sizing_rejects_unowned_pressure_controls(tmp_path, location):
    raw = dict(problem=dict(base_configuration=str(CONFIGURATION)),
        variables=[dict(parameter='cold_ua', lower_bound=1., upper_bound=10., initial_value=5.)],
        objective=dict(type='minimize_total_ua'), constraints=[],
        optimizer=dict(objective_scale=1., unavailable_objective_penalty=1e6,
            unavailable_constraint_margin=-1., maximum_iterations=10,
            function_tolerance=1e-6, constraint_scales={}))
    if location == 'constraints':
        raw['constraints'].append(dict(type='maximum_pressure_equalization_error', limit=.1))
    else:
        raw['optimizer']['constraint_scales']['maximum_pressure_equalization_error'] = .1
    path = tmp_path / 'sizing.toml'
    path.write_text(dumps(raw))
    with pytest.raises(ValueError, match='Unsupported sizing constraint|Unknown constraint scales'):
        load_sizing_problem(path)


@pytest.mark.parametrize('version', [1, 2])
def test_research_rejects_superseded_schema(tmp_path, version):
    path = initialize_v2(tmp_path / 'study.toml')
    raw = load_study(path).data
    raw['schema_version'] = version
    path.write_text(dumps(raw))
    with pytest.raises(ValueError, match='Only study schema_version = 3'):
        load_study(path)


@pytest.mark.parametrize('module,name', [
    ('dada_solver.kinematics', 'VolumeKinematics'),
    ('dada_solver.four_bar', 'FourBarKinematics'),
    ('dada_solver.six_bar', 'load_six_bar_mechanism'),
    ('dada_solver.exchangers.air_wall', 'AirWallMotor'),
    ('dada_solver.exchangers.wall_cycle', 'solve_periodic_wall_motor'),
    ('dada_solver.exchangers.hardware', 'connect_hardware'),
])
def test_public_modules_expose_no_superseded_facades(module, name):
    assert not hasattr(importlib.import_module(module), name)


def test_research_cli_rejects_reference_selector(tmp_path):
    with pytest.raises(SystemExit) as caught:
        main(['evaluate', 'study.toml', '--output', str(tmp_path / 'result.json'), '--reference'])
    assert caught.value.code == 2


def test_uncompressed_journal_is_rejected(tmp_path):
    path = tmp_path / 'history.jsonl'
    path.write_text('')
    with pytest.raises(ValueError, match='Uncompressed'):
        CampaignHistory(tmp_path)
    with pytest.raises(ValueError, match='compressed format'):
        read_journal(path)


def test_per_candidate_recovery_layout_is_rejected(tmp_path):
    path = tmp_path / 'candidates'
    path.mkdir()
    (path / 'record.json').write_text('{}')
    with pytest.raises(ValueError, match='Per-candidate recovery'):
        CampaignHistory(tmp_path)


@pytest.mark.parametrize('kind', ['research_v1', 'research_v2'])
def test_resume_rejects_superseded_research_definitions(tmp_path, kind):
    from dada_solver.campaign.definition import CampaignDefinition
    (tmp_path / 'definition.json').write_text(json.dumps(dict(definition_kind=kind)))
    with pytest.raises(ValueError, match='current research_v3'):
        CampaignDefinition.resume(tmp_path)
