"""Validation, ownership and exact integer identities for the first study protocol."""
import hashlib
import json
from pathlib import Path

import pytest

from dada_solver.campaign.parameters import IntegerParameter
from dada_solver.research.cli import initialize
from dada_solver.research.schema import load_study, compile_study, candidate_for_values


@pytest.fixture
def study_path(tmp_path):
    return initialize(tmp_path/'study.toml')


def test_integer_decoding_ties_and_endpoints():
    p = IntegerParameter('count', 2, 6, 4)
    assert [p.decode(u) for u in (0,.125,.375,.625,.875,1)] == [2,2,4,4,6,6]
    for k in range(2,7): assert p.decode(p.encode(k)) == k
    for value in (True,3.0,0,7):
        with pytest.raises(ValueError): p.encode(value)
    for value in (-.01,1.01,float('nan')):
        with pytest.raises(ValueError): p.decode(value)


@pytest.mark.parametrize('values', [(True,6,4),(2.,6,4),(0,6,4),(2,2,2),(2,6,7)])
def test_invalid_integer_parameter(values):
    with pytest.raises(ValueError): IntegerParameter('count',*values)


def test_initialization_is_portable_and_does_not_overwrite(study_path):
    source = study_path.read_bytes()
    with pytest.raises(ValueError,match='already exists'): initialize(study_path)
    assert study_path.read_bytes() == source
    study = load_study(study_path)
    assert len(study.space.parameters) == 5
    assert study.basis.heat_in.inputs.air_inlet_temperature_k == 558.15
    assert study.basis.configuration.hot_reservoir_temperature == 598.15


@pytest.mark.parametrize('before,after,message', [
    ('schema_version = 1','schema_version = 2','schema_version'),
    ('unit = "m"','unit = "mm"','requires'),
    ('lower = 2400','lower = 2400.0','integers'),
    ('name = "volume.swept_ratio"','name = "common.hot_ua"','derived parameter'),
    ('hot_air_inlet_K = 558.15','hot_air_inlet_K = 598.15','differs'),
    ('mechanical_losses = "unknown"','mechanical_losses = "zero"','policies'),
    ('domain = "fixed_global_bounds"','domain = "local_refinement"','bounded Sobol'),
    ('unit = "W"','unit = "kW"','requires unit'),
    ('maximum_cycles = 100','maximum_cycles = true','positive integer'),
    ('schema_version = 1','schema_version = 1\nsurprise = true','unknown keys'),
])
def test_rejects_ambiguous_or_unsupported_studies(study_path,before,after,message):
    study_path.write_text(study_path.read_text().replace(before,after,1))
    with pytest.raises(ValueError,match=message): load_study(study_path)


def test_hash_changes_fail_before_geometry_or_integration(study_path):
    basis = study_path.with_suffix('.basis.json')
    basis.write_text(basis.read_text()+' ')
    with pytest.raises(ValueError,match='SHA-256'): load_study(study_path)


def test_invalid_geometry_is_identified_before_solver(study_path,monkeypatch):
    basis = study_path.with_suffix('.basis.json')
    old_hash = hashlib.sha256(basis.read_bytes()).hexdigest()
    raw = json.loads(basis.read_text()); raw['small']['primary_ground'] = 100
    basis.write_text(json.dumps(raw))
    study_path.write_text(study_path.read_text().replace(old_hash,hashlib.sha256(basis.read_bytes()).hexdigest()))
    monkeypatch.setattr('dada_solver.campaign.evaluator.solve_periodic_wall_motor',lambda *a,**k:pytest.fail('integration started'))
    with pytest.raises(ValueError,match='six-bar geometry'): load_study(study_path)


def test_presentation_path_and_invocation_budget_do_not_change_study_identity(study_path,tmp_path):
    first = load_study(study_path)
    other = initialize(tmp_path/'other.toml')
    other.write_text(other.read_text().replace('name = "Dada-Engine Research"','name = "A different label"').replace('default_budget = "2m"','default_budget = "5m"'))
    second = load_study(other)
    assert first.study_id == second.study_id
    assert compile_study(first).definition_id == compile_study(second).definition_id
    other.write_text(other.read_text().replace('required_power = 25.0','required_power = 30.0'))
    assert load_study(other).study_id != first.study_id


def test_exact_configuration_keeps_requested_floats_and_integer_types(study_path):
    definition = compile_study(load_study(study_path))
    physical = {p.name:p.initial for p in definition.space.parameters}
    candidate = candidate_for_values(definition,physical)
    assert candidate.payload['physical'] == physical
    assert type(candidate.payload['physical']['microtube.heat_in.tube_count']) is int
    definition.study.space.encode(physical)
    from dada_solver.campaign.candidate import Candidate
    a = [0.5]*5; b = a.copy(); b[1] += 1e-8
    def make(u): return Candidate.create(definition.space,u,families=definition.families,
                                        numerical_settings=definition.numerical_settings,definition_id=definition.definition_id)
    assert make(a).payload['physical'] == make(b).payload['physical']
    assert make(a).candidate_id != make(b).candidate_id
