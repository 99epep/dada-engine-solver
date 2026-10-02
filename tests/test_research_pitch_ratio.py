"""Explicit ratio ownership changes packing without changing legacy studies."""
import tomllib
import pytest
from dada_solver.research.presets import initialize_v2
from dada_solver.research.schema import load_study, compile_study
from dada_solver.research.study_io import dumps


def test_ratio_parameter_with_legacy_basis(tmp_path):
    path = initialize_v2(tmp_path/'study.toml', 'harmonic', 'harmonic')
    raw = tomllib.loads(path.read_text())
    raw['parameters'] = [r for r in raw['parameters'] if not r['name'].endswith('.pitch_m')]
    for side in ('heat_in','heat_out'):
        raw['parameters'].append(dict(name=f'microtube.{side}.pitch_ratio', unit='1',
            kind='continuous', initial=1.2, lower=1.01, upper=2., transform='linear'))
    path.write_text(dumps(raw))
    study = load_study(path)
    definition = compile_study(study)
    physical = dict(definition.fixed_parameters, **{p.name:p.initial for p in definition.space.parameters})
    for diameter in (.00008, .0018):
        physical['microtube.heat_in.inner_diameter_m'] = diameter
        bank = definition.adapter.build(physical).heat_in.bank
        assert bank.pitch_ratio == 1.2
        assert bank.effective_pitch_m == 1.2*(diameter+2*bank.wall_thickness_m)
    raw['parameters'][-1]['lower'] = 1.
    path.write_text(dumps(raw))
    with pytest.raises(ValueError, match='pitch_ratio'): load_study(path)
