"""Retain current exchanger evidence without obsolete trial redirects."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REMOVED = (
    'DOUBLED_EXCHANGER_TRIAL.md', 'EXCESS_AIR_TRIAL.md',
    'HIGHER_TEMPERATURE_TRIAL.md', 'PARALLEL_EXCHANGER_TRIAL.md',
    'history/EXCHANGER_TRIALS.md', 'history/DOCUMENTATION_MAP.md',
)


def test_obsolete_trial_layer_is_absent():
    for name in REMOVED:
        assert not (ROOT / 'docs' / name).exists()
    instructions = (ROOT / 'AGENTS.md').read_text()
    for name in REMOVED[:4]:
        assert name not in instructions


def test_scientific_references_remain():
    for name in ('DOTY_SCREENING', 'EXCHANGER_VALIDATION',
                 'MICROTUBE_GAS_MODEL', 'EXTERNAL_STREAM_THERMAL_MODEL', 'HEAT_EXCHANGERS'):
        assert (ROOT / 'docs' / (name + '.md')).is_file()


def test_exchanger_validation_is_independent_of_application_artifacts():
    text = (ROOT / 'docs' / 'EXCHANGER_VALIDATION.md').read_text()
    for stale_reference in ('examples/', 'outputs/', 'recorded September 2026'):
        assert stale_reference not in text


def test_sizing_references_are_independent_of_application_artifacts():
    for name in ('SIZING.md', 'HEAT_EXCHANGERS.md'):
        text = (ROOT / 'docs' / name).read_text()
        for stale_reference in ('examples/', 'outputs/'):
            assert stale_reference not in text
