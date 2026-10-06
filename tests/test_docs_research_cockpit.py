"""Durable cockpit documentation boundaries, without frozen report snapshots."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / 'docs/DADA_ENGINE_RESEARCH_COCKPIT.md'


def test_current_cockpit_role_and_links():
    text = DOC.read_text()
    assert text.startswith('# Research cockpit and execution UX')
    for target in ('DADA_ENGINE_RESEARCH.md', 'DADA_ENGINE_RESEARCH_REFERENCE.md'):
        assert f']({target})' in text
    for token in ('TTY', 'Non-TTY', 'derived', '`bounds`', '`filter`'):
        assert token in text
    assert 'never executes' in text


def test_journal_documentation_matches_current_path_selection(tmp_path):
    from dada_solver.campaign.history import journal_path
    import pytest

    text = DOC.read_text()
    assert journal_path(tmp_path).name == 'history.jsonl.gz'
    (tmp_path / 'history.jsonl').touch()
    with pytest.raises(ValueError):
        journal_path(tmp_path)
    (tmp_path / 'history.jsonl.gz').touch()
    with pytest.raises(ValueError, match='Uncompressed'):
        journal_path(tmp_path)
    for token in ('history.jsonl.gz', 'recovery.json', 'rejected',
                  'CRC', 'Read-only inspection', 'Execution resume'):
        assert token in text
    assert 'Append the complete result to `history.jsonl`' not in text


def test_no_historical_artifact_dependencies():
    text = DOC.read_text()
    for token in ('../outputs/', 'examples/', 'human_cell_stage2c_multi'):
        assert token not in text
