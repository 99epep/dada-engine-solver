"""Narrative Research archives stay separate from current references."""
from pathlib import Path
import re

DOCS = Path(__file__).resolve().parents[1] / 'docs'
NAMES = tuple(f'DADA_ENGINE_RESEARCH_{part}.md' for part in (
    'ARCHITECTURE_AUDIT', 'IMPLEMENTATION_PLAN', 'MIGRATION_MATRIX',
))


def test_archives_live_in_history_and_are_indexed():
    index = (DOCS / 'history/README.md').read_text()
    for name in NAMES:
        assert not (DOCS / name).exists()
        assert (DOCS / 'history' / name).is_file()
        assert f']({name})' in index
    assert '](RESEARCH_VALIDATION_LEDGER.md)' in index


def test_root_document_links_use_history_paths():
    for path in DOCS.glob('*.md'):
        targets = re.findall(r'\]\(([^)]+)\)', path.read_text())
        assert not set(NAMES).intersection(targets), path


def test_audit_material_has_not_moved():
    assert (DOCS / 'DADA_ENGINE_RESEARCH_MIGRATION_MATRIX.csv').is_file()
    assert (DOCS / 'research_audit/generate_inventory.py').is_file()
    matrix = (DOCS / 'history' / NAMES[-1]).read_text()
    assert '](../DADA_ENGINE_RESEARCH_MIGRATION_MATRIX.csv)' in matrix
    assert '](../research_audit/example_inventory.json)' in matrix
