"""Narrative Research archives stay separate from current references."""
from pathlib import Path
import re

DOCS = Path(__file__).resolve().parents[1] / 'docs'
RESEARCH_ARCHIVES = tuple(f'DADA_ENGINE_RESEARCH_{part}.md' for part in (
    'ARCHITECTURE_AUDIT', 'IMPLEMENTATION_PLAN', 'MIGRATION_MATRIX',
))
NAMES = RESEARCH_ARCHIVES + ('RESEARCH_LIMIT_OWNERSHIP_AUDIT.md',)


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


def test_intermediate_audit_material_is_removed():
    assert not (DOCS / 'research_audit').exists()
    assert not (DOCS / 'DADA_ENGINE_RESEARCH_MIGRATION_MATRIX.csv').exists()


def test_current_ownership_authority():
    index = (DOCS / 'README.md').read_text()
    row = next(line for line in index.splitlines()
               if 'Model domain versus design requirement' in line)
    assert '](PHYSICS_DECISIONS.md)' in row
    for name in ('DADA_ENGINE_RESEARCH_REFERENCE.md',
                 'DADA_ENGINE_RESEARCH_KINEMATICS.md'):
        text = (DOCS / name).read_text()
        assert '[current ownership rule](PHYSICS_DECISIONS.md)' in text
        assert '](RESEARCH_LIMIT_OWNERSHIP_AUDIT.md)' not in text
    audit = (DOCS / 'history/RESEARCH_LIMIT_OWNERSHIP_AUDIT.md').read_text()
    assert '](../PHYSICS_DECISIONS.md)' in audit
