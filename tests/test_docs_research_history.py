"""Consumed Research migration material stays removed."""
from pathlib import Path

DOCS = Path(__file__).resolve().parents[1] / 'docs'
RESEARCH_ARCHIVES = tuple(f'DADA_ENGINE_RESEARCH_{part}.md' for part in (
    'ARCHITECTURE_AUDIT', 'IMPLEMENTATION_PLAN', 'MIGRATION_MATRIX',
))
NAMES = RESEARCH_ARCHIVES + ('RESEARCH_LIMIT_OWNERSHIP_AUDIT.md',)


def test_consumed_research_archives_are_removed():
    names = NAMES + ('RESEARCH_VALIDATION_LEDGER.md',)
    for name in names:
        assert not (DOCS / name).exists()
        assert not (DOCS / 'history' / name).exists()
    for path in DOCS.rglob('*.md'):
        text = path.read_text()
        assert not any(name in text for name in names), path


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
