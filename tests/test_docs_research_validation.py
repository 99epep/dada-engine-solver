from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CURRENT = ROOT / "docs" / "validation.md"
HISTORICAL = ROOT / "docs" / "history" / "RESEARCH_VALIDATION_LEDGER.md"
OLD = ROOT / "docs" / "DADA_ENGINE_RESEARCH_VALIDATION.md"


def read(path):
    return path.read_text(encoding="utf-8")


def test_validation_roles_are_separated():
    assert CURRENT.is_file()
    assert HISTORICAL.is_file()
    assert not OLD.exists()
    current = read(CURRENT)
    historical = read(HISTORICAL)
    assert current.startswith("# Validation and evidence")
    assert historical.startswith("# Research implementation validation ledger")
    assert "current validation status" in historical.lower()
    assert "../validation.md" in historical


def test_historical_ledger_has_no_runtime_output_dependencies():
    historical = read(HISTORICAL)
    forbidden = (
        "../outputs/", "outputs/", "examples/", "candidate 3952", "candidate 501",
        "ranks 1/4/12/50", "595 passed", "679 passed", "744 passed", "760 passed",
        "782 passed", "912 passed",
    )
    for token in forbidden:
        assert token not in historical


def test_historical_ledger_keeps_validation_boundaries():
    historical = read(HISTORICAL).lower()
    for concept in (
        "software parity", "physical validation", "scientific identity", "periodic",
        "mechanism artifact", "external", "capacity", "journal", "inventory", "limit ownership",
    ):
        assert concept in historical


def test_current_validation_remains_authoritative():
    current = read(CURRENT)
    assert "## Current check record" in current
    assert "## Physical limits of the evidence" in current
    assert "not validate" in current.lower()


def test_no_docs_reference_the_old_filename():
    old_name = "DADA_ENGINE_RESEARCH_VALIDATION.md"
    offenders = []
    for path in (ROOT / "docs").rglob("*.md"):
        if old_name in read(path):
            offenders.append(path.relative_to(ROOT).as_posix())
    assert not offenders, offenders
