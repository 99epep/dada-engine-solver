from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CURRENT = ROOT / "docs/validation.md"


def test_current_validation_remains_authoritative():
    text = CURRENT.read_text()
    assert text.startswith("# Validation and evidence")
    assert "## Current regression evidence" in text
    assert "## Physical limits of the evidence" in text
    assert "not validate" in text.lower()


def test_parity_boundaries_remain_in_current_reference():
    text = CURRENT.read_text()
    assert "## Software parity and regression tolerances" in text
    for concept in ("adaptive integration", "physical acceptance", "old model",
                    "current-physics regressions"):
        assert concept in text
    assert (ROOT / "tests/test_research_v2_thermodynamics.py").is_file()


def test_old_validation_page_is_not_reintroduced():
    name = "DADA_ENGINE_RESEARCH_VALIDATION.md"
    assert not (ROOT / "docs" / name).exists()
    for path in (ROOT / "docs").rglob("*.md"):
        assert name not in path.read_text(), path
