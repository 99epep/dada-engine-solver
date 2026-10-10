"""Keep the current user guide focused on usable commands and durable links."""
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "DADA_ENGINE_RESEARCH.md"


def test_current_guide_and_command_workflow():
    text = DOC.read_text(encoding="utf-8")
    assert text.startswith("# Dada Engine Research user guide")
    assert "dada-research" in text
    for command in (
        "init", "validate", "evaluate", "run", "resume", "status",
        "report", "compare", "refine", "rescale",
    ):
        assert re.search(rf"^research {command} ", text, re.MULTILINE)
    for preset in ("kinematics", "external-stream-refrigeration", "external-stream-motor"):
        assert f"research init {preset}" in text


def test_guide_preserves_report_anchor_and_reference_ownership():
    text = DOC.read_text(encoding="utf-8")
    assert "## Reports, curves and animations" in text
    for name in ("REFERENCE", "COCKPIT", "CAPACITY", "KINEMATICS"):
        assert f"(DADA_ENGINE_RESEARCH_{name}.md)" in text


def test_guide_has_no_historical_dependencies_or_legacy_workflow():
    text = DOC.read_text(encoding="utf-8")
    for token in (
        "../outputs/", "examples/", "Historical V1 walkthrough", "sixbar-thermo5d",
        "research init structured-c2-3952",
        "candidate 3952", "candidate 501", "research_kinematics_v2/comparison",
        "DADA_ENGINE_RESEARCH_MIGRATION_MATRIX.md",
    ):
        assert token not in text
