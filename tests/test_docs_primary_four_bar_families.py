from pathlib import Path
import subprocess
import sys
import tomllib


ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "PRIMARY_FOUR_BAR_FAMILIES.md"
REPRO = ROOT / "docs" / "repro" / "primary_four_bar_families"
DATA = REPRO / "families.toml"
VERIFY = REPRO / "verify.py"

EXPECTED_IDS = [
    "compact_balanced",
    "compact_offset",
    "long_ground",
    "long_coupler_short_rocker",
]


def test_primary_four_bar_reproduction_check():
    completed = subprocess.run(
        [sys.executable, str(VERIFY), "--check"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr


def test_primary_four_bar_catalogue_uses_descriptive_ids():
    data = tomllib.loads(DATA.read_text(encoding="utf-8"))
    assert [family["id"] for family in data["family"]] == EXPECTED_IDS
    assert all(not family_id.isdigit() for family_id in EXPECTED_IDS)


def test_primary_four_bar_doc_has_no_historical_output_dependencies():
    text = DOC.read_text(encoding="utf-8")

    forbidden = (
        "outputs/",
        "examples/",
        "3952",
        "Historical ID",
        "historical ID",
        "rank-1",
        "rank-4",
        "rank-12",
        "rank-50",
        "Family 1",
        "Family 4",
        "Family 12",
        "Family 50",
        "V4.1 score",
    )

    for token in forbidden:
        assert token not in text


def test_primary_four_bar_repro_has_no_output_dependency():
    for path in (DATA, VERIFY):
        text = path.read_text(encoding="utf-8")
        assert "outputs/" not in text
        assert "examples/" not in text
        assert "3952" not in text
