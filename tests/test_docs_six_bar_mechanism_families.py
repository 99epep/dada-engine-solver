from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tomllib

import pytest

from dada_solver.six_bar import SixBarCylinderMechanism
from dada_solver.mechanism_diagnostics import six_bar_metrics


ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "SIX_BAR_MECHANISM_FAMILIES.md"
REPRO = ROOT / "docs" / "repro" / "six_bar_mechanism_families"
FAMILIES = REPRO / "families.toml"
SCREEN = REPRO / "design_screen.toml"
VERIFY = REPRO / "verify.py"

EXPECTED_IDS = (
    "compact_balanced",
    "compact_offset",
    "long_ground",
    "long_coupler_short_rocker",
)


def load_families():
    with FAMILIES.open("rb") as stream:
        return tomllib.load(stream)


def load_verifier():
    spec = importlib.util.spec_from_file_location("six_bar_catalogue_verifier", VERIFY)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def catalogue():
    verifier = load_verifier()
    screen = verifier.load_toml(SCREEN)
    return verifier, verifier.evaluate_catalogue(load_families(), screen)


def test_reproduction_files_exist_without_json():
    assert DOC.is_file()
    assert FAMILIES.is_file()
    assert SCREEN.is_file()
    assert VERIFY.is_file()
    assert not list(REPRO.glob("*.json"))


def test_stable_family_ids_and_eight_mechanisms():
    data = load_families()
    assert tuple(data["family_order"]) == EXPECTED_IDS
    assert set(data["families"]) == set(EXPECTED_IDS)
    count = 0
    for family_id in EXPECTED_IDS:
        family = data["families"][family_id]
        assert family["primary_seed_family"] == family_id
        assert "small" in family
        assert "large" in family
        count += 2
    assert count == 8


def test_exact_geometry_sentinels():
    data = load_families()["families"]
    assert data["compact_balanced"]["small"]["primary_ground"] == 1.6171985844849184
    assert data["compact_offset"]["large"]["link_ef"] == 1.5033623712897515
    assert data["long_ground"]["small"]["piston_rod"] == 8.671128855660843
    assert data["long_coupler_short_rocker"]["large"]["link_gf"] == 9.769509440556504


def test_all_eight_construct_with_production_mechanism():
    data = load_families()["families"]
    for family_id in EXPECTED_IDS:
        for side in ("small", "large"):
            mechanism = SixBarCylinderMechanism(**data[family_id][side])
            metrics = six_bar_metrics(mechanism, samples=1440)
            assert metrics["zero_crossing_count"] == 2
            assert metrics["stroke_over_crank"] > 0.0


def test_documentary_verifier_check_mode():
    completed = subprocess.run(
        [sys.executable, str(VERIFY), "--check"],
        cwd=ROOT,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert completed.stdout.startswith("OK six_bar_mechanism_families")
    assert "families=4" in completed.stdout
    assert "mechanisms=8" in completed.stdout
    assert "samples=1440/5760" in completed.stdout


def test_clean_material_has_no_historical_dependencies():
    forbidden = (
        "outputs/", "examples/", "3952", "candidate_id",
        "rank_01", "rank_04", "rank_12", "rank_50",
        "Family 1", "Family 4", "Family 12", "Family 50", "thermo5d",
        "2048", "tolerance_reference.toml",
    )
    for path in (DOC, FAMILIES, SCREEN, VERIFY):
        text = path.read_text(encoding="utf-8")
        for token in forbidden:
            assert token not in text, f"{token!r} remains in {path}"


def test_screen_preserves_all_standard_limits():
    with SCREEN.open("rb") as stream:
        screen = tomllib.load(stream)
    assert screen["canonical_samples"] == 1440
    assert screen["dense_cross_check_samples"] == 5760
    assert [(r["metric"], r["relation"], r["limit"], r["unit"])
            for r in screen["constraints"]] == [
        ("stroke_over_crank", "minimum", 1.0, "crank_radius"),
        ("stroke_over_crank", "maximum", 3.0, "crank_radius"),
        ("minimum_primary_transmission_sine", "minimum", 0.30, "1"),
        ("minimum_secondary_transmission_sine", "minimum", 0.30, "1"),
        ("minimum_rod_axis_cosine", "minimum", 0.95, "1"),
        ("EH_over_crank", "maximum", 7.0, "crank_radius"),
        ("crank_axis_to_EFH_clearance_over_crank", "minimum", 0.5, "crank_radius"),
        ("H_axis_lateral_rms_over_stroke", "maximum", 0.25, "1"),
        ("H_axis_lateral_span_over_stroke", "maximum", 0.65, "1"),
        ("zero_crossing_count", "equal", 2, "1"),
    ]


def test_both_screens_and_every_individual_margin(catalogue):
    verifier, results = catalogue
    assert verifier.screen_failures(results) == []
    for sides in results.values():
        for record in sides.values():
            assert len(record["constraints"]) == 10
            for row in record["constraints"]:
                assert row["canonical_pass"] and row["dense_pass"]
                if row["relation"] == "equal":
                    assert row["margin"] is None
                    assert row["dense_margin"] is None
                else:
                    for value, margin in (("canonical", "margin"), ("dense", "dense_margin")):
                        expected = (row[value] - row["limit"] if row["relation"] == "minimum"
                                    else row["limit"] - row[value])
                        assert row[margin] == expected
                        assert row[margin] >= 0


def test_markdown_tables_reproduce_document(catalogue):
    verifier, results = catalogue
    doc = DOC.read_text(encoding="utf-8")
    for section in verifier.markdown_tables(results).split("### `")[1:]:
        _, tables = section.split("`\n", 1)
        assert tables.strip() in doc
    assert "INSERT_" not in doc


def test_dense_only_failure_is_not_hidden(monkeypatch, capsys):
    verifier = load_verifier()
    row = dict(metric="H_axis_lateral_span_over_stroke", relation="maximum", limit=0.65,
               unit="1", canonical=0.64999, dense=0.65001,
               margin=0.00001, dense_margin=-0.00001,
               canonical_pass=True, dense_pass=False)
    results = {"compact_balanced": {"large": {"constraints": [row]}}}
    monkeypatch.setattr(verifier, "evaluate_catalogue", lambda *args: results)
    monkeypatch.setattr(sys, "argv", [str(VERIFY), "--check"])
    with pytest.raises(SystemExit) as error:
        verifier.main()
    message = str(error.value)
    assert "compact_balanced/large" in message
    assert "canonical=0.64999 pass=True" in message
    assert "dense=0.65001 pass=False" in message
    assert "OK" not in capsys.readouterr().out
