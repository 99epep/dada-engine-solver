from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tomllib


ROOT = Path(__file__).resolve().parents[1]

DOC = ROOT / "docs" / "MECHANISM_SYNTHESIS_SEARCH.md"
REPRO = ROOT / "docs" / "repro" / "mechanism_synthesis_search"
POLICY = REPRO / "policy.toml"
SATURATION = REPRO / "saturation_reference.toml"
VERIFY = REPRO / "verify.py"

OLD_DOC = ROOT / "docs" / ("MECHANISM_SYNTHESIS_SEARCH_" + "THEORY.md")
OLD_NAME = "MECHANISM_SYNTHESIS_SEARCH_" + "THEORY.md"


def test_document_was_renamed_and_reproduction_exists():
    assert DOC.is_file()
    assert not OLD_DOC.exists()

    assert POLICY.is_file()
    assert SATURATION.is_file()
    assert VERIFY.is_file()

    assert not list(REPRO.glob("*.json"))


def test_documentary_verifier_check_mode():
    completed = subprocess.run(
        [sys.executable, str(VERIFY), "--check"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert completed.stdout.startswith("OK mechanism_synthesis_search")


def test_policy_keeps_reference_search_definition():
    with POLICY.open("rb") as stream:
        policy = tomllib.load(stream)

    assert policy["schema_version"] == 1
    assert policy["reference_motion_id"] == "structured_c2_260k_reference_motor_motion"

    assert policy["primary_geometry"]["branches"] == [-1, 1]
    assert policy["hard_constraints"]["minimum_primary_transmission_sine"] == 0.30
    assert policy["soft_preferences"]["preferred_primary_transmission_sine"] == 0.35

    score = policy["score_weights"]
    assert score["monotonicity"] == 2.0
    assert score["endpoint_zero"] == 0.40
    assert score["turnaround"] == 0.65
    assert score["fast_slow_speed_ratio"] == 0.55
    assert score["fast_displacement_fraction"] == 0.80
    assert score["long_branch_mirror_asymmetry"] == 0.45

    long_branch = policy["long_branch"]
    assert long_branch["criterion"] == "mirror_symmetry"
    assert long_branch["target_profile_matching"] is False
    assert long_branch["smoothness_penalty"] is False
    assert long_branch["curvature_penalty"] is False
    assert long_branch["local_extrema_penalty"] is False

    saturation = policy["fresh_island_saturation"]
    assert saturation["generations"] == 260
    assert saturation["population_size"] == 16
    assert saturation["historical_seeds"] == 0
    assert saturation["expected_initial_population_per_island"] == 96


def test_minimal_saturation_reference():
    with SATURATION.open("rb") as stream:
        reference = tomllib.load(stream)

    counts = reference["family_counts_descending"]

    assert sum(counts) == 512
    assert len(counts) == 16
    assert sum(value == 1 for value in counts) == 8
    assert sum(value == 2 for value in counts) == 3

    assert counts[:4] == [244, 228, 19, 4]

    assert reference["expected_good_turing_unseen_capture_mass"] == 8 / 512
    assert reference["expected_top1_capture_share"] == 244 / 512
    assert reference["expected_top4_capture_share"] == 495 / 512
    assert reference["expected_top10_capture_share"] == 506 / 512


def test_clean_material_has_no_historical_artifact_dependencies():
    files = [DOC, POLICY, SATURATION, VERIFY]

    forbidden = (
        "outputs/",
        "examples/",
        "3952",
        "V1",
        "V2",
        "V3",
        "V4",
        "V4.1",
        "Family 1",
        "Family 4",
        "Family 12",
        "Family 50",
        "rank-1",
        "rank-4",
        "rank-12",
        "rank-50",
    )

    for path in files:
        text = path.read_text(encoding="utf-8")
        for token in forbidden:
            assert token not in text, f"{token!r} remains in {path}"


def test_old_document_name_has_no_remaining_consumers():
    checked_roots = (
        ROOT / "docs",
        ROOT / "src",
        ROOT / "tests",
    )

    suffixes = {".md", ".py", ".toml", ".txt", ".rst"}

    offenders = []
    for base in checked_roots:
        for path in base.rglob("*"):
            if not path.is_file() or path.suffix not in suffixes:
                continue
            if path == Path(__file__).resolve():
                continue
            text = path.read_text(encoding="utf-8")
            if OLD_NAME in text:
                offenders.append(path.relative_to(ROOT))

    assert not offenders, offenders


def test_method_describes_current_handoff_without_campaign_records():
    from dada_solver.research.synthesis import STAGES

    text = DOC.read_text()
    for stage in STAGES:
        assert f'`{stage}`' in text
    assert 'declarative stage and ownership contracts' in text
    for contract in ('GeometrySearch', 'direct_geometry_islands_v1', 'full_local_polish',
                     'category round-robin', 'not a certified continuous root count',
                     'standalone sortable HTML', 'max-evaluations'):
        assert contract in text
    assert 'no injection of design candidates' in text
    assert 'declare the\nweights' in text
    for artifact in ('policy.toml', 'saturation_reference.toml',
                     '## 23. Reproducing', 'completed reference saturation experiment'):
        assert artifact not in text


def test_six_bar_stages_are_real_hierarchical_operators():
    text=DOC.read_text()
    for token in ('hierarchical_six_bar_geometry_v1','SixBarPrimaryMechanism','component = "primary"',
                  'E is not P','deferred','primary_topology_cadence_v2','missing_cadence','PrimaryCadencePolicy',
                  'pivot_envelope_radius','minimum_trajectory_extent','parent lineages never merge',
                  'q_mirror(theta) = q_source(-theta)','--stage opposite_local_adaptation'):
        assert token in text
    assert 'six-bar search remains an explicit implementation boundary' not in text
    assert 'Normalized projection position MSE dominates' not in text
    assert 'primary_linearity_weight' not in text
    assert 'primary_topology_weight' not in text
    assert 'Position and\nvelocity-profile RMS remain diagnostics only' in text


def test_six_bar_design_screen_is_explicit_and_separate_from_fit():
    text=DOC.read_text()
    for token in ('six_bar_design_v1','--mechanical-screen none','1 440 samples',
                  'maximum_position_rms','Each selected parent','Rejections retain aggregate counters'):
        assert token in text


def test_primary_hard_screen_and_soft_symmetry_have_separate_ownership():
    text = DOC.read_text()
    for token in ('six_bar_primary_design_v1', '[primary_mechanical]',
                  'minimum_primary_transmission_sine >= 0.30',
                  'E_span_over_crank >= 0.75', 'max(ptp(E_x), ptp(E_y))',
                  'long_symmetry_weight = 0.45', 'primary_constraints',
                  'primary_transmission_below_minimum', 'primary_e_span_below_minimum'):
        assert token in text
    assert 'Long-branch symmetry is disabled unless requested explicitly' not in text
