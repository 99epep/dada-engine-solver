from pathlib import Path

from dada_solver.research.families import (
    FAMILIES,
    PHYSICAL_FAMILIES,
    SIXBAR_CONTINUOUS,
    PRIMARY_COORDINATES,
    DOWNSTREAM_COORDINATES,
    available_metrics,
    parameter_specs,
)
from dada_solver.research.synthesis import STAGES, release_coordinates


ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "DADA_ENGINE_RESEARCH_KINEMATICS.md"


def text():
    return DOC.read_text(encoding="utf-8")


def test_document_is_current_unversioned_reference():
    content = text()
    assert content.startswith("# Research kinematics and mechanism reference")
    assert "Research V2 — kinematics" not in content
    assert "V2 reference" not in content


def test_current_family_registry_is_documented():
    content = text()
    assert len(FAMILIES) == 11
    assert set(PHYSICAL_FAMILIES) == {"slider_crank", "four_bar", "six_bar"}
    for family in FAMILIES:
        assert f"`{family}`" in content


def test_six_bar_coordinate_contract_is_documented():
    content = text()
    assert len(SIXBAR_CONTINUOUS) == 15
    assert len(PRIMARY_COORDINATES) == 6
    assert len(DOWNSTREAM_COORDINATES) == 9
    for name in SIXBAR_CONTINUOUS:
        assert f"`{name}`" in content


def test_synthesis_release_contract_is_current():
    assert len(release_coordinates("primary_discovery", ("large",))) == 6
    assert len(release_coordinates("downstream_fit", ("large",))) == 9
    assert len(release_coordinates("full_local_polish", ("large",))) == 15
    paired = release_coordinates("paired_thermodynamic", ("small", "large"))
    assert len(paired) == 30
    assert not any("branch" in name for name in paired)
    assert release_coordinates("mirror_initialization") == ()
    assert release_coordinates("hardware_retuning") == ()
    for stage in STAGES:
        assert f"`{stage}`" in text()


def test_documented_physical_metrics_are_current():
    content = text()
    for family in PHYSICAL_FAMILIES:
        for metric in available_metrics(family):
            if metric in {"maximum_absolute_first_derivative", "maximum_absolute_second_derivative"}:
                continue
            assert f"`{metric}`" in content


def test_current_public_interfaces_are_documented():
    content = text()
    for symbol in (
        "MechanismArtifact", "MechanismLibrary", "side_metrics", "available_metrics",
        "SynthesisRequest", "FreshIslandPolicy", "release_coordinates", "sample_motion",
        "sample_report_volumes", "parameter_specs", "margin_record", "validate_mechanical_constraint",
    ):
        assert f"`{symbol}" in content


def test_historical_campaign_material_is_absent():
    content = text()
    forbidden = (
        "candidate 501", "candidate 3952", "candidate 3335", "ranks 1, 4, 12 and 50",
        "DD5", "../outputs/", "examples/", "tests/fixtures/hybrid_compact",
        "tests/fixtures/research_v2", "research_kinematics_v2/comparison.html",
        "pitch_ratio", "collector_half_angle_deg", "conduit_area_ratio",
        "additional_internal_volume_m3", "geometry_conduit_area_v1",
    )
    for token in forbidden:
        assert token not in content


def test_named_family_coordinates_follow_production_ownership():
    content = text()
    for family in FAMILIES:
        if family in ("free_spline", "fourier_c2"):
            continue  # Variable-size coordinates are described by ranges.
        settings = dict(family=family)
        if family == "four_bar":
            settings["output"] = "rocker"
        for side in ("small", "large"):
            for name in parameter_specs(settings, side):
                assert f"`{name}`" in content
    assert len(parameter_specs(dict(family="structured_c2_15p"), "small")) == 8
    assert len(parameter_specs(dict(family="structured_c2_15p"), "large")) == 7
    assert len(parameter_specs(dict(family="hybrid_compact"), "small")) == 5
    assert len(parameter_specs(dict(family="hybrid_compact"), "large")) == 4


def test_second_derivative_availability_matches_reference():
    assert {family for family in FAMILIES
            if "maximum_absolute_second_derivative" in available_metrics(family)} == {
        "harmonic", "slider_crank", "free_spline", "fourier_c2", "structured_c2_15p",
    }


def test_motion_target_family_protocols_and_interactive_boundary_are_documented():
    content = text()
    for symbol in ('MotionTarget', 'SynthesisPlan', 'MotionRefitRequest',
                   'synthesis_protocol', 'mechanism_state', 'animate_mechanism',
                   'plot_motion_comparison', 'mechanism_catalogue'):
        assert symbol in content
    assert '15 structured coordinates on the pair' in content
    assert 'missing' in content.lower() or 'absent derivatives' in content
    assert 'never finite-differenced' in content
    assert 'not implemented' in content
    assert 'FuncAnimation' in content
    assert 'global_discovery' in content


def test_paired_thermodynamic_and_hardware_generation_use_real_research():
    content=text()
    for token in ('mechanism adapt', 'paired_thermodynamic()', 'local_regions_v1',
                  '6 for slider-crank, 22 for four-bar, 30 for six-bar',
                  'source_exact', 'hardware_retuning()', 'mechanism retune',
                  'source-active', 'zero mechanism coordinates', 'choice_scope'):
        assert token in content
    assert 'paired thermodynamics and hardware retuning are\nnot implemented' not in content


def test_refit_is_current_geometric_initialization_and_study_handoff():
    content = text()
    for token in ('MotionRefitResult', 'structured_c2_15p', 'exactly 15 active parameters',
                  't = -theta/(2*pi)', 'reference-pressure', '--report', '--radius'):
        assert token in content
    assert 'No thermodynamic integration' in content


def test_refit_search_defaults_are_structural_local_bounds():
    content = text()
    assert '0.1 and lies in (0, 0.5]' in content
    assert 'half to 1.5 times the fitted value' in content
    assert 'not a physical displacement metric' in content
    assert 'Existing studies and campaign identities are not rewritten' in content


def test_executable_six_bar_hierarchy_and_independent_pair_are_documented():
    content=text()
    for token in ('SixBarPrimaryMechanism','component = "primary"','primary_discovery',
                  'downstream_fit','opposite_local_adaptation','first side untouched',
                  'no implicit symmetry coupling','No global fifteen-dimensional fit'):
        assert token in content
    assert 'Six-bar hierarchy, fresh-primary' not in content
