from pathlib import Path

from dada_solver.research.study_schema import POLICIES_V3
from dada_solver.research.charge import POLICY as CHARGE_POLICY, METHOD as CHARGE_METHOD
from dada_solver.research.margins import PHYSICAL_CONSTRAINTS
from dada_solver.campaign.parameters import ContinuousParameter, IntegerParameter, ChoiceParameter


ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "DADA_ENGINE_RESEARCH_REFERENCE.md"
GUIDE = ROOT / "docs" / "DADA_ENGINE_RESEARCH.md"


def text(path=DOC):
    return path.read_text(encoding="utf-8")


def test_reference_has_current_role():
    content = text()
    assert content.startswith("# Research technical reference")
    assert "Research reference and historical notes" not in content
    assert "## Historical V1 walkthrough" not in content


def test_supported_schema_generations_are_documented():
    content = text()
    assert "`schema_version = 3`" in content
    assert "performs no automatic schema conversion" in content


def test_current_policy_values_are_documented():
    content = text()
    for value in set(POLICIES_V3.values()):
        assert f"`{value}`" in content
    assert f"`{CHARGE_POLICY}`" in content
    assert f"`{CHARGE_METHOD}`" in content


def test_current_constraint_vocabulary_is_documented():
    content = text()
    for name, (field, relation, unit) in PHYSICAL_CONSTRAINTS.items():
        assert f"`{name}`" in content
        row = f"| `{name}` | {'None' if field is None else f'`{field}`'} | `{unit}` |"
        assert row in content


def test_parameter_kinds_and_identity_levels_are_documented():
    content = text()
    assert ContinuousParameter
    assert IntegerParameter
    assert ChoiceParameter
    for token in ("continuous", "integer", "choice", "nearest_even_v1",
                  "study_id", "definition_id", "candidate_id"):
        assert f"`{token}`" in content or token in content


def test_search_domains_and_scheduling_are_documented():
    content = text()
    for token in ("fixed_global_bounds", "local_regions_v1", "evaluate_initial",
                  "round_robin", "evaluate_centers", "radius_fraction"):
        assert f"`{token}`" in content or token in content


def test_required_external_anchors_remain():
    content = text()
    assert "## Local refinement and explicit initial evaluations" in content
    assert "## Filling at a reference pressure and maximum total gas volume" in content


def test_reference_has_no_historical_output_dependency():
    content = text()
    forbidden = ("../outputs/", "examples/", "## Historical V1 walkthrough",
                 "bounded Human Cell check", "Scope of this reference machine",
                 "Acceptance and the next review", "298.15/558.15")
    for token in forbidden:
        assert token not in content


def test_user_guide_describes_reference_as_current():
    content = text(GUIDE)
    assert "Reference and historical details" not in content
    assert "historical presets" not in content


def test_microtube_productivity_objective_contract():
    content = text()
    assert '| `maximize_cooling_power_per_total_microtube` | `W/microtube` | Refrigeration |' in content
    assert 'cooling_power_per_total_microtube_w' in content
    assert 'microtube.heat_in.tube_count + microtube.heat_out.tube_count' in content
    assert '| `maximize_cooling_cop_times_power_per_total_microtube` | `W/microtube` | Refrigeration |' in content
    assert ('`maximize_cooling_cop_times_power_per_total_microtube` maximizes\n'
            '`cooling_cop * cooling_power_per_total_microtube_w`.') in content
    assert 'does not impose minimum COP or cooling power by itself' in content


def test_study_editor_documents_domains_and_scheduler_separately():
    content = text().split("## Study editor", 1)[1].split("\n## ", 1)[0]
    for command in ("study edit", "study release", "study freeze"):
        assert command in content
    for contract in ("prompt_toolkit", "reference_boxes", "--bounds", "--choices",
                     "fixed_global_bounds", "local_regions_v1", "--recenter"):
        assert contract in content
    assert "all active parameters" in content
    assert "unevaluated" in content
    assert "artifact" in content
    assert "integration" in content
