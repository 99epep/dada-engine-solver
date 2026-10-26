"""Unit tests for expression objectives, exchanger length-to-diameter ratio, and CLI --set."""
import math
from types import SimpleNamespace
from pathlib import Path
import pytest

from dada_solver.campaign.parameters import ContinuousParameter, IntegerParameter, ChoiceParameter
from dada_solver.exchangers.microtube_geometry import MicrotubeBank
from dada_solver.research.configuration import ConstraintEntry, ExchangerLengthToDiameter
from dada_solver.research.geometry_helpers import exchanger_length_to_diameter
from dada_solver.research.margins import enrich_constraints, PHYSICAL_CONSTRAINTS
from dada_solver.research.objectives import ExpressionObjective
from dada_solver.research.schema import load_study, compile_study, candidate_for_values
from dada_solver.sizing.objectives import ObjectiveValue
from dada_solver.sizing.constraints import ConstraintValue


# =============================================================================
# 1. Tests for exchanger_length_to_diameter helper
# =============================================================================

def test_exchanger_length_to_diameter_direct_attributes():
    geom = SimpleNamespace(microtube_length=0.4, microtube_diameter=0.08)
    assert exchanger_length_to_diameter(geom) == pytest.approx(5.0)


def test_exchanger_length_to_diameter_mapping():
    geom = {"microtube_length": 0.5, "microtube_diameter": 0.1}
    assert exchanger_length_to_diameter(geom) == pytest.approx(5.0)


def test_exchanger_length_to_diameter_microtube_bank():
    # MicrotubeBank with circular collectors computes bundle_diameter_m in dimensions()
    bank = MicrotubeBank(
        tube_count=100,
        tube_length_m=0.3,
        inner_diameter_m=0.002,
        wall_thickness_m=0.0001,
        pitch_ratio=1.25,
        collector_half_angle_deg=30.0,
    )
    bundle_d = bank.dimensions()["bundle_diameter_m"]
    assert bundle_d > 0
    ratio = exchanger_length_to_diameter(bank)
    assert ratio == pytest.approx(0.3 / bundle_d)


def test_exchanger_length_to_diameter_with_exchanger_models():
    bank_in = MicrotubeBank(
        tube_count=100,
        tube_length_m=0.4,
        inner_diameter_m=0.002,
        wall_thickness_m=0.0001,
        pitch_ratio=1.25,
        collector_half_angle_deg=30.0,
    )
    bank_out = MicrotubeBank(
        tube_count=200,
        tube_length_m=0.6,
        inner_diameter_m=0.002,
        wall_thickness_m=0.0001,
        pitch_ratio=1.25,
        collector_half_angle_deg=30.0,
    )
    machine = SimpleNamespace(
        heat_in=SimpleNamespace(bank=bank_in),
        heat_out=SimpleNamespace(bank=bank_out),
    )
    r_in = 0.4 / bank_in.dimensions()["bundle_diameter_m"]
    r_out = 0.6 / bank_out.dimensions()["bundle_diameter_m"]

    assert exchanger_length_to_diameter(machine, side="heat_in") == pytest.approx(r_in)
    assert exchanger_length_to_diameter(machine, side="heat_out") == pytest.approx(r_out)
    # Default without side returns min(r_in, r_out)
    assert exchanger_length_to_diameter(machine) == pytest.approx(min(r_in, r_out))


def test_exchanger_length_to_diameter_invalid_values():
    with pytest.raises(AttributeError):
        exchanger_length_to_diameter(None)

    with pytest.raises(AttributeError):
        exchanger_length_to_diameter(SimpleNamespace())

    with pytest.raises(ValueError, match="positive"):
        exchanger_length_to_diameter(SimpleNamespace(microtube_length=0.4, microtube_diameter=0.0))

    with pytest.raises(ValueError, match="positive"):
        exchanger_length_to_diameter(SimpleNamespace(microtube_length=0.4, microtube_diameter=-0.05))

    with pytest.raises(ValueError, match="nonnegative"):
        exchanger_length_to_diameter(SimpleNamespace(microtube_length=-0.1, microtube_diameter=0.05))


# =============================================================================
# 2. Tests for ExchangerLengthToDiameter constraint
# =============================================================================

def test_exchanger_length_to_diameter_constraint_class():
    c = ExchangerLengthToDiameter(value=3.0)
    assert c.name == "exchanger_length_to_diameter"
    assert c.type == "exchanger_length_to_diameter"
    assert c.value == 3.0

    with pytest.raises(ValueError, match="greater than 0"):
        ExchangerLengthToDiameter(value=0.0)

    with pytest.raises(ValueError, match="greater than 0"):
        ExchangerLengthToDiameter(value=-1.5)


def test_exchanger_length_to_diameter_constraint_evaluation():
    constraint = ExchangerLengthToDiameter(value=4.0)

    # Satisfied evaluation: ratio = 5.0 >= 4.0
    evaluation_sat = SimpleNamespace(
        design=SimpleNamespace(microtube_length=0.5, microtube_diameter=0.1)
    )
    res_sat = constraint.evaluate(evaluation_sat)
    assert isinstance(res_sat, ConstraintValue)
    assert res_sat.name == "exchanger_length_to_diameter"
    assert res_sat.available is True
    assert res_sat.satisfied is True
    assert res_sat.margin == pytest.approx(1.0)

    # Violated evaluation: ratio = 3.0 < 4.0
    evaluation_viol = SimpleNamespace(
        design=SimpleNamespace(microtube_length=0.3, microtube_diameter=0.1)
    )
    res_viol = constraint.evaluate(evaluation_viol)
    assert res_viol.available is True
    assert res_viol.satisfied is False
    assert res_viol.margin == pytest.approx(-1.0)

    # Unavailable evaluation (missing geometry)
    res_unavail = constraint.evaluate(SimpleNamespace(design=SimpleNamespace()))
    assert res_unavail.available is False
    assert res_unavail.satisfied is False
    assert res_unavail.margin is None


def test_enrich_constraints_with_exchanger_length_to_diameter():
    declarations = [{"type": "exchanger_length_to_diameter", "value": 3.0, "unit": "1"}]
    records = [{"name": "exchanger_length_to_diameter", "margin": 1.5, "satisfied": True, "available": True}]
    enriched = enrich_constraints(records, declarations)
    assert len(enriched) == 1
    item = enriched[0]
    assert item["name"] == "exchanger_length_to_diameter"
    assert item["limit"] == 3.0
    assert item["value"] == pytest.approx(4.5)  # limit + margin for minimum relation
    assert item["satisfied"] is True
    assert item["state"] == "satisfied"


# =============================================================================
# 3. Tests for ExpressionObjective
# =============================================================================

def test_expression_objective_standalone_evaluation():
    obj = ExpressionObjective(name="test_obj", expr="power * 2 + 10")
    val = obj.evaluate({"power": 50.0})
    assert val == pytest.approx(110.0)


def test_expression_objective_caret_exponentiation():
    obj = ExpressionObjective(name="test_obj", expr="(ratio - 1)^2")
    val = obj.evaluate({"ratio": 4.0})
    assert val == pytest.approx(9.0)


def test_expression_objective_math_functions():
    obj = ExpressionObjective(name="test_obj", expr="sqrt(power) * abs(ratio)")
    val = obj.evaluate({"power": 100.0, "ratio": -2.5})
    assert val == pytest.approx(25.0)


def test_expression_objective_ratio_from_geometry():
    obj = ExpressionObjective(name="test_obj", expr="cooling_power * ratio")
    geom = SimpleNamespace(microtube_length=0.6, microtube_diameter=0.2)
    val = obj.evaluate({"cooling_power": 150.0}, geometry=geom)
    assert val == pytest.approx(450.0)


def test_expression_objective_standalone_errors():
    obj = ExpressionObjective(name="test_obj", expr="unknown_var + 1")
    with pytest.raises(ValueError, match="Error evaluating expression"):
        obj.evaluate({"power": 50.0})

    obj_div0 = ExpressionObjective(name="test_obj", expr="1 / 0")
    with pytest.raises(ValueError, match="Error evaluating expression"):
        obj_div0.evaluate({})


def test_expression_objective_campaign_evaluation():
    perf = SimpleNamespace(
        gas_power=120.0,
        heat_in_power=200.0,
        heat_out_power=-320.0,
        thermal_efficiency=0.375,
        mechanical_input_power=120.0,
        motor_power=120.0,
        cooling_power=None,
        cooling_cop=None,
    )
    geom = SimpleNamespace(microtube_length=0.4, microtube_diameter=0.1)
    evaluation = SimpleNamespace(usable=True, performance=perf, design=geom)

    # Maximize direction: returns negated score
    obj_max = ExpressionObjective(name="obj_max", expr="power * ratio", direction="maximize")
    res_max = obj_max.evaluate(evaluation)
    assert isinstance(res_max, ObjectiveValue)
    assert res_max.available is True
    # 120 * 4.0 = 480 -> negated for maximization: -480
    assert res_max.value == pytest.approx(-480.0)

    # Minimize direction: returns positive score
    obj_min = ExpressionObjective(name="obj_min", expr="power * ratio", direction="minimize")
    res_min = obj_min.evaluate(evaluation)
    assert res_min.available is True
    assert res_min.value == pytest.approx(480.0)


def test_expression_objective_unusable_evaluation():
    obj = ExpressionObjective(name="test_obj", expr="power")
    res_unusable = obj.evaluate(SimpleNamespace(usable=False, performance=None))
    assert isinstance(res_unusable, ObjectiveValue)
    assert res_unusable.available is False
    assert res_unusable.value is None


# =============================================================================
# 4. Tests for Study TOML loading & compilation with new features
# =============================================================================

def test_load_and_compile_study_with_expression_objective_and_exchanger_constraint(tmp_path):
    slider_toml = Path("examples/pedal_cooling_cell/slider/slider_coolcell.toml").read_text()
    p_expr = slider_toml.replace(
        '"type" = "maximize_cooling_cop"',
        '"type" = "cooling_power_w * ratio"\n"direction" = "maximize"',
    )
    extra = """
[[ "constraints" ]]
"type" = "exchanger_length_to_diameter"
"value" = 2.5
"unit" = "1"
"""
    custom_study_path = tmp_path / "study.toml"
    # Copy required artifact and basis files
    src_dir = Path("examples/pedal_cooling_cell/slider")
    for file in src_dir.glob("*.json"):
        (tmp_path / file.name).write_bytes(file.read_bytes())
    custom_study_path.write_text(p_expr + extra)

    study = load_study(custom_study_path)
    assert study.data["objective"]["type"] == "cooling_power_w * ratio"
    assert study.data["objective"]["direction"] == "maximize"

    defn = compile_study(study)
    assert isinstance(defn.objective, ExpressionObjective)
    assert defn.objective.expr == "cooling_power_w * ratio"
    assert defn.objective.direction == "maximize"

    constraint_names = [c.name for c in defn.constraints]
    assert "exchanger_length_to_diameter" in constraint_names
    ex_c = next(c for c in defn.constraints if c.name == "exchanger_length_to_diameter")
    assert isinstance(ex_c, ExchangerLengthToDiameter)
    assert ex_c.value == 2.5


def test_invalid_expression_syntax_rejected_at_load_time(tmp_path):
    slider_toml = Path("examples/pedal_cooling_cell/slider/slider_coolcell.toml").read_text()
    p_invalid = slider_toml.replace(
        '"type" = "maximize_cooling_cop"',
        '"type" = "power + * 2"',
    )
    custom_study_path = tmp_path / "invalid_study.toml"
    src_dir = Path("examples/pedal_cooling_cell/slider")
    for file in src_dir.glob("*.json"):
        (tmp_path / file.name).write_bytes(file.read_bytes())
    custom_study_path.write_text(p_invalid)

    with pytest.raises(Exception, match="Invalid expression objective syntax"):
        load_study(custom_study_path)


# =============================================================================
# 5. Tests for CLI evaluate --set option
# =============================================================================

def test_cli_evaluate_with_set_assignment(tmp_path):
    from dada_solver.research.cli import evaluate

    study_path = Path("examples/pedal_cooling_cell/slider/slider_coolcell.toml")
    output_path = tmp_path / "eval_out.json"

    # Override an active parameter with --set
    record = evaluate(
        study_path,
        output_path,
        assignments=["microtube.heat_in.tube_count=220"],
    )
    assert record["status"] == "feasible"
    assert record["physical"]["microtube.heat_in.tube_count"] == 220
    assert output_path.exists()
