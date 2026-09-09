from pathlib import Path

import pytest

from dada_solver.exchangers.doty import load_references, scale_bank

DATA = Path(__file__).resolve().parents[1] / "examples/data/doty_1991_nitrogen_reference.csv"


def estimate(reference, banks=1, flow=None, **overrides):
    inputs = dict(parallel_banks=banks,
                  total_mass_flow_kg_s=flow or reference.mass_flow_kg_s,
                  ua_flow_exponent=0.0, ua_condition_factor=1.0,
                  viscosity_ratio=1.0, density_ratio=1.0,
                  additional_pressure_drop_pa=0.0)
    inputs.update(overrides)
    return scale_bank(reference, **inputs)


def test_measured_anchors_and_parallel_replication():
    for reference in load_references(DATA):
        result = estimate(reference, banks=4, flow=4 * reference.mass_flow_kg_s)
        assert result["estimated_ua_w_k"] == pytest.approx(4 * reference.ua_w_k)
        assert result["estimated_tube_pressure_drop_pa"] == pytest.approx(reference.tube_pressure_drop_pa)
        assert result["sum_bank_envelope_m3"] == pytest.approx(0.000798)


def test_density_viscosity_and_flow_scaling():
    ref = load_references(DATA)[0]
    result = estimate(ref, flow=2 * ref.mass_flow_kg_s, viscosity_ratio=1.5,
                      density_ratio=3, additional_pressure_drop_pa=100)
    assert result["estimated_tube_and_additional_pressure_drop_pa"] == pytest.approx(5200)
    assert ref.tube_mean_temperature_k == pytest.approx((296.15 + 366.15) / 2)


@pytest.mark.parametrize("overrides", [dict(total_mass_flow_kg_s=0),
    dict(total_mass_flow_kg_s=-1), dict(density_ratio=float("nan")),
    dict(parallel_banks=1.5), dict(parallel_banks=True), dict(ua_flow_exponent=2)])
def test_unsupported_inputs_are_explicit(overrides):
    with pytest.raises(ValueError):
        estimate(load_references(DATA)[0], **overrides)
