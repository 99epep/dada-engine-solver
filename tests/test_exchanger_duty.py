import numpy as np
import pytest

from dada_solver.exchangers.duty import extract_exchanger_duty, summarize_port
from tests.test_periodic import create_static_cycle


def test_signed_mass_and_activity_use_time_not_sample_counts():
    # One symmetric reversal from -2 to +2 kg/s over four seconds.
    coarse = summarize_port([0, 4], [-2, 2], 0.5)
    refined = summarize_port([0, 0.2, 1, 2, 3, 4], [-2, -1.8, -1, 0, 1, 2], 0.5)
    assert coarse == pytest.approx(refined)
    assert coarse["forward_mass_per_cycle_kg"] == pytest.approx(2.0)
    assert coarse["reverse_mass_per_cycle_kg"] == pytest.approx(2.0)
    assert coarse["duration_above_activity_threshold_s"] == pytest.approx(1.0)
    assert coarse["duty_fraction_above_activity_threshold"] == pytest.approx(0.25)


def test_interrupted_forward_pulse_preserves_zero_flow_and_duplicate_event_samples():
    summary = summarize_port([0, 1, 1, 2, 3, 4], [0, 0, 0, 2, 0, 0], 0.5)
    assert summary["forward_mass_per_cycle_kg"] == pytest.approx(2.0)
    assert summary["reverse_mass_per_cycle_kg"] == 0.0
    assert summary["duration_above_activity_threshold_s"] == pytest.approx(1.0)


def test_zero_flow_is_available_without_inventing_a_pulse():
    summary = summarize_port([0, 1, 2], [0, 0, 0])
    assert summary["forward_mass_per_cycle_kg"] == 0.0
    assert summary["duration_above_activity_threshold_s"] == 0.0
    assert summary["reverse_to_forward_mass_ratio"] is None


def test_duty_export_preserves_simultaneous_state_and_port_columns(ideal_gas):
    model, _, cycle = create_static_cycle(ideal_gas)
    summary, columns = extract_exchanger_duty(cycle, model)
    assert summary["period_s"] == pytest.approx(2*np.pi/model.angular_speed)
    assert summary["intended_flow"] == "intermittent_unidirectional"
    assert columns["H_i_pressure_Pa"] == pytest.approx(2e5)
    assert columns["H_o_temperature_K"] == pytest.approx(300.0)
    assert columns["H_i_to_L_mass_flow_kg_s"] == pytest.approx(0.0)
    assert summary["exchangers"]["H_i"]["mean_heat_received_W"] == pytest.approx(0.0)


@pytest.mark.parametrize("time,flow", [([0, 0], [1, 1]), ([1, 0], [1, 1]), ([0, 1], [1]), ([0, 1], [1, np.nan])])
def test_invalid_histories_are_rejected(time, flow):
    with pytest.raises(ValueError):
        summarize_port(time, flow)
