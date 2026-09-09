import pytest

from dada_solver.heat_transfer import ReservoirHeatTransfer


def test_zero_heat_transfer_at_reservoir_temperature() -> None:
    exchanger = ReservoirHeatTransfer(
        conductance=125.0,
        reservoir_temperature=290.0,
    )

    assert exchanger.heat_rate(290.0) == pytest.approx(0.0)


def test_heat_rate_sign_is_heat_received_by_gas() -> None:
    exchanger = ReservoirHeatTransfer(
        conductance=10.0,
        reservoir_temperature=300.0,
    )

    assert exchanger.heat_rate(280.0) > 0.0
    assert exchanger.heat_rate(320.0) < 0.0

