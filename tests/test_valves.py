import pytest

from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.hydraulics import CompressibleOrifice
from dada_solver.valves import PassiveCheckValve, ValveState


@pytest.fixture
def valve() -> PassiveCheckValve:
    return PassiveCheckValve(
        flow_model=CompressibleOrifice(effective_flow_area=1.0e-6),
        opening_pressure_difference=2.0e3,
        closing_pressure_difference=0.5e3,
    )


def test_valve_hysteresis_preserves_state_inside_threshold_band(
    valve: PassiveCheckValve,
) -> None:
    closed = valve.updated_state(ValveState.CLOSED, 101.0e3, 100.0e3)
    opened = valve.updated_state(ValveState.OPEN, 101.0e3, 100.0e3)

    assert closed.current_state is ValveState.CLOSED
    assert opened.current_state is ValveState.OPEN


def test_valve_opens_and_closes_only_from_pressure_thresholds(
    valve: PassiveCheckValve,
) -> None:
    opening = valve.updated_state(ValveState.CLOSED, 102.0e3, 100.0e3)
    closing = valve.updated_state(ValveState.OPEN, 100.5e3, 100.0e3)

    assert opening.current_state is ValveState.OPEN
    assert opening.changed
    assert closing.current_state is ValveState.CLOSED
    assert closing.changed


def test_check_valve_blocks_reverse_flow_even_while_open(
    valve: PassiveCheckValve,
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    result = valve.flow(
        ValveState.OPEN,
        upstream_pressure=1.0e5,
        downstream_pressure=1.2e5,
        upstream_temperature=300.0,
        gas=ideal_gas,
    )

    assert result.mass_flow_rate == 0.0


def test_zero_hysteresis_valve_is_closed_at_exact_pressure_equality() -> None:
    ideal = PassiveCheckValve(CompressibleOrifice(1.0e-6), 0.0, 0.0)

    equal = ideal.updated_state(ValveState.CLOSED, 1.0e5, 1.0e5)
    forward = ideal.updated_state(ValveState.CLOSED, 1.0e5 + 1.0, 1.0e5)
    reseated = ideal.updated_state(ValveState.OPEN, 1.0e5, 1.0e5)

    assert equal.current_state is ValveState.CLOSED
    assert forward.current_state is ValveState.OPEN
    assert reseated.current_state is ValveState.CLOSED
