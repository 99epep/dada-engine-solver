import math

import pytest

from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.hydraulics import (
    BlendedHydraulicFlow,
    CompressibleOrifice,
    QuasiSteadyCompressibleDuct,
    RectangularDuct,
    SeriesDuctOrifice,
)


def test_orifice_has_zero_flow_at_equal_pressure(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    orifice = CompressibleOrifice(effective_flow_area=2.0e-6)

    result = orifice.directed_flow(2.0e5, 2.0e5, 300.0, ideal_gas)

    assert result.mass_flow_rate == 0.0
    assert not result.is_choked


def test_orifice_matches_unchoked_equation(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    cda = 1.5e-6
    upstream_pressure = 2.0e5
    downstream_pressure = 1.8e5
    temperature = 310.0
    orifice = CompressibleOrifice(cda)
    gamma = ideal_gas.heat_capacity_ratio
    ratio = downstream_pressure / upstream_pressure
    expected = cda * upstream_pressure * math.sqrt(
        2.0
        * gamma
        / (ideal_gas.gas_constant * temperature * (gamma - 1.0))
        * (ratio ** (2.0 / gamma) - ratio ** ((gamma + 1.0) / gamma))
    )

    result = orifice.directed_flow(
        upstream_pressure, downstream_pressure, temperature, ideal_gas
    )

    assert not result.is_choked
    assert result.mass_flow_rate == pytest.approx(expected)


def test_regularized_orifice_is_linear_near_equal_pressure(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    orifice = CompressibleOrifice(1.0e-4, pressure_regularization=1.0)
    pressure = 1.0e5
    first = orifice.directed_flow(pressure, pressure - 1.0e-4, 300.0, ideal_gas)
    second = orifice.directed_flow(pressure, pressure - 2.0e-4, 300.0, ideal_gas)

    assert second.mass_flow_rate / first.mass_flow_rate == pytest.approx(
        2.0, rel=2.0e-4
    )


def test_zero_regularization_retains_published_orifice_law(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    implicit_default = CompressibleOrifice(2.0e-6)
    explicit_zero = CompressibleOrifice(2.0e-6, pressure_regularization=0.0)

    assert implicit_default.directed_flow(1.0e5, 99999.0, 300.0, ideal_gas) == (
        explicit_zero.directed_flow(1.0e5, 99999.0, 300.0, ideal_gas)
    )


def test_choked_flow_is_independent_of_downstream_pressure(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    orifice = CompressibleOrifice(effective_flow_area=1.0e-6)

    first = orifice.directed_flow(4.0e5, 1.0e5, 300.0, ideal_gas)
    second = orifice.directed_flow(4.0e5, 0.5e5, 300.0, ideal_gas)

    assert first.is_choked
    assert second.is_choked
    assert first.mass_flow_rate == pytest.approx(second.mass_flow_rate)


def test_bidirectional_orifice_uses_actual_upstream_temperature(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    orifice = CompressibleOrifice(effective_flow_area=1.0e-6)

    signed = orifice.bidirectional_flow(
        first_pressure=1.0e5,
        second_pressure=2.0e5,
        first_temperature=250.0,
        second_temperature=350.0,
        gas=ideal_gas,
    )
    expected = orifice.directed_flow(2.0e5, 1.0e5, 350.0, ideal_gas)

    assert signed.mass_flow_rate == pytest.approx(-expected.mass_flow_rate)


def _duct() -> QuasiSteadyCompressibleDuct:
    return QuasiSteadyCompressibleDuct(
        RectangularDuct(
            flow_area=1.0e-4,
            hydraulic_diameter=2.0e-3,
            length=0.1,
            aspect_ratio=0.2,
            dynamic_viscosity=1.8e-5,
            minor_loss_coefficient=1.0,
        )
    )


def test_geometric_duct_flow_is_zero_at_equal_pressure(ideal_gas) -> None:
    result = _duct().directed_flow(1.0e5, 1.0e5, 300.0, ideal_gas)

    assert result.mass_flow_rate == 0.0
    assert not result.is_choked


def test_geometric_duct_flow_increases_with_pressure_difference(ideal_gas) -> None:
    low = _duct().directed_flow(1.01e5, 1.0e5, 300.0, ideal_gas)
    high = _duct().directed_flow(1.10e5, 1.0e5, 300.0, ideal_gas)

    assert 0.0 < low.mass_flow_rate < high.mass_flow_rate


def test_geometric_duct_bidirectional_flow_is_antisymmetric(ideal_gas) -> None:
    forward = _duct().bidirectional_flow(
        1.1e5, 1.0e5, 300.0, 300.0, ideal_gas
    )
    reverse = _duct().bidirectional_flow(
        1.0e5, 1.1e5, 300.0, 300.0, ideal_gas
    )

    assert reverse.mass_flow_rate == pytest.approx(-forward.mass_flow_rate)


def test_series_orifice_reduces_duct_flow(ideal_gas) -> None:
    duct = _duct()
    series = SeriesDuctOrifice(duct, CompressibleOrifice(5.0e-5))

    duct_only = duct.directed_flow(1.2e5, 1.0e5, 300.0, ideal_gas)
    combined = series.directed_flow(1.2e5, 1.0e5, 300.0, ideal_gas)

    assert 0.0 < combined.mass_flow_rate < duct_only.mass_flow_rate


def test_hydraulic_continuation_interpolates_flow_without_changing_endpoints(
    ideal_gas,
) -> None:
    initial = CompressibleOrifice(1.0e-6)
    final = CompressibleOrifice(3.0e-6)
    initial_flow = initial.directed_flow(1.2e5, 1.0e5, 300.0, ideal_gas)
    final_flow = final.directed_flow(1.2e5, 1.0e5, 300.0, ideal_gas)

    start = BlendedHydraulicFlow(initial, final, 0.0).directed_flow(
        1.2e5, 1.0e5, 300.0, ideal_gas
    )
    middle = BlendedHydraulicFlow(initial, final, 0.25).directed_flow(
        1.2e5, 1.0e5, 300.0, ideal_gas
    )
    end = BlendedHydraulicFlow(initial, final, 1.0).directed_flow(
        1.2e5, 1.0e5, 300.0, ideal_gas
    )

    assert start.mass_flow_rate == pytest.approx(initial_flow.mass_flow_rate)
    assert end.mass_flow_rate == pytest.approx(final_flow.mass_flow_rate)
    assert middle.mass_flow_rate == pytest.approx(
        0.75 * initial_flow.mass_flow_rate + 0.25 * final_flow.mass_flow_rate
    )
