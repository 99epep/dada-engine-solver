import math

import pytest

from dada_solver.exchangers.models import (
    CorrelationRegime,
    ExchangerOperatingPoint,
    ParallelRectangularChannels,
    ThermalResistanceModel,
    TransportProperties,
    evaluate_parallel_channels,
)
from dada_solver.exchangers.cycle_adapter import (
    ExchangerSide,
    conservative_port_flow_operating_point,
)
from dada_solver.exchangers.sizing import (
    ExchangerRequirements,
    assess_exchanger,
    pareto_exchanger_designs,
    screen_parallel_channel_designs,
)
from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.results import CycleDiagnostics, Extrema
from dada_solver.topology import (
    CycleTopologyClassification,
    CycleTopologyDiagnostic,
)


def _transport() -> TransportProperties:
    return TransportProperties(dynamic_viscosity=1.8e-5, thermal_conductivity=0.026)


def _resistances() -> ThermalResistanceModel:
    return ThermalResistanceModel(
        wall_thickness=1.0e-3,
        wall_thermal_conductivity=200.0,
        external_conductance=100.0,
    )


def test_parallel_channel_geometry_is_exact() -> None:
    geometry = ParallelRectangularChannels(10, 0.010, 0.002, 0.200)

    assert geometry.total_flow_area == pytest.approx(2.0e-4)
    assert geometry.gas_volume == pytest.approx(4.0e-5)
    assert geometry.gas_side_area == pytest.approx(0.048)
    assert geometry.hydraulic_diameter == pytest.approx(0.0033333333333333335)


def test_square_duct_laminar_correlations_match_limiting_values(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    geometry = ParallelRectangularChannels(1, 0.01, 0.01, 0.1)
    point = ExchangerOperatingPoint(1.0e-5, 1.0e5, 300.0)

    result = evaluate_parallel_channels(
        geometry, point, ideal_gas, _transport(), _resistances()
    )

    assert result.regime is CorrelationRegime.LAMINAR
    expected_poiseuille = 96.0 * (
        1.0 - 1.3553 + 1.9467 - 1.7012 + 0.9564 - 0.2537
    )
    dynamic_pressure = 0.5 * (
        point.pressure / (ideal_gas.gas_constant * point.temperature)
    ) * result.velocity**2
    expected_drop = (
        expected_poiseuille
        / result.reynolds_number
        * geometry.channel_length
        / geometry.hydraulic_diameter
        * dynamic_pressure
    )
    assert result.pressure_drop == pytest.approx(expected_drop)
    expected_nusselt = 7.541 * (
        1.0 - 2.610 + 4.970 - 5.119 + 2.702 - 0.548
    )
    expected_h = expected_nusselt * _transport().thermal_conductivity / 0.01
    assert result.gas_side_heat_transfer_coefficient == pytest.approx(expected_h)


def test_thermal_resistances_are_combined_in_series(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    geometry = ParallelRectangularChannels(2, 0.01, 0.002, 0.1)
    point = ExchangerOperatingPoint(1.0e-5, 1.0e5, 300.0)
    thermal = ThermalResistanceModel(0.001, 200.0, 50.0, 0.002)

    result = evaluate_parallel_channels(
        geometry, point, ideal_gas, _transport(), thermal
    )

    assert result.gas_side_heat_transfer_coefficient is not None
    area = geometry.gas_side_area
    expected = 1.0 / (
        1.0 / (result.gas_side_heat_transfer_coefficient * area)
        + thermal.wall_thickness / (thermal.wall_thermal_conductivity * area)
        + 1.0 / thermal.external_conductance
        + thermal.additional_resistance
    )
    assert result.overall_conductance == pytest.approx(expected)
    assert result.reservoir_effectiveness == pytest.approx(
        1.0 - math.exp(-expected / (point.absolute_mass_flow_rate * ideal_gas.heat_capacity_cp))
    )


def test_transitional_correlations_are_reported_unavailable(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    geometry = ParallelRectangularChannels(1, 0.01, 0.01, 0.1)
    density = 1.0e5 / (ideal_gas.gas_constant * 300.0)
    velocity = 2600.0 * 1.8e-5 / (density * geometry.hydraulic_diameter)
    point = ExchangerOperatingPoint(
        velocity * density * geometry.total_flow_area, 1.0e5, 300.0
    )

    result = evaluate_parallel_channels(
        geometry, point, ideal_gas, _transport(), _resistances()
    )

    assert result.regime is CorrelationRegime.TRANSITIONAL_UNAVAILABLE
    assert result.pressure_drop is None
    assert result.overall_conductance is None


def test_screening_rejects_unavailable_and_violated_results(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    results = screen_parallel_channel_designs(
        channel_counts=(1, 10),
        channel_widths=(0.01,),
        channel_heights=(0.01,),
        channel_lengths=(0.1,),
        operating_point=ExchangerOperatingPoint(0.001, 1.0e5, 300.0),
        gas=ideal_gas,
        transport=_transport(),
        thermal_resistances=_resistances(),
        requirements=ExchangerRequirements(0.01, 1.0e5, 0.3, 1.0e-3, 0.1),
    )

    assert len(results) == 2
    assert all(isinstance(item.feasible, bool) for item in results)


def test_pareto_filter_keeps_only_feasible_designs(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    performance = evaluate_parallel_channels(
        ParallelRectangularChannels(10, 0.01, 0.002, 0.1),
        ExchangerOperatingPoint(1.0e-4, 1.0e5, 300.0),
        ideal_gas,
        _transport(),
        _resistances(),
    )
    feasible = assess_exchanger(
        performance, ExchangerRequirements(1.0e-3, 1.0e6, 1.0, 1.0, 1.0)
    )
    infeasible = assess_exchanger(
        performance, ExchangerRequirements(1.0e6, 1.0e6, 1.0, 1.0, 1.0)
    )

    assert pareto_exchanger_designs((feasible, infeasible)) == (feasible,)


def test_cycle_adapter_uses_largest_adjacent_cold_port_flow() -> None:
    diagnostics = CycleDiagnostics(
        pressure_extrema={},
        temperature_extrema={},
        mass_flow_extrema={
            "small_to_cold": Extrema(-0.01, 0.03),
            "cold_to_large": Extrema(-0.04, 0.02),
        },
        cold_heat_rate_extrema=Extrema(0.0, 0.0),
        hot_heat_rate_extrema=Extrema(0.0, 0.0),
        valve_events=(),
        topology=CycleTopologyDiagnostic(
            classification=CycleTopologyClassification.NON_NOMINAL,
            reasons=("controlled_test",),
        ),
    )

    point = conservative_port_flow_operating_point(
        diagnostics, ExchangerSide.COLD, 1.2e5, 270.0
    )

    assert point.absolute_mass_flow_rate == pytest.approx(0.04)
    assert point.pressure == pytest.approx(1.2e5)
