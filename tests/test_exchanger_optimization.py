import pytest

from dada_solver.exchangers.models import (
    ExchangerOperatingPoint,
    ThermalResistanceModel,
    TransportProperties,
)
from dada_solver.exchangers.optimization import (
    ContinuousGeometryBounds,
    ExchangerObjective,
    optimize_parallel_channel_geometry,
)
from dada_solver.exchangers.sizing import ExchangerRequirements
from dada_solver.fluids import CaloricallyPerfectGas


def test_optimizer_returns_feasible_minimum_volume_candidate(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    requirements = ExchangerRequirements(
        minimum_overall_conductance=5.0,
        maximum_pressure_drop=1.0e4,
        maximum_mach_number=0.3,
        maximum_gas_volume=1.0e-3,
        maximum_pressure_drop_fraction=0.1,
    )

    result = optimize_parallel_channel_geometry(
        channel_counts=(20, 40),
        bounds=ContinuousGeometryBounds(0.01, 0.02, 0.001, 0.002, 0.1, 0.2),
        seed_widths=(0.01, 0.02),
        seed_heights=(0.001, 0.002),
        seed_lengths=(0.1, 0.2),
        operating_point=ExchangerOperatingPoint(0.001, 1.0e5, 300.0),
        gas=ideal_gas,
        transport=TransportProperties(1.8e-5, 0.026),
        thermal_resistances=ThermalResistanceModel(0.001, 200.0, 100.0),
        requirements=requirements,
        objective=ExchangerObjective.MINIMUM_GAS_VOLUME,
        maximum_iterations=50,
    )

    assert result.success
    assert result.best is not None
    assert result.best.feasible
    assert result.best.performance.overall_conductance is not None
    assert result.best.performance.overall_conductance >= 5.0


def test_optimizer_reports_absence_of_feasible_seed(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    result = optimize_parallel_channel_geometry(
        channel_counts=(1,),
        bounds=ContinuousGeometryBounds(0.01, 0.02, 0.001, 0.002, 0.1, 0.2),
        seed_widths=(0.01, 0.02),
        seed_heights=(0.001, 0.002),
        seed_lengths=(0.1, 0.2),
        operating_point=ExchangerOperatingPoint(0.1, 1.0e5, 300.0),
        gas=ideal_gas,
        transport=TransportProperties(1.8e-5, 0.026),
        thermal_resistances=ThermalResistanceModel(0.001, 200.0, 100.0),
        requirements=ExchangerRequirements(1.0e6, 1.0, 0.01, 1.0e-6, 1.0e-5),
    )

    assert not result.success
    assert result.best is None
    assert result.local_optimization_count == 0
