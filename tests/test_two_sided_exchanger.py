import pytest

from dada_solver.exchangers.models import CorrelationRegime, ParallelRectangularChannels
from dada_solver.exchangers.two_sided import (
    ChannelSidePerformance,
    IncompressibleFluidProperties,
    LiquidLoopOperatingPoint,
    combine_channel_sides,
    evaluate_incompressible_channel_side,
)


def test_incompressible_channel_side_conserves_declared_flow() -> None:
    geometry = ParallelRectangularChannels(10, 0.01, 0.002, 3.0)
    fluid = IncompressibleFluidProperties(1000.0, 4200.0, 1.0e-3, 0.6)
    result = evaluate_incompressible_channel_side(geometry, 0.1, fluid)

    assert result.velocity == pytest.approx(0.5)
    assert result.reynolds_number > 0.0
    assert result.pressure_drop is not None
    assert result.heat_transfer_coefficient is not None


def test_two_sided_resistances_and_pump_power_are_explicit() -> None:
    side = ChannelSidePerformance(
        reynolds_number=1000.0,
        prandtl_number=7.0,
        velocity=1.0,
        pressure_drop=2000.0,
        heat_transfer_coefficient=1000.0,
        regime=CorrelationRegime.LAMINAR,
        fully_developed_assumption_satisfied=True,
    )
    operating = LiquidLoopOperatingPoint(0.2, 0.5)
    result = combine_channel_sides(
        side,
        side,
        operating,
        liquid_density=1000.0,
        shared_heat_transfer_area=2.0,
        wall_thickness=0.0,
        wall_thermal_conductivity=200.0,
    )

    assert result.overall_conductance == pytest.approx(1000.0)
    assert result.liquid_pump_power == pytest.approx(0.8)
