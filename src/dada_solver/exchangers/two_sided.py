"""Two-sided channel heat-exchanger screening with an explicit liquid loop."""

from __future__ import annotations

from dataclasses import dataclass
import math

from dada_solver.exchangers.models import (
    CorrelationRegime,
    ExchangerPerformance,
    ParallelRectangularChannels,
    _darcy_friction_factor,
    _fully_developed_available,
    _nusselt_number,
    _regime,
)


@dataclass(frozen=True, slots=True)
class IncompressibleFluidProperties:
    density: float
    specific_heat: float
    dynamic_viscosity: float
    thermal_conductivity: float

    def __post_init__(self) -> None:
        for name, value in (
            ("Density", self.density),
            ("Specific heat", self.specific_heat),
            ("Dynamic viscosity", self.dynamic_viscosity),
            ("Thermal conductivity", self.thermal_conductivity),
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")


@dataclass(frozen=True, slots=True)
class LiquidLoopOperatingPoint:
    mass_flow_rate: float
    pump_efficiency: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.mass_flow_rate) or self.mass_flow_rate <= 0.0:
            raise ValueError("Liquid mass-flow rate must be finite and positive.")
        if (
            not math.isfinite(self.pump_efficiency)
            or not 0.0 < self.pump_efficiency <= 1.0
        ):
            raise ValueError("Pump efficiency must lie in (0, 1].")


@dataclass(frozen=True, slots=True)
class ChannelSidePerformance:
    reynolds_number: float
    prandtl_number: float
    velocity: float
    pressure_drop: float | None
    heat_transfer_coefficient: float | None
    regime: CorrelationRegime
    fully_developed_assumption_satisfied: bool


@dataclass(frozen=True, slots=True)
class TwoSidedExchangerPerformance:
    gas_side: ChannelSidePerformance
    liquid_side: ChannelSidePerformance
    shared_heat_transfer_area: float
    overall_conductance: float | None
    liquid_pump_power: float | None
    available: bool
    unavailable_reasons: tuple[str, ...]


def channel_side_from_gas_performance(
    performance: ExchangerPerformance,
) -> ChannelSidePerformance:
    """Extract the gas-side terms needed by the explicit two-sided model."""

    return ChannelSidePerformance(
        reynolds_number=performance.reynolds_number,
        prandtl_number=performance.prandtl_number,
        velocity=performance.velocity,
        pressure_drop=performance.pressure_drop,
        heat_transfer_coefficient=performance.gas_side_heat_transfer_coefficient,
        regime=performance.regime,
        fully_developed_assumption_satisfied=(
            performance.fully_developed_assumption_satisfied
        ),
    )


def evaluate_incompressible_channel_side(
    geometry: ParallelRectangularChannels,
    mass_flow_rate: float,
    fluid: IncompressibleFluidProperties,
) -> ChannelSidePerformance:
    """Evaluate one steady incompressible side with the channel correlations."""

    if not math.isfinite(mass_flow_rate) or mass_flow_rate <= 0.0:
        raise ValueError("Mass-flow rate must be finite and positive.")
    velocity = mass_flow_rate / (fluid.density * geometry.total_flow_area)
    diameter = geometry.hydraulic_diameter
    reynolds = fluid.density * velocity * diameter / fluid.dynamic_viscosity
    prandtl = fluid.specific_heat * fluid.dynamic_viscosity / fluid.thermal_conductivity
    regime = _regime(reynolds, prandtl)
    developed = _fully_developed_available(
        geometry.channel_length, diameter, reynolds, prandtl, regime
    )
    friction = _darcy_friction_factor(
        reynolds,
        geometry.aspect_ratio,
        geometry.surface_roughness / diameter,
        regime,
    )
    pressure_drop = None
    if friction is not None and developed:
        pressure_drop = (
            friction * geometry.channel_length / diameter
            + geometry.minor_loss_coefficient
        ) * 0.5 * fluid.density * velocity**2
    nusselt = _nusselt_number(
        reynolds, prandtl, geometry.aspect_ratio, friction, regime
    )
    coefficient = (
        None
        if nusselt is None or not developed
        else nusselt * fluid.thermal_conductivity / diameter
    )
    return ChannelSidePerformance(
        reynolds_number=reynolds,
        prandtl_number=prandtl,
        velocity=velocity,
        pressure_drop=pressure_drop,
        heat_transfer_coefficient=coefficient,
        regime=regime,
        fully_developed_assumption_satisfied=developed,
    )


def combine_channel_sides(
    gas_side: ChannelSidePerformance,
    liquid_side: ChannelSidePerformance,
    liquid_operating_point: LiquidLoopOperatingPoint,
    *,
    liquid_density: float,
    shared_heat_transfer_area: float,
    wall_thickness: float,
    wall_thermal_conductivity: float,
    additional_resistance: float = 0.0,
) -> TwoSidedExchangerPerformance:
    """Combine already evaluated sides without inventing external conductance."""

    for name, value in (
        ("Shared heat-transfer area", shared_heat_transfer_area),
        ("Wall thermal conductivity", wall_thermal_conductivity),
        ("Liquid density", liquid_density),
    ):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive.")
    if not math.isfinite(wall_thickness) or wall_thickness < 0.0:
        raise ValueError("Wall thickness must be finite and non-negative.")
    if not math.isfinite(additional_resistance) or additional_resistance < 0.0:
        raise ValueError("Additional resistance must be finite and non-negative.")
    unavailable = []
    if gas_side.heat_transfer_coefficient is None:
        unavailable.append("gas_side_heat_transfer")
    if liquid_side.heat_transfer_coefficient is None:
        unavailable.append("liquid_side_heat_transfer")
    if liquid_side.pressure_drop is None:
        unavailable.append("liquid_side_pressure_drop")
    if unavailable:
        return TwoSidedExchangerPerformance(
            gas_side, liquid_side, shared_heat_transfer_area, None, None,
            False, tuple(unavailable)
        )
    assert gas_side.heat_transfer_coefficient is not None
    assert liquid_side.heat_transfer_coefficient is not None
    assert liquid_side.pressure_drop is not None
    resistance = (
        1.0 / (gas_side.heat_transfer_coefficient * shared_heat_transfer_area)
        + wall_thickness / (wall_thermal_conductivity * shared_heat_transfer_area)
        + 1.0 / (liquid_side.heat_transfer_coefficient * shared_heat_transfer_area)
        + additional_resistance
    )
    conductance = 1.0 / resistance
    volumetric_flow = liquid_operating_point.mass_flow_rate / liquid_density
    return TwoSidedExchangerPerformance(
        gas_side, liquid_side, shared_heat_transfer_area, conductance,
        volumetric_flow * liquid_side.pressure_drop / liquid_operating_point.pump_efficiency,
        True, ()
    )


def evaluate_two_sided_exchanger(
    gas_side: ChannelSidePerformance,
    liquid_geometry: ParallelRectangularChannels,
    liquid_fluid: IncompressibleFluidProperties,
    liquid_operating_point: LiquidLoopOperatingPoint,
    *,
    shared_heat_transfer_area: float,
    wall_thickness: float,
    wall_thermal_conductivity: float,
    additional_resistance: float = 0.0,
) -> TwoSidedExchangerPerformance:
    """Evaluate the liquid side, overall UA and liquid hydraulic input."""

    liquid_side = evaluate_incompressible_channel_side(
        liquid_geometry, liquid_operating_point.mass_flow_rate, liquid_fluid
    )
    combined = combine_channel_sides(
        gas_side,
        liquid_side,
        liquid_operating_point,
        liquid_density=liquid_fluid.density,
        shared_heat_transfer_area=shared_heat_transfer_area,
        wall_thickness=wall_thickness,
        wall_thermal_conductivity=wall_thermal_conductivity,
        additional_resistance=additional_resistance,
    )
    return combined
