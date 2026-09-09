"""Low-order geometric models for parallel-channel heat exchangers.

The closures in this module are screening models. They do not replace the
four-control-volume cycle model and they do not resolve manifolds, entrance
regions, conjugate conduction, pulsatile-flow phase lag, or axial wall storage.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math

from dada_solver.fluids import CaloricallyPerfectGas


def _positive(name: str, value: float) -> None:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive.")


class CorrelationRegime(Enum):
    """Flow regime used by the channel correlations."""

    LAMINAR = "laminar"
    TRANSITIONAL_UNAVAILABLE = "transitional_unavailable"
    TURBULENT = "turbulent"
    OUTSIDE_VALIDITY = "outside_validity"


@dataclass(frozen=True, slots=True)
class TransportProperties:
    """Gas transport properties evaluated at the declared mean state."""

    dynamic_viscosity: float
    thermal_conductivity: float

    def __post_init__(self) -> None:
        _positive("Dynamic viscosity", self.dynamic_viscosity)
        _positive("Thermal conductivity", self.thermal_conductivity)


@dataclass(frozen=True, slots=True)
class ParallelRectangularChannels:
    """Identical straight rectangular gas passages connected in parallel."""

    channel_count: int
    channel_width: float
    channel_height: float
    channel_length: float
    surface_roughness: float = 0.0
    minor_loss_coefficient: float = 0.0

    def __post_init__(self) -> None:
        if (
            isinstance(self.channel_count, bool)
            or not isinstance(self.channel_count, int)
            or self.channel_count <= 0
        ):
            raise ValueError("Channel count must be a positive integer.")
        for name, value in (
            ("Channel width", self.channel_width),
            ("Channel height", self.channel_height),
            ("Channel length", self.channel_length),
        ):
            _positive(name, value)
        for name, value in (
            ("Surface roughness", self.surface_roughness),
            ("Minor-loss coefficient", self.minor_loss_coefficient),
        ):
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and non-negative.")

    @property
    def single_channel_area(self) -> float:
        return self.channel_width * self.channel_height

    @property
    def total_flow_area(self) -> float:
        return self.channel_count * self.single_channel_area

    @property
    def wetted_perimeter(self) -> float:
        return 2.0 * (self.channel_width + self.channel_height)

    @property
    def hydraulic_diameter(self) -> float:
        return 4.0 * self.single_channel_area / self.wetted_perimeter

    @property
    def gas_volume(self) -> float:
        return self.total_flow_area * self.channel_length

    @property
    def gas_side_area(self) -> float:
        return self.channel_count * self.wetted_perimeter * self.channel_length

    @property
    def aspect_ratio(self) -> float:
        return min(self.channel_width, self.channel_height) / max(
            self.channel_width, self.channel_height
        )


@dataclass(frozen=True, slots=True)
class ThermalResistanceModel:
    """Wall and reservoir-side terms in the overall thermal resistance."""

    wall_thickness: float
    wall_thermal_conductivity: float
    external_conductance: float
    additional_resistance: float = 0.0

    def __post_init__(self) -> None:
        if not math.isfinite(self.wall_thickness) or self.wall_thickness < 0.0:
            raise ValueError("Wall thickness must be finite and non-negative.")
        _positive("Wall thermal conductivity", self.wall_thermal_conductivity)
        _positive("External conductance", self.external_conductance)
        if not math.isfinite(self.additional_resistance) or self.additional_resistance < 0.0:
            raise ValueError("Additional resistance must be finite and non-negative.")


@dataclass(frozen=True, slots=True)
class ExchangerOperatingPoint:
    """Mean-state operating point for a conservative screening calculation."""

    absolute_mass_flow_rate: float
    pressure: float
    temperature: float

    def __post_init__(self) -> None:
        for name, value in (
            ("Absolute mass-flow rate", self.absolute_mass_flow_rate),
            ("Pressure", self.pressure),
            ("Temperature", self.temperature),
        ):
            _positive(name, value)


@dataclass(frozen=True, slots=True)
class ExchangerPerformance:
    """Predicted performance and explicit correlation-availability state."""

    geometry: ParallelRectangularChannels
    reynolds_number: float
    prandtl_number: float
    mach_number: float
    velocity: float
    pressure_drop: float | None
    gas_side_heat_transfer_coefficient: float | None
    overall_conductance: float | None
    reservoir_effectiveness: float | None
    regime: CorrelationRegime
    thermal_correlation_available: bool
    hydraulic_correlation_available: bool
    pressure_drop_fraction: float | None
    fully_developed_assumption_satisfied: bool
    assumptions: tuple[str, ...]


def evaluate_parallel_channels(
    geometry: ParallelRectangularChannels,
    operating_point: ExchangerOperatingPoint,
    gas: CaloricallyPerfectGas,
    transport: TransportProperties,
    thermal_resistances: ThermalResistanceModel,
) -> ExchangerPerformance:
    """Evaluate one straight-channel design at one mean operating point.

    Laminar rectangular-duct correlations are fully developed values. The
    turbulent closure uses the smooth-duct Gnielinski correlation. The
    transitional interval is deliberately reported as unavailable.
    """

    density = operating_point.pressure / (
        gas.gas_constant * operating_point.temperature
    )
    velocity = operating_point.absolute_mass_flow_rate / (
        density * geometry.total_flow_area
    )
    diameter = geometry.hydraulic_diameter
    reynolds = density * velocity * diameter / transport.dynamic_viscosity
    prandtl = (
        gas.heat_capacity_cp
        * transport.dynamic_viscosity
        / transport.thermal_conductivity
    )
    sound_speed = math.sqrt(
        gas.heat_capacity_ratio * gas.gas_constant * operating_point.temperature
    )
    mach = velocity / sound_speed
    regime = _regime(reynolds, prandtl)
    fully_developed = _fully_developed_available(
        geometry.channel_length, diameter, reynolds, prandtl, regime
    )

    friction_factor = _darcy_friction_factor(
        reynolds,
        geometry.aspect_ratio,
        geometry.surface_roughness / diameter,
        regime,
    )
    pressure_drop = None
    pressure_fraction = None
    if friction_factor is not None and fully_developed:
        dynamic_pressure = 0.5 * density * velocity * velocity
        pressure_drop = (
            friction_factor * geometry.channel_length / diameter
            + geometry.minor_loss_coefficient
        ) * dynamic_pressure
        pressure_fraction = pressure_drop / operating_point.pressure

    nusselt = _nusselt_number(
        reynolds, prandtl, geometry.aspect_ratio, friction_factor, regime
    )
    heat_coefficient = None
    overall_conductance = None
    effectiveness = None
    if nusselt is not None and fully_developed:
        heat_coefficient = nusselt * transport.thermal_conductivity / diameter
        area = geometry.gas_side_area
        total_resistance = (
            1.0 / (heat_coefficient * area)
            + thermal_resistances.wall_thickness
            / (thermal_resistances.wall_thermal_conductivity * area)
            + 1.0 / thermal_resistances.external_conductance
            + thermal_resistances.additional_resistance
        )
        overall_conductance = 1.0 / total_resistance
        capacity_rate = (
            operating_point.absolute_mass_flow_rate * gas.heat_capacity_cp
        )
        effectiveness = -math.expm1(-overall_conductance / capacity_rate)

    return ExchangerPerformance(
        geometry=geometry,
        reynolds_number=reynolds,
        prandtl_number=prandtl,
        mach_number=mach,
        velocity=velocity,
        pressure_drop=pressure_drop,
        gas_side_heat_transfer_coefficient=heat_coefficient,
        overall_conductance=overall_conductance,
        reservoir_effectiveness=effectiveness,
        regime=regime,
        thermal_correlation_available=nusselt is not None and fully_developed,
        hydraulic_correlation_available=friction_factor is not None and fully_developed,
        pressure_drop_fraction=pressure_fraction,
        fully_developed_assumption_satisfied=fully_developed,
        assumptions=(
            "steady mean-state flow",
            "uniform flow distribution among channels",
            "fully developed internal flow",
            "constant properties",
            "negligible axial wall conduction and wall heat capacity",
            "external conductance supplied independently",
        ),
    )


def _regime(reynolds: float, prandtl: float) -> CorrelationRegime:
    if reynolds < 2300.0:
        return CorrelationRegime.LAMINAR
    if reynolds < 4000.0:
        return CorrelationRegime.TRANSITIONAL_UNAVAILABLE
    if reynolds <= 5.0e6 and 0.5 <= prandtl <= 2000.0:
        return CorrelationRegime.TURBULENT
    return CorrelationRegime.OUTSIDE_VALIDITY


def _fully_developed_available(
    length: float,
    diameter: float,
    reynolds: float,
    prandtl: float,
    regime: CorrelationRegime,
) -> bool:
    if regime is CorrelationRegime.LAMINAR:
        hydrodynamic_entrance_length = 0.05 * reynolds * diameter
        thermal_entrance_length = 0.05 * reynolds * prandtl * diameter
        return length >= max(hydrodynamic_entrance_length, thermal_entrance_length)
    if regime is CorrelationRegime.TURBULENT:
        return length >= 10.0 * diameter
    return False


def _rectangular_polynomial(aspect_ratio: float, coefficients: tuple[float, ...]) -> float:
    return sum(value * aspect_ratio**power for power, value in enumerate(coefficients))


def _darcy_friction_factor(
    reynolds: float,
    aspect_ratio: float,
    relative_roughness: float,
    regime: CorrelationRegime,
) -> float | None:
    if regime is CorrelationRegime.LAMINAR:
        poiseuille = 96.0 * _rectangular_polynomial(
            aspect_ratio,
            (1.0, -1.3553, 1.9467, -1.7012, 0.9564, -0.2537),
        )
        return poiseuille / reynolds
    if regime is CorrelationRegime.TURBULENT:
        return (
            -1.8
            * math.log10((relative_roughness / 3.7) ** 1.11 + 6.9 / reynolds)
        ) ** -2
    return None


def _nusselt_number(
    reynolds: float,
    prandtl: float,
    aspect_ratio: float,
    friction_factor: float | None,
    regime: CorrelationRegime,
) -> float | None:
    if regime is CorrelationRegime.LAMINAR:
        return 7.541 * _rectangular_polynomial(
            aspect_ratio,
            (1.0, -2.610, 4.970, -5.119, 2.702, -0.548),
        )
    if regime is CorrelationRegime.TURBULENT:
        assert friction_factor is not None
        friction = friction_factor
        numerator = (friction / 8.0) * (reynolds - 1000.0) * prandtl
        denominator = 1.0 + 12.7 * math.sqrt(friction / 8.0) * (
            prandtl ** (2.0 / 3.0) - 1.0
        )
        return numerator / denominator
    return None
