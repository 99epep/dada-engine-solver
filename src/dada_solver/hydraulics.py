"""Replaceable hydraulic-flow closures for compressible gas links."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Protocol

from dada_solver.fluids import CaloricallyPerfectGas


@dataclass(frozen=True, slots=True)
class FlowResult:
    """Mass-flow result, positive in the requested upstream-to-downstream direction."""

    mass_flow_rate: float
    is_choked: bool


class HydraulicFlowModel(Protocol):
    """Interface for a directed hydraulic-flow calculation."""

    def directed_flow(
        self,
        upstream_pressure: float,
        downstream_pressure: float,
        upstream_temperature: float,
        gas: CaloricallyPerfectGas,
    ) -> FlowResult:
        """Return a non-negative directed mass flow."""

    def bidirectional_flow(
        self,
        first_pressure: float,
        second_pressure: float,
        first_temperature: float,
        second_temperature: float,
        gas: CaloricallyPerfectGas,
    ) -> FlowResult:
        """Return signed flow, positive from the first node to the second."""


@dataclass(frozen=True, slots=True)
class HydraulicNetworkModels:
    """Replaceable flow closures for all four thermodynamic network links."""

    large_to_hot: HydraulicFlowModel
    small_to_cold: HydraulicFlowModel
    hot_to_small_valve: HydraulicFlowModel
    cold_to_large_valve: HydraulicFlowModel


@dataclass(frozen=True, slots=True)
class BlendedHydraulicFlow:
    """Numerical continuation between two closures; endpoints retain meaning."""

    initial: HydraulicFlowModel
    final: HydraulicFlowModel
    final_fraction: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.final_fraction) or not 0.0 <= self.final_fraction <= 1.0:
            raise ValueError("Final hydraulic fraction must lie in [0, 1].")

    def directed_flow(
        self, upstream_pressure, downstream_pressure, upstream_temperature, gas
    ) -> FlowResult:
        initial = self.initial.directed_flow(
            upstream_pressure, downstream_pressure, upstream_temperature, gas
        )
        final = self.final.directed_flow(
            upstream_pressure, downstream_pressure, upstream_temperature, gas
        )
        fraction = self.final_fraction
        return FlowResult(
            (1.0 - fraction) * initial.mass_flow_rate
            + fraction * final.mass_flow_rate,
            initial.is_choked or final.is_choked,
        )

    def bidirectional_flow(
        self,
        first_pressure,
        second_pressure,
        first_temperature,
        second_temperature,
        gas,
    ) -> FlowResult:
        initial = self.initial.bidirectional_flow(
            first_pressure, second_pressure, first_temperature, second_temperature, gas
        )
        final = self.final.bidirectional_flow(
            first_pressure, second_pressure, first_temperature, second_temperature, gas
        )
        fraction = self.final_fraction
        return FlowResult(
            (1.0 - fraction) * initial.mass_flow_rate
            + fraction * final.mass_flow_rate,
            initial.is_choked or final.is_choked,
        )


@dataclass(frozen=True, slots=True)
class RectangularDuct:
    """Straight rectangular hydraulic passage with explicit minor losses."""

    flow_area: float
    hydraulic_diameter: float
    length: float
    aspect_ratio: float
    dynamic_viscosity: float
    relative_roughness: float = 0.0
    minor_loss_coefficient: float = 0.0

    def __post_init__(self) -> None:
        for name, value in (
            ("Flow area", self.flow_area),
            ("Hydraulic diameter", self.hydraulic_diameter),
            ("Length", self.length),
            ("Dynamic viscosity", self.dynamic_viscosity),
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")
        if not 0.0 < self.aspect_ratio <= 1.0:
            raise ValueError("Aspect ratio must lie in (0, 1].")
        for name, value in (
            ("Relative roughness", self.relative_roughness),
            ("Minor-loss coefficient", self.minor_loss_coefficient),
        ):
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and non-negative.")


@dataclass(frozen=True, slots=True)
class QuasiSteadyCompressibleDuct:
    """Mean-density Darcy closure with an isentropic sonic mass-flow cap."""

    duct: RectangularDuct
    iteration_count: int = 48

    def directed_flow(
        self,
        upstream_pressure: float,
        downstream_pressure: float,
        upstream_temperature: float,
        gas: CaloricallyPerfectGas,
    ) -> FlowResult:
        for name, value in (
            ("upstream pressure", upstream_pressure),
            ("downstream pressure", downstream_pressure),
            ("upstream temperature", upstream_temperature),
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")
        if downstream_pressure >= upstream_pressure:
            return FlowResult(0.0, False)
        available_drop = upstream_pressure - downstream_pressure
        density = 0.5 * (upstream_pressure + downstream_pressure) / (
            gas.gas_constant * upstream_temperature
        )
        sonic_cap = CompressibleOrifice(self.duct.flow_area).directed_flow(
            upstream_pressure, downstream_pressure, upstream_temperature, gas
        )

        aspect = self.duct.aspect_ratio
        poiseuille = 96.0 * (
            1.0
            - 1.3553 * aspect
            + 1.9467 * aspect**2
            - 1.7012 * aspect**3
            + 0.9564 * aspect**4
            - 0.2537 * aspect**5
        )
        linear = (
            poiseuille
            * self.duct.dynamic_viscosity
            * self.duct.length
            / (2.0 * density * self.duct.flow_area * self.duct.hydraulic_diameter**2)
        )
        quadratic = self.duct.minor_loss_coefficient / (
            2.0 * density * self.duct.flow_area**2
        )
        if quadratic == 0.0:
            laminar_flow = available_drop / linear
        else:
            laminar_flow = (
                -linear + math.sqrt(linear * linear + 4.0 * quadratic * available_drop)
            ) / (2.0 * quadratic)
        laminar_reynolds = (
            laminar_flow
            * self.duct.hydraulic_diameter
            / (self.duct.flow_area * self.duct.dynamic_viscosity)
        )
        if laminar_reynolds < 2300.0:
            mass_flow = min(laminar_flow, sonic_cap.mass_flow_rate)
            return FlowResult(
                mass_flow,
                mass_flow == sonic_cap.mass_flow_rate and sonic_cap.is_choked,
            )

        def predicted_drop(mass_flow: float) -> float:
            reynolds = (
                mass_flow
                * self.duct.hydraulic_diameter
                / (self.duct.flow_area * self.duct.dynamic_viscosity)
            )
            friction = _rectangular_darcy_factor(
                reynolds, self.duct.aspect_ratio, self.duct.relative_roughness
            )
            loss = (
                friction * self.duct.length / self.duct.hydraulic_diameter
                + self.duct.minor_loss_coefficient
            )
            velocity = mass_flow / (density * self.duct.flow_area)
            return loss * 0.5 * density * velocity * velocity

        upper = sonic_cap.mass_flow_rate
        if predicted_drop(upper) <= available_drop:
            return FlowResult(upper, True)
        lower = 0.0
        for _ in range(self.iteration_count):
            middle = 0.5 * (lower + upper)
            if predicted_drop(middle) < available_drop:
                lower = middle
            else:
                upper = middle
        return FlowResult(0.5 * (lower + upper), False)

    def bidirectional_flow(
        self,
        first_pressure: float,
        second_pressure: float,
        first_temperature: float,
        second_temperature: float,
        gas: CaloricallyPerfectGas,
    ) -> FlowResult:
        if first_pressure >= second_pressure:
            return self.directed_flow(
                first_pressure, second_pressure, first_temperature, gas
            )
        reverse = self.directed_flow(
            second_pressure, first_pressure, second_temperature, gas
        )
        return FlowResult(-reverse.mass_flow_rate, reverse.is_choked)


@dataclass(frozen=True, slots=True)
class SeriesDuctOrifice:
    """Quasi-steady duct in series with an optional localized orifice."""

    duct: QuasiSteadyCompressibleDuct
    orifice: CompressibleOrifice | None = None
    iteration_count: int = 48

    def directed_flow(
        self,
        upstream_pressure: float,
        downstream_pressure: float,
        upstream_temperature: float,
        gas: CaloricallyPerfectGas,
    ) -> FlowResult:
        if self.orifice is None:
            return self.duct.directed_flow(
                upstream_pressure, downstream_pressure, upstream_temperature, gas
            )
        if downstream_pressure >= upstream_pressure:
            return FlowResult(0.0, False)
        base = self.duct.duct
        orifice_loss_coefficient = (
            base.flow_area / self.orifice.effective_flow_area
        ) ** 2
        equivalent = QuasiSteadyCompressibleDuct(
            RectangularDuct(
                flow_area=base.flow_area,
                hydraulic_diameter=base.hydraulic_diameter,
                length=base.length,
                aspect_ratio=base.aspect_ratio,
                dynamic_viscosity=base.dynamic_viscosity,
                relative_roughness=base.relative_roughness,
                minor_loss_coefficient=(
                    base.minor_loss_coefficient + orifice_loss_coefficient
                ),
            ),
            iteration_count=self.iteration_count,
        )
        duct_result = equivalent.directed_flow(
            upstream_pressure, downstream_pressure, upstream_temperature, gas
        )
        orifice_cap = self.orifice.directed_flow(
            upstream_pressure, downstream_pressure, upstream_temperature, gas
        )
        mass_flow = min(duct_result.mass_flow_rate, orifice_cap.mass_flow_rate)
        return FlowResult(
            mass_flow,
            duct_result.is_choked or mass_flow == orifice_cap.mass_flow_rate,
        )

    def bidirectional_flow(
        self,
        first_pressure: float,
        second_pressure: float,
        first_temperature: float,
        second_temperature: float,
        gas: CaloricallyPerfectGas,
    ) -> FlowResult:
        if first_pressure >= second_pressure:
            return self.directed_flow(
                first_pressure, second_pressure, first_temperature, gas
            )
        reverse = self.directed_flow(
            second_pressure, first_pressure, second_temperature, gas
        )
        return FlowResult(-reverse.mass_flow_rate, reverse.is_choked)


def _rectangular_darcy_factor(
    reynolds: float, aspect_ratio: float, relative_roughness: float
) -> float:
    if reynolds <= 0.0:
        return math.inf
    if reynolds < 2300.0:
        polynomial = (
            1.0
            - 1.3553 * aspect_ratio
            + 1.9467 * aspect_ratio**2
            - 1.7012 * aspect_ratio**3
            + 0.9564 * aspect_ratio**4
            - 0.2537 * aspect_ratio**5
        )
        return 96.0 * polynomial / reynolds
    return (
        -1.8
        * math.log10((relative_roughness / 3.7) ** 1.11 + 6.9 / reynolds)
    ) ** -2


@dataclass(frozen=True, slots=True)
class CompressibleOrifice:
    """First-level isentropic orifice with optional near-zero regularization."""

    effective_flow_area: float
    pressure_regularization: float = 0.0

    def __post_init__(self) -> None:
        if not math.isfinite(self.effective_flow_area) or self.effective_flow_area <= 0.0:
            raise ValueError("Effective flow area CdA must be finite and positive.")
        if (
            not math.isfinite(self.pressure_regularization)
            or self.pressure_regularization < 0.0
        ):
            raise ValueError("Pressure regularization must be finite and non-negative.")

    def critical_pressure_ratio(self, gas: CaloricallyPerfectGas) -> float:
        gamma = gas.heat_capacity_ratio
        return (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))

    def directed_flow(
        self,
        upstream_pressure: float,
        downstream_pressure: float,
        upstream_temperature: float,
        gas: CaloricallyPerfectGas,
    ) -> FlowResult:
        """Return flow from a prescribed upstream state, with reverse flow excluded."""

        for name, value in (
            ("upstream pressure", upstream_pressure),
            ("downstream pressure", downstream_pressure),
            ("upstream temperature", upstream_temperature),
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")

        if downstream_pressure >= upstream_pressure:
            return FlowResult(mass_flow_rate=0.0, is_choked=False)

        gamma = gas.heat_capacity_ratio
        pressure_ratio = downstream_pressure / upstream_pressure
        critical_ratio = self.critical_pressure_ratio(gas)

        if pressure_ratio <= critical_ratio:
            factor = math.sqrt(gamma / (gas.gas_constant * upstream_temperature))
            factor *= (2.0 / (gamma + 1.0)) ** (
                (gamma + 1.0) / (2.0 * (gamma - 1.0))
            )
            return FlowResult(
                mass_flow_rate=self._regularize(
                    self.effective_flow_area * upstream_pressure * factor,
                    upstream_pressure - downstream_pressure,
                ),
                is_choked=True,
            )

        radicand = (
            2.0
            * gamma
            / (gas.gas_constant * upstream_temperature * (gamma - 1.0))
            * (
                pressure_ratio ** (2.0 / gamma)
                - pressure_ratio ** ((gamma + 1.0) / gamma)
            )
        )
        return FlowResult(
            mass_flow_rate=self._regularize(
                self.effective_flow_area
                * upstream_pressure
                * math.sqrt(max(0.0, radicand)),
                upstream_pressure - downstream_pressure,
            ),
            is_choked=False,
        )

    def _regularize(self, mass_flow_rate: float, pressure_difference: float) -> float:
        """Linearize the square-root law only near zero pressure difference."""

        if self.pressure_regularization == 0.0:
            return mass_flow_rate
        return mass_flow_rate * math.sqrt(
            pressure_difference
            / (pressure_difference + self.pressure_regularization)
        )

    def bidirectional_flow(
        self,
        first_pressure: float,
        second_pressure: float,
        first_temperature: float,
        second_temperature: float,
        gas: CaloricallyPerfectGas,
    ) -> FlowResult:
        """Return signed flow, positive from the first node to the second node."""

        if first_pressure >= second_pressure:
            return self.directed_flow(
                first_pressure, second_pressure, first_temperature, gas
            )
        reverse = self.directed_flow(
            second_pressure, first_pressure, second_temperature, gas
        )
        return FlowResult(-reverse.mass_flow_rate, reverse.is_choked)
