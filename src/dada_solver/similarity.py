"""Exact extensive scaling transformation for the first-level 0D model."""

from __future__ import annotations

from dataclasses import replace
import math

from dada_solver.configuration import (
    ChargeConfiguration,
    HydraulicNetworkConfiguration,
    SimulationConfiguration,
    ValveThresholdConfiguration,
)
from dada_solver.geometry import CylinderVolumeLimits, MachineVolumes
from dada_solver.performance import CyclePerformance


def required_power_similarity_factor(
    reference: CyclePerformance, target_cooling_power: float
) -> float:
    """Return the exact first-level similarity multiplier for a cooling target."""

    if not math.isfinite(target_cooling_power) or target_cooling_power <= 0.0:
        raise ValueError("Target cooling power must be finite and positive.")
    if reference.cooling_power <= 0.0 or reference.cooling_cop is None:
        raise ValueError("Reference point must operate as a refrigerator.")
    return target_cooling_power / reference.cooling_power


def scale_extensive_machine(
    configuration: SimulationConfiguration,
    factor: float,
) -> SimulationConfiguration:
    """Scale all extensive first-level quantities at unchanged cycle frequency.

    Volumes, gas inventory, UA and effective flow areas use the same factor.
    Dimensionless kinematics and intensive filling conditions are preserved.
    A physical four-bar length scale uses the cube root so piston geometric
    similarity is retained. This mathematical similarity does not establish
    that a real heat exchanger's UA scales linearly with its gas volume.
    """

    return scale_capacity_and_speed(configuration, factor, 1.0)


def scale_capacity_and_speed(
    configuration: SimulationConfiguration,
    capacity_factor: float,
    speed_factor: float,
) -> SimulationConfiguration:
    """Scale machine capacity and cycle frequency as independent factors.

    Pressure and temperature histories remain similar when volume and inventory
    scale with ``capacity_factor``, angular speed scales with ``speed_factor``,
    and both UA and CdA scale with their product. Cooling and input powers scale
    with that same product.
    """

    for name, value in (
        ("Capacity factor", capacity_factor),
        ("Speed factor", speed_factor),
    ):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive.")
    factor = capacity_factor
    rate_factor = capacity_factor * speed_factor
    old = configuration.machine_volumes
    volumes = MachineVolumes(
        small_cylinder=_scale_cylinder(old.small_cylinder, factor),
        large_cylinder=_scale_cylinder(old.large_cylinder, factor),
        cold_heat_exchanger=old.cold_heat_exchanger * factor,
        hot_heat_exchanger=old.hot_heat_exchanger * factor,
    )
    hydraulics = configuration.hydraulics
    scaled_hydraulics = HydraulicNetworkConfiguration(
        large_to_hot_cda=hydraulics.large_to_hot_cda * rate_factor,
        small_to_cold_cda=hydraulics.small_to_cold_cda * rate_factor,
        hot_to_small_valve_cda=hydraulics.hot_to_small_valve_cda * rate_factor,
        cold_to_large_valve_cda=hydraulics.cold_to_large_valve_cda * rate_factor,
        orifice_pressure_regularization=hydraulics.orifice_pressure_regularization,
    )
    charge = configuration.charge
    scaled_charge = (
        ChargeConfiguration(temperature=charge.temperature, pressure=charge.pressure)
        if charge.pressure is not None
        else ChargeConfiguration(
            temperature=charge.temperature,
            total_mass=charge.total_mass * factor,  # type: ignore[operator]
        )
    )
    ground_distance = configuration.four_bar_ground_distance
    return replace(
        configuration,
        machine_volumes=volumes,
        hydraulics=scaled_hydraulics,
        cold_thermal_conductance=configuration.cold_thermal_conductance * rate_factor,
        hot_thermal_conductance=configuration.hot_thermal_conductance * rate_factor,
        angular_speed=configuration.angular_speed * speed_factor,
        charge=scaled_charge,
        four_bar_ground_distance=(
            None if ground_distance is None else ground_distance * factor ** (1.0 / 3.0)
        ),
    )


def scale_volume_at_constant_inventory(
    configuration: SimulationConfiguration,
    volume_factor: float,
) -> SimulationConfiguration:
    """Trade volume against pressure while preserving the ideal 0D cycle.

    Every control volume and CdA is multiplied by ``volume_factor`` and the
    filling pressure and absolute valve thresholds are divided by it. UA,
    angular speed and gas inventory remain unchanged. Consequently `P*V`, mass
    flow, heat rate, work, temperatures and ideal COP are invariant.

    This is an exact similarity only for the first-level ideal-gas/orifice
    model. Higher pressure still requires real-gas, seal, force, exchanger and
    structural checks.
    """

    if not math.isfinite(volume_factor) or volume_factor <= 0.0:
        raise ValueError("Volume factor must be finite and positive.")
    charge = configuration.charge
    if charge.pressure is None:
        raise ValueError("Pressure-volume scaling requires a pressure-based charge.")
    old = configuration.machine_volumes
    hydraulics = configuration.hydraulics
    ground_distance = configuration.four_bar_ground_distance

    def scale_thresholds(
        thresholds: ValveThresholdConfiguration,
    ) -> ValveThresholdConfiguration:
        return ValveThresholdConfiguration(
            opening_pressure_difference=(
                thresholds.opening_pressure_difference / volume_factor
            ),
            closing_pressure_difference=(
                thresholds.closing_pressure_difference / volume_factor
            ),
        )

    return replace(
        configuration,
        machine_volumes=MachineVolumes(
            small_cylinder=_scale_cylinder(old.small_cylinder, volume_factor),
            large_cylinder=_scale_cylinder(old.large_cylinder, volume_factor),
            cold_heat_exchanger=old.cold_heat_exchanger * volume_factor,
            hot_heat_exchanger=old.hot_heat_exchanger * volume_factor,
        ),
        hydraulics=HydraulicNetworkConfiguration(
            large_to_hot_cda=hydraulics.large_to_hot_cda * volume_factor,
            small_to_cold_cda=hydraulics.small_to_cold_cda * volume_factor,
            hot_to_small_valve_cda=(
                hydraulics.hot_to_small_valve_cda * volume_factor
            ),
            cold_to_large_valve_cda=(
                hydraulics.cold_to_large_valve_cda * volume_factor
            ),
            orifice_pressure_regularization=(
                hydraulics.orifice_pressure_regularization / volume_factor
            ),
        ),
        charge=ChargeConfiguration(
            temperature=charge.temperature,
            pressure=charge.pressure / volume_factor,
        ),
        hot_to_small_valve=scale_thresholds(configuration.hot_to_small_valve),
        cold_to_large_valve=scale_thresholds(configuration.cold_to_large_valve),
        four_bar_ground_distance=(
            None
            if ground_distance is None
            else ground_distance * volume_factor ** (1.0 / 3.0)
        ),
    )


def _scale_cylinder(
    cylinder: CylinderVolumeLimits, factor: float
) -> CylinderVolumeLimits:
    return CylinderVolumeLimits(cylinder.minimum * factor, cylinder.maximum * factor)
