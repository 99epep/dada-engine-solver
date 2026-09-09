"""Bounded design variables and explicit configuration transformations."""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
import math
from typing import Mapping

from dada_solver.configuration import (
    ChargeConfiguration,
    HydraulicNetworkConfiguration,
    SimulationConfiguration,
)
from dada_solver.geometry import CylinderVolumeLimits, MachineVolumes


class DesignParameter(Enum):
    SMALL_SWEPT_VOLUME = "small_swept_volume"
    LARGE_SWEPT_VOLUME = "large_swept_volume"
    SMALL_CLEARANCE_VOLUME = "small_clearance_volume"
    LARGE_CLEARANCE_VOLUME = "large_clearance_volume"
    SMALL_CLEARANCE_RATIO = "small_clearance_ratio"
    LARGE_CLEARANCE_RATIO = "large_clearance_ratio"
    COLD_HEAT_EXCHANGER_VOLUME = "cold_heat_exchanger_volume"
    HOT_HEAT_EXCHANGER_VOLUME = "hot_heat_exchanger_volume"
    COLD_UA = "cold_ua"
    HOT_UA = "hot_ua"
    COLD_RESERVOIR_TEMPERATURE = "cold_reservoir_temperature"
    HOT_RESERVOIR_TEMPERATURE = "hot_reservoir_temperature"
    LARGE_TO_HOT_CDA = "large_to_hot_cda"
    SMALL_TO_COLD_CDA = "small_to_cold_cda"
    HOT_TO_SMALL_VALVE_CDA = "hot_to_small_valve_cda"
    COLD_TO_LARGE_VALVE_CDA = "cold_to_large_valve_cda"
    CHARGE_PRESSURE = "charge_pressure"
    TOTAL_GAS_MASS = "total_gas_mass"
    ANGULAR_SPEED = "angular_speed"
    SMALL_PHASE_OFFSET = "small_phase_offset"
    COMMON_LAMBDA_TARGET = "common_lambda_target"
    SMALL_LAMBDA_TARGET = "small_lambda_target"
    LARGE_LAMBDA_TARGET = "large_lambda_target"
    ADIABATIC_SECTOR_FRACTION = "adiabatic_sector_fraction"


@dataclass(frozen=True, slots=True)
class DesignVariable:
    parameter: DesignParameter
    lower_bound: float
    upper_bound: float
    initial_value: float

    def __post_init__(self) -> None:
        if self.parameter is DesignParameter.ANGULAR_SPEED:
            if not all(math.isfinite(value) for value in (
                self.lower_bound, self.upper_bound, self.initial_value
            )):
                raise ValueError("Angular-speed bounds and initial value must be finite.")
            if self.lower_bound <= 0.0 <= self.upper_bound:
                raise ValueError("Angular-speed bounds must not include zero or cross operating modes.")
            if self.upper_bound <= self.lower_bound:
                raise ValueError("Design-variable upper bound must exceed lower bound.")
            if not self.lower_bound <= self.initial_value <= self.upper_bound:
                raise ValueError("Initial design value must lie within its bounds.")
            return
        allow_zero = self.parameter in {
            DesignParameter.COLD_UA,
            DesignParameter.HOT_UA,
            DesignParameter.SMALL_PHASE_OFFSET,
        }
        if not math.isfinite(self.lower_bound) or (
            self.lower_bound < 0.0 if allow_zero else self.lower_bound <= 0.0
        ):
            qualifier = "non-negative" if allow_zero else "positive"
            raise ValueError(f"Design-variable lower bound must be finite and {qualifier}.")
        for name, value in (
            ("upper bound", self.upper_bound),
            ("initial value", self.initial_value),
        ):
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"Design-variable {name} must be finite and non-negative.")
        if self.upper_bound <= self.lower_bound:
            raise ValueError("Design-variable upper bound must exceed lower bound.")
        if not self.lower_bound <= self.initial_value <= self.upper_bound:
            raise ValueError("Initial design value must lie within its bounds.")


@dataclass(frozen=True, slots=True)
class DesignPoint:
    values: Mapping[DesignParameter, float]

    def __post_init__(self) -> None:
        for parameter, value in self.values.items():
            if not isinstance(parameter, DesignParameter):
                raise TypeError("Design-point keys must be DesignParameter values.")
            if parameter is DesignParameter.ANGULAR_SPEED:
                if not math.isfinite(value) or value == 0.0:
                    raise ValueError("Angular speed must be finite and non-zero.")
                continue
            allow_zero = parameter in {
                DesignParameter.COLD_UA,
                DesignParameter.HOT_UA,
                DesignParameter.SMALL_PHASE_OFFSET,
            }
            invalid = value < 0.0 if allow_zero else value <= 0.0
            if not math.isfinite(value) or invalid:
                qualifier = "non-negative" if allow_zero else "positive"
                raise ValueError(
                    f"Design value {parameter.value} must be finite and {qualifier}."
                )

    def value_or(self, parameter: DesignParameter, fallback: float) -> float:
        return float(self.values.get(parameter, fallback))


def apply_design_point(
    base: SimulationConfiguration,
    point: DesignPoint,
) -> SimulationConfiguration:
    """Return a new configuration with only explicitly supported changes."""

    harmonic_parameters = {DesignParameter.SMALL_PHASE_OFFSET}
    piecewise_parameters = {
        DesignParameter.COMMON_LAMBDA_TARGET,
        DesignParameter.SMALL_LAMBDA_TARGET,
        DesignParameter.LARGE_LAMBDA_TARGET,
        DesignParameter.ADIABATIC_SECTOR_FRACTION,
    }
    if base.kinematics_type != "harmonic_example" and harmonic_parameters & point.values.keys():
        raise ValueError("Small phase offset applies only to harmonic_example kinematics.")
    if base.kinematics_type != "ideal_piecewise_linear" and piecewise_parameters & point.values.keys():
        raise ValueError("Lambda and sector-duration parameters apply only to ideal_piecewise_linear kinematics.")
    common_lambda = point.values.get(DesignParameter.COMMON_LAMBDA_TARGET)
    if common_lambda is not None and (
        DesignParameter.SMALL_LAMBDA_TARGET in point.values
        or DesignParameter.LARGE_LAMBDA_TARGET in point.values
    ):
        raise ValueError(
            "Common Lambda target cannot be combined with an individual Lambda target."
        )

    base_volumes = base.machine_volumes
    small_swept = point.value_or(
        DesignParameter.SMALL_SWEPT_VOLUME,
        base_volumes.small_cylinder.swept,
    )
    large_swept = point.value_or(
        DesignParameter.LARGE_SWEPT_VOLUME,
        base_volumes.large_cylinder.swept,
    )
    small_clearance = _resolve_clearance(
        point,
        DesignParameter.SMALL_CLEARANCE_VOLUME,
        DesignParameter.SMALL_CLEARANCE_RATIO,
        base_volumes.small_cylinder.minimum,
        small_swept,
    )
    large_clearance = _resolve_clearance(
        point,
        DesignParameter.LARGE_CLEARANCE_VOLUME,
        DesignParameter.LARGE_CLEARANCE_RATIO,
        base_volumes.large_cylinder.minimum,
        large_swept,
    )
    volumes = MachineVolumes(
        small_cylinder=CylinderVolumeLimits(
            small_clearance, small_clearance + small_swept
        ),
        large_cylinder=CylinderVolumeLimits(
            large_clearance, large_clearance + large_swept
        ),
        cold_heat_exchanger=point.value_or(
            DesignParameter.COLD_HEAT_EXCHANGER_VOLUME,
            base_volumes.cold_heat_exchanger,
        ),
        hot_heat_exchanger=point.value_or(
            DesignParameter.HOT_HEAT_EXCHANGER_VOLUME,
            base_volumes.hot_heat_exchanger,
        ),
    )
    old_hydraulics = base.hydraulics
    hydraulics = HydraulicNetworkConfiguration(
        large_to_hot_cda=point.value_or(
            DesignParameter.LARGE_TO_HOT_CDA, old_hydraulics.large_to_hot_cda
        ),
        small_to_cold_cda=point.value_or(
            DesignParameter.SMALL_TO_COLD_CDA, old_hydraulics.small_to_cold_cda
        ),
        hot_to_small_valve_cda=point.value_or(
            DesignParameter.HOT_TO_SMALL_VALVE_CDA,
            old_hydraulics.hot_to_small_valve_cda,
        ),
        cold_to_large_valve_cda=point.value_or(
            DesignParameter.COLD_TO_LARGE_VALVE_CDA,
            old_hydraulics.cold_to_large_valve_cda,
        ),
    )
    charge = _updated_charge(base.charge, point)
    return replace(
        base,
        machine_volumes=volumes,
        cold_thermal_conductance=point.value_or(
            DesignParameter.COLD_UA, base.cold_thermal_conductance
        ),
        hot_thermal_conductance=point.value_or(
            DesignParameter.HOT_UA, base.hot_thermal_conductance
        ),
        cold_reservoir_temperature=point.value_or(
            DesignParameter.COLD_RESERVOIR_TEMPERATURE,
            base.cold_reservoir_temperature,
        ),
        hot_reservoir_temperature=point.value_or(
            DesignParameter.HOT_RESERVOIR_TEMPERATURE,
            base.hot_reservoir_temperature,
        ),
        hydraulics=hydraulics,
        charge=charge,
        angular_speed=point.value_or(DesignParameter.ANGULAR_SPEED, base.angular_speed),
        small_phase_offset_degrees=math.degrees(
            point.value_or(
                DesignParameter.SMALL_PHASE_OFFSET,
                math.radians(base.small_phase_offset_degrees or 0.0),
            )
        ),
        small_lambda_target=point.value_or(
            DesignParameter.SMALL_LAMBDA_TARGET,
            float(common_lambda) if common_lambda is not None else base.small_lambda_target or 0.7,
        ),
        large_lambda_target=point.value_or(
            DesignParameter.LARGE_LAMBDA_TARGET,
            float(common_lambda) if common_lambda is not None else base.large_lambda_target or 0.7,
        ),
        adiabatic_sector_fraction=point.value_or(
            DesignParameter.ADIABATIC_SECTOR_FRACTION,
            base.adiabatic_sector_fraction or 0.25,
        ),
    )


def _updated_charge(
    base: ChargeConfiguration,
    point: DesignPoint,
) -> ChargeConfiguration:
    pressure = point.values.get(DesignParameter.CHARGE_PRESSURE)
    total_mass = point.values.get(DesignParameter.TOTAL_GAS_MASS)
    if pressure is not None and total_mass is not None:
        raise ValueError("Charge pressure and total gas mass cannot both be designed.")
    if pressure is not None:
        return ChargeConfiguration(temperature=base.temperature, pressure=float(pressure))
    if total_mass is not None:
        return ChargeConfiguration(
            temperature=base.temperature, total_mass=float(total_mass)
        )
    return base


def _resolve_clearance(
    point: DesignPoint,
    absolute_parameter: DesignParameter,
    ratio_parameter: DesignParameter,
    base_clearance: float,
    swept_volume: float,
) -> float:
    absolute = point.values.get(absolute_parameter)
    ratio = point.values.get(ratio_parameter)
    if absolute is not None and ratio is not None:
        raise ValueError(
            f"{absolute_parameter.value} and {ratio_parameter.value} cannot both be designed."
        )
    if ratio is not None:
        return float(ratio) * swept_volume
    if absolute is not None:
        return float(absolute)
    return base_clearance
