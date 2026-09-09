"""Validated physical and numerical configuration structures."""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import tomllib

from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.four_bar import SharedCrankRockerDesign
from dada_solver.geometry import CylinderVolumeLimits, MachineVolumes
from dada_solver.hydraulics import HydraulicNetworkModels
from dada_solver.humidity import HumidityScreeningConfiguration


def _require_positive(name: str, value: float) -> None:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive.")


@dataclass(frozen=True, slots=True)
class HydraulicNetworkConfiguration:
    """Required effective flow areas for every first-level hydraulic link."""

    large_to_hot_cda: float
    small_to_cold_cda: float
    hot_to_small_valve_cda: float
    cold_to_large_valve_cda: float
    orifice_pressure_regularization: float = 0.0

    def __post_init__(self) -> None:
        for name, value in (
            ("large-to-hot CdA", self.large_to_hot_cda),
            ("small-to-cold CdA", self.small_to_cold_cda),
            ("hot-to-small valve CdA", self.hot_to_small_valve_cda),
            ("cold-to-large valve CdA", self.cold_to_large_valve_cda),
        ):
            _require_positive(name, value)
        if (
            not math.isfinite(self.orifice_pressure_regularization)
            or self.orifice_pressure_regularization < 0.0
        ):
            raise ValueError(
                "Orifice pressure regularization must be finite and non-negative."
            )


@dataclass(frozen=True, slots=True)
class ValidityThresholds:
    """User-selected limits required for a physical-validity verdict."""

    maximum_pressure_equalization_error: float
    maximum_mach_number: float
    maximum_isothermality_error: float
    maximum_compressibility_deviation: float
    maximum_cp_variation: float

    def __post_init__(self) -> None:
        for name, value in (
            ("pressure equalization error", self.maximum_pressure_equalization_error),
            ("Mach number", self.maximum_mach_number),
            ("isothermality error", self.maximum_isothermality_error),
            ("compressibility deviation", self.maximum_compressibility_deviation),
            ("Cp variation", self.maximum_cp_variation),
        ):
            _require_positive(f"maximum {name}", value)


@dataclass(frozen=True, slots=True)
class ValveThresholdConfiguration:
    opening_pressure_difference: float
    closing_pressure_difference: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.opening_pressure_difference):
            raise ValueError("Valve opening pressure difference must be finite.")
        if not math.isfinite(self.closing_pressure_difference):
            raise ValueError("Valve closing pressure difference must be finite.")
        if self.closing_pressure_difference > self.opening_pressure_difference:
            raise ValueError("Valve thresholds must satisfy delta_P_close <= delta_P_open.")


@dataclass(frozen=True, slots=True)
class ChargeConfiguration:
    temperature: float
    pressure: float | None = None
    total_mass: float | None = None

    def __post_init__(self) -> None:
        _require_positive("charge temperature", self.temperature)
        if (self.pressure is None) == (self.total_mass is None):
            raise ValueError("Specify exactly one of charge pressure or total mass.")
        if self.pressure is not None:
            _require_positive("charge pressure", self.pressure)
        if self.total_mass is not None:
            _require_positive("total gas mass", self.total_mass)

    def resolved_pressure(
        self,
        gas: CaloricallyPerfectGas,
        total_filling_volume: float,
    ) -> float:
        if self.pressure is not None:
            return self.pressure
        assert self.total_mass is not None
        return self.total_mass * gas.gas_constant * self.temperature / total_filling_volume


@dataclass(frozen=True, slots=True)
class NumericalConfiguration:
    integration_method: str
    integration_relative_tolerance: float
    integration_absolute_tolerance: float
    maximum_step_angle_degrees: float
    maximum_valve_events_per_cycle: int
    maximum_cycles: int
    periodic_relative_tolerance: float
    periodic_mass_absolute_tolerance: float
    periodic_energy_absolute_tolerance: float

    def __post_init__(self) -> None:
        if self.integration_method not in {"RK45", "Radau", "BDF", "LSODA"}:
            raise ValueError("Unsupported integration method.")
        for name, value in (
            ("integration relative tolerance", self.integration_relative_tolerance),
            ("integration absolute tolerance", self.integration_absolute_tolerance),
            ("maximum step angle", self.maximum_step_angle_degrees),
            ("periodic relative tolerance", self.periodic_relative_tolerance),
            ("periodic mass absolute tolerance", self.periodic_mass_absolute_tolerance),
            ("periodic energy absolute tolerance", self.periodic_energy_absolute_tolerance),
        ):
            _require_positive(name, value)
        if self.maximum_valve_events_per_cycle <= 0 or self.maximum_cycles <= 0:
            raise ValueError("Numerical count limits must be positive integers.")


@dataclass(frozen=True, slots=True)
class SimulationConfiguration:
    gas: CaloricallyPerfectGas
    machine_volumes: MachineVolumes
    cold_reservoir_temperature: float
    hot_reservoir_temperature: float
    cold_thermal_conductance: float
    hot_thermal_conductance: float
    angular_speed: float
    charge: ChargeConfiguration
    hydraulics: HydraulicNetworkConfiguration
    hot_to_small_valve: ValveThresholdConfiguration
    cold_to_large_valve: ValveThresholdConfiguration
    validity: ValidityThresholds
    numerical: NumericalConfiguration
    kinematics_type: str
    small_phase_offset_degrees: float | None
    small_lambda_target: float | None
    large_lambda_target: float | None
    adiabatic_sector_fraction: float | None
    example_data: bool
    valve_model: str = "discrete_hysteretic"
    four_bar_ground_distance: float | None = None
    four_bar_connecting_rod_ratio: float | None = None
    four_bar_crank_angle_offset_degrees: float | None = None
    four_bar_crank_direction: int | None = None
    hydraulic_flow_models: HydraulicNetworkModels | None = None
    humidity_screening: HumidityScreeningConfiguration | None = None
    shared_four_bar_design: SharedCrankRockerDesign | None = None

    @property
    def motor_operation(self) -> bool:
        """Negative study crank speed selects reversed motor operation."""

        return self.angular_speed < 0.0

    @property
    def heat_in_reservoir_temperature(self) -> float:
        return (
            self.hot_reservoir_temperature if self.motor_operation
            else self.cold_reservoir_temperature
        )

    @property
    def heat_out_reservoir_temperature(self) -> float:
        return (
            self.cold_reservoir_temperature if self.motor_operation
            else self.hot_reservoir_temperature
        )

    def __post_init__(self) -> None:
        if self.valve_model not in {
            "discrete_hysteretic",
            "continuous_ideal_diode",
        }:
            raise ValueError(f"Unsupported valve model: {self.valve_model}.")
        for name, value in (
            ("cold reservoir temperature", self.cold_reservoir_temperature),
            ("hot reservoir temperature", self.hot_reservoir_temperature),
        ):
            _require_positive(name, value)
        if not math.isfinite(self.angular_speed) or self.angular_speed == 0.0:
            raise ValueError("Angular speed must be finite and non-zero.")
        for name, value in (
            ("cold thermal conductance", self.cold_thermal_conductance),
            ("hot thermal conductance", self.hot_thermal_conductance),
        ):
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and non-negative.")
        if self.kinematics_type == "harmonic_example":
            if self.small_phase_offset_degrees is None or not math.isfinite(self.small_phase_offset_degrees):
                raise ValueError("Small-cylinder phase offset must be finite.")
            if not self.example_data:
                raise ValueError("harmonic_example must be explicitly marked as example_data.")
        elif self.kinematics_type == "ideal_piecewise_linear":
            for name, value in (("small Lambda target", self.small_lambda_target), ("large Lambda target", self.large_lambda_target)):
                if value is None or not math.isfinite(value) or not 0.0 < value < 1.0:
                    raise ValueError(f"{name} must lie strictly between zero and one.")
            if self.adiabatic_sector_fraction is None or not 0.0 < self.adiabatic_sector_fraction < 0.5:
                raise ValueError("Adiabatic sector fraction must lie strictly between zero and 0.5.")
        elif self.kinematics_type in {
            "published_f65_opposed",
            "published_e0_opposed",
        }:
            for name, value in (
                ("four-bar ground distance", self.four_bar_ground_distance),
                ("four-bar connecting-rod ratio", self.four_bar_connecting_rod_ratio),
            ):
                if value is None or not math.isfinite(value) or value <= 0.0:
                    raise ValueError(f"{name} must be finite and positive.")
            if (
                self.four_bar_crank_angle_offset_degrees is None
                or not math.isfinite(self.four_bar_crank_angle_offset_degrees)
            ):
                raise ValueError("Four-bar crank angle offset must be finite.")
            if self.four_bar_crank_direction not in (-1, 1):
                raise ValueError("Four-bar crank direction must be -1 or 1.")
        elif self.kinematics_type == "shared_crank_rocker":
            if self.shared_four_bar_design is None:
                raise ValueError("Shared-crank rocker design is required.")
            if self.four_bar_ground_distance is None or self.four_bar_ground_distance <= 0.0:
                raise ValueError("Four-bar ground distance must be finite and positive.")
            if self.four_bar_crank_angle_offset_degrees is None:
                raise ValueError("Four-bar crank angle offset must be finite.")
            if self.four_bar_crank_direction not in (-1, 1):
                raise ValueError("Four-bar crank direction must be -1 or 1.")
        else:
            raise ValueError(f"Unsupported kinematics type: {self.kinematics_type}.")

def load_simulation_configuration(path: str | Path) -> SimulationConfiguration:
    """Load a complete SI-unit simulation configuration from TOML."""

    config_path = Path(path)
    with config_path.open("rb") as stream:
        data = tomllib.load(stream)
    try:
        gas_data = data["working_gas"]
        geometry = data["geometry"]
        reservoirs = data["reservoirs"]
        heat_transfer = data["heat_transfer"]
        operation = data["operation"]
        charge_data = data["charge"]
        hydraulic_data = data["hydraulics"]
        valve_data = data["valves"]
        validity_data = data["validity"]
        numerical_data = data["numerical"]
        kinematics_data = data["kinematics"]
        humidity_data = data.get("humidity")

        gas = CaloricallyPerfectGas(
            gas_constant=float(gas_data["gas_constant"]),
            heat_capacity_cp=float(gas_data["heat_capacity_cp"]),
            heat_capacity_cv=float(gas_data["heat_capacity_cv"]),
        )
        machine_volumes = MachineVolumes(
            small_cylinder=CylinderVolumeLimits(
                float(geometry["small_cylinder_minimum_volume"]),
                float(geometry["small_cylinder_maximum_volume"]),
            ),
            large_cylinder=CylinderVolumeLimits(
                float(geometry["large_cylinder_minimum_volume"]),
                float(geometry["large_cylinder_maximum_volume"]),
            ),
            cold_heat_exchanger=float(geometry["cold_heat_exchanger_volume"]),
            hot_heat_exchanger=float(geometry["hot_heat_exchanger_volume"]),
        )
        charge = ChargeConfiguration(
            temperature=float(charge_data["temperature"]),
            pressure=_optional_float(charge_data, "pressure"),
            total_mass=_optional_float(charge_data, "total_mass"),
        )
        hydraulics = HydraulicNetworkConfiguration(
            large_to_hot_cda=float(hydraulic_data["large_to_hot_cda"]),
            small_to_cold_cda=float(hydraulic_data["small_to_cold_cda"]),
            hot_to_small_valve_cda=float(hydraulic_data["hot_to_small_valve_cda"]),
            cold_to_large_valve_cda=float(hydraulic_data["cold_to_large_valve_cda"]),
            orifice_pressure_regularization=float(
                hydraulic_data.get("orifice_pressure_regularization", 0.0)
            ),
        )
        numerical = NumericalConfiguration(
            integration_method=str(numerical_data.get("integration_method", "RK45")),
            integration_relative_tolerance=float(
                numerical_data["integration_relative_tolerance"]
            ),
            integration_absolute_tolerance=float(
                numerical_data["integration_absolute_tolerance"]
            ),
            maximum_step_angle_degrees=float(
                numerical_data["maximum_step_angle_degrees"]
            ),
            maximum_valve_events_per_cycle=int(
                numerical_data["maximum_valve_events_per_cycle"]
            ),
            maximum_cycles=int(numerical_data["maximum_cycles"]),
            periodic_relative_tolerance=float(
                numerical_data["periodic_relative_tolerance"]
            ),
            periodic_mass_absolute_tolerance=float(
                numerical_data["periodic_mass_absolute_tolerance"]
            ),
            periodic_energy_absolute_tolerance=float(
                numerical_data["periodic_energy_absolute_tolerance"]
            ),
        )
        return SimulationConfiguration(
            gas=gas,
            machine_volumes=machine_volumes,
            cold_reservoir_temperature=float(reservoirs["cold_temperature"]),
            hot_reservoir_temperature=float(reservoirs["hot_temperature"]),
            cold_thermal_conductance=float(heat_transfer["cold_ua"]),
            hot_thermal_conductance=float(heat_transfer["hot_ua"]),
            angular_speed=float(operation["angular_speed"]),
            charge=charge,
            hydraulics=hydraulics,
            hot_to_small_valve=_load_valve_thresholds(valve_data["hot_to_small"]),
            cold_to_large_valve=_load_valve_thresholds(valve_data["cold_to_large"]),
            validity=ValidityThresholds(
                maximum_pressure_equalization_error=float(
                    validity_data["maximum_pressure_equalization_error"]
                ),
                maximum_mach_number=float(validity_data["maximum_mach_number"]),
                maximum_isothermality_error=float(
                    validity_data["maximum_isothermality_error"]
                ),
                maximum_compressibility_deviation=float(
                    validity_data["maximum_compressibility_deviation"]
                ),
                maximum_cp_variation=float(validity_data["maximum_cp_variation"]),
            ),
            numerical=numerical,
            kinematics_type=str(kinematics_data["type"]),
            small_phase_offset_degrees=_optional_float(kinematics_data, "small_phase_offset_degrees"),
            small_lambda_target=_optional_float(kinematics_data, "small_lambda_target"),
            large_lambda_target=_optional_float(kinematics_data, "large_lambda_target"),
            adiabatic_sector_fraction=_optional_float(kinematics_data, "adiabatic_sector_fraction"),
            example_data=bool(data["metadata"]["example_data"]),
            valve_model=str(valve_data.get("model", "discrete_hysteretic")),
            four_bar_ground_distance=_optional_float(kinematics_data, "ground_distance"),
            four_bar_connecting_rod_ratio=_optional_float(
                kinematics_data, "connecting_rod_to_projected_stroke_ratio"
            ),
            four_bar_crank_angle_offset_degrees=_optional_float(
                kinematics_data, "crank_angle_offset_degrees"
            ),
            four_bar_crank_direction=(
                None
                if "crank_direction" not in kinematics_data
                else int(kinematics_data["crank_direction"])
            ),
            humidity_screening=(
                None
                if humidity_data is None
                else HumidityScreeningConfiguration(
                    relative_humidity=float(humidity_data["relative_humidity"])
                )
            ),
            shared_four_bar_design=_load_shared_four_bar_design(kinematics_data),
        )
    except KeyError as error:
        raise ValueError(f"Missing required configuration field: {error.args[0]}") from error
    except (TypeError, ValueError) as error:
        if isinstance(error, ValueError) and not str(error).startswith("could not convert"):
            raise
        raise ValueError(f"Invalid configuration value: {error}") from error


def _optional_float(data: dict[str, object], key: str) -> float | None:
    value = data.get(key)
    return None if value is None else float(value)


def _load_valve_thresholds(data: dict[str, object]) -> ValveThresholdConfiguration:
    return ValveThresholdConfiguration(
        opening_pressure_difference=float(data["opening_pressure_difference"]),
        closing_pressure_difference=float(data["closing_pressure_difference"]),
    )


def _load_shared_four_bar_design(
    data: dict[str, object],
) -> SharedCrankRockerDesign | None:
    if data.get("type") != "shared_crank_rocker":
        return None
    small = data["small_four_bar"]
    large = data["large_four_bar"]
    if not isinstance(small, dict) or not isinstance(large, dict):
        raise ValueError("Four-bar loop configurations must be tables.")
    return SharedCrankRockerDesign(
        crank_ratio=float(data["crank_ratio"]),
        small_coupler_ratio=float(small["coupler_ratio"]),
        small_rocker_ratio=float(small["rocker_ratio"]),
        small_output_along_ratio=float(small["output_along_ratio"]),
        small_output_normal_ratio=float(small["output_normal_ratio"]),
        large_coupler_ratio=float(large["coupler_ratio"]),
        large_rocker_ratio=float(large["rocker_ratio"]),
        large_output_along_ratio=float(large["output_along_ratio"]),
        large_output_normal_ratio=float(large["output_normal_ratio"]),
        small_pivot_x_ratio=float(small.get("pivot_x_ratio", 1.0)),
        small_pivot_y_ratio=float(small.get("pivot_y_ratio", 0.0)),
        large_pivot_x_ratio=float(large.get("pivot_x_ratio", 1.0)),
        large_pivot_y_ratio=float(large.get("pivot_y_ratio", 0.0)),
        small_four_bar_branch=int(small.get("assembly_branch", 1)),
        large_four_bar_branch=int(large.get("assembly_branch", -1)),
        small_slider_rod_ratio=float(small.get("slider_rod_ratio", 5.0)),
        large_slider_rod_ratio=float(large.get("slider_rod_ratio", 5.0)),
    )
