"""TOML input for finite-grid heat-exchanger screening studies."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tomllib

from dada_solver.exchangers.models import (
    ExchangerOperatingPoint,
    ThermalResistanceModel,
    TransportProperties,
)
from dada_solver.exchangers.sizing import ExchangerRequirements
from dada_solver.fluids import CaloricallyPerfectGas


@dataclass(frozen=True, slots=True)
class ExchangerScreeningConfiguration:
    gas: CaloricallyPerfectGas
    transport: TransportProperties
    operating_point: ExchangerOperatingPoint
    thermal_resistances: ThermalResistanceModel
    requirements: ExchangerRequirements
    channel_counts: tuple[int, ...]
    channel_widths: tuple[float, ...]
    channel_heights: tuple[float, ...]
    channel_lengths: tuple[float, ...]
    surface_roughness: float
    minor_loss_coefficient: float
    example_data: bool


def load_exchanger_screening_configuration(
    path: str | Path,
) -> ExchangerScreeningConfiguration:
    """Load an explicitly bounded exchanger study in SI units."""

    with Path(path).open("rb") as stream:
        data = tomllib.load(stream)
    try:
        gas_data = data["working_gas"]
        transport = data["transport"]
        operating = data["operating_point"]
        thermal = data["thermal_resistances"]
        requirements = data["requirements"]
        grid = data["parallel_channel_grid"]
        return ExchangerScreeningConfiguration(
            gas=CaloricallyPerfectGas(
                gas_constant=float(gas_data["gas_constant"]),
                heat_capacity_cp=float(gas_data["heat_capacity_cp"]),
                heat_capacity_cv=float(gas_data["heat_capacity_cv"]),
            ),
            transport=TransportProperties(
                dynamic_viscosity=float(transport["dynamic_viscosity"]),
                thermal_conductivity=float(transport["thermal_conductivity"]),
            ),
            operating_point=ExchangerOperatingPoint(
                absolute_mass_flow_rate=float(operating["absolute_mass_flow_rate"]),
                pressure=float(operating["pressure"]),
                temperature=float(operating["temperature"]),
            ),
            thermal_resistances=ThermalResistanceModel(
                wall_thickness=float(thermal["wall_thickness"]),
                wall_thermal_conductivity=float(
                    thermal["wall_thermal_conductivity"]
                ),
                external_conductance=float(thermal["external_conductance"]),
                additional_resistance=float(
                    thermal.get("additional_resistance", 0.0)
                ),
            ),
            requirements=ExchangerRequirements(
                minimum_overall_conductance=float(
                    requirements["minimum_overall_conductance"]
                ),
                maximum_pressure_drop=float(requirements["maximum_pressure_drop"]),
                maximum_mach_number=float(requirements["maximum_mach_number"]),
                maximum_gas_volume=float(requirements["maximum_gas_volume"]),
                maximum_pressure_drop_fraction=float(
                    requirements["maximum_pressure_drop_fraction"]
                ),
                minimum_effectiveness=_optional_float(
                    requirements, "minimum_effectiveness"
                ),
            ),
            channel_counts=tuple(int(value) for value in grid["channel_counts"]),
            channel_widths=tuple(float(value) for value in grid["channel_widths"]),
            channel_heights=tuple(float(value) for value in grid["channel_heights"]),
            channel_lengths=tuple(float(value) for value in grid["channel_lengths"]),
            surface_roughness=float(grid.get("surface_roughness", 0.0)),
            minor_loss_coefficient=float(grid.get("minor_loss_coefficient", 0.0)),
            example_data=bool(data["metadata"]["example_data"]),
        )
    except KeyError as error:
        raise ValueError(
            f"Missing exchanger configuration field: {error.args[0]}"
        ) from error


def _optional_float(data: dict[str, object], key: str) -> float | None:
    value = data.get(key)
    return None if value is None else float(value)
