"""Diagnostics for deciding whether quasi-steady duct hydraulics are adequate."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from numpy.typing import ArrayLike

from dada_solver.exchangers.models import ParallelRectangularChannels
from dada_solver.fluids import CaloricallyPerfectGas


@dataclass(frozen=True, slots=True)
class HydraulicValidityThresholds:
    maximum_mach_number: float
    maximum_pressure_drop_fraction: float
    maximum_acoustic_time_ratio: float
    maximum_inertial_pressure_fraction: float

    def __post_init__(self) -> None:
        for name, value in (
            ("Maximum Mach number", self.maximum_mach_number),
            ("Maximum pressure-drop fraction", self.maximum_pressure_drop_fraction),
            ("Maximum acoustic-time ratio", self.maximum_acoustic_time_ratio),
            ("Maximum inertial-pressure fraction", self.maximum_inertial_pressure_fraction),
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")


@dataclass(frozen=True, slots=True)
class HydraulicQuasiSteadyValidity:
    maximum_mach_number: float
    maximum_pressure_drop_fraction: float
    acoustic_time_ratio: float
    maximum_inertial_pressure_fraction: float
    quasi_steady_valid: bool
    inertia_model_recommended: bool
    failed_criteria: tuple[str, ...]


def assess_hydraulic_quasi_steady_validity(
    geometry: ParallelRectangularChannels,
    mass_flow_history: ArrayLike,
    period: float,
    reference_pressure: float,
    reference_temperature: float,
    maximum_pressure_drop: float,
    gas: CaloricallyPerfectGas,
    thresholds: HydraulicValidityThresholds,
) -> HydraulicQuasiSteadyValidity:
    """Estimate acoustic and lumped inertial terms from a periodic flow history."""

    for name, value in (
        ("Period", period),
        ("Reference pressure", reference_pressure),
        ("Reference temperature", reference_temperature),
    ):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive.")
    if not math.isfinite(maximum_pressure_drop) or maximum_pressure_drop < 0.0:
        raise ValueError("Maximum pressure drop must be finite and non-negative.")
    flow = np.asarray(mass_flow_history, dtype=float)
    if flow.ndim != 1 or flow.size < 3 or not np.all(np.isfinite(flow)):
        raise ValueError("Mass-flow history must contain at least three finite samples.")

    density = reference_pressure / (gas.gas_constant * reference_temperature)
    sound_speed = math.sqrt(
        gas.heat_capacity_ratio * gas.gas_constant * reference_temperature
    )
    velocity = np.abs(flow) / (density * geometry.total_flow_area)
    maximum_mach = float(np.max(velocity) / sound_speed)
    pressure_fraction = maximum_pressure_drop / reference_pressure
    acoustic_ratio = geometry.channel_length / sound_speed / period
    sample_interval = period / flow.size
    derivative = np.gradient(flow, sample_interval, edge_order=2)
    inertial_pressure = geometry.channel_length / geometry.total_flow_area * derivative
    inertial_fraction = float(np.max(np.abs(inertial_pressure)) / reference_pressure)

    failed = []
    if maximum_mach > thresholds.maximum_mach_number:
        failed.append("maximum_mach_number")
    if pressure_fraction > thresholds.maximum_pressure_drop_fraction:
        failed.append("maximum_pressure_drop_fraction")
    if acoustic_ratio > thresholds.maximum_acoustic_time_ratio:
        failed.append("maximum_acoustic_time_ratio")
    if inertial_fraction > thresholds.maximum_inertial_pressure_fraction:
        failed.append("maximum_inertial_pressure_fraction")
    inertia_failed = bool(
        {"maximum_acoustic_time_ratio", "maximum_inertial_pressure_fraction"}
        & set(failed)
    )
    return HydraulicQuasiSteadyValidity(
        maximum_mach,
        pressure_fraction,
        acoustic_ratio,
        inertial_fraction,
        not failed,
        inertia_failed,
        tuple(failed),
    )
