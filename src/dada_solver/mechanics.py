"""Thermodynamic load boundary for future mechanical models."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Protocol

import numpy as np
from numpy.typing import NDArray

from dada_solver.dynamics import ThermodynamicModel
from dada_solver.integration import CycleIntegrationResult
from dada_solver.state import ThermodynamicState


@dataclass(frozen=True, slots=True)
class PistonFaceAreas:
    """Gas-loaded piston face areas in m^2."""

    small: float
    large: float

    def __post_init__(self) -> None:
        for name, value in (("small", self.small), ("large", self.large)):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} piston face area must be finite and positive.")


@dataclass(frozen=True, slots=True)
class PistonExternalPressures:
    """Absolute pressures acting on the non-gas side of each piston in Pa."""

    small: float
    large: float

    def __post_init__(self) -> None:
        for name, value in (("small", self.small), ("large", self.large)):
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(
                    f"{name} external pressure must be finite and non-negative."
                )


@dataclass(frozen=True, slots=True)
class ThermodynamicLoadHistory:
    """Gas-side loads only; no inertia, friction or structural reactions."""

    angles: NDArray[np.float64]
    small_pressure: NDArray[np.float64]
    large_pressure: NDArray[np.float64]
    small_gas_force: NDArray[np.float64]
    large_gas_force: NDArray[np.float64]
    generalized_gas_torque: NDArray[np.float64]
    generalized_net_pressure_torque: NDArray[np.float64] | None
    small_net_pressure_force: NDArray[np.float64] | None
    large_net_pressure_force: NDArray[np.float64] | None


class MechanicalResponseModel(Protocol):
    """Future replaceable linkage, inertia, friction and load model."""

    def evaluate(self, loads: ThermodynamicLoadHistory) -> object:
        """Return a model-specific mechanical response."""


def extract_thermodynamic_loads(
    cycle: CycleIntegrationResult,
    model: ThermodynamicModel,
    piston_areas: PistonFaceAreas,
    external_pressures: PistonExternalPressures | None = None,
) -> ThermodynamicLoadHistory:
    """Extract absolute gas loads and optional net pressure forces."""

    if not cycle.completed:
        raise ValueError("Mechanical load extraction requires a completed cycle.")
    sample_count = cycle.angles.size
    small_pressure = np.empty(sample_count)
    large_pressure = np.empty(sample_count)
    generalized_torque = extract_generalized_gas_torque(cycle, model)
    for index, (angle, values) in enumerate(
        zip(cycle.angles, cycle.states.T, strict=True)
    ):
        state = ThermodynamicState.from_array(values)
        pressures = state.pressures(model.gas, model.volumes(float(angle)))
        small_pressure[index] = pressures[0]
        large_pressure[index] = pressures[1]
    return ThermodynamicLoadHistory(
        angles=cycle.angles.copy(),
        small_pressure=small_pressure,
        large_pressure=large_pressure,
        small_gas_force=small_pressure * piston_areas.small,
        large_gas_force=large_pressure * piston_areas.large,
        generalized_gas_torque=generalized_torque,
        generalized_net_pressure_torque=(
            None
            if external_pressures is None
            else generalized_torque
            - np.array(
                [
                    external_pressures.small
                    * model.kinematics.small_cylinder_volume_derivative(float(angle))
                    + external_pressures.large
                    * model.kinematics.large_cylinder_volume_derivative(float(angle))
                    for angle in cycle.angles
                ]
            )
        ),
        small_net_pressure_force=(
            None
            if external_pressures is None
            else (small_pressure - external_pressures.small) * piston_areas.small
        ),
        large_net_pressure_force=(
            None
            if external_pressures is None
            else (large_pressure - external_pressures.large) * piston_areas.large
        ),
    )


def extract_generalized_gas_torque(
    cycle: CycleIntegrationResult,
    model: ThermodynamicModel,
) -> NDArray[np.float64]:
    """Return sum(P_i*dV_i/dtheta), without a linkage or shaft-load model."""

    if not cycle.completed:
        raise ValueError("Gas torque extraction requires a completed cycle.")
    torque = np.empty(cycle.angles.size)
    for index, (angle, values) in enumerate(
        zip(cycle.angles, cycle.states.T, strict=True)
    ):
        state = ThermodynamicState.from_array(values)
        pressures = state.pressures(model.gas, model.volumes(float(angle)))
        torque[index] = (
            pressures[0]
            * model.kinematics.small_cylinder_volume_derivative(float(angle))
            + pressures[1]
            * model.kinematics.large_cylinder_volume_derivative(float(angle))
        )
    return torque
