"""Thermodynamic evaluation of shared-crank four-bar design candidates."""

from __future__ import annotations

from dataclasses import dataclass, replace
import math

from dada_solver.configuration import ChargeConfiguration, SimulationConfiguration
from dada_solver.dynamics import ThermodynamicModel
from dada_solver.factory import (
    build_initial_state,
    build_model,
    build_periodic_solver,
    initial_valve_topology,
)
from dada_solver.four_bar import (
    SharedCrankFourBarVolumeKinematics,
    SharedCrankRockerDesign,
    shared_crank_rocker_kinematics,
)
from dada_solver.performance import CyclePerformance, calculate_cycle_performance
from dada_solver.periodic import PeriodicResult, PeriodicStatus
from dada_solver.results import CycleDiagnostics, extract_cycle_diagnostics
from dada_solver.state import ThermodynamicState
from dada_solver.validity import ValidityReport, assess_cycle_validity


@dataclass(frozen=True, slots=True)
class ReflectedRockerParameters:
    """Compact reflected design used for the first local E0 search."""

    crank_ratio: float
    coupler_ratio: float
    rocker_ratio: float
    output_angle_degrees: float
    crank_angle_offset_degrees: float
    connecting_rod_ratio: float = 5.0

    def design(self) -> SharedCrankRockerDesign:
        angle = math.radians(self.output_angle_degrees)
        along = self.rocker_ratio * math.cos(angle)
        normal = self.rocker_ratio * math.sin(angle)
        return SharedCrankRockerDesign(
            crank_ratio=self.crank_ratio,
            small_coupler_ratio=self.coupler_ratio,
            small_rocker_ratio=self.rocker_ratio,
            small_output_along_ratio=along,
            small_output_normal_ratio=normal,
            large_coupler_ratio=self.coupler_ratio,
            large_rocker_ratio=self.rocker_ratio,
            large_output_along_ratio=along,
            large_output_normal_ratio=-normal,
            small_slider_rod_ratio=self.connecting_rod_ratio,
            large_slider_rod_ratio=self.connecting_rod_ratio,
        )


@dataclass(frozen=True, slots=True)
class MechanicalDesignEvaluation:
    parameters: ReflectedRockerParameters
    periodic: PeriodicResult | None
    performance: CyclePerformance | None
    diagnostics: CycleDiagnostics | None
    validity: ValidityReport | None
    model: ThermodynamicModel | None
    kinematics: SharedCrankFourBarVolumeKinematics | None
    total_gas_mass: float
    error: str | None = None

    @property
    def converged(self) -> bool:
        return (
            self.periodic is not None
            and self.periodic.status is PeriodicStatus.CONVERGED
            and self.performance is not None
        )

    @property
    def satisfies_reference_constraints(self) -> bool:
        """Apply only the two currently selected cooling-cell constraints."""

        return bool(
            self.converged
            and self.performance is not None
            and self.performance.cooling_power >= 100.0
            and self.performance.mechanical_input_power <= 150.0
            and self.performance.cooling_cop is not None
        )


@dataclass(frozen=True, slots=True)
class ScaledPistonGeometry:
    """Physical scale implied by swept volumes and a selected large-piston ratio."""

    ground_distance: float
    small_stroke: float
    large_stroke: float
    small_bore: float
    large_bore: float
    small_bore_to_stroke_ratio: float
    large_bore_to_stroke_ratio: float


def scale_for_large_bore_to_stroke_ratio(
    unit_scale_kinematics: SharedCrankFourBarVolumeKinematics,
    *,
    current_ground_distance: float,
    target_large_bore_to_stroke_ratio: float,
) -> ScaledPistonGeometry:
    """Scale a dimensionless mechanism without changing its volume waveforms.

    Cylinder swept volumes remain fixed. Link lengths and physical strokes scale
    together, while piston areas are reconstructed from swept volume divided by
    stroke. The target ratio therefore determines a unique mechanism scale.
    """

    if current_ground_distance <= 0.0 or not math.isfinite(current_ground_distance):
        raise ValueError("Current ground distance must be finite and positive.")
    if (
        target_large_bore_to_stroke_ratio <= 0.0
        or not math.isfinite(target_large_bore_to_stroke_ratio)
    ):
        raise ValueError("Target bore-to-stroke ratio must be finite and positive.")
    large_swept = unit_scale_kinematics.large_volume_limits.swept
    stroke_per_ground = (
        unit_scale_kinematics.large_physical_stroke / current_ground_distance
    )
    target_large_stroke = (
        4.0 * large_swept
        / (math.pi * target_large_bore_to_stroke_ratio**2)
    ) ** (1.0 / 3.0)
    ground_distance = target_large_stroke / stroke_per_ground
    scale = ground_distance / current_ground_distance
    small_stroke = unit_scale_kinematics.small_physical_stroke * scale
    large_stroke = unit_scale_kinematics.large_physical_stroke * scale
    small_bore = math.sqrt(
        4.0 * unit_scale_kinematics.small_volume_limits.swept
        / (math.pi * small_stroke)
    )
    large_bore = math.sqrt(4.0 * large_swept / (math.pi * large_stroke))
    return ScaledPistonGeometry(
        ground_distance=ground_distance,
        small_stroke=small_stroke,
        large_stroke=large_stroke,
        small_bore=small_bore,
        large_bore=large_bore,
        small_bore_to_stroke_ratio=small_bore / small_stroke,
        large_bore_to_stroke_ratio=large_bore / large_stroke,
    )


class FixedInventoryMechanicalEvaluator:
    """Evaluate kinematics while holding thermodynamic hardware and mass fixed."""

    def __init__(self, configuration: SimulationConfiguration) -> None:
        if configuration.four_bar_ground_distance is None:
            raise ValueError("Mechanical evaluation requires a four-bar configuration.")
        self.configuration = configuration
        reference_model = build_model(configuration)
        self.total_gas_mass = build_initial_state(
            configuration, reference_model
        ).total_mass
        self._fixed_mass_configuration = replace(
            configuration,
            charge=ChargeConfiguration(
                temperature=configuration.charge.temperature,
                total_mass=self.total_gas_mass,
            ),
        )

    def build_candidate_model(
        self, parameters: ReflectedRockerParameters
    ) -> tuple[ThermodynamicModel, SharedCrankFourBarVolumeKinematics]:
        configuration = self._fixed_mass_configuration
        assert configuration.four_bar_ground_distance is not None
        assert configuration.four_bar_crank_direction is not None
        kinematics = shared_crank_rocker_kinematics(
            parameters.design(),
            configuration.machine_volumes.small_cylinder,
            configuration.machine_volumes.large_cylinder,
            ground_distance=configuration.four_bar_ground_distance,
            crank_angle_offset=math.radians(parameters.crank_angle_offset_degrees),
            crank_direction=configuration.four_bar_crank_direction,
        )
        return replace(build_model(configuration), kinematics=kinematics), kinematics

    def evaluate(
        self,
        parameters: ReflectedRockerParameters,
        warm_start: ThermodynamicState | None = None,
    ) -> MechanicalDesignEvaluation:
        try:
            model, kinematics = self.build_candidate_model(parameters)
            filling_state = build_initial_state(self._fixed_mass_configuration, model)
            initial_state = filling_state if warm_start is None else _fixed_mass_state(
                warm_start, self.total_gas_mass
            )
            periodic = build_periodic_solver(
                self._fixed_mass_configuration, model
            ).solve(initial_state, initial_valve_topology())
            if periodic.status is not PeriodicStatus.CONVERGED or periodic.final_cycle is None:
                return MechanicalDesignEvaluation(
                    parameters, periodic, None, None, None, model, kinematics,
                    self.total_gas_mass
                )
            cycle = periodic.final_cycle
            cycle_start = ThermodynamicState.from_array(cycle.states[:, 0])
            performance = calculate_cycle_performance(
                cycle, cycle_start, model.signed_angular_speed
            )
            return MechanicalDesignEvaluation(
                parameters,
                periodic,
                performance,
                extract_cycle_diagnostics(cycle, model),
                assess_cycle_validity(
                    cycle, model, self._fixed_mass_configuration.validity
                ),
                model,
                kinematics,
                self.total_gas_mass,
            )
        except (ValueError, FloatingPointError) as error:
            return MechanicalDesignEvaluation(
                parameters, None, None, None, None, None, None,
                self.total_gas_mass, str(error)
            )


def _fixed_mass_state(
    state: ThermodynamicState, target_total_mass: float
) -> ThermodynamicState:
    factor = target_total_mass / state.total_mass
    return ThermodynamicState.from_array(state.as_array() * factor)
