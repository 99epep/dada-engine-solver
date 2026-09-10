"""Construction of simulation objects from validated configuration."""

from __future__ import annotations

import math
from dataclasses import replace

from dada_solver.configuration import SimulationConfiguration
from dada_solver.free_kinematics import FreeKinematics
from dada_solver.dynamics import ThermodynamicModel, ValveTopology
from dada_solver.heat_transfer import ReservoirHeatTransfer
from dada_solver.hydraulics import CompressibleOrifice
from dada_solver.integration import CycleIntegrator, SegmentIntegrationProgress
from typing import Callable
from dada_solver.kinematics import KinematicsModel, HarmonicVolumeKinematics, IdealPiecewiseLinearVolumeKinematics, ReversedVolumeKinematics
from dada_solver.four_bar import (
    SharedCrankFourBarVolumeKinematics,
    published_e0_opposed_kinematics,
    published_f65_opposed_kinematics,
    shared_crank_rocker_kinematics,
)
from dada_solver.periodic import SuccessiveCycleSolver
from dada_solver.state import ThermodynamicState, UniformCharge
from dada_solver.valves import PassiveCheckValve, ValveState


def build_model(configuration: SimulationConfiguration, *, kinematics: KinematicsModel | None = None) -> ThermodynamicModel:
    """Build the first-level model without adding unconfigured physical data."""

    if kinematics is not None:
        # Injected implementations use study angle, before the operation reversal.
        if isinstance(kinematics, ReversedVolumeKinematics):
            raise ValueError('Inject study-angle kinematics, not an already reversed wrapper.')
        if (kinematics.small_volume_limits != configuration.machine_volumes.small_cylinder
                or kinematics.large_volume_limits != configuration.machine_volumes.large_cylinder):
            raise ValueError('Injected kinematics and configuration volume limits must agree.')
        validator = getattr(kinematics, 'require_feasible', None)
        if validator is not None:
            validator()
    elif configuration.kinematics_type == "free":
        kinematics = FreeKinematics(configuration.free_kinematics)
        kinematics.require_feasible()
    elif configuration.kinematics_type == "harmonic_example":
        assert configuration.small_phase_offset_degrees is not None
        kinematics = HarmonicVolumeKinematics(
            configuration.machine_volumes.small_cylinder,
            configuration.machine_volumes.large_cylinder,
            math.radians(configuration.small_phase_offset_degrees),
        )
    elif configuration.kinematics_type == "ideal_piecewise_linear":
        assert configuration.small_lambda_target is not None
        assert configuration.large_lambda_target is not None
        assert configuration.adiabatic_sector_fraction is not None
        kinematics = IdealPiecewiseLinearVolumeKinematics(
            configuration.machine_volumes.small_cylinder,
            configuration.machine_volumes.large_cylinder,
            configuration.small_lambda_target,
            configuration.large_lambda_target,
            configuration.adiabatic_sector_fraction,
        )
    elif configuration.kinematics_type == "published_f65_opposed":
        assert configuration.four_bar_ground_distance is not None
        assert configuration.four_bar_connecting_rod_ratio is not None
        assert configuration.four_bar_crank_angle_offset_degrees is not None
        assert configuration.four_bar_crank_direction is not None
        kinematics = published_f65_opposed_kinematics(
            configuration.machine_volumes.small_cylinder,
            configuration.machine_volumes.large_cylinder,
            ground_distance=configuration.four_bar_ground_distance,
            connecting_rod_to_projected_stroke_ratio=(
                configuration.four_bar_connecting_rod_ratio
            ),
            crank_angle_offset=math.radians(
                configuration.four_bar_crank_angle_offset_degrees
            ),
            crank_direction=configuration.four_bar_crank_direction,
        )
    elif configuration.kinematics_type == "published_e0_opposed":
        assert configuration.four_bar_ground_distance is not None
        assert configuration.four_bar_connecting_rod_ratio is not None
        assert configuration.four_bar_crank_angle_offset_degrees is not None
        assert configuration.four_bar_crank_direction is not None
        kinematics = published_e0_opposed_kinematics(
            configuration.machine_volumes.small_cylinder,
            configuration.machine_volumes.large_cylinder,
            ground_distance=configuration.four_bar_ground_distance,
            connecting_rod_to_projected_stroke_ratio=(
                configuration.four_bar_connecting_rod_ratio
            ),
            crank_angle_offset=math.radians(
                configuration.four_bar_crank_angle_offset_degrees
            ),
            crank_direction=configuration.four_bar_crank_direction,
        )
    else:
        assert configuration.shared_four_bar_design is not None
        assert configuration.four_bar_ground_distance is not None
        assert configuration.four_bar_crank_angle_offset_degrees is not None
        assert configuration.four_bar_crank_direction is not None
        kinematics = shared_crank_rocker_kinematics(
            configuration.shared_four_bar_design,
            configuration.machine_volumes.small_cylinder,
            configuration.machine_volumes.large_cylinder,
            ground_distance=configuration.four_bar_ground_distance,
            crank_angle_offset=math.radians(
                configuration.four_bar_crank_angle_offset_degrees
            ),
            crank_direction=configuration.four_bar_crank_direction,
        )
    if configuration.motor_operation:
        if isinstance(kinematics, SharedCrankFourBarVolumeKinematics):
            # Keep the mechanical adapter available to animation and sizing.
            kinematics = replace(kinematics, crank_direction=-kinematics.crank_direction)
        else:
            kinematics = ReversedVolumeKinematics(kinematics)
    hydraulics = configuration.hydraulics
    flow_models = configuration.hydraulic_flow_models
    large_hot_model = (
        CompressibleOrifice(
            hydraulics.large_to_hot_cda,
            hydraulics.orifice_pressure_regularization,
        )
        if flow_models is None else flow_models.large_to_hot
    )
    small_cold_model = (
        CompressibleOrifice(
            hydraulics.small_to_cold_cda,
            hydraulics.orifice_pressure_regularization,
        )
        if flow_models is None else flow_models.small_to_cold
    )
    hot_small_model = (
        CompressibleOrifice(
            hydraulics.hot_to_small_valve_cda,
            hydraulics.orifice_pressure_regularization,
        )
        if flow_models is None else flow_models.hot_to_small_valve
    )
    cold_large_model = (
        CompressibleOrifice(
            hydraulics.cold_to_large_valve_cda,
            hydraulics.orifice_pressure_regularization,
        )
        if flow_models is None else flow_models.cold_to_large_valve
    )
    return ThermodynamicModel(
        gas=configuration.gas,
        machine_volumes=configuration.machine_volumes,
        kinematics=kinematics,
        angular_speed=abs(configuration.angular_speed),
        study_crank_direction=-1 if configuration.motor_operation else 1,
        cold_heat_transfer=ReservoirHeatTransfer(
            configuration.cold_thermal_conductance,
            configuration.heat_in_reservoir_temperature,
        ),
        hot_heat_transfer=ReservoirHeatTransfer(
            configuration.hot_thermal_conductance,
            configuration.heat_out_reservoir_temperature,
        ),
        large_hot_link=large_hot_model,
        small_cold_link=small_cold_model,
        hot_small_valve=PassiveCheckValve(
            hot_small_model,
            configuration.hot_to_small_valve.opening_pressure_difference,
            configuration.hot_to_small_valve.closing_pressure_difference,
        ),
        cold_large_valve=PassiveCheckValve(
            cold_large_model,
            configuration.cold_to_large_valve.opening_pressure_difference,
            configuration.cold_to_large_valve.closing_pressure_difference,
        ),
        continuous_ideal_diodes=(
            configuration.valve_model == "continuous_ideal_diode"
        ),
    )


def build_initial_state(
    configuration: SimulationConfiguration,
    model: ThermodynamicModel,
) -> ThermodynamicState:
    """Build the uniform filling state at theta=0, not a periodic state."""

    volumes = model.volumes(0.0)
    pressure = configuration.charge.resolved_pressure(
        configuration.gas, volumes.total
    )
    return UniformCharge(pressure, configuration.charge.temperature).create_state(
        configuration.gas, volumes
    )


def build_periodic_solver(
    configuration: SimulationConfiguration,
    model: ThermodynamicModel,
    integration_progress_callback: Callable[[SegmentIntegrationProgress], None] | None = None,
) -> SuccessiveCycleSolver:
    numerical = configuration.numerical
    integrator = CycleIntegrator(
        model=model,
        relative_tolerance=numerical.integration_relative_tolerance,
        absolute_tolerance=numerical.integration_absolute_tolerance,
        maximum_step_angle=math.radians(numerical.maximum_step_angle_degrees),
        maximum_events=numerical.maximum_valve_events_per_cycle,
        integration_method=numerical.integration_method,
        progress_callback=integration_progress_callback,
    )
    return SuccessiveCycleSolver(
        integrator=integrator,
        maximum_cycles=numerical.maximum_cycles,
        relative_tolerance=numerical.periodic_relative_tolerance,
        mass_absolute_tolerance=numerical.periodic_mass_absolute_tolerance,
        energy_absolute_tolerance=numerical.periodic_energy_absolute_tolerance,
    )


def initial_valve_topology() -> ValveTopology:
    """Return the study's theta=0 nominal initial guess; passive checks follow."""

    return ValveTopology(ValveState.CLOSED, ValveState.CLOSED)
