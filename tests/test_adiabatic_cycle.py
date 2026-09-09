import math

import numpy as np
import pytest

from dada_solver.dynamics import ThermodynamicModel, ValveTopology
from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.geometry import CylinderVolumeLimits, MachineVolumes
from dada_solver.heat_transfer import ReservoirHeatTransfer
from dada_solver.hydraulics import CompressibleOrifice, FlowResult
from dada_solver.integration import CycleIntegrator
from dada_solver.kinematics import HarmonicVolumeKinematics
from dada_solver.state import ThermodynamicState, UniformCharge
from dada_solver.valves import PassiveCheckValve, ValveState
from dada_solver.verification import closed_adiabatic_invariant


class ZeroInternalFlow:
    def bidirectional_flow(self, *args, **kwargs) -> FlowResult:
        return FlowResult(0.0, False)


def create_isolated_adiabatic_model(
    gas: CaloricallyPerfectGas,
) -> ThermodynamicModel:
    small = CylinderVolumeLimits(1.0e-4, 3.0e-4)
    large = CylinderVolumeLimits(2.0e-4, 6.0e-4)
    machine_volumes = MachineVolumes(small, large, 1.0e-4, 1.0e-4)
    blocked_valve = PassiveCheckValve(
        CompressibleOrifice(1.0e-6),
        opening_pressure_difference=1.0e12,
        closing_pressure_difference=0.0,
    )
    model = ThermodynamicModel(
        gas=gas,
        machine_volumes=machine_volumes,
        kinematics=HarmonicVolumeKinematics(small, large, math.pi),
        angular_speed=10.0,
        cold_heat_transfer=ReservoirHeatTransfer(0.0, 280.0),
        hot_heat_transfer=ReservoirHeatTransfer(0.0, 320.0),
        large_hot_link=CompressibleOrifice(1.0e-6),
        small_cold_link=CompressibleOrifice(1.0e-6),
        hot_small_valve=blocked_valve,
        cold_large_valve=blocked_valve,
    )
    object.__setattr__(model, "large_hot_link", ZeroInternalFlow())
    object.__setattr__(model, "small_cold_link", ZeroInternalFlow())
    return model


def maximum_relative_invariant_error(
    model: ThermodynamicModel,
    angles: np.ndarray,
    states: np.ndarray,
) -> float:
    gamma = model.gas.heat_capacity_ratio
    initial_state = ThermodynamicState.from_array(states[:, 0])
    initial_pressures = initial_state.pressures(model.gas, model.volumes(angles[0]))
    initial_volumes = model.volumes(angles[0])
    references = (
        closed_adiabatic_invariant(
            initial_pressures[0], initial_volumes.small_cylinder, gamma
        ),
        closed_adiabatic_invariant(
            initial_pressures[1], initial_volumes.large_cylinder, gamma
        ),
    )
    errors = []
    for angle, values in zip(angles, states.T, strict=True):
        state = ThermodynamicState.from_array(values)
        volumes = model.volumes(float(angle))
        pressures = state.pressures(model.gas, volumes)
        invariants = (
            closed_adiabatic_invariant(pressures[0], volumes.small_cylinder, gamma),
            closed_adiabatic_invariant(pressures[1], volumes.large_cylinder, gamma),
        )
        errors.extend(
            abs(value - reference) / reference
            for value, reference in zip(invariants, references, strict=True)
        )
    return max(errors)


def test_mobile_closed_adiabatic_cycle_preserves_analytical_invariants(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    model = create_isolated_adiabatic_model(ideal_gas)
    initial_state = UniformCharge(2.0e5, 300.0).create_state(
        ideal_gas, model.volumes(0.0)
    )
    topology = ValveTopology(ValveState.CLOSED, ValveState.CLOSED)
    cycle = CycleIntegrator(
        model,
        relative_tolerance=1.0e-10,
        absolute_tolerance=1.0e-12,
        maximum_step_angle=math.radians(0.5),
    ).integrate_cycle(initial_state, topology)

    assert cycle.completed
    assert cycle.events == ()
    assert maximum_relative_invariant_error(
        model, cycle.angles, cycle.states
    ) < 2.0e-10
    np.testing.assert_allclose(
        cycle.final_state.as_array(),
        initial_state.as_array(),
        rtol=2.0e-10,
        atol=1.0e-12,
    )


def test_adiabatic_invariant_error_decreases_with_numerical_refinement(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    model = create_isolated_adiabatic_model(ideal_gas)
    initial_state = UniformCharge(2.0e5, 300.0).create_state(
        ideal_gas, model.volumes(0.0)
    )
    topology = ValveTopology(ValveState.CLOSED, ValveState.CLOSED)
    coarse = CycleIntegrator(
        model,
        relative_tolerance=1.0e-4,
        absolute_tolerance=1.0e-7,
        maximum_step_angle=math.radians(20.0),
    ).integrate_cycle(initial_state, topology)
    fine = CycleIntegrator(
        model,
        relative_tolerance=1.0e-9,
        absolute_tolerance=1.0e-12,
        maximum_step_angle=math.radians(1.0),
    ).integrate_cycle(initial_state, topology)

    coarse_error = maximum_relative_invariant_error(
        model, coarse.angles, coarse.states
    )
    fine_error = maximum_relative_invariant_error(model, fine.angles, fine.states)

    assert fine_error < coarse_error
    assert fine_error < 1.0e-8
