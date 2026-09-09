from dataclasses import dataclass

from dada_solver.dynamics import ThermodynamicModel, ValveTopology
from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.geometry import CylinderVolumeLimits, MachineVolumes
from dada_solver.heat_transfer import ReservoirHeatTransfer
from dada_solver.hydraulics import CompressibleOrifice
from dada_solver.integration import CycleIntegrator, IntegrationStatus
from dada_solver.periodic import PeriodicStatus, SuccessiveCycleSolver
from dada_solver.state import UniformCharge
from dada_solver.valves import PassiveCheckValve, ValveState


@dataclass(frozen=True)
class ConstantKinematics:
    small_volume: float
    large_volume: float

    def small_cylinder_volume(self, theta: float) -> float:
        return self.small_volume

    def small_cylinder_volume_derivative(self, theta: float) -> float:
        return 0.0

    def large_cylinder_volume(self, theta: float) -> float:
        return self.large_volume

    def large_cylinder_volume_derivative(self, theta: float) -> float:
        return 0.0


def test_periodic_state_converges_on_static_controlled_case(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    machine_volumes = MachineVolumes(
        small_cylinder=CylinderVolumeLimits(1.0e-4, 3.0e-4),
        large_cylinder=CylinderVolumeLimits(2.0e-4, 6.0e-4),
        cold_heat_exchanger=1.0e-4,
        hot_heat_exchanger=1.0e-4,
    )
    kinematics = ConstantKinematics(2.0e-4, 5.0e-4)
    valve = PassiveCheckValve(CompressibleOrifice(1.0e-6), 1000.0, 0.0)
    model = ThermodynamicModel(
        gas=ideal_gas,
        machine_volumes=machine_volumes,
        kinematics=kinematics,
        angular_speed=10.0,
        cold_heat_transfer=ReservoirHeatTransfer(10.0, 300.0),
        hot_heat_transfer=ReservoirHeatTransfer(10.0, 300.0),
        large_hot_link=CompressibleOrifice(1.0e-6),
        small_cold_link=CompressibleOrifice(1.0e-6),
        hot_small_valve=valve,
        cold_large_valve=valve,
    )
    initial_state = UniformCharge(2.0e5, 300.0).create_state(
        ideal_gas, model.volumes(0.0)
    )
    initial_topology = ValveTopology(ValveState.CLOSED, ValveState.CLOSED)
    integrator = CycleIntegrator(model)

    cycle = integrator.integrate_cycle(initial_state, initial_topology)
    result = SuccessiveCycleSolver(integrator, maximum_cycles=2).solve(
        initial_state, initial_topology
    )

    assert cycle.status is IntegrationStatus.COMPLETED
    assert result.status is PeriodicStatus.CONVERGED
    assert len(result.history) == 1
    assert result.history[0].normalized_state_error == 0.0
    assert result.history[0].absolute_mass_residual == 0.0
    assert cycle.statistics.elapsed_seconds >= 0.0
    assert cycle.statistics.solve_segments == 1
    assert cycle.statistics.right_hand_side_evaluations > 0
    assert cycle.statistics.event_function_evaluations > 0
    assert cycle.statistics.stored_samples == cycle.angles.size


def create_static_cycle(ideal_gas: CaloricallyPerfectGas):
    """Return a completed equilibrium cycle and its model for diagnostics tests."""
    machine_volumes = MachineVolumes(
        small_cylinder=CylinderVolumeLimits(1.0e-4, 3.0e-4),
        large_cylinder=CylinderVolumeLimits(2.0e-4, 6.0e-4),
        cold_heat_exchanger=1.0e-4,
        hot_heat_exchanger=1.0e-4,
    )
    model = ThermodynamicModel(
        gas=ideal_gas,
        machine_volumes=machine_volumes,
        kinematics=ConstantKinematics(2.0e-4, 5.0e-4),
        angular_speed=10.0,
        cold_heat_transfer=ReservoirHeatTransfer(10.0, 300.0),
        hot_heat_transfer=ReservoirHeatTransfer(10.0, 300.0),
        large_hot_link=CompressibleOrifice(1.0e-6),
        small_cold_link=CompressibleOrifice(1.0e-6),
        hot_small_valve=PassiveCheckValve(
            CompressibleOrifice(1.0e-6), 1000.0, 0.0
        ),
        cold_large_valve=PassiveCheckValve(
            CompressibleOrifice(1.0e-6), 1000.0, 0.0
        ),
    )
    state = UniformCharge(2.0e5, 300.0).create_state(
        ideal_gas, model.volumes(0.0)
    )
    topology = ValveTopology(ValveState.CLOSED, ValveState.CLOSED)
    return model, state, CycleIntegrator(model).integrate_cycle(state, topology)


def test_state_is_continuous_when_valves_open_at_initial_angle(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    model, state, _ = create_static_cycle(ideal_gas)
    values = state.as_array()
    values[5] *= 1.02
    values[7] *= 1.02
    from dada_solver.state import ThermodynamicState

    initial_state = ThermodynamicState.from_array(values)
    cycle = CycleIntegrator(model).integrate_cycle(
        initial_state,
        ValveTopology(ValveState.CLOSED, ValveState.CLOSED),
    )

    assert cycle.completed
    assert len(cycle.events) >= 2
    assert all(event.angle == 0.0 for event in cycle.events[:2])
    # A topology event performs no state reset before integration starts.
    import numpy as np

    np.testing.assert_array_equal(cycle.states[:, 0], initial_state.as_array())
