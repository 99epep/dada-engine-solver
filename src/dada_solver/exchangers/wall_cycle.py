"""Reusable periodic evaluation for the existing ten-state air-wall motor."""
from __future__ import annotations

from dataclasses import dataclass
import math
from types import SimpleNamespace
import numpy as np

from dada_solver.dynamics import ValveTopology
from dada_solver.exchangers.air_wall import AirWallMotor
from dada_solver.integration import IntegrationInterrupted
from dada_solver.performance import ConservationReport, CyclePerformance, OperatingMode
from dada_solver.state import ThermodynamicState
from dada_solver.valves import ValveState


@dataclass(frozen=True)
class WallCycleResult:
    status: str
    message: str
    history: tuple[dict, ...]
    angles: np.ndarray | None
    trajectory: np.ndarray | None
    last_complete_state: np.ndarray | None

    @property
    def converged(self):
        return self.status == 'converged'


def solve_periodic_wall_motor(wrapper: AirWallMotor, initial_state, *, maximum_cycles,
        progress_callback=None, accelerate_walls=False, cycle_callback=None) -> WallCycleResult:
    """Repeat unmodified complete cycles using the screening convergence rule."""
    from dada_solver.exchangers.wall_iteration import extrapolate_wall_energies
    state = np.asarray(initial_state, dtype=float)
    history, wall_history = [], []
    last_angles = last_trajectory = None
    atol = np.array([1e-13, 1e-8]*4 + [1e-7]*2 + [1e-8]*5)
    for cycle in range(1, maximum_cycles+1):
        try:
            angles, trajectory = wrapper.integrate_cycle(state, rtol=1e-8,
                atol=atol, maximum_step_angle=math.pi/360,
                progress_callback=progress_callback, progress_interval_seconds=1.0)
        except IntegrationInterrupted as error:
            return WallCycleResult('interrupted', str(error), tuple(history),
                last_angles, last_trajectory, state.copy() if history else None)
        end = trajectory[:10, -1]
        error = float(np.max(np.abs(end-state) /
            (1e-12 + 1e-6*np.maximum(np.abs(end), np.abs(state)))))
        history.append(dict(cycle=cycle, normalized_state_error=error))
        last_angles, last_trajectory = angles, trajectory
        state = end
        if cycle_callback is not None:
            cycle_callback(cycle, end.copy(), error, history[-1])
        if error <= 1:
            return WallCycleResult('converged', 'Periodic steady state converged.',
                tuple(history), angles, trajectory, end.copy())
        wall_history.append(end[8:10].copy())
        if accelerate_walls and cycle % 10 == 0 and cycle < maximum_cycles and len(wall_history) >= 3:
            proposed = extrapolate_wall_energies(*wall_history[-3:], np.array([
                wrapper.heat_in.wall_capacity_j_k, wrapper.heat_out.wall_capacity_j_k]))
            if proposed is not None:
                state = state.copy(); state[8:10] = proposed
                history[-1]['wall_initial_guess_extrapolated'] = True
            wall_history = []
    return WallCycleResult('maximum_cycles',
        'Maximum cycle count reached before periodic convergence.', tuple(history),
        last_angles, last_trajectory, state.copy())


class WallDiagnosticCycle:
    """Completed-cycle view consumed by existing diagnostics and mechanics."""
    def __init__(self, angles, trajectory):
        self.angles = angles
        self.states = trajectory[:8]
        self.cold_heat = trajectory[10]
        self.hot_heat = trajectory[11]
        self.gas_work = trajectory[14]
        topology = ValveTopology(ValveState.CLOSED, ValveState.CLOSED)
        self.topologies = tuple(topology for _ in angles)
        self.final_topology = topology
        self.events = ()
        self.completed = True

    @property
    def final_state(self):
        return ThermodynamicState.from_array(self.states[:, -1])


def wall_cycle_performance(wrapper, trajectory):
    """Use external-air heat, including wall storage, at the machine boundary."""
    frequency = abs(wrapper.model.signed_angular_speed)/(2*math.pi)
    heat_in = float(trajectory[10, -1]-trajectory[10, 0])
    heat_out = float(trajectory[11, -1]-trajectory[11, 0])
    work = float(trajectory[14, -1]-trajectory[14, 0])
    gas0, gas1 = ThermodynamicState.from_array(trajectory[:8, 0]), ThermodynamicState.from_array(trajectory[:8, -1])
    mass = gas1.total_mass-gas0.total_mass
    stored = gas1.total_internal_energy-gas0.total_internal_energy + float(np.sum(trajectory[8:10, -1]-trajectory[8:10, 0]))
    residual = stored-heat_in-heat_out+work
    scale = max(abs(gas0.total_internal_energy)+abs(float(np.sum(trajectory[8:10, 0]))),
                abs(heat_in)+abs(heat_out)+abs(work), np.finfo(float).tiny)
    motor = wrapper.model.signed_angular_speed < 0 and heat_in > 0 and heat_out < 0 and work > 0
    return CyclePerformance(heat_in, heat_out, work, -work, None, None,
        heat_in*frequency, -heat_out*frequency, -work*frequency,
        OperatingMode.MOTOR if motor else OperatingMode.NON_REFRIGERATION,
        ConservationReport(mass, mass/gas0.total_mass, residual, residual/scale))


def convergence_summary(history):
    errors = [float(item['normalized_state_error']) for item in history]
    return dict(cycles_completed=len(errors), first_normalized_periodic_error=errors[0] if errors else None,
        last_normalized_periodic_error=errors[-1] if errors else None,
        normalized_periodic_error_ratio=(errors[-1]/errors[0] if errors and errors[0] else None),
        improving=bool(len(errors) >= 2 and errors[-1] < errors[0]), history=list(history))
