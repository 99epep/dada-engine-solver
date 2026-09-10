"""Event-driven integration of one imposed-kinematics cycle."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math
import time
from typing import Callable

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from dada_solver.dynamics import ThermodynamicModel, ValveTopology
from dada_solver.state import ThermodynamicState
from dada_solver.valves import ValveState


class IntegrationStatus(Enum):
    COMPLETED = "completed"
    INVALID_PHYSICAL_STATE = "invalid_physical_state"
    NUMERICAL_FAILURE = "numerical_integration_failure"
    MAXIMUM_EVENT_COUNT_REACHED = "maximum_event_count_reached"
    INTERRUPTED = "interrupted"


class IntegrationInterrupted(RuntimeError):
    """Cooperative interruption requested at an integration progress boundary."""


@dataclass(frozen=True, slots=True)
class ValveEvent:
    angle: float
    valve_name: str
    previous_state: ValveState
    current_state: ValveState


@dataclass(frozen=True, slots=True)
class IntegrationStatistics:
    elapsed_seconds: float
    solve_segments: int
    right_hand_side_evaluations: int
    jacobian_evaluations: int
    linear_decompositions: int
    event_function_evaluations: int
    stored_samples: int
    segments: tuple["SegmentIntegrationStatistics", ...]


@dataclass(frozen=True, slots=True)
class SegmentIntegrationStatistics:
    start_angle: float
    end_angle: float
    hot_to_small_state: ValveState
    cold_to_large_state: ValveState
    elapsed_seconds: float
    right_hand_side_evaluations: int
    jacobian_evaluations: int
    linear_decompositions: int
    stored_samples: int
    ended_at_event: bool


@dataclass(frozen=True, slots=True)
class SegmentIntegrationProgress:
    """Optional live progress from a potentially stiff solve segment."""

    phase: str
    start_angle: float
    current_angle: float
    target_angle: float
    hot_to_small_state: ValveState
    cold_to_large_state: ValveState
    elapsed_seconds: float
    right_hand_side_evaluations: int


@dataclass(slots=True)
class _IntegrationCounters:
    start_time: float
    solve_segments: int = 0
    right_hand_side_evaluations: int = 0
    jacobian_evaluations: int = 0
    linear_decompositions: int = 0
    event_function_evaluations: int = 0
    segments: list[SegmentIntegrationStatistics] | None = None

    def __post_init__(self) -> None:
        if self.segments is None:
            self.segments = []


@dataclass(frozen=True, slots=True)
class CycleIntegrationResult:
    status: IntegrationStatus
    message: str
    angles: NDArray[np.float64]
    states: NDArray[np.float64]
    cold_heat: NDArray[np.float64]
    hot_heat: NDArray[np.float64]
    gas_work: NDArray[np.float64]
    topologies: tuple[ValveTopology, ...]
    events: tuple[ValveEvent, ...]
    final_topology: ValveTopology
    statistics: IntegrationStatistics

    @property
    def completed(self) -> bool:
        return self.status is IntegrationStatus.COMPLETED

    @property
    def final_state(self) -> ThermodynamicState:
        return ThermodynamicState.from_array(self.states[:, -1])


@dataclass(frozen=True, slots=True)
class CycleIntegrator:
    """Integrate a cycle while locating every passive-valve threshold crossing."""

    model: ThermodynamicModel
    relative_tolerance: float = 1.0e-8
    absolute_tolerance: float = 1.0e-10
    maximum_step_angle: float = math.radians(1.0)
    maximum_events: int = 100
    integration_method: str = "RK45"
    use_jacobian_sparsity: bool = True
    progress_callback: Callable[[SegmentIntegrationProgress], None] | None = None
    progress_interval_seconds: float = 10.0

    def __post_init__(self) -> None:
        if self.relative_tolerance <= 0.0 or self.absolute_tolerance <= 0.0:
            raise ValueError("Integration tolerances must be positive.")
        if self.maximum_step_angle <= 0.0:
            raise ValueError("Maximum angular step must be positive.")
        if self.maximum_events <= 0:
            raise ValueError("Maximum event count must be positive.")
        if self.integration_method not in {"RK45", "Radau", "BDF", "LSODA"}:
            raise ValueError("Unsupported integration method.")
        if self.progress_interval_seconds <= 0.0:
            raise ValueError("Progress interval must be positive.")

    def integrate_cycle(
        self,
        initial_state: ThermodynamicState,
        initial_topology: ValveTopology,
    ) -> CycleIntegrationResult:
        """Integrate from theta=0 to theta=2*pi with continuous conserved state."""

        theta = 0.0
        end_angle = 2.0 * math.pi
        topology = initial_topology
        augmented_state = np.concatenate((initial_state.as_array(), np.zeros(3)))
        angle_parts: list[NDArray[np.float64]] = []
        state_parts: list[NDArray[np.float64]] = []
        topology_samples: list[ValveTopology] = []
        events: list[ValveEvent] = []
        counters = _IntegrationCounters(time.perf_counter())

        if not self.model.continuous_ideal_diodes:
            topology = self._apply_transitions(
                theta, initial_state, topology, events
            )
        breakpoint_provider = getattr(self.model.kinematics, "breakpoint_angles", None)
        kinematic_breakpoints = (
            tuple(breakpoint_provider()) if breakpoint_provider is not None else ()
        )

        try:
            while theta < end_angle:
                segment_start = theta
                segment_start_time = time.perf_counter()
                segment_end = next(
                    (angle for angle in kinematic_breakpoints if angle > theta + 1.0e-12),
                    end_angle,
                )
                event_functions = (
                    None
                    if self.model.continuous_ideal_diodes
                    else self._event_functions(topology, counters)
                )
                jacobian_options = (
                    {"jac_sparsity": self._jacobian_sparsity()}
                    if self.use_jacobian_sparsity
                    and self.integration_method in {"Radau", "BDF"}
                    else {}
                )
                segment_rhs_evaluations = 0
                last_progress_time = segment_start_time

                def right_hand_side(angle, values):
                    nonlocal segment_rhs_evaluations, last_progress_time
                    segment_rhs_evaluations += 1
                    now = time.perf_counter()
                    if (
                        self.progress_callback is not None
                        and now - last_progress_time >= self.progress_interval_seconds
                    ):
                        self.progress_callback(
                            SegmentIntegrationProgress(
                                phase="running",
                                start_angle=segment_start,
                                current_angle=float(angle),
                                target_angle=segment_end,
                                hot_to_small_state=topology.hot_to_small,
                                cold_to_large_state=topology.cold_to_large,
                                elapsed_seconds=now - segment_start_time,
                                right_hand_side_evaluations=segment_rhs_evaluations,
                            )
                        )
                        last_progress_time = now
                    return self._right_hand_side(angle, values, topology)

                solution = solve_ivp(
                    right_hand_side,
                    (theta, segment_end),
                    augmented_state,
                    rtol=self.relative_tolerance,
                    atol=self.absolute_tolerance,
                    max_step=self.maximum_step_angle,
                    method=self.integration_method,
                    events=event_functions,
                    dense_output=self.model.continuous_ideal_diodes,
                    **jacobian_options,
                )
                counters.solve_segments += 1
                counters.right_hand_side_evaluations += solution.nfev
                counters.jacobian_evaluations += solution.njev
                counters.linear_decompositions += solution.nlu
                root_observed = not self.model.continuous_ideal_diodes and any(
                    event_times.size for event_times in solution.t_events
                )
                assert counters.segments is not None
                counters.segments.append(
                    SegmentIntegrationStatistics(
                        start_angle=segment_start,
                        end_angle=float(solution.t[-1]),
                        hot_to_small_state=topology.hot_to_small,
                        cold_to_large_state=topology.cold_to_large,
                        elapsed_seconds=time.perf_counter() - segment_start_time,
                        right_hand_side_evaluations=solution.nfev,
                        jacobian_evaluations=solution.njev,
                        linear_decompositions=solution.nlu,
                        stored_samples=solution.t.size,
                        ended_at_event=root_observed,
                    )
                )
                if not solution.success:
                    return self._result(
                        IntegrationStatus.NUMERICAL_FAILURE,
                        solution.message,
                        angle_parts,
                        state_parts,
                        topology_samples,
                        events,
                        topology,
                        counters,
                    )

                # Keep both samples at a valve event: the preceding segment ends
                # with the old topology and the next one starts with the new
                # topology. Their conserved states must be exactly identical.
                start_index = 0
                angle_parts.append(solution.t[start_index:])
                state_parts.append(solution.y[:, start_index:])
                if self.model.continuous_ideal_diodes:
                    self._record_continuous_diode_events(
                        solution,
                        topology,
                        events,
                    )
                    segment_topologies = [
                        self.model.effective_topology(
                            float(angle),
                            ThermodynamicState.from_array(values[:8]),
                            topology,
                        )
                        for angle, values in zip(
                            solution.t[start_index:],
                            solution.y[:, start_index:].T,
                            strict=True,
                        )
                    ]
                    topology_samples.extend(segment_topologies)
                    topology = segment_topologies[-1]
                else:
                    topology_samples.extend(
                        [topology] * solution.t[start_index:].size
                    )
                augmented_state = solution.y[:, -1].copy()
                theta = float(solution.t[-1])

                if theta >= end_angle:
                    break
                if not root_observed:
                    # The solver stopped exactly at an imposed-kinematics slope
                    # discontinuity. State and valve topology remain continuous.
                    continue
                if len(events) >= self.maximum_events:
                    return self._result(
                        IntegrationStatus.MAXIMUM_EVENT_COUNT_REACHED,
                        "Maximum valve event count reached.",
                        angle_parts,
                        state_parts,
                        topology_samples,
                        events,
                        topology,
                        counters,
                    )

                event_state = ThermodynamicState.from_array(augmented_state[:8])
                previous_topology = topology
                triggered_indices = tuple(
                    index
                    for index, event_times in enumerate(solution.t_events)
                    if event_times.size
                    and math.isclose(
                        float(event_times[-1]), theta, rel_tol=0.0, abs_tol=1.0e-12
                    )
                )
                topology = self._apply_root_transitions(
                    theta,
                    event_state,
                    topology,
                    events,
                    triggered_indices,
                )
                if topology == previous_topology:
                    return self._result(
                        IntegrationStatus.NUMERICAL_FAILURE,
                        "A valve event root did not change the hydraulic topology.",
                        angle_parts,
                        state_parts,
                        topology_samples,
                        events,
                        topology,
                        counters,
                    )

        except IntegrationInterrupted:
            raise
        except ValueError as error:
            return self._result(
                IntegrationStatus.INVALID_PHYSICAL_STATE,
                str(error),
                angle_parts,
                state_parts,
                topology_samples,
                events,
                topology,
                counters,
            )
        except (FloatingPointError, OverflowError) as error:
            return self._result(
                IntegrationStatus.NUMERICAL_FAILURE,
                str(error),
                angle_parts,
                state_parts,
                topology_samples,
                events,
                topology,
                counters,
            )

        return self._result(
            IntegrationStatus.COMPLETED,
            "Cycle integration completed.",
            angle_parts,
            state_parts,
            topology_samples,
            events,
            topology,
            counters,
        )

    def _right_hand_side(
        self,
        theta: float,
        augmented_state: NDArray[np.float64],
        topology: ValveTopology,
    ) -> NDArray[np.float64]:
        state = ThermodynamicState.from_array(augmented_state[:8])
        rates = self.model.evaluate(theta, state, topology)
        derivative = np.empty(11, dtype=float)
        derivative[:8] = rates.state_derivative / self.model.angular_speed
        derivative[8] = rates.cold_heat_rate / self.model.angular_speed
        derivative[9] = rates.hot_heat_rate / self.model.angular_speed
        derivative[10] = rates.gas_work_rate / self.model.angular_speed
        return derivative

    @staticmethod
    def _jacobian_sparsity() -> NDArray[np.bool_]:
        """Return a topology-independent superset of exact state dependencies."""

        pattern = np.zeros((11, 11), dtype=bool)
        # Each mass/energy pair depends on its own state and its two neighboring
        # control volumes in the four-node hydraulic cycle.
        dependencies = (
            (0, 2, 3),  # S depends on S, C and H
            (1, 2, 3),  # L depends on L, C and H
            (0, 1, 2),  # C depends on C, S and L
            (0, 1, 3),  # H depends on H, S and L
        )
        for volume, neighbors in enumerate(dependencies):
            rows = (2 * volume, 2 * volume + 1)
            columns = tuple(
                component
                for neighbor in neighbors
                for component in (2 * neighbor, 2 * neighbor + 1)
            )
            pattern[np.ix_(rows, columns)] = True
        pattern[8, 4:6] = True  # cold heat integral depends on C
        pattern[9, 6:8] = True  # hot heat integral depends on H
        pattern[10, 0:4] = True  # piston work depends on S and L
        return pattern

    def _event_functions(
        self,
        topology: ValveTopology,
        counters: _IntegrationCounters,
        *,
        continuous: bool = False,
    ) -> tuple[Callable[[float, NDArray[np.float64]], float], ...]:
        def hot_small(theta: float, values: NDArray[np.float64]) -> float:
            counters.event_function_evaluations += 1
            state = ThermodynamicState.from_array(values[:8])
            pressures = state.pressures(self.model.gas, self.model.volumes(theta))
            threshold = 0.0 if continuous else (
                self.model.hot_small_valve.opening_pressure_difference
                if topology.hot_to_small is ValveState.CLOSED
                else self.model.hot_small_valve.closing_pressure_difference
            )
            return float(pressures[3] - pressures[0] - threshold)

        def cold_large(theta: float, values: NDArray[np.float64]) -> float:
            counters.event_function_evaluations += 1
            state = ThermodynamicState.from_array(values[:8])
            pressures = state.pressures(self.model.gas, self.model.volumes(theta))
            threshold = 0.0 if continuous else (
                self.model.cold_large_valve.opening_pressure_difference
                if topology.cold_to_large is ValveState.CLOSED
                else self.model.cold_large_valve.closing_pressure_difference
            )
            return float(pressures[2] - pressures[1] - threshold)

        hot_small.terminal = not continuous  # type: ignore[attr-defined]
        cold_large.terminal = not continuous  # type: ignore[attr-defined]
        hot_small.direction = (0.0 if continuous else (  # type: ignore[attr-defined]
            1.0 if topology.hot_to_small is ValveState.CLOSED else -1.0
        ))
        cold_large.direction = (0.0 if continuous else (  # type: ignore[attr-defined]
            1.0 if topology.cold_to_large is ValveState.CLOSED else -1.0
        ))
        return hot_small, cold_large

    def _record_continuous_diode_events(
        self,
        solution,
        initial_topology: ValveTopology,
        events: list[ValveEvent],
    ) -> None:
        """Record non-terminal zero-pressure crossings from continuous diodes."""

        roots = self._continuous_diode_roots(solution)
        topology = initial_topology
        assert solution.sol is not None
        if solution.t.size > 1:
            probe = min(float(solution.t[1]), float(solution.t[0]) + 1.0e-8)
            values = solution.sol(probe)
            state = ThermodynamicState.from_array(values[:8])
            after_start = self.model.effective_topology(probe, state, topology)
            for name, previous, current in (
                (
                    "hot_to_small",
                    topology.hot_to_small,
                    after_start.hot_to_small,
                ),
                (
                    "cold_to_large",
                    topology.cold_to_large,
                    after_start.cold_to_large,
                ),
            ):
                if previous is not current:
                    events.append(
                        ValveEvent(float(solution.t[0]), name, previous, current)
                    )
            topology = after_start
        for angle, index in roots:
            if math.isclose(angle, float(solution.t[0]), abs_tol=1.0e-10):
                continue
            epsilon = min(1.0e-8, 0.25 * max(solution.t[-1] - solution.t[0], 1.0e-12))
            probe = min(solution.t[-1], angle + epsilon)
            if probe == angle:
                probe = max(solution.t[0], angle - epsilon)
            values = solution.sol(probe)
            state = ThermodynamicState.from_array(values[:8])
            after = self.model.effective_topology(probe, state, topology)
            previous = (
                topology.hot_to_small if index == 0 else topology.cold_to_large
            )
            current = after.hot_to_small if index == 0 else after.cold_to_large
            if current is previous:
                current = (
                    ValveState.OPEN
                    if previous is ValveState.CLOSED
                    else ValveState.CLOSED
                )
            if events and events[-1].valve_name == (
                "hot_to_small" if index == 0 else "cold_to_large"
            ) and math.isclose(events[-1].angle, angle, abs_tol=1.0e-10):
                continue
            name = "hot_to_small" if index == 0 else "cold_to_large"
            events.append(ValveEvent(angle, name, previous, current))
            topology = (
                ValveTopology(current, topology.cold_to_large)
                if index == 0
                else ValveTopology(topology.hot_to_small, current)
            )

    def _continuous_diode_roots(self, solution) -> list[tuple[float, int]]:
        """Locate diode pressure crossings after a topology-free integration."""

        assert solution.sol is not None

        def pressure_difference(angle: float, index: int) -> float:
            state = ThermodynamicState.from_array(solution.sol(angle)[:8])
            pressures = state.pressures(
                self.model.gas, self.model.volumes(angle)
            )
            return float(
                pressures[3] - pressures[0]
                if index == 0
                else pressures[2] - pressures[1]
            )

        roots: list[tuple[float, int]] = []
        for valve_index in (0, 1):
            values = np.array(
                [
                    pressure_difference(float(angle), valve_index)
                    for angle in solution.t
                ]
            )
            for left_index in range(solution.t.size - 1):
                left_angle = float(solution.t[left_index])
                right_angle = float(solution.t[left_index + 1])
                left = values[left_index]
                right = values[left_index + 1]
                if left == 0.0:
                    root = left_angle
                elif left * right < 0.0:
                    root = brentq(
                        lambda angle: pressure_difference(angle, valve_index),
                        left_angle,
                        right_angle,
                    )
                else:
                    continue
                if not roots or not any(
                    index == valve_index
                    and math.isclose(existing, root, abs_tol=1.0e-10)
                    for existing, index in roots[-4:]
                ):
                    roots.append((root, valve_index))
        return sorted(roots)

    def _apply_transitions(
        self,
        theta: float,
        state: ThermodynamicState,
        topology: ValveTopology,
        events: list[ValveEvent],
    ) -> ValveTopology:
        updated = self.model.valve_transitions(theta, state, topology)
        for name, previous, current in (
            ("hot_to_small", topology.hot_to_small, updated.hot_to_small),
            ("cold_to_large", topology.cold_to_large, updated.cold_to_large),
        ):
            if previous is not current:
                events.append(ValveEvent(theta, name, previous, current))
        return updated

    def _apply_root_transitions(
        self,
        theta: float,
        state: ThermodynamicState,
        topology: ValveTopology,
        events: list[ValveEvent],
        triggered_indices: tuple[int, ...],
    ) -> ValveTopology:
        """Apply located roots exactly, then detect any simultaneous passive event."""

        hot_state = topology.hot_to_small
        cold_state = topology.cold_to_large
        if 0 in triggered_indices:
            hot_state = (
                ValveState.OPEN
                if hot_state is ValveState.CLOSED
                else ValveState.CLOSED
            )
        if 1 in triggered_indices:
            cold_state = (
                ValveState.OPEN
                if cold_state is ValveState.CLOSED
                else ValveState.CLOSED
            )
        rooted = ValveTopology(hot_state, cold_state)

        # Evaluate only valves that did not supply the terminal root. This captures
        # a simultaneous threshold condition without undoing a located transition
        # because of floating-point evaluation a few ulps away from its root.
        passively_updated = self.model.valve_transitions(theta, state, rooted)
        if 0 in triggered_indices:
            passively_updated = ValveTopology(
                rooted.hot_to_small, passively_updated.cold_to_large
            )
        if 1 in triggered_indices:
            passively_updated = ValveTopology(
                passively_updated.hot_to_small, rooted.cold_to_large
            )

        for name, previous, current in (
            ("hot_to_small", topology.hot_to_small, passively_updated.hot_to_small),
            ("cold_to_large", topology.cold_to_large, passively_updated.cold_to_large),
        ):
            if previous is not current:
                events.append(ValveEvent(theta, name, previous, current))
        return passively_updated

    def _result(
        self,
        status: IntegrationStatus,
        message: str,
        angle_parts: list[NDArray[np.float64]],
        state_parts: list[NDArray[np.float64]],
        topologies: list[ValveTopology],
        events: list[ValveEvent],
        final_topology: ValveTopology,
        counters: _IntegrationCounters,
    ) -> CycleIntegrationResult:
        if angle_parts:
            angles = np.concatenate(angle_parts)
            augmented = np.concatenate(state_parts, axis=1)
        else:
            angles = np.empty(0, dtype=float)
            augmented = np.empty((11, 0), dtype=float)
        return CycleIntegrationResult(
            status=status,
            message=message,
            angles=angles,
            states=augmented[:8],
            cold_heat=augmented[8],
            hot_heat=augmented[9],
            gas_work=augmented[10],
            topologies=tuple(topologies),
            events=tuple(events),
            final_topology=final_topology,
            statistics=IntegrationStatistics(
                elapsed_seconds=time.perf_counter() - counters.start_time,
                solve_segments=counters.solve_segments,
                right_hand_side_evaluations=counters.right_hand_side_evaluations,
                jacobian_evaluations=counters.jacobian_evaluations,
                linear_decompositions=counters.linear_decompositions,
                event_function_evaluations=counters.event_function_evaluations,
                stored_samples=angles.size,
                segments=tuple(counters.segments or ()),
            ),
        )
