"""Reusable periodic evaluation for the existing ten-state air-wall motor."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
import math
import time
from types import SimpleNamespace
import numpy as np

from dada_solver.dynamics import ValveTopology
from dada_solver.exchangers.air_wall import AirWallMotor
from dada_solver.integration import IntegrationInterrupted
from dada_solver.wall_backend import WallBackendSettings, WallRHS
from dada_solver.performance import ConservationReport, CyclePerformance, OperatingMode
from dada_solver.state import ThermodynamicState
from dada_solver.valves import ValveState


@dataclass(frozen=True)
class WallCycleNumericalSettings:
    """Numerical controls specific to the ten-state wall-cycle calculation."""
    integration_method: str = 'LSODA'
    integration_relative_tolerance: float = 1e-8
    integration_absolute_tolerances: tuple[float, ...] = tuple([1e-13, 1e-8]*4 + [1e-7]*2 + [1e-8]*5)
    maximum_step_angle_radians: float = math.pi/360
    periodic_relative_tolerance: float = 1e-6
    periodic_absolute_tolerance: float = 1e-12
    progress_interval_seconds: float = 1.0
    accelerate_walls: bool = False

    def __post_init__(self):
        if self.integration_method not in {'RK45','Radau','BDF','LSODA'}:
            raise ValueError('Unsupported wall integration method.')
        if len(self.integration_absolute_tolerances) != 15:
            raise ValueError('Wall integration requires 15 absolute tolerances.')
        values = (self.integration_relative_tolerance, *self.integration_absolute_tolerances,
            self.maximum_step_angle_radians, self.periodic_relative_tolerance,
            self.periodic_absolute_tolerance, self.progress_interval_seconds)
        if any(not math.isfinite(x) or x <= 0 for x in values):
            raise ValueError('Wall numerical tolerances, step and progress interval must be positive.')


@dataclass(frozen=True)
class WallCycleResult:
    status: str
    message: str
    history: tuple[dict, ...]
    angles: np.ndarray | None
    trajectory: np.ndarray | None
    last_complete_state: np.ndarray | None
    solver_statistics: tuple[dict, ...] = ()
    backend_statistics: dict = field(default_factory=dict)

    anderson_statistics: dict = field(default_factory=dict)

    @property
    def converged(self):
        return self.status == 'converged'


def solve_periodic_wall_motor(wrapper: AirWallMotor, initial_state, *, maximum_cycles,
        progress_callback=None, settings=WallCycleNumericalSettings(), cycle_callback=None,
        statistics_callback=None, measure_rhs_time=False,
        adaptive_acceleration=None, backend=WallBackendSettings(), exact_kinematics_cache=True,
        anderson_acceleration=None) -> WallCycleResult:
    """Repeat complete physical cycles at the original periodic tolerances.

    Optional adaptive wall guesses replace the fixed ten-cycle schedule. Failed
    trials consume the attempt budget and are reported, then roll back to the
    retained physical endpoint. Only retained cycles notify cycle_callback.
    Campaign defaults and their serialized numerical settings remain unchanged.
    Experimental Anderson mixing instead uses recent full physical map pairs;
    it suppresses fixed wall extrapolation and cannot certify convergence itself.
    """
    from dada_solver.exchangers.wall_iteration import (
        extrapolate_wall_energies, adaptive_wall_proposal, AdaptiveWallAccelerationSettings,
        AndersonAccelerationSettings, PeriodicMapPair, anderson_proposal, mass_compatible, MASS_INDICES, _physical_state)
    if adaptive_acceleration is not None and not isinstance(adaptive_acceleration, AdaptiveWallAccelerationSettings):
        raise TypeError('Expected AdaptiveWallAccelerationSettings or None.')
    if anderson_acceleration is not None and not isinstance(anderson_acceleration,AndersonAccelerationSettings):
        raise TypeError('Expected AndersonAccelerationSettings or None.')
    if adaptive_acceleration is not None and anderson_acceleration is not None:
        raise ValueError('Anderson and adaptive wall acceleration are mutually exclusive.')
    if not isinstance(backend,WallBackendSettings):
        raise TypeError('Expected WallBackendSettings.')
    from dada_solver.kinematics_cache import prepare_exact_kinematics, ExactAngleKinematics
    cache = None
    if exact_kinematics_cache and isinstance(wrapper, AirWallMotor):
        prepared = prepare_exact_kinematics(wrapper.model.kinematics)
        if isinstance(prepared, ExactAngleKinematics):
            cache = prepared
            wrapper = replace(wrapper, model=replace(wrapper.model, kinematics=prepared))
    state = np.asarray(initial_state, dtype=float)
    history, wall_history = [], []
    statistics = []
    last_angles = last_trajectory = None
    pending = None
    anderson_pending = None
    pairs = []
    anderson_counts = (dict(attempted_cycles=0,retained_cycles=0,proposal_count=0,
        accepted_trial_count=0,rejected_trial_count=0,rollback_count=0,domain_failure_count=0)
        if anderson_acceleration is not None else {})
    charge = float(_physical_state(state)[MASS_INDICES].sum()) if anderson_acceleration is not None else None
    recovery_remaining = 0
    capacities = (np.array([wrapper.heat_in.wall_capacity_j_k, wrapper.heat_out.wall_capacity_j_k])
                  if adaptive_acceleration is not None else None)
    rhs = None
    backend_options = {}
    def backend_snapshot():
        snapshot = rhs.snapshot() if rhs is not None else dict(requested_backend=backend.name,
            actual_backend='python' if backend.name=='python' else 'not_started')
        if cache is not None: snapshot['exact_kinematics_cache'] = cache.snapshot()
        return snapshot
    atol = np.asarray(settings.integration_absolute_tolerances)
    for cycle in range(1, maximum_cycles+1):
        input_state = state.copy() if anderson_acceleration is not None else None
        if anderson_acceleration is not None: anderson_counts['attempted_cycles'] += 1
        def record_segment(record):
            record = dict(record, cycle=cycle)
            statistics.append(record)
            if statistics_callback is not None:
                statistics_callback(record)

        def restore_anchor(reason, **details):
            nonlocal state, pending, recovery_remaining, wall_history, last_angles, last_trajectory
            state = pending['state'].copy()
            last_angles, last_trajectory = pending['angles'], pending['trajectory']
            record_segment(dict(phase='wall_acceleration_rollback', status='rejected',
                reason=reason, retained_cycle=pending['cycle'],
                retained_normalized_error=pending['error'], **details))
            pending = None
            wall_history = []
            recovery_remaining = adaptive_acceleration.recovery_cycles

        def rollback_anderson(reason, **details):
            nonlocal state, anderson_pending, pairs
            state=anderson_pending['state'].copy()
            anderson_counts['rejected_trial_count']+=1
            anderson_counts['rollback_count']+=1
            record_segment(dict(phase='anderson_trial',status='rolled_back',anderson_trial=True,
                rollback_reason=reason,previous_normalized_error=anderson_pending['error'],**details))
            anderson_pending=None;pairs=[]

        try:
            if rhs is None and (backend.name!='python' or backend.profile):
                # Keep Python's preflight-before-first-progress ordering.
                rhs = WallRHS(wrapper,backend)
                backend_options['rhs'] = rhs
            angles, trajectory = wrapper.integrate_cycle(state,
                integration_method=settings.integration_method,
                rtol=settings.integration_relative_tolerance, atol=atol,
                maximum_step_angle=settings.maximum_step_angle_radians,
                progress_callback=progress_callback,
                progress_interval_seconds=settings.progress_interval_seconds,
                statistics_callback=record_segment, measure_rhs_time=measure_rhs_time, **backend_options)
            if anderson_acceleration is not None:
                if not np.all(np.isfinite(trajectory)) or np.any(trajectory[:10]<=0):
                    raise ValueError('Physical cycle returned nonfinite or nonpositive states.')
                if not mass_compatible(trajectory[:10,-1],charge):
                    raise ValueError('Physical cycle changed the fixed gas charge.')
        except IntegrationInterrupted as error:
            if anderson_pending is not None:
                record_segment(dict(phase='anderson_trial',status='interrupted',anderson_trial=True,
                    previous_normalized_error=anderson_pending['error']))
            return WallCycleResult('interrupted', str(error), tuple(history),
                last_angles, last_trajectory,
                last_trajectory[:10, -1].copy() if last_trajectory is not None else None, tuple(statistics), backend_snapshot(), dict(anderson_counts))
        except (ValueError, RuntimeError, ArithmeticError) as error:
            if anderson_pending is not None:
                from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
                if isinstance(error,MicrotubeDomainError):anderson_counts['domain_failure_count']+=1
                rollback_anderson('integration_failed',exception=type(error).__name__,message=str(error))
                continue
            if pending is None:
                raise
            restore_anchor('accelerated_integration_failed', exception=type(error).__name__, message=str(error))
            continue
        finally:
            if rhs is not None and statistics_callback is not None:
                statistics_callback(dict(phase='rhs_backend',cycle=cycle,**backend_snapshot()))
        end = trajectory[:10, -1]
        normalized = np.abs(end-state) / (
            settings.periodic_absolute_tolerance + settings.periodic_relative_tolerance*np.maximum(np.abs(end), np.abs(state)))
        error = float(np.max(normalized))
        history.append(dict(cycle=cycle, normalized_state_error=error))
        if anderson_acceleration is not None:
            history[-1].update(anderson_trial=anderson_pending is not None,retained=False)
            if anderson_pending is not None:
                history[-1]['trial_normalized_error']=error
                if error>1 and error>anderson_acceleration.maximum_residual_growth*anderson_pending['error']:
                    history[-1].update(anderson_rolled_back=True,rollback_reason='residual_growth')
                    rollback_anderson('residual_growth',trial_normalized_error=error)
                    continue
                anderson_counts['accepted_trial_count']+=1
                history[-1]['anderson_accepted']=True
                record_segment(dict(phase='anderson_trial',status='accepted',anderson_trial=True,trial_normalized_error=error))
                anderson_pending=None
            history[-1]['retained']=True
            anderson_counts['retained_cycles']+=1
            pairs.append(PeriodicMapPair(input_state,end,error));pairs=pairs[-anderson_acceleration.memory:]
        if adaptive_acceleration is not None:
            history[-1].update(normalized_gas_error=float(np.max(normalized[:8])),
                               normalized_wall_error=float(np.max(normalized[8:10])))
        if pending is not None:
            if 'trial_error' in pending:
                safe = (error <= 1 or (error <= adaptive_acceleration.maximum_residual_growth*pending['error']
                        and error < .5*pending['trial_error']))
                if not safe:
                    history[-1]['wall_acceleration_rolled_back'] = True
                    restore_anchor('transient_failed_to_contract', trial_normalized_error=error)
                    continue
                wall_history = [pending['trial_wall']]
            else:
                allowed = adaptive_acceleration.maximum_residual_growth*max(pending['error'],pending['normalized_jump'])
                if error > allowed:
                    history[-1]['wall_acceleration_rolled_back'] = True
                    restore_anchor('normalized_residual_growth', trial_normalized_error=error,
                                   allowed_normalized_error=allowed)
                    continue
                if error > pending['error'] and error > 1:
                    # A wall-only jump leaves the gas to adjust physically. A
                    # plausible transient is provisional until the next cycle
                    # demonstrably contracts; it is not a convergence result.
                    history[-1]['wall_acceleration_pending_validation'] = True
                    pending['trial_error'] = error
                    pending['trial_wall'] = end[8:10].copy()
                    state = end.copy()
                    continue
                wall_history = []
            history[-1]['wall_acceleration_accepted'] = True
            pending = None
        last_angles, last_trajectory = angles, trajectory
        state = end
        if cycle_callback is not None:
            cycle_callback(cycle, end.copy(), error, history[-1])
        if error <= 1:
            return WallCycleResult('converged', 'Periodic steady state converged.',
                tuple(history), angles, trajectory, end.copy(), tuple(statistics), backend_snapshot(), dict(anderson_counts))
        if anderson_acceleration is not None:
            if cycle<maximum_cycles:
                started=time.perf_counter()
                proposal=anderson_proposal(pairs,settings=anderson_acceleration)
                record_segment(dict(phase='anderson_decision',**proposal.statistics,
                    previous_normalized_error=error,elapsed_seconds=time.perf_counter()-started))
                if proposal.state is not None:
                    anderson_counts['proposal_count']+=1
                    anderson_pending=dict(state=end.copy(),error=error)
                    state=proposal.state.copy()
            continue
        wall_history.append(end[8:10].copy())
        if adaptive_acceleration is not None:
            wall_history = wall_history[-4:]
            recovery_remaining = max(0, recovery_remaining-1)
            if len(wall_history)==4 and recovery_remaining==0 and cycle < maximum_cycles:
                decision_started = time.perf_counter()
                proposal = adaptive_wall_proposal(wall_history, capacities, settings=adaptive_acceleration)
                record_segment(dict(phase='wall_acceleration_decision',
                    status='proposed' if proposal is not None else 'skipped',
                    elapsed_seconds=time.perf_counter()-decision_started))
                if proposal is not None:
                    pending = dict(state=end.copy(), error=error, cycle=cycle,
                        angles=angles, trajectory=trajectory)
                    state = end.copy()
                    state[8:10] = proposal.wall_energies
                    pending['normalized_jump'] = float(np.max(np.abs(state[8:10]-end[8:10]) / (
                        settings.periodic_absolute_tolerance + settings.periodic_relative_tolerance*
                        np.maximum(np.abs(state[8:10]),np.abs(end[8:10])))))
                    history[-1]['wall_initial_guess_extrapolated'] = True
                    history[-1]['adaptive_wall_proposal'] = asdict(proposal)
                    wall_history = []
        elif settings.accelerate_walls and cycle % 10 == 0 and cycle < maximum_cycles and len(wall_history) >= 3:
            proposed = extrapolate_wall_energies(*wall_history[-3:], np.array([
                wrapper.heat_in.wall_capacity_j_k, wrapper.heat_out.wall_capacity_j_k]))
            if proposed is not None:
                state = state.copy(); state[8:10] = proposed
                history[-1]['wall_initial_guess_extrapolated'] = True
            wall_history = []
    return WallCycleResult('maximum_cycles',
        'Maximum cycle count reached before periodic convergence.', tuple(history),
        last_angles, last_trajectory,
        last_trajectory[:10, -1].copy() if last_trajectory is not None else None, tuple(statistics), backend_snapshot(), dict(anderson_counts))


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
