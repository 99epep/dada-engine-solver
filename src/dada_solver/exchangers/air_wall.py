"""Conservative lumped wall storage with a finite-capacity external air stream.

The air passage is quasi-steady against a uniform wall. Wall conductances and
capacity are explicit inputs; Doty's overall UA does not identify their split.
"""
from dataclasses import dataclass
import math
from dada_solver import numerical_primitives as numeric
import time

import numpy as np
from scipy.integrate import solve_ivp

from dada_solver.dynamics import ThermodynamicModel
from dada_solver.exchangers.base import LumpedWallThermalModel
from dada_solver.state import ThermodynamicState
from dada_solver.dynamics import ValveTopology
from dada_solver.valves import ValveState
from dada_solver.integration import IntegrationInterrupted


@dataclass(frozen=True)
class AirWallExchanger:
    gas_wall_conductance_w_k: float
    air_wall_conductance_w_k: float
    wall_capacity_j_k: float
    air_mass_flow_kg_s: float
    air_cp_j_kg_k: float
    air_inlet_temperature_k: float
    gas_film: object | None = None

    def __post_init__(self):
        for value in (self.wall_capacity_j_k, self.air_cp_j_kg_k, self.air_inlet_temperature_k):
            if not math.isfinite(value) or value <= 0:
                raise ValueError('Wall capacity, air cp and inlet temperature must be positive.')
        for value in (self.gas_wall_conductance_w_k, self.air_wall_conductance_w_k, self.air_mass_flow_kg_s):
            if not math.isfinite(value) or value < 0:
                raise ValueError('Conductances and air flow must be finite and nonnegative.')
        capacity_rate = self.air_mass_flow_kg_s*self.air_cp_j_kg_k
        effective = (capacity_rate * -math.expm1(-self.air_wall_conductance_w_k/capacity_rate)
                     if capacity_rate > 0 else 0.0)
        object.__setattr__(self, '_air_capacity_rate', capacity_rate)
        object.__setattr__(self, '_effective_air_conductance', effective)

    @property
    def requires_flow_context(self):
        return self.gas_film is not None

    def rates(self, gas_temperature_k, wall_energy_j, *, context=None):
        wall_temperature = wall_energy_j / self.wall_capacity_j_k
        if any(not math.isfinite(v) or v <= 0 for v in (wall_temperature, gas_temperature_k)):
            raise ValueError('Gas and wall temperatures must be positive.')
        capacity_rate = self._air_capacity_rate
        effective = self._effective_air_conductance
        conductance = self.gas_wall_conductance_w_k
        if self.gas_film is not None:
            if context is None: raise ValueError('Variable gas film requires instantaneous flow context.')
            conductance, _ = self.gas_film.evaluate(gas_temperature_k,wall_temperature,context)
        air_heat,gas_heat,wall_rate = numeric.wall_heat_rates(effective,self.air_inlet_temperature_k,
            wall_temperature,conductance,gas_temperature_k)
        return dict(gas_heat_w=gas_heat, air_heat_w=air_heat,
                    wall_energy_rate_w=wall_rate,
                    air_outlet_temperature_k=(self.air_inlet_temperature_k-air_heat/capacity_rate
                                              if capacity_rate > 0 else None))


@dataclass(frozen=True)
class AirWallMotor:
    """Opt-in 10-state motor wrapper; existing eight-state solver is unchanged.

    State: eight conservative gas values, then wall energies for H_i and H_o.
    Integration appends air heats, gas heats and gas work (five quadratures).
    Only continuous ideal diodes are supported; hysteretic memory is not hidden.
    """
    model: ThermodynamicModel
    heat_in: LumpedWallThermalModel
    heat_out: LumpedWallThermalModel

    def __post_init__(self):
        if not self.model.continuous_ideal_diodes or self.model.study_crank_direction != -1:
            raise ValueError('Air-wall wrapper requires motor operation with continuous ideal diodes.')

    def derivative(self, angle, values):
        gas = ThermodynamicState.from_array(np.asarray(values[:8]))
        point = self.model.instantaneous_point(
            angle, gas, ValveTopology(ValveState.CLOSED, ValveState.CLOSED))
        incoming, outgoing = self._thermal_rates_from_point(point, values)
        rates = self.model.assemble_rates(point, incoming['gas_heat_w'], outgoing['gas_heat_w'])
        return np.r_[rates.state_derivative,
                     incoming['wall_energy_rate_w'], outgoing['wall_energy_rate_w'],
                     incoming['air_heat_w'], outgoing['air_heat_w'],
                     incoming['gas_heat_w'], outgoing['gas_heat_w'], rates.gas_work_rate] / self.model.angular_speed

    def flow_contexts(self, angle, values):
        """Geometry-independent pressure/flow state for contextual wall closures."""
        gas = ThermodynamicState.from_array(np.asarray(values[:8]))
        point = self.model.instantaneous_point(
            angle, gas, ValveTopology(ValveState.CLOSED, ValveState.CLOSED))
        return self._flow_contexts_from_point(point)

    def _flow_contexts_from_point(self, point):
        pressures, flows = point.pressures, point.flows
        frequency = self.model.angular_speed/(2*math.pi)
        contexts = (
            dict(frequency_hz=frequency,passages=((flows.small_to_cold,pressures[0],pressures[2]),
                                                 (flows.cold_to_large,pressures[2],pressures[1]))),
            dict(frequency_hz=frequency,passages=((flows.large_to_hot,pressures[1],pressures[3]),
                                                 (flows.hot_to_small,pressures[3],pressures[0]))))
        for index,context in zip((2,3),contexts):
            context['port_pressure_ratios'] = tuple(max(p1,p2)/min(p1,p2) for _,p1,p2 in context['passages'])
            # A closed valve supports a pressure jump; it is not an axial tube
            # pressure gradient at zero flow. Retain the port ratio separately.
            context['passages'] = tuple((flow,p1,p2) if flow else
                (flow,pressures[index],pressures[index]) for flow,p1,p2 in context['passages'])
        return contexts

    def thermal_rates(self, angle, values):
        """Single diagnostic/integration access path; never substitutes static UA."""
        gas = ThermodynamicState.from_array(np.asarray(values[:8]))
        contextual = any(getattr(x,'requires_flow_context',False) for x in (self.heat_in,self.heat_out))
        if contextual:
            point = self.model.instantaneous_point(
                angle, gas, ValveTopology(ValveState.CLOSED, ValveState.CLOSED))
            return self._thermal_rates_from_point(point, values)
        return self._thermal_rates(gas.temperatures(self.model.gas), values, (None, None))

    def _thermal_rates_from_point(self, point, values):
        contextual = any(getattr(x,'requires_flow_context',False) for x in (self.heat_in,self.heat_out))
        contexts = self._flow_contexts_from_point(point) if contextual else (None, None)
        return self._thermal_rates(point.temperatures, values, contexts)

    def _thermal_rates(self, temperatures, values, contexts):
        return tuple(exchanger.rates(temperatures[index],values[index+6],context=context)
                     if getattr(exchanger,'requires_flow_context',False)
                     else exchanger.rates(temperatures[index],values[index+6])
                     for index,exchanger,context in zip((2,3),(self.heat_in,self.heat_out),contexts))

    def integrate_cycle(self, state, *, integration_method='LSODA', rtol=1e-7, atol=1e-10,
                        maximum_step_angle=math.pi/360, progress_callback=None,
                        progress_interval_seconds=10.0, statistics_callback=None,
                        measure_rhs_time=False, rhs=None):
        initial = np.asarray(state, dtype=float)
        if initial.shape != (10,):
            raise ValueError('Expected eight gas states and two wall energies.')
        evaluate_rhs = self.derivative if rhs is None else rhs
        preflight_started = time.perf_counter()
        try:
            evaluate_rhs(0, np.r_[initial, np.zeros(5)])
        except Exception as error:
            if statistics_callback is not None:
                statistics_callback(dict(phase='integration_preflight', status='failed',
                    elapsed_seconds=time.perf_counter()-preflight_started,
                    exception=type(error).__name__, message=str(error)))
            raise
        if statistics_callback is not None:
            statistics_callback(dict(phase='integration_preflight', status='completed',
                elapsed_seconds=time.perf_counter()-preflight_started, rhs_calls=1))
        if not math.isfinite(maximum_step_angle) or maximum_step_angle <= 0:
            raise ValueError('Maximum angular step must be positive.')
        provider = getattr(self.model.kinematics, 'breakpoint_angles', lambda: ())
        edges = sorted({0.0, 2*math.pi, *(float(x) for x in provider() if 0 < x < 2*math.pi)})
        values = np.r_[initial, np.zeros(5)]
        angles, histories = [], []
        for lower, upper in zip(edges[:-1], edges[1:]):
            last_progress = time.perf_counter()
            evaluations = 0
            rhs_seconds = 0.0
            started = time.perf_counter()
            def derivative(angle, state):
                nonlocal last_progress, evaluations, rhs_seconds
                evaluations += 1
                now = time.perf_counter()
                if progress_callback is not None and now-last_progress >= progress_interval_seconds:
                    progress_callback(dict(phase='running', start_angle=lower,
                        current_angle=float(angle), target_angle=upper,
                        elapsed_seconds=now-last_progress,
                        right_hand_side_evaluations=evaluations))
                    last_progress = now
                if not measure_rhs_time:
                    return evaluate_rhs(angle, state)
                before = time.perf_counter()
                try:
                    return evaluate_rhs(angle, state)
                finally:
                    rhs_seconds += time.perf_counter()-before

            def report(solution=None, error=None):
                if statistics_callback is None:
                    return
                steps = np.diff(solution.t) if solution is not None else np.array([])
                statistics_callback(dict(
                    phase='solver_segment', method=integration_method,
                    start_angle=lower, target_angle=upper,
                    ends_at_kinematic_breakpoint=upper in edges[1:-1],
                    status=('interrupted' if isinstance(error, IntegrationInterrupted) else
                            'failed' if error is not None or not solution.success else 'completed'),
                    message=str(error) if error is not None else solution.message,
                    nfev=int(solution.nfev) if solution is not None else None,
                    njev=int(solution.njev) if solution is not None else None,
                    nlu=int(solution.nlu) if solution is not None else None,
                    rhs_calls=evaluations,
                    accepted_samples=len(solution.t) if solution is not None else None,
                    accepted_steps=len(steps) if solution is not None else None,
                    minimum_step_angle=float(np.min(steps)) if len(steps) else None,
                    median_step_angle=float(np.median(steps)) if len(steps) else None,
                    maximum_step_angle=float(np.max(steps)) if len(steps) else None,
                    integration_elapsed_seconds=time.perf_counter()-started,
                    rhs_elapsed_seconds=rhs_seconds if measure_rhs_time else None,
                ))

            try:
                # Fast segments may all finish before the timed RHS callback.
                # Always check cancellation/deadlines at physical boundaries.
                if progress_callback is not None:
                    progress_callback(dict(phase='running', start_angle=lower,
                        current_angle=lower, target_angle=upper,
                        elapsed_seconds=0.0, right_hand_side_evaluations=0))
                solution = solve_ivp(derivative, (lower, upper), values, method=integration_method,
                                     rtol=rtol, atol=atol, max_step=maximum_step_angle)
            except Exception as error:
                report(error=error)
                raise
            report(solution)
            if not solution.success:
                raise RuntimeError(solution.message)
            angles.extend(solution.t)
            histories.extend(solution.y.T)
            values = solution.y[:, -1]
        return np.asarray(angles), np.asarray(histories).T
