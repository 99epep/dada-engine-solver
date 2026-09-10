"""Conservative lumped wall storage with a finite-capacity external air stream.

The air passage is quasi-steady against a uniform wall. Wall conductances and
capacity are explicit inputs; Doty's overall UA does not identify their split.
"""
from dataclasses import dataclass, replace
import math
import time

import numpy as np
from scipy.integrate import solve_ivp

from dada_solver.dynamics import ThermodynamicModel
from dada_solver.heat_transfer import PrescribedHeatRate
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

    def __post_init__(self):
        for value in (self.wall_capacity_j_k, self.air_cp_j_kg_k, self.air_inlet_temperature_k):
            if not math.isfinite(value) or value <= 0:
                raise ValueError('Wall capacity, air cp and inlet temperature must be positive.')
        for value in (self.gas_wall_conductance_w_k, self.air_wall_conductance_w_k, self.air_mass_flow_kg_s):
            if not math.isfinite(value) or value < 0:
                raise ValueError('Conductances and air flow must be finite and nonnegative.')

    def rates(self, gas_temperature_k, wall_energy_j):
        wall_temperature = wall_energy_j / self.wall_capacity_j_k
        if any(not math.isfinite(v) or v <= 0 for v in (wall_temperature, gas_temperature_k)):
            raise ValueError('Gas and wall temperatures must be positive.')
        capacity_rate = self.air_mass_flow_kg_s*self.air_cp_j_kg_k
        effective = (capacity_rate * -math.expm1(-self.air_wall_conductance_w_k/capacity_rate)
                     if capacity_rate > 0 else 0.0)
        air_heat = effective*(self.air_inlet_temperature_k-wall_temperature)
        gas_heat = self.gas_wall_conductance_w_k*(wall_temperature-gas_temperature_k)
        return dict(gas_heat_w=gas_heat, air_heat_w=air_heat,
                    wall_energy_rate_w=air_heat-gas_heat,
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
        temperatures = gas.temperatures(self.model.gas)
        incoming = self.heat_in.rates(temperatures[2], values[8])
        outgoing = self.heat_out.rates(temperatures[3], values[9])
        # Reuse all existing conservative transport and passive-valve equations.
        instantaneous = replace(self.model,
            cold_heat_transfer=PrescribedHeatRate(incoming['gas_heat_w']),
            hot_heat_transfer=PrescribedHeatRate(outgoing['gas_heat_w']))
        rates = instantaneous.evaluate(angle, gas, ValveTopology(ValveState.CLOSED, ValveState.CLOSED))
        return np.r_[rates.state_derivative,
                     incoming['wall_energy_rate_w'], outgoing['wall_energy_rate_w'],
                     incoming['air_heat_w'], outgoing['air_heat_w'],
                     incoming['gas_heat_w'], outgoing['gas_heat_w'], rates.gas_work_rate] / self.model.angular_speed

    def integrate_cycle(self, state, *, rtol=1e-7, atol=1e-10,
                        maximum_step_angle=math.pi/360, progress_callback=None,
                        progress_interval_seconds=10.0):
        initial = np.asarray(state, dtype=float)
        if initial.shape != (10,):
            raise ValueError('Expected eight gas states and two wall energies.')
        self.derivative(0, np.r_[initial, np.zeros(5)])
        if not math.isfinite(maximum_step_angle) or maximum_step_angle <= 0:
            raise ValueError('Maximum angular step must be positive.')
        provider = getattr(self.model.kinematics, 'breakpoint_angles', lambda: ())
        edges = sorted({0.0, 2*math.pi, *(float(x) for x in provider() if 0 < x < 2*math.pi)})
        values = np.r_[initial, np.zeros(5)]
        angles, histories = [], []
        for lower, upper in zip(edges[:-1], edges[1:]):
            last_progress = time.perf_counter()
            evaluations = 0
            def derivative(angle, state):
                nonlocal last_progress, evaluations
                evaluations += 1
                now = time.perf_counter()
                if progress_callback is not None and now-last_progress >= progress_interval_seconds:
                    progress_callback(dict(phase='running', start_angle=lower,
                        current_angle=float(angle), target_angle=upper,
                        elapsed_seconds=now-last_progress,
                        right_hand_side_evaluations=evaluations))
                    last_progress = now
                return self.derivative(angle, state)
            solution = solve_ivp(derivative, (lower, upper), values, method='LSODA',
                                 rtol=rtol, atol=atol, max_step=maximum_step_angle)
            if not solution.success:
                raise RuntimeError(solution.message)
            angles.extend(solution.t)
            histories.extend(solution.y.T)
            values = solution.y[:, -1]
        return np.asarray(angles), np.asarray(histories).T
