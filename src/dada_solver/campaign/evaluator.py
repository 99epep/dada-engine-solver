"""Campaign evaluation for reservoir and lumped-wall machine models."""
from dataclasses import asdict, dataclass, replace
from enum import Enum
import math
import time
from types import SimpleNamespace
import numpy as np

from dada_solver.campaign.adapters import PreflightRejection
from dada_solver.exchangers.air_wall import AirWallMotor
from dada_solver.exchangers.wall_cycle import (WallDiagnosticCycle,
    convergence_summary, solve_periodic_wall_motor, wall_cycle_performance)
from dada_solver.factory import build_initial_state, initial_valve_topology
from dada_solver.integration import IntegrationInterrupted
from dada_solver.kinematics import KinematicConstraintViolation
from dada_solver.results import extract_cycle_diagnostics
from dada_solver.sizing.evaluator import evaluate_configuration, EvaluationStatus
from dada_solver.state import ThermodynamicState, UniformCharge
from dada_solver.validity import assess_cycle_validity, ValidityVerdict


@dataclass(frozen=True)
class EvaluationControl:
    deadline: float | None = None
    clock: object = time.monotonic
    previous_records: tuple = ()
    def check(self, _progress=None):
        if self.deadline is not None and self.clock() >= self.deadline:
            raise IntegrationInterrupted('Campaign wall-clock deadline reached.')


def json_values(value):
    if isinstance(value, Enum): return value.value
    if isinstance(value, dict): return {str(k): json_values(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)): return [json_values(v) for v in value]
    if isinstance(value, np.generic): return value.item()
    return value


def rejected(status, reason, diagnostics=None):
    return dict(status=status, reason=reason, integrated=False, converged=False,
        objective=None, constraints=[], metrics={}, derived={}, periodic_cycle_count=0,
        periodic_convergence=convergence_summary(()), warm_start_source=None,
        warm_start_normalized_distance=None, warm_start_source_status=None,
        final_periodic_state=None, initial_guess_state=None, preflight_diagnostics=diagnostics)


def select_warm_start(records, normalized, layout, family, direction):
    choices = []
    for record in records:
        state = record.get('initial_guess_state') or record.get('final_periodic_state')
        compatible = (state and state.get('state_layout') == layout and
            state.get('evaluator_family') == family and state.get('operating_direction') == direction and
            record.get('status') not in {'integration_failure','invalid_parameterization','invalid_kinematics','invalid_exchanger'})
        if compatible:
            distance = float(np.linalg.norm(np.asarray(normalized)-np.asarray(record['normalized'])))
            choices.append((not record.get('converged', False), distance, record['candidate_id'], record))
    if not choices: return None, None
    selected = min(choices)[-1]
    return selected, float(np.linalg.norm(np.asarray(normalized)-np.asarray(selected['normalized'])))


def _state_record(values, layout, family, direction, *, periodic, wall_capacities=None):
    return dict(state_layout=layout, evaluator_family=family, operating_direction=direction,
        values=np.asarray(values).tolist(), total_mass_kg=float(np.asarray(values)[:8:2].sum()),
        wall_capacities_j_k=wall_capacities, periodic_solution=periodic,
        reusable_only_as_initial_guess=not periodic)


def rescale_wall_state(values, target_gas_mass, old_wall_capacities, new_wall_capacities):
    """Preserve gas specific energies and wall temperatures for a new design."""
    state = np.asarray(values, dtype=float).copy()
    state[:8] *= target_gas_mass/state[:8:2].sum()
    state[8:10] *= np.asarray(new_wall_capacities)/np.asarray(old_wall_capacities)
    return state


def finalize_microtube_validity(validity, maximum_reynolds, maximum_mach, mach_limit,
                               *, requires_laminar=True, domain_failures=()):
    """Replace generic Mach unavailability and recompute the global verdict."""
    failed = list(validity.failed_criteria)
    if requires_laminar and maximum_reynolds >= 2300: failed.append('microtube_internal_laminar_reynolds')
    failed.extend(domain_failures)
    if maximum_mach > mach_limit: failed.append('microtube_mach_number')
    unavailable = tuple(x for x in validity.unavailable_criteria if x != 'mach_number')
    verdict = (ValidityVerdict.INVALID if failed else
        ValidityVerdict.INDETERMINATE if unavailable else ValidityVerdict.VALID)
    return replace(validity, verdict=verdict, maximum_mach_number=maximum_mach,
        mach_number_status='available', failed_criteria=tuple(failed),
        unavailable_criteria=unavailable)


def _gas_trend(history):
    return convergence_summary([dict(cycle=x.cycle_number,
        normalized_state_error=x.normalized_state_error) for x in history])


class MachineEvaluator:
    def __init__(self, definition): self.definition = definition
    def evaluate(self, candidate): return self.evaluate_with_control(candidate, EvaluationControl())

    def _unavailable_constraints(self):
        return [dict(name=c.name, margin=None, satisfied=False, available=False)
                for c in self.definition.constraints]

    def evaluate_with_control(self, candidate, control):
        try:
            physical = dict(self.definition.fixed_parameters)
            physical.update(candidate.payload['physical'])
            design = self.definition.adapter.build(physical)
        except PreflightRejection as error: return rejected(error.status, str(error), error.diagnostics)
        try: built = design.build()
        except KinematicConstraintViolation as error:
            return rejected('invalid_kinematics', str(error), [asdict(d) for d in error.diagnostics])
        except (ValueError, ArithmeticError) as error: return rejected('invalid_exchanger', str(error))
        model = built.model if isinstance(built, AirWallMotor) else built
        derived = dict(small_volume_limits=asdict(model.machine_volumes.small_cylinder),
            large_volume_limits=asdict(model.machine_volumes.large_cylinder),
            heat_in_gas_volume_m3=model.machine_volumes.cold_heat_exchanger,
            heat_out_gas_volume_m3=model.machine_volumes.hot_heat_exchanger,
            small_physical_stroke_m=model.kinematics.small_physical_stroke,
            large_physical_stroke_m=model.kinematics.large_physical_stroke)
        direction = 'motor' if design.configuration.motor_operation else 'receiver'
        if isinstance(built, AirWallMotor):
            return self._wall(candidate, design, built, derived, direction, control)
        return self._reservoir(candidate, design, built, derived, direction, control)

    def _reservoir(self, candidate, design, model, derived, direction, control):
        source, distance = select_warm_start(control.previous_records, candidate.payload['normalized'],
            'four_gas_volumes_m_U', 'reservoir_8_state', direction)
        initial = topology = None
        if source:
            saved = source.get('initial_guess_state') or source['final_periodic_state']
            initial = ThermodynamicState.from_array(saved['values'])
            topology_data = saved.get('topology')
            if topology_data:
                from dada_solver.valves import ValveState
                from dada_solver.dynamics import ValveTopology
                topology = ValveTopology(ValveState(topology_data['hot_to_small']), ValveState(topology_data['cold_to_large']))
        try:
            evaluation = evaluate_configuration(design.configuration, initial, topology, control.check, model=model)
        except IntegrationInterrupted as error:
            result = rejected('budget_exhausted', str(error)); result.update(integrated=True, derived=derived,
                constraints=self._unavailable_constraints())
            return result
        except (ValueError, RuntimeError, ArithmeticError) as error:
            result = rejected('integration_failure', f'{type(error).__name__}: {error}'); result.update(integrated=True, derived=derived)
            return result
        periodic_status = getattr(evaluation.periodic, 'status', None)
        if getattr(periodic_status, 'value', None) == 'interrupted':
            cycle = evaluation.periodic.final_cycle
            guess = None
            if cycle is not None and cycle.completed:
                guess = _state_record(cycle.states[:, -1], 'four_gas_volumes_m_U', 'reservoir_8_state', direction, periodic=False)
                guess['topology'] = json_values(asdict(cycle.final_topology))
            result = rejected('budget_exhausted', evaluation.periodic.message)
            result.update(integrated=True, derived=derived, periodic_cycle_count=len(evaluation.periodic.history),
                periodic_convergence=_gas_trend(evaluation.periodic.history), initial_guess_state=guess,
                constraints=self._unavailable_constraints(),
                warm_start_source=source['candidate_id'] if source else None,
                warm_start_normalized_distance=distance, warm_start_source_status=source['status'] if source else None)
            return result
        cycle = getattr(evaluation.periodic, 'final_cycle', None)
        guess = None
        if cycle is not None and cycle.completed:
            guess = _state_record(cycle.states[:, -1], 'four_gas_volumes_m_U', 'reservoir_8_state', direction, periodic=evaluation.usable)
            guess['topology'] = json_values(asdict(cycle.final_topology))
        return self._assessment(evaluation, derived, source, distance, _gas_trend(evaluation.periodic.history), guess)

    def _wall_initial(self, config, wrapper):
        base_initial = build_initial_state(config, wrapper.model); volumes = wrapper.model.volumes(0)
        pressure = base_initial.total_mass*config.gas.gas_constant*config.charge.temperature/volumes.total
        gas = UniformCharge(pressure, config.charge.temperature).create_state(config.gas, volumes).as_array()
        return np.r_[gas, wrapper.heat_in.wall_capacity_j_k*wrapper.heat_in.air_inlet_temperature_k,
            wrapper.heat_out.wall_capacity_j_k*wrapper.heat_out.air_inlet_temperature_k]

    def _wall(self, candidate, design, wrapper, derived, direction, control):
        layout, family = 'four_gas_volumes_m_U_plus_H_i_H_o_wall_energy', 'microtube_wall_10_state'
        source, distance = select_warm_start(control.previous_records, candidate.payload['normalized'], layout, family, direction)
        target = self._wall_initial(design.configuration, wrapper); state = target.copy()
        if source:
            saved = source.get('initial_guess_state') or source['final_periodic_state']; old = np.asarray(saved['values'], dtype=float)
            new_caps = np.array([wrapper.heat_in.wall_capacity_j_k, wrapper.heat_out.wall_capacity_j_k])
            state = rescale_wall_state(old, target[:8:2].sum(), saved['wall_capacities_j_k'], new_caps)
        try:
            periodic = solve_periodic_wall_motor(wrapper, state,
                maximum_cycles=design.configuration.numerical.maximum_cycles,
                progress_callback=control.check,
                settings=self.definition.wall_numerical_settings)
        except (ValueError, RuntimeError, ArithmeticError) as error:
            from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
            status = 'invalid_exchanger' if isinstance(error, MicrotubeDomainError) else 'integration_failure'
            result = rejected(status, f'{type(error).__name__}: {error}'); result.update(integrated=True, derived=derived)
            return result
        caps = [wrapper.heat_in.wall_capacity_j_k, wrapper.heat_out.wall_capacity_j_k]
        guess = (_state_record(periodic.last_complete_state, layout, family, direction,
            periodic=periodic.converged, wall_capacities=caps) if periodic.last_complete_state is not None else None)
        convergence = convergence_summary(periodic.history)
        warm = dict(warm_start_source=source['candidate_id'] if source else None,
            warm_start_normalized_distance=distance, warm_start_source_status=source['status'] if source else None)
        if periodic.status == 'interrupted':
            result = rejected('budget_exhausted', periodic.message)
            result.update(integrated=True, derived=derived, periodic_cycle_count=len(periodic.history),
                periodic_convergence=convergence, initial_guess_state=guess, **warm)
            result['constraints'] = self._unavailable_constraints()
            return result
        cycle = WallDiagnosticCycle(periodic.angles, periodic.trajectory)
        performance = wall_cycle_performance(wrapper, periodic.trajectory) if periodic.converged else None
        def wall_heat_rates(index, _angle, gas_state):
            incoming, outgoing = wrapper.thermal_rates(_angle, periodic.trajectory[:,index])
            return incoming['gas_heat_w'], outgoing['gas_heat_w']
        diagnostics = (extract_cycle_diagnostics(cycle, wrapper.model,
            heat_rate_provider=wall_heat_rates) if periodic.converged else None)
        validity = assess_cycle_validity(cycle, wrapper.model, design.configuration.validity) if periodic.converged else None
        max_re, max_mach = self._tube_validity(wrapper, periodic.angles, periodic.trajectory)
        from dada_solver.exchangers.gas_diagnostics import cycle_microtube_diagnostics
        gas_domains = cycle_microtube_diagnostics(wrapper, periodic.angles, periodic.trajectory)
        if validity is not None:
            validity = finalize_microtube_validity(validity, max_re, max_mach,
                design.configuration.validity.maximum_mach_number,
                requires_laminar=gas_domains is None,
                domain_failures=gas_domains['failed_criteria'] if gas_domains else ())
        evaluation = SimpleNamespace(configuration=design.configuration, usable=periodic.converged,
            status=EvaluationStatus.CONVERGED if periodic.converged else EvaluationStatus.NOT_CONVERGED,
            periodic=SimpleNamespace(message=periodic.message), performance=performance,
            diagnostics=diagnostics, validity=validity, model=wrapper.model, cycle=cycle)
        domain = dict(name='microtube_model_domain', margin=float(min(2300-max_re,
            design.configuration.validity.maximum_mach_number-max_mach)),
            satisfied=bool(max_re < 2300 and max_mach <= design.configuration.validity.maximum_mach_number), available=True)
        if gas_domains is not None:
            satisfied = gas_domains['model_validity']=='valid' and max_mach<=design.configuration.validity.maximum_mach_number
            domain.update(satisfied=satisfied, margin=1. if satisfied else -1.)
        return self._assessment(evaluation, dict(derived, maximum_tube_reynolds=max_re,
            maximum_tube_mach_number=max_mach, microtube_gas_domains=gas_domains), source, distance, convergence, guess, [domain])

    def _tube_validity(self, wrapper, angles, trajectory):
        if any(getattr(w,'requires_flow_context',False) for w in (wrapper.heat_in,wrapper.heat_out)):
            from dada_solver.exchangers.gas_diagnostics import cycle_microtube_diagnostics
            domains = cycle_microtube_diagnostics(wrapper,angles,trajectory)['passages'].values()
            return (max(x['hydraulic_upstream_ranges']['reynolds']['maximum'] for x in domains),
                    max(x['hydraulic_upstream_ranges']['mach']['maximum'] for x in domains))
        max_re = max_mach = 0.; gas = wrapper.model.gas
        for angle, values in zip(angles, trajectory.T):
            state = ThermodynamicState.from_array(values[:8]); temps = state.temperatures(gas)
            pressures = state.pressures(gas, wrapper.model.volumes(float(angle)))
            flows = wrapper.model.evaluate(float(angle), state, initial_valve_topology()).flows
            samples = ((flows.small_to_cold, 2, wrapper.model.small_cold_link),
                (flows.cold_to_large, 2, wrapper.model.cold_large_valve.flow_model),
                (flows.large_to_hot, 3, wrapper.model.large_hot_link),
                (flows.hot_to_small, 3, wrapper.model.hot_small_valve.flow_model))
            for flow, index, link in samples:
                bank = link.bank; area = bank.dimensions()['tube_flow_area_m2']
                max_re = max(max_re, abs(flow)*bank.inner_diameter_m/(area*link.viscosity_pa_s))
                density = pressures[index]/(gas.gas_constant*temps[index])
                speed = math.sqrt(gas.heat_capacity_cp/gas.heat_capacity_cv*gas.gas_constant*temps[index])
                max_mach = max(max_mach, abs(flow)/(density*area*speed))
        return float(max_re), float(max_mach)

    def _assessment(self, evaluation, derived, source, distance, convergence, guess, extra=()):
        objective = self.definition.objective.evaluate(evaluation)
        constraints = [asdict(c.evaluate(evaluation)) for c in self.definition.constraints] + list(extra)
        objective_valid = objective.available and objective.value is not None and np.isfinite(objective.value)
        feasible = evaluation.usable and objective_valid and all(c['available'] and c['satisfied'] for c in constraints)
        status = ('feasible' if feasible else 'converged_infeasible' if evaluation.usable else
            'periodic_non_convergence' if evaluation.status is EvaluationStatus.NOT_CONVERGED else
            'integration_failure')
        metrics = {}
        if evaluation.performance is not None:
            p = evaluation.performance
            metrics = dict(indicated_power_w=p.gas_power, heat_input_w=p.heat_in_power,
                heat_out_w=p.heat_out_power, indicated_thermal_efficiency=p.thermal_efficiency,
                useful_mechanical_power_w=None, mechanical_losses='unknown', conservation=asdict(p.conservation),
                validity=json_values(asdict(evaluation.validity)))
        reasons = [c['name']+(': unavailable' if not c['available'] else ': violated') for c in constraints if not c['available'] or not c['satisfied']]
        if not objective_valid: reasons.append('objective unavailable or nonfinite')
        return dict(status=status, integrated=True, converged=evaluation.usable,
            reason='; '.join(reasons) if evaluation.usable else evaluation.periodic.message,
            objective=asdict(objective), constraints=constraints, metrics=metrics, derived=derived,
            periodic_cycle_count=convergence['cycles_completed'], periodic_convergence=convergence,
            warm_start_source=source['candidate_id'] if source else None,
            warm_start_normalized_distance=distance, warm_start_source_status=source['status'] if source else None,
            final_periodic_state=guess if evaluation.usable else None, initial_guess_state=guess,
            preflight_diagnostics=None)
