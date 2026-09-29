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
    statistics_callback: object = None
    progress_callback: object = None
    def check(self, _progress=None):
        if self.progress_callback is not None: self.progress_callback(_progress)
        if self.deadline is not None and self.clock() >= self.deadline:
            raise IntegrationInterrupted('Campaign wall-clock deadline reached.')

    def measure(self, phase, function, *args, **kwargs):
        """Optional side-channel timings; persisted campaign records stay unchanged."""
        if self.statistics_callback is None:
            return function(*args, **kwargs)
        before = time.perf_counter()
        status = 'completed'
        try:
            return function(*args, **kwargs)
        except Exception:
            status = 'failed'
            raise
        finally:
            self.statistics_callback(dict(phase=phase, status=status,
                elapsed_seconds=time.perf_counter()-before))


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
        timings = {}
        if hasattr(self.definition, 'study'):
            original_callback = control.statistics_callback
            def observe(record):
                if 'elapsed_seconds' in record:
                    phase = record['phase']
                    timings[phase] = timings.get(phase, 0.0) + record['elapsed_seconds']
                if original_callback is not None: original_callback(record)
            control = replace(control, statistics_callback=observe)
        candidate_budget = getattr(self.definition, 'candidate_budget_seconds', None)
        if candidate_budget is not None:
            deadline = control.clock() + candidate_budget
            control = replace(control, deadline=min(deadline, control.deadline) if control.deadline is not None else deadline)
        result = control.measure('candidate_evaluation', self._evaluate_with_control, candidate, control)
        if hasattr(self.definition, 'study'):
            from dada_solver.research.margins import enrich_constraints
            result['timing_seconds'] = timings
            if not result['integrated'] and result.get('preflight_diagnostics'):
                diagnostics=result['preflight_diagnostics']
                if all(isinstance(d,dict) and 'name' in d and 'margin' in d for d in diagnostics):
                    result['constraints']=list(diagnostics)+self._unavailable_constraints()
            result['constraints']=enrich_constraints(result['constraints'],self.definition.study.data['constraints'])
        if control.statistics_callback is not None:
            control.statistics_callback(dict(phase='candidate_result', status=result['status'],
                warm_start_source=result.get('warm_start_source'),
                warm_start_normalized_distance=result.get('warm_start_normalized_distance'),
                periodic_cycle_count=result['periodic_cycle_count']))
        return json_values(result) if hasattr(self.definition, 'study') else result

    def _evaluate_with_control(self, candidate, control):
        try:
            physical = dict(self.definition.fixed_parameters)
            physical.update(candidate.payload['physical'])
            design = control.measure('parameter_preflight', self.definition.adapter.build, physical)
        except PreflightRejection as error: return rejected(error.status, str(error), error.diagnostics)
        try: built = control.measure('build', design.build)
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
        mechanical=getattr(design.kinematics,'mechanical_diagnostics',())
        if mechanical:
            derived['mechanical_constraints']=list(mechanical)
            derived['kinematic_metrics']=list(getattr(design.kinematics,'mechanical_metrics',()))
        if isinstance(built, AirWallMotor):
            result=self._wall(candidate, design, built, derived, direction, control)
        else:
            result=self._reservoir(candidate, design, built, derived, direction, control)
        if mechanical: result['constraints']=list(mechanical)+result['constraints']
        return result

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
        from dada_solver.fluids import CaloricallyPerfectGas
        pressure = (base_initial.total_mass*config.gas.gas_constant*config.charge.temperature/volumes.total
            if type(config.gas) is CaloricallyPerfectGas else
            config.gas.state_from_rho_t(base_initial.total_mass/volumes.total,config.charge.temperature).pressure)
        gas = UniformCharge(pressure, config.charge.temperature).create_state(config.gas, volumes).as_array()
        return np.r_[gas, wrapper.heat_in.wall_capacity_j_k*wrapper.heat_in.external_inlet_temperature_k,
            wrapper.heat_out.wall_capacity_j_k*wrapper.heat_out.external_inlet_temperature_k]

    def _wall(self, candidate, design, wrapper, derived, direction, control):
        layout, family = 'four_gas_volumes_m_U_plus_H_i_H_o_wall_energy', 'microtube_wall_10_state'
        source, distance = select_warm_start(control.previous_records, candidate.payload['normalized'], layout, family, direction)
        from dada_solver.tabulated_fluid import FluidDomainError
        try: target = self._wall_initial(design.configuration, wrapper)
        except FluidDomainError as error:
            result=rejected('invalid_fluid_domain',str(error))
            result.update(derived=derived,constraints=self._unavailable_constraints())
            return result
        state = target.copy()
        external = getattr(self.definition, 'initial_wall_state', None) if source is None else None
        caps = [wrapper.heat_in.wall_capacity_j_k, wrapper.heat_out.wall_capacity_j_k]
        if source:
            saved = source.get('initial_guess_state') or source['final_periodic_state']; old = np.asarray(saved['values'], dtype=float)
            new_caps = np.array([wrapper.heat_in.wall_capacity_j_k, wrapper.heat_out.wall_capacity_j_k])
            state = rescale_wall_state(old, target[:8:2].sum(), saved['wall_capacities_j_k'], new_caps)
        elif external:
            # The fixed-inventory reference preserves the original gas state
            # bit for bit. Only wall energies scale with changed wall capacity.
            state = np.asarray(external['values'], dtype=float).copy()
            state[8:10] *= np.asarray(caps)/np.asarray(external['wall_capacities_j_k'])
        backend_last = {}
        first_microtube_failure = None
        def retain_failure(result):
            if first_microtube_failure is not None:
                result.setdefault('diagnostics', {})['first_microtube_failure'] = first_microtube_failure
            return result
        def record_backend(record):
            nonlocal first_microtube_failure
            if record['phase'] == 'microtube_failure' and first_microtube_failure is None:
                first_microtube_failure = record['snapshot']
            if record['phase']=='rhs_backend': backend_last.update(record)
            if control.statistics_callback is not None: control.statistics_callback(record)
        def solve(initial):
            return control.measure('periodic_integration', solve_periodic_wall_motor, wrapper, initial,
                maximum_cycles=design.configuration.numerical.maximum_cycles,
                progress_callback=control.check,
                settings=self.definition.wall_numerical_settings,
                backend=self.definition.wall_backend,
                adaptive_acceleration=getattr(self.definition, "adaptive_wall_acceleration", None),
                statistics_callback=record_backend)
        safe_retry = False
        try:
            from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
            try:
                periodic = solve(state)
            except MicrotubeDomainError as error:
                from dada_solver.exchangers.failure_diagnostics import microtube_failure_snapshot
                if first_microtube_failure is None:
                    first_microtube_failure = microtube_failure_snapshot(error)
                if not getattr(self.definition, 'safe_domain_retry', False): raise
                control.check()
                safe_retry = True
                safe = target.copy()
                safe[8:10] = state[8:10]
                periodic = solve(safe)
        except IntegrationInterrupted as error:
            result = rejected('budget_exhausted', str(error))
            result.update(integrated=True, derived=derived, constraints=self._unavailable_constraints())
            return retain_failure(result)
        except (ValueError, RuntimeError, ArithmeticError) as error:
            from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
            status = ('invalid_fluid_domain' if isinstance(error,FluidDomainError) else
                'invalid_exchanger' if isinstance(error, MicrotubeDomainError) else 'integration_failure')
            result = rejected(status, f'{type(error).__name__}: {error}'); result.update(integrated=True, derived=derived)
            if getattr(error, 'native_solver_stderr', None):
                result['technical_diagnostics'] = dict(native_solver_stderr=error.native_solver_stderr)
            if backend_last: result['rhs_backend'] = backend_last
            if hasattr(self.definition, 'safe_domain_retry'):
                result.update(safe_retry_used=safe_retry,
                    warm_start_source=source['candidate_id'] if source else external['source_candidate_id'] if external else None)
            return retain_failure(result)
        caps = [wrapper.heat_in.wall_capacity_j_k, wrapper.heat_out.wall_capacity_j_k]
        guess = (_state_record(periodic.last_complete_state, layout, family, direction,
            periodic=periodic.converged, wall_capacities=caps) if periodic.last_complete_state is not None else None)
        convergence = convergence_summary(periodic.history)
        if self.definition.wall_backend.name != "python":
            convergence["rhs_backend"] = periodic.backend_statistics
        warm = dict(warm_start_source=source['candidate_id'] if source else None,
            warm_start_normalized_distance=distance, warm_start_source_status=source['status'] if source else None)
        if external:
            warm.update(warm_start_source=external['source_candidate_id'],
                        warm_start_source_status='external_reference_initial_guess')
        if hasattr(self.definition, 'safe_domain_retry'):
            warm['safe_retry_used'] = safe_retry
        if periodic.status == 'interrupted':
            result = rejected('budget_exhausted', periodic.message)
            result.update(integrated=True, derived=derived, periodic_cycle_count=len(periodic.history),
                periodic_convergence=convergence, initial_guess_state=guess, **warm)
            result['constraints'] = self._unavailable_constraints()
            return retain_failure(result)
        if not periodic.converged:
            result = rejected('periodic_non_convergence', periodic.message)
            result.update(integrated=True, derived=derived, periodic_cycle_count=len(periodic.history),
                          periodic_convergence=convergence, initial_guess_state=guess,
                          constraints=self._unavailable_constraints(), **warm)
            return retain_failure(result)
        cycle = WallDiagnosticCycle(periodic.angles, periodic.trajectory)
        performance = wall_cycle_performance(wrapper, periodic.trajectory) if periodic.converged else None
        from dada_solver.diagnostic_replay import replay_wall_trajectory
        replay = control.measure('shared_replay', replay_wall_trajectory, wrapper, cycle, periodic.trajectory) if periodic.converged else None
        diagnostics = (control.measure('generic_diagnostics', extract_cycle_diagnostics, cycle, wrapper.model,
            replay=replay) if periodic.converged else None)
        validity = control.measure('generic_validity', assess_cycle_validity, cycle, wrapper.model, design.configuration.validity, replay=replay) if periodic.converged else None
        from dada_solver.exchangers.gas_diagnostics import cycle_microtube_diagnostics
        gas_domains = control.measure('microtube_diagnostics', cycle_microtube_diagnostics, wrapper, periodic.angles, periodic.trajectory, replay=replay)
        max_re, max_mach = self._tube_validity(wrapper, periodic.angles, periodic.trajectory, gas_domains=gas_domains,replay=replay)
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
        result = self._assessment(evaluation, dict(derived, maximum_tube_reynolds=max_re,
            maximum_tube_mach_number=max_mach, microtube_gas_domains=gas_domains), source, distance, convergence, guess, [domain])
        if hasattr(self.definition, 'study'):
            result.update(warm)
            result['metrics'].update(
                maximum_pressure_pa=max(v.maximum for v in diagnostics.pressure_extrema.values()),
                maximum_temperature_k=max(v.maximum for v in diagnostics.temperature_extrema.values()),
                maximum_absolute_mass_flow_kg_s=max(max(abs(v.minimum),abs(v.maximum)) for v in diagnostics.mass_flow_extrema.values()),
                total_mass_kg=float(periodic.trajectory[:8:2,-1].sum()))
            result['diagnostics'] = json_values(asdict(diagnostics))
            result['derived']['hardware'] = dict(heat_in=dict(design.heat_in.build().metadata), heat_out=dict(design.heat_out.build().metadata))
            if self.definition.study.data['schema_version']==3:
                streams={}
                for i,(side,exchanger) in enumerate((('heat_in',wrapper.heat_in),('heat_out',wrapper.heat_out))):
                    stream=exchanger.external_stream
                    outlets=[s.walls[i].external_outlet_temperature_k for s in replay.samples]
                    streams[side]=dict(asdict(stream),capacity_rate_w_k=stream.capacity_rate_w_k,
                        outlet_minimum_k=min(outlets) if all(v is not None for v in outlets) else None,
                        outlet_maximum_k=max(outlets) if all(v is not None for v in outlets) else None,
                        heat_into_machine_per_cycle_j=performance.heat_in_per_cycle if i==0 else performance.heat_out_per_cycle,
                        mean_heat_into_machine_w=performance.heat_in_power if i==0 else performance.heat_out_power,
                        external_loop_losses='excluded; hydraulics and pump/fan consumption unmodelled')
                result['derived']['external_streams']=streams
                result['metrics'].update(operating_mode=performance.operating_mode.value,
                    heating_power_w=performance.heating_power,heating_cop=performance.heating_cop)
                result['rhs_backend']=periodic.backend_statistics
            else:
                result['derived']['air_inlet_temperatures_k'] = dict(heat_in=wrapper.heat_in.external_inlet_temperature_k, heat_out=wrapper.heat_out.external_inlet_temperature_k)
                air = {}
                for i, (side, exchanger, ports) in enumerate((
                        ('heat_in',wrapper.heat_in,('small_to_cold','cold_to_large')),
                        ('heat_out',wrapper.heat_out,('large_to_hot','hot_to_small')))):
                    peak = max(max(abs(diagnostics.mass_flow_extrema[p].minimum),abs(diagnostics.mass_flow_extrema[p].maximum)) for p in ports)
                    capacity = exchanger.air_mass_flow_kg_s * exchanger.air_cp_j_kg_k
                    air[side] = dict(peak_internal_mass_flow_kg_s=peak,
                        external_to_peak_internal_capacity_rate_ratio=capacity/(peak*design.configuration.gas.heat_capacity_cp) if peak else None,
                        maximum_external_air_temperature_change_k=max(abs(s.walls[i].air_heat_w)/capacity for s in replay.samples),
                        scope='sampled_cycle; fixed_external_flow; finite_film_resistance_retained')
                result['derived']['external_air_capacity_diagnostics'] = air
            result['derived']['local_reflux'] = dict(
                threshold_kg_s=1e-8,
                detected=any(v.minimum < -1e-8 for v in diagnostics.mass_flow_extrema.values()),
                minimum_signed_flows_kg_s={k:v.minimum for k,v in diagnostics.mass_flow_extrema.items()})
        return retain_failure(result)

    def _tube_validity(self, wrapper, angles, trajectory, *, gas_domains=None, replay=None):
        if replay is not None: replay.require(wrapper=wrapper, angles=angles, trajectory=trajectory)
        if any(getattr(w,'requires_flow_context',False) for w in (wrapper.heat_in,wrapper.heat_out)):
            from dada_solver.exchangers.gas_diagnostics import cycle_microtube_diagnostics
            if gas_domains is None:
                gas_domains = cycle_microtube_diagnostics(wrapper,angles,trajectory)
            domains = gas_domains['passages'].values()
            return (max(x['hydraulic_upstream_ranges']['reynolds']['maximum'] for x in domains),
                    max(x['hydraulic_upstream_ranges']['mach']['maximum'] for x in domains))
        max_re = max_mach = 0.; gas = wrapper.model.gas
        for index,(angle, values) in enumerate(zip(angles, trajectory.T)):
            if replay is None:
                state = ThermodynamicState.from_array(values[:8]); temps = state.temperatures(gas, wrapper.model.volumes(float(angle)))
                pressures = state.pressures(gas, wrapper.model.volumes(float(angle)))
                flows = wrapper.model.evaluate(float(angle), state, initial_valve_topology()).flows
            else:
                point = replay.samples[index].point
                temps, pressures, flows = point.temperatures, point.pressures, point.flows
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
            if getattr(self.definition,'identity',{}).get('definition_kind') in ('research_v2','research_v3'):
                cooling=self.definition.objective.name in ('maximize_cooling_cop','maximize_cooling_power')
                metrics.update(cooling_power_w=p.cooling_power if cooling else None, cooling_cop=p.cooling_cop if cooling else None,
                               indicated_mechanical_input_power_w=p.mechanical_input_power if cooling else None)
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
