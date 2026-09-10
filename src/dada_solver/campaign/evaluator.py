"""Campaign assessment through the existing physical evaluator and sizing objects."""
from dataclasses import asdict
from enum import Enum
import numpy as np
from dada_solver.campaign.adapters import PreflightRejection
from dada_solver.kinematics import KinematicConstraintViolation
from dada_solver.sizing.evaluator import evaluate_configuration, EvaluationStatus
from dada_solver.sizing.problem import SizingAssessment


def json_values(value):
    if isinstance(value, Enum): return value.value
    if isinstance(value, dict): return {str(k): json_values(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)): return [json_values(v) for v in value]
    if isinstance(value, np.generic): return value.item()
    return value


def rejected(status, reason, diagnostics=None):
    return dict(status=status, reason=reason, integrated=False, converged=False,
        objective=None, constraints=[], metrics={}, derived={}, periodic_cycle_count=0,
        warm_start_source=None, final_periodic_state=None, preflight_diagnostics=diagnostics)


class MachineEvaluator:
    """Cold-start evaluation; final states are retained for future compatible reuse."""
    def __init__(self, definition):
        self.definition = definition

    def evaluate(self, candidate):
        try:
            design = self.definition.adapter.build(candidate.payload['physical'])
        except PreflightRejection as error:
            return rejected(error.status, str(error), error.diagnostics)
        try:
            model = design.build()
        except KinematicConstraintViolation as error:
            return rejected('invalid_kinematics', str(error), [asdict(d) for d in error.diagnostics])
        except (ValueError, ArithmeticError) as error:
            return rejected('invalid_exchanger', str(error))
        derived = dict(small_volume_limits=asdict(model.machine_volumes.small_cylinder),
            large_volume_limits=asdict(model.machine_volumes.large_cylinder),
            heat_in_gas_volume_m3=model.machine_volumes.cold_heat_exchanger,
            heat_out_gas_volume_m3=model.machine_volumes.hot_heat_exchanger,
            small_physical_stroke_m=model.kinematics.small_physical_stroke,
            large_physical_stroke_m=model.kinematics.large_physical_stroke)
        try:
            evaluation = evaluate_configuration(design.configuration, model=model)
        except (ValueError, RuntimeError, ArithmeticError) as error:
            result = rejected('integration_failure', f'{type(error).__name__}: {error}')
            result.update(integrated=True, derived=derived)
            return result
        objective = self.definition.objective.evaluate(evaluation)
        constraints = tuple(c.evaluate(evaluation) for c in self.definition.constraints)
        assessment = SizingAssessment(evaluation, objective, constraints)
        converged = evaluation.usable
        objective_valid = objective.available and objective.value is not None and np.isfinite(objective.value)
        if converged:
            status = 'feasible' if assessment.feasible and objective_valid else 'converged_infeasible'
        elif evaluation.status is EvaluationStatus.NOT_CONVERGED:
            status = 'periodic_non_convergence'
        else:
            status = 'integration_failure'
        metrics, final_state = {}, None
        if evaluation.performance is not None:
            performance = evaluation.performance
            metrics = dict(indicated_power_w=performance.gas_power,
                heat_input_w=performance.heat_in_power, heat_out_w=performance.heat_out_power,
                indicated_thermal_efficiency=performance.thermal_efficiency,
                useful_mechanical_power_w=None, mechanical_losses='unknown',
                conservation=asdict(performance.conservation),
                validity=json_values(asdict(evaluation.validity)))
            cycle = evaluation.cycle
            final_state = dict(state_layout='four_gas_volumes_m_U',
                values=cycle.states[:, -1].tolist(),
                topology=json_values(asdict(cycle.final_topology)),
                total_mass_kg=float(cycle.states[::2, -1].sum()),
                reusable_only_as_initial_guess=True)
        reasons = [c.name+(': unavailable' if not c.available else ': violated')
                   for c in constraints if not c.available or not c.satisfied]
        if not objective_valid: reasons.append('objective unavailable or nonfinite')
        return dict(status=status, integrated=True, converged=converged,
            reason='; '.join(reasons) if converged else evaluation.periodic.message,
            objective=asdict(objective), constraints=[asdict(c) for c in constraints],
            metrics=metrics, derived=derived, periodic_cycle_count=len(evaluation.periodic.history),
            warm_start_source=None, final_periodic_state=final_state, preflight_diagnostics=None)
