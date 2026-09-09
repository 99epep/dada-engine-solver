"""TOML loading for explicit constrained-sizing problems."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tomllib

from dada_solver.configuration import load_simulation_configuration
from dada_solver.mechanics import PistonFaceAreas
from dada_solver.thermal_load_configuration import (
    CoolingCellScenario,
    load_cooling_cell_scenario,
)
from dada_solver.sizing.constraints import (
    MaximumAbsoluteMassFlow,
    MaximumMachNumber,
    MaximumIsothermalityError,
    MaximumAbsoluteGeneralizedGasTorque,
    MaximumPistonGasForce,
    MaximumPressure,
    MaximumPressureEqualizationError,
    MaximumTemperature,
    MinimumCoolingPower,
    RequirePeriodicConvergence,
    RequireNominalCycleTopology,
    RequireValidThermodynamicModel,
    SizingConstraint,
    CompleteCoolingTaskWithinTime,
    MaximumMechanicalInputPower,
)
from dada_solver.sizing.design import DesignParameter, DesignVariable
from dada_solver.sizing.evaluator import ThermodynamicSizingEvaluator
from dada_solver.sizing.objectives import (
    MaximizeCoolingCop,
    MaximizeThermalEfficiency,
    MaximizeMotorPower,
    MinimizeChargingPressure,
    MinimizeTotalSweptVolume,
    MinimizeTotalUa,
    SizingObjective,
)
from dada_solver.sizing.optimizer import OptimizationSettings
from dada_solver.sizing.problem import SizingProblem


@dataclass(frozen=True, slots=True)
class LoadedSizingProblem:
    problem: SizingProblem
    optimization_settings: OptimizationSettings
    base_configuration_path: Path
    cooling_cell_scenario: CoolingCellScenario | None


def load_sizing_problem(path: str | Path) -> LoadedSizingProblem:
    """Load variables, objective, constraints and numerical scales from TOML."""

    sizing_path = Path(path)
    with sizing_path.open("rb") as stream:
        data = tomllib.load(stream)
    try:
        base_path = Path(data["problem"]["base_configuration"])
        if not base_path.is_absolute():
            base_path = sizing_path.parent / base_path
        base = load_simulation_configuration(base_path)
        scenario = _load_optional_cooling_scenario(data["problem"], sizing_path)
        variables = tuple(_load_variable(item) for item in data["variables"])
        objective = _load_objective(data["objective"])
        constraints = tuple(
            _load_constraint(item, scenario) for item in data.get("constraints", [])
        )
        optimizer = data["optimizer"]
        scales = {
            str(name): float(value)
            for name, value in optimizer["constraint_scales"].items()
        }
        settings = OptimizationSettings(
            objective_scale=float(optimizer["objective_scale"]),
            constraint_scales=scales,
            unavailable_objective_penalty=float(
                optimizer["unavailable_objective_penalty"]
            ),
            unavailable_constraint_margin=float(
                optimizer["unavailable_constraint_margin"]
            ),
            maximum_iterations=int(optimizer["maximum_iterations"]),
            function_tolerance=float(optimizer["function_tolerance"]),
        )
    except KeyError as error:
        raise ValueError(f"Missing sizing configuration field: {error.args[0]}") from error
    problem = SizingProblem(
        variables=variables,
        objective=objective,
        constraints=constraints,
        evaluator=ThermodynamicSizingEvaluator(base),
    )
    return LoadedSizingProblem(problem, settings, base_path, scenario)


def _load_variable(data: dict[str, object]) -> DesignVariable:
    try:
        parameter = DesignParameter(str(data["parameter"]))
    except ValueError as error:
        raise ValueError(f"Unsupported design parameter: {data['parameter']}") from error
    return DesignVariable(
        parameter=parameter,
        lower_bound=float(data["lower_bound"]),
        upper_bound=float(data["upper_bound"]),
        initial_value=float(data["initial_value"]),
    )


def _load_objective(data: dict[str, object]) -> SizingObjective:
    objective_type = str(data["type"])
    objectives: dict[str, SizingObjective] = {
        "maximize_cooling_cop": MaximizeCoolingCop(),
        "maximize_thermal_efficiency": MaximizeThermalEfficiency(),
        "maximize_motor_power": MaximizeMotorPower(),
        "minimize_total_swept_volume": MinimizeTotalSweptVolume(),
        "minimize_charging_pressure": MinimizeChargingPressure(),
        "minimize_total_ua": MinimizeTotalUa(),
    }
    try:
        return objectives[objective_type]
    except KeyError as error:
        raise ValueError(f"Unsupported sizing objective: {objective_type}") from error


def _load_constraint(
    data: dict[str, object],
    scenario: CoolingCellScenario | None,
) -> SizingConstraint:
    constraint_type = str(data["type"])
    if constraint_type == "minimum_motor_power":
        from dada_solver.sizing.constraints import MinimumMotorPower
        return MinimumMotorPower(float(data["required_power"]))
    if constraint_type == "minimum_cooling_power":
        return MinimumCoolingPower(float(data["required_power"]))
    if constraint_type == "complete_reference_cooling_task":
        return CompleteCoolingTaskWithinTime(
            _require_cooling_scenario(scenario).reference_task
        )
    if constraint_type == "complete_ambitious_cooling_task":
        return CompleteCoolingTaskWithinTime(
            _require_cooling_scenario(scenario).ambitious_task
        )
    if constraint_type == "maximum_nominal_human_power":
        return MaximumMechanicalInputPower(
            _require_cooling_scenario(scenario).human_power.nominal_power
        )
    if constraint_type == "maximum_pressure":
        return MaximumPressure(float(data["limit"]))
    if constraint_type == "maximum_temperature":
        return MaximumTemperature(float(data["limit"]))
    if constraint_type == "maximum_absolute_mass_flow":
        return MaximumAbsoluteMassFlow(float(data["limit"]))
    if constraint_type == "maximum_mach_number":
        return MaximumMachNumber(float(data["limit"]))
    if constraint_type == "maximum_pressure_equalization_error":
        return MaximumPressureEqualizationError(float(data["limit"]))
    if constraint_type == "maximum_isothermality_error":
        return MaximumIsothermalityError(float(data["limit"]))
    if constraint_type == "nominal_cycle_topology":
        return RequireNominalCycleTopology()
    if constraint_type == "valid_thermodynamic_model":
        return RequireValidThermodynamicModel()
    if constraint_type == "periodic_convergence":
        return RequirePeriodicConvergence()
    if constraint_type == "maximum_absolute_generalized_gas_torque":
        return MaximumAbsoluteGeneralizedGasTorque(float(data["limit"]))
    if constraint_type == "maximum_piston_gas_force":
        return MaximumPistonGasForce(
            limit=float(data["limit"]),
            piston_areas=PistonFaceAreas(
                small=float(data["small_piston_area"]),
                large=float(data["large_piston_area"]),
            ),
        )
    raise ValueError(f"Unsupported sizing constraint: {constraint_type}")


def _load_optional_cooling_scenario(
    problem_data: dict[str, object],
    sizing_path: Path,
) -> CoolingCellScenario | None:
    value = problem_data.get("cooling_load_configuration")
    if value is None:
        return None
    path = Path(str(value))
    if not path.is_absolute():
        path = sizing_path.parent / path
    return load_cooling_cell_scenario(path)


def _require_cooling_scenario(
    scenario: CoolingCellScenario | None,
) -> CoolingCellScenario:
    if scenario is None:
        raise ValueError(
            "Cooling-task constraint requires problem.cooling_load_configuration."
        )
    return scenario
