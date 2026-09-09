"""Command-line entry point for one configured constrained-sizing run."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from dada_solver.sizing.configuration import load_sizing_problem
from dada_solver.sizing.export import write_feasibility_csv
from dada_solver.sizing.feasibility import (
    FeasibilityStudySettings,
    run_feasibility_study,
)
from dada_solver.sizing.optimizer import SlsqpSizingOptimizer
from dada_solver.sizing.design import DesignParameter
from dada_solver.sizing.sensitivity import SamplingScale, run_one_at_a_time_sensitivity
from dada_solver.sizing.reporting import (
    format_assessment_report,
    format_feasibility_report,
    format_optimization_report,
)


def main(arguments: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run constrained DADA sizing.")
    parser.add_argument("configuration", type=Path, help="Sizing TOML file")
    parser.add_argument(
        "--initial-only",
        action="store_true",
        help="Evaluate and report the initial point without optimization",
    )
    parser.add_argument(
        "--feasibility-samples",
        type=int,
        help="Run this many Latin-hypercube samples instead of optimization",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Non-negative random seed for feasibility sampling",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        help="Write feasibility assessments to a CSV file",
    )
    parser.add_argument(
        "--sensitivity",
        choices=[item.value for item in DesignParameter],
        help="Sweep one configured variable while holding the others at their initial values",
    )
    parser.add_argument(
        "--sensitivity-points", type=int, default=9,
        help="Number of points in a one-at-a-time sensitivity sweep",
    )
    parser.add_argument(
        "--sensitivity-scale", choices=[item.value for item in SamplingScale],
        default=SamplingScale.LINEAR.value,
        help="Sampling scale for a one-at-a-time sensitivity sweep",
    )
    parser.add_argument(
        "--sensitivity-lower", type=float,
        help="Optional local lower bound for the sensitivity sweep",
    )
    parser.add_argument(
        "--sensitivity-upper", type=float,
        help="Optional local upper bound for the sensitivity sweep",
    )
    options = parser.parse_args(arguments)
    try:
        loaded = load_sizing_problem(options.configuration)
        if options.sensitivity is not None:
            result = run_one_at_a_time_sensitivity(
                loaded.problem,
                DesignParameter(options.sensitivity),
                options.sensitivity_points,
                SamplingScale(options.sensitivity_scale),
                options.sensitivity_lower,
                options.sensitivity_upper,
            )
            print("DADA one-at-a-time sensitivity study")
            print(f"parameter = {result.parameter.value}")
            print(f"scale = {result.scale.value}")
            print(
                "value,evaluation_status,feasible,cooling_power,"
                "mechanical_input_power,cooling_cop,minimum_H_i_temperature,"
                "maximum_temperature,theta_zero_H_i_L_pressure_difference,"
                "H_i_to_L_closing_angle_deg,H_o_to_S_opening_angle_deg,"
                "phase_one_both_closed_angle_deg,relative_mass_residual,"
                "relative_energy_residual,maximum_pressure_equalization_error,"
                "cycle_topology,validity"
            )
            for assessment in result.assessments:
                value = assessment.evaluation.point.values[result.parameter]
                evaluation = assessment.evaluation
                performance = evaluation.performance
                diagnostics = evaluation.diagnostics
                validity = evaluation.validity
                cooling_power = "unavailable" if performance is None else f"{performance.cooling_power:.12e}"
                mechanical_power = "unavailable" if performance is None else f"{performance.mechanical_input_power:.12e}"
                cop = None if performance is None else performance.cooling_cop
                cop_text = "unavailable" if cop is None else f"{cop:.12e}"
                topology = "unavailable" if diagnostics is None else diagnostics.topology.classification.value
                validity_text = "unavailable" if validity is None else validity.verdict.value
                cold_minimum = (
                    "unavailable" if diagnostics is None
                    else f"{diagnostics.temperature_extrema['C'].minimum:.12e}"
                )
                maximum_temperature = (
                    "unavailable" if diagnostics is None
                    else f"{max(item.maximum for item in diagnostics.temperature_extrema.values()):.12e}"
                )
                pressure_difference = "unavailable"
                cold_closing = "unavailable"
                hot_opening = "unavailable"
                both_closed = "unavailable"
                mass_residual = "unavailable"
                energy_residual = "unavailable"
                pressure_equalization = "unavailable"
                if performance is not None:
                    mass_residual = (
                        f"{performance.conservation.relative_mass_residual:.12e}"
                    )
                    energy_residual = (
                        f"{performance.conservation.relative_energy_residual:.12e}"
                    )
                if validity is not None:
                    pressure_equalization = (
                        f"{validity.maximum_pressure_equalization_error:.12e}"
                    )
                if evaluation.cycle is not None and evaluation.model is not None:
                    from dada_solver.state import ThermodynamicState

                    initial_state = ThermodynamicState.from_array(
                        evaluation.cycle.states[:, 0]
                    )
                    pressures = initial_state.pressures(
                        evaluation.model.gas, evaluation.model.volumes(0.0)
                    )
                    pressure_difference = f"{pressures[2] - pressures[1]:.12e}"
                if diagnostics is not None:
                    cold_event = next(
                        (
                            item for item in diagnostics.valve_events
                            if item.valve_name == "cold_to_large"
                            and item.transition == "closed"
                        ),
                        None,
                    )
                    hot_event = next(
                        (
                            item for item in diagnostics.valve_events
                            if item.valve_name == "hot_to_small"
                            and item.transition == "open"
                        ),
                        None,
                    )
                    if cold_event is not None:
                        cold_closing = f"{cold_event.angle_degrees:.12e}"
                    if hot_event is not None:
                        hot_opening = f"{hot_event.angle_degrees:.12e}"
                    if cold_event is not None and hot_event is not None:
                        duration = hot_event.angle_degrees - cold_event.angle_degrees
                        if duration >= -math.sqrt(math.ulp(1.0)):
                            both_closed = f"{max(0.0, duration):.12e}"
                print(
                    f"{value:.12e},{evaluation.status.value},"
                    f"{str(assessment.feasible).lower()},{cooling_power},"
                    f"{mechanical_power},{cop_text},{cold_minimum},"
                    f"{maximum_temperature},{pressure_difference},{cold_closing},"
                    f"{hot_opening},{both_closed},{mass_residual},{energy_residual},"
                    f"{pressure_equalization},{topology},{validity_text}",
                    flush=True,
                )
            return 0
        if options.initial_only:
            assessment = loaded.problem.assess(loaded.problem.initial_point)
            print(format_assessment_report(assessment))
            return 0 if assessment.feasible else 2
        if options.feasibility_samples is not None:
            settings = FeasibilityStudySettings(
                sample_count=options.feasibility_samples,
                random_seed=options.seed,
                include_initial_point=True,
            )

            def progress(index, total, assessment):
                print(
                    f"feasibility_progress = {index}/{total}, "
                    f"status={assessment.evaluation.status.value}, "
                    f"feasible={str(assessment.feasible).lower()}",
                    flush=True,
                )

            result = run_feasibility_study(
                loaded.problem,
                settings,
                progress_callback=progress,
            )
            print(format_feasibility_report(result))
            if options.csv is not None:
                write_feasibility_csv(result, options.csv)
                print(f"feasibility_csv = {options.csv}")
            return 0 if result.feasible_count else 2
        result = SlsqpSizingOptimizer(
            loaded.problem, loaded.optimization_settings
        ).optimize()
        print(format_optimization_report(result))
        return 0 if result.success else 2
    except (OSError, ValueError) as error:
        print(f"sizing_error = {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
