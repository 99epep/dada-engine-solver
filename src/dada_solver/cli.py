"""Command-line entry point for periodic DADA thermodynamic simulations."""

from __future__ import annotations

import argparse
from pathlib import Path

from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import (
    build_initial_state,
    build_model,
    build_periodic_solver,
    initial_valve_topology,
)
from dada_solver.performance import calculate_cycle_performance
from dada_solver.periodic import PeriodicStatus
from dada_solver.plotting import plot_cycle_diagnostics
from dada_solver.reporting import format_simulation_report
from dada_solver.results import extract_cycle_diagnostics
from dada_solver.state import ThermodynamicState
from dada_solver.validity import assess_cycle_validity
from dada_solver.humidity import assess_moisture_phase_change_risk


def main(arguments: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Simulate the first-level DADA thermodynamic model."
    )
    parser.add_argument("configuration", type=Path, help="TOML configuration file")
    parser.add_argument("--plot", type=Path, help="Write a diagnostic plot to this path")
    parser.add_argument(
        "--integration-profile",
        action="store_true",
        help="Include final-cycle solver statistics for every integration segment",
    )
    options = parser.parse_args(arguments)

    try:
        configuration = load_simulation_configuration(options.configuration)
        model = build_model(configuration)
        print(f"signed_crank_angular_speed_rad_per_s = {configuration.angular_speed:.12e}")
        print(f"requested_operation = {'motor' if configuration.motor_operation else 'refrigeration'}")
        print(f"heat_in_reservoir_temperature_K = {model.cold_heat_transfer.reservoir_temperature:.12e}")
        print(f"heat_out_reservoir_temperature_K = {model.hot_heat_transfer.reservoir_temperature:.12e}")
        initial_state = build_initial_state(configuration, model)
        periodic = build_periodic_solver(configuration, model).solve(
            initial_state, initial_valve_topology()
        )
        if periodic.status is not PeriodicStatus.CONVERGED or periodic.final_cycle is None:
            print(f"periodic_status = {periodic.status.value}")
            print(f"message = {periodic.message}")
            return 2

        cycle = periodic.final_cycle
        cycle_initial_state = ThermodynamicState.from_array(cycle.states[:, 0])
        performance = calculate_cycle_performance(
            cycle, cycle_initial_state, model.signed_angular_speed
        )
        diagnostics = extract_cycle_diagnostics(cycle, model)
        validity = assess_cycle_validity(cycle, model, configuration.validity)
        moisture = assess_moisture_phase_change_risk(
            cycle,
            model,
            configuration.humidity_screening,
            configuration.charge.temperature,
            configuration.charge.resolved_pressure(
                configuration.gas, model.volumes(0.0).total
            ),
        )
        print(
            format_simulation_report(
                periodic,
                performance,
                diagnostics,
                validity,
                initial_state.total_mass,
                include_integration_segments=options.integration_profile,
                moisture=moisture,
            )
        )
        if options.plot is not None:
            plot_cycle_diagnostics(
                cycle,
                model,
                options.plot,
                humidity=configuration.humidity_screening,
                charge_temperature=configuration.charge.temperature,
                charge_pressure=configuration.charge.resolved_pressure(
                    configuration.gas, model.volumes(0.0).total
                ),
            )
            print(f"diagnostic_plot = {options.plot}")
        return 0
    except (OSError, ValueError) as error:
        print(f"simulation_error = {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
