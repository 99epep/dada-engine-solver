"""Configuration loading for coupled cycle/exchanger fixed-point sizing."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tomllib

from dada_solver.configuration import SimulationConfiguration, load_simulation_configuration
from dada_solver.exchangers.configuration import (
    ExchangerScreeningConfiguration,
    load_exchanger_screening_configuration,
)
from dada_solver.exchangers.coupling import (
    CoupledSizingSettings,
    ExchangerGeometrySearch,
)
from dada_solver.exchangers.optimization import (
    ContinuousGeometryBounds,
    ExchangerObjective,
)


@dataclass(frozen=True, slots=True)
class LoadedCoupledExchangerProblem:
    simulation: SimulationConfiguration
    cold_search: ExchangerGeometrySearch
    hot_search: ExchangerGeometrySearch
    settings: CoupledSizingSettings
    example_data: bool


def load_coupled_exchanger_problem(
    path: str | Path,
) -> LoadedCoupledExchangerProblem:
    """Load references to one cycle and two independent exchanger searches."""

    problem_path = Path(path)
    with problem_path.open("rb") as stream:
        data = tomllib.load(stream)
    try:
        problem = data["problem"]
        simulation_path = _relative_path(
            problem_path, str(problem["simulation_configuration"])
        )
        cold_path = _relative_path(
            problem_path, str(problem["cold_exchanger_configuration"])
        )
        hot_path = _relative_path(
            problem_path, str(problem["hot_exchanger_configuration"])
        )
        simulation = load_simulation_configuration(simulation_path)
        cold = load_exchanger_screening_configuration(cold_path)
        hot = load_exchanger_screening_configuration(hot_path)
        _require_same_gas(simulation, cold, "cold")
        _require_same_gas(simulation, hot, "hot")
        settings_data = data["coupling"]
        settings = CoupledSizingSettings(
            maximum_iterations=int(settings_data["maximum_iterations"]),
            relative_tolerance=float(settings_data["relative_tolerance"]),
            flow_under_relaxation=float(settings_data["flow_under_relaxation"]),
            exchanger_maximum_iterations=int(
                settings_data["exchanger_maximum_iterations"]
            ),
            exchanger_function_tolerance=float(
                settings_data["exchanger_function_tolerance"]
            ),
            hydraulic_continuation_steps=int(
                settings_data.get("hydraulic_continuation_steps", 1)
            ),
        )
        objective_data = data.get("objectives", {})
        partition = data["hydraulic_partition"]
        cold_objective = ExchangerObjective(
            str(objective_data.get("cold", "minimum_gas_volume"))
        )
        hot_objective = ExchangerObjective(
            str(objective_data.get("hot", "minimum_gas_volume"))
        )
        return LoadedCoupledExchangerProblem(
            simulation=simulation,
            cold_search=_search(
                cold,
                cold_objective,
                float(partition["cold_inlet_core_fraction"]),
                float(partition["cold_inlet_collector_loss_coefficient"]),
                float(partition["cold_outlet_collector_loss_coefficient"]),
            ),
            hot_search=_search(
                hot,
                hot_objective,
                float(partition["hot_inlet_core_fraction"]),
                float(partition["hot_inlet_collector_loss_coefficient"]),
                float(partition["hot_outlet_collector_loss_coefficient"]),
            ),
            settings=settings,
            example_data=bool(data["metadata"]["example_data"]),
        )
    except KeyError as error:
        raise ValueError(
            f"Missing coupled exchanger configuration field: {error.args[0]}"
        ) from error


def _search(
    config: ExchangerScreeningConfiguration,
    objective: ExchangerObjective,
    inlet_fraction: float,
    inlet_collector_loss: float,
    outlet_collector_loss: float,
) -> ExchangerGeometrySearch:
    return ExchangerGeometrySearch(
        channel_counts=config.channel_counts,
        bounds=ContinuousGeometryBounds(
            min(config.channel_widths),
            max(config.channel_widths),
            min(config.channel_heights),
            max(config.channel_heights),
            min(config.channel_lengths),
            max(config.channel_lengths),
        ),
        seed_widths=config.channel_widths,
        seed_heights=config.channel_heights,
        seed_lengths=config.channel_lengths,
        transport=config.transport,
        thermal_resistances=config.thermal_resistances,
        requirements=config.requirements,
        inlet_core_resistance_fraction=inlet_fraction,
        inlet_collector_loss_coefficient=inlet_collector_loss,
        outlet_collector_loss_coefficient=outlet_collector_loss,
        objective=objective,
        surface_roughness=config.surface_roughness,
        minor_loss_coefficient=config.minor_loss_coefficient,
    )


def _relative_path(parent: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else parent.parent / path


def _require_same_gas(
    simulation: SimulationConfiguration,
    exchanger: ExchangerScreeningConfiguration,
    side: str,
) -> None:
    if exchanger.gas != simulation.gas:
        raise ValueError(
            f"The {side} exchanger working gas must exactly match the cycle gas."
        )
