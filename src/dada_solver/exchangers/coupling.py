"""Fixed-point coupling between periodic cycles and exchanger geometry sizing."""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
import math
from typing import Callable

from dada_solver.configuration import SimulationConfiguration
from dada_solver.exchangers.cycle_adapter import (
    ExchangerSide,
    conservative_port_flow_operating_point,
)
from dada_solver.exchangers.models import ThermalResistanceModel, TransportProperties
from dada_solver.exchangers.optimization import (
    ContinuousGeometryBounds,
    ExchangerObjective,
    ExchangerOptimizationResult,
    optimize_parallel_channel_geometry,
)
from dada_solver.exchangers.sizing import ExchangerRequirements
from dada_solver.geometry import MachineVolumes
from dada_solver.hydraulics import (
    BlendedHydraulicFlow,
    CompressibleOrifice,
    HydraulicNetworkModels,
    QuasiSteadyCompressibleDuct,
    RectangularDuct,
    SeriesDuctOrifice,
)
from dada_solver.sizing.evaluator import (
    DesignEvaluation,
    EvaluationStatus,
    evaluate_configuration,
)
from dada_solver.state import ThermodynamicState
from dada_solver.integration import SegmentIntegrationProgress


@dataclass(frozen=True, slots=True)
class ExchangerGeometrySearch:
    """All explicit physical and numerical inputs for one exchanger search."""

    channel_counts: tuple[int, ...]
    bounds: ContinuousGeometryBounds
    seed_widths: tuple[float, ...]
    seed_heights: tuple[float, ...]
    seed_lengths: tuple[float, ...]
    transport: TransportProperties
    thermal_resistances: ThermalResistanceModel
    requirements: ExchangerRequirements
    inlet_core_resistance_fraction: float
    inlet_collector_loss_coefficient: float
    outlet_collector_loss_coefficient: float
    objective: ExchangerObjective = ExchangerObjective.MINIMUM_GAS_VOLUME
    surface_roughness: float = 0.0
    minor_loss_coefficient: float = 0.0

    def __post_init__(self) -> None:
        if not 0.0 < self.inlet_core_resistance_fraction < 1.0:
            raise ValueError("Inlet core resistance fraction must lie in (0, 1).")
        for name, value in (
            ("Inlet collector loss coefficient", self.inlet_collector_loss_coefficient),
            ("Outlet collector loss coefficient", self.outlet_collector_loss_coefficient),
        ):
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and non-negative.")


@dataclass(frozen=True, slots=True)
class CoupledSizingSettings:
    """Termination settings for the cycle/exchanger fixed-point iteration."""

    maximum_iterations: int = 12
    relative_tolerance: float = 1.0e-3
    flow_under_relaxation: float = 0.5
    exchanger_maximum_iterations: int = 200
    exchanger_function_tolerance: float = 1.0e-9
    hydraulic_continuation_steps: int = 1

    def __post_init__(self) -> None:
        if self.maximum_iterations <= 0 or self.exchanger_maximum_iterations <= 0:
            raise ValueError("Iteration limits must be positive.")
        if not math.isfinite(self.relative_tolerance) or self.relative_tolerance <= 0.0:
            raise ValueError("Relative tolerance must be finite and positive.")
        if not 0.0 < self.flow_under_relaxation <= 1.0:
            raise ValueError("Flow under-relaxation must lie in (0, 1].")
        if (
            not math.isfinite(self.exchanger_function_tolerance)
            or self.exchanger_function_tolerance <= 0.0
        ):
            raise ValueError(
                "Exchanger function tolerance must be finite and positive."
            )
        if self.hydraulic_continuation_steps <= 0:
            raise ValueError("Hydraulic continuation step count must be positive.")


class CoupledSizingStatus(Enum):
    CONVERGED = "converged"
    MAXIMUM_ITERATIONS_REACHED = "maximum_iterations_reached"
    CYCLE_EVALUATION_FAILURE = "cycle_evaluation_failure"
    COLD_EXCHANGER_INFEASIBLE = "cold_exchanger_infeasible"
    HOT_EXCHANGER_INFEASIBLE = "hot_exchanger_infeasible"


@dataclass(frozen=True, slots=True)
class CoupledSizingIteration:
    iteration: int
    cold_screening_mass_flow: float
    hot_screening_mass_flow: float
    cold_ua: float
    hot_ua: float
    cold_gas_volume: float
    hot_gas_volume: float
    normalized_change: float


@dataclass(frozen=True, slots=True)
class CoupledExchangerSizingResult:
    status: CoupledSizingStatus
    message: str
    configuration: SimulationConfiguration
    evaluation: DesignEvaluation | None
    cold_optimization: ExchangerOptimizationResult | None
    hot_optimization: ExchangerOptimizationResult | None
    history: tuple[CoupledSizingIteration, ...]
    reference_evaluation: DesignEvaluation | None


CycleEvaluationFunction = Callable[[SimulationConfiguration], DesignEvaluation]


class CoupledExchangerSizer:
    """Iterate the periodic cycle and two geometry searches to a fixed point."""

    def __init__(
        self,
        cold_search: ExchangerGeometrySearch,
        hot_search: ExchangerGeometrySearch,
        settings: CoupledSizingSettings = CoupledSizingSettings(),
        cycle_evaluator: CycleEvaluationFunction | None = None,
        progress_callback: Callable[[str], None] | None = None,
    ) -> None:
        self.cold_search = cold_search
        self.hot_search = hot_search
        self.settings = settings
        self._cycle_evaluator = cycle_evaluator
        self._warm_state: ThermodynamicState | None = None
        self._warm_topology = None
        self._progress_callback = progress_callback

    def _progress(self, message: str) -> None:
        if self._progress_callback is not None:
            self._progress_callback(message)

    def _integration_progress(self, item: SegmentIntegrationProgress) -> None:
        self._progress(
            "integration_progress: "
            f"angle_deg={math.degrees(item.current_angle):.6f}, "
            f"segment_deg=[{math.degrees(item.start_angle):.6f}, "
            f"{math.degrees(item.target_angle):.6f}], "
            f"hot_to_small={item.hot_to_small_state.value}, "
            f"cold_to_large={item.cold_to_large_state.value}, "
            f"elapsed_s={item.elapsed_seconds:.3f}, "
            f"rhs={item.right_hand_side_evaluations}"
        )

    def _evaluate(self, configuration: SimulationConfiguration) -> DesignEvaluation:
        if self._cycle_evaluator is not None:
            return self._cycle_evaluator(configuration)
        evaluation = evaluate_configuration(
            configuration,
            self._warm_state,
            self._warm_topology,
            self._integration_progress if self._progress_callback is not None else None,
        )
        if evaluation.cycle is not None and evaluation.cycle.completed:
            self._warm_state = evaluation.cycle.final_state
            self._warm_topology = evaluation.cycle.final_topology
        return evaluation

    def solve(
        self, configuration: SimulationConfiguration
    ) -> CoupledExchangerSizingResult:
        self._warm_state = None
        self._warm_topology = None
        current = configuration
        history: list[CoupledSizingIteration] = []
        previous_flows: tuple[float, float] | None = None
        last_evaluation = None
        reference_evaluation = None
        cold_result = None
        hot_result = None
        prefetched_evaluation: DesignEvaluation | None = None

        for iteration in range(1, self.settings.maximum_iterations + 1):
            evaluation = prefetched_evaluation or self._evaluate(current)
            prefetched_evaluation = None
            last_evaluation = evaluation
            if reference_evaluation is None:
                reference_evaluation = evaluation
            if (
                evaluation.status is not EvaluationStatus.CONVERGED
                or evaluation.diagnostics is None
            ):
                return self._result(
                    CoupledSizingStatus.CYCLE_EVALUATION_FAILURE,
                    f"Periodic cycle evaluation failed with status {evaluation.status.value}.",
                    current,
                    evaluation,
                    cold_result,
                    hot_result,
                    history,
                    reference_evaluation,
                )

            diagnostics = evaluation.diagnostics
            cold_point = conservative_port_flow_operating_point(
                diagnostics,
                ExchangerSide.COLD,
                representative_pressure=diagnostics.pressure_extrema["C"].minimum,
                representative_temperature=diagnostics.temperature_extrema["C"].maximum,
            )
            hot_point = conservative_port_flow_operating_point(
                diagnostics,
                ExchangerSide.HOT,
                representative_pressure=diagnostics.pressure_extrema["H"].minimum,
                representative_temperature=diagnostics.temperature_extrema["H"].maximum,
            )
            flows = (cold_point.absolute_mass_flow_rate, hot_point.absolute_mass_flow_rate)
            flow_change = math.inf
            if previous_flows is not None:
                factor = self.settings.flow_under_relaxation
                flows = tuple(
                    old + factor * (new - old)
                    for old, new in zip(previous_flows, flows, strict=True)
                )
                flow_change = max(
                    abs(new - old) / max(abs(old), abs(new))
                    for old, new in zip(previous_flows, flows, strict=True)
                )
                cold_point = replace(cold_point, absolute_mass_flow_rate=flows[0])
                hot_point = replace(hot_point, absolute_mass_flow_rate=flows[1])
            previous_flows = flows

            cold_result = self._optimize(self.cold_search, cold_point, current)
            cold_result = self._select_hydraulically_feasible(
                cold_result,
                self.cold_search,
                cold_point,
                current.hydraulics.cold_to_large_valve_cda,
                current,
            )
            if not cold_result.success or cold_result.best is None:
                return self._result(
                    CoupledSizingStatus.COLD_EXCHANGER_INFEASIBLE,
                    cold_result.message,
                    current,
                    evaluation,
                    cold_result,
                    hot_result,
                    history,
                    reference_evaluation,
                )
            hot_result = self._optimize(self.hot_search, hot_point, current)
            hot_result = self._select_hydraulically_feasible(
                hot_result,
                self.hot_search,
                hot_point,
                current.hydraulics.hot_to_small_valve_cda,
                current,
            )
            if not hot_result.success or hot_result.best is None:
                return self._result(
                    CoupledSizingStatus.HOT_EXCHANGER_INFEASIBLE,
                    hot_result.message,
                    current,
                    evaluation,
                    cold_result,
                    hot_result,
                    history,
                    reference_evaluation,
                )

            cold_performance = cold_result.best.performance
            hot_performance = hot_result.best.performance
            assert cold_performance.overall_conductance is not None
            assert hot_performance.overall_conductance is not None
            updated = _replace_exchanger_properties(
                current,
                cold_ua=cold_performance.overall_conductance,
                hot_ua=hot_performance.overall_conductance,
                cold_volume=cold_performance.geometry.gas_volume,
                hot_volume=hot_performance.geometry.gas_volume,
                cold_geometry=cold_performance.geometry,
                hot_geometry=hot_performance.geometry,
                cold_search=self.cold_search,
                hot_search=self.hot_search,
            )
            change = max(_configuration_change(current, updated), flow_change)
            history.append(
                CoupledSizingIteration(
                    iteration=iteration,
                    cold_screening_mass_flow=flows[0],
                    hot_screening_mass_flow=flows[1],
                    cold_ua=cold_performance.overall_conductance,
                    hot_ua=hot_performance.overall_conductance,
                    cold_gas_volume=cold_performance.geometry.gas_volume,
                    hot_gas_volume=hot_performance.geometry.gas_volume,
                    normalized_change=change,
                )
            )
            if change <= self.settings.relative_tolerance:
                verified = self._evaluate(updated)
                if (
                    verified.status is not EvaluationStatus.CONVERGED
                    or verified.diagnostics is None
                ):
                    return self._result(
                        CoupledSizingStatus.CYCLE_EVALUATION_FAILURE,
                        "Final updated exchanger configuration failed periodic verification.",
                        updated,
                        verified,
                        cold_result,
                        hot_result,
                        history,
                        reference_evaluation,
                    )
                return self._result(
                    CoupledSizingStatus.CONVERGED,
                    "Cycle and exchanger UA/volume iteration converged.",
                    updated,
                    verified,
                    cold_result,
                    hot_result,
                    history,
                    reference_evaluation,
                )
            continuation_evaluation = None
            continuation_steps = (
                1
                if current.hydraulic_flow_models == updated.hydraulic_flow_models
                else self.settings.hydraulic_continuation_steps
            )
            for step in range(1, continuation_steps + 1):
                fraction = step / continuation_steps
                self._progress(
                    f"hydraulic_continuation: fraction={fraction:.6f}, "
                    f"step={step}/{continuation_steps}"
                )
                continuation_configuration = _interpolate_configuration(
                    current, updated, fraction
                )
                continuation_evaluation = self._evaluate(continuation_configuration)
                if continuation_evaluation.status is not EvaluationStatus.CONVERGED:
                    return self._result(
                        CoupledSizingStatus.CYCLE_EVALUATION_FAILURE,
                        "Hydraulic continuation failed at final-model fraction "
                        f"{fraction:.6f} with status "
                        f"{continuation_evaluation.status.value}.",
                        continuation_configuration,
                        continuation_evaluation,
                        cold_result,
                        hot_result,
                        history,
                        reference_evaluation,
                    )
            current = updated
            prefetched_evaluation = continuation_evaluation

        if prefetched_evaluation is not None:
            last_evaluation = prefetched_evaluation
        return self._result(
            CoupledSizingStatus.MAXIMUM_ITERATIONS_REACHED,
            "Maximum cycle/exchanger iteration count reached.",
            current,
            last_evaluation,
            cold_result,
            hot_result,
            history,
            reference_evaluation,
        )

    def _optimize(
        self,
        search: ExchangerGeometrySearch,
        operating_point,
        configuration: SimulationConfiguration,
    ) -> ExchangerOptimizationResult:
        return optimize_parallel_channel_geometry(
            channel_counts=search.channel_counts,
            bounds=search.bounds,
            seed_widths=search.seed_widths,
            seed_heights=search.seed_heights,
            seed_lengths=search.seed_lengths,
            operating_point=operating_point,
            gas=configuration.gas,
            transport=search.transport,
            thermal_resistances=search.thermal_resistances,
            requirements=search.requirements,
            objective=search.objective,
            surface_roughness=search.surface_roughness,
            minor_loss_coefficient=search.minor_loss_coefficient,
            maximum_iterations=self.settings.exchanger_maximum_iterations,
            function_tolerance=self.settings.exchanger_function_tolerance,
        )

    def _select_hydraulically_feasible(
        self, result, search, point, valve_cda, configuration
    ):
        candidates = []
        allowed_drop = min(
            search.requirements.maximum_pressure_drop,
            point.pressure * search.requirements.maximum_pressure_drop_fraction,
            point.pressure * (1.0 - 1.0e-12),
        )
        downstream = point.pressure - allowed_drop
        for candidate in result.feasible_candidates:
            inlet, outlet = _port_models(
                candidate.performance.geometry, search, valve_cda
            )
            capacities = (
                inlet.directed_flow(
                    point.pressure, downstream, point.temperature, configuration.gas
                ).mass_flow_rate,
                outlet.directed_flow(
                    point.pressure, downstream, point.temperature, configuration.gas
                ).mass_flow_rate,
            )
            if min(capacities) >= point.absolute_mass_flow_rate:
                candidates.append(candidate)
        if not candidates:
            return replace(
                result,
                success=False,
                message=(
                    "No thermally feasible candidate supplies the required flow "
                    "through both geometric ports within the pressure-drop limit."
                ),
                best=None,
                feasible_candidates=(),
                pareto_candidates=(),
            )
        best = min(
            candidates,
            key=lambda item: (
                item.performance.geometry.gas_volume
                if search.objective is ExchangerObjective.MINIMUM_GAS_VOLUME
                else item.performance.pressure_drop
                if search.objective is ExchangerObjective.MINIMUM_PRESSURE_DROP
                else item.performance.geometry.gas_side_area
            ),
        )
        return replace(result, best=best, feasible_candidates=tuple(candidates))

    @staticmethod
    def _result(
        status,
        message,
        configuration,
        evaluation,
        cold_optimization,
        hot_optimization,
        history,
        reference_evaluation,
    ) -> CoupledExchangerSizingResult:
        return CoupledExchangerSizingResult(
            status,
            message,
            configuration,
            evaluation,
            cold_optimization,
            hot_optimization,
            tuple(history),
            reference_evaluation,
        )


def _replace_exchanger_properties(
    configuration: SimulationConfiguration,
    *,
    cold_ua: float,
    hot_ua: float,
    cold_volume: float,
    hot_volume: float,
    cold_geometry,
    hot_geometry,
    cold_search: ExchangerGeometrySearch,
    hot_search: ExchangerGeometrySearch,
) -> SimulationConfiguration:
    old = configuration.machine_volumes
    volumes = MachineVolumes(
        small_cylinder=old.small_cylinder,
        large_cylinder=old.large_cylinder,
        cold_heat_exchanger=cold_volume,
        hot_heat_exchanger=hot_volume,
    )
    cold_inlet, cold_outlet = _port_models(
        cold_geometry, cold_search, configuration.hydraulics.cold_to_large_valve_cda
    )
    hot_inlet, hot_outlet = _port_models(
        hot_geometry, hot_search, configuration.hydraulics.hot_to_small_valve_cda
    )
    return replace(
        configuration,
        machine_volumes=volumes,
        cold_thermal_conductance=cold_ua,
        hot_thermal_conductance=hot_ua,
        hydraulic_flow_models=HydraulicNetworkModels(
            large_to_hot=hot_inlet,
            small_to_cold=cold_inlet,
            hot_to_small_valve=hot_outlet,
            cold_to_large_valve=cold_outlet,
        ),
    )


def _port_models(geometry, search, valve_cda):
    fraction = search.inlet_core_resistance_fraction

    def duct(part: float, collector_loss: float) -> QuasiSteadyCompressibleDuct:
        return QuasiSteadyCompressibleDuct(
            RectangularDuct(
                flow_area=geometry.total_flow_area,
                hydraulic_diameter=geometry.hydraulic_diameter,
                length=geometry.channel_length * part,
                aspect_ratio=geometry.aspect_ratio,
                dynamic_viscosity=search.transport.dynamic_viscosity,
                relative_roughness=(
                    geometry.surface_roughness / geometry.hydraulic_diameter
                ),
                minor_loss_coefficient=(
                    geometry.minor_loss_coefficient * part + collector_loss
                ),
            )
        )

    return (
        SeriesDuctOrifice(
            duct(fraction, search.inlet_collector_loss_coefficient)
        ),
        SeriesDuctOrifice(
            duct(1.0 - fraction, search.outlet_collector_loss_coefficient),
            CompressibleOrifice(valve_cda),
        ),
    )


def _configuration_change(
    old: SimulationConfiguration,
    new: SimulationConfiguration,
) -> float:
    pairs = (
        (old.cold_thermal_conductance, new.cold_thermal_conductance),
        (old.hot_thermal_conductance, new.hot_thermal_conductance),
        (old.machine_volumes.cold_heat_exchanger, new.machine_volumes.cold_heat_exchanger),
        (old.machine_volumes.hot_heat_exchanger, new.machine_volumes.hot_heat_exchanger),
    )
    return max(
        abs(updated - previous) / max(abs(previous), abs(updated))
        for previous, updated in pairs
    )


def _interpolate_configuration(
    initial: SimulationConfiguration,
    final: SimulationConfiguration,
    fraction: float,
) -> SimulationConfiguration:
    """Create one conservative homotopy state; fraction one is exactly final."""

    if fraction >= 1.0:
        return final
    if not 0.0 < fraction < 1.0:
        raise ValueError("Continuation fraction must lie in (0, 1].")
    initial_network = initial.hydraulic_flow_models or _reference_network(initial)
    final_network = final.hydraulic_flow_models or _reference_network(final)

    def blend(first: float, second: float) -> float:
        return first + fraction * (second - first)

    initial_volumes = initial.machine_volumes
    final_volumes = final.machine_volumes
    volumes = MachineVolumes(
        small_cylinder=final_volumes.small_cylinder,
        large_cylinder=final_volumes.large_cylinder,
        cold_heat_exchanger=blend(
            initial_volumes.cold_heat_exchanger,
            final_volumes.cold_heat_exchanger,
        ),
        hot_heat_exchanger=blend(
            initial_volumes.hot_heat_exchanger,
            final_volumes.hot_heat_exchanger,
        ),
    )

    def flow(first, second):
        return BlendedHydraulicFlow(first, second, fraction)

    return replace(
        final,
        machine_volumes=volumes,
        cold_thermal_conductance=blend(
            initial.cold_thermal_conductance, final.cold_thermal_conductance
        ),
        hot_thermal_conductance=blend(
            initial.hot_thermal_conductance, final.hot_thermal_conductance
        ),
        hydraulic_flow_models=HydraulicNetworkModels(
            large_to_hot=flow(initial_network.large_to_hot, final_network.large_to_hot),
            small_to_cold=flow(initial_network.small_to_cold, final_network.small_to_cold),
            hot_to_small_valve=flow(
                initial_network.hot_to_small_valve, final_network.hot_to_small_valve
            ),
            cold_to_large_valve=flow(
                initial_network.cold_to_large_valve, final_network.cold_to_large_valve
            ),
        ),
    )


def _reference_network(configuration: SimulationConfiguration) -> HydraulicNetworkModels:
    hydraulics = configuration.hydraulics
    return HydraulicNetworkModels(
        large_to_hot=CompressibleOrifice(hydraulics.large_to_hot_cda),
        small_to_cold=CompressibleOrifice(hydraulics.small_to_cold_cda),
        hot_to_small_valve=CompressibleOrifice(hydraulics.hot_to_small_valve_cda),
        cold_to_large_valve=CompressibleOrifice(hydraulics.cold_to_large_valve_cda),
    )
