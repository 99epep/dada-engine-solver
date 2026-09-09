from pathlib import Path

import numpy as np
import pytest

from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_model
from dada_solver.exchangers.coupling import (
    CoupledExchangerSizer,
    CoupledSizingSettings,
    CoupledSizingStatus,
    ExchangerGeometrySearch,
)
from dada_solver.exchangers.models import ThermalResistanceModel, TransportProperties
from dada_solver.exchangers.optimization import ContinuousGeometryBounds
from dada_solver.exchangers.sizing import ExchangerRequirements
from dada_solver.results import CycleDiagnostics, Extrema
from dada_solver.sizing.design import DesignPoint
from dada_solver.sizing.evaluator import DesignEvaluation, EvaluationStatus
from dada_solver.sizing.evaluator import _rescale_state_inventory
from dada_solver.state import UniformCharge
from dada_solver.topology import (
    CycleTopologyClassification,
    CycleTopologyDiagnostic,
)


def _base_configuration():
    return load_simulation_configuration(
        Path(__file__).parents[1] / "examples" / "harmonic_controlled_example.toml"
    )


def _diagnostics() -> CycleDiagnostics:
    return CycleDiagnostics(
        pressure_extrema={
            name: Extrema(1.0e5, 1.2e5) for name in ("S", "L", "C", "H")
        },
        temperature_extrema={
            name: Extrema(280.0, 300.0) for name in ("S", "L", "C", "H")
        },
        mass_flow_extrema={
            "small_to_cold": Extrema(-5.0e-5, 1.0e-4),
            "cold_to_large": Extrema(0.0, 8.0e-5),
            "large_to_hot": Extrema(-5.0e-5, 8.0e-5),
            "hot_to_small": Extrema(0.0, 1.0e-4),
        },
        cold_heat_rate_extrema=Extrema(-1.0, 1.0),
        hot_heat_rate_extrema=Extrema(-1.0, 1.0),
        valve_events=(),
        topology=CycleTopologyDiagnostic(CycleTopologyClassification.NOMINAL, ()),
    )


def _evaluation(configuration, status=EvaluationStatus.CONVERGED):
    return DesignEvaluation(
        point=DesignPoint({}),
        configuration=configuration,
        status=status,
        periodic=None,  # type: ignore[arg-type]
        performance=None,
        diagnostics=_diagnostics() if status is EvaluationStatus.CONVERGED else None,
        validity=None,
        model=None,
        cycle=None,
    )


def _search() -> ExchangerGeometrySearch:
    return ExchangerGeometrySearch(
        channel_counts=(20, 40),
        bounds=ContinuousGeometryBounds(0.01, 0.02, 0.001, 0.002, 0.1, 0.2),
        seed_widths=(0.01, 0.02),
        seed_heights=(0.001, 0.002),
        seed_lengths=(0.1, 0.2),
        transport=TransportProperties(1.8e-5, 0.026),
        thermal_resistances=ThermalResistanceModel(0.001, 200.0, 100.0),
        requirements=ExchangerRequirements(5.0, 1.0e4, 0.3, 1.0e-3, 0.1),
        inlet_core_resistance_fraction=0.5,
        inlet_collector_loss_coefficient=0.0,
        outlet_collector_loss_coefficient=0.0,
    )


def test_coupled_sizer_reinjects_geometry_until_fixed_point() -> None:
    calls = []

    def evaluate(configuration):
        calls.append(configuration)
        return _evaluation(configuration)

    result = CoupledExchangerSizer(
        _search(),
        _search(),
        CoupledSizingSettings(maximum_iterations=4, relative_tolerance=1.0e-8),
        cycle_evaluator=evaluate,
    ).solve(_base_configuration())

    assert result.status is CoupledSizingStatus.CONVERGED
    assert len(calls) == 3
    assert len(result.history) == 2
    assert result.cold_optimization is not None
    assert result.cold_optimization.best is not None
    cold = result.cold_optimization.best.performance
    assert result.configuration.machine_volumes.cold_heat_exchanger == cold.geometry.gas_volume
    assert result.configuration.cold_thermal_conductance == cold.overall_conductance


def test_coupled_sizer_propagates_cycle_failure() -> None:
    result = CoupledExchangerSizer(
        _search(),
        _search(),
        cycle_evaluator=lambda configuration: _evaluation(
            configuration, EvaluationStatus.NUMERICAL_FAILURE
        ),
    ).solve(_base_configuration())

    assert result.status is CoupledSizingStatus.CYCLE_EVALUATION_FAILURE
    assert not result.history


def test_warm_start_rescaling_preserves_specific_energy_and_target_mass() -> None:
    configuration = _base_configuration()
    gas = configuration.gas
    state = UniformCharge(1.0e5, 290.0).create_state(
        gas,
        build_model(configuration).volumes(0.0),
    )

    scaled = _rescale_state_inventory(state, 1.7 * state.total_mass)

    assert scaled.total_mass == pytest.approx(1.7 * state.total_mass)
    np.testing.assert_allclose(scaled.energies / scaled.masses, state.energies / state.masses)
