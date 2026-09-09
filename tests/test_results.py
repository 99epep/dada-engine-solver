import pytest

from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.results import extract_cycle_diagnostics
from dada_solver.topology import CycleTopologyClassification
from tests.test_periodic import create_static_cycle


def test_static_cycle_diagnostics_report_uniform_state_and_non_nominal_topology(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    model, _, cycle = create_static_cycle(ideal_gas)

    diagnostics = extract_cycle_diagnostics(cycle, model)

    for extrema in diagnostics.pressure_extrema.values():
        assert extrema.minimum == pytest.approx(2.0e5)
        assert extrema.maximum == pytest.approx(2.0e5)
    assert diagnostics.valve_events == ()
    assert (
        diagnostics.topology.classification
        is CycleTopologyClassification.NON_NOMINAL
    )

