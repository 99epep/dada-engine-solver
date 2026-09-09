import numpy as np

from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_initial_state, build_model, initial_valve_topology
from dada_solver.integration import CycleIntegrator


def test_conserved_state_is_exactly_continuous_at_nonzero_valve_events() -> None:
    configuration = load_simulation_configuration(
        "examples/harmonic_controlled_example.toml"
    )
    model = build_model(configuration)
    initial_state = build_initial_state(configuration, model)
    cycle = CycleIntegrator(
        model,
        relative_tolerance=configuration.numerical.integration_relative_tolerance,
        absolute_tolerance=configuration.numerical.integration_absolute_tolerance,
        maximum_step_angle=0.02,
    ).integrate_cycle(initial_state, initial_valve_topology())

    assert cycle.completed
    nonzero_events = [event for event in cycle.events if event.angle > 0.0]
    assert nonzero_events
    for event in nonzero_events:
        indices = np.flatnonzero(cycle.angles == event.angle)
        assert indices.size == 2
        before, after = indices
        assert cycle.topologies[before] != cycle.topologies[after]
        np.testing.assert_array_equal(
            cycle.states[:, before], cycle.states[:, after]
        )
        assert cycle.cold_heat[before] == cycle.cold_heat[after]
        assert cycle.hot_heat[before] == cycle.hot_heat[after]
        assert cycle.gas_work[before] == cycle.gas_work[after]
