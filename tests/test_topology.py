import math
from dataclasses import replace
import pytest

from dada_solver.integration import ValveEvent
from dada_solver.topology import CycleTopologyClassification, classify_cycle_topology
from dada_solver.valves import ValveState


def test_nominal_event_sequence_is_recognized() -> None:
    events = (
        ValveEvent(0.5, "cold_to_large", ValveState.OPEN, ValveState.CLOSED),
        ValveEvent(1.5, "hot_to_small", ValveState.CLOSED, ValveState.OPEN),
        ValveEvent(3.5, "hot_to_small", ValveState.OPEN, ValveState.CLOSED),
        ValveEvent(5.5, "cold_to_large", ValveState.CLOSED, ValveState.OPEN),
    )

    diagnostic = classify_cycle_topology(events)

    assert diagnostic.classification is CycleTopologyClassification.NOMINAL
    assert diagnostic.reasons == ()


def test_simultaneous_or_repeated_events_are_reported_not_rejected() -> None:
    events = (
        ValveEvent(0.5, "hot_to_small", ValveState.CLOSED, ValveState.OPEN),
        ValveEvent(0.5, "cold_to_large", ValveState.CLOSED, ValveState.OPEN),
        ValveEvent(1.5, "hot_to_small", ValveState.OPEN, ValveState.CLOSED),
        ValveEvent(2.0, "hot_to_small", ValveState.CLOSED, ValveState.OPEN),
    )

    diagnostic = classify_cycle_topology(events)

    assert diagnostic.classification is CycleTopologyClassification.NON_NOMINAL
    assert "simultaneous_valve_events_observed" in diagnostic.reasons
    assert "repeated_valve_transition_observed" in diagnostic.reasons


@pytest.mark.parametrize("offset", range(4))
def test_nominal_classification_is_invariant_to_cycle_origin(offset):
    events = (
        ValveEvent(0.5, "cold_to_large", ValveState.OPEN, ValveState.CLOSED),
        ValveEvent(1.5, "hot_to_small", ValveState.CLOSED, ValveState.OPEN),
        ValveEvent(3.5, "hot_to_small", ValveState.OPEN, ValveState.CLOSED),
        ValveEvent(5.5, "cold_to_large", ValveState.CLOSED, ValveState.OPEN),
    )
    rotated = events[offset:] + events[:offset]
    rotated = tuple(replace(event, angle=0.5+i) for i, event in enumerate(rotated))
    assert classify_cycle_topology(rotated).classification is CycleTopologyClassification.NOMINAL
