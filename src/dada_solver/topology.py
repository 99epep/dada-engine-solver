"""Post-processing diagnostics for the observed passive-valve topology."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math

from dada_solver.integration import ValveEvent
from dada_solver.valves import ValveState


class CycleTopologyClassification(Enum):
    NOMINAL = "nominal"
    NON_NOMINAL = "non_nominal"


@dataclass(frozen=True, slots=True)
class CycleTopologyDiagnostic:
    classification: CycleTopologyClassification
    reasons: tuple[str, ...]


def classify_cycle_topology(
    events: tuple[ValveEvent, ...],
    simultaneous_angle_tolerance: float = 1.0e-10,
) -> CycleTopologyDiagnostic:
    """Classify an observed event sequence without constraining it."""

    expected = (
        ("cold_to_large", ValveState.CLOSED),
        ("hot_to_small", ValveState.OPEN),
        ("hot_to_small", ValveState.CLOSED),
        ("cold_to_large", ValveState.OPEN),
    )
    observed = tuple((event.valve_name, event.current_state) for event in events)
    reasons: list[str] = []
    # A periodic orbit has no preferred first event. Changing its angular
    # origin (including reversing the prescribed kinematics) may rotate this
    # sequence, but must not permit overlapping or repeated valve openings.
    if not any(observed == expected[i:] + expected[:i] for i in range(len(expected))):
        reasons.append("observed_event_sequence_differs_from_nominal_sequence")

    for first, second in zip(events, events[1:]):
        if math.isclose(
            first.angle,
            second.angle,
            rel_tol=0.0,
            abs_tol=simultaneous_angle_tolerance,
        ):
            reasons.append("simultaneous_valve_events_observed")
            break

    counts = {
        (name, state): observed.count((name, state))
        for name in ("hot_to_small", "cold_to_large")
        for state in (ValveState.OPEN, ValveState.CLOSED)
    }
    if any(count > 1 for count in counts.values()):
        reasons.append("repeated_valve_transition_observed")

    if reasons:
        return CycleTopologyDiagnostic(
            CycleTopologyClassification.NON_NOMINAL,
            tuple(reasons),
        )
    return CycleTopologyDiagnostic(CycleTopologyClassification.NOMINAL, ())
