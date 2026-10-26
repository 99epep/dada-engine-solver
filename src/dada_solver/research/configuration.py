'''Configuration models for research constraints.

Only a minimal model for the new ``exchanger_length_to_diameter`` constraint is
provided here.  The existing study loading logic (in ``study_schema.py``) reads
the raw TOML entries and validates them against ``PHYSICAL_CONSTRAINTS`` from
``margins.py``.  Adding a Pydantic/Dataclass model makes the constraint
available for other parts of the codebase (e.g., UI editors) without affecting
runtime behaviour.
'''

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Any

from dada_solver.sizing.constraints import ConstraintValue
from .geometry_helpers import exchanger_length_to_diameter


@dataclass(frozen=True, slots=True)
class ConstraintEntry:
    """Base class for research constraints.

    ``type`` identifies the constraint name, and ``value`` holds the numeric
    value (for constraints that store a single number).  Sub‑classes can add
    further validation logic if required.
    """
    value: float
    type: str = "constraint"


# Specific constraint for the micro‑tube bundle length‑to‑diameter ratio.
@dataclass(frozen=True, slots=True)
class ExchangerLengthToDiameter(ConstraintEntry):
    value: float  # length / diameter, must be > 0
    type: Literal["exchanger_length_to_diameter"] = "exchanger_length_to_diameter"

    def __post_init__(self) -> None:
        if self.value <= 0:
            raise ValueError("exchanger_length_to_diameter must be greater than 0")

    @property
    def name(self) -> str:
        return self.type

    def evaluate(self, evaluation: Any) -> ConstraintValue:
        """Evaluate the length-to-diameter constraint against an evaluation or design."""
        target = getattr(evaluation, "design", None)
        if target is None:
            target = getattr(evaluation, "model", None) or evaluation
        try:
            ratio = exchanger_length_to_diameter(target)
            margin = float(ratio - self.value)
            return ConstraintValue(
                name=self.name,
                margin=margin,
                satisfied=bool(margin >= 0.0),
                available=True,
            )
        except Exception:
            return ConstraintValue(
                name=self.name,
                margin=None,
                satisfied=False,
                available=False,
            )
