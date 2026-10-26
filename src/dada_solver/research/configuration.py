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
from typing import Literal


@dataclass(frozen=True, slots=True)
class ConstraintEntry:
    """Base class for research constraints.

    ``type`` identifies the constraint name, and ``value`` holds the numeric
    value (for constraints that store a single number).  Sub‑classes can add
    further validation logic if required.
    """
    type: str
    value: float


# Specific constraint for the micro‑tube bundle length‑to‑diameter ratio.
@dataclass(frozen=True, slots=True)
class ExchangerLengthToDiameter(ConstraintEntry):
    type: Literal["exchanger_length_to_diameter"] = "exchanger_length_to_diameter"
    value: float  # length / diameter, must be > 0

    def __post_init__(self) -> None:
        if self.value <= 0:
            raise ValueError("exchanger_length_to_diameter must be greater than 0")
