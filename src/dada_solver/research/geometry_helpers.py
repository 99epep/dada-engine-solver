'''Utility functions for research geometry calculations.

This module provides helpers that are used by research code to compute
derived geometric quantities from the underlying machine geometry.
''' 

from __future__ import annotations

from typing import Any

# The geometry model used by the research code is defined in
# ``src/dada_solver/machine_basis.py`` as ``MachineGeometry`` (or a similar
# dataclass).  Importing the full class here would create a circular import in
# some contexts, therefore the function accepts a generic ``Any`` and accesses
# the required attributes dynamically.

def exchanger_length_to_diameter(geometry: Any) -> float:
    """Return the micro‑tube bundle length‑to‑diameter ratio.

    The ratio is defined as ``microtube_length / microtube_diameter`` and
    excludes any collector geometry.  ``geometry`` is expected to provide the
    attributes ``microtube_length`` and ``microtube_diameter`` which are stored
    in the machine basis.  A ``ValueError`` is raised if the denominator is
    zero or negative.
    """
    length = getattr(geometry, "microtube_length", None)
    diameter = getattr(geometry, "microtube_diameter", None)
    if length is None or diameter is None:
        raise AttributeError(
            "Geometry object must have 'microtube_length' and 'microtube_diameter' attributes"
        )
    if diameter <= 0:
        raise ValueError("microtube_diameter must be positive for ratio calculation")
    return length / diameter
