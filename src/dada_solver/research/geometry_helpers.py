'''Utility functions for research geometry calculations.

This module provides helpers that are used by research code to compute
derived geometric quantities from the underlying machine geometry.
'''

from __future__ import annotations

import math
from typing import Any, Mapping


def _bank_ratio(bank: Any) -> float:
    """Compute length / diameter for a single MicrotubeBank."""
    length = getattr(bank, "tube_length_m", None)
    if length is None:
        length = getattr(bank, "microtube_length", None)
    diameter = getattr(bank, "bundle_diameter_m", None)
    if diameter is None and hasattr(bank, "dimensions"):
        dims = bank.dimensions()
        diameter = dims.get("bundle_diameter_m")
        if diameter is None and "bundle_face_area_m2" in dims:
            diameter = math.sqrt(4.0 * dims["bundle_face_area_m2"] / math.pi)
        elif diameter is None and "core_width_m" in dims and "core_height_m" in dims:
            diameter = math.sqrt(4.0 * dims["core_width_m"] * dims["core_height_m"] / math.pi)
    if length is None or diameter is None:
        raise AttributeError(
            "Geometry object must have 'microtube_length' and 'microtube_diameter' attributes"
        )
    if not math.isfinite(diameter) or diameter <= 0:
        raise ValueError("microtube_diameter must be positive for ratio calculation")
    if not math.isfinite(length) or length < 0:
        raise ValueError("microtube_length must be nonnegative")
    return float(length / diameter)


def exchanger_length_to_diameter(geometry: Any, *, side: str | None = None) -> float:
    """Return the micro‑tube bundle length‑to‑diameter ratio.

    The ratio is defined as ``microtube_length / microtube_diameter`` and
    excludes any collector geometry. The function accepts:
    - Objects providing ``microtube_length`` and ``microtube_diameter`` directly.
    - :class:`~dada_solver.exchangers.microtube_geometry.MicrotubeBank` instances.
    - Exchanger models containing a ``.bank`` attribute.
    - Complete machine designs/bases containing ``heat_in`` and/or ``heat_out`` exchangers.
      When multiple exchangers are present and no ``side`` is chosen, the minimum ratio
      is returned (representing the critical constraint bound).
    - Mappings with matching keys.

    A ``ValueError`` is raised if the denominator is zero or negative.
    """
    if geometry is None:
        raise AttributeError("Geometry object cannot be None")

    # 1. Direct attribute access (e.g. mock objects or dedicated geometry wrappers)
    length = getattr(geometry, "microtube_length", None)
    diameter = getattr(geometry, "microtube_diameter", None)
    if length is not None and diameter is not None:
        if not math.isfinite(diameter) or diameter <= 0:
            raise ValueError("microtube_diameter must be positive for ratio calculation")
        if not math.isfinite(length) or length < 0:
            raise ValueError("microtube_length must be nonnegative")
        return float(length / diameter)

    # 2. Mapping access
    if isinstance(geometry, Mapping):
        if "microtube_length" in geometry and "microtube_diameter" in geometry:
            l_val = geometry["microtube_length"]
            d_val = geometry["microtube_diameter"]
            if not math.isfinite(d_val) or d_val <= 0:
                raise ValueError("microtube_diameter must be positive for ratio calculation")
            if not math.isfinite(l_val) or l_val < 0:
                raise ValueError("microtube_length must be nonnegative")
            return float(l_val / d_val)

    # 3. Direct MicrotubeBank instance or object with .tube_length_m and .dimensions()
    if hasattr(geometry, "tube_length_m"):
        return _bank_ratio(geometry)

    # 4. ExchangerModel with .bank attribute
    if hasattr(geometry, "bank") and geometry.bank is not None:
        return _bank_ratio(geometry.bank)

    # 5. MachineDesign or MachineBasis or namespace with heat_in / heat_out
    hi = getattr(geometry, "heat_in", None)
    ho = getattr(geometry, "heat_out", None)
    if side == "heat_in" and hi is not None:
        return _bank_ratio(getattr(hi, "bank", hi))
    if side == "heat_out" and ho is not None:
        return _bank_ratio(getattr(ho, "bank", ho))
    if hi is not None and ho is not None:
        r_in = _bank_ratio(getattr(hi, "bank", hi))
        r_out = _bank_ratio(getattr(ho, "bank", ho))
        return min(r_in, r_out)
    if hi is not None:
        return _bank_ratio(getattr(hi, "bank", hi))
    if ho is not None:
        return _bank_ratio(getattr(ho, "bank", ho))

    raise AttributeError(
        "Geometry object must have 'microtube_length' and 'microtube_diameter' attributes"
    )
