"""Explicit whole-bank scaling of Doty's steady fluid measurements.

This screening model is not a DADA thermal closure. UA is gas/gas overall
conductance; pressure drop is tube-side only. Property and pulse corrections
are caller-supplied scenarios, not experimentally established coefficients.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from math import isfinite
from pathlib import Path


@dataclass(frozen=True)
class DotyReference:
    mass_flow_kg_s: float
    ua_w_k: float
    tube_pressure_drop_pa: float
    pressure_pa: float
    tube_mean_temperature_k: float


def load_references(path: str | Path) -> tuple[DotyReference, ...]:
    """Load measured anchors from homogeneous nitrogen or helium reference file.
    
    Raises ValueError if the file is empty, contains mixed fluids, or
    lists an unsupported fluid.
    """
    with Path(path).open(newline="") as stream:
        rows = tuple(csv.DictReader(stream))
    
    if not rows:
        raise ValueError("Expected non-empty reference measurements.")
    
    fluids = {row["fluid"] for row in rows}
    if len(fluids) > 1:
        raise ValueError("Expected homogeneous fluid (all nitrogen or all helium).")
    
    fluid = fluids.pop()
    if fluid not in ("nitrogen", "helium"):
        raise ValueError(f"Unsupported fluid: {fluid}. Expected nitrogen or helium.")
    
    return tuple(DotyReference(
        float(row["mass_flow_kg_s"]), float(row["UA_W_K"]),
        float(row["pressure_drop_Pa"]), float(row["pressure_Pa"]),
        (float(row["T3_K"]) + float(row["T4_K"])) / 2,
    ) for row in rows)


def scale_bank(reference: DotyReference, *, parallel_banks: int,
               total_mass_flow_kg_s: float, ua_flow_exponent: float,
               ua_condition_factor: float, viscosity_ratio: float,
               density_ratio: float, additional_pressure_drop_pa: float) -> dict:
    """Scale identical banks with equal flow splitting on both gas sides.

    UA = N UA_ref (mass_flow / (N mass_flow_ref))**exponent * factor.
    Tube pressure loss follows anchored laminar mu*mass_flow/rho scaling.
    Density ratio refers to representative tube-side density, not pressure alone.
    Zero flow is excluded: steady effectiveness cannot predict idle heat transfer.
    """
    positive = (reference.mass_flow_kg_s, reference.ua_w_k,
                reference.tube_pressure_drop_pa, reference.pressure_pa,
                reference.tube_mean_temperature_k, total_mass_flow_kg_s,
                ua_condition_factor, viscosity_ratio, density_ratio)
    if (type(parallel_banks) is not int or parallel_banks < 1
            or any(not isfinite(v) or v <= 0 for v in positive)
            or not isfinite(ua_flow_exponent) or not 0 <= ua_flow_exponent <= 1
            or not isfinite(additional_pressure_drop_pa)
            or additional_pressure_drop_pa < 0):
        raise ValueError("Expected positive finite inputs, integer bank count, and exponent in [0, 1].")
    ratio = total_mass_flow_kg_s / (parallel_banks * reference.mass_flow_kg_s)
    return {
        "parallel_banks": parallel_banks,
        "per_bank_mass_flow_ratio": ratio,
        "estimated_ua_w_k": parallel_banks * reference.ua_w_k * ratio**ua_flow_exponent * ua_condition_factor,
        "estimated_tube_pressure_drop_pa": reference.tube_pressure_drop_pa * ratio * viscosity_ratio / density_ratio,
        "additional_pressure_drop_pa": additional_pressure_drop_pa,
        "estimated_tube_and_additional_pressure_drop_pa": reference.tube_pressure_drop_pa * ratio * viscosity_ratio / density_ratio + additional_pressure_drop_pa,
        "sum_bank_envelope_m3": parallel_banks * 0.190 * 0.030 * 0.035,
        "sum_bank_mass_kg": parallel_banks * 0.3,
        "status": "unvalidated_steady_scaling",
    }
