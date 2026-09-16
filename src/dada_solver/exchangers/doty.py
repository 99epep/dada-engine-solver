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


def validate_gas_hydraulics(path: str | Path) -> list[dict]:
    """Uncalibrated absolute Doty checks with declared property conventions.

    Reported pressure is treated as mean tube pressure; T=(T3+T4)/2.
    Experimental output uncertainty is not inflated to hide model discrepancy.
    Thermal measurements are retained, but no missing shell-side closure is
    invented to turn this tube-hydraulic test into a complete UA prediction.
    """
    from .gas_transport import DiluteGasTransport
    from .microtube_geometry import MicrotubeBank
    from .gas_correlations import MicrotubeGasModel, compressible_poiseuille
    import math
    with Path(path).open(newline='') as stream:
        rows=list(csv.DictReader(stream))
    bank=MicrotubeBank(309,.127,.00033,.0001524,.00125,.001)
    output=[]
    for index,row in enumerate(rows):
        tr=DiluteGasTransport(row['fluid'])
        t=(float(row['T3_K'])+float(row['T4_K']))/2
        p=float(row['pressure_Pa']); flow=float(row['mass_flow_kg_s'])
        measured=float(row['pressure_drop_Pa']);uncertainty=float(row['pressure_drop_uncertainty_Pa'])
        for model_id,mu in [('legacy_constant_300K',tr.viscosity(300)),('compressible_variable_transport',tr.viscosity(t))]:
            predicted=128*mu*bank.tube_length_m*flow*tr.gas_constant*t/(p*bank.tube_count*math.pi*bank.inner_diameter_m**4)
            recovered=compressible_poiseuille(p+predicted/2,p-predicted/2,t,mu,bank.tube_length_m,
                bank.inner_diameter_m,bank.tube_count,tr.gas_constant)
            diag=MicrotubeGasModel(tr).diagnose(bank,flow,p+predicted/2,p-predicted/2,t)
            output.append(dict(fluid=row['fluid'],row=index,model_id=model_id,
                quantity='tube_pressure_drop_Pa',measured=measured,predicted=predicted,
                absolute_error=abs(predicted-measured),signed_error=predicted-measured,
                relative_error=(predicted-measured)/measured,
                inside_experimental_uncertainty=abs(predicted-measured)<=uncertainty,
                experimental_output_uncertainty=uncertainty,
                model_validity=dict(thermal_and_hydraulic_screen=diag.model_validity,issues=diag.issues,
                    reynolds=diag.reynolds,knudsen=diag.knudsen,mach=diag.mach),
                inverse_flow_residual_kg_s=recovered-flow,
                source='Doty et al. 1991 Tables 1-2; '+str(path),
                property_provenance=tr.provenance,
                assumptions='Mean reported pressure; mean T3/T4; 309 tubes, 0.33 mm ID, 127 mm length; no fitted multiplier',
                thermal_validation=dict(measured_UA_W_K=float(row['UA_W_K']),measured_effectiveness=float(row['effectiveness']),
                    predicted_UA_W_K=None,inside_experimental_uncertainty=None,
                    status='unavailable_without_independent_shell_side_and_axial_temperature_model')))
    return output
