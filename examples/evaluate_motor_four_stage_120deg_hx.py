"""Evaluate the first 120/60/120/60 motor sizing point with a redesigned H_o.

This is a controlled 2x2 comparison:

    A. current H_i/H_o + retained champion timing
    B. current H_i/H_o + 120/60/120/60 timing
    C. redesigned H_o   + retained champion timing
    D. redesigned H_o   + 120/60/120/60 timing

The cylinder intermediate volume levels remain those of the retained 7D champion.
Only the shared timing knots change in the target-timing cases.

Redesigned H_o:
    - 6050 tubes
    - 0.33 mm tube ID (unchanged)
    - retained K2 H_o tube length (~25.19 mm)
    - wall thickness, pitch and header depth unchanged
    - outlet-valve CdA multiplied by 194.2288/120 so the valve is not sized
      around the old 194-degree low-pressure exchange.

H_i remains the retained K2 exchanger.

Every case uses air, 2 Hz, and a uniform 100 kPa absolute filling condition at
theta=0. Gas inventory is therefore derived independently for each physical
geometry; it is not held fixed when H_o dead volume changes.

Outputs:
    outputs/motor_four_stage_120deg_hx/report.json

Run:
    PYTHONPATH=src python3 examples/evaluate_motor_four_stage_120deg_hx.py
"""

from __future__ import annotations

from dataclasses import replace
import json
import math
from pathlib import Path
import time

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.configuration import ChargeConfiguration
from dada_solver.four_stage_kinematics import FourStageVolumeKinematics

from compare_motor_motion_laws_stage7A5 import (
    A5_CAMPAIGN,
    ROOT,
    _candidate_design_and_mass,
    _evaluate,
)
from optimize_motor_exchanger_asymmetry_stageK1 import (
    _hardware_metrics,
    _scaled_exchanger,
)
from optimize_motor_piecewise_stageP3 import _feasibility


THERMO3D_REPORT = ROOT / "outputs" / "motor_four_stage_thermo3d" / "report.json"
OUTPUT_DIRECTORY = ROOT / "outputs" / "motor_four_stage_120deg_hx"
OUTPUT_REPORT = OUTPUT_DIRECTORY / "report.json"

CHARGE_PRESSURE_PA = 100_000.0
NEW_HO_TUBE_COUNT = 6050

TARGET_TIMING = {
    "t1": 120.0 / 360.0,
    "t2": 180.0 / 360.0,
    "t3": 300.0 / 360.0,
}


def _load_control() -> dict:
    report = json.loads(THERMO3D_REPORT.read_text())
    control = report.get("control_1bar")
    if control is None or not control.get("feasible"):
        raise RuntimeError("Missing feasible 1-bar thermo-3D control.")
    return control


def _geometry_metrics(exchanger) -> dict:
    built = exchanger.build()
    metadata = dict(built.metadata)
    dimensions = exchanger.bank.dimensions()
    return {
        "tube_count": exchanger.bank.tube_count,
        "tube_length_m": exchanger.bank.tube_length_m,
        "inner_diameter_m": exchanger.bank.inner_diameter_m,
        "wall_thickness_m": exchanger.bank.wall_thickness_m,
        "pitch_m": exchanger.bank.pitch_m,
        "header_depth_m": exchanger.bank.header_depth_m,
        "outlet_valve_cda_m2": exchanger.outlet_valve_cda_m2,
        "tube_flow_area_m2": dimensions["tube_flow_area_m2"],
        "tube_internal_area_m2": dimensions["tube_internal_area_m2"],
        "tube_gas_volume_m3": dimensions["tube_gas_volume_m3"],
        "header_gas_volume_m3": dimensions["header_gas_volume_m3"],
        "working_gas_volume_m3": dimensions["working_gas_volume_m3"],
        "core_width_m": dimensions["core_width_m"],
        "core_height_m": dimensions["core_height_m"],
        "fluid_envelope_length_m": dimensions["fluid_envelope_length_m"],
        "wall_capacity_j_k": float(metadata["wall_capacity_j_k"]),
        "overall_static_conductance_w_k": float(
            metadata["overall_static_conductance_w_k"]
        ),
        "gas_film_resistance_k_w": float(metadata["gas_film_resistance_k_w"]),
        "air_film_resistance_k_w": float(metadata["air_film_resistance_k_w"]),
        "metal_resistance_k_w": float(metadata["metal_resistance_k_w"]),
        "air_velocity_m_s": float(metadata["air_velocity_m_s"]),
        "air_reynolds": float(metadata["air_reynolds"]),
        "air_pressure_drop_pa": float(metadata["air_pressure_drop_pa"]),
        "fan_electrical_power_w": float(metadata["fan_electrical_power_w"]),
    }


def _motion_parameters(control: dict, target: bool) -> dict:
    retained = control["motion_parameters"]
    result = {
        "a_l": float(retained["a_l"]),
        "b_l": float(retained["b_l"]),
        "a_s": float(retained["a_s"]),
        "b_s": float(retained["b_s"]),
    }
    if target:
        result.update(TARGET_TIMING)
    else:
        result.update(
            {
                "t1": float(retained["t1"]),
                "t2": float(retained["t2"]),
                "t3": float(retained["t3"]),
            }
        )
    return result


def _phase_degrees(motion: dict) -> dict:
    t1 = motion["t1"]
    t2 = motion["t2"]
    t3 = motion["t3"]
    return {
        "low_pressure_exchange_deg": 360.0 * t1,
        "compression_deg": 360.0 * (t2 - t1),
        "high_pressure_exchange_deg": 360.0 * (t3 - t2),
        "expansion_deg": 360.0 * (1.0 - t3),
    }


def _make_design(base, control, *, redesigned_ho: bool, target_timing: bool):
    heat_in = _scaled_exchanger(
        base.heat_in,
        float(control["parameters"]["k_i"]),
    )
    heat_out = _scaled_exchanger(
        base.heat_out,
        float(control["parameters"]["k_o"]),
    )

    old_low_pressure_deg = (
        360.0 * float(control["motion_parameters"]["t1"])
    )
    valve_scale = old_low_pressure_deg / 120.0

    if redesigned_ho:
        new_bank = replace(
            heat_out.bank,
            tube_count=NEW_HO_TUBE_COUNT,
        )
        heat_out = replace(
            heat_out,
            bank=new_bank,
            outlet_valve_cda_m2=(
                heat_out.outlet_valve_cda_m2 * valve_scale
            ),
        )

    charge = ChargeConfiguration(
        temperature=base.configuration.charge.temperature,
        pressure=CHARGE_PRESSURE_PA,
    )
    configuration = replace(base.configuration, charge=charge)

    motion = _motion_parameters(control, target_timing)
    limits = configuration.machine_volumes
    kinematics = FourStageVolumeKinematics(
        limits.small_cylinder,
        limits.large_cylinder,
        **motion,
    )

    design = replace(
        base,
        configuration=configuration,
        heat_in=heat_in,
        heat_out=heat_out,
        kinematics=kinematics,
    )
    return design, motion, heat_in, heat_out, valve_scale


def _run_case(label, design, motion, heat_in, heat_out, definition):
    before = time.monotonic()
    result, state = _evaluate(label, design, definition)
    feasible, reasons = (
        _feasibility(result)
        if result.get("status") == "converged"
        else (False, [])
    )
    return {
        "label": label,
        "motion_parameters": motion,
        "phase_degrees": _phase_degrees(motion),
        "H_i": _geometry_metrics(heat_in),
        "H_o": _geometry_metrics(heat_out),
        "result": result,
        "feasible": feasible,
        "physical_constraint_failures": reasons,
        "elapsed_seconds": time.monotonic() - before,
        "last_complete_state": None if state is None else state.tolist(),
    }


def main():
    control = _load_control()
    definition = CampaignDefinition(A5_CAMPAIGN)
    base, _legacy_mass, _saved, _best, _candidate = (
        _candidate_design_and_mass(definition)
    )

    scenarios = (
        ("A_current_hardware_retained_timing", False, False),
        ("B_current_hardware_target_120_60_120_60", False, True),
        ("C_new_Ho_retained_timing", True, False),
        ("D_new_Ho_target_120_60_120_60", True, True),
    )

    records = []
    valve_scale = None
    for label, redesigned_ho, target_timing in scenarios:
        design, motion, heat_in, heat_out, valve_scale = _make_design(
            base,
            control,
            redesigned_ho=redesigned_ho,
            target_timing=target_timing,
        )
        record = _run_case(
            label,
            design,
            motion,
            heat_in,
            heat_out,
            definition,
        )
        records.append(record)

        result = record["result"]
        eta = result.get("indicated_thermal_efficiency")
        power = result.get("indicated_power_w")
        qin = result.get("heat_input_w")
        eta_text = "n/a" if eta is None else f"{100.0*eta:.6f}%"
        p_text = "n/a" if power is None else f"{power:.3f} W"
        q_text = "n/a" if qin is None else f"{qin:.3f} W"
        print(
            f"{label:45s}  "
            f"{result.get('status'):>12s}  "
            f"eta={eta_text:>11s}  P={p_text:>10s}  Qin={q_text:>11s}"
        )

    baseline = records[0]
    eta0 = baseline["result"].get("indicated_thermal_efficiency")
    p0 = baseline["result"].get("indicated_power_w")

    for record in records:
        eta = record["result"].get("indicated_thermal_efficiency")
        power = record["result"].get("indicated_power_w")
        record["delta_vs_case_A"] = {
            "efficiency_percentage_point": (
                None if eta is None or eta0 is None
                else 100.0 * (eta - eta0)
            ),
            "power_w": (
                None if power is None or p0 is None
                else power - p0
            ),
        }

    report = {
        "experiment": (
            "First physical exchanger/timing screening for the provisional "
            "120/60/120/60 sizing cycle."
        ),
        "design_decisions": {
            "charge_pressure_pa": CHARGE_PRESSURE_PA,
            "frequency_hz": abs(base.configuration.angular_speed) / (2.0 * math.pi),
            "H_i": "Retained K2 geometry.",
            "H_o_new_tube_count": NEW_HO_TUBE_COUNT,
            "H_o_new_inner_diameter_m": float(
                records[2]["H_o"]["inner_diameter_m"]
            ),
            "H_o_new_tube_length_m": float(
                records[2]["H_o"]["tube_length_m"]
            ),
            "H_o_valve_cda_scale": valve_scale,
            "H_o_valve_reason": (
                "Scale retained H_o outlet-valve CdA by retained LP duration / "
                "120 deg so the valve is not sized around the old 194-deg exchange."
            ),
            "target_phase_degrees": {
                "low_pressure_exchange": 120.0,
                "compression": 60.0,
                "high_pressure_exchange": 120.0,
                "expansion": 60.0,
            },
            "volume_levels": {
                key: float(control["motion_parameters"][key])
                for key in ("a_l", "b_l", "a_s", "b_s")
            },
            "mass_policy": (
                "Uniform 100 kPa fill at theta=0 for each physical geometry; "
                "gas inventory changes when exchanger dead volume changes."
            ),
        },
        "reference_control_from_thermo3d": {
            "efficiency": control["result"]["indicated_thermal_efficiency"],
            "power_w": control["result"]["indicated_power_w"],
            "heat_input_w": control["result"]["heat_input_w"],
        },
        "cases": records,
    }

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    OUTPUT_REPORT.write_text(json.dumps(report, indent=2) + "\n")

    print("\nH_o geometry comparison")
    for record in (records[0], records[2]):
        ho = record["H_o"]
        print(
            f"{record['label']:45s} "
            f"N={ho['tube_count']:4d}  "
            f"Aflow={1e6*ho['tube_flow_area_m2']:.1f} mm^2  "
            f"UA={ho['overall_static_conductance_w_k']:.3f} W/K  "
            f"Vgas={1e6*ho['working_gas_volume_m3']:.2f} mL  "
            f"face={1000*ho['core_width_m']:.1f}x"
            f"{1000*ho['core_height_m']:.1f} mm  "
            f"CdA={1e6*ho['outlet_valve_cda_m2']:.1f} mm^2"
        )

    print(f"\nWrote {OUTPUT_REPORT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
