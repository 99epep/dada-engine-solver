"""Evaluate only the new 150/30/150/30 timing cases.

This reuses exactly the two physical exchanger configurations from
evaluate_motor_four_stage_120deg_hx.py:

- current K2 H_i / current K2 H_o
- current K2 H_i / redesigned 6050-tube H_o

The redesigned H_o is NOT resized for the 150-degree exchange. In particular,
its outlet-valve CdA remains the one selected in the earlier 120-degree design
exercise. Only the four-stage timing changes.

Target chronology:
    low-pressure exchange   150 deg
    compression              30 deg
    high-pressure exchange  150 deg
    expansion                30 deg

The intermediate volume levels a_l, b_l, a_s, b_s remain those of the retained
7D champion.

For comparison, the report also copies the four already-computed cases from:
    outputs/motor_four_stage_120deg_hx/report.json

Run:
    PYTHONPATH=src python3 examples/evaluate_motor_four_stage_150_30_hx.py

Output:
    outputs/motor_four_stage_150_30_hx/report.json
"""

from __future__ import annotations

import json
from pathlib import Path

from dada_solver.campaign.definition import CampaignDefinition

from compare_motor_motion_laws_stage7A5 import (
    A5_CAMPAIGN,
    ROOT,
    _candidate_design_and_mass,
)

import evaluate_motor_four_stage_120deg_hx as previous


PREVIOUS_REPORT = (
    ROOT / "outputs" / "motor_four_stage_120deg_hx" / "report.json"
)
OUTPUT_DIRECTORY = ROOT / "outputs" / "motor_four_stage_150_30_hx"
OUTPUT_REPORT = OUTPUT_DIRECTORY / "report.json"

TARGET_TIMING = {
    "t1": 150.0 / 360.0,
    "t2": 180.0 / 360.0,
    "t3": 330.0 / 360.0,
}


def _summary(record: dict) -> dict:
    result = record["result"]
    return {
        "label": record["label"],
        "phase_degrees": record["phase_degrees"],
        "efficiency": result.get("indicated_thermal_efficiency"),
        "power_w": result.get("indicated_power_w"),
        "heat_input_w": result.get("heat_input_w"),
        "heat_out_w": result.get("heat_out_w"),
        "total_mass_kg": result.get("total_mass_kg"),
        "maximum_pressure_pa": result.get("maximum_pressure_pa"),
        "maximum_temperature_k": result.get("maximum_temperature_k"),
        "maximum_absolute_mass_flow_kg_s": result.get(
            "maximum_absolute_mass_flow_kg_s"
        ),
        "validity": result.get("validity"),
    }


def main():
    if not PREVIOUS_REPORT.exists():
        raise FileNotFoundError(PREVIOUS_REPORT)

    old_report = json.loads(PREVIOUS_REPORT.read_text())

    control = previous._load_control()
    definition = CampaignDefinition(A5_CAMPAIGN)
    base, _legacy_mass, _saved, _best, _candidate = (
        _candidate_design_and_mass(definition)
    )

    # Deliberately change timing only. The redesigned H_o construction,
    # including its 194.23/120 outlet-valve CdA scale, remains untouched.
    previous.TARGET_TIMING = dict(TARGET_TIMING)

    scenarios = (
        (
            "E_current_hardware_target_150_30_150_30",
            False,
        ),
        (
            "F_new_Ho_target_150_30_150_30",
            True,
        ),
    )

    records = []
    for label, redesigned_ho in scenarios:
        design, motion, heat_in, heat_out, _valve_scale = previous._make_design(
            base,
            control,
            redesigned_ho=redesigned_ho,
            target_timing=True,
        )
        record = previous._run_case(
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
        eta_text = "n/a" if eta is None else f"{100.0 * eta:.6f}%"
        p_text = "n/a" if power is None else f"{power:.3f} W"
        q_text = "n/a" if qin is None else f"{qin:.3f} W"
        print(
            f"{label:45s}  "
            f"{result.get('status'):>12s}  "
            f"eta={eta_text:>11s}  P={p_text:>10s}  Qin={q_text:>11s}"
        )

    old_by_label = {record["label"]: record for record in old_report["cases"]}
    comparison = {
        "A_current_hardware_retained_timing": _summary(
            old_by_label["A_current_hardware_retained_timing"]
        ),
        "B_current_hardware_120_60_120_60": _summary(
            old_by_label["B_current_hardware_target_120_60_120_60"]
        ),
        "C_new_Ho_retained_timing": _summary(
            old_by_label["C_new_Ho_retained_timing"]
        ),
        "D_new_Ho_120_60_120_60": _summary(
            old_by_label["D_new_Ho_target_120_60_120_60"]
        ),
        "E_current_hardware_150_30_150_30": _summary(records[0]),
        "F_new_Ho_150_30_150_30": _summary(records[1]),
    }

    # Explicit deltas answer the immediate hypothesis.
    b = old_by_label["B_current_hardware_target_120_60_120_60"]["result"]
    d = old_by_label["D_new_Ho_target_120_60_120_60"]["result"]
    e = records[0]["result"]
    f = records[1]["result"]

    deltas = {
        "current_hardware_150_30_vs_120_60": {
            "efficiency_percentage_point": 100.0 * (
                e["indicated_thermal_efficiency"]
                - b["indicated_thermal_efficiency"]
            ),
            "power_w": e["indicated_power_w"] - b["indicated_power_w"],
            "heat_input_w": e["heat_input_w"] - b["heat_input_w"],
        },
        "new_Ho_150_30_vs_120_60": {
            "efficiency_percentage_point": 100.0 * (
                f["indicated_thermal_efficiency"]
                - d["indicated_thermal_efficiency"]
            ),
            "power_w": f["indicated_power_w"] - d["indicated_power_w"],
            "heat_input_w": f["heat_input_w"] - d["heat_input_w"],
        },
    }

    report = {
        "experiment": (
            "Timing-only 150/30/150/30 test using exactly the two hardware "
            "configurations from the prior 120/60/120/60 comparison."
        ),
        "hypothesis": (
            "Shorter compression/expansion may reduce thermally parasitic "
            "exchange outside the two principal exchange phases."
        ),
        "target_phase_degrees": {
            "low_pressure_exchange": 150.0,
            "compression": 30.0,
            "high_pressure_exchange": 150.0,
            "expansion": 30.0,
        },
        "hardware_policy": (
            "No hardware resizing. The 6050-tube H_o and its outlet-valve CdA "
            "are identical to the previous 120-degree design case."
        ),
        "new_cases": records,
        "comparison": comparison,
        "deltas_vs_120_60_120_60": deltas,
    }

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    OUTPUT_REPORT.write_text(json.dumps(report, indent=2) + "\n")

    print("\nDelta 150/30/150/30 versus 120/60/120/60")
    for name, delta in deltas.items():
        print(
            f"{name:38s} "
            f"deta={delta['efficiency_percentage_point']:+.6f} point  "
            f"dP={delta['power_w']:+.3f} W  "
            f"dQin={delta['heat_input_w']:+.3f} W"
        )

    print(f"\nWrote {OUTPUT_REPORT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
