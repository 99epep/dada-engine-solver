"""Causal phase-heat ablation for the retained motor and redesigned H_o.

Purpose
-------
Measure whether gas <-> exchanger-wall heat transfer during compression and/or
expansion is beneficial or parasitic, WITHOUT changing kinematics, hydraulics,
geometry, charge policy, or external-air boundary conditions.

Two physical configurations are tested:
    A: retained K2 H_i + retained K2 H_o, retained champion timing
    C: retained K2 H_i + redesigned 6050-tube H_o, retained champion timing

The already-computed A and C results are copied from:
    outputs/motor_four_stage_120deg_hx/report.json

For each physical configuration, four ablations are evaluated:
    1. suppress H_o wall <-> gas heat only during compression;
    2. suppress both H_i and H_o wall <-> gas heat during compression;
    3. suppress both H_i and H_o wall <-> gas heat during expansion;
    4. suppress both H_i and H_o wall <-> gas heat during compression+expansion.

External air <-> wall heat exchange remains active at all times. When a
gas-wall path is suppressed, the corresponding wall still exchanges heat with
its external-air stream and stores/releases wall energy normally. Only the
wall <-> working-gas term is set to zero.

This is deliberately a diagnostic counterfactual, not proposed hardware and
not an optimized motor.

Run:
    PYTHONPATH=src python3 examples/analyze_motor_phase_heat_ablation.py

Output:
    outputs/motor_phase_heat_ablation/report.json
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import json
import math
from pathlib import Path
import time

import numpy as np

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.dynamics import ValveTopology
from dada_solver.exchangers.air_wall import AirWallMotor
from dada_solver.exchangers.wall_cycle import (
    convergence_summary,
    solve_periodic_wall_motor,
    wall_cycle_performance,
)
from dada_solver.heat_transfer import PrescribedHeatRate
from dada_solver.state import ThermodynamicState
from dada_solver.valves import ValveState

from compare_motor_motion_laws_stage7A5 import (
    A5_CAMPAIGN,
    ROOT,
    _candidate_design_and_mass,
)
import evaluate_motor_four_stage_120deg_hx as hx120


SOURCE_REPORT = ROOT / "outputs" / "motor_four_stage_120deg_hx" / "report.json"
OUTPUT_DIRECTORY = ROOT / "outputs" / "motor_phase_heat_ablation"
OUTPUT_REPORT = OUTPUT_DIRECTORY / "report.json"


@dataclass(frozen=True)
class PhaseHeatGateAirWallMotor(AirWallMotor):
    """AirWallMotor with selected wall-to-gas paths disabled by motor phase."""

    t1: float
    t2: float
    t3: float
    block_hi_compression: bool = False
    block_ho_compression: bool = False
    block_hi_expansion: bool = False
    block_ho_expansion: bool = False

    def _phase(self, angle: float) -> str:
        # The built motor model integrates forward motor-cycle angle from
        # 0 to 2*pi. This is the same convention used by the existing
        # phase-diagnostics script.
        u = (float(angle) / (2.0 * math.pi)) % 1.0
        if u < self.t1:
            return "low_pressure_exchange"
        if u < self.t2:
            return "compression"
        if u < self.t3:
            return "high_pressure_exchange"
        return "expansion"

    def derivative(self, angle, values):
        gas = ThermodynamicState.from_array(np.asarray(values[:8]))
        temperatures = gas.temperatures(self.model.gas)

        incoming = self.heat_in.rates(temperatures[2], values[8])
        outgoing = self.heat_out.rates(temperatures[3], values[9])

        hi_gas_heat = float(incoming["gas_heat_w"])
        ho_gas_heat = float(outgoing["gas_heat_w"])
        phase = self._phase(float(angle))

        if phase == "compression":
            if self.block_hi_compression:
                hi_gas_heat = 0.0
            if self.block_ho_compression:
                ho_gas_heat = 0.0
        elif phase == "expansion":
            if self.block_hi_expansion:
                hi_gas_heat = 0.0
            if self.block_ho_expansion:
                ho_gas_heat = 0.0

        instantaneous = replace(
            self.model,
            cold_heat_transfer=PrescribedHeatRate(hi_gas_heat),
            hot_heat_transfer=PrescribedHeatRate(ho_gas_heat),
        )
        rates = instantaneous.evaluate(
            angle,
            gas,
            ValveTopology(ValveState.CLOSED, ValveState.CLOSED),
        )

        # External-air <-> wall exchange remains unchanged. Only the gas-wall
        # leg is ablated, so wall energy follows Q_air_to_wall - Q_wall_to_gas.
        hi_wall_rate = float(incoming["air_heat_w"]) - hi_gas_heat
        ho_wall_rate = float(outgoing["air_heat_w"]) - ho_gas_heat

        return np.r_[
            rates.state_derivative,
            hi_wall_rate,
            ho_wall_rate,
            incoming["air_heat_w"],
            outgoing["air_heat_w"],
            hi_gas_heat,
            ho_gas_heat,
            rates.gas_work_rate,
        ] / self.model.angular_speed


def _build_physical_case(
    *,
    control: dict,
    base,
    redesigned_ho: bool,
):
    design, motion, heat_in, heat_out, _valve_scale = hx120._make_design(
        base,
        control,
        redesigned_ho=redesigned_ho,
        target_timing=False,
    )
    wrapper = design.build()
    if not isinstance(wrapper, AirWallMotor):
        raise TypeError("Expected AirWallMotor.")
    return design, motion, wrapper


def _gate_wrapper(
    wrapper: AirWallMotor,
    motion: dict,
    *,
    block_hi_compression=False,
    block_ho_compression=False,
    block_hi_expansion=False,
    block_ho_expansion=False,
):
    return PhaseHeatGateAirWallMotor(
        model=wrapper.model,
        heat_in=wrapper.heat_in,
        heat_out=wrapper.heat_out,
        t1=float(motion["t1"]),
        t2=float(motion["t2"]),
        t3=float(motion["t3"]),
        block_hi_compression=block_hi_compression,
        block_ho_compression=block_ho_compression,
        block_hi_expansion=block_hi_expansion,
        block_ho_expansion=block_ho_expansion,
    )


def _evaluate_gated(
    *,
    label: str,
    wrapper: PhaseHeatGateAirWallMotor,
    initial_state: np.ndarray,
    design,
    definition: CampaignDefinition,
):
    before = time.monotonic()
    periodic = solve_periodic_wall_motor(
        wrapper,
        initial_state,
        maximum_cycles=design.configuration.numerical.maximum_cycles,
        settings=definition.wall_numerical_settings,
    )

    result = {
        "label": label,
        "status": periodic.status,
        "message": periodic.message,
        "convergence": convergence_summary(periodic.history),
        "elapsed_seconds": time.monotonic() - before,
    }
    if not periodic.converged:
        result["last_complete_state"] = (
            None
            if periodic.last_complete_state is None
            else periodic.last_complete_state.tolist()
        )
        return result

    assert periodic.trajectory is not None
    performance = wall_cycle_performance(wrapper, periodic.trajectory)
    result.update(
        {
            "indicated_thermal_efficiency": performance.thermal_efficiency,
            "indicated_power_w": performance.gas_power,
            "heat_input_w": performance.heat_in_power,
            "heat_out_w": performance.heat_out_power,
            "gas_work_per_cycle_j": performance.gas_work_per_cycle,
            "heat_input_per_cycle_j": performance.heat_in_per_cycle,
            "heat_out_per_cycle_j": performance.heat_out_per_cycle,
            "last_complete_state": periodic.last_complete_state.tolist(),
        }
    )
    return result


def _baseline_summary(record: dict) -> dict:
    r = record["result"]
    return {
        "label": record["label"],
        "status": r["status"],
        "indicated_thermal_efficiency": r["indicated_thermal_efficiency"],
        "indicated_power_w": r["indicated_power_w"],
        "heat_input_w": r["heat_input_w"],
        "heat_out_w": r["heat_out_w"],
        "gas_work_per_cycle_j": r["indicated_power_w"] / 2.0,
        "heat_input_per_cycle_j": r["heat_input_w"] / 2.0,
        "heat_out_per_cycle_j": r["heat_out_w"] / 2.0,
        "last_complete_state": record["last_complete_state"],
    }


def _delta(result: dict, baseline: dict) -> dict:
    if result.get("indicated_thermal_efficiency") is None:
        return {}
    return {
        "efficiency_percentage_point": 100.0 * (
            result["indicated_thermal_efficiency"]
            - baseline["indicated_thermal_efficiency"]
        ),
        "power_w": result["indicated_power_w"] - baseline["indicated_power_w"],
        "heat_input_w": result["heat_input_w"] - baseline["heat_input_w"],
        "heat_out_w": result["heat_out_w"] - baseline["heat_out_w"],
        "gas_work_per_cycle_j": (
            result["gas_work_per_cycle_j"]
            - baseline["gas_work_per_cycle_j"]
        ),
        "heat_input_per_cycle_j": (
            result["heat_input_per_cycle_j"]
            - baseline["heat_input_per_cycle_j"]
        ),
    }


def main():
    if not SOURCE_REPORT.exists():
        raise FileNotFoundError(SOURCE_REPORT)

    source = json.loads(SOURCE_REPORT.read_text())
    by_label = {case["label"]: case for case in source["cases"]}

    baseline_A = _baseline_summary(
        by_label["A_current_hardware_retained_timing"]
    )
    baseline_C = _baseline_summary(
        by_label["C_new_Ho_retained_timing"]
    )

    control = hx120._load_control()
    definition = CampaignDefinition(A5_CAMPAIGN)
    base, _legacy_mass, _saved, _best, _candidate = (
        _candidate_design_and_mass(definition)
    )

    cases = {
        "A_current_hardware": {
            "redesigned_ho": False,
            "baseline": baseline_A,
        },
        "C_new_Ho": {
            "redesigned_ho": True,
            "baseline": baseline_C,
        },
    }

    ablations = (
        (
            "block_Ho_compression",
            dict(block_ho_compression=True),
        ),
        (
            "block_both_compression",
            dict(
                block_hi_compression=True,
                block_ho_compression=True,
            ),
        ),
        (
            "block_both_expansion",
            dict(
                block_hi_expansion=True,
                block_ho_expansion=True,
            ),
        ),
        (
            "block_both_compression_and_expansion",
            dict(
                block_hi_compression=True,
                block_ho_compression=True,
                block_hi_expansion=True,
                block_ho_expansion=True,
            ),
        ),
    )

    all_results = {}
    for case_name, case in cases.items():
        design, motion, physical_wrapper = _build_physical_case(
            control=control,
            base=base,
            redesigned_ho=case["redesigned_ho"],
        )
        baseline = case["baseline"]
        initial = np.asarray(baseline["last_complete_state"], dtype=float)

        case_results = {
            "baseline": baseline,
            "motion_parameters": motion,
            "ablations": [],
        }

        for ablation_name, flags in ablations:
            wrapper = _gate_wrapper(
                physical_wrapper,
                motion,
                **flags,
            )
            label = f"{case_name}__{ablation_name}"
            result = _evaluate_gated(
                label=label,
                wrapper=wrapper,
                initial_state=initial,
                design=design,
                definition=definition,
            )
            result["ablation"] = ablation_name
            result["gate_flags"] = flags
            result["delta_vs_same_hardware_baseline"] = _delta(
                result, baseline
            )
            case_results["ablations"].append(result)

            eta = result.get("indicated_thermal_efficiency")
            power = result.get("indicated_power_w")
            qin = result.get("heat_input_w")
            deta = result["delta_vs_same_hardware_baseline"].get(
                "efficiency_percentage_point"
            )
            eta_text = "n/a" if eta is None else f"{100.0*eta:.6f}%"
            p_text = "n/a" if power is None else f"{power:.3f} W"
            q_text = "n/a" if qin is None else f"{qin:.3f} W"
            d_text = "n/a" if deta is None else f"{deta:+.6f} pt"
            print(
                f"{label:58s} "
                f"{result['status']:>12s}  "
                f"eta={eta_text:>11s}  "
                f"P={p_text:>10s}  "
                f"Qin={q_text:>11s}  "
                f"deta={d_text}"
            )

        all_results[case_name] = case_results

    report = {
        "experiment": (
            "Causal gas-wall heat-transfer ablation during compression and "
            "expansion, with kinematics and hardware frozen."
        ),
        "interpretation_guardrail": (
            "These are counterfactual diagnostics only. No kinematic or "
            "exchanger optimization is performed, and the heat gates are not "
            "proposed physical hardware."
        ),
        "source_report": str(SOURCE_REPORT.relative_to(ROOT)),
        "phase_definition": (
            "Retained champion timing; phase boundaries are the retained "
            "t1/t2/t3 values for each physical case."
        ),
        "external_air_policy": (
            "External-air to wall heat exchange remains active during every "
            "phase; only selected wall-to-working-gas paths are suppressed."
        ),
        "cases": all_results,
    }

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    OUTPUT_REPORT.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nWrote {OUTPUT_REPORT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
