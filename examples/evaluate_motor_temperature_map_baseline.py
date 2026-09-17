"""Frozen-design temperature sweep for the final production-gas 9D champion.

The exact final periodic state of the selected 300 K-delta-T champion is reused
for the 598.15 K reference point. Neighboring temperature points are then warm
started sequentially from the preceding converged periodic state.

Cold source: 298.15 K.
Hot-source grid: 508.15, 538.15, 568.15, 598.15, 628.15, 658.15 K.

Only the hot external-air source condition changes. Geometry, motion, gas
inventory convention and external-air mass flow remain fixed. Hot-air density
is scaled inversely with absolute temperature at constant pressure.

Run:
    PYTHONPATH=src python3 examples/evaluate_motor_temperature_map_baseline.py

Output:
    outputs/motor_temperature_map/frozen_reference.json
"""

from __future__ import annotations

from dataclasses import replace
import argparse
import json
from pathlib import Path
import time

import numpy as np

from dada_solver.integration import IntegrationInterrupted

from compare_motor_motion_laws_stage7A5 import ROOT, _evaluate
from refine_motor_four_stage_hx9d_variable_gas import (
    CHARGE_PRESSURE_PA,
    _build_design,
    _load_basis,
    _uniform_state_and_hardware,
)


REFERENCE_REPORT = (
    ROOT / "outputs" / "motor_four_stage_hx9d_variable_gas" / "report.json"
)
DEFAULT_OUTPUT = (
    ROOT / "outputs" / "motor_temperature_map" / "frozen_reference.json"
)
COLD_SOURCE_K = 298.15
HOT_SOURCE_GRID_K = (508.15, 538.15, 568.15, 598.15, 628.15, 658.15)
REFERENCE_POWER_FLOOR_W = 40.0


def _load_reference() -> tuple[dict, np.ndarray]:
    report = json.loads(REFERENCE_REPORT.read_text())
    best = report.get("best_feasible")
    if best is None:
        raise RuntimeError("Missing best_feasible in final production-gas 9D report.")

    p = best["parameters"]
    parameters = {
        **{
            name: float(p[name])
            for name in ("t1", "t2", "t3", "a_l", "b_l", "a_s", "b_s")
        },
        "n_i": int(p["n_i"]),
        "n_o": int(p["n_o"]),
    }

    saved = best.get("last_complete_state")
    if saved is None:
        raise RuntimeError(
            "Final 9D champion has no saved periodic state; cannot use exact warm start."
        )
    state = np.asarray(saved, dtype=float)
    if state.shape != (10,):
        raise RuntimeError(f"Expected 10-state champion state, got {state.shape}.")
    return parameters, state


def _with_hot_source(design, hot_temperature_k: float):
    reference_inputs = design.heat_in.inputs
    density = (
        reference_inputs.air_density_kg_m3
        * reference_inputs.air_inlet_temperature_k
        / hot_temperature_k
    )
    hot_inputs = replace(
        reference_inputs,
        air_inlet_temperature_k=hot_temperature_k,
        air_density_kg_m3=density,
    )
    return replace(design, heat_in=replace(design.heat_in, inputs=hot_inputs))


def _summary(result: dict, hot_temperature_k: float, total_mass_kg: float) -> dict:
    eta = result.get("indicated_thermal_efficiency")
    eta_c = 1.0 - COLD_SOURCE_K / hot_temperature_k
    power = result.get("indicated_power_w")
    return {
        "hot_source_temperature_k": hot_temperature_k,
        "hot_source_temperature_deg_c": hot_temperature_k - 273.15,
        "cold_source_temperature_k": COLD_SOURCE_K,
        "source_delta_t_k": hot_temperature_k - COLD_SOURCE_K,
        "carnot_efficiency": eta_c,
        "indicated_thermal_efficiency": eta,
        "fraction_of_carnot": (
            float(eta) / eta_c if eta is not None and eta_c > 0.0 else None
        ),
        "indicated_power_w": power,
        "heat_input_w": result.get("heat_input_w"),
        "heat_out_w": result.get("heat_out_w"),
        "total_gas_mass_kg": total_mass_kg,
        "passes_reference_40w_guard": (
            power is not None and float(power) >= REFERENCE_POWER_FLOOR_W
        ),
        "model_validity": result.get("validity"),
        "maximum_tube_reynolds": result.get("maximum_tube_reynolds"),
        "maximum_tube_mach_number": result.get("maximum_tube_mach_number"),
        "maximum_absolute_mass_flow_kg_s": result.get(
            "maximum_absolute_mass_flow_kg_s"
        ),
        "diagnostics": result.get("diagnostics"),
        "status": result.get("status"),
        "message": result.get("message"),
    }


def _print_row(row: dict) -> None:
    eta = row["indicated_thermal_efficiency"]
    frac = row["fraction_of_carnot"]
    print(
        json.dumps(
            {
                "delta_t_k": row["source_delta_t_k"],
                "hot_k": row["hot_source_temperature_k"],
                "status": row["status"],
                "eta_pct": 100.0 * eta if eta is not None else None,
                "carnot_pct": 100.0 * row["carnot_efficiency"],
                "fraction_of_carnot_pct": (
                    100.0 * frac if frac is not None else None
                ),
                "power_w": row["indicated_power_w"],
                "heat_input_w": row["heat_input_w"],
                "passes_40w": row["passes_reference_40w_guard"],
                "message": row["message"],
            }
        ),
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidate-seconds",
        type=float,
        default=0.0,
        help=(
            "Optional per-temperature wall-clock limit. "
            "0 disables the limit (recommended for this six-point sweep)."
        ),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.candidate_seconds < 0.0:
        raise ValueError("--candidate-seconds must be non-negative.")

    if not REFERENCE_REPORT.exists():
        raise FileNotFoundError(
            f"Missing reference report: {REFERENCE_REPORT.relative_to(ROOT)}"
        )

    definition, base = _load_basis()
    parameters, champion_state = _load_reference()
    reference_design = _build_design(base, parameters)

    cold_inputs = replace(
        reference_design.heat_out.inputs,
        air_inlet_temperature_k=COLD_SOURCE_K,
    )
    reference_design = replace(
        reference_design,
        heat_out=replace(reference_design.heat_out, inputs=cold_inputs),
    )

    # Validate that the rebuilt geometry still matches the saved champion state.
    uniform_state, hardware = _uniform_state_and_hardware(reference_design)
    uniform_mass = float(np.sum(uniform_state[:8:2]))
    champion_mass = float(np.sum(champion_state[:8:2]))
    if not np.isclose(champion_mass, uniform_mass, rtol=1e-9, atol=1e-12):
        raise RuntimeError(
            "Rebuilt champion gas inventory does not match saved final 9D state: "
            f"{champion_mass:.12g} vs {uniform_mass:.12g} kg."
        )

    raw_results: dict[float, dict] = {}

    def run_point(hot_temperature_k: float, initial_state: np.ndarray) -> np.ndarray:
        design = _with_hot_source(reference_design, hot_temperature_k)

        progress_callback = None
        if args.candidate_seconds > 0.0:
            deadline = time.monotonic() + args.candidate_seconds

            def progress_callback(_: object) -> None:
                if time.monotonic() >= deadline:
                    raise IntegrationInterrupted(
                        "Temperature point wall-clock limit exhausted."
                    )

        print(
            f"Starting deltaT={hot_temperature_k-COLD_SOURCE_K:.0f} K "
            f"(Th={hot_temperature_k:.2f} K)...",
            flush=True,
        )

        try:
            result, state = _evaluate(
                f"temperature_map_frozen_{int(round(hot_temperature_k * 100))}",
                design,
                definition,
                initial_state=np.asarray(initial_state, dtype=float),
                progress_callback=progress_callback,
            )
        except IntegrationInterrupted as exc:
            result = {"status": "interrupted", "message": str(exc)}
            state = None
        except (ValueError, RuntimeError) as exc:
            result = {"status": "integration_failure", "message": str(exc)}
            state = None

        raw_results[hot_temperature_k] = result
        row = _summary(result, hot_temperature_k, champion_mass)
        _print_row(row)

        if state is None or result.get("status") != "converged":
            print(
                "No converged periodic state available for warm-start propagation; "
                "reusing previous state.",
                flush=True,
            )
            return np.asarray(initial_state, dtype=float)
        return np.asarray(state, dtype=float)

    # The 598.15 K point starts from the exact converged state saved by the final
    # 9D campaign, rather than from a uniform initial condition.
    state_300 = run_point(598.15, champion_state)

    previous = state_300
    for hot in (568.15, 538.15, 508.15):
        previous = run_point(hot, previous)

    previous = state_300
    for hot in (628.15, 658.15):
        previous = run_point(hot, previous)

    rows = [
        _summary(raw_results[hot], hot, champion_mass)
        for hot in HOT_SOURCE_GRID_K
    ]

    reference = json.loads(REFERENCE_REPORT.read_text())["best_feasible"]
    report = {
        "study": "frozen 300 K-delta-T champion across six source temperatures",
        "reference_report": str(REFERENCE_REPORT.relative_to(ROOT)),
        "reference_candidate_index": reference["index"],
        "reference_parameters": parameters,
        "reference_periodic_state_reused": True,
        "charge_pressure_pa": CHARGE_PRESSURE_PA,
        "charge_temperature_k": base.configuration.charge.temperature,
        "cold_source_temperature_k": COLD_SOURCE_K,
        "hot_air_density_policy": (
            "reference density scaled inversely with absolute hot-air temperature "
            "at constant pressure"
        ),
        "candidate_seconds": args.candidate_seconds,
        "reference_power_floor_w": REFERENCE_POWER_FLOOR_W,
        "reference_power_floor_role": (
            "reported guard only; does not suppress frozen-design points"
        ),
        "hardware": hardware,
        "results": rows,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Saved {args.output}", flush=True)


if __name__ == "__main__":
    main()
