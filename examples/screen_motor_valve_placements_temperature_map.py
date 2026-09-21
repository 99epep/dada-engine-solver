#!/usr/bin/env python3
"""Replay UD/DU/UU valve placements on every available temperature-map champion.

The historical temperature-map champions use the DD topology
(downstream/downstream). Their saved results are used directly as the DD
reference; only the three alternative passive-valve placements are integrated.

Output:
    outputs/motor_valve_placement_temperature_map.json

Run:
    PYTHONPATH=src python3 examples/screen_motor_valve_placements_temperature_map.py
"""

from __future__ import annotations

from dataclasses import replace
import hashlib
import json
from pathlib import Path
import re
from types import SimpleNamespace
import time

import numpy as np

from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.wall_backend import WallBackendSettings

from compare_motor_motion_laws_stage7A5 import (
    A5_CAMPAIGN,
    ROOT,
    _evaluate,
)
from optimize_motor_temperature_point_fixed import (
    COLD_K,
    MAX_FLOW,
    MAX_P,
    MAX_T,
    base_geometry,
    build_design,
    feasibility,
)
from refine_motor_four_stage_hx9d_variable_gas import _load_basis


MAP_ROOT = ROOT / "outputs" / "motor_temperature_map"
OUTPUT = ROOT / "outputs" / "motor_valve_placement_temperature_map.json"

# Historical DD is read from the champion and is NOT recomputed.
ALTERNATIVE_PLACEMENTS = (
    ("UD", "upstream", "downstream"),
    ("DU", "downstream", "upstream"),
    ("UU", "upstream", "upstream"),
)

TOPOLOGY_PLACEMENTS = {
    "DD": ("downstream", "downstream"),
    "UD": ("upstream", "downstream"),
    "DU": ("downstream", "upstream"),
    "UU": ("upstream", "upstream"),
}


def discover_reports() -> list[tuple[float, Path]]:
    """Return temperature-map reports having an exploitable champion."""
    found: list[tuple[float, Path]] = []

    pattern = re.compile(r"^dT_([0-9]+(?:\.[0-9]+)?)K$")

    for path in MAP_ROOT.glob("dT_*K/report.json"):
        match = pattern.match(path.parent.name)
        if match is None:
            continue

        delta_t = float(match.group(1))

        try:
            report = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue

        champion = report.get("champion")
        if not champion:
            continue

        state = champion.get("last_complete_state")
        if state is None:
            continue

        values = np.asarray(state, dtype=float)
        if values.shape != (10,) or not np.all(np.isfinite(values)):
            continue

        if not champion.get("parameters"):
            continue

        found.append((delta_t, path))

    return sorted(found)


def extrema(result: dict, group: str) -> float | None:
    items = result.get("diagnostics", {}).get(group, {})
    if not isinstance(items, dict):
        return None

    values = []
    for item in items.values():
        try:
            values.append(float(item["maximum"]))
        except (KeyError, TypeError, ValueError):
            pass

    return max(values) if values else None


def result_max_pressure(result: dict) -> float | None:
    value = result.get("maximum_pressure_pa")
    return float(value) if value is not None else extrema(result, "pressure_extrema")


def result_max_temperature(result: dict) -> float | None:
    value = result.get("maximum_temperature_k")
    return float(value) if value is not None else extrema(
        result, "temperature_extrema"
    )


def result_max_abs_flow(result: dict) -> float | None:
    value = result.get("maximum_absolute_mass_flow_kg_s")
    if value is not None:
        return float(value)

    items = result.get("diagnostics", {}).get("mass_flow_extrema", {})
    if not isinstance(items, dict):
        return None

    values = []
    for item in items.values():
        try:
            values.extend(
                (
                    abs(float(item["minimum"])),
                    abs(float(item["maximum"])),
                )
            )
        except (KeyError, TypeError, ValueError):
            pass

    return max(values) if values else None


def safe_float(value):
    return float(value) if value is not None else None


def carnot_efficiency(delta_t_k: float) -> float:
    hot_k = COLD_K + delta_t_k
    return 1.0 - COLD_K / hot_k


def normalize_result_fields(result: dict) -> dict:
    validity = result.get("validity") or {}

    return dict(
        status=result.get("status"),
        indicated_thermal_efficiency=safe_float(
            result.get("indicated_thermal_efficiency")
        ),
        indicated_power_w=safe_float(result.get("indicated_power_w")),
        heat_input_w=safe_float(result.get("heat_input_w")),
        heat_out_w=safe_float(result.get("heat_out_w")),
        maximum_pressure_pa=result_max_pressure(result),
        maximum_temperature_k=result_max_temperature(result),
        maximum_absolute_mass_flow_kg_s=result_max_abs_flow(result),
        validity_verdict=validity.get("verdict"),
        failed_validity_criteria=list(validity.get("failed_criteria") or []),
        maximum_tube_reynolds=safe_float(result.get("maximum_tube_reynolds")),
        maximum_tube_mach=safe_float(
            result.get("maximum_tube_mach_number")
        ),
    )


def historical_dd(
    source: dict,
    delta_t_k: float,
    power_floor_w: float,
) -> dict:
    """Convert the saved historical DD champion to the common screen schema."""
    champion = source["champion"]
    result = champion["result"]

    try:
        feasible, reasons = feasibility(result, power_floor_w)
    except Exception:
        # A temperature-map champion was selected under these constraints.
        # Keep the historical feasibility marker if the helper cannot consume
        # an older report schema, but do not fabricate failed criteria.
        feasible = bool(champion.get("feasible", True))
        reasons = []

    normalized = normalize_result_fields(result)
    eta = normalized["indicated_thermal_efficiency"]
    carnot = carnot_efficiency(delta_t_k)

    return dict(
        topology_code="DD",
        heat_in_valve_placement="downstream",
        heat_out_valve_placement="downstream",
        evaluation_origin="historical_champion",
        converged=result.get("status") == "converged",
        feasible_under_original_constraints=bool(feasible),
        failed_constraints=list(reasons),
        fraction_of_carnot=(
            eta / carnot if eta is not None and carnot > 0.0 else None
        ),
        periodic_cycle_count=None,
        elapsed_seconds=None,
        requested_rhs_backend=None,
        actual_rhs_backend=None,
        fallback_count=None,
        fallback_reason=None,
        final_periodic_state=champion.get("last_complete_state"),
        **normalized,
    )


def delta_fields(row: dict, dd: dict) -> None:
    eta = row.get("indicated_thermal_efficiency")
    eta_dd = dd.get("indicated_thermal_efficiency")

    power = row.get("indicated_power_w")
    power_dd = dd.get("indicated_power_w")

    if eta is not None and eta_dd is not None:
        row["efficiency_absolute_gain"] = eta - eta_dd
        row["efficiency_gain_percentage_points"] = 100.0 * (eta - eta_dd)
        row["efficiency_relative_gain"] = (
            (eta / eta_dd - 1.0) if eta_dd != 0.0 else None
        )
    else:
        row["efficiency_absolute_gain"] = None
        row["efficiency_gain_percentage_points"] = None
        row["efficiency_relative_gain"] = None

    if power is not None and power_dd is not None:
        row["power_absolute_gain_w"] = power - power_dd
        row["power_relative_gain"] = (
            (power / power_dd - 1.0) if power_dd != 0.0 else None
        )
    else:
        row["power_absolute_gain_w"] = None
        row["power_relative_gain"] = None


def best_topology(
    results: dict[str, dict],
    *,
    feasible_only: bool,
) -> str | None:
    candidates = []

    for code, row in results.items():
        eta = row.get("indicated_thermal_efficiency")

        if not row.get("converged"):
            continue
        if eta is None:
            continue
        if feasible_only and not row.get("feasible_under_original_constraints"):
            continue

        candidates.append((float(eta), code))

    return max(candidates)[1] if candidates else None


def write_report(payload: dict) -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")


def main() -> None:
    reports = discover_reports()

    if not reports:
        raise RuntimeError(
            f"No usable temperature-map champions found below {MAP_ROOT}"
        )

    print(
        "Found champions:",
        ", ".join(f"{delta:g} K" for delta, _ in reports),
        flush=True,
    )

    # Same numerical settings used by the validated 300 K valve screen.
    campaign = CampaignDefinition(A5_CAMPAIGN)
    _, base = _load_basis()
    geometry = base_geometry(base)

    backend = WallBackendSettings("numba")

    definition = SimpleNamespace(
        wall_numerical_settings=campaign.wall_numerical_settings,
        wall_backend=backend,
        exact_kinematics_cache=True,
        shared_replay=True,
    )

    payload = dict(
        study=(
            "passive check-valve placement replay on frozen "
            "temperature-map champions"
        ),
        nominal_flow="S -> H_i -> L -> H_o -> S",
        cold_source_temperature_k=COLD_K,
        topology_definitions={
            code: {
                "heat_in_valve_placement": placements[0],
                "heat_out_valve_placement": placements[1],
            }
            for code, placements in TOPOLOGY_PLACEMENTS.items()
        },
        global_constraints=dict(
            maximum_pressure_pa=MAX_P,
            maximum_temperature_k=MAX_T,
            maximum_absolute_mass_flow_kg_s=MAX_FLOW,
            valid_thermodynamic_exchanger_model=True,
        ),
        temperatures=[],
        complete=False,
    )

    total_started = time.perf_counter()

    print()
    print(
        " ΔT  topo   conv  feasible      eta %      Δη pp"
        "       P W       ΔP W  cycles    time s",
        flush=True,
    )
    print("-" * 88, flush=True)

    for delta_t_k, source_path in reports:
        source = json.loads(source_path.read_text())
        champion = source["champion"]

        hot_k = COLD_K + delta_t_k
        power_floor_w = float(source["power_floor_w"])
        source_complete = bool(source.get("campaign_complete", False))
        source_sha = hashlib.sha256(source_path.read_bytes()).hexdigest()

        initial = np.asarray(
            champion["last_complete_state"],
            dtype=float,
        )

        parameters = champion["parameters"]

        dd = historical_dd(source, delta_t_k, power_floor_w)
        delta_fields(dd, dd)

        results: dict[str, dict] = {"DD": dd}

        print(
            f"{delta_t_k:4.0f}  DD    "
            f"{str(dd['converged']):>5}  "
            f"{str(dd['feasible_under_original_constraints']):>8}  "
            f"{100*dd['indicated_thermal_efficiency']:9.4f}  "
            f"{0.0:9.4f}  "
            f"{dd['indicated_power_w']:8.3f}  "
            f"{0.0:9.3f}  "
            f"{'-':>6}  {'-':>8}",
            flush=True,
        )

        for code, heat_in, heat_out in ALTERNATIVE_PLACEMENTS:
            design = build_design(
                base,
                geometry,
                parameters,
                hot_k,
            )

            design = replace(
                design,
                configuration=replace(
                    design.configuration,
                    heat_in_valve_placement=heat_in,
                    heat_out_valve_placement=heat_out,
                ),
            )

            captures = []
            started = time.perf_counter()

        try:
            result, state = _evaluate(
                f"temperature_map_{delta_t_k:g}K_valve_{code}",
                design,
                definition,
                initial_state=initial.copy(),
                periodic_observer=captures.append,
            )
        except MicrotubeDomainError as error:
            elapsed = time.perf_counter() - started

            row = dict(
                topology_code=code,
                heat_in_valve_placement=heat_in,
                heat_out_valve_placement=heat_out,
                evaluation_origin="valve_placement_replay",
                status="invalid_exchanger",
                converged=False,
                feasible_under_original_constraints=False,
                failed_constraints=["invalid_exchanger"],
                failure_reason=str(error),
                fraction_of_carnot=None,
                periodic_cycle_count=None,
                elapsed_seconds=elapsed,
                requested_rhs_backend="numba",
                actual_rhs_backend=None,
                fallback_count=None,
                fallback_reason=None,
                final_periodic_state=None,
                backend_statistics=None,
                indicated_thermal_efficiency=None,
                indicated_power_w=None,
                heat_input_w=None,
                heat_out_w=None,
                maximum_pressure_pa=None,
                maximum_temperature_k=None,
                maximum_absolute_mass_flow_kg_s=None,
                validity_verdict=None,
                failed_validity_criteria=[],
                maximum_tube_reynolds=None,
                maximum_tube_mach=None,
            )

            delta_fields(row, dd)
            results[code] = row

            print(
                f"{delta_t_k:4.0f}  {code:2}    "
                f"{'False':>5}  "
                f"{'False':>8}  "
                f"{'nan':>9}  "
                f"{'nan':>9}  "
                f"{'nan':>8}  "
                f"{'nan':>9}  "
                f"{'-':>6}  "
                f"{elapsed:8.3f}  "
                f"INVALID: {error}",
                flush=True,
            )

            continue

            elapsed = time.perf_counter() - started

            if not captures:
                raise RuntimeError(
                    f"{delta_t_k:g} K {code}: periodic observer "
                    "did not receive a result."
                )

            periodic = captures[0]
            stats = periodic.backend_statistics or {}

            actual_backend = stats.get("actual_backend")
            fallback_calls = int(stats.get("fallback_calls", 0) or 0)

            if actual_backend != "numba" or fallback_calls:
                raise RuntimeError(
                    f"{delta_t_k:g} K {code} did not run entirely "
                    f"on Numba: {stats}"
                )

            converged = result.get("status") == "converged"

            if converged:
                feasible, reasons = feasibility(
                    result,
                    power_floor_w,
                )
            else:
                feasible, reasons = False, []

            normalized = normalize_result_fields(result)
            eta = normalized["indicated_thermal_efficiency"]
            carnot = carnot_efficiency(delta_t_k)

            row = dict(
                topology_code=code,
                heat_in_valve_placement=heat_in,
                heat_out_valve_placement=heat_out,
                evaluation_origin="valve_placement_replay",
                converged=converged,
                feasible_under_original_constraints=bool(feasible),
                failed_constraints=list(reasons),
                fraction_of_carnot=(
                    eta / carnot
                    if eta is not None and carnot > 0.0
                    else None
                ),
                periodic_cycle_count=len(periodic.history),
                elapsed_seconds=elapsed,
                requested_rhs_backend=stats.get("requested_backend"),
                actual_rhs_backend=actual_backend,
                fallback_count=fallback_calls,
                fallback_reason=stats.get("fallback_reason"),
                final_periodic_state=(
                    state.tolist() if state is not None else None
                ),
                backend_statistics=stats,
                **normalized,
            )

            delta_fields(row, dd)
            results[code] = row

            eta_display = (
                100.0 * eta if eta is not None else float("nan")
            )
            gain_display = (
                row["efficiency_gain_percentage_points"]
                if row["efficiency_gain_percentage_points"] is not None
                else float("nan")
            )
            power_display = (
                row["indicated_power_w"]
                if row["indicated_power_w"] is not None
                else float("nan")
            )
            power_gain_display = (
                row["power_absolute_gain_w"]
                if row["power_absolute_gain_w"] is not None
                else float("nan")
            )

            print(
                f"{delta_t_k:4.0f}  {code:2}    "
                f"{str(converged):>5}  "
                f"{str(feasible):>8}  "
                f"{eta_display:9.4f}  "
                f"{gain_display:+9.4f}  "
                f"{power_display:8.3f}  "
                f"{power_gain_display:+9.3f}  "
                f"{len(periodic.history):6d}  "
                f"{elapsed:8.3f}",
                flush=True,
            )

        best_converged = best_topology(
            results,
            feasible_only=False,
        )
        best_feasible = best_topology(
            results,
            feasible_only=True,
        )

        entry = dict(
            delta_t_k=delta_t_k,
            hot_source_temperature_k=hot_k,
            power_floor_w=power_floor_w,
            source_campaign_complete=source_complete,
            source_report=str(source_path.relative_to(ROOT)),
            source_sha256=source_sha,
            results=results,
            best_converged_efficiency=best_converged,
            best_feasible_efficiency=best_feasible,
        )

        payload["temperatures"].append(entry)

        # Checkpoint after every completed temperature.
        payload["elapsed_seconds"] = (
            time.perf_counter() - total_started
        )
        write_report(payload)

        provisional = " [PROVISIONAL SOURCE]" if not source_complete else ""
        print(
            f"      -> best converged={best_converged}, "
            f"best feasible={best_feasible}{provisional}",
            flush=True,
        )
        print(flush=True)

    payload["complete"] = True
    payload["elapsed_seconds"] = time.perf_counter() - total_started

    feasible_winners = [
        entry["best_feasible_efficiency"]
        for entry in payload["temperatures"]
        if entry["best_feasible_efficiency"] is not None
    ]

    if feasible_winners and len(set(feasible_winners)) == 1:
        payload["single_feasible_winner_across_map"] = feasible_winners[0]
    else:
        payload["single_feasible_winner_across_map"] = None

    write_report(payload)

    print("=" * 88)
    print("SUMMARY")
    print("=" * 88)

    for entry in payload["temperatures"]:
        marker = (
            ""
            if entry["source_campaign_complete"]
            else " (source campaign incomplete)"
        )
        print(
            f"ΔT={entry['delta_t_k']:g} K -> "
            f"best feasible: {entry['best_feasible_efficiency']}"
            f"{marker}"
        )

    if payload["single_feasible_winner_across_map"] is not None:
        print(
            "\nSame feasible efficiency winner at every screened point:",
            payload["single_feasible_winner_across_map"],
        )
    else:
        print(
            "\nThe best feasible valve placement depends on temperature "
            "or some points have no common winner."
        )

    print(
        f"\nTotal elapsed: {payload['elapsed_seconds']:.3f} s"
    )
    print(f"Report: {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
