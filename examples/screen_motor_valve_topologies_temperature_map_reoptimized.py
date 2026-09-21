#!/usr/bin/env python3
"""Short re-optimization of passive-valve topologies across the existing DD temperature map.

Purpose
-------
Estimate whether the best valve topology changes with source delta-T.

Existing DD champions are used as references and are NOT recomputed.
For each requested temperature, the three alternative topologies UD/DU/UU are
independently re-optimized with a deliberately short two-pass campaign:

  * thermo 5D
  * motion 7D

Defaults are 24 thermo + 24 motion evaluations, with 6 evaluations per radius.
This is a topology-map screening, not a final optimization.

Default temperatures:
    160, 200, 240, 270, 300, 330, 360 K

Default topologies:
    UD, DU, UU

The script is resumable. Each temperature/topology has its own output directory.

Run:
    PYTHONPATH=src python3 examples/screen_motor_valve_topologies_temperature_map_reoptimized.py

Outputs:
    outputs/motor_valve_topology_temperature_screen/
        dT_160K/UD/
        dT_160K/DU/
        ...
        summary.json
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
import os
from pathlib import Path
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]

MAP_ROOT = ROOT / "outputs/motor_temperature_map"
OUTPUT_ROOT = ROOT / "outputs/motor_valve_topology_temperature_screen"
SUMMARY = OUTPUT_ROOT / "summary.json"

DEFAULT_TEMPERATURES = (160.0, 200.0, 240.0, 270.0, 300.0, 330.0, 360.0)

TOPOLOGIES = {
    "DD": ("downstream", "downstream"),
    "UD": ("upstream", "downstream"),
    "DU": ("downstream", "upstream"),
    "UU": ("upstream", "upstream"),
}
DEFAULT_ALTERNATIVES = ("UD", "DU", "UU")


def temp_tag(delta_t_k: float) -> str:
    return f"dT_{delta_t_k:g}K".replace(".", "p")


def source_report_path(delta_t_k: float) -> Path:
    return MAP_ROOT / temp_tag(delta_t_k) / "report.json"


def parse_temperatures(text: str) -> tuple[float, ...]:
    values = tuple(float(x.strip()) for x in text.split(",") if x.strip())
    if not values:
        raise argparse.ArgumentTypeError("At least one delta-T is required.")
    if any(x <= 0 for x in values):
        raise argparse.ArgumentTypeError("All delta-T values must be positive.")
    return values


def parse_topologies(text: str) -> tuple[str, ...]:
    values = tuple(x.strip().upper() for x in text.split(",") if x.strip())
    if not values:
        raise argparse.ArgumentTypeError("At least one topology is required.")
    unknown = [x for x in values if x not in DEFAULT_ALTERNATIVES]
    if unknown:
        raise argparse.ArgumentTypeError(
            f"Only UD, DU and UU are screened here; invalid: {unknown}"
        )
    return values


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def champion_row(
    code: str,
    delta_t_k: float,
    report: dict,
    origin: str,
    power_floor_w: float,
) -> dict:
    champion = report.get("champion")
    if champion is None:
        return {
            "delta_t_k": delta_t_k,
            "topology_code": code,
            "campaign_complete": bool(report.get("campaign_complete", False)),
            "origin": origin,
            "power_floor_w": power_floor_w,
            "status": "no_champion",
            "feasible": False,
        }

    result = champion.get("result", {})
    return {
        "delta_t_k": delta_t_k,
        "topology_code": code,
        "campaign_complete": bool(report.get("campaign_complete", False)),
        "origin": origin,
        "power_floor_w": power_floor_w,
        "phase": champion.get("phase"),
        "status": result.get("status"),
        "indicated_thermal_efficiency": result.get("indicated_thermal_efficiency"),
        "fraction_of_carnot": champion.get("fraction_of_carnot"),
        "indicated_power_w": result.get("indicated_power_w"),
        "heat_input_w": result.get("heat_input_w"),
        "heat_out_w": result.get("heat_out_w"),
        "parameters": champion.get("parameters"),
        "physical_constraint_failures": champion.get(
            "physical_constraint_failures", []
        ),
        "feasible": bool(champion.get("feasible", False)),
        "thermo_complete": bool(
            report.get("thermo_pass", {}).get("complete", False)
        ),
        "motion_complete": bool(
            report.get("motion_pass", {}).get("complete", False)
        ),
        "thermo_nonseed_evaluations": report.get(
            "thermo_pass", {}
        ).get("completed_nonseed_evaluations"),
        "motion_nonseed_evaluations": report.get(
            "motion_pass", {}
        ).get("completed_nonseed_evaluations"),
    }


def enrich_temperature_rows(rows: list[dict]) -> list[dict]:
    dd = next((row for row in rows if row["topology_code"] == "DD"), None)
    dd_eta = dd.get("indicated_thermal_efficiency") if dd else None
    dd_ratio = dd.get("fraction_of_carnot") if dd else None
    dd_power = dd.get("indicated_power_w") if dd else None

    enriched = []
    for original in rows:
        row = dict(original)

        eta = row.get("indicated_thermal_efficiency")
        ratio = row.get("fraction_of_carnot")
        power = row.get("indicated_power_w")

        row["efficiency_gain_vs_DD_percentage_points"] = (
            100.0 * (eta - dd_eta)
            if eta is not None and dd_eta is not None
            else None
        )
        row["fraction_of_carnot_gain_vs_DD_percentage_points"] = (
            100.0 * (ratio - dd_ratio)
            if ratio is not None and dd_ratio is not None
            else None
        )
        row["power_gain_vs_DD_w"] = (
            power - dd_power
            if power is not None and dd_power is not None
            else None
        )

        enriched.append(row)

    return enriched


def temperature_summary(delta_t_k: float, rows: list[dict]) -> dict:
    enriched = enrich_temperature_rows(rows)

    eligible = [
        row
        for row in enriched
        if row.get("campaign_complete")
        and row.get("feasible")
        and row.get("indicated_thermal_efficiency") is not None
    ]

    best_eta = (
        max(eligible, key=lambda row: row["indicated_thermal_efficiency"])
        if eligible
        else None
    )
    best_ratio = (
        max(eligible, key=lambda row: row["fraction_of_carnot"])
        if eligible
        else None
    )

    return {
        "delta_t_k": delta_t_k,
        "rows": enriched,
        "best_feasible_efficiency": (
            best_eta["topology_code"] if best_eta else None
        ),
        "best_feasible_fraction_of_carnot": (
            best_ratio["topology_code"] if best_ratio else None
        ),
    }


def save_summary(
    studies: list[dict],
    elapsed_seconds: float | None = None,
) -> dict:
    winners = {
        f"{study['delta_t_k']:g}": study["best_feasible_efficiency"]
        for study in studies
    }

    ratio_winners = {
        f"{study['delta_t_k']:g}": study[
            "best_feasible_fraction_of_carnot"
        ]
        for study in studies
    }

    completed_winners = [
        winner for winner in winners.values() if winner is not None
    ]

    payload = {
        "study": "short re-optimized passive-valve topology temperature screen",
        "temperatures_k": [study["delta_t_k"] for study in studies],
        "topology_definitions": {
            code: {
                "heat_in_valve_placement": placements[0],
                "heat_out_valve_placement": placements[1],
            }
            for code, placements in TOPOLOGIES.items()
        },
        "studies": studies,
        "best_efficiency_topology_by_temperature": winners,
        "best_fraction_of_carnot_topology_by_temperature": ratio_winners,
        "same_efficiency_winner_at_all_completed_temperatures": (
            completed_winners[0]
            if completed_winners
            and len(set(completed_winners)) == 1
            else None
        ),
        "elapsed_seconds": elapsed_seconds,
    }

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    SUMMARY.write_text(json.dumps(payload, indent=2) + "\n")
    return payload


def fmt(value, width: int, digits: int) -> str:
    if value is None:
        return f"{'-':>{width}}"
    return f"{value:{width}.{digits}f}"


def print_temperature(study: dict) -> None:
    print()
    print(f"ΔT = {study['delta_t_k']:g} K")
    print(
        "topo | complete | feasible | eta % | eta/Carnot % | "
        "P W | Δη pp | Δ(eta/C) pp | ΔP W"
    )
    print("-" * 103)

    for row in study["rows"]:
        eta = row.get("indicated_thermal_efficiency")
        ratio = row.get("fraction_of_carnot")
        power = row.get("indicated_power_w")

        print(
            f"{row['topology_code']:>4} | "
            f"{str(row.get('campaign_complete', False)):>8} | "
            f"{str(row.get('feasible', False)):>8} | "
            f"{fmt(100.0 * eta if eta is not None else None, 7, 4)} | "
            f"{fmt(100.0 * ratio if ratio is not None else None, 12, 4)} | "
            f"{fmt(power, 7, 3)} | "
            f"{fmt(row.get('efficiency_gain_vs_DD_percentage_points'), 6, 4)} | "
            f"{fmt(row.get('fraction_of_carnot_gain_vs_DD_percentage_points'), 11, 4)} | "
            f"{fmt(row.get('power_gain_vs_DD_w'), 6, 3)}"
        )

    print(
        "Best feasible:",
        study["best_feasible_efficiency"],
        "(efficiency),",
        study["best_feasible_fraction_of_carnot"],
        "(eta/Carnot)",
    )


def write_topology_manifest(
    output_dir: Path,
    code: str,
    heat_in: str,
    heat_out: str,
    delta_t_k: float,
    source: Path,
    power_floor_w: float,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "valve_topology.json"

    expected = {
        "topology_code": code,
        "heat_in_valve_placement": heat_in,
        "heat_out_valve_placement": heat_out,
        "delta_t_k": delta_t_k,
        "power_floor_w": power_floor_w,
        "source_report": str(source.relative_to(ROOT)),
    }

    if path.exists():
        existing = load_json(path)
        if existing != expected:
            raise RuntimeError(
                f"Topology manifest changed in {output_dir}; "
                "use a fresh output directory."
            )
    else:
        path.write_text(json.dumps(expected, indent=2) + "\n")


def worker(args: argparse.Namespace) -> None:
    delta_t_k = args.worker_delta_t
    code = args.worker_topology

    source = source_report_path(delta_t_k)
    if not source.exists():
        raise FileNotFoundError(source)

    source_data = load_json(source)
    power_floor_w = float(source_data["power_floor_w"])

    heat_in, heat_out = TOPOLOGIES[code]
    output_dir = OUTPUT_ROOT / temp_tag(delta_t_k) / code

    write_topology_manifest(
        output_dir,
        code,
        heat_in,
        heat_out,
        delta_t_k,
        source,
        power_floor_w,
    )

    examples_dir = str(ROOT / "examples")
    if examples_dir not in sys.path:
        sys.path.insert(0, examples_dir)

    import optimize_motor_temperature_point as optimizer

    original_build_design = optimizer.build_design

    def build_design_with_topology(base, geometry, parameters, hot_k):
        design = original_build_design(
            base,
            geometry,
            parameters,
            hot_k,
        )
        return replace(
            design,
            configuration=replace(
                design.configuration,
                heat_in_valve_placement=heat_in,
                heat_out_valve_placement=heat_out,
            ),
        )

    optimizer.build_design = build_design_with_topology

    optimizer_argv = [
        str(Path(optimizer.__file__)),
        "--delta-t-k",
        f"{delta_t_k:g}",
        "--power-floor-w",
        f"{power_floor_w:g}",
        "--source-report",
        str(source),
        "--output-directory",
        str(output_dir),
        "--thermo-evaluations",
        str(args.thermo_evaluations),
        "--motion-evaluations",
        str(args.motion_evaluations),
        "--thermo-evaluations-per-radius",
        str(args.thermo_evaluations_per_radius),
        "--motion-evaluations-per-radius",
        str(args.motion_evaluations_per_radius),
        "--thermo-budget-seconds",
        str(args.thermo_budget_seconds),
        "--motion-budget-seconds",
        str(args.motion_budget_seconds),
        "--candidate-seconds",
        str(args.candidate_seconds),
        "--seed",
        str(args.seed),
    ]

    print(
        f"=== ΔT={delta_t_k:g} K {code}: "
        f"H_i={heat_in}, H_o={heat_out}, "
        f"power floor={power_floor_w:g} W ===",
        flush=True,
    )

    old_argv = sys.argv
    try:
        sys.argv = optimizer_argv
        optimizer.main()
    finally:
        sys.argv = old_argv


def run_worker_process(
    args: argparse.Namespace,
    delta_t_k: float,
    code: str,
) -> int:
    cmd = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker-topology",
        code,
        "--worker-delta-t",
        f"{delta_t_k:g}",
        "--thermo-evaluations",
        str(args.thermo_evaluations),
        "--motion-evaluations",
        str(args.motion_evaluations),
        "--thermo-evaluations-per-radius",
        str(args.thermo_evaluations_per_radius),
        "--motion-evaluations-per-radius",
        str(args.motion_evaluations_per_radius),
        "--thermo-budget-seconds",
        str(args.thermo_budget_seconds),
        "--motion-budget-seconds",
        str(args.motion_budget_seconds),
        "--candidate-seconds",
        str(args.candidate_seconds),
        "--seed",
        str(args.seed),
    ]

    env = os.environ.copy()
    src = str(ROOT / "src")
    env["PYTHONPATH"] = (
        src
        + (
            os.pathsep + env["PYTHONPATH"]
            if env.get("PYTHONPATH")
            else ""
        )
    )

    return subprocess.run(cmd, cwd=ROOT, env=env).returncode


def collect_temperature(
    delta_t_k: float,
    requested_topologies: tuple[str, ...],
) -> dict:
    source = source_report_path(delta_t_k)
    source_data = load_json(source)
    power_floor_w = float(source_data["power_floor_w"])

    rows = [
        champion_row(
            "DD",
            delta_t_k,
            source_data,
            "existing_DD_reference",
            power_floor_w,
        )
    ]

    for code in requested_topologies:
        report_path = (
            OUTPUT_ROOT
            / temp_tag(delta_t_k)
            / code
            / "report.json"
        )
        if report_path.exists():
            rows.append(
                champion_row(
                    code,
                    delta_t_k,
                    load_json(report_path),
                    "short_independent_reoptimization",
                    power_floor_w,
                )
            )

    return temperature_summary(delta_t_k, rows)


def run_parent(args: argparse.Namespace) -> None:
    for delta_t_k in args.temperatures:
        source = source_report_path(delta_t_k)
        if not source.exists():
            raise FileNotFoundError(source)

        source_data = load_json(source)
        if not source_data.get("campaign_complete"):
            raise RuntimeError(
                f"DD source campaign is incomplete at ΔT={delta_t_k:g} K."
            )
        if source_data.get("champion") is None:
            raise RuntimeError(
                f"DD source has no champion at ΔT={delta_t_k:g} K."
            )

    started = time.perf_counter()
    studies: list[dict] = []

    for delta_t_k in args.temperatures:
        print()
        print("#" * 90)
        print(f"SCREEN ΔT={delta_t_k:g} K")
        print("#" * 90)

        for code in args.topologies:
            returncode = run_worker_process(args, delta_t_k, code)
            if returncode != 0:
                print(
                    f"WARNING: ΔT={delta_t_k:g} K {code} "
                    f"returned {returncode}; continuing.",
                    flush=True,
                )

        study = collect_temperature(delta_t_k, args.topologies)
        studies.append(study)
        payload = save_summary(
            studies,
            time.perf_counter() - started,
        )
        print_temperature(study)

    payload = save_summary(
        studies,
        time.perf_counter() - started,
    )

    print()
    print("=" * 90)
    print("TEMPERATURE / TOPOLOGY SUMMARY")
    print("=" * 90)
    for study in payload["studies"]:
        print(
            f"ΔT={study['delta_t_k']:g} K -> "
            f"best eta={study['best_feasible_efficiency']}, "
            f"best eta/Carnot="
            f"{study['best_feasible_fraction_of_carnot']}"
        )

    common = payload["same_efficiency_winner_at_all_completed_temperatures"]
    if common is not None:
        print(
            f"Same efficiency winner at every completed temperature: {common}"
        )
    else:
        print("Efficiency winner changes across the screened temperatures.")

    print(f"Total elapsed: {payload['elapsed_seconds']:.3f} s")
    print("Summary:", SUMMARY.relative_to(ROOT))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--temperatures",
        type=parse_temperatures,
        default=DEFAULT_TEMPERATURES,
        help="Comma-separated delta-T values.",
    )
    parser.add_argument(
        "--topologies",
        type=parse_topologies,
        default=DEFAULT_ALTERNATIVES,
        help="Comma-separated alternative topologies; default UD,DU,UU.",
    )

    parser.add_argument(
        "--worker-topology",
        choices=DEFAULT_ALTERNATIVES,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--worker-delta-t",
        type=float,
        help=argparse.SUPPRESS,
    )

    parser.add_argument("--thermo-evaluations", type=int, default=24)
    parser.add_argument("--motion-evaluations", type=int, default=24)
    parser.add_argument(
        "--thermo-evaluations-per-radius",
        type=int,
        default=6,
    )
    parser.add_argument(
        "--motion-evaluations-per-radius",
        type=int,
        default=6,
    )
    parser.add_argument(
        "--thermo-budget-seconds",
        type=float,
        default=1800.0,
    )
    parser.add_argument(
        "--motion-budget-seconds",
        type=float,
        default=1800.0,
    )
    parser.add_argument(
        "--candidate-seconds",
        type=float,
        default=600.0,
    )
    parser.add_argument("--seed", type=int, default=26091730)

    args = parser.parse_args()

    counts = (
        args.thermo_evaluations,
        args.motion_evaluations,
        args.thermo_evaluations_per_radius,
        args.motion_evaluations_per_radius,
    )
    if min(counts) <= 0:
        raise ValueError("Evaluation counts must be positive.")

    if min(
        args.thermo_budget_seconds,
        args.motion_budget_seconds,
        args.candidate_seconds,
    ) <= 0.0:
        raise ValueError("Budgets must be positive.")

    if args.worker_topology is not None:
        if args.worker_delta_t is None:
            raise ValueError("--worker-delta-t is required in worker mode.")
        worker(args)
    else:
        run_parent(args)


if __name__ == "__main__":
    main()
