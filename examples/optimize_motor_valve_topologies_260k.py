#!/usr/bin/env python3
"""Re-optimize the three alternative passive-valve topologies at ΔT = 260 K.

DD is not recomputed: the refined 260 K DD champion is used as the reference.

Each requested alternative topology gets its own independent two-pass campaign:
  1. thermo 5D: swept_ratio, n_i, length_i_m, n_o, length_o_m
  2. motion 7D: t1, t2, t3, a_l, b_l, a_s, b_s

The existing optimize_motor_temperature_point.py implementation is reused
unchanged. This launcher injects the selected valve placement into every
candidate before evaluation.

Default run:
    PYTHONPATH=src python3 examples/optimize_motor_valve_topologies_260k.py

Outputs:
    outputs/motor_valve_optimization_260K/UD/
    outputs/motor_valve_optimization_260K/DU/
    outputs/motor_valve_optimization_260K/UU/
    outputs/motor_valve_optimization_260K/summary.json

The script is resumable: rerunning it resumes each topology from its existing
history files.
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

SOURCE = ROOT / "outputs/motor_temperature_peak_refinement/dT_260K/report.json"
OUTPUT_ROOT = ROOT / "outputs/motor_valve_optimization_260K"
SUMMARY = OUTPUT_ROOT / "summary.json"

DELTA_T_K = 260.0
POWER_FLOOR_W = 25.0

TOPOLOGIES = {
    "DD": ("downstream", "downstream"),
    "UD": ("upstream", "downstream"),
    "DU": ("downstream", "upstream"),
    "UU": ("upstream", "upstream"),
}

DEFAULT_ALTERNATIVES = ("UD", "DU", "UU")


def parse_topologies(text: str) -> tuple[str, ...]:
    values = tuple(x.strip().upper() for x in text.split(",") if x.strip())
    if not values:
        raise argparse.ArgumentTypeError("At least one topology is required.")
    unknown = [x for x in values if x not in DEFAULT_ALTERNATIVES]
    if unknown:
        raise argparse.ArgumentTypeError(
            "Only the alternative topologies UD, DU and UU are optimized here; "
            f"unknown/unsupported selection: {unknown}"
        )
    return values


def champion_row(code: str, report: dict, origin: str) -> dict:
    champion = report.get("champion")
    if champion is None:
        return {
            "topology_code": code,
            "campaign_complete": bool(report.get("campaign_complete", False)),
            "origin": origin,
            "status": "no_champion",
        }

    result = champion.get("result", {})
    return {
        "topology_code": code,
        "campaign_complete": bool(report.get("campaign_complete", False)),
        "origin": origin,
        "phase": champion.get("phase"),
        "indicated_thermal_efficiency": result.get(
            "indicated_thermal_efficiency"
        ),
        "fraction_of_carnot": champion.get("fraction_of_carnot"),
        "indicated_power_w": result.get("indicated_power_w"),
        "parameters": champion.get("parameters"),
        "physical_constraint_failures": champion.get(
            "physical_constraint_failures", []
        ),
        "feasible": bool(champion.get("feasible", False)),
        "last_complete_state": champion.get("last_complete_state"),
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


def save_summary(rows: list[dict], elapsed: float | None = None) -> dict:
    valid = [
        row
        for row in rows
        if row.get("campaign_complete")
        and row.get("indicated_thermal_efficiency") is not None
        and row.get("feasible", True)
    ]

    best_eta = (
        max(valid, key=lambda row: row["indicated_thermal_efficiency"])
        if valid
        else None
    )
    best_carnot = (
        max(valid, key=lambda row: row["fraction_of_carnot"])
        if valid
        else None
    )

    dd = next((row for row in rows if row["topology_code"] == "DD"), None)
    dd_eta = (
        dd.get("indicated_thermal_efficiency")
        if dd is not None
        else None
    )
    dd_power = dd.get("indicated_power_w") if dd is not None else None
    dd_ratio = dd.get("fraction_of_carnot") if dd is not None else None

    enriched = []
    for original in rows:
        row = dict(original)
        eta = row.get("indicated_thermal_efficiency")
        power = row.get("indicated_power_w")
        ratio = row.get("fraction_of_carnot")

        row["efficiency_gain_vs_DD"] = (
            eta - dd_eta
            if eta is not None and dd_eta is not None
            else None
        )
        row["efficiency_gain_vs_DD_percentage_points"] = (
            100.0 * (eta - dd_eta)
            if eta is not None and dd_eta is not None
            else None
        )
        row["power_gain_vs_DD_w"] = (
            power - dd_power
            if power is not None and dd_power is not None
            else None
        )
        row["fraction_of_carnot_gain_vs_DD"] = (
            ratio - dd_ratio
            if ratio is not None and dd_ratio is not None
            else None
        )
        row["fraction_of_carnot_gain_vs_DD_percentage_points"] = (
            100.0 * (ratio - dd_ratio)
            if ratio is not None and dd_ratio is not None
            else None
        )
        enriched.append(row)

    payload = {
        "study": "260 K independent passive-valve topology re-optimization",
        "delta_t_k": DELTA_T_K,
        "power_floor_w": POWER_FLOOR_W,
        "source_DD_report": str(SOURCE.relative_to(ROOT)),
        "rows": enriched,
        "best_feasible_indicated_efficiency": (
            best_eta["topology_code"] if best_eta else None
        ),
        "best_feasible_fraction_of_carnot": (
            best_carnot["topology_code"] if best_carnot else None
        ),
        "elapsed_seconds": elapsed,
    }

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    SUMMARY.write_text(json.dumps(payload, indent=2) + "\n")
    return payload


def print_summary(payload: dict) -> None:
    print()
    print(
        "topo | complete | feasible | eta % | eta/Carnot % | "
        "P W | Δη pp | Δ(eta/C) pp"
    )
    print("-" * 86)

    for row in payload["rows"]:
        eta = row.get("indicated_thermal_efficiency")
        ratio = row.get("fraction_of_carnot")
        power = row.get("indicated_power_w")
        de = row.get("efficiency_gain_vs_DD_percentage_points")
        dc = row.get("fraction_of_carnot_gain_vs_DD_percentage_points")

        def f(value, width, digits):
            if value is None:
                return f"{'-':>{width}}"
            return f"{value:{width}.{digits}f}"

        print(
            f"{row['topology_code']:>4} | "
            f"{str(row.get('campaign_complete', False)):>8} | "
            f"{str(row.get('feasible', False)):>8} | "
            f"{f(100.0 * eta if eta is not None else None, 7, 4)} | "
            f"{f(100.0 * ratio if ratio is not None else None, 12, 4)} | "
            f"{f(power, 7, 3)} | "
            f"{f(de, 6, 4)} | "
            f"{f(dc, 11, 4)}"
        )

    print()
    print(
        "Best feasible indicated efficiency:",
        payload["best_feasible_indicated_efficiency"],
    )
    print(
        "Best feasible eta/Carnot:",
        payload["best_feasible_fraction_of_carnot"],
    )
    if payload.get("elapsed_seconds") is not None:
        print(f"Total elapsed: {payload['elapsed_seconds']:.3f} s")
    print("Summary:", SUMMARY.relative_to(ROOT))


def write_topology_manifest(
    output_dir: Path,
    code: str,
    heat_in: str,
    heat_out: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "valve_topology.json"
    expected = {
        "topology_code": code,
        "heat_in_valve_placement": heat_in,
        "heat_out_valve_placement": heat_out,
        "delta_t_k": DELTA_T_K,
        "source_report": str(SOURCE.relative_to(ROOT)),
    }

    if path.exists():
        existing = json.loads(path.read_text())
        if existing != expected:
            raise RuntimeError(
                f"Topology manifest changed in {output_dir}; "
                "use a fresh output directory."
            )
    else:
        path.write_text(json.dumps(expected, indent=2) + "\n")


def worker(args: argparse.Namespace) -> None:
    """Run one topology in an isolated Python process."""
    code = args.worker
    heat_in, heat_out = TOPOLOGIES[code]
    output_dir = OUTPUT_ROOT / code

    write_topology_manifest(output_dir, code, heat_in, heat_out)

    # Import only inside the worker: the optimizer module is patched locally
    # and the parent process stays untouched.
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

    # The current optimizer already selects WallBackendSettings("numba").
    optimizer_argv = [
        str(Path(optimizer.__file__)),
        "--delta-t-k",
        f"{DELTA_T_K:g}",
        "--power-floor-w",
        f"{POWER_FLOOR_W:g}",
        "--source-report",
        str(SOURCE),
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
        f"=== {code}: H_i={heat_in}, H_o={heat_out} ===",
        flush=True,
    )
    old_argv = sys.argv
    try:
        sys.argv = optimizer_argv
        optimizer.main()
    finally:
        sys.argv = old_argv


def run_parent(args: argparse.Namespace) -> None:
    if not SOURCE.exists():
        raise FileNotFoundError(SOURCE)

    source_report = json.loads(SOURCE.read_text())
    if not source_report.get("campaign_complete"):
        raise RuntimeError("The refined 260 K DD source campaign is incomplete.")
    if source_report.get("champion") is None:
        raise RuntimeError("The refined 260 K DD source has no champion.")

    started = time.perf_counter()

    # DD is the already-refined reference and is never recomputed.
    rows = [champion_row("DD", source_report, "refined_DD_reference")]
    save_summary(rows, 0.0)

    for code in args.topologies:
        cmd = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--worker",
            code,
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

        completed = subprocess.run(cmd, cwd=ROOT, env=env)

        report_path = OUTPUT_ROOT / code / "report.json"
        if report_path.exists():
            report = json.loads(report_path.read_text())
            rows.append(
                champion_row(code, report, "independent_reoptimization")
            )
        else:
            heat_in, heat_out = TOPOLOGIES[code]
            rows.append(
                {
                    "topology_code": code,
                    "campaign_complete": False,
                    "origin": "independent_reoptimization",
                    "status": "worker_failed_without_report",
                    "returncode": completed.returncode,
                    "heat_in_valve_placement": heat_in,
                    "heat_out_valve_placement": heat_out,
                }
            )

        payload = save_summary(
            rows,
            time.perf_counter() - started,
        )
        print_summary(payload)

        if completed.returncode != 0:
            print(
                f"\nWARNING: {code} worker returned "
                f"{completed.returncode}; continuing with the next topology.",
                flush=True,
            )

    payload = save_summary(
        rows,
        time.perf_counter() - started,
    )
    print_summary(payload)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--topologies",
        type=parse_topologies,
        default=DEFAULT_ALTERNATIVES,
        help="Comma-separated alternative topologies; default: UD,DU,UU.",
    )
    parser.add_argument(
        "--worker",
        choices=DEFAULT_ALTERNATIVES,
        help=argparse.SUPPRESS,
    )

    parser.add_argument("--thermo-evaluations", type=int, default=96)
    parser.add_argument("--motion-evaluations", type=int, default=96)
    parser.add_argument(
        "--thermo-evaluations-per-radius",
        type=int,
        default=24,
    )
    parser.add_argument(
        "--motion-evaluations-per-radius",
        type=int,
        default=24,
    )
    parser.add_argument(
        "--thermo-budget-seconds",
        type=float,
        default=3600.0,
    )
    parser.add_argument(
        "--motion-budget-seconds",
        type=float,
        default=3600.0,
    )
    parser.add_argument(
        "--candidate-seconds",
        type=float,
        default=900.0,
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

    if args.worker:
        worker(args)
    else:
        run_parent(args)


if __name__ == "__main__":
    main()
