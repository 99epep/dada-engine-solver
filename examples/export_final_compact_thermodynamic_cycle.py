#!/usr/bin/env python3
"""Export the final compact DADA-engine champion for SVG animation.

Mechanical source:
  outputs/sixbar_rank01_compact_rods_hlat90/best_pair.json

Thermodynamic source:
  outputs/sixbar_rank01_compact_thermo5d_deep/rank_01/report.json

This is intentionally a thin wrapper around the already validated
examples/export_rank01_thermodynamic_cycle.py exporter. It only redirects
that exporter to the final compact mechanism and final deep thermo-5D
champion, preserving its CSV schema, periodic-cycle reconstruction, angular
convention and animation-model format.

Outputs:
  outputs/final_compact_thermodynamic_cycle.csv
  outputs/final_compact_thermodynamic_cycle_raw.csv
  outputs/final_compact_thermodynamic_cycle_metadata.json
  outputs/final_compact_animation_model.json

Run from the repository root:

  PYTHONPATH=src python3 examples/export_final_compact_thermodynamic_cycle.py

Then animate with:

  python3 examples/animate_sixbar_graphs_svg.py \\
    --source outputs/final_compact_thermodynamic_cycle.csv \\
    --model outputs/final_compact_animation_model.json
"""

from __future__ import annotations

from pathlib import Path
import json

import export_rank01_thermodynamic_cycle as exporter


ROOT = Path.cwd()

THERMO_REPORT = (
    ROOT
    / "outputs"
    / "sixbar_rank01_compact_thermo5d_deep"
    / "rank_01"
    / "report.json"
)

PAIR_PATH = (
    ROOT
    / "outputs"
    / "sixbar_rank01_compact_rods_hlat90"
    / "best_pair.json"
)

PAIR_REPORT = (
    ROOT
    / "outputs"
    / "sixbar_rank01_compact_rods_hlat90"
    / "report.json"
)

OUTPUT_PREFIX = ROOT / "outputs" / "final_compact_thermodynamic_cycle"
MODEL_OUTPUT = ROOT / "outputs" / "final_compact_animation_model.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def validate_sources() -> None:
    for path in (THERMO_REPORT, PAIR_PATH, PAIR_REPORT):
        if not path.exists():
            raise FileNotFoundError(path)

    thermo = load_json(THERMO_REPORT)
    pair = load_json(PAIR_PATH)
    pair_report = load_json(PAIR_REPORT)

    if int(thermo.get("family_rank", -1)) != 1:
        raise RuntimeError(
            f"Expected family_rank=1 in {THERMO_REPORT}, "
            f"got {thermo.get('family_rank')!r}."
        )

    thermo_best = thermo.get("best_feasible")
    if not isinstance(thermo_best, dict):
        raise RuntimeError(f"No best_feasible in {THERMO_REPORT}.")
    if thermo_best.get("result", {}).get("status") != "converged":
        raise RuntimeError("Final thermo champion is not converged.")
    if thermo_best.get("last_complete_state") is None:
        raise RuntimeError("Final thermo champion has no periodic state.")

    pair_id = pair.get("candidate_id")
    report_best = pair_report.get("best_feasible")
    if not isinstance(report_best, dict):
        raise RuntimeError(f"No best_feasible in {PAIR_REPORT}.")
    if report_best.get("candidate_id") != pair_id:
        raise RuntimeError(
            "Compact mechanism best_pair.json and report.json disagree:\n"
            f"  pair:   {pair_id}\n"
            f"  report: {report_best.get('candidate_id')}"
        )

    expected_pair_id = thermo.get("fixed_pair", {}).get("candidate_id")
    thermo_pair_id = thermo_best.get("pair_candidate_id")

    if expected_pair_id != pair_id:
        raise RuntimeError(
            "Thermo report fixed-pair ID does not match hlat90 champion:\n"
            f"  thermo fixed pair: {expected_pair_id}\n"
            f"  hlat90 pair:       {pair_id}"
        )

    if thermo_pair_id != pair_id:
        raise RuntimeError(
            "Thermo champion pair ID does not match hlat90 champion:\n"
            f"  thermo champion: {thermo_pair_id}\n"
            f"  hlat90 pair:     {pair_id}"
        )

    params = thermo_best.get("parameters", {})
    required = ("swept_ratio", "n_i", "length_i_m", "n_o", "length_o_m")
    missing = [name for name in required if name not in params]
    if missing:
        raise RuntimeError(
            f"Final thermo champion is missing parameters: {missing}"
        )


def main() -> None:
    validate_sources()

    # Redirect the thermo-5D helper used internally by the validated exporter.
    exporter.thermo5d.PAIR_SOURCES = {
        1: (PAIR_PATH, PAIR_REPORT),
    }

    # Redirect the exporter's defaults to the final deep thermo champion.
    exporter.DEFAULT_REPORT = THERMO_REPORT
    exporter.DEFAULT_PREFIX = OUTPUT_PREFIX
    exporter.DEFAULT_MODEL = MODEL_OUTPUT

    exporter.main()


if __name__ == "__main__":
    main()
