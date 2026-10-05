#!/usr/bin/env python3
"""Overnight frequency sweep for Pedal Cell Retro.

Runs local Research campaigns at 0.25, 0.35, 0.45 and 0.50 Hz, seeded from
nearby completed frequency campaigns. The script is restartable and refreshes
the COP-Qcold Pareto comparison after each completed frequency.

Run from the repository root:
    python sweep_frequency_overnight.py

Useful overrides:
    python sweep_frequency_overnight.py --budget 90m
    python sweep_frequency_overnight.py --budget 3h --max-candidates 512
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tomllib


ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from dada_solver.research.report import inspect, select_records
from dada_solver.research.study_io import dumps


FREQUENCY_ROOT = Path("outputs/pedal_cell_retro/compact/frequency")

# Target -> existing campaign used as the local-search seed.
TARGETS = (
    (0.25, FREQUENCY_ROOT / "0p3Hz" / "campaign"),
    (0.35, FREQUENCY_ROOT / "0p3Hz" / "campaign"),
    (0.45, FREQUENCY_ROOT / "0p4Hz" / "campaign"),
    (0.50, FREQUENCY_ROOT / "0p4Hz" / "campaign"),
)


def research(*args: str) -> None:
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["PYTHONPATH"] = str(SRC)
    cmd = [sys.executable, "-m", "dada_solver.research", *map(str, args)]
    print("\n+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, cwd=ROOT, env=env)


def best_feasible(campaign: Path) -> dict:
    data = inspect(ROOT / campaign)
    rows = select_records(data, ["best"])
    if not rows or rows[0].get("status") != "feasible":
        raise RuntimeError(f"No feasible seed candidate in {campaign}")
    return rows[0]


def frequency_label(freq: float) -> str:
    return f"{freq:g}".replace(".", "p") + "Hz"


def patch_frequency(study_path: Path, freq: float, seed_campaign: Path) -> None:
    raw = tomllib.loads(study_path.read_text(encoding="utf-8"))

    found = False
    for row in raw.get("parameters", []):
        if row.get("name") == "operation.frequency_hz":
            # Frequency must remain fixed; remove any active-coordinate fields
            # defensively in case the source study changes later.
            unit = row.get("unit", "Hz")
            row.clear()
            row.update(name="operation.frequency_hz", unit=unit, value=float(freq))
            found = True
            break
    if not found:
        raise RuntimeError("operation.frequency_hz not found in generated study")

    study = raw.setdefault("study", {})
    study["name"] = f"Pedal Cell Retro — frequency sweep — {freq:g} Hz"
    study["purpose"] = "cop_qcold_frequency_sweep"

    study_path.write_text(dumps(raw), encoding="utf-8")


def prepare_study(freq: float, seed_campaign: Path, radius: float) -> Path:
    folder = ROOT / FREQUENCY_ROOT / frequency_label(freq)
    folder.mkdir(parents=True, exist_ok=True)
    study = folder / "study.toml"

    if not study.exists():
        seed = best_feasible(seed_campaign)
        print(
            f"\nPreparing {freq:g} Hz from {seed_campaign} "
            f"best={seed['candidate_id'][:12]} COP="
            f"{(seed.get('metrics') or {}).get('cooling_cop')}",
            flush=True,
        )
        research(
            "refine", str(seed_campaign),
            "--candidate", seed["candidate_id"],
            "--radius", str(radius),
            "--output", str(study.relative_to(ROOT)),
        )
    else:
        print(f"\nUsing existing study: {study.relative_to(ROOT)}", flush=True)

    patch_frequency(study, freq, seed_campaign)
    research("validate", str(study.relative_to(ROOT)))
    return study


def compact_best(campaign: Path, freq: float) -> dict:
    best = best_feasible(campaign)
    metrics = best.get("metrics") or {}
    return {
        "frequency_hz": freq,
        "candidate_id": best["candidate_id"],
        "cooling_cop": metrics.get("cooling_cop"),
        "cooling_power_w": metrics.get("cooling_power_w"),
        "campaign": str(campaign),
    }


def run_target(freq: float, seed_campaign: Path, args: argparse.Namespace) -> dict:
    folder = ROOT / FREQUENCY_ROOT / frequency_label(freq)
    campaign_rel = FREQUENCY_ROOT / frequency_label(freq) / "campaign"
    campaign = ROOT / campaign_rel
    marker = folder / "overnight_sweep_complete.json"

    if marker.exists() and not args.force_resume:
        result = json.loads(marker.read_text(encoding="utf-8"))
        print(
            f"\n=== {freq:g} Hz already complete: "
            f"COP={result.get('cooling_cop')} "
            f"Qcold={result.get('cooling_power_w')} W ===",
            flush=True,
        )
        return result

    study = prepare_study(freq, seed_campaign, args.radius)
    study_rel = study.relative_to(ROOT)

    print(f"\n=== RUN {freq:g} Hz ===", flush=True)
    if campaign.exists() and any(campaign.iterdir()):
        print(f"Resuming {campaign_rel}", flush=True)
        research(
            "resume", str(campaign_rel),
            "--budget", args.budget,
            "--max-candidates", str(args.max_candidates),
        )
    else:
        research(
            "run", str(study_rel),
            "--budget", args.budget,
            "--max-candidates", str(args.max_candidates),
            "--directory", str(campaign_rel),
        )

    result = compact_best(campaign_rel, freq)
    marker.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        f"DONE {freq:g} Hz — COP={result['cooling_cop']} "
        f"Qcold={result['cooling_power_w']} W "
        f"best={result['candidate_id'][:12]}",
        flush=True,
    )
    return result


def refresh_pareto() -> None:
    comparator = ROOT / "compare_frequency_pareto.py"
    if not comparator.exists():
        print("Pareto comparator not found; skipping refresh.", flush=True)
        return
    cmd = [sys.executable, str(comparator), "--cloud"]
    print("\n+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, cwd=ROOT)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run the overnight 0.25/0.35/0.45/0.50 Hz sweep.")
    p.add_argument(
        "--budget", default="2h",
        help="Research wall-clock budget per frequency (default: 2h; four steps ~= 8h total).",
    )
    p.add_argument("--max-candidates", type=int, default=512)
    p.add_argument("--radius", type=float, default=0.20)
    p.add_argument(
        "--force-resume", action="store_true",
        help="Add another budget slice even when this script already marked the frequency complete.",
    )
    p.add_argument(
        "--no-pareto", action="store_true",
        help="Do not regenerate compare_frequency_pareto.py output after each step.",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if args.max_candidates <= 0:
        raise SystemExit("--max-candidates must be positive")
    if not (0 < args.radius <= 1):
        raise SystemExit("--radius must be in (0, 1]")

    print("Overnight COP-Qcold frequency sweep")
    print("  frequencies:    0.25  0.35  0.45  0.50 Hz")
    print("  budget:        ", args.budget, "per frequency")
    print("  max candidates:", args.max_candidates, "per frequency")
    print("  local radius:  ", args.radius)
    print("  restartable:    yes")

    # Fail early if the two required seed campaigns are unavailable.
    for _, source in TARGETS:
        if not (ROOT / source).exists():
            raise SystemExit(f"Missing seed campaign: {source}")

    results = []
    for freq, source in TARGETS:
        result = run_target(freq, source, args)
        results.append(result)
        if not args.no_pareto:
            refresh_pareto()

    summary = ROOT / FREQUENCY_ROOT / "overnight_sweep_summary.json"
    summary.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")

    print("\n=== SWEEP COMPLETE ===")
    for row in results:
        print(
            f"{row['frequency_hz']:>4g} Hz  "
            f"COP={row.get('cooling_cop')}  Qcold={row.get('cooling_power_w')} W  "
            f"{row['candidate_id'][:12]}"
        )
    print("summary:", summary.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as exc:
        print(
            f"\nCommand failed with exit code {exc.returncode}. "
            "Existing campaign state is preserved; rerun the same command to resume.",
            file=sys.stderr,
        )
        raise SystemExit(exc.returncode)
