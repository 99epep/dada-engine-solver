#!/usr/bin/env python3
"""Raise DADA refrigeration frequency step by step while reoptimizing only the 6 microtube geometry coordinates."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tomllib


ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))


from dada_solver.research.study_io import dumps
from dada_solver.research.report import inspect, select_records

HX_ACTIVE = {
    "microtube.heat_in.tube_count",
    "microtube.heat_in.tube_length_m",
    "microtube.heat_in.inner_diameter_m",
    "microtube.heat_out.tube_count",
    "microtube.heat_out.tube_length_m",
    "microtube.heat_out.inner_diameter_m",
}

DEFAULT_SOURCE = Path("outputs/pedal_cell_retro/opt8h_local/campaign")
DEFAULT_CANDIDATE = "bf589558b59a4f2b3cf0a4c208ccf03ffb5340bfe9cd437635f90c0eb23ac689"
DEFAULT_OUTPUT = Path("outputs/pedal_cell_retro/frequency_hx_continuation")
DEFAULT_FREQUENCIES = [0.25, 0.40, 0.60, 0.80, 1.00, 1.40, 2.00, 2.80]

def run_research(args: list[str]) -> None:
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["PYTHONPATH"] = str(SRC)

    cmd = [sys.executable, "-m", "dada_solver.research", *args]
    print("\n+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, env=env)
    
def run(cmd: list[str]) -> None:
    print("\n+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def make_fixed(row: dict, value) -> dict:
    return {"name": row["name"], "unit": row["unit"], "value": value}


def patch_study(study_path: Path, frequency_hz: float, step_index: int) -> None:
    """Freeze all active coordinates except the 6 HX coordinates; fix frequency to this step."""
    raw = tomllib.loads(study_path.read_text(encoding="utf-8"))
    params = []
    seen_frequency = False

    for row in raw["parameters"]:
        name = row["name"]
        if name == "operation.frequency_hz":
            params.append(make_fixed(row, float(frequency_hz)))
            seen_frequency = True
        elif "initial" in row and name not in HX_ACTIVE:
            # refine() has already replaced initial by the selected source candidate value.
            params.append(make_fixed(row, row["initial"]))
        else:
            params.append(row)

    if not seen_frequency:
        raise RuntimeError("operation.frequency_hz not found in study")

    raw["parameters"] = params
    active = {p["name"] for p in params if "initial" in p}
    if active != HX_ACTIVE:
        raise RuntimeError(
            "Expected exactly the six HX coordinates active; "
            f"missing={sorted(HX_ACTIVE-active)}, extra={sorted(active-HX_ACTIVE)}"
        )

    search = raw["search"]
    if search.get("domain") != "local_regions_v1":
        raise RuntimeError(f"Expected local_regions_v1, got {search.get('domain')!r}")

    # local_regions_v1 requires center keys to match the active parameter space exactly.
    for region in search["regions"]:
        center = region["center"]
        region["center"] = {name: center[name] for name in sorted(HX_ACTIVE)}

    raw["study"]["name"] = f"Pedal Cell Retro — HX continuation {frequency_hz:g} Hz"
    raw["study"]["purpose"] = (
        "progressive_frequency_continuation_with_fixed_kinematics_"
        "and_local_microtube_reoptimization"
    )
    study_path.write_text(dumps(raw), encoding="utf-8")


def best_record(campaign: Path) -> dict:
    data = inspect(campaign)
    chosen = select_records(data, ["best"])
    if not chosen:
        raise RuntimeError(f"No feasible candidate found in {campaign}")
    record = chosen[0]
    if record["status"] != "feasible":
        raise RuntimeError(f"Best candidate is not feasible: {record['status']}")
    return record


def compact_result(freq: float, record: dict, campaign: Path) -> dict:
    metrics = record.get("metrics") or {}
    derived = record.get("derived") or {}
    return {
        "frequency_hz": freq,
        "candidate_id": record["candidate_id"],
        "campaign": str(campaign),
        "status": record["status"],
        "cooling_cop": metrics.get("cooling_cop"),
        "cooling_power_w": metrics.get("cooling_power_w"),
        "indicated_mechanical_input_power_w": metrics.get("indicated_mechanical_input_power_w"),
        "maximum_tube_mach_number": derived.get("maximum_tube_mach_number"),
        "maximum_mach_number": (metrics.get("validity") or {}).get("maximum_mach_number"),
        "maximum_absolute_mass_flow_kg_s": metrics.get("maximum_absolute_mass_flow_kg_s"),
        "maximum_pressure_pa": metrics.get("maximum_pressure_pa"),
        "maximum_temperature_k": metrics.get("maximum_temperature_k"),
        "physical": record.get("physical"),
    }


def write_progress(root: Path) -> None:
    rows = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(root.glob("step_*/step_result.json"))]
    (root / "progress.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")

    fields = [
        "frequency_hz", "candidate_id", "cooling_cop", "cooling_power_w",
        "indicated_mechanical_input_power_w", "maximum_tube_mach_number",
        "maximum_mach_number", "maximum_absolute_mass_flow_kg_s",
        "maximum_pressure_pa", "maximum_temperature_k", "campaign",
    ]
    with (root / "progress.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k) for k in fields})


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Progressively increase frequency while locally reoptimizing only the two microtube exchangers."
    )
    p.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    p.add_argument("--candidate", default=DEFAULT_CANDIDATE)
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--frequencies", nargs="+", type=float, default=DEFAULT_FREQUENCIES, metavar="HZ")
    p.add_argument("--radius", type=float, default=0.15)
    p.add_argument("--budget", default="45m", help="Research wall-clock budget per step")
    p.add_argument("--max-candidates", type=int, default=512)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not args.source.exists():
        raise SystemExit(f"Source does not exist: {args.source}")
    if not (0 < args.radius <= 1):
        raise SystemExit("--radius must be in (0,1]")
    if args.max_candidates <= 0:
        raise SystemExit("--max-candidates must be positive")
    if not args.frequencies or any(f <= 0 for f in args.frequencies):
        raise SystemExit("All frequencies must be positive")
    if any(b <= a for a, b in zip(args.frequencies, args.frequencies[1:])):
        raise SystemExit("Frequencies must be strictly increasing")

    root = args.output
    root.mkdir(parents=True, exist_ok=True)
    source = args.source
    candidate = args.candidate

    print("Progressive frequency / HX continuation")
    print("  source:     ", source)
    print("  candidate:  ", candidate)
    print("  frequencies:", " ".join(f"{x:g}" for x in args.frequencies), "Hz")
    print("  radius:     ", args.radius)
    print("  budget:     ", args.budget, "per step")
    print("  output:     ", root)

    for index, freq in enumerate(args.frequencies, 1):
        label = f"{freq:g}".replace(".", "p")
        step = root / f"step_{index:02d}_{label}Hz"
        study = step / "study.toml"
        campaign = step / "campaign"
        marker = step / "step_result.json"
        step.mkdir(parents=True, exist_ok=True)

        if marker.exists():
            result = json.loads(marker.read_text(encoding="utf-8"))
            candidate = result["candidate_id"]
            source = Path(result["campaign"])
            print(f"\n[{index:02d}] {freq:g} Hz already complete: COP={result.get('cooling_cop')} best={candidate[:12]}")
            continue

        print(f"\n=== STEP {index:02d}: {freq:g} Hz — source {candidate[:12]} ===")

        if not study.exists():
            run_research([
                "refine", str(source), "--candidate", candidate,
                "--radius", str(args.radius), "--output", str(study),
            ])
        else:
            print("Using existing study:", study)

        # Idempotent: also repairs a study left between refine() and patching
        # if the script was interrupted at exactly that point.
        patch_study(study, freq, index)
        run_research(["validate", str(study)])

        if campaign.exists() and any(campaign.iterdir()):
            print("Resuming existing campaign:", campaign)
            run_research([
                "resume", str(campaign), "--budget", args.budget,
                "--max-candidates", str(args.max_candidates),
            ])
        else:
            run_research([
                "run", str(study), "--budget", args.budget,
                "--max-candidates", str(args.max_candidates), "--directory", str(campaign),
            ])

        record = best_record(campaign)
        result = compact_result(freq, record, campaign)
        marker.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        write_progress(root)

        print(f"\nSTEP {index:02d} COMPLETE")
        print(f"  frequency: {freq:g} Hz")
        print(f"  best:      {record['candidate_id']}")
        print(f"  COP:       {result['cooling_cop']}")
        print(f"  Qcold:     {result['cooling_power_w']} W")
        print(f"  Pin:       {result['indicated_mechanical_input_power_w']} W")
        print(f"  tube Mach: {result['maximum_tube_mach_number']}")

        # Next step starts from this optimized exchanger geometry.
        source = campaign
        candidate = record["candidate_id"]

    write_progress(root)
    print("\nContinuation complete")
    print("  JSON:", root / "progress.json")
    print("  CSV: ", root / "progress.csv")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as exc:
        print(f"\nCommand failed with exit code {exc.returncode}.", file=sys.stderr)
        print("Campaign state is preserved; rerun the script to resume the current step.", file=sys.stderr)
        raise SystemExit(exc.returncode)
