#!/usr/bin/env python3
"""Extract candidates near a common COP-Qcold point and compare only optimized variables."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from dada_solver.research.report import inspect

DEFAULT_ROOT = Path("outputs/pedal_cell_retro/compact/frequency")


def parse_frequency(name: str) -> float | None:
    m = re.fullmatch(r"(\d+)p(\d+)Hz", name)
    if not m:
        return None
    return float(f"{m.group(1)}.{m.group(2)}")


def campaign_dirs(freq_dir: Path) -> list[Path]:
    return [
        p for p in sorted(freq_dir.iterdir())
        if p.is_dir() and (p / "definition.json").exists() and (p / "study.json").exists()
    ]


def flatten_dict(d: dict, prefix: str = "") -> dict:
    out = {}
    for k, v in d.items():
        key = f"{prefix}.{k}" if prefix else str(k)
        if isinstance(v, dict):
            out.update(flatten_dict(v, key))
        elif isinstance(v, (str, int, float, bool)) or v is None:
            out[key] = v
    return out


def distance(q, cop, tq, tc, qs, cs):
    return math.hypot((q - tq) / qs, (cop - tc) / cs)


def rows_from_campaign(campaign, freq, args):
    data = inspect(campaign)

    # Research stores fixed values separately from the coordinates that were
    # actually optimized.  Keep the latter so the comparison CSV contains
    # only variables that were free during the search.
    scientific = data.get("scientific") or {}
    active_parameter_names = {
        p["name"]
        for p in scientific.get("parameters", [])
        if isinstance(p, dict) and p.get("name")
    }

    unique = {r["candidate_id"]: r for r in data["records"]}
    rows = []
    for r in unique.values():
        if r.get("status") != "feasible":
            continue
        m = r.get("metrics") or {}
        q = m.get("cooling_power_w")
        cop = m.get("cooling_cop")
        if q is None or cop is None:
            continue
        q, cop = float(q), float(cop)
        if not (math.isfinite(q) and math.isfinite(cop)):
            continue
        if args.qcold_tol is not None and abs(q - args.qcold) > args.qcold_tol:
            continue
        if args.cop_tol is not None and abs(cop - args.cop) > args.cop_tol:
            continue
        rows.append({
            "frequency_hz": freq,
            "campaign": str(campaign),
            "candidate_id": r["candidate_id"],
            "cooling_power_w": q,
            "cooling_cop": cop,
            "delta_qcold_w": q - args.qcold,
            "delta_cop": cop - args.cop,
            "distance": distance(q, cop, args.qcold, args.cop, args.qcold_scale, args.cop_scale),
            "physical": r.get("physical") or {},
            "active_parameter_names": sorted(active_parameter_names),
            "metrics": m,
        })
    return rows


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract candidates nearest a common COP-Qcold point across frequencies."
    )
    p.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    p.add_argument("--qcold", type=float, default=13.0)
    p.add_argument("--cop", type=float, default=3.3)
    p.add_argument("--qcold-scale", type=float, default=1.0,
                   help="Distance normalization for Qcold (default 1 W).")
    p.add_argument("--cop-scale", type=float, default=0.1,
                   help="Distance normalization for COP (default 0.1).")
    p.add_argument("--qcold-tol", type=float, default=None)
    p.add_argument("--cop-tol", type=float, default=None)
    p.add_argument("--output-prefix", type=Path,
                   default=DEFAULT_ROOT / "common_point_candidates")
    return p.parse_args()


def main():
    args = parse_args()
    if args.qcold_scale <= 0 or args.cop_scale <= 0:
        raise SystemExit("Scales must be positive.")

    root = ROOT / args.root
    all_selected = []

    freq_dirs = []
    for d in root.iterdir():
        if d.is_dir():
            f = parse_frequency(d.name)
            if f is not None:
                freq_dirs.append((f, d))

    for freq, freq_dir in sorted(freq_dirs):
        rows = []
        for campaign in campaign_dirs(freq_dir):
            try:
                rows.extend(rows_from_campaign(campaign, freq, args))
            except Exception as exc:
                print(f"WARNING {campaign}: {exc}", file=sys.stderr)

        if not rows:
            continue

        chosen = min(
            rows,
            key=lambda r: (
                r["distance"],
                abs(r["delta_qcold_w"]),
                abs(r["delta_cop"]),
            ),
        )
        all_selected.append(chosen)

    if not all_selected:
        raise SystemExit("No matching feasible candidates found.")

    prefix = ROOT / args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)

    summary_csv = prefix.with_suffix(".csv")
    params_csv = prefix.with_name(prefix.name + "_parameters").with_suffix(".csv")
    json_path = prefix.with_suffix(".json")

    with summary_csv.open("w", newline="", encoding="utf-8") as f:
        fields = [
            "frequency_hz", "campaign", "candidate_id",
            "cooling_power_w", "cooling_cop",
            "delta_qcold_w", "delta_cop", "distance",
        ]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in all_selected:
            w.writerow({k: r[k] for k in fields})

    # Detailed comparison: only coordinates that were actually optimized.
    # No fixed inputs and no derived diagnostics are included here.
    flat_rows = []
    active_union = set()

    for r in all_selected:
        physical = r["physical"]
        active_names = set(r.get("active_parameter_names") or [])
        active_union.update(active_names)

        row = {
            "frequency_hz": r["frequency_hz"],
            "campaign": r["campaign"],
            "candidate_id": r["candidate_id"],
            "cooling_power_w": r["cooling_power_w"],
            "cooling_cop": r["cooling_cop"],
        }

        for name in active_names:
            if name in physical:
                row[name] = physical[name]

        flat_rows.append(row)

    fixed = [
        "frequency_hz", "campaign", "candidate_id",
        "cooling_power_w", "cooling_cop",
    ]
    fields = fixed + sorted(active_union)

    with params_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(flat_rows)

    json_path.write_text(
        json.dumps({
            "target": {
                "cooling_power_w": args.qcold,
                "cooling_cop": args.cop,
                "qcold_scale": args.qcold_scale,
                "cop_scale": args.cop_scale,
                "qcold_tolerance": args.qcold_tol,
                "cop_tolerance": args.cop_tol,
            },
            "selected": all_selected,
        }, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    print(f"Target: Qcold={args.qcold:g} W, COP={args.cop:g}")
    for r in all_selected:
        print(
            f"{r['frequency_hz']:>5g} Hz  "
            f"Qcold={r['cooling_power_w']:.4f} W  "
            f"COP={r['cooling_cop']:.4f}  "
            f"d={r['distance']:.3f}  "
            f"{r['candidate_id'][:12]}"
        )
    print("\nWrote:")
    print(" ", summary_csv.relative_to(ROOT))
    print(" ", params_csv.relative_to(ROOT))
    print(" ", json_path.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
