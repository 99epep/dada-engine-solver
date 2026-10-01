#!/usr/bin/env python3

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import tomllib

from dada_solver.research.cli import evaluate
from dada_solver.research.report import inspect, select_records
from dada_solver.research.study_io import dumps


CAMPAIGN = Path(
    "outputs/pedal_cell_retro/5L_2p8Hz_coupled/campaign"
)
SOURCE_STUDY = Path(
    "outputs/pedal_cell_retro/5L_2p8Hz_coupled/study.toml"
)
SOURCE_BASIS = SOURCE_STUDY.with_suffix(".basis.json")

OUTPUT = Path(
    "outputs/pedal_cell_retro/5L_2p8Hz_valve_screen"
)

PLACEMENTS = {
    "DD": ("downstream", "downstream"),
    "UD": ("upstream",   "downstream"),
    "DU": ("downstream", "upstream"),
    "UU": ("upstream",   "upstream"),
}


def select_candidates():
    data = inspect(CAMPAIGN)

    if len(data["best"]) < 2:
        raise RuntimeError("Campaign has fewer than two ranked candidates.")

    selected = [
        ("best1", data["best"][0]),
        ("best2", data["best"][1]),
        ("22107e", select_records(data, ["22107e"])[0]),
    ]

    # Avoid evaluating the same candidate twice if 22107e is already top-2.
    unique = []
    seen = set()
    for label, record in selected:
        if record["candidate_id"] in seen:
            continue
        seen.add(record["candidate_id"])
        unique.append((label, record))

    return unique


def make_study(label, record, code, heat_in, heat_out):
    target = OUTPUT / f"{label}_{record['candidate_id'][:12]}" / code
    target.mkdir(parents=True, exist_ok=True)

    study_path = target / "study.toml"
    basis_path = target / "study.basis.json"

    raw = tomllib.loads(SOURCE_STUDY.read_text())
    basis = json.loads(SOURCE_BASIS.read_text())

    # Exact physical coordinates of the selected campaign candidate.
    physical = record["physical"]

    for row in raw["parameters"]:
        if "initial" in row:
            row["initial"] = physical[row["name"]]

    # Change ONLY passive check-valve placement.
    basis["configuration"]["heat_in_valve_placement"] = heat_in
    basis["configuration"]["heat_out_valve_placement"] = heat_out

    # Keep exchanger-local descriptions consistent with configuration.
    basis["heat_in"]["valve_placement"] = heat_in
    basis["heat_out"]["valve_placement"] = heat_out

    basis_text = json.dumps(basis, indent=2) + "\n"
    basis_sha = hashlib.sha256(basis_text.encode()).hexdigest()

    raw["sources"]["machine"]["path"] = basis_path.name
    raw["sources"]["machine"]["sha256"] = basis_sha

    raw["study"]["name"] = (
        f"5 L / 2.8 Hz valve screen — {label} — {code}"
    )
    raw["study"]["purpose"] = (
        "frozen_candidate_passive_check_valve_placement_screen"
    )

    # Explicitly keep the same neutral starting policy for every topology.
    raw["warm_start"] = {"initial_source": "uniform"}

    basis_path.write_text(basis_text)
    study_path.write_text(dumps(raw))

    return study_path, target / "evaluation.json"


def fmt(x, digits=4):
    if x is None:
        return "—"
    return f"{x:.{digits}f}"


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)

    candidates = select_candidates()
    results = []

    print("Selected candidates:")
    for label, record in candidates:
        m = record.get("metrics", {})
        print(
            f"  {label:7s} {record['candidate_id']} "
            f"COP={m.get('cooling_cop')} "
            f"Qcold={m.get('cooling_power_w')} W"
        )

    print()
    print("Running valve-placement screen...")

    for label, source in candidates:
        for code, (heat_in, heat_out) in PLACEMENTS.items():
            study_path, evaluation_path = make_study(
                label, source, code, heat_in, heat_out
            )

            if evaluation_path.exists():
                artifact = json.loads(evaluation_path.read_text())
                record = artifact["record"]
            else:
                print(
                    f"{label:7s} {source['candidate_id'][:12]} "
                    f"{code} ...",
                    flush=True,
                )
                record = evaluate(
                    study_path,
                    evaluation_path,
                    budget="5m",
                )

            metrics = record.get("metrics", {})
            validity = metrics.get("validity") or {}

            row = {
                "label": label,
                "source_candidate_id": source["candidate_id"],
                "topology": code,
                "heat_in_valve_placement": heat_in,
                "heat_out_valve_placement": heat_out,
                "status": record.get("status"),
                "cooling_cop": metrics.get("cooling_cop"),
                "cooling_power_w": metrics.get("cooling_power_w"),
                "mechanical_input_w":
                    metrics.get("indicated_mechanical_input_power_w"),
                "maximum_pressure_pa":
                    metrics.get("maximum_pressure_pa"),
                "maximum_temperature_k":
                    metrics.get("maximum_temperature_k"),
                "maximum_absolute_mass_flow_kg_s":
                    metrics.get("maximum_absolute_mass_flow_kg_s"),
                "maximum_mach_number":
                    validity.get("maximum_mach_number"),
            }

            results.append(row)

    (OUTPUT / "summary.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )

    print()
    print(
        "candidate  topo  status                  "
        "COP       Qcold W    Pin W     "
        "Pmax bar   mdot kg/s   Mach"
    )
    print("-" * 113)

    for r in results:
        print(
            f"{r['label']:9s} "
            f"{r['topology']:4s} "
            f"{r['status'][:22]:22s} "
            f"{fmt(r['cooling_cop'], 4):>8s} "
            f"{fmt(r['cooling_power_w'], 2):>10s} "
            f"{fmt(r['mechanical_input_w'], 2):>9s} "
            f"{fmt(r['maximum_pressure_pa'] / 1e5 if r['maximum_pressure_pa'] else None, 3):>10s} "
            f"{fmt(r['maximum_absolute_mass_flow_kg_s'], 5):>11s} "
            f"{fmt(r['maximum_mach_number'], 4):>8s}"
        )

    print()
    print("Best topology per candidate:")

    for label, source in candidates:
        rows = [
            r for r in results
            if r["source_candidate_id"] == source["candidate_id"]
            and r["cooling_cop"] is not None
        ]

        if not rows:
            print(f"  {label}: no refrigeration solution")
            continue

        best = max(rows, key=lambda r: r["cooling_cop"])
        print(
            f"  {label:7s}: {best['topology']} "
            f"COP={best['cooling_cop']:.5f}, "
            f"Qcold={best['cooling_power_w']:.2f} W"
        )

    print()
    print(f"Full results: {OUTPUT / 'summary.json'}")


if __name__ == "__main__":
    main()
