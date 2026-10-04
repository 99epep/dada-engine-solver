#!/usr/bin/env python3

from __future__ import annotations

import copy
import hashlib
import json
import tomllib
from pathlib import Path

from dada_solver.research.cli import evaluate
from dada_solver.research.study_io import dumps


SOURCE_DIR = Path("outputs/pedal_cell_retro/compact/final")
OUTPUT_DIR = Path("outputs/pedal_cell_retro/compact/final_valve_screen")

CHAMPIONS = ("DD", "DU", "UD", "UU")

# Code = placement Hi, puis placement Ho.
PLACEMENTS = {
    "DD": ("downstream", "downstream"),
    "DU": ("downstream", "upstream"),
    "UD": ("upstream", "downstream"),
    "UU": ("upstream", "upstream"),
}


def load_toml(path: Path) -> dict:
    with path.open("rb") as fh:
        return tomllib.load(fh)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def make_study(
    champion: str,
    topology: str,
) -> tuple[Path, Path]:
    source_study_path = SOURCE_DIR / f"{champion}.toml"
    source_basis_path = SOURCE_DIR / f"{champion}.basis.json"

    if not source_study_path.exists():
        raise FileNotFoundError(source_study_path)
    if not source_basis_path.exists():
        raise FileNotFoundError(source_basis_path)

    raw = copy.deepcopy(load_toml(source_study_path))
    basis = copy.deepcopy(load_json(source_basis_path))

    heat_in, heat_out = PLACEMENTS[topology]

    # Ne modifier QUE le placement des clapets.
    basis["configuration"]["heat_in_valve_placement"] = heat_in
    basis["configuration"]["heat_out_valve_placement"] = heat_out

    basis["heat_in"]["valve_placement"] = heat_in
    basis["heat_out"]["valve_placement"] = heat_out

    case_dir = OUTPUT_DIR / champion / topology
    case_dir.mkdir(parents=True, exist_ok=True)

    basis_path = case_dir / "machine.basis.json"
    study_path = case_dir / "study.toml"

    basis_bytes = (
        json.dumps(
            basis,
            indent=2,
            sort_keys=False,
        )
        + "\n"
    ).encode("utf-8")

    basis_path.write_bytes(basis_bytes)

    # Le study généré doit pointer sur sa propre basis.
    raw["sources"]["machine"]["path"] = basis_path.name
    raw["sources"]["machine"]["sha256"] = sha256_bytes(basis_bytes)

    # Même politique de démarrage pour les 16 cas.
    raw["warm_start"] = {
        "initial_source": "uniform",
    }

    study_path.write_text(dumps(raw), encoding="utf-8")

    return study_path, case_dir / "evaluation.json"


def metric(record: dict, key: str):
    value = record.get("metrics", {}).get(key)
    if value is not None:
        return value

    value = record.get("derived", {}).get(key)
    if value is not None:
        return value

    return record.get(key)


def run_case(champion: str, topology: str) -> dict:
    study_path, evaluation_path = make_study(champion, topology)

    if evaluation_path.exists():
        print(f"{champion} × {topology}: reuse {evaluation_path}")
        record = load_json(evaluation_path)
    else:
        print(f"{champion} × {topology}: evaluate")
        record = evaluate(
            study_path,
            evaluation_path,
            budget="5m",
        )

    validity = record.get("validity", {})
    if not isinstance(validity, dict):
        validity = {}

    return {
        "champion": champion,
        "topology": topology,
        "status": record.get("status"),
        "cooling_cop": metric(record, "cooling_cop"),
        "cooling_power_w": metric(record, "cooling_power_w"),
        "mechanical_input_power_w": metric(
            record,
            "indicated_mechanical_input_power_w",
        ),
        "maximum_pressure_pa": metric(
            record,
            "maximum_pressure_pa",
        ),
        "maximum_temperature_k": metric(
            record,
            "maximum_temperature_k",
        ),
        "maximum_absolute_mass_flow_kg_s": metric(
            record,
            "maximum_absolute_mass_flow_kg_s",
        ),
        "maximum_mach_number": validity.get("maximum_mach_number"),
        "evaluation": str(evaluation_path),
    }


def fmt(value, digits=4):
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def print_cop_matrix(results: list[dict]) -> None:
    by_case = {
        (r["champion"], r["topology"]): r
        for r in results
    }

    print()
    print("COP matrix")
    print("source \\ valves | " + " | ".join(f"{x:>8}" for x in PLACEMENTS))
    print("-" * 58)

    for champion in CHAMPIONS:
        values = []
        for topology in PLACEMENTS:
            r = by_case[(champion, topology)]
            cop = r["cooling_cop"]
            values.append(
                f"{cop:8.4f}" if isinstance(cop, (int, float)) else f"{'-':>8}"
            )
        print(f"{champion:15} | " + " | ".join(values))


def print_power_matrix(results: list[dict]) -> None:
    by_case = {
        (r["champion"], r["topology"]): r
        for r in results
    }

    print()
    print("Cooling power matrix [W]")
    print("source \\ valves | " + " | ".join(f"{x:>8}" for x in PLACEMENTS))
    print("-" * 58)

    for champion in CHAMPIONS:
        values = []
        for topology in PLACEMENTS:
            r = by_case[(champion, topology)]
            power = r["cooling_power_w"]
            values.append(
                f"{power:8.2f}"
                if isinstance(power, (int, float))
                else f"{'-':>8}"
            )
        print(f"{champion:15} | " + " | ".join(values))


def print_row_winners(results: list[dict]) -> None:
    print()
    print("Best valve topology for each frozen champion")

    for champion in CHAMPIONS:
        rows = [
            r
            for r in results
            if r["champion"] == champion
            and r["status"] == "feasible"
            and isinstance(r["cooling_cop"], (int, float))
        ]

        if not rows:
            print(f"{champion}: no feasible result")
            continue

        best = max(rows, key=lambda r: r["cooling_cop"])

        source = next(
            (
                r
                for r in rows
                if r["topology"] == champion
            ),
            None,
        )

        delta = None
        if source and source["cooling_cop"]:
            delta = (
                100.0
                * (
                    best["cooling_cop"]
                    - source["cooling_cop"]
                )
                / source["cooling_cop"]
            )

        text = (
            f"{champion}: {best['topology']} "
            f"COP={fmt(best['cooling_cop'])}, "
            f"Qcold={fmt(best['cooling_power_w'], 2)} W"
        )

        if delta is not None:
            text += f", ΔCOP vs original={delta:+.3f}%"

        print(text)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    results = []

    for champion in CHAMPIONS:
        for topology in PLACEMENTS:
            result = run_case(champion, topology)
            results.append(result)

            print(
                f"  status={result['status']} "
                f"COP={fmt(result['cooling_cop'])} "
                f"Qcold={fmt(result['cooling_power_w'], 2)} W "
                f"Pin={fmt(result['mechanical_input_power_w'], 2)} W"
            )

    summary_path = OUTPUT_DIR / "summary.json"
    summary_path.write_text(
        json.dumps(results, indent=2) + "\n",
        encoding="utf-8",
    )

    print_cop_matrix(results)
    print_power_matrix(results)
    print_row_winners(results)

    feasible = [
        r
        for r in results
        if r["status"] == "feasible"
        and isinstance(r["cooling_cop"], (int, float))
    ]

    if feasible:
        best = max(feasible, key=lambda r: r["cooling_cop"])

        print()
        print(
            "Global best: "
            f"champion {best['champion']} × "
            f"valves {best['topology']} | "
            f"COP={fmt(best['cooling_cop'])} | "
            f"Qcold={fmt(best['cooling_power_w'], 2)} W | "
            f"Pin={fmt(best['mechanical_input_power_w'], 2)} W"
        )

    print()
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
