#!/usr/bin/env python3
"""Verify the retained paired six-bar mechanism catalogue."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import tomllib

from dada_solver.six_bar import SixBarCylinderMechanism
from dada_solver.mechanism_diagnostics import six_bar_metrics
from dada_solver.research.margins import validate_mechanical_constraint


HERE = Path(__file__).resolve().parent
FAMILIES_PATH = HERE / "families.toml"
SCREEN_PATH = HERE / "design_screen.toml"

MECHANISM_FIELDS = (
    "primary_ground",
    "primary_coupler",
    "primary_rocker",
    "primary_e_along",
    "primary_e_normal",
    "primary_phase",
    "second_pivot_x",
    "second_pivot_y",
    "link_ef",
    "link_gf",
    "h_along_over_ef",
    "h_normal_over_ef",
    "piston_rod",
    "slider_axis_offset",
    "slider_axis_angle",
    "primary_branch",
    "second_branch",
)


def load_toml(path: Path) -> dict:
    with path.open("rb") as stream:
        return tomllib.load(stream)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def build_mechanism(data: dict) -> SixBarCylinderMechanism:
    missing = [name for name in MECHANISM_FIELDS if name not in data]
    extra = sorted(set(data) - set(MECHANISM_FIELDS))
    require(not missing, f"Missing mechanism fields: {missing}")
    require(not extra, f"Unexpected mechanism fields: {extra}")
    for name in ("primary_branch", "second_branch"):
        require(type(data[name]) is int and data[name] in (-1, 1),
                f"{name}: expected an integer assembly branch -1 or +1.")
    return SixBarCylinderMechanism(
        **{
            name: int(data[name])
            if name in ("primary_branch", "second_branch")
            else float(data[name])
            for name in MECHANISM_FIELDS
        }
    )


def constraint_margin(value: float, relation: str, limit: float):
    if relation == "minimum":
        return float(value) - float(limit)
    if relation == "maximum":
        return float(limit) - float(value)
    if relation == "equal":
        return None
    raise ValueError(f"Unknown relation {relation!r}")


def constraint_passes(value: float, relation: str, limit: float) -> bool:
    if relation == "minimum":
        return float(value) >= float(limit)
    if relation == "maximum":
        return float(value) <= float(limit)
    if relation == "equal":
        return math.isclose(float(value), float(limit), rel_tol=0.0, abs_tol=1e-12)
    raise ValueError(f"Unknown relation {relation!r}")


def validate_screen(screen: dict) -> None:
    require(screen.get("schema_version") == 1, "Unsupported screen schema.")
    require(int(screen["canonical_samples"]) >= 360, "Canonical sample count is too small.")
    require(
        int(screen["dense_cross_check_samples"]) > int(screen["canonical_samples"]),
        "Dense cross-check must use more samples than the canonical screen.",
    )
    require(bool(screen["constraints"]), "The screen must contain constraints.")
    for row in screen["constraints"]:
        validate_mechanical_constraint(row, "six_bar", scoped=False)


def evaluate_catalogue(families_data: dict, screen: dict):
    require(families_data.get("schema_version") == 1, "Unsupported family schema.")
    require(families_data.get("length_unit") == "crank_radius", "Expected crank-radius lengths.")
    require(families_data.get("angle_unit") == "rad", "Expected radian angles.")
    validate_screen(screen)
    order = families_data["family_order"]
    families = families_data["families"]

    require(len(order) == 4, "Expected four retained family lineages.")
    require(len(set(order)) == len(order), "Family identifiers must be unique.")
    require(set(order) == set(families), "family_order and family tables disagree.")

    canonical_samples = int(screen["canonical_samples"])
    dense_samples = int(screen["dense_cross_check_samples"])
    results = {}

    for family_id in order:
        family = families[family_id]
        results[family_id] = {}
        require(
            family.get("primary_seed_family") == family_id,
            f"{family_id}: lineage identifier mismatch.",
        )
        for side in ("small", "large"):
            mechanism = build_mechanism(family[side])
            canonical = six_bar_metrics(mechanism, samples=canonical_samples)
            dense = six_bar_metrics(mechanism, samples=dense_samples)

            rows = []
            for constraint in screen["constraints"]:
                metric = constraint["metric"]
                relation = constraint["relation"]
                limit = constraint["limit"]
                require(metric in canonical, f"{family_id}/{side}: missing metric {metric}.")
                require(metric in dense, f"{family_id}/{side}: dense metric {metric} missing.")
                c_value = canonical[metric]
                d_value = dense[metric]
                rows.append(
                    {
                        **constraint,
                        "canonical": c_value,
                        "dense": d_value,
                        "margin": constraint_margin(c_value, relation, limit),
                        "dense_margin": constraint_margin(d_value, relation, limit),
                        "canonical_pass": constraint_passes(c_value, relation, limit),
                        "dense_pass": constraint_passes(d_value, relation, limit),
                    }
                )

            results[family_id][side] = {
                "mechanism": mechanism,
                "canonical": canonical,
                "dense": dense,
                "constraints": rows,
            }
    return results


def screen_failures(results: dict) -> list[str]:
    """Collect all failures, retaining both grid values rather than stopping early."""
    failures = []
    for family_id, sides in results.items():
        for side, record in sides.items():
            for row in record["constraints"]:
                if not row["canonical_pass"] or not row["dense_pass"]:
                    failures.append(
                        f"{family_id}/{side}: {row['metric']} {row['relation']} {row['limit']!r}; "
                        f"canonical={row['canonical']!r} pass={row['canonical_pass']}; "
                        f"dense={row['dense']!r} pass={row['dense_pass']}"
                    )
    return failures


def format_number(value) -> str:
    if isinstance(value, int):
        return str(value)
    return f"{float(value):.6f}"


def format_margin(value) -> str:
    return "—" if value is None else f"{value:+.9g}"


def markdown_tables(results: dict) -> str:
    wanted = (
        "stroke_over_crank",
        "minimum_primary_transmission_sine",
        "minimum_secondary_transmission_sine",
        "minimum_rod_axis_cosine",
        "EH_over_crank",
        "H_axis_lateral_rms_over_stroke",
        "H_axis_lateral_span_over_stroke",
        "crank_axis_to_EFH_clearance_over_crank",
        "zero_crossing_count",
    )
    lines = []
    for family_id, sides in results.items():
        lines.extend((f"### `{family_id}`", "", "| Metric | SMALL | LARGE |", "|---|---:|---:|"))
        for metric in wanted:
            small = sides["small"]["canonical"][metric]
            large = sides["large"]["canonical"][metric]
            lines.append(f"| `{metric}` | {format_number(small)} | {format_number(large)} |")
        lines.extend(("", "Signed margins (canonical / dense); equality rows report pass/fail.", "",
                      "| Constraint | Unit | SMALL margin | LARGE margin |",
                      "|---|---|---:|---:|"))
        for small, large in zip(sides["small"]["constraints"], sides["large"]["constraints"]):
            symbol = {"minimum": "≥", "maximum": "≤", "equal": "="}[small["relation"]]
            cells = []
            for row in (small, large):
                if row["relation"] == "equal":
                    cells.append(" / ".join("pass" if row[k] else "FAIL"
                                            for k in ("canonical_pass", "dense_pass")))
                else:
                    cells.append(f"{format_margin(row['margin'])} / {format_margin(row['dense_margin'])}")
            lines.append(f"| `{small['metric']}` {symbol} {small['limit']} | {small['unit']} | "
                         f"{cells[0]} | {cells[1]} |")
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--check", action="store_true", help="Fail if either nominal screen fails.")
    modes.add_argument("--markdown", action="store_true", help="Print metric and individual-margin tables.")
    args = parser.parse_args()

    families_data = load_toml(FAMILIES_PATH)
    screen = load_toml(SCREEN_PATH)
    results = evaluate_catalogue(families_data, screen)
    failures = screen_failures(results)

    if args.markdown:
        print(markdown_tables(results))
    elif not args.check:
        print("Paired six-bar mechanism catalogue")
        print(f"families={len(results)} mechanisms={2 * len(results)} "
              f"canonical_samples={screen['canonical_samples']} "
              f"dense_samples={screen['dense_cross_check_samples']}")
        for family_id, sides in results.items():
            for side, record in sides.items():
                print(f"\n{family_id}/{side}")
                for row in record["constraints"]:
                    print(f"  {row['metric']} {row['relation']} {row['limit']} [{row['unit']}]: "
                          f"canonical={row['canonical']:.12g} margin={format_margin(row['margin'])} "
                          f"pass={row['canonical_pass']}; "
                          f"dense={row['dense']:.12g} margin={format_margin(row['dense_margin'])} "
                          f"pass={row['dense_pass']}")

    if failures:
        raise SystemExit("Screen failures:\n" + "\n".join(failures))
    if args.check:
        print("OK six_bar_mechanism_families "
              f"families={len(results)} mechanisms={2 * len(results)} "
              f"samples={screen['canonical_samples']}/{screen['dense_cross_check_samples']}")


if __name__ == "__main__":
    main()
