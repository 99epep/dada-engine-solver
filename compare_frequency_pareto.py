#!/usr/bin/env python3
"""Compare COP-Qcold Pareto fronts across DADA Research frequency campaigns.

Expected layout (default):
    outputs/pedal_cell_retro/compact/frequency/
        0p3Hz/campaign/
        0p4Hz/campaign/
        0p6Hz/campaign2/   # explicit built-in override
        ...

New frequency directories are discovered automatically.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path


DEFAULT_ROOT = Path("outputs/pedal_cell_retro/compact/frequency")
DEFAULT_OVERRIDES = {"0p6Hz": "campaign2"}
FREQUENCY_RE = re.compile(r"^(?P<value>[0-9]+(?:p[0-9]+)?)Hz$")


def find_repo_root(start: Path) -> Path:
    """Find a checkout containing src/dada_solver, starting from script and cwd."""
    candidates = [start.resolve(), Path.cwd().resolve()]
    seen: set[Path] = set()
    for base in candidates:
        for p in (base, *base.parents):
            if p in seen:
                continue
            seen.add(p)
            if (p / "src" / "dada_solver").is_dir():
                return p
    raise SystemExit(
        "Cannot find repository root (expected src/dada_solver). "
        "Run this script from the dada-engine-solver checkout."
    )


REPO_ROOT = find_repo_root(Path(__file__).resolve().parent)
sys.path.insert(0, str(REPO_ROOT / "src"))

from dada_solver.research.report import inspect  # noqa: E402


def frequency_from_name(name: str) -> float | None:
    match = FREQUENCY_RE.fullmatch(name)
    if not match:
        return None
    return float(match.group("value").replace("p", "."))


def parse_override(text: str) -> tuple[str, str]:
    try:
        frequency_dir, campaign = text.split("=", 1)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("expected FREQUENCY_DIR=CAMPAIGN, e.g. 0p6Hz=campaign2") from exc
    if frequency_from_name(frequency_dir) is None or not campaign:
        raise argparse.ArgumentTypeError("expected e.g. 0p6Hz=campaign2")
    return frequency_dir, campaign


def choose_campaign(frequency_dir: Path, overrides: dict[str, str]) -> Path | None:
    """Choose the campaign to compare for one frequency directory."""
    if frequency_dir.name in overrides:
        chosen = frequency_dir / overrides[frequency_dir.name]
        if not chosen.is_dir():
            raise SystemExit(f"Configured campaign does not exist: {chosen}")
        return chosen

    # Normal convention: campaign/.  This deliberately wins over campaign2, etc.
    default = frequency_dir / "campaign"
    if default.is_dir():
        return default

    # Fallback for future directories that contain only one campaign-like directory.
    campaigns = sorted(
        p for p in frequency_dir.iterdir()
        if p.is_dir() and p.name.startswith("campaign")
    )
    if len(campaigns) == 1:
        return campaigns[0]
    if not campaigns:
        return None

    names = ", ".join(p.name for p in campaigns)
    raise SystemExit(
        f"Several campaigns found in {frequency_dir} ({names}) but no 'campaign/' default. "
        f"Select one with --campaign {frequency_dir.name}=NAME."
    )


def finite_number(value) -> float | None:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def feasible_points(campaign: Path) -> list[dict]:
    """Return distinct feasible candidates carrying finite COP and Qcold."""
    data = inspect(campaign)
    distinct: dict[str, dict] = {}
    for record in data["records"]:
        if record.get("status") != "feasible":
            continue
        metrics = record.get("metrics") or {}
        cop = finite_number(metrics.get("cooling_cop"))
        qcold = finite_number(metrics.get("cooling_power_w"))
        if cop is None or qcold is None:
            continue
        distinct[record["candidate_id"]] = {
            "candidate_id": record["candidate_id"],
            "cop": cop,
            "qcold_w": qcold,
        }
    return list(distinct.values())


def pareto_front(points: list[dict]) -> list[dict]:
    """Return non-dominated points when COP and Qcold are both maximized."""
    # Descending Qcold: a point is non-dominated iff its COP exceeds every point
    # already seen at equal or higher Qcold.  Equal COP keeps only the larger Qcold.
    ordered = sorted(points, key=lambda p: (-p["qcold_w"], -p["cop"]))
    front: list[dict] = []
    best_cop = -math.inf
    for point in ordered:
        if point["cop"] > best_cop:
            front.append(point)
            best_cop = point["cop"]
    return sorted(front, key=lambda p: p["qcold_w"])


def discover(root: Path, overrides: dict[str, str]) -> list[dict]:
    series = []
    if not root.is_dir():
        raise SystemExit(f"Frequency root does not exist: {root}")

    for frequency_dir in root.iterdir():
        if not frequency_dir.is_dir():
            continue
        frequency_hz = frequency_from_name(frequency_dir.name)
        if frequency_hz is None:
            continue
        campaign = choose_campaign(frequency_dir, overrides)
        if campaign is None:
            print(f"warning: no campaign found in {frequency_dir}; skipped", file=sys.stderr)
            continue
        points = feasible_points(campaign)
        if not points:
            print(f"warning: no feasible COP/Qcold points in {campaign}; skipped", file=sys.stderr)
            continue
        series.append(
            {
                "frequency_hz": frequency_hz,
                "frequency_dir": frequency_dir.name,
                "campaign": campaign,
                "points": points,
                "front": pareto_front(points),
            }
        )

    series.sort(key=lambda s: s["frequency_hz"])
    return series


def write_csv(series: list[dict], path: Path) -> None:
    import csv

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["frequency_hz", "campaign", "candidate_id", "cooling_power_w", "cooling_cop"])
        for s in series:
            for p in s["front"]:
                writer.writerow([
                    s["frequency_hz"],
                    s["campaign"],
                    p["candidate_id"],
                    p["qcold_w"],
                    p["cop"],
                ])


def plot(series: list[dict], output: Path, show_cloud: bool) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit("matplotlib is required to draw the comparison") from exc

    fig, ax = plt.subplots(figsize=(11, 7))

    for s in series:
        label = f"{s['frequency_hz']:g} Hz"
        front = s["front"]

        if show_cloud:
            ax.scatter(
                [p["qcold_w"] for p in s["points"]],
                [p["cop"] for p in s["points"]],
                s=9,
                alpha=0.10,
            )

        ax.plot(
            [p["qcold_w"] for p in front],
            [p["cop"] for p in front],
            marker="o",
            markersize=4,
            linewidth=1.8,
            label=label,
        )

    ax.set_xlabel("Qcold [W]")
    ax.set_ylabel("Cooling COP [-]")
    ax.set_title("COP–Qcold Pareto fronts by operating frequency")
    ax.grid(True, alpha=0.25)
    ax.legend(title="Frequency", ncols=2)
    fig.tight_layout()

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    print(f"plot: {output}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help=f"frequency directory (default: {DEFAULT_ROOT})",
    )
    parser.add_argument(
        "--campaign",
        action="append",
        default=[],
        type=parse_override,
        metavar="FREQUENCY_DIR=CAMPAIGN",
        help="override campaign selection; may be repeated (e.g. 0p6Hz=campaign2)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_ROOT / "pareto_cop_qcold_by_frequency.png",
        help="output image (.png, .svg, .pdf, ...)",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=DEFAULT_ROOT / "pareto_cop_qcold_by_frequency.csv",
        help="CSV containing every Pareto point",
    )
    parser.add_argument(
        "--cloud",
        action="store_true",
        help="also draw all feasible candidates faintly behind each Pareto front",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    overrides = dict(DEFAULT_OVERRIDES)
    overrides.update(dict(args.campaign))

    root = args.root if args.root.is_absolute() else REPO_ROOT / args.root
    output = args.output if args.output.is_absolute() else REPO_ROOT / args.output
    csv_path = args.csv if args.csv.is_absolute() else REPO_ROOT / args.csv

    series = discover(root, overrides)
    if not series:
        raise SystemExit(f"No usable frequency campaigns found under {root}")

    print("Compared campaigns:")
    for s in series:
        print(
            f"  {s['frequency_hz']:g} Hz: {s['campaign']}  "
            f"({len(s['points'])} feasible points, {len(s['front'])} Pareto points)"
        )

    write_csv(series, csv_path)
    print(f"csv:  {csv_path}")
    plot(series, output, args.cloud)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
