"""Plot Stage F1 FreeKinematics motion laws over the harmonic reference.

Defaults: every feasible F1 candidate that beats the FreeKinematics harmonic
baseline, plus the exact 250-degree harmonic and the 12-control F1 baseline.

Outputs:
    outputs/motor_free_local_stageF1_motion_overlay.png
    outputs/motor_free_local_stageF1_motion_delta.png
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import tempfile

import numpy as np

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.kinematics import HarmonicVolumeKinematics

from compare_motor_motion_laws_stage7A5 import (
    ROOT,
    A5_CAMPAIGN,
    _candidate_design_and_mass,
)
from optimize_motor_free_local_stageF1 import (
    ANGLES_DEG,
    BASELINE,
    PHASE_DEG,
    _make_kinematics,
)

F1_CSV = ROOT / "outputs" / "motor_free_local_stageF1.csv"
OUTPUT_OVERLAY = ROOT / "outputs" / "motor_free_local_stageF1_motion_overlay.png"
OUTPUT_DELTA = ROOT / "outputs" / "motor_free_local_stageF1_motion_delta.png"


def setup_matplotlib():
    cache = Path(tempfile.gettempdir()) / "dada_solver_matplotlib"
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib.pyplot as plt
    return plt


def load_rows():
    with F1_CSV.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    numeric = (
        "index", "small_60_deg", "small_90_deg",
        "large_150_deg", "large_210_deg",
        "indicated_thermal_efficiency", "indicated_power_w",
    )
    for row in rows:
        for key in numeric:
            row[key] = float(row[key])
        row["index"] = int(row["index"])
        row["feasible"] = row["feasible"].strip().lower() == "true"
    return rows


def candidate_values(row):
    return {
        "small_60_deg": row["small_60_deg"],
        "small_90_deg": row["small_90_deg"],
        "large_150_deg": row["large_150_deg"],
        "large_210_deg": row["large_210_deg"],
    }


def motion(kinematics, angles):
    small = np.array([
        kinematics.small_cylinder_volume(float(theta))
        for theta in angles
    ])
    large = np.array([
        kinematics.large_cylinder_volume(float(theta))
        for theta in angles
    ])
    return small, large


def legend_below(axis, ncol=2):
    axis.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=ncol,
        frameon=False,
    )


def decorate(axis):
    axis.set_xlim(0.0, 360.0)
    axis.set_xticks(np.arange(0.0, 361.0, 45.0))
    axis.grid(True, alpha=0.3)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--top",
        type=int,
        default=None,
        help=(
            "Show N best feasible candidates instead of only candidates "
            "that beat the F1 harmonic baseline."
        ),
    )
    args = parser.parse_args()

    definition = CampaignDefinition(A5_CAMPAIGN)
    base_design, _mass, _saved, _best, _candidate_id = (
        _candidate_design_and_mass(definition)
    )
    small_limits = base_design.configuration.machine_volumes.small_cylinder
    large_limits = base_design.configuration.machine_volumes.large_cylinder

    rows = load_rows()
    baseline_row = next(
        row for row in rows if row["label"] == "free_harmonic_baseline"
    )
    baseline_eta = baseline_row["indicated_thermal_efficiency"]

    feasible = sorted(
        (
            row for row in rows
            if row["feasible"]
            and math.isfinite(row["indicated_thermal_efficiency"])
        ),
        key=lambda row: row["indicated_thermal_efficiency"],
        reverse=True,
    )

    if args.top is None:
        selected = [
            row for row in feasible
            if row["indicated_thermal_efficiency"] > baseline_eta
        ]
    else:
        if args.top < 1:
            raise ValueError("--top must be >= 1")
        selected = feasible[:args.top]

    exact_harmonic = HarmonicVolumeKinematics(
        small_limits,
        large_limits,
        math.radians(PHASE_DEG),
    )
    free_baseline = _make_kinematics(
        small_limits,
        large_limits,
        dict(BASELINE),
    )

    angle_deg = np.linspace(0.0, 360.0, 1441)
    angle_rad = np.radians(angle_deg)

    exact_s, exact_l = motion(exact_harmonic, angle_rad)
    base_s, base_l = motion(free_baseline, angle_rad)

    candidates = []
    for rank, row in enumerate(selected, 1):
        kin = _make_kinematics(
            small_limits,
            large_limits,
            candidate_values(row),
        )
        small, large = motion(kin, angle_rad)
        candidates.append((rank, row, small, large))

    plt = setup_matplotlib()

    # Physical-volume overlay.
    fig, axes = plt.subplots(2, 1, figsize=(12, 11), sharex=True)

    axes[0].plot(
        angle_deg, exact_s * 1000.0,
        linestyle="--", linewidth=2.0,
        label=f"Exact harmonic {PHASE_DEG:.0f}°",
    )
    axes[0].plot(
        angle_deg, base_s * 1000.0,
        linestyle=":", linewidth=2.0,
        label=f"F1 spline baseline ({100*baseline_eta:.4f}%)",
    )
    for rank, row, small, _large in candidates:
        axes[0].plot(
            angle_deg, small * 1000.0,
            label=(
                f"#{rank} index {row['index']} — "
                f"{100*row['indicated_thermal_efficiency']:.4f}%"
            ),
        )
    axes[0].set_ylabel("Small-cylinder volume (L)")
    legend_below(axes[0], ncol=2)

    axes[1].plot(
        angle_deg, exact_l * 1000.0,
        linestyle="--", linewidth=2.0,
        label=f"Exact harmonic {PHASE_DEG:.0f}°",
    )
    axes[1].plot(
        angle_deg, base_l * 1000.0,
        linestyle=":", linewidth=2.0,
        label=f"F1 spline baseline ({100*baseline_eta:.4f}%)",
    )
    for rank, row, _small, large in candidates:
        axes[1].plot(
            angle_deg, large * 1000.0,
            label=(
                f"#{rank} index {row['index']} — "
                f"{100*row['indicated_thermal_efficiency']:.4f}%"
            ),
        )
    axes[1].set_ylabel("Large-cylinder volume (L)")
    axes[1].set_xlabel("Cycle angle (deg)")
    legend_below(axes[1], ncol=2)

    for axis in axes:
        for angle in ANGLES_DEG:
            axis.axvline(angle, linewidth=0.45, alpha=0.12)
        decorate(axis)

    for angle in (60.0, 90.0):
        axes[0].axvline(angle, linestyle="--", linewidth=0.8, alpha=0.5)
    for angle in (150.0, 210.0):
        axes[1].axvline(angle, linestyle="--", linewidth=0.8, alpha=0.5)

    fig.suptitle(
        "Stage F1 — FreeKinematics motions over the 250° harmonic reference\n"
        f"{len(selected)} improving candidate(s)",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.97), h_pad=3.1)
    OUTPUT_OVERLAY.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_OVERLAY, dpi=160)
    plt.close(fig)

    # Difference from exact harmonic in mL.
    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    axes[0].axhline(0.0, linewidth=1.0)
    axes[1].axhline(0.0, linewidth=1.0)

    axes[0].plot(
        angle_deg, (base_s - exact_s) * 1e6,
        linestyle=":", linewidth=2.0,
        label="F1 spline baseline - exact harmonic",
    )
    axes[1].plot(
        angle_deg, (base_l - exact_l) * 1e6,
        linestyle=":", linewidth=2.0,
        label="F1 spline baseline - exact harmonic",
    )

    for rank, row, small, large in candidates:
        label = (
            f"#{rank} index {row['index']} — "
            f"{100*row['indicated_thermal_efficiency']:.4f}%"
        )
        axes[0].plot(angle_deg, (small - exact_s) * 1e6, label=label)
        axes[1].plot(angle_deg, (large - exact_l) * 1e6, label=label)

    axes[0].set_ylabel("Small V - harmonic (mL)")
    axes[1].set_ylabel("Large V - harmonic (mL)")
    axes[1].set_xlabel("Cycle angle (deg)")

    for axis in axes:
        for angle in ANGLES_DEG:
            axis.axvline(angle, linewidth=0.45, alpha=0.12)
        decorate(axis)
        legend_below(axis, ncol=2)

    for angle in (60.0, 90.0):
        axes[0].axvline(angle, linestyle="--", linewidth=0.8, alpha=0.5)
    for angle in (150.0, 210.0):
        axes[1].axvline(angle, linestyle="--", linewidth=0.8, alpha=0.5)

    fig.suptitle(
        "Stage F1 — departure from exact harmonic motion",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.97), h_pad=3.1)
    fig.savefig(OUTPUT_DELTA, dpi=160)
    plt.close(fig)

    print(
        f"F1 spline baseline: {100*baseline_eta:.6f}% "
        f"({baseline_row['indicated_power_w']:.3f} W)"
    )
    print(f"Candidates shown: {len(selected)}")
    for rank, row, _small, _large in candidates:
        print(
            f"  #{rank} index={row['index']} "
            f"eta={100*row['indicated_thermal_efficiency']:.6f}% "
            f"P={row['indicated_power_w']:.3f} W"
        )

    print(f"\nWrote {OUTPUT_OVERLAY.relative_to(ROOT)}")
    print(f"Wrote {OUTPUT_DELTA.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
