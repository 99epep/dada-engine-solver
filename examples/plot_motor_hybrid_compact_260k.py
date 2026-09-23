#!/usr/bin/env python3
"""Plot the 260 K hybrid-compact champion against the aligned 16-control source."""

from pathlib import Path
import argparse
import json

import matplotlib.pyplot as plt
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--report",
        type=Path,
        default=Path("outputs/motor_hybrid_compact_260k/report.json"),
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/motor_hybrid_compact_260k/comparison.png"),
    )
    args = ap.parse_args()

    report = json.loads(args.report.read_text())
    best = report["best_feasible"]
    p = best["parameters"]

    csv_path = args.report.parent / "best_motion.csv"
    d = np.genfromtxt(csv_path, delimiter=",", names=True)

    deg = d["motor_angle_deg"]
    ss = d["source16_aligned_small"]
    sl = d["source16_aligned_large"]
    hs = d["hybrid_small"]
    hl = d["hybrid_large"]
    dss = d["source16_aligned_dsmall_dt"]
    dsl = d["source16_aligned_dlarge_dt"]
    dhs = d["hybrid_dsmall_dt"]
    dhl = d["hybrid_dlarge_dt"]

    smax = p["small_max_deg"] % 360.0
    smin = (smax + p["small_down_duration_deg"]) % 360.0
    lmax = 0.0
    lmin = p["large_down_duration_deg"] % 360.0

    skink = (
        smax
        + p["small_down_kink_u"] * p["small_down_duration_deg"]
    ) % 360.0
    lup = 360.0 - p["large_down_duration_deg"]
    lkink = (lmin + p["large_up_kink_u"] * lup) % 360.0

    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)

    ax = axes[0, 0]
    ax.plot(deg, ss, "--", label="Source16 aligned")
    ax.plot(deg, hs, label="Hybrid champion")
    for x in (smax, smin):
        ax.axvline(x, linewidth=0.9, linestyle=":")
    ax.axvline(skink, linewidth=1.0, linestyle="-.")
    ax.set_title("Small piston — normalized volume")
    ax.set_ylabel("q")
    ax.legend()
    ax.grid(True, alpha=0.25)

    ax = axes[0, 1]
    ax.plot(deg, sl, "--", label="Source16 aligned")
    ax.plot(deg, hl, label="Hybrid champion")
    for x in (lmax, lmin):
        ax.axvline(x, linewidth=0.9, linestyle=":")
    ax.axvline(lkink, linewidth=1.0, linestyle="-.")
    ax.set_title("Large piston — normalized volume")
    ax.set_ylabel("q")
    ax.legend()
    ax.grid(True, alpha=0.25)

    ax = axes[1, 0]
    ax.plot(deg, dss, "--", label="Source16 aligned")
    ax.plot(deg, dhs, label="Hybrid champion")
    ax.axhline(0.0, linewidth=0.8)
    for x in (smax, smin):
        ax.axvline(x, linewidth=0.9, linestyle=":")
    ax.axvline(skink, linewidth=1.0, linestyle="-.")
    ax.set_title("Small piston — normalized velocity")
    ax.set_ylabel("dq/dt")
    ax.set_xlabel("Motor angle [deg]")
    ax.legend()
    ax.grid(True, alpha=0.25)

    ax = axes[1, 1]
    ax.plot(deg, dsl, "--", label="Source16 aligned")
    ax.plot(deg, dhl, label="Hybrid champion")
    ax.axhline(0.0, linewidth=0.8)
    for x in (lmax, lmin):
        ax.axvline(x, linewidth=0.9, linestyle=":")
    ax.axvline(lkink, linewidth=1.0, linestyle="-.")
    ax.set_title("Large piston — normalized velocity")
    ax.set_ylabel("dq/dt")
    ax.set_xlabel("Motor angle [deg]")
    ax.legend()
    ax.grid(True, alpha=0.25)

    for ax in axes.ravel():
        ax.set_xlim(0, 360)
        ax.set_xticks(np.arange(0, 361, 45))

    eta = 100.0 * best["result"]["indicated_thermal_efficiency"]
    power = best["result"]["indicated_power_w"]
    source_eta = 100.0 * report["source16"]["efficiency"]
    source_power = report["source16"]["power_W"]

    fig.suptitle(
        "Hybrid compact champion vs source16\n"
        f"Hybrid: η={eta:.3f}%  P={power:.3f} W   |   "
        f"Source16: η={source_eta:.3f}%  P={source_power:.3f} W\n"
        f"S kink={skink:.1f}°   L kink={lkink:.1f}°   "
        f"BP roundings: L↓={100*p['large_down_rounding']:.1f}%  "
        f"S↑={100*p['small_up_rounding']:.1f}%"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.91))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    print(f"Saved {args.output}")
    print(json.dumps({
        "small_max_deg": smax,
        "small_min_deg": smin,
        "small_kink_deg": skink,
        "large_max_deg": lmax,
        "large_min_deg": lmin,
        "large_kink_deg": lkink,
    }, indent=2))


if __name__ == "__main__":
    main()
