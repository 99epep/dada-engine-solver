#!/usr/bin/env python3
"""Plot the 7-mobile-knot spline result and overlay the knot positions."""

from pathlib import Path
import argparse
import json

import matplotlib.pyplot as plt
import numpy as np


def _sorted_periodic(x_deg, y):
    x = np.asarray(x_deg, dtype=float) % 360.0
    y = np.asarray(y, dtype=float)
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    return x, y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--report",
        type=Path,
        default=Path("outputs/motor_mobile_spline_7k_260k/report.json"),
    )
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()

    report = json.loads(args.report.read_text())
    best = report.get("best_feasible")
    if best is None:
        raise RuntimeError("No best_feasible in report.")

    csv_path = args.report.parent / "best_motion.csv"
    data = np.genfromtxt(csv_path, delimiter=",", names=True)

    deg = data["motor_angle_deg"]
    src_s = data["source16_small"]
    src_l = data["source16_large"]
    mob_s = data["mobile7_small"]
    mob_l = data["mobile7_large"]
    dsrc_s = data["source16_dsmall_dt"]
    dsrc_l = data["source16_dlarge_dt"]
    dmob_s = data["mobile7_dsmall_dt"]
    dmob_l = data["mobile7_dlarge_dt"]

    kdiag = best["knot_diagnostics"]
    small_knots = np.asarray(kdiag["small_motor_knot_deg"], dtype=float) % 360.0
    large_knots = np.asarray(kdiag["large_motor_knot_deg"], dtype=float) % 360.0

    small_knot_y = np.interp(small_knots, deg, mob_s)
    large_knot_y = np.interp(large_knots, deg, mob_l)
    small_knot_v = np.interp(small_knots, deg, dmob_s)
    large_knot_v = np.interp(large_knots, deg, dmob_l)

    fig, axes = plt.subplots(3, 1, figsize=(12, 11), sharex=True)

    ax = axes[0]
    ax.plot(deg, src_s, "--", label="S — source 16 controls")
    ax.plot(deg, src_l, "--", label="L — source 16 controls")
    ax.plot(deg, mob_s, label="S — mobile 7 knots")
    ax.plot(deg, mob_l, label="L — mobile 7 knots")
    ax.plot(small_knots, small_knot_y, "o", linestyle="None", label="S knots")
    ax.plot(large_knots, large_knot_y, "s", linestyle="None", label="L knots")
    ax.set_ylabel("Normalized volume")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=3)

    ax = axes[1]
    ax.plot(deg, mob_s - src_s, label="S mobile7 − source16")
    ax.plot(deg, mob_l - src_l, label="L mobile7 − source16")
    ax.axhline(0.0, linewidth=0.8)
    for x in small_knots:
        ax.axvline(x, linewidth=0.8, alpha=0.20)
    for x in large_knots:
        ax.axvline(x, linewidth=0.8, alpha=0.20)
    ax.set_ylabel("Δ normalized volume")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=2)

    ax = axes[2]
    ax.plot(deg, dsrc_s, "--", label="S — source 16 controls")
    ax.plot(deg, dsrc_l, "--", label="L — source 16 controls")
    ax.plot(deg, dmob_s, label="S — mobile 7 knots")
    ax.plot(deg, dmob_l, label="L — mobile 7 knots")
    ax.plot(small_knots, small_knot_v, "o", linestyle="None", label="S knots")
    ax.plot(large_knots, large_knot_v, "s", linestyle="None", label="L knots")
    ax.axhline(0.0, linewidth=0.8)
    ax.set_ylabel("dq/dt")
    ax.set_xlabel("Motor angle [deg]")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=3)

    for ax in axes:
        ax.set_xlim(0, 360)
        ax.set_xticks(np.arange(0, 361, 45))

    eta = best["result"]["indicated_thermal_efficiency"]
    power = best["result"]["indicated_power_w"]
    src = report.get("source16", {})
    title = f"7 mobile knots vs 16-control source — η={100*eta:.5f}%, P={power:.3f} W"
    if "efficiency" in src:
        title += f" — source η={100*src['efficiency']:.5f}%"
    fig.suptitle(title)

    sk, _ = _sorted_periodic(small_knots, small_knot_y)
    lk, _ = _sorted_periodic(large_knots, large_knot_y)
    fig.text(0.01, 0.955, "S knots [deg]: " + ", ".join(f"{x:.1f}" for x in sk), fontsize=9, va="top")
    fig.text(0.01, 0.937, "L knots [deg]: " + ", ".join(f"{x:.1f}" for x in lk), fontsize=9, va="top")

    fig.tight_layout(rect=(0, 0, 1, 0.92))

    output = args.output or (args.report.parent / "comparison_with_knots.png")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    print(f"Saved {output}")


if __name__ == "__main__":
    main()
