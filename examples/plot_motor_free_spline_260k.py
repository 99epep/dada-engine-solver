#!/usr/bin/env python3
"""Plot Fourier-8H versus best free-spline motion from a campaign report.

Works with either:
  outputs/motor_free_spline_260k_v3/report.json
or:
  outputs/motor_free_spline_260k_refine/report.json

The corresponding best_motion.csv must be in the same directory.
"""

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
        default=Path("outputs/motor_free_spline_260k_v3/report.json"),
    )
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()

    report = json.loads(args.report.read_text())
    best = report["best_feasible"]
    if best is None:
        raise RuntimeError("No best_feasible in report.")

    csv_path = args.report.parent / "best_motion.csv"
    data = np.genfromtxt(csv_path, delimiter=",", names=True)

    deg = data["motor_angle_deg"]
    fs = data["source_fourier_small"]
    fl = data["source_fourier_large"]
    ss = data["spline_small"]
    sl = data["spline_large"]
    dfs = data["source_fourier_dsmall_dt"]
    dfl = data["source_fourier_dlarge_dt"]
    dss = data["spline_dsmall_dt"]
    dsl = data["spline_dlarge_dt"]

    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

    ax = axes[0]
    ax.plot(deg, fs, "--", label="S — Fourier 8H")
    ax.plot(deg, fl, "--", label="L — Fourier 8H")
    ax.plot(deg, ss, label="S — spline")
    ax.plot(deg, sl, label="L — spline")
    ax.set_ylabel("Normalized volume")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=2)

    ax = axes[1]
    ax.plot(deg, ss - fs, label="S spline − 8H")
    ax.plot(deg, sl - fl, label="L spline − 8H")
    ax.axhline(0.0, linewidth=0.8)
    ax.set_ylabel("Δ normalized volume")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=2)

    ax = axes[2]
    ax.plot(deg, dfs, "--", label="S — Fourier 8H")
    ax.plot(deg, dfl, "--", label="L — Fourier 8H")
    ax.plot(deg, dss, label="S — spline")
    ax.plot(deg, dsl, label="L — spline")
    ax.axhline(0.0, linewidth=0.8)
    ax.set_ylabel("dq/dt")
    ax.set_xlabel("Motor angle [deg]")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=2)

    for ax in axes:
        ax.set_xlim(0, 360)
        ax.set_xticks(np.arange(0, 361, 45))

    source = report.get("source_fourier", {})
    eta8 = source.get("efficiency")
    p8 = source.get("power_W")
    eta = best["result"]["indicated_thermal_efficiency"]
    power = best["result"]["indicated_power_w"]

    title = f"Free spline versus Fourier 8H — spline η={100*eta:.5f}%"
    if eta8 is not None:
        title += f", 8H η={100*eta8:.5f}%"
    if p8 is not None:
        title += f" — P={power:.3f} W vs {p8:.3f} W"
    fig.suptitle(title)

    fig.tight_layout(rect=(0, 0, 1, 0.965))

    output = args.output
    if output is None:
        output = args.report.parent / "comparison_vs_fourier.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    print(f"Saved {output}")


if __name__ == "__main__":
    main()
