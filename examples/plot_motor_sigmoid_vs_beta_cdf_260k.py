#!/usr/bin/env python3
"""Compare the power-sigmoid and beta-CDF champions against the aligned source16 motion."""

from pathlib import Path
import argparse
import json

import matplotlib.pyplot as plt
import numpy as np


def _load_report(path: Path):
    report = json.loads(path.read_text())
    best = report.get("best_feasible")
    if best is None:
        raise RuntimeError(f"No best_feasible in {path}")
    return report, best


def _load_motion(csv_path: Path, small_name: str, large_name: str,
                 dsmall_name: str, dlarge_name: str):
    data = np.genfromtxt(csv_path, delimiter=",", names=True)
    return {
        "deg": data["motor_angle_deg"],
        "source_small": data["source16_aligned_small"],
        "source_large": data["source16_aligned_large"],
        "small": data[small_name],
        "large": data[large_name],
        "dsource_small": data["source16_aligned_dsmall_dt"],
        "dsource_large": data["source16_aligned_dlarge_dt"],
        "dsmall": data[dsmall_name],
        "dlarge": data[dlarge_name],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--sigmoid-report",
        type=Path,
        default=Path("outputs/motor_power_sigmoid_260k_fast/report.json"),
    )
    ap.add_argument(
        "--beta-report",
        type=Path,
        default=Path("outputs/motor_beta_cdf_260k/report.json"),
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/compare_sigmoid_beta_cdf_260k.png"),
    )
    args = ap.parse_args()

    s_report, s_best = _load_report(args.sigmoid_report)
    b_report, b_best = _load_report(args.beta_report)

    s_motion = _load_motion(
        args.sigmoid_report.parent / "best_motion.csv",
        "power_sigmoid_small",
        "power_sigmoid_large",
        "power_sigmoid_dsmall_dt",
        "power_sigmoid_dlarge_dt",
    )
    b_motion = _load_motion(
        args.beta_report.parent / "best_motion.csv",
        "beta_cdf_small",
        "beta_cdf_large",
        "beta_cdf_dsmall_dt",
        "beta_cdf_dlarge_dt",
    )

    deg = s_motion["deg"]

    fig, axes = plt.subplots(3, 2, figsize=(13, 12), sharex=True)

    ax = axes[0, 0]
    ax.plot(deg, s_motion["source_small"], "--", label="Source16 aligned")
    ax.plot(deg, s_motion["small"], label="Power-sigmoid champion")
    ax.plot(deg, b_motion["small"], label="Beta-CDF champion")
    ax.set_title("Small piston — normalized volume")
    ax.set_ylabel("q")
    ax.grid(True, alpha=0.25)
    ax.legend()

    ax = axes[0, 1]
    ax.plot(deg, s_motion["source_large"], "--", label="Source16 aligned")
    ax.plot(deg, s_motion["large"], label="Power-sigmoid champion")
    ax.plot(deg, b_motion["large"], label="Beta-CDF champion")
    ax.set_title("Large piston — normalized volume")
    ax.set_ylabel("q")
    ax.grid(True, alpha=0.25)
    ax.legend()

    ax = axes[1, 0]
    ax.plot(deg, s_motion["dsource_small"], "--", label="Source16 aligned")
    ax.plot(deg, s_motion["dsmall"], label="Power-sigmoid champion")
    ax.plot(deg, b_motion["dsmall"], label="Beta-CDF champion")
    ax.axhline(0.0, linewidth=0.8)
    ax.set_title("Small piston — normalized velocity")
    ax.set_ylabel("dq/dt")
    ax.grid(True, alpha=0.25)
    ax.legend()

    ax = axes[1, 1]
    ax.plot(deg, s_motion["dsource_large"], "--", label="Source16 aligned")
    ax.plot(deg, s_motion["dlarge"], label="Power-sigmoid champion")
    ax.plot(deg, b_motion["dlarge"], label="Beta-CDF champion")
    ax.axhline(0.0, linewidth=0.8)
    ax.set_title("Large piston — normalized velocity")
    ax.set_ylabel("dq/dt")
    ax.grid(True, alpha=0.25)
    ax.legend()

    ax = axes[2, 0]
    ax.plot(deg, s_motion["small"] - s_motion["source_small"],
            label="Power-sigmoid − source16")
    ax.plot(deg, b_motion["small"] - b_motion["source_small"],
            label="Beta-CDF − source16")
    ax.axhline(0.0, linewidth=0.8)
    ax.set_title("Small piston — deviation from source16")
    ax.set_ylabel("Δq")
    ax.set_xlabel("Motor angle [deg]")
    ax.grid(True, alpha=0.25)
    ax.legend()

    ax = axes[2, 1]
    ax.plot(deg, s_motion["large"] - s_motion["source_large"],
            label="Power-sigmoid − source16")
    ax.plot(deg, b_motion["large"] - b_motion["source_large"],
            label="Beta-CDF − source16")
    ax.axhline(0.0, linewidth=0.8)
    ax.set_title("Large piston — deviation from source16")
    ax.set_ylabel("Δq")
    ax.set_xlabel("Motor angle [deg]")
    ax.grid(True, alpha=0.25)
    ax.legend()

    for ax in axes.ravel():
        ax.set_xlim(0, 360)
        ax.set_xticks(np.arange(0, 361, 45))

    s_eta = 100.0 * s_best["result"]["indicated_thermal_efficiency"]
    s_pow = s_best["result"]["indicated_power_w"]
    b_eta = 100.0 * b_best["result"]["indicated_thermal_efficiency"]
    b_pow = b_best["result"]["indicated_power_w"]
    src_eta = 100.0 * s_report["source16"]["efficiency"]
    src_pow = s_report["source16"]["power_W"]

    fig.suptitle(
        "Compact-law champions vs source16\n"
        f"Source16: η={src_eta:.3f}%  P={src_pow:.3f} W   |   "
        f"Power-sigmoid: η={s_eta:.3f}%  P={s_pow:.3f} W   |   "
        f"Beta-CDF: η={b_eta:.3f}%  P={b_pow:.3f} W"
    )

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
