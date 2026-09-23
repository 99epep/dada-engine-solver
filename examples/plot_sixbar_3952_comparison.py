#!/usr/bin/env python3
"""Plot candidate-3952 target motion against the synthesized six-bar pair."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path.cwd()
DEFAULT_CSV = ROOT / "outputs" / "sixbar_3952_motion_comparison.csv"
DEFAULT_OUT = ROOT / "outputs" / "sixbar_3952_plots"


def cyclic_mask(x_deg: np.ndarray, center_deg: float, half_width_deg: float) -> np.ndarray:
    distance = np.abs((x_deg - center_deg + 180.0) % 360.0 - 180.0)
    return distance <= half_width_deg


def shade_mask(ax, x: np.ndarray, mask: np.ndarray) -> None:
    if not np.any(mask):
        return
    padded = np.r_[False, mask, False]
    changes = np.flatnonzero(padded[1:] != padded[:-1])
    for start, stop in changes.reshape(-1, 2):
        x0 = x[max(0, start)]
        x1 = x[min(len(x) - 1, stop - 1)]
        ax.axvspan(x0, x1, alpha=0.14)


def save_overlay(x, target, mechanism, ylabel, title, path, hp_mask):
    fig, ax = plt.subplots(figsize=(11.5, 5.2))
    ax.plot(x, target, label="3952 target")
    ax.plot(x, mechanism, label="six-bar")
    shade_mask(ax, x, hp_mask)
    ax.set_xlim(0.0, 360.0)
    ax.set_xlabel("Crank angle θ (deg)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def save_error(x, error, title, path, hp_mask):
    fig, ax = plt.subplots(figsize=(11.5, 4.6))
    ax.plot(x, error)
    shade_mask(ax, x, hp_mask)
    ax.axhline(0.0, linewidth=0.8)
    ax.set_xlim(0.0, 360.0)
    ax.set_xlabel("Crank angle θ (deg)")
    ax.set_ylabel("six-bar − target")
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--exclude-half-width-deg", type=float, default=8.0)
    ap.add_argument(
        "--small-hp-center-t-deg", type=float, default=216.73421234322282,
        help="Small HP feature center in normalized structured-motion t degrees.",
    )
    ap.add_argument(
        "--large-hp-center-t-deg", type=float, default=320.90381952231166,
        help="Large HP feature center in normalized structured-motion t degrees.",
    )
    args = ap.parse_args()

    raw = np.genfromtxt(args.csv, delimiter=",", names=True)
    x = np.asarray(raw["theta_deg"], dtype=float)
    stat_mask = x < 360.0 - 1e-9

    sides = {
        "small": {
            "center": (-args.small_hp_center_t_deg) % 360.0,
            "q_t": "small_target_q",
            "q_m": "small_mechanism_q",
            "dq_t": "small_target_dq",
            "dq_m": "small_mechanism_dq",
            "ddq_t": "small_target_ddq",
            "ddq_m": "small_mechanism_ddq",
        },
        "large": {
            "center": (-args.large_hp_center_t_deg) % 360.0,
            "q_t": "large_target_q",
            "q_m": "large_mechanism_q",
            "dq_t": "large_target_dq",
            "dq_m": "large_mechanism_dq",
            "ddq_t": "large_target_ddq",
            "ddq_m": "large_mechanism_ddq",
        },
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)

    for side, cfg in sides.items():
        hp_mask = cyclic_mask(x, cfg["center"], args.exclude_half_width_deg)
        fit_mask = stat_mask & ~hp_mask

        qt = np.asarray(raw[cfg["q_t"]], dtype=float)
        qm = np.asarray(raw[cfg["q_m"]], dtype=float)
        dqt = np.asarray(raw[cfg["dq_t"]], dtype=float)
        dqm = np.asarray(raw[cfg["dq_m"]], dtype=float)
        ddqt = np.asarray(raw[cfg["ddq_t"]], dtype=float)
        ddqm = np.asarray(raw[cfg["ddq_m"]], dtype=float)

        rms_full = float(np.sqrt(np.mean((qm[stat_mask] - qt[stat_mask]) ** 2)))
        rms_fit = float(np.sqrt(np.mean((qm[fit_mask] - qt[fit_mask]) ** 2)))

        side_title = side.capitalize()
        save_overlay(
            x, qt, qm, "Normalized piston position q",
            f"{side_title} cylinder — position | full RMS={100*rms_full:.3f}% | "
            f"RMS outside HP ±{args.exclude_half_width_deg:g}°={100*rms_fit:.3f}%",
            args.output_dir / f"{side}_position.png", hp_mask,
        )
        save_error(
            x, qm - qt,
            f"{side_title} cylinder — position error",
            args.output_dir / f"{side}_position_error.png", hp_mask,
        )
        save_overlay(
            x, dqt, dqm, "dq/dθ (rad⁻¹)",
            f"{side_title} cylinder — velocity",
            args.output_dir / f"{side}_velocity.png", hp_mask,
        )
        save_overlay(
            x, ddqt, ddqm, "d²q/dθ² (rad⁻²)",
            f"{side_title} cylinder — acceleration",
            args.output_dir / f"{side}_acceleration.png", hp_mask,
        )

        print(
            f"{side}: full position RMS={100*rms_full:.6f}% ; "
            f"outside HP ±{args.exclude_half_width_deg:g}°={100*rms_fit:.6f}%"
        )

    print(f"Wrote plots to {args.output_dir}")


if __name__ == "__main__":
    main()
