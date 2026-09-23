#!/usr/bin/env python3
"""Plot candidate-3952 target motion against the synthesized six-bar pair.

The shaded regions are detected automatically from the six dominant acceleration
fronts in the target; no fixed angle or width is used.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path.cwd()
DEFAULT_CSV = ROOT / "outputs" / "sixbar_3952_motion_comparison.csv"
DEFAULT_OUT = ROOT / "outputs" / "sixbar_3952_plots"


def local_extrema_mask(y: np.ndarray, maxima: bool) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    if maxima:
        return (y > np.roll(y, 1)) & (y >= np.roll(y, -1))
    return (y < np.roll(y, 1)) & (y <= np.roll(y, -1))


def detect_band(theta_deg: np.ndarray, ddq: np.ndarray, front_count: int = 6):
    jerk = np.abs(np.roll(ddq, -1) - np.roll(ddq, 1))
    peaks = np.flatnonzero(local_extrema_mask(jerk, True))
    ranked = peaks[np.argsort(jerk[peaks])[::-1]]
    if len(ranked) < front_count + 1:
        raise ValueError("not enough acceleration fronts")

    selected = np.sort(ranked[:front_count])
    isolation = float(jerk[ranked[front_count - 1]] / max(jerk[ranked[front_count]], 1e-30))

    n = len(ddq)
    gaps = np.diff(np.r_[selected, selected[0] + n])
    cut = int(np.argmax(gaps))
    first = int(selected[(cut + 1) % front_count])
    offsets = np.sort((selected - first) % n)
    last = int((first + int(offsets[-1])) % n)

    minima = local_extrema_mask(jerk, False)

    def prev_min(i):
        for k in range(1, n):
            j = (i - k) % n
            if minima[j]:
                return j
        return (i - 1) % n

    def next_min(i):
        for k in range(1, n):
            j = (i + k) % n
            if minima[j]:
                return j
        return (i + 1) % n

    start = prev_min(first)
    end = next_min(last)
    idx = np.arange(n)
    mask = ((idx - start) % n) <= ((end - start) % n)
    fronts = [int((first + int(x)) % n) for x in offsets]
    return mask, {
        "start": float(theta_deg[start]),
        "end": float(theta_deg[end]),
        "fronts": [float(theta_deg[i]) for i in fronts],
        "isolation": isolation,
    }


def shade_mask(ax, x, mask):
    padded = np.r_[False, mask, False]
    changes = np.flatnonzero(padded[1:] != padded[:-1])
    for start, stop in changes.reshape(-1, 2):
        ax.axvspan(x[start], x[stop - 1], alpha=0.14)


def overlay(x, target, mechanism, ylabel, title, path, mask):
    fig, ax = plt.subplots(figsize=(11.5, 5.2))
    ax.plot(x, target, label="3952 target")
    ax.plot(x, mechanism, label="six-bar")
    shade_mask(ax, x, mask)
    ax.set_xlim(0, 360)
    ax.set_xlabel("Crank angle θ (deg)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def error_plot(x, error, title, path, mask):
    fig, ax = plt.subplots(figsize=(11.5, 4.6))
    ax.plot(x, error)
    shade_mask(ax, x, mask)
    ax.axhline(0, linewidth=0.8)
    ax.set_xlim(0, 360)
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
    ap.add_argument("--front-count", type=int, default=6)
    args = ap.parse_args()

    raw = np.genfromtxt(args.csv, delimiter=",", names=True)
    if raw["theta_deg"][-1] >= 360.0 - 1e-9:
        raw = raw[:-1]
    x = np.asarray(raw["theta_deg"], dtype=float)

    sides = {
        "small": {
            "q_t": "small_target_q", "q_m": "small_mechanism_q",
            "dq_t": "small_target_dq", "dq_m": "small_mechanism_dq",
            "ddq_t": "small_target_ddq", "ddq_m": "small_mechanism_ddq",
        },
        "large": {
            "q_t": "large_target_q", "q_m": "large_mechanism_q",
            "dq_t": "large_target_dq", "dq_m": "large_mechanism_dq",
            "ddq_t": "large_target_ddq", "ddq_m": "large_mechanism_ddq",
        },
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)

    for side, c in sides.items():
        qt = np.asarray(raw[c["q_t"]], dtype=float)
        qm = np.asarray(raw[c["q_m"]], dtype=float)
        dqt = np.asarray(raw[c["dq_t"]], dtype=float)
        dqm = np.asarray(raw[c["dq_m"]], dtype=float)
        ddqt = np.asarray(raw[c["ddq_t"]], dtype=float)
        ddqm = np.asarray(raw[c["ddq_m"]], dtype=float)

        excluded, band = detect_band(x, ddqt, args.front_count)
        fit = ~excluded
        full_rms = float(np.sqrt(np.mean((qm - qt) ** 2)))
        fit_rms = float(np.sqrt(np.mean((qm[fit] - qt[fit]) ** 2)))

        title = side.capitalize()
        overlay(
            x, qt, qm, "Normalized piston position q",
            f"{title} position | full RMS={100*full_rms:.3f}% | "
            f"RMS outside noise={100*fit_rms:.3f}%",
            args.output_dir / f"{side}_position.png", excluded,
        )
        error_plot(
            x, qm - qt, f"{title} position error",
            args.output_dir / f"{side}_position_error.png", excluded,
        )
        overlay(
            x, dqt, dqm, "dq/dθ (rad⁻¹)", f"{title} velocity",
            args.output_dir / f"{side}_velocity.png", excluded,
        )
        overlay(
            x, ddqt, ddqm, "d²q/dθ² (rad⁻²)", f"{title} acceleration",
            args.output_dir / f"{side}_acceleration.png", excluded,
        )

        print(
            f"{side}: excluded {100*np.mean(excluded):.2f}% of cycle; "
            f"fronts={','.join(f'{a:.3f}' for a in band['fronts'])}; "
            f"isolation={band['isolation']:.2f}; "
            f"full RMS={100*full_rms:.6f}%; fit RMS={100*fit_rms:.6f}%"
        )

    print(f"Wrote plots to {args.output_dir}")


if __name__ == "__main__":
    main()
