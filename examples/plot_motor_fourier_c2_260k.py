#!/usr/bin/env python3
"""Plot the retained 260 K linear UU motion against the best smooth candidate.

Run:
    PYTHONPATH=src python3 examples/plot_motor_fourier_c2_260k.py
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from compare_motor_motion_laws_stage7A5 import ROOT
from optimize_motor_fourier_c2_260k import (
    FourierVolumeKinematics,
    SOURCE_REPORT,
    _linear_fraction,
)
from dada_solver.geometry import CylinderVolumeLimits


REPORT = ROOT / "outputs" / "motor_fourier_c2_260k" / "report.json"
OUTPUT = ROOT / "outputs" / "motor_fourier_c2_260k" / "comparison.png"


def main():
    if not REPORT.exists():
        raise FileNotFoundError(REPORT)

    study = json.loads(REPORT.read_text())
    best = study.get("best_feasible")
    if best is None:
        raise RuntimeError("No feasible smooth candidate yet.")

    source = json.loads(SOURCE_REPORT.read_text())
    champion = source.get("champion")
    if champion is None:
        raise RuntimeError("UU source report has no champion.")
    p = champion["parameters"]

    lim = CylinderVolumeLimits(minimum=1e-6, maximum=2e-6)
    kin = FourierVolumeKinematics(
        lim,
        lim,
        tuple(best["small_coefficients"]),
        tuple(best["large_coefficients"]),
    )

    t = np.linspace(0.0, 1.0, 2001)
    angle = 360.0 * t
    s, l = kin.normalized_fractions(t)
    sd, ld = kin.normalized_fraction_derivatives(t)

    sl = _linear_fraction(t, p, "small")
    ll = _linear_fraction(t, p, "large")

    fig = plt.figure(figsize=(11, 8))
    ax1 = fig.add_axes([0.09, 0.56, 0.86, 0.37])
    ax2 = fig.add_axes([0.09, 0.10, 0.86, 0.34])

    ax1.plot(angle, sl, "--", label="Small — retained linear UU")
    ax1.plot(angle, ll, "--", label="Large — retained linear UU")
    ax1.plot(angle, s, label="Small — best smooth")
    ax1.plot(angle, l, label="Large — best smooth")
    ax1.set_xlim(0, 360)
    ax1.set_ylim(-0.03, 1.03)
    ax1.set_xticks(np.arange(0, 361, 45))
    ax1.set_ylabel("Swept-volume fraction")
    ax1.grid(True, alpha=0.25)
    ax1.legend(ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.32))

    ax2.plot(angle, sd, label="dS/dt — best smooth")
    ax2.plot(angle, ld, label="dL/dt — best smooth")
    ax2.axhline(0.0, linewidth=0.8)
    ax2.set_xlim(0, 360)
    ax2.set_xticks(np.arange(0, 361, 45))
    ax2.set_xlabel("Forward motor angle (deg)")
    ax2.set_ylabel("Fraction / cycle")
    ax2.grid(True, alpha=0.25)
    ax2.legend(ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.31))

    eta_ref = study["reference_linear_UU"].get("efficiency")
    eta_best = best["result"].get("indicated_thermal_efficiency")
    power = best["result"].get("indicated_power_w")
    fig.suptitle(
        f"260 K UU: linear reference {100*eta_ref:.4f}%  |  "
        f"smooth best {100*eta_best:.4f}%  |  {power:.2f} W"
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=160)
    print(OUTPUT)


if __name__ == "__main__":
    main()
