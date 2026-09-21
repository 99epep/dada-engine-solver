#!/usr/bin/env python3
"""Plot the best 6-harmonic smooth motion against the retained linear UU motion.

Run:
    PYTHONPATH=src python3 examples/plot_motor_fourier_c2_6h_refine.py
"""
from __future__ import annotations

import json
import numpy as np
import matplotlib.pyplot as plt

from compare_motor_motion_laws_stage7A5 import ROOT
from dada_solver.geometry import CylinderVolumeLimits
from optimize_motor_fourier_c2_260k import FourierVolumeKinematics


REPORT = ROOT / "outputs" / "motor_fourier_c2_6h_refine" / "report.json"
OUTPUT_PNG = ROOT / "outputs" / "motor_fourier_c2_6h_refine" / "comparison_vs_linear.png"
OUTPUT_CSV = ROOT / "outputs" / "motor_fourier_c2_6h_refine" / "comparison_vs_linear.csv"


def linear_fraction(t, p, piston):
    knots = np.asarray([0.0, p["t1"], p["t2"], p["t3"], 1.0])
    if piston == "small":
        values = np.asarray([p["b_s"], 1.0, p["a_s"], 0.0, p["b_s"]])
    else:
        values = np.asarray([1.0, p["b_l"], 0.0, p["a_l"], 1.0])
    return np.interp(np.asarray(t), knots, values)


def main():
    if not REPORT.exists():
        raise FileNotFoundError(REPORT)

    report = json.loads(REPORT.read_text())
    best = report.get("best_feasible")
    if best is None:
        raise RuntimeError("No best_feasible candidate in report.")
    ref = report.get("reference_linear_UU")
    if ref is None or ref.get("parameters") is None:
        raise RuntimeError("Missing reference_linear_UU in report.")

    p = ref["parameters"]

    lim = CylinderVolumeLimits(minimum=1e-6, maximum=2e-6)
    kin = FourierVolumeKinematics(
        lim,
        lim,
        tuple(best["small_coefficients"]),
        tuple(best["large_coefficients"]),
        harmonics=best["harmonics"],
    )

    t = np.linspace(0.0, 1.0, 2401)
    angle = 360.0 * t

    s6, l6 = kin.normalized_fractions(t)
    ds6, dl6 = kin.normalized_fraction_derivatives(t)
    sl = linear_fraction(t, p, "small")
    ll = linear_fraction(t, p, "large")

    delta_s = s6 - sl
    delta_l = l6 - ll

    arr = np.column_stack(
        (t, angle, sl, ll, s6, l6, delta_s, delta_l, ds6, dl6)
    )
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        OUTPUT_CSV,
        arr,
        delimiter=",",
        header=(
            "motor_time_fraction,motor_angle_deg,"
            "small_linear,large_linear,small_6h,large_6h,"
            "delta_small,delta_large,dsmall6h_dt,dlarge6h_dt"
        ),
        comments="",
    )

    fig = plt.figure(figsize=(11, 10))
    ax1 = fig.add_axes([0.09, 0.71, 0.86, 0.22])
    ax2 = fig.add_axes([0.09, 0.40, 0.86, 0.20])
    ax3 = fig.add_axes([0.09, 0.10, 0.86, 0.20])

    ax1.plot(angle, sl, "--", label="Small — linear UU")
    ax1.plot(angle, ll, "--", label="Large — linear UU")
    ax1.plot(angle, s6, label="Small — best 6H")
    ax1.plot(angle, l6, label="Large — best 6H")
    ax1.set_xlim(0, 360)
    ax1.set_ylim(-0.03, 1.03)
    ax1.set_xticks(np.arange(0, 361, 45))
    ax1.set_ylabel("Swept-volume fraction")
    ax1.grid(True, alpha=0.25)
    ax1.legend(ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.40))

    ax2.plot(angle, delta_s, label="6H - linear (small)")
    ax2.plot(angle, delta_l, label="6H - linear (large)")
    ax2.axhline(0.0, linewidth=0.8)
    ax2.set_xlim(0, 360)
    ax2.set_xticks(np.arange(0, 361, 45))
    ax2.set_ylabel("Fraction difference")
    ax2.grid(True, alpha=0.25)
    ax2.legend(ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.42))

    ax3.plot(angle, ds6, label="dS/dt - best 6H")
    ax3.plot(angle, dl6, label="dL/dt - best 6H")
    ax3.axhline(0.0, linewidth=0.8)
    ax3.set_xlim(0, 360)
    ax3.set_xticks(np.arange(0, 361, 45))
    ax3.set_xlabel("Forward motor angle (deg)")
    ax3.set_ylabel("Fraction / cycle")
    ax3.grid(True, alpha=0.25)
    ax3.legend(ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.42))

    eta_ref = ref["efficiency"]
    eta_best = best["result"]["indicated_thermal_efficiency"]
    power_ref = ref["power_W"]
    power_best = best["result"]["indicated_power_w"]

    fig.suptitle(
        "260 K UU - retained linear motion vs best 6-harmonic smooth motion\n"
        f"eta: {100*eta_ref:.4f}% -> {100*eta_best:.4f}%    "
        f"power: {power_ref:.2f} W -> {power_best:.2f} W"
    )

    fig.savefig(OUTPUT_PNG, dpi=160)
    print(OUTPUT_PNG)
    print(OUTPUT_CSV)


if __name__ == "__main__":
    main()
