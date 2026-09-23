#!/usr/bin/env python3
"""Study plots for the final 15-parameter structured C2 260 K campaign.

Reads:
  outputs/motor_hybrid_c2_15p_260k/report.json
  outputs/motor_hybrid_c2_15p_260k/best_motion.csv
  outputs/motor_hybrid_c2_15p_260k/history.jsonl

Produces:
  kinematics_detail.png
  kinematics_phase_plane.png
  optimization_convergence.png
"""

from __future__ import annotations

from pathlib import Path
import argparse
import json

import matplotlib.pyplot as plt
import numpy as np


def load_history(path):
    rows = []
    if not path.exists():
        return rows
    with path.open() as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def kink_geometry(p):
    smax = p["small_max_deg"] % 360.0
    sdown = p["small_down_duration_deg"]
    smin = (smax + sdown) % 360.0

    sku = p["small_down_kink_u"]
    skw_rel = p["small_down_kink_width_rel"]
    sk_center = (smax + sku*sdown) % 360.0
    sk_full_width = skw_rel * min(sku, 1.0-sku) * sdown

    lmax = 0.0
    lmin = p["large_down_duration_deg"] % 360.0
    lup = 360.0 - p["large_down_duration_deg"]

    lku = p["large_up_kink_u"]
    lkw_rel = p["large_up_kink_width_rel"]
    lk_center = (lmin + lku*lup) % 360.0
    lk_full_width = lkw_rel * min(lku, 1.0-lku) * lup

    return {
        "small_max": smax,
        "small_min": smin,
        "large_max": lmax,
        "large_min": lmin,
        "small_kink_center": sk_center,
        "small_kink_width": sk_full_width,
        "large_kink_center": lk_center,
        "large_kink_width": lk_full_width,
    }


def shade_periodic(ax, center, width):
    lo = center - 0.5*width
    hi = center + 0.5*width
    if lo >= 0 and hi <= 360:
        ax.axvspan(lo, hi, alpha=0.12)
    elif lo < 0:
        ax.axvspan(0, hi, alpha=0.12)
        ax.axvspan(360+lo, 360, alpha=0.12)
    else:
        ax.axvspan(lo, 360, alpha=0.12)
        ax.axvspan(0, hi-360, alpha=0.12)


def plot_kinematics(report, motion, outdir):
    best = report["best_feasible"]
    p = best["parameters"]
    g = kink_geometry(p)

    deg = motion["motor_angle_deg"]

    src_q = (
        motion["source16_aligned_small"],
        motion["source16_aligned_large"],
    )
    cur_q = (
        motion["structured15_small"],
        motion["structured15_large"],
    )
    src_v = (
        motion["source16_aligned_dsmall_dt"],
        motion["source16_aligned_dlarge_dt"],
    )
    cur_v = (
        motion["structured15_dsmall_dt"],
        motion["structured15_dlarge_dt"],
    )
    src_a = (
        motion["source16_aligned_ddsmall_dt2"],
        motion["source16_aligned_ddlarge_dt2"],
    )
    cur_a = (
        motion["structured15_ddsmall_dt2"],
        motion["structured15_ddlarge_dt2"],
    )

    fig, axes = plt.subplots(3, 2, figsize=(13, 11), sharex=True)
    names = ("Small piston", "Large piston")

    extrema = (
        (g["small_max"], g["small_min"]),
        (g["large_max"], g["large_min"]),
    )
    kinks = (
        (g["small_kink_center"], g["small_kink_width"]),
        (g["large_kink_center"], g["large_kink_width"]),
    )

    for j in range(2):
        for i, (source, current, ylabel, title) in enumerate((
            (src_q[j], cur_q[j], "q", "normalized position"),
            (src_v[j], cur_v[j], "dq/dt", "normalized velocity"),
            (src_a[j], cur_a[j], "d²q/dt²", "normalized acceleration"),
        )):
            ax = axes[i, j]
            ax.plot(deg, source, "--", label="Source16 aligned")
            ax.plot(deg, current, label="Structured 15p")
            for x in extrema[j]:
                ax.axvline(x, linestyle=":", linewidth=0.9)
            shade_periodic(ax, *kinks[j])
            if i > 0:
                ax.axhline(0.0, linewidth=0.8)
            ax.set_title(f"{names[j]} — {title}")
            ax.set_ylabel(ylabel)
            ax.grid(True, alpha=0.25)
            if i == 0:
                ax.legend()

    for ax in axes.ravel():
        ax.set_xlim(0, 360)
        ax.set_xticks(np.arange(0, 361, 45))

    axes[2, 0].set_xlabel("Motor angle [deg]")
    axes[2, 1].set_xlabel("Motor angle [deg]")

    eta = 100*best["result"]["indicated_thermal_efficiency"]
    source_eta = 100*report["source16"]["efficiency"]

    fig.suptitle(
        "Structured 15p champion vs 16-control source\n"
        f"η15p={eta:.5f}%   η16c={source_eta:.5f}%   "
        f"Δ={eta-source_eta:+.5f} percentage point\n"
        f"BP mid-q: S↑={p['small_up_bp_mid_q']:.4f}, "
        f"L↓={p['large_down_bp_mid_q']:.4f}   |   "
        f"HP kink widths: S↓={p['small_down_kink_width_rel']:.4f}, "
        f"L↑={p['large_up_kink_width_rel']:.4f}"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(outdir/"kinematics_detail.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for j, ax in enumerate(axes):
        ax.plot(src_q[j], src_v[j], "--", label="Source16 aligned")
        ax.plot(cur_q[j], cur_v[j], label="Structured 15p")
        ax.axhline(0.0, linewidth=0.8)
        ax.set_title(f"{names[j]} — phase plane")
        ax.set_xlabel("Normalized position q")
        ax.set_ylabel("dq/dt")
        ax.grid(True, alpha=0.25)
        ax.legend()

    fig.suptitle(
        "Position–velocity phase plane\n"
        "Useful for seeing extremum sharpness and the HP-side change of slope"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(outdir/"kinematics_phase_plane.png", dpi=180)
    plt.close(fig)


def plot_history(report, history, outdir):
    feasible = []
    for r in history:
        res = r.get("result") or {}
        eta = res.get("indicated_thermal_efficiency")
        if (
            r.get("feasible")
            and res.get("status") == "converged"
            and eta is not None
        ):
            feasible.append((
                int(r.get("proposal_index", r.get("index", 0))),
                float(eta),
                int(r.get("radius_index", -1)),
            ))

    if not feasible:
        return

    feasible.sort()
    x = np.array([r[0] for r in feasible])
    eta = 100*np.array([r[1] for r in feasible])
    radius = np.array([r[2] for r in feasible])

    best_so_far = np.maximum.accumulate(eta)

    fig, ax = plt.subplots(figsize=(13, 6))
    ax.scatter(x, eta, s=9, alpha=0.25, label="Feasible candidates")
    ax.plot(x, best_so_far, linewidth=1.7, label="Best so far")

    source_eta = 100*report["source16"]["efficiency"]
    ax.axhline(source_eta, linestyle="--", linewidth=1.1, label="16-control champion")

    changes = np.flatnonzero(radius[1:] != radius[:-1]) + 1
    for k in changes:
        ax.axvline(x[k], linestyle=":", linewidth=0.7)

    best = report["best_feasible"]
    bx = best["proposal_index"]
    by = 100*best["result"]["indicated_thermal_efficiency"]
    ax.scatter([bx], [by], s=55, label="Final 15p champion")

    ax.set_title(
        "15p optimization convergence — "
        f"{report['attempted_total']} evaluated candidates"
    )
    ax.set_xlabel("Proposal index")
    ax.set_ylabel("Indicated thermal efficiency [%]")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir/"optimization_convergence.png", dpi=180)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--directory",
        type=Path,
        default=Path("outputs/motor_hybrid_c2_15p_260k"),
    )
    args = ap.parse_args()

    d = args.directory
    report = json.loads((d/"report.json").read_text())
    motion = np.genfromtxt(d/"best_motion.csv", delimiter=",", names=True)
    history = load_history(d/"history.jsonl")

    plot_kinematics(report, motion, d)
    plot_history(report, history, d)

    p = report["best_feasible"]["parameters"]
    g = kink_geometry(p)

    print(json.dumps({
        "efficiency_percent":
            100*report["best_feasible"]["result"]["indicated_thermal_efficiency"],
        "power_W":
            report["best_feasible"]["result"]["indicated_power_w"],
        "source16_efficiency_percent":
            100*report["source16"]["efficiency"],
        "small_max_deg": g["small_max"],
        "small_min_deg": g["small_min"],
        "large_max_deg": g["large_max"],
        "large_min_deg": g["large_min"],
        "small_HP_kink_center_deg": g["small_kink_center"],
        "small_HP_kink_full_width_deg": g["small_kink_width"],
        "large_HP_kink_center_deg": g["large_kink_center"],
        "large_HP_kink_full_width_deg": g["large_kink_width"],
        "small_BP_mid_q": p["small_up_bp_mid_q"],
        "large_BP_mid_q": p["large_down_bp_mid_q"],
    }, indent=2))

    print(f"Saved {d/'kinematics_detail.png'}")
    print(f"Saved {d/'kinematics_phase_plane.png'}")
    if history:
        print(f"Saved {d/'optimization_convergence.png'}")


if __name__ == "__main__":
    main()
