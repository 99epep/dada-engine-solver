#!/usr/bin/env python3
"""Plot exact linear, independently discovered spline, and established spline."""
from pathlib import Path
import argparse, json
import matplotlib.pyplot as plt
import numpy as np

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--report",type=Path,
        default=Path("outputs/motor_spline_from_linear_260k/report.json"))
    ap.add_argument("--output",type=Path,default=None)
    args=ap.parse_args()
    report=json.loads(args.report.read_text())
    best=report.get("best_feasible")
    if best is None: raise RuntimeError("No best_feasible in report.")
    data=np.genfromtxt(args.report.parent/"best_motion.csv",delimiter=",",names=True)
    deg=data["motor_angle_deg"]
    fig,axes=plt.subplots(3,1,figsize=(12,10),sharex=True)

    ax=axes[0]
    ax.plot(deg,data["linear_small"],"--",label="S — exact linear")
    ax.plot(deg,data["linear_large"],"--",label="L — exact linear")
    ax.plot(deg,data["candidate_small"],label="S — spline from linear")
    ax.plot(deg,data["candidate_large"],label="L — spline from linear")
    if "reference_spline_small" in data.dtype.names:
        ax.plot(deg,data["reference_spline_small"],":",label="S — established spline")
        ax.plot(deg,data["reference_spline_large"],":",label="L — established spline")
    ax.set_ylabel("Normalized volume");ax.grid(True,alpha=.25);ax.legend(ncol=3)

    ax=axes[1]
    ax.plot(deg,data["candidate_small"]-data["linear_small"],label="S candidate − linear")
    ax.plot(deg,data["candidate_large"]-data["linear_large"],label="L candidate − linear")
    if "reference_spline_small" in data.dtype.names:
        ax.plot(deg,data["candidate_small"]-data["reference_spline_small"],":",
                label="S candidate − established spline")
        ax.plot(deg,data["candidate_large"]-data["reference_spline_large"],":",
                label="L candidate − established spline")
    ax.axhline(0,lw=.8);ax.set_ylabel("Δ normalized volume")
    ax.grid(True,alpha=.25);ax.legend(ncol=2)

    ax=axes[2]
    ax.plot(deg,data["linear_dsmall_dt"],"--",label="S — exact linear")
    ax.plot(deg,data["linear_dlarge_dt"],"--",label="L — exact linear")
    ax.plot(deg,data["candidate_dsmall_dt"],label="S — spline from linear")
    ax.plot(deg,data["candidate_dlarge_dt"],label="L — spline from linear")
    if "reference_spline_dsmall_dt" in data.dtype.names:
        ax.plot(deg,data["reference_spline_dsmall_dt"],":",label="S — established spline")
        ax.plot(deg,data["reference_spline_dlarge_dt"],":",label="L — established spline")
    ax.axhline(0,lw=.8);ax.set_ylabel("dq/dt");ax.set_xlabel("Motor angle [deg]")
    ax.grid(True,alpha=.25);ax.legend(ncol=3)

    for ax in axes:
        ax.set_xlim(0,360);ax.set_xticks(np.arange(0,361,45))

    eta=best["result"]["indicated_thermal_efficiency"]
    power=best["result"]["indicated_power_w"]
    known=report.get("known_refined_spline_diagnostic_only")
    title=f"Spline independently seeded from linear — η={100*eta:.5f}%, P={power:.3f} W"
    if known:title+=f" — established spline η={100*known['efficiency']:.5f}%"
    fig.suptitle(title);fig.tight_layout(rect=(0,0,1,.965))
    out=args.output or args.report.parent/"comparison_linear_to_splines.png"
    out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out,dpi=180)
    print(f"Saved {out}")

if __name__=="__main__":
    main()
