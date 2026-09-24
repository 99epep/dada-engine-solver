#!/usr/bin/env python3
"""Visualize selected LARGE second-dyad candidates for 3952.

Reads:
  - outputs/large_second_dyad_v41_panel_3952.json
  - outputs/large_primary_cadence_3952_v4_1_rescore.json
  - outputs/motor_hybrid_c2_15p_260k/candidate_3952_motion_target.csv

For each requested primary rank, uses the family's best downstream candidate and writes:
  - a summary plot (target vs candidate position/velocity/error)
  - a mechanism strip with several snapshots over one cycle
  - optionally a GIF animation for one chosen rank

Typical uses:
  PYTHONPATH=src python3 examples/visualize_large_second_dyad_v41_panel_3952.py \
    --ranks 1,4,12,50

  PYTHONPATH=src python3 examples/visualize_large_second_dyad_v41_panel_3952.py \
    --ranks 1 --animate-rank 1
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

import search_large_sixbar_stageL1 as stage
import search_motor_hybrid_3952_sixbar as motor3952

ROOT = Path.cwd()
DOWNSTREAM = ROOT / "outputs" / "large_second_dyad_v41_panel_3952.json"
PRIMARY_LIB = ROOT / "outputs" / "large_primary_cadence_3952_v4_1_rescore.json"
TARGET = ROOT / "outputs" / "motor_hybrid_c2_15p_260k" / "candidate_3952_motion_target.csv"
OUTDIR = ROOT / "outputs" / "large_second_dyad_v41_visuals"

DOWNSTREAM_NAMES = [
    "second_pivot_x",
    "second_pivot_y",
    "link_EF",
    "link_GF",
    "H_along_over_EF",
    "H_normal_over_EF",
    "piston_rod",
    "slider_axis_offset",
    "slider_axis_angle",
]


def load_json_relaxed(path: Path):
    return json.loads(path.read_text(encoding="utf-8").replace("NaN", "null"))


def bridge_mask(y: np.ndarray, mask: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    if not np.any(mask):
        return y.copy()
    n = len(y)
    starts = [i for i in range(n) if mask[i] and not mask[(i - 1) % n]]
    ends = [i for i in range(n) if mask[i] and not mask[(i + 1) % n]]
    if len(starts) != 1 or len(ends) != 1:
        raise ValueError("expected one contiguous cyclic band")
    start, end = starts[0], ends[0]
    left = (start - 1) % n
    right = (end + 1) % n
    total = (right - left) % n
    out = y.copy()
    for k in range(1, total):
        i = (left + k) % n
        f = k / total
        out[i] = (1.0 - f) * y[left] + f * y[right]
    return out


def primary_vector(item: dict) -> np.ndarray:
    p = item["primary"]
    return np.array(
        [
            p["ground"],
            p["coupler"],
            p["rocker"],
            p["E_along"],
            p["E_normal"],
            p["phase"],
        ],
        dtype=float,
    )


def reconstruct(primary_item: dict, best: dict, theta: np.ndarray):
    pv = primary_vector(primary_item)
    pbranch = int(primary_item["branch"])
    E, Ed, psine, pdata = stage.primary(theta, pv, pbranch)

    params = best["parameters"]
    gx = float(params["second_pivot_x"])
    gy = float(params["second_pivot_y"])
    lef = float(params["link_EF"])
    lgf = float(params["link_GF"])
    ha = float(params["H_along_over_EF"])
    hn = float(params["H_normal_over_EF"])
    rod = float(params["piston_rod"])
    axis_offset = float(params["slider_axis_offset"])
    axis_angle = float(params["slider_axis_angle"])
    second_branch = int(best["second_branch"])

    G = np.array((gx, gy))
    delta = G - E
    dist = np.linalg.norm(delta, axis=1)
    u0 = delta / dist[:, None]
    along = (lef * lef - lgf * lgf + dist * dist) / (2.0 * dist)
    h = np.sqrt(lef * lef - along * along)
    n0 = np.column_stack((-u0[:, 1], u0[:, 0]))
    F = E + along[:, None] * u0 + second_branch * h[:, None] * n0

    ef = F - E
    gf = F - G
    det = ef[:, 0] * gf[:, 1] - ef[:, 1] * gf[:, 0]
    rhs = np.sum(ef * Ed, axis=1)
    Fd = np.column_stack((rhs * gf[:, 1] / det, -rhs * gf[:, 0] / det))

    u = ef / lef
    ud = (Fd - Ed) / lef
    nvec = np.column_stack((-u[:, 1], u[:, 0]))
    nd = np.column_stack((-ud[:, 1], ud[:, 0]))
    H = E + ha * lef * u + hn * lef * nvec
    Hd = Ed + ha * lef * ud + hn * lef * nd

    axis = np.array((math.cos(axis_angle), math.sin(axis_angle)))
    normal = np.array((-axis[1], axis[0]))
    origin = axis_offset * normal

    rel = H - origin
    longitudinal = rel @ axis
    transverse = rel @ normal
    longitudinal_d = Hd @ axis
    transverse_d = Hd @ normal
    margin = np.sqrt(rod * rod - transverse * transverse)
    slider = longitudinal + margin
    slider_d = longitudinal_d - transverse * transverse_d / margin
    stroke = float(np.ptp(slider))
    q = 1.0 - (slider - float(np.min(slider))) / stroke
    dq = -slider_d / stroke

    return {
        "pv": pv,
        "primary_branch": pbranch,
        "primary_sine": psine,
        "primary_data": pdata,
        "G": G,
        "E": E,
        "F": F,
        "H": H,
        "axis": axis,
        "normal": normal,
        "origin": origin,
        "rod": rod,
        "stroke": stroke,
        "slider": slider,
        "q": q,
        "dq": dq,
    }


def save_summary_plot(rank: int, best: dict, rec: dict, theta_deg: np.ndarray,
                      target_q: np.ndarray, target_dq: np.ndarray,
                      target_dq_bridge: np.ndarray, hp_mask: np.ndarray,
                      outdir: Path):
    q = rec["q"]
    dq = rec["dq"]
    err = q - target_q

    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

    axes[0].plot(theta_deg, target_q, label="target")
    axes[0].plot(theta_deg, q, label="candidate")
    axes[0].set_ylabel("q")
    axes[0].legend()

    axes[1].plot(theta_deg, err)
    axes[1].set_ylabel("q error")

    axes[2].plot(theta_deg, target_dq, label="target raw")
    axes[2].plot(theta_deg, target_dq_bridge, label="target bridged")
    axes[2].plot(theta_deg, dq, label="candidate")
    axes[2].set_ylabel("dq/dθ")
    axes[2].legend()

    axes[3].plot(theta_deg, np.where(hp_mask, 1.0, 0.0))
    axes[3].set_ylabel("HP band")
    axes[3].set_xlabel("crank angle (deg)")

    for ax in axes:
        ax.grid(True)

    fig.suptitle(
        f"Large second dyad, primary rank {rank} | full={100*best['position_rms_full_cycle']:.3f}% "
        f"HP={100*best['position_rms_hp_band']:.3f}% nonHP={100*best['position_rms_non_hp']:.3f}% "
        f"zeros={best['zero_crossing_count']} rod/stroke={best['piston_rod_over_stroke']:.3f}"
    )
    fig.tight_layout()
    path = outdir / f"rank_{rank:02d}_summary.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def save_mechanism_strip(rank: int, best: dict, rec: dict, theta_deg: np.ndarray,
                         outdir: Path, nframes: int = 8):
    E = rec["E"]
    F = rec["F"]
    H = rec["H"]
    G = rec["G"]
    origin = rec["origin"]
    axis = rec["axis"]
    rod = rec["rod"]
    slider = rec["slider"]

    idx = np.linspace(0, len(theta_deg) - 1, nframes, dtype=int)
    pts = np.vstack([E, F, H, G[None, :]])
    xmin, ymin = np.min(pts, axis=0) - 1.5
    xmax, ymax = np.max(pts, axis=0) + 1.5

    fig, axes = plt.subplots(2, math.ceil(nframes / 2), figsize=(14, 7))
    axes = np.array(axes).reshape(-1)

    for ax, i in zip(axes, idx):
        A = np.array((0.0, 0.0))
        D = np.array((1.0, 0.0))
        Ei = E[i]
        Fi = F[i]
        Hi = H[i]
        si = slider[i]
        slider_pt = origin + si * axis

        ax.plot([A[0], Ei[0]], [A[1], Ei[1]])
        ax.plot([Ei[0], D[0]], [Ei[1], D[1]])
        ax.plot([D[0], A[0]], [D[1], A[1]])
        ax.plot([Ei[0], Fi[0]], [Ei[1], Fi[1]])
        ax.plot([G[0], Fi[0]], [G[1], Fi[1]])
        ax.plot([Ei[0], Hi[0]], [Ei[1], Hi[1]])
        ax.plot([Hi[0], slider_pt[0]], [Hi[1], slider_pt[1]])
        ax.plot([G[0]], [G[1]], marker='o')
        ax.plot([A[0], D[0], Ei[0], Fi[0], Hi[0], slider_pt[0]],
                [A[1], D[1], Ei[1], Fi[1], Hi[1], slider_pt[1]],
                marker='o', linestyle='None')
        axis_line = np.vstack([origin - 2 * axis, origin + (max(slider) + 2) * axis])
        ax.plot(axis_line[:, 0], axis_line[:, 1], linestyle='--')
        ax.set_title(f"{theta_deg[i]:.1f}°")
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.grid(True)

    for ax in axes[len(idx):]:
        ax.axis('off')

    fig.suptitle(
        f"Mechanism snapshots, primary rank {rank} | second branch {best['second_branch']:+d}"
    )
    fig.tight_layout()
    path = outdir / f"rank_{rank:02d}_mechanism_strip.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def save_animation(rank: int, best: dict, rec: dict, theta_deg: np.ndarray,
                   outdir: Path):
    E = rec["E"]
    F = rec["F"]
    H = rec["H"]
    G = rec["G"]
    origin = rec["origin"]
    axis = rec["axis"]
    slider = rec["slider"]

    pts = np.vstack([E, F, H, G[None, :]])
    xmin, ymin = np.min(pts, axis=0) - 1.5
    xmax, ymax = np.max(pts, axis=0) + 1.5

    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.grid(True)

    crank, = ax.plot([], [])
    coupler, = ax.plot([], [])
    rocker, = ax.plot([], [])
    ef_link, = ax.plot([], [])
    gf_link, = ax.plot([], [])
    eh_link, = ax.plot([], [])
    rod_link, = ax.plot([], [])
    points, = ax.plot([], [], marker='o', linestyle='None')
    axis_line, = ax.plot([], [], linestyle='--')
    title = ax.set_title("")

    A = np.array((0.0, 0.0))
    D = np.array((1.0, 0.0))
    slider_axis_line = np.vstack([origin - 2 * axis, origin + (max(slider) + 2) * axis])
    axis_line.set_data(slider_axis_line[:, 0], slider_axis_line[:, 1])

    def update(i):
        Ei = E[i]
        Fi = F[i]
        Hi = H[i]
        slider_pt = origin + slider[i] * axis
        crank.set_data([A[0], Ei[0]], [A[1], Ei[1]])
        coupler.set_data([Ei[0], D[0]], [Ei[1], D[1]])
        rocker.set_data([D[0], A[0]], [D[1], A[1]])
        ef_link.set_data([Ei[0], Fi[0]], [Ei[1], Fi[1]])
        gf_link.set_data([G[0], Fi[0]], [G[1], Fi[1]])
        eh_link.set_data([Ei[0], Hi[0]], [Ei[1], Hi[1]])
        rod_link.set_data([Hi[0], slider_pt[0]], [Hi[1], slider_pt[1]])
        points.set_data(
            [A[0], D[0], G[0], Ei[0], Fi[0], Hi[0], slider_pt[0]],
            [A[1], D[1], G[1], Ei[1], Fi[1], Hi[1], slider_pt[1]],
        )
        title.set_text(f"rank {rank} | θ={theta_deg[i]:.1f}°")
        return crank, coupler, rocker, ef_link, gf_link, eh_link, rod_link, points, axis_line, title

    ani = animation.FuncAnimation(fig, update, frames=len(theta_deg), interval=30, blit=False)
    path = outdir / f"rank_{rank:02d}_animation.gif"
    writer = animation.PillowWriter(fps=25)
    ani.save(path, writer=writer)
    plt.close(fig)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--downstream", type=Path, default=DOWNSTREAM)
    ap.add_argument("--primary-library", type=Path, default=PRIMARY_LIB)
    ap.add_argument("--target", type=Path, default=TARGET)
    ap.add_argument("--outdir", type=Path, default=OUTDIR)
    ap.add_argument("--ranks", default="1,4,12,50")
    ap.add_argument("--animate-rank", type=int, default=None)
    args = ap.parse_args()

    data = load_json_relaxed(args.downstream)
    prim = load_json_relaxed(args.primary_library)
    primary_by_rank = {int(x["library_rank"]): x for x in prim.get("candidates", [])}
    families = {int(f["primary_rank"]): f for f in data.get("families", [])}

    raw = np.genfromtxt(args.target, delimiter=",", names=True)
    if raw["theta_rad"][-1] >= 2.0 * math.pi - 1e-9:
        raw = raw[:-1]
    theta = np.asarray(raw["theta_rad"], dtype=float)
    theta_deg = np.mod(np.degrees(theta), 360.0)
    target_q = np.asarray(raw["large_fraction_0_1"], dtype=float)
    target_dq = np.asarray(raw["large_dq_dtheta_per_rad"], dtype=float)

    band = motor3952._detect_noise_bands(args.target, 6)["large"]
    hp_mask = motor3952._theta_band_mask(theta, band)
    target_dq_bridge = bridge_mask(target_dq, hp_mask)

    ranks = [int(x) for x in args.ranks.split(",") if x.strip()]
    args.outdir.mkdir(parents=True, exist_ok=True)

    written = []
    for rank in ranks:
        fam = families.get(rank)
        if fam is None or not fam.get("best"):
            print(f"rank {rank}: no best candidate available")
            continue
        primary_item = primary_by_rank.get(rank)
        if primary_item is None:
            print(f"rank {rank}: primary not found in library")
            continue

        best = fam["best"]
        rec = reconstruct(primary_item, best, theta)

        p1 = save_summary_plot(rank, best, rec, theta_deg, target_q, target_dq,
                               target_dq_bridge, hp_mask, args.outdir)
        p2 = save_mechanism_strip(rank, best, rec, theta_deg, args.outdir)
        written.extend([p1, p2])
        print(f"rank {rank}: wrote {p1.name}, {p2.name}")

        if args.animate_rank is not None and rank == args.animate_rank:
            p3 = save_animation(rank, best, rec, theta_deg, args.outdir)
            written.append(p3)
            print(f"rank {rank}: wrote {p3.name}")

    print("\nOutput directory:", args.outdir)
    for path in written:
        print(path)


if __name__ == "__main__":
    main()
