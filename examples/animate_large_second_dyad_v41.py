#!/usr/bin/env python3
"""Interactively animate a LARGE six-bar candidate from the V4.1 downstream panel.

This reuses the exact mechanism representation of
examples/animate_small_sixbar_candidate_2c.py:

- primary four-bar A-B-C-D;
- rigid primary coupler plate B-C-E;
- secondary dyad E-F-G;
- rigid secondary plate E-F-H;
- finite piston rod H-P;
- slider axis.

By default this script opens an interactive matplotlib window and does NOT
encode a GIF. GIF export is opt-in with --gif.

Typical use:

    PYTHONPATH=src python3 \
      examples/animate_large_second_dyad_v41.py \
      --rank 1

Optional GIF, deliberately low-cost defaults:

    PYTHONPATH=src python3 \
      examples/animate_large_second_dyad_v41.py \
      --rank 1 --gif outputs/rank1.gif --frames 90 --fps 15 --dpi 90
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter, FuncAnimation
import numpy as np

ROOT = Path.cwd()
DEFAULT_RESULT = ROOT / "outputs" / "large_second_dyad_v41_panel_3952.json"


def load_json_relaxed(path: Path):
    return json.loads(path.read_text(encoding="utf-8").replace("NaN", "null"))


def load_rank_candidate(path: Path, rank: int) -> dict:
    data = load_json_relaxed(path)
    for family in data.get("families", []):
        if int(family.get("primary_rank", -999)) != rank:
            continue
        best = family.get("best")
        if best is None:
            raise ValueError(f"Primary rank {rank} exists but has no feasible downstream candidate.")
        return best
    raise ValueError(f"Primary rank {rank} not found in {path}")


def get_primary(best: dict) -> dict:
    p = best["primary"]
    return {
        "ground": float(p["ground"]),
        "coupler": float(p["coupler"]),
        "rocker": float(p["rocker"]),
        "assembly_branch": int(p["assembly_branch"]),
        "E_along": float(p["E_along"]),
        "E_normal": float(p["E_normal"]),
        "phase": float(p["phase"]),
    }


def primary_state(theta: float, primary: dict) -> dict[str, np.ndarray]:
    g = primary["ground"]
    c = primary["coupler"]
    r = primary["rocker"]
    branch = primary["assembly_branch"]

    angle = theta + primary["phase"]
    A = np.array((0.0, 0.0))
    B = np.array((math.cos(angle), math.sin(angle)))
    D = np.array((g, 0.0))

    delta = D - B
    distance = float(np.hypot(delta[0], delta[1]))
    if distance <= 1e-12:
        raise ValueError("Primary four-bar circle centers coincide.")

    direction = delta / distance
    along = (c*c - r*r + distance*distance) / (2.0*distance)
    height2 = c*c - along*along
    if height2 <= 1e-12:
        raise ValueError("Primary four-bar does not close.")

    height = math.sqrt(height2)
    left = np.array((-direction[1], direction[0]))
    C = B + along*direction + branch*height*left

    coupler_unit = (C - B) / c
    coupler_normal = np.array((-coupler_unit[1], coupler_unit[0]))
    E = B + primary["E_along"]*coupler_unit + primary["E_normal"]*coupler_normal

    return {"A": A, "B": B, "C": C, "D": D, "E": E}


def sixbar_state(theta: float, best: dict, primary: dict) -> dict[str, np.ndarray]:
    state = primary_state(theta, primary)
    E = state["E"]

    params = best["parameters"]
    second_branch = int(best["second_branch"])

    G = np.array((float(params["second_pivot_x"]), float(params["second_pivot_y"])))
    link_ef = float(params["link_EF"])
    link_gf = float(params["link_GF"])

    delta = G - E
    distance = float(np.hypot(delta[0], delta[1]))
    if distance <= 1e-12:
        raise ValueError("Second dyad circle centers coincide.")

    direction = delta / distance
    along = (link_ef*link_ef - link_gf*link_gf + distance*distance) / (2.0*distance)
    height2 = link_ef*link_ef - along*along
    if height2 <= 1e-12:
        raise ValueError("Second dyad does not close.")

    height = math.sqrt(height2)
    left = np.array((-direction[1], direction[0]))
    F = E + along*direction + second_branch*height*left

    ef_unit = (F - E) / link_ef
    ef_normal = np.array((-ef_unit[1], ef_unit[0]))
    H = (
        E
        + float(params["H_along_over_EF"])*link_ef*ef_unit
        + float(params["H_normal_over_EF"])*link_ef*ef_normal
    )

    rod = float(params["piston_rod"])
    axis_offset = float(params["slider_axis_offset"])
    axis_angle = float(params["slider_axis_angle"])

    axis = np.array((math.cos(axis_angle), math.sin(axis_angle)))
    normal = np.array((-axis[1], axis[0]))
    origin = axis_offset*normal

    relative = H - origin
    longitudinal = float(relative @ axis)
    transverse = float(relative @ normal)
    margin2 = rod*rod - transverse*transverse
    if margin2 <= 1e-12:
        raise ValueError("Finite piston rod cannot reach slider.")

    margin = math.sqrt(margin2)
    slider = longitudinal + margin
    P = origin + slider*axis
    Q = origin + longitudinal*axis

    state.update({
        "G": G,
        "F": F,
        "H": H,
        "P": P,
        "Q": Q,
        "axis_origin": origin,
        "axis_direction": axis,
    })
    return state


def build_samples(best: dict, frame_count: int):
    primary = get_primary(best)
    thetas = np.linspace(0.0, 2.0*math.pi, frame_count, endpoint=False)
    samples = [sixbar_state(float(theta), best, primary) for theta in thetas]
    return primary, thetas, samples


def bounds_from_samples(samples):
    points = []
    for state in samples:
        for name in ("A", "B", "C", "D", "E", "F", "G", "H", "P", "Q", "axis_origin"):
            points.append(state[name])
    cloud = np.array(points)
    xmin = float(np.min(cloud[:, 0])); xmax = float(np.max(cloud[:, 0]))
    ymin = float(np.min(cloud[:, 1])); ymax = float(np.max(cloud[:, 1]))
    span = max(xmax-xmin, ymax-ymin)
    margin = 0.12*span + 0.4
    return xmin-margin, xmax+margin, ymin-margin, ymax+margin, span+2*margin


def build_animation(result_json: Path, rank: int, frames: int, fps: int):
    best = load_rank_candidate(result_json, rank)
    primary, thetas, samples = build_samples(best, frames)
    xmin, xmax, ymin, ymax, full_span = bounds_from_samples(samples)

    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    ax.set_aspect("equal")
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
    ax.grid(True, alpha=0.25)
    ax.set_title(f"LARGE V4.1 second dyad — primary rank {rank}")

    primary_plate_color = "C0"
    second_plate_color = "C1"

    primary_ground_line, = ax.plot([], [], "-", lw=2.0)
    crank_line, = ax.plot([], [], "-", lw=2.5)
    coupler_line, = ax.plot([], [], "-", lw=2.0, color=primary_plate_color)
    coupler_plate_be, = ax.plot([], [], "-", lw=1.5, color=primary_plate_color)
    coupler_plate_ce, = ax.plot([], [], "-", lw=1.5, color=primary_plate_color)
    rocker_line, = ax.plot([], [], "-", lw=2.0)
    second_input_line, = ax.plot([], [], "-", lw=2.0, color=second_plate_color)
    second_plate_eh, = ax.plot([], [], "-", lw=1.5, color=second_plate_color)
    second_plate_fh, = ax.plot([], [], "-", lw=1.5, color=second_plate_color)
    second_ground_line, = ax.plot([], [], "-", lw=2.0)
    piston_rod_line, = ax.plot([], [], "-", lw=2.4)
    perpendicular_hint, = ax.plot([], [], "--", lw=1.0, alpha=0.5)
    slider_axis_line, = ax.plot([], [], "--", lw=1.6, alpha=0.7)
    h_trail_line, = ax.plot([], [], "-", lw=1.2, alpha=0.6)
    p_trail_line, = ax.plot([], [], "-", lw=1.2, alpha=0.6)

    joint_points = ax.scatter([], [])
    fixed_points = ax.scatter([], [], marker="s")
    special_points = ax.scatter([], [], marker="o")

    text = ax.text(
        0.02, 0.98, "", transform=ax.transAxes, va="top", ha="left",
        family="monospace", fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    h_trail = []
    p_trail = []

    def artists():
        return (
            primary_ground_line, crank_line,
            coupler_line, coupler_plate_be, coupler_plate_ce, rocker_line,
            second_input_line, second_plate_eh, second_plate_fh,
            second_ground_line, piston_rod_line, perpendicular_hint,
            slider_axis_line, h_trail_line, p_trail_line,
            joint_points, fixed_points, special_points, text,
        )

    def init():
        return artists()

    def update(i):
        s = samples[i]
        A, B, C, D, E = s["A"], s["B"], s["C"], s["D"], s["E"]
        G, F, H, P, Q = s["G"], s["F"], s["H"], s["P"], s["Q"]

        axis_origin = s["axis_origin"]
        axis_dir = s["axis_direction"]
        axis_half = 0.65*full_span
        axis_start = axis_origin - axis_half*axis_dir
        axis_end = axis_origin + axis_half*axis_dir

        primary_ground_line.set_data([A[0], D[0]], [A[1], D[1]])
        crank_line.set_data([A[0], B[0]], [A[1], B[1]])
        coupler_line.set_data([B[0], C[0]], [B[1], C[1]])
        coupler_plate_be.set_data([B[0], E[0]], [B[1], E[1]])
        coupler_plate_ce.set_data([C[0], E[0]], [C[1], E[1]])
        rocker_line.set_data([C[0], D[0]], [C[1], D[1]])
        second_input_line.set_data([E[0], F[0]], [E[1], F[1]])
        second_plate_eh.set_data([E[0], H[0]], [E[1], H[1]])
        second_plate_fh.set_data([F[0], H[0]], [F[1], H[1]])
        second_ground_line.set_data([G[0], F[0]], [G[1], F[1]])
        piston_rod_line.set_data([H[0], P[0]], [H[1], P[1]])
        perpendicular_hint.set_data([H[0], Q[0]], [H[1], Q[1]])
        slider_axis_line.set_data([axis_start[0], axis_end[0]], [axis_start[1], axis_end[1]])

        h_trail.append(H.copy()); p_trail.append(P.copy())
        trail_len = min(60, len(h_trail))
        hh = np.array(h_trail[-trail_len:]); pp = np.array(p_trail[-trail_len:])
        h_trail_line.set_data(hh[:, 0], hh[:, 1]); p_trail_line.set_data(pp[:, 0], pp[:, 1])

        joint_points.set_offsets(np.array([B, C, E, F]))
        fixed_points.set_offsets(np.array([A, D, G]))
        special_points.set_offsets(np.array([H, P]))

        theta_deg = math.degrees(float(thetas[i])) % 360.0
        text.set_text("\n".join([
            f"theta            {theta_deg:8.2f} deg",
            f"frame            {i+1:4d} / {len(samples):4d}",
            f"position RMS     {100*best['position_rms']:.4f} %",
            f"HP RMS           {100*best.get('position_rms_hp_band', float('nan')):.4f} %",
            f"non-HP RMS       {100*best.get('position_rms_non_hp', float('nan')):.4f} %",
            f"stroke / crank   {best['stroke_over_crank']:.4f}",
            f"rod / stroke     {best['piston_rod_over_stroke']:.4f}",
            f"primary sine     {best.get('minimum_primary_transmission_sine', float('nan')):.4f}",
            f"secondary sine   {best.get('minimum_secondary_transmission_sine', float('nan')):.4f}",
            "",
            f"AD               {primary['ground']:.4f}",
            f"BC               {primary['coupler']:.4f}",
            f"CD               {primary['rocker']:.4f}",
        ]))
        return artists()

    ani = FuncAnimation(
        fig, update, frames=len(samples), init_func=init,
        interval=1000.0/fps, blit=False, repeat=True,
    )
    return fig, ani


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    ap.add_argument("--rank", type=int, default=1)
    ap.add_argument("--frames", type=int, default=120)
    ap.add_argument("--fps", type=int, default=20)
    ap.add_argument("--gif", type=Path, default=None,
                    help="Optional GIF output. If omitted: interactive display only.")
    ap.add_argument("--dpi", type=int, default=90)
    args = ap.parse_args()

    fig, ani = build_animation(args.result, args.rank, args.frames, args.fps)

    if args.gif is None:
        print(f"Interactive animation: rank {args.rank}. Close the matplotlib window to exit.")
        plt.show()
        return

    args.gif.parent.mkdir(parents=True, exist_ok=True)
    ani.save(args.gif, writer=PillowWriter(fps=args.fps), dpi=args.dpi)
    plt.close(fig)
    print(f"Wrote {args.gif}")


if __name__ == "__main__":
    main()
