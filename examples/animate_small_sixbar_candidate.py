#!/usr/bin/env python3
"""Animate a synthesized small-cylinder mechanism and save it as a GIF.

Supports:
- Stage 2A JSON (H implicitly equal to F)
- Stage 2B JSON (free point H on EF)

Typical usage from dada-engine-solver:

    PYTHONPATH=src python3 examples/animate_small_sixbar_candidate.py \
        outputs/small_sixbar_stage2b.json \
        --output outputs/small_sixbar_stage2b.gif
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter, FuncAnimation
import numpy as np


PRIMARY = {
    "ground": 2.045193376592419,
    "coupler": 2.1152136168766678,
    "rocker": 1.8202134329709354,
    "branch": -1,
    "e_along": 1.873315626010214,
    "e_normal": -1.5362965718510164,
    "phase": -3.0863675667472057,
}


def load_best_candidate(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))

    if isinstance(data, dict) and data.get("best") is not None:
        return data["best"]

    if isinstance(data, dict) and "runs" in data:
        candidates = [
            run.get("best_dense")
            for run in data["runs"]
            if run.get("best_dense") is not None
        ]
        if candidates:
            return min(candidates, key=lambda item: item["score"])

    raise ValueError("Could not find a best candidate in the JSON result.")


def primary_state(theta: float) -> dict[str, np.ndarray]:
    g = PRIMARY["ground"]
    c = PRIMARY["coupler"]
    r = PRIMARY["rocker"]
    branch = PRIMARY["branch"]

    angle = theta + PRIMARY["phase"]
    A = np.array((0.0, 0.0))
    B = np.array((math.cos(angle), math.sin(angle)))
    D = np.array((g, 0.0))

    delta = D - B
    distance = float(np.hypot(delta[0], delta[1]))
    if distance <= 1e-12:
        raise ValueError("Primary four-bar circle centers coincide.")

    direction = delta / distance
    along = (c * c - r * r + distance * distance) / (2.0 * distance)
    height2 = c * c - along * along
    if height2 <= 1e-12:
        raise ValueError("Primary four-bar does not close.")
    height = math.sqrt(height2)

    left = np.array((-direction[1], direction[0]))
    C = B + along * direction + branch * height * left

    coupler_unit = (C - B) / c
    coupler_normal = np.array((-coupler_unit[1], coupler_unit[0]))
    E = (
        B
        + PRIMARY["e_along"] * coupler_unit
        + PRIMARY["e_normal"] * coupler_normal
    )

    return {"A": A, "B": B, "C": C, "D": D, "E": E}


def sixbar_state(theta: float, best: dict) -> dict[str, np.ndarray]:
    primary = primary_state(theta)
    E = primary["E"]

    params = best["parameters"]
    second_branch = int(best["second_branch"])

    G = np.array((params["second_pivot_x"], params["second_pivot_y"]))
    link_ef = float(params["link_EF"])
    link_gf = float(params["link_GF"])

    delta = G - E
    distance = float(np.hypot(delta[0], delta[1]))
    if distance <= 1e-12:
        raise ValueError("Second dyad circle centers coincide.")

    direction = delta / distance
    along = (link_ef * link_ef - link_gf * link_gf + distance * distance) / (
        2.0 * distance
    )
    height2 = link_ef * link_ef - along * along
    if height2 <= 1e-12:
        raise ValueError("Second dyad does not close.")
    height = math.sqrt(height2)

    left = np.array((-direction[1], direction[0]))
    F = E + along * direction + second_branch * height * left

    if "H_along_over_EF" in params and "H_normal_over_EF" in params:
        ef_unit = (F - E) / link_ef
        ef_normal = np.array((-ef_unit[1], ef_unit[0]))
        H = (
            E
            + params["H_along_over_EF"] * link_ef * ef_unit
            + params["H_normal_over_EF"] * link_ef * ef_normal
        )
    else:
        H = F.copy()

    rod = float(params["piston_rod"])
    axis_offset = float(params["slider_axis_offset"])
    axis_angle = float(params["slider_axis_angle"])

    axis = np.array((math.cos(axis_angle), math.sin(axis_angle)))
    normal = np.array((-axis[1], axis[0]))
    origin = axis_offset * normal

    relative = H - origin
    longitudinal = float(relative @ axis)
    transverse = float(relative @ normal)

    margin2 = rod * rod - transverse * transverse
    if margin2 <= 1e-12:
        raise ValueError("Finite piston rod cannot reach slider.")
    margin = math.sqrt(margin2)

    slider = longitudinal + margin
    P = origin + slider * axis
    Q = origin + longitudinal * axis

    out = dict(primary)
    out.update(
        {
            "G": G,
            "F": F,
            "H": H,
            "P": P,
            "Q": Q,
            "axis_origin": origin,
            "axis_direction": axis,
            "axis_normal": normal,
        }
    )
    return out


def build_samples(best: dict, frame_count: int) -> list[dict]:
    thetas = np.linspace(0.0, 2.0 * math.pi, frame_count, endpoint=False)
    return [sixbar_state(float(theta), best) for theta in thetas]


def bounds_from_samples(samples: list[dict]) -> tuple[float, float, float, float, float]:
    points = []
    for sample in samples:
        for name in ("A", "B", "C", "D", "E", "F", "G", "H", "P", "Q", "axis_origin"):
            points.append(sample[name])

    cloud = np.array(points)
    xmin = float(np.min(cloud[:, 0]))
    xmax = float(np.max(cloud[:, 0]))
    ymin = float(np.min(cloud[:, 1]))
    ymax = float(np.max(cloud[:, 1]))

    span = max(xmax - xmin, ymax - ymin)
    margin = 0.15 * span + 0.5
    return xmin - margin, xmax + margin, ymin - margin, ymax + margin, span + 2 * margin


def make_animation(result_json: Path, output_gif: Path, frames: int, fps: int, dpi: int) -> None:
    best = load_best_candidate(result_json)
    samples = build_samples(best, frames)

    xmin, xmax, ymin, ymax, full_span = bounds_from_samples(samples)

    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    ax.set_aspect("equal")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.grid(True, alpha=0.25)
    ax.set_title(result_json.stem)

    primary_ground_line, = ax.plot([], [], "-", lw=2.0)
    crank_line, = ax.plot([], [], "-", lw=2.5)
    coupler_line, = ax.plot([], [], "-", lw=2.0)
    rocker_line, = ax.plot([], [], "-", lw=2.0)

    second_input_line, = ax.plot([], [], "-", lw=2.0)
    second_ground_line, = ax.plot([], [], "-", lw=2.0)

    piston_rod_line, = ax.plot([], [], "-", lw=2.4)
    perpendicular_hint, = ax.plot([], [], "--", lw=1.0, alpha=0.5)
    slider_axis_line, = ax.plot([], [], "--", lw=1.6, alpha=0.7)

    coupler_plate_be, = ax.plot([], [], "-", lw=1.4)
    coupler_plate_ce, = ax.plot([], [], "-", lw=1.4)

    second_plate_eh, = ax.plot([], [], "-", lw=1.4)
    second_plate_fh, = ax.plot([], [], "-", lw=1.4)

    h_trail_line, = ax.plot([], [], "-", lw=1.2, alpha=0.6)
    p_trail_line, = ax.plot([], [], "-", lw=1.2, alpha=0.6)

    joint_points = ax.scatter([], [])
    fixed_points = ax.scatter([], [], marker="s")
    special_points = ax.scatter([], [], marker="o")

    text = ax.text(
        0.02,
        0.98,
        "",
        transform=ax.transAxes,
        va="top",
        ha="left",
        family="monospace",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    h_trail = []
    p_trail = []

    def init():
        return (
            primary_ground_line,
            crank_line,
            coupler_line,
            coupler_plate_be,
            coupler_plate_ce,
            rocker_line,
            second_input_line,
            second_ground_line,
            second_plate_eh,
            second_plate_fh,
            piston_rod_line,
            perpendicular_hint,
            slider_axis_line,
            h_trail_line,
            p_trail_line,
            joint_points,
            fixed_points,
            special_points,
            text,
        )

    def update(i: int):
        sample = samples[i]

        A, B, C, D, E = sample["A"], sample["B"], sample["C"], sample["D"], sample["E"]
        G, F, H, P, Q = sample["G"], sample["F"], sample["H"], sample["P"], sample["Q"]
        axis_origin = sample["axis_origin"]
        axis_dir = sample["axis_direction"]

        axis_half = 0.65 * full_span
        axis_start = axis_origin - axis_half * axis_dir
        axis_end = axis_origin + axis_half * axis_dir

        primary_ground_line.set_data([A[0], D[0]], [A[1], D[1]])
        crank_line.set_data([A[0], B[0]], [A[1], B[1]])
        coupler_line.set_data([B[0], C[0]], [B[1], C[1]])
        rocker_line.set_data([C[0], D[0]], [C[1], D[1]])

        second_input_line.set_data([E[0], F[0]], [E[1], F[1]])
        second_ground_line.set_data([G[0], F[0]], [G[1], F[1]])

        piston_rod_line.set_data([H[0], P[0]], [H[1], P[1]])
        perpendicular_hint.set_data([H[0], Q[0]], [H[1], Q[1]])
        slider_axis_line.set_data([axis_start[0], axis_end[0]], [axis_start[1], axis_end[1]])

        coupler_plate_be.set_data([B[0], E[0]], [B[1], E[1]])
        coupler_plate_ce.set_data([C[0], E[0]], [C[1], E[1]])

        second_plate_eh.set_data([E[0], H[0]], [E[1], H[1]])
        second_plate_fh.set_data([F[0], H[0]], [F[1], H[1]])

        h_trail.append(H.copy())
        p_trail.append(P.copy())
        trail_len = min(60, len(h_trail))
        h_hist = np.array(h_trail[-trail_len:])
        p_hist = np.array(p_trail[-trail_len:])
        h_trail_line.set_data(h_hist[:, 0], h_hist[:, 1])
        p_trail_line.set_data(p_hist[:, 0], p_hist[:, 1])

        joints = np.array([B, C, E, F])
        fixed = np.array([A, D, G])
        special = np.array([H, P])

        joint_points.set_offsets(joints)
        fixed_points.set_offsets(fixed)
        special_points.set_offsets(special)

        text.set_text(
            "\n".join(
                [
                    f"frame            {i+1:4d} / {len(samples):4d}",
                    f"position RMS     {100*best['position_rms']:.4f} %",
                    f"dRMS             {best['derivative_rms']:.6f}",
                    f"stroke / crank   {best['stroke_over_crank']:.4f}",
                    f"Hlat / stroke    {best.get('H_lateral_span_over_piston_stroke', float('nan')):.4f}",
                    f"2nd trans. sine  {best.get('minimum_secondary_transmission_sine', float('nan')):.4f}",
                    f"rod-axis cosine  {best.get('minimum_rod_axis_cosine', float('nan')):.4f}",
                ]
            )
        )

        return (
            primary_ground_line,
            crank_line,
            coupler_line,
            rocker_line,
            coupler_plate_be,
            coupler_plate_ce,
            second_input_line,
            second_ground_line,
            second_plate_eh,
            second_plate_fh,
            piston_rod_line,
            perpendicular_hint,
            slider_axis_line,
            h_trail_line,
            p_trail_line,
            joint_points,
            fixed_points,
            special_points,
            text,
        )

    animation = FuncAnimation(
        fig,
        update,
        frames=len(samples),
        init_func=init,
        interval=1000 / fps,
        blit=False,
        repeat=True,
    )

    output_gif.parent.mkdir(parents=True, exist_ok=True)
    animation.save(output_gif, writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("result_json", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--frames", type=int, default=120)
    parser.add_argument("--fps", type=int, default=20)
    parser.add_argument("--dpi", type=int, default=110)
    args = parser.parse_args()

    if args.frames < 16:
        parser.error("--frames must be at least 16.")
    if args.fps < 1:
        parser.error("--fps must be positive.")

    output = args.output if args.output is not None else args.result_json.with_suffix(".gif")
    make_animation(args.result_json, output, args.frames, args.fps, args.dpi)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
