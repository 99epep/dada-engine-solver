#!/usr/bin/env python3
"""Animate the small and large six-bar mechanisms together, face-to-face.

Layout:
    [mechanism S] -- [cylinder S]   exchangers + check valves   [cylinder L] -- [mechanism L]

Features:
- both mechanisms animated on the same cycle;
- cylinders aligned on a common horizontal axis;
- cylinder back-heads face each other;
- mechanisms kept outside the cylinders;
- no axes, no grid, no text;
- central transfer block with two exchangers and two opposite check valves;
- piston/cylinder width of S = 0.8 * width of L.

Example:

PYTHONPATH=src python3 examples/animate_sixbar_pair.py \
    outputs/small_sixbar_stage2f_r4_freeH_tightaxis.json \
    outputs/large_sixbar_stageL1.json \
    --small-restart 0 \
    --small-branch +1 \
    --output outputs/sixbar_pair.gif
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Circle, Polygon, Rectangle
import numpy as np


DEFAULT_PRIMARY = {
    "ground": 2.045193376592419,
    "coupler": 2.1152136168766678,
    "rocker": 1.8202134329709354,
    "assembly_branch": -1,
    "E_along": 1.873315626010214,
    "E_normal": -1.5362965718510164,
    "phase": -3.0863675667472057,
}


def load_candidate(path: Path, restart: int | None, branch: int | None) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))

    if restart is not None or branch is not None:
        if "runs" not in data:
            raise ValueError(f"{path} has no runs[] section for run selection.")
        target = None
        for run in data["runs"]:
            ok = True
            if restart is not None:
                ok = ok and (int(run.get("restart")) == int(restart))
            if branch is not None:
                ok = ok and (int(run.get("second_branch")) == int(branch))
            if ok and run.get("best_dense") is not None:
                target = run["best_dense"]
                break
        if target is None:
            raise ValueError(
                f"Could not find restart={restart!r}, branch={branch!r} in {path}"
            )
        return target

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

    raise ValueError(f"Could not find a candidate in {path}")


def get_primary(best: dict) -> dict:
    p = best.get("primary")
    if p is None:
        return DEFAULT_PRIMARY.copy()

    return {
        "ground": float(p["ground"]),
        "coupler": float(p["coupler"]),
        "rocker": float(p["rocker"]),
        "assembly_branch": int(p.get("assembly_branch", -1)),
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

    E = (
        B
        + primary["E_along"] * coupler_unit
        + primary["E_normal"] * coupler_normal
    )

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
    along = (
        link_ef*link_ef - link_gf*link_gf + distance*distance
    ) / (2.0*distance)

    height2 = link_ef*link_ef - along*along
    if height2 <= 1e-12:
        raise ValueError("Second dyad does not close.")

    height = math.sqrt(height2)
    left = np.array((-direction[1], direction[0]))
    F = E + along*direction + second_branch*height*left

    if "H_along_over_EF" in params and "H_normal_over_EF" in params:
        ef_unit = (F - E) / link_ef
        ef_normal = np.array((-ef_unit[1], ef_unit[0]))
        H = (
            E
            + float(params["H_along_over_EF"]) * link_ef * ef_unit
            + float(params["H_normal_over_EF"]) * link_ef * ef_normal
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

    margin2 = rod*rod - transverse*transverse
    if margin2 <= 1e-12:
        raise ValueError("Finite piston rod cannot reach slider.")

    margin = math.sqrt(margin2)
    slider = longitudinal + margin
    P = origin + slider*axis

    state.update({
        "G": G,
        "F": F,
        "H": H,
        "P": P,
        "axis_origin": origin,
        "axis_direction": axis,
        "axis_angle": axis_angle,
        "rod_length": rod,
    })
    return state


POINT_KEYS = ["A", "B", "C", "D", "E", "F", "G", "H", "P", "axis_origin"]


def rotate_point(p: np.ndarray, angle: float) -> np.ndarray:
    c = math.cos(angle)
    s = math.sin(angle)
    return np.array((c*p[0] - s*p[1], s*p[0] + c*p[1]))


def align_state(state: dict[str, np.ndarray], mirror_x: bool) -> dict[str, np.ndarray]:
    angle = -float(state["axis_angle"])
    out = {}
    for key in POINT_KEYS:
        out[key] = rotate_point(state[key], angle)

    axis_dir = rotate_point(state["axis_direction"], angle)
    if mirror_x:
        for key in POINT_KEYS:
            out[key] = np.array((-out[key][0], out[key][1]))
        axis_dir = np.array((-axis_dir[0], axis_dir[1]))

    out["axis_direction"] = axis_dir
    out["rod_length"] = state["rod_length"]
    return out


def build_aligned_samples(best: dict, frame_count: int, side: str):
    primary = get_primary(best)
    thetas = np.linspace(0.0, 2.0*math.pi, frame_count, endpoint=False)
    mirror_x = (side == "right")
    samples = [align_state(sixbar_state(float(theta), best, primary), mirror_x) for theta in thetas]
    return samples


def place_samples(samples, side: str, inner_head_x: float, cyl_width: float):
    p_values = np.array([s["P"][0] for s in samples], dtype=float)
    pmin = float(np.min(p_values))
    pmax = float(np.max(p_values))
    stroke = max(1e-9, pmax - pmin)

    piston_len = 1.15 * cyl_width
    wall = 0.12 * cyl_width
    chamber_clear = 0.35 * cyl_width
    rod_clear = 0.50 * cyl_width

    if side == "left":
        tx = inner_head_x - chamber_clear - piston_len/2.0 - pmax
    else:
        tx = inner_head_x + chamber_clear + piston_len/2.0 - pmin

    placed = []
    for s in samples:
        d = {}
        for key, value in s.items():
            if isinstance(value, np.ndarray):
                d[key] = value + np.array((tx, 0.0))
            else:
                d[key] = value
        placed.append(d)

    p_world = np.array([s["P"][0] for s in placed], dtype=float)
    pmin_w = float(np.min(p_world))
    pmax_w = float(np.max(p_world))

    if side == "left":
        x_inner = inner_head_x
        x_outer = pmin_w - piston_len/2.0 - rod_clear
    else:
        x_inner = inner_head_x
        x_outer = pmax_w + piston_len/2.0 + rod_clear

    cylinder = {
        "side": side,
        "x_inner": x_inner,
        "x_outer": x_outer,
        "width": cyl_width,
        "piston_length": piston_len,
        "wall": wall,
        "stroke": stroke,
    }
    return placed, cylinder


def draw_check_valve(ax, center, size, direction, lw=3.0):
    x, y = center
    r = 0.24 * size

    circ = Circle((x, y), r, fill=False, lw=lw, color="black")
    ax.add_patch(circ)

    jaw = 0.40 * size
    if direction == "right":
        ax.plot([x - 0.45*size, x - 0.05*size], [y + jaw, y], color="black", lw=lw)
        ax.plot([x - 0.45*size, x - 0.05*size], [y - jaw, y], color="black", lw=lw)
        tri = Polygon(
            [(x + 0.02*size, y - 0.10*size), (x + 0.20*size, y), (x + 0.02*size, y + 0.10*size)],
            closed=True, fill=False, edgecolor="black", lw=lw,
        )
    else:
        ax.plot([x + 0.45*size, x + 0.05*size], [y + jaw, y], color="black", lw=lw)
        ax.plot([x + 0.45*size, x + 0.05*size], [y - jaw, y], color="black", lw=lw)
        tri = Polygon(
            [(x - 0.02*size, y - 0.10*size), (x - 0.20*size, y), (x - 0.02*size, y + 0.10*size)],
            closed=True, fill=False, edgecolor="black", lw=lw,
        )

    ax.add_patch(tri)


def draw_exchanger(ax, x0, x1, y, amplitude, turns, lw=3.0):
    xs = np.linspace(x0, x1, 2*turns + 1)
    ys = np.full_like(xs, y)
    for i in range(1, len(xs) - 1):
        ys[i] = y + amplitude if i % 2 else y - amplitude
    ax.plot(xs, ys, color="black", lw=lw)


def draw_transfer_block(ax, x_left, x_right, large_width, lw=3.2):
    h = 0.72 * large_width
    y_top = 0.80 * large_width
    y_bot = -0.80 * large_width

    manifold = 0.16 * (x_right - x_left)
    xl = x_left + 0.05 * (x_right - x_left)
    xr = x_right - 0.05 * (x_right - x_left)

    # Manifolds connected to cylinder heads.
    ax.plot([x_left, xl], [0, 0], color="black", lw=lw)
    ax.plot([xl, xl], [y_bot, y_top], color="black", lw=lw)
    ax.plot([x_right, xr], [0, 0], color="black", lw=lw)
    ax.plot([xr, xr], [y_bot, y_top], color="black", lw=lw)

    # Branches.
    top_y = y_top
    bot_y = y_bot
    span = xr - xl

    # Top branch: valve then exchanger, left -> right.
    vt = xl + 0.20 * span
    et0 = xl + 0.46 * span
    et1 = xl + 0.78 * span

    ax.plot([xl, vt - 0.32*large_width], [top_y, top_y], color="black", lw=lw)
    draw_check_valve(ax, (vt, top_y), 0.78*large_width, "right", lw=lw)
    ax.plot([vt + 0.25*large_width, et0], [top_y, top_y], color="black", lw=lw)
    draw_exchanger(ax, et0, et1, top_y, 0.42*large_width, 6, lw=lw)
    ax.plot([et1, xr], [top_y, top_y], color="black", lw=lw)

    # Bottom branch: exchanger then valve, right -> left.
    vb = xr - 0.20 * span
    eb0 = xl + 0.22 * span
    eb1 = xl + 0.54 * span

    ax.plot([xl, eb0], [bot_y, bot_y], color="black", lw=lw)
    draw_exchanger(ax, eb0, eb1, bot_y, 0.42*large_width, 6, lw=lw)
    ax.plot([eb1, vb - 0.25*large_width], [bot_y, bot_y], color="black", lw=lw)
    draw_check_valve(ax, (vb, bot_y), 0.78*large_width, "left", lw=lw)
    ax.plot([vb + 0.32*large_width, xr], [bot_y, bot_y], color="black", lw=lw)


def draw_cylinder(ax, cyl, piston_center_x, piston_facecolor, body_facecolor="#f1f1f1", lw=3.0):
    side = cyl["side"]
    x0 = min(cyl["x_inner"], cyl["x_outer"])
    x1 = max(cyl["x_inner"], cyl["x_outer"])
    w = cyl["width"]
    piston_len = cyl["piston_length"]

    # Cylinder body
    body = Rectangle((x0, -w/2.0), x1 - x0, w, facecolor=body_facecolor, edgecolor="black", lw=lw)
    ax.add_patch(body)

    # Back head (facing the exchanger)
    if side == "left":
        ax.plot([cyl["x_inner"], cyl["x_inner"]], [-w/2.0, w/2.0], color="black", lw=lw)
    else:
        ax.plot([cyl["x_inner"], cyl["x_inner"]], [-w/2.0, w/2.0], color="black", lw=lw)

    # Rod-side guide / packing
    gland_x = cyl["x_outer"] + (0.14*w if side == "left" else -0.14*w)
    for yy in np.linspace(-0.36*w, 0.36*w, 6):
        ax.plot([gland_x - 0.03*w, gland_x + 0.03*w], [yy, yy], color="black", lw=lw*0.55)

    # Piston
    piston = Rectangle(
        (piston_center_x - piston_len/2.0, -0.46*w),
        piston_len,
        0.92*w,
        facecolor=piston_facecolor,
        edgecolor="none",
    )
    ax.add_patch(piston)


def draw_mechanism(ax, s, side: str, color="black", lw=2.8):
    A, B, C, D, E = s["A"], s["B"], s["C"], s["D"], s["E"]
    G, F, H, P = s["G"], s["F"], s["H"], s["P"]

    # Primary
    ax.plot([A[0], D[0]], [A[1], D[1]], color=color, lw=lw)
    ax.plot([A[0], B[0]], [A[1], B[1]], color=color, lw=lw)
    ax.plot([B[0], C[0]], [B[1], C[1]], color=color, lw=lw)
    ax.plot([C[0], D[0]], [C[1], D[1]], color=color, lw=lw)
    ax.plot([B[0], E[0]], [B[1], E[1]], color=color, lw=lw*0.9)
    ax.plot([C[0], E[0]], [C[1], E[1]], color=color, lw=lw*0.9)

    # Secondary
    ax.plot([E[0], F[0]], [E[1], F[1]], color=color, lw=lw)
    ax.plot([G[0], F[0]], [G[1], F[1]], color=color, lw=lw)
    ax.plot([E[0], H[0]], [E[1], H[1]], color=color, lw=lw*0.9)
    ax.plot([F[0], H[0]], [F[1], H[1]], color=color, lw=lw*0.9)

    # Piston rod
    ax.plot([H[0], P[0]], [H[1], P[1]], color=color, lw=lw)

    # Small joint dots
    for pt in (A, B, C, D, E, F, G, H):
        circ = Circle((pt[0], pt[1]), 0.055, facecolor="white", edgecolor="black", lw=lw*0.45)
        ax.add_patch(circ)


def collect_cloud(samples):
    pts = []
    for s in samples:
        for key in ("A", "B", "C", "D", "E", "F", "G", "H", "P"):
            pts.append(s[key])
    return np.array(pts)


def make_animation(
    small_json: Path,
    large_json: Path,
    output_gif: Path,
    small_restart: int | None,
    small_branch: int | None,
    large_restart: int | None,
    large_branch: int | None,
    frames: int,
    fps: int,
    dpi: int,
    inter_cylinder_gap: float,
    large_cylinder_width: float,
    small_width_ratio: float,
):
    small_best = load_candidate(small_json, small_restart, small_branch)
    large_best = load_candidate(large_json, large_restart, large_branch)

    small_samples_local = build_aligned_samples(small_best, frames, "left")
    large_samples_local = build_aligned_samples(large_best, frames, "right")

    small_width = small_width_ratio * large_cylinder_width
    x_left_head = -inter_cylinder_gap / 2.0
    x_right_head = inter_cylinder_gap / 2.0

    small_samples, small_cyl = place_samples(small_samples_local, "left", x_left_head, small_width)
    large_samples, large_cyl = place_samples(large_samples_local, "right", x_right_head, large_cylinder_width)

    cloud = np.vstack((collect_cloud(small_samples), collect_cloud(large_samples)))

    x_extra = [
        small_cyl["x_outer"], small_cyl["x_inner"],
        large_cyl["x_outer"], large_cyl["x_inner"],
        -inter_cylinder_gap/2.0, inter_cylinder_gap/2.0
    ]
    y_extra = [
        -large_cylinder_width, large_cylinder_width,
        -large_cylinder_width, large_cylinder_width
    ]

    xmin = min(float(np.min(cloud[:, 0])), min(x_extra))
    xmax = max(float(np.max(cloud[:, 0])), max(x_extra))
    ymin = min(float(np.min(cloud[:, 1])), min(y_extra))
    ymax = max(float(np.max(cloud[:, 1])), max(y_extra))

    spanx = xmax - xmin
    spany = ymax - ymin
    span = max(spanx, spany)
    margin = 0.08 * span + 0.8

    fig, ax = plt.subplots(figsize=(12.8, 5.4))
    ax.set_aspect("equal")
    ax.set_xlim(xmin - margin, xmax + margin)
    ax.set_ylim(ymin - margin, ymax + margin)
    ax.axis("off")

    def update(i):
        ax.clear()
        ax.set_aspect("equal")
        ax.set_xlim(xmin - margin, xmax + margin)
        ax.set_ylim(ymin - margin, ymax + margin)
        ax.axis("off")

        draw_transfer_block(ax, small_cyl["x_inner"], large_cyl["x_inner"], large_cylinder_width, lw=3.0)

        sS = small_samples[i]
        sL = large_samples[i]

        draw_cylinder(ax, small_cyl, sS["P"][0], piston_facecolor="#e02b2b", body_facecolor="#eeeeee", lw=3.0)
        draw_cylinder(ax, large_cyl, sL["P"][0], piston_facecolor="#c71515", body_facecolor="#eeeeee", lw=3.0)

        draw_mechanism(ax, sS, "left", lw=2.5)
        draw_mechanism(ax, sL, "right", lw=2.5)

        return ()

    anim = FuncAnimation(
        fig,
        update,
        frames=frames,
        interval=1000/fps,
        blit=False,
        repeat=True,
    )

    output_gif.parent.mkdir(parents=True, exist_ok=True)
    anim.save(output_gif, writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("small_json", type=Path)
    ap.add_argument("large_json", type=Path)

    ap.add_argument("--small-restart", type=int, default=None)
    ap.add_argument("--small-branch", type=int, choices=[-1, 1], default=None)
    ap.add_argument("--large-restart", type=int, default=None)
    ap.add_argument("--large-branch", type=int, choices=[-1, 1], default=None)

    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--frames", type=int, default=120)
    ap.add_argument("--fps", type=int, default=20)
    ap.add_argument("--dpi", type=int, default=110)

    ap.add_argument("--inter-cylinder-gap", type=float, default=6.0)
    ap.add_argument("--large-cylinder-width", type=float, default=1.35)
    ap.add_argument("--small-width-ratio", type=float, default=0.8)

    args = ap.parse_args()

    output = args.output if args.output else Path("outputs/sixbar_pair.gif")
    make_animation(
        small_json=args.small_json,
        large_json=args.large_json,
        output_gif=output,
        small_restart=args.small_restart,
        small_branch=args.small_branch,
        large_restart=args.large_restart,
        large_branch=args.large_branch,
        frames=args.frames,
        fps=args.fps,
        dpi=args.dpi,
        inter_cylinder_gap=args.inter_cylinder_gap,
        large_cylinder_width=args.large_cylinder_width,
        small_width_ratio=args.small_width_ratio,
    )
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
