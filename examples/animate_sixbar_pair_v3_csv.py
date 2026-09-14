#!/usr/bin/env python3
"""Animate the small and large six-bar mechanisms together, face-to-face (V3 CSV).

Differences vs V3:
- thermodynamic temperatures are read from a cycle CSV;
- optional metadata JSON can be provided;
- automatic or manual column selection.

Example:

PYTHONPATH=src python3 examples/animate_sixbar_pair_v3_csv.py \
    outputs/champion-6bar_thermodynamic_cycle.csv \
    outputs/champion-6bar_thermodynamic_cycle_metadata.json \
    outputs/small_sixbar_stage2f_r4_freeH_tightaxis.json \
    outputs/large_sixbar_stageL1.json \
    --small-restart 0 \
    --small-branch +1 \
    --output outputs/sixbar_pair_v3.gif
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
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

COLORS = {
    "crank": "#1f77b4",
    "coupler": "#ff7f0e",
    "rocker": "#2ca02c",
    "secondary_input": "#9467bd",
    "secondary_plate": "#d62728",
    "rod": "#444444",
    "ground": "#111111",
    "joint_fill": "#ffffff",
    "piston": "#555555",
    "pipe": "#111111",
}


def load_candidate(path: Path, restart: int | None, branch: int | None) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))

    if restart is not None or branch is not None:
        if "runs" not in data:
            raise ValueError(f"{path} has no runs[] section for run selection.")
        for run in data["runs"]:
            ok = True
            if restart is not None:
                ok = ok and (int(run.get("restart")) == int(restart))
            if branch is not None:
                ok = ok and (int(run.get("second_branch")) == int(branch))
            if ok and run.get("best_dense") is not None:
                return run["best_dense"]
        raise ValueError(f"Could not find restart={restart!r}, branch={branch!r} in {path}")

    if isinstance(data, dict) and data.get("best") is not None:
        return data["best"]

    if isinstance(data, dict) and "runs" in data:
        candidates = [run.get("best_dense") for run in data["runs"] if run.get("best_dense") is not None]
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
    direction = delta / distance
    along = (c*c - r*r + distance*distance) / (2.0*distance)
    height = math.sqrt(max(1e-12, c*c - along*along))
    left = np.array((-direction[1], direction[0]))
    C = B + along*direction + branch*height*left

    coupler_unit = (C - B) / c
    coupler_normal = np.array((-coupler_unit[1], coupler_unit[0]))
    E = B + primary["E_along"] * coupler_unit + primary["E_normal"] * coupler_normal
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
    direction = delta / distance
    along = (link_ef*link_ef - link_gf*link_gf + distance*distance) / (2.0*distance)
    height = math.sqrt(max(1e-12, link_ef*link_ef - along*along))
    left = np.array((-direction[1], direction[0]))
    F = E + along*direction + second_branch*height*left

    ef_unit = (F - E) / link_ef
    ef_normal = np.array((-ef_unit[1], ef_unit[0]))
    H = (
        E
        + float(params["H_along_over_EF"]) * link_ef * ef_unit
        + float(params["H_normal_over_EF"]) * link_ef * ef_normal
    )

    rod = float(params["piston_rod"])
    axis_offset = float(params["slider_axis_offset"])
    axis_angle = float(params["slider_axis_angle"])

    axis = np.array((math.cos(axis_angle), math.sin(axis_angle)))
    normal = np.array((-axis[1], axis[0]))
    origin = axis_offset * normal

    relative = H - origin
    longitudinal = float(relative @ axis)
    transverse = float(relative @ normal)

    margin = math.sqrt(max(1e-12, rod*rod - transverse*transverse))
    slider = longitudinal + margin
    P = origin + slider*axis

    state.update({
        "G": G, "F": F, "H": H, "P": P,
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

    y_shift = -out["axis_origin"][1]
    for key in POINT_KEYS:
        out[key] = out[key] + np.array((0.0, y_shift))

    if mirror_x:
        for key in POINT_KEYS:
            out[key] = np.array((-out[key][0], out[key][1]))
        axis_dir = np.array((-axis_dir[0], axis_dir[1]))

    out["axis_direction"] = axis_dir
    out["rod_length"] = state["rod_length"]
    return out


def build_aligned_samples(best: dict, frame_count: int, side: str):
    primary = get_primary(best)
    base_thetas = np.linspace(0.0, 2.0*math.pi, frame_count, endpoint=False)
    mirror_x = (side == "right")
    eval_thetas = (-base_thetas) % (2.0*math.pi) if side == "right" else base_thetas
    return [align_state(sixbar_state(float(theta), best, primary), mirror_x) for theta in eval_thetas]


def place_samples(samples, side: str, inner_head_x: float, cyl_width: float):
    p_values = np.array([s["P"][0] for s in samples], dtype=float)
    pmin = float(np.min(p_values))
    pmax = float(np.max(p_values))
    rod_clear = 0.70

    tx = (inner_head_x - pmax) if side == "left" else (inner_head_x - pmin)

    placed = []
    for s in samples:
        d = {}
        for key, value in s.items():
            d[key] = value + np.array((tx, 0.0)) if isinstance(value, np.ndarray) else value
        placed.append(d)

    p_world = np.array([s["P"][0] for s in placed], dtype=float)
    pmin_w = float(np.min(p_world))
    pmax_w = float(np.max(p_world))

    if side == "left":
        x_inner = inner_head_x
        x_outer = pmin_w - rod_clear
    else:
        x_inner = inner_head_x
        x_outer = pmax_w + rod_clear

    return placed, {
        "side": side,
        "x_inner": x_inner,
        "x_outer": x_outer,
        "width": cyl_width,
    }


def load_cycle_csv(csv_path: Path):
    with csv_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        columns = {name: [] for name in fieldnames}
        for row in reader:
            for name in fieldnames:
                val = row.get(name, "")
                try:
                    columns[name].append(float(val))
                except Exception:
                    columns[name].append(np.nan)
    arrays = {k: np.asarray(v, dtype=float) for k, v in columns.items()}
    return arrays


def score_name(name: str, patterns: Iterable[str]) -> int:
    n = name.lower()
    return sum(10 for p in patterns if p in n)


def pick_column(columns: dict[str, np.ndarray], patterns: list[str]):
    scored = []
    for name, arr in columns.items():
        if np.all(~np.isfinite(arr)):
            continue
        s = score_name(name, patterns)
        if s > 0:
            scored.append((s, name))
    if not scored:
        return None
    scored.sort(key=lambda x: (-x[0], len(x[1])))
    return scored[0][1]


def load_cycle_series(
    csv_path: Path,
    metadata_path: Path | None,
    theta_column: str | None,
    small_temp_column: str | None,
    large_temp_column: str | None,
    print_columns: bool = False,
):
    columns = load_cycle_csv(csv_path)

    if print_columns:
        print("Detected CSV columns:")
        for name, arr in columns.items():
            finite = np.isfinite(arr)
            n = int(np.sum(finite))
            print(f"  {name}  [n_valid={n}]")
        if metadata_path and metadata_path.exists():
            print(f"\nMetadata file: {metadata_path}")
        return None

    if theta_column is None:
        theta_column = pick_column(columns, [
            "theta_rad", "theta", "angle_rad", "crank_angle", "cycle_angle", "angle"
        ])
    if small_temp_column is None:
        small_temp_column = pick_column(columns, [
            "small_temperature", "small_gas_temperature", "small_temp", "temperature_small",
            "small_cylinder_temperature", "tsmall", "small_t"
        ])
    if large_temp_column is None:
        large_temp_column = pick_column(columns, [
            "large_temperature", "large_gas_temperature", "large_temp", "temperature_large",
            "large_cylinder_temperature", "tlarge", "large_t"
        ])

    if theta_column is None or small_temp_column is None or large_temp_column is None:
        names = ", ".join(columns.keys())
        raise ValueError(
            "Could not determine CSV columns automatically.\n"
            "Run with --print-cycle-columns and then set:\n"
            "  --theta-column ... --small-temp-column ... --large-temp-column ...\n"
            f"Available columns: {names}"
        )

    theta = columns[theta_column]
    small_t = columns[small_temp_column]
    large_t = columns[large_temp_column]

    mask = np.isfinite(theta) & np.isfinite(small_t) & np.isfinite(large_t)
    theta = theta[mask]
    small_t = small_t[mask]
    large_t = large_t[mask]

    if len(theta) >= 2 and abs(theta[-1] - theta[0] - 2.0*math.pi) < 1e-6:
        theta = theta[:-1]
        small_t = small_t[:-1]
        large_t = large_t[:-1]

    return {
        "theta": theta,
        "small_temp": small_t,
        "large_temp": large_t,
        "theta_column": theta_column,
        "small_temp_column": small_temp_column,
        "large_temp_column": large_temp_column,
    }


def resample_periodic(theta_ref, values_ref, theta_query):
    theta_ref = np.asarray(theta_ref, dtype=float)
    values_ref = np.asarray(values_ref, dtype=float)
    theta_query = np.asarray(theta_query, dtype=float)
    order = np.argsort(theta_ref)
    theta_ref = theta_ref[order]
    values_ref = values_ref[order]
    theta_ext = np.concatenate([theta_ref, [theta_ref[0] + 2.0*math.pi]])
    values_ext = np.concatenate([values_ref, [values_ref[0]]])
    tq = np.mod(theta_query, 2.0*math.pi)
    return np.interp(tq, theta_ext, values_ext)


def temp_to_rgb(temp_k: float, tmin: float, tmax: float):
    alpha = 0.5 if tmax <= tmin else float((temp_k - tmin) / (tmax - tmin))
    alpha = max(0.0, min(1.0, alpha))
    return (alpha, 0.0, 1.0 - alpha)


def draw_ground_symbol(ax, p, scale, lw=2.0):
    x, y = p
    ax.plot([x, x], [y - 0.08*scale, y - 0.22*scale], color=COLORS["ground"], lw=lw, solid_capstyle="round")
    base_y = y - 0.24*scale
    half = 0.18*scale
    ax.plot([x - half, x + half], [base_y, base_y], color=COLORS["ground"], lw=lw, solid_capstyle="round")
    for k in range(5):
        xi = x - half + k*(2*half/4)
        ax.plot([xi - 0.04*scale, xi + 0.04*scale], [base_y - 0.07*scale, base_y], color=COLORS["ground"], lw=lw*0.7, solid_capstyle="round")


def draw_joint(ax, p, r=0.075, lw=1.2):
    ax.add_patch(Circle((p[0], p[1]), r, facecolor=COLORS["joint_fill"], edgecolor="black", lw=lw))


def _polyline_collection(points, values, tmin, tmax, lw=3.0, z=3):
    segs = np.stack([points[:-1], points[1:]], axis=1)
    mids = 0.5 * (values[:-1] + values[1:])
    colors = [temp_to_rgb(v, tmin, tmax) for v in mids]
    return LineCollection(segs, colors=colors, linewidths=lw, capstyle="round", joinstyle="round", zorder=z)


def spring_points(x0, x1, y, amplitude, turns, n=200):
    s = np.linspace(0.0, 1.0, n)
    x = x0 + (x1 - x0) * s
    yv = y + amplitude * np.sin(2.0*math.pi*turns*s)
    return np.column_stack([x, yv])


def draw_check_valve(ax, center, size, direction, lw=3.0):
    x, y = center
    half_h = 0.18 * size
    tri_len = 0.32 * size
    gap = 0.06 * size
    bar_h = 0.42 * size

    if direction == "right":
        tri = np.array([[x - tri_len/2, y - half_h], [x - tri_len/2, y + half_h], [x + tri_len/2 - gap, y]])
        bar_x = x + tri_len/2 + gap
    else:
        tri = np.array([[x + tri_len/2, y - half_h], [x + tri_len/2, y + half_h], [x - tri_len/2 + gap, y]])
        bar_x = x - tri_len/2 - gap

    ax.add_patch(Polygon(tri, closed=True, fill=False, edgecolor="black", lw=lw, joinstyle="round"))
    ax.plot([bar_x, bar_x], [y - bar_h/2, y + bar_h/2], color="black", lw=lw, solid_capstyle="round")


def draw_transfer_block(ax, x_left, x_right, large_width, t_small, t_large, tmin, tmax, lw_pipe=2.6):
    y_top = 1.15 * large_width
    y_bot = -1.15 * large_width
    xl = x_left + 0.06 * (x_right - x_left)
    xr = x_right - 0.06 * (x_right - x_left)
    span = xr - xl

    for x_start, x_end, yy in [(x_left, xl, 0.0), (x_right, xr, 0.0)]:
        ax.plot([x_start, x_end], [yy, yy], color=COLORS["pipe"], lw=lw_pipe, solid_capstyle="round")
    ax.plot([xl, xl], [y_bot, y_top], color=COLORS["pipe"], lw=lw_pipe, solid_capstyle="round")
    ax.plot([xr, xr], [y_bot, y_top], color=COLORS["pipe"], lw=lw_pipe, solid_capstyle="round")

    vt = xl + 0.20 * span
    et0 = xl + 0.42 * span
    et1 = xl + 0.80 * span
    ax.plot([xl, vt - 0.26*large_width], [y_top, y_top], color=COLORS["pipe"], lw=lw_pipe, solid_capstyle="round")
    draw_check_valve(ax, (vt, y_top), 0.72*large_width, "right", lw=lw_pipe*0.9)
    ax.plot([vt + 0.20*large_width, et0], [y_top, y_top], color=COLORS["pipe"], lw=lw_pipe, solid_capstyle="round")
    pts_top = spring_points(et0, et1, y_top, 0.32*large_width, 6)
    ax.plot(pts_top[:, 0], pts_top[:, 1], color=COLORS["pipe"], lw=lw_pipe, solid_capstyle="round")
    ax.plot([et1, xr], [y_top, y_top], color=COLORS["pipe"], lw=lw_pipe, solid_capstyle="round")

    vb = xr - 0.20 * span
    eb0 = xl + 0.20 * span
    eb1 = xl + 0.58 * span
    ax.plot([xl, eb0], [y_bot, y_bot], color=COLORS["pipe"], lw=lw_pipe, solid_capstyle="round")
    pts_bot = spring_points(eb0, eb1, y_bot, 0.32*large_width, 6)
    ax.plot(pts_bot[:, 0], pts_bot[:, 1], color=COLORS["pipe"], lw=lw_pipe, solid_capstyle="round")
    ax.plot([eb1, vb - 0.20*large_width], [y_bot, y_bot], color=COLORS["pipe"], lw=lw_pipe, solid_capstyle="round")
    draw_check_valve(ax, (vb, y_bot), 0.72*large_width, "left", lw=lw_pipe*0.9)
    ax.plot([vb + 0.26*large_width, xr], [y_bot, y_bot], color=COLORS["pipe"], lw=lw_pipe, solid_capstyle="round")

    overlay_lw = 1.35
    cL = temp_to_rgb(t_small, tmin, tmax)
    cR = temp_to_rgb(t_large, tmin, tmax)

    ax.plot([x_left, xl], [0, 0], color=cL, lw=overlay_lw, solid_capstyle="round")
    ax.plot([x_right, xr], [0, 0], color=cR, lw=overlay_lw, solid_capstyle="round")
    ax.plot([xl, xl], [0, y_top], color=cL, lw=overlay_lw, solid_capstyle="round")
    ax.plot([xr, xr], [0, y_top], color=cR, lw=overlay_lw, solid_capstyle="round")
    ax.plot([xl, xl], [0, y_bot], color=cL, lw=overlay_lw, solid_capstyle="round")
    ax.plot([xr, xr], [0, y_bot], color=cR, lw=overlay_lw, solid_capstyle="round")

    def line_points(a, b, y):
        return np.array([[a, y], [b, y]], dtype=float)

    ax.add_collection(_polyline_collection(line_points(xl, et0, y_top), np.array([t_small, 0.5*(t_small+t_large)]), tmin, tmax, lw=overlay_lw, z=4))
    ax.add_collection(_polyline_collection(pts_top, np.linspace(t_small, t_large, len(pts_top)), tmin, tmax, lw=overlay_lw, z=4))
    ax.add_collection(_polyline_collection(line_points(et1, xr, y_top), np.array([0.5*(t_small+t_large), t_large]), tmin, tmax, lw=overlay_lw, z=4))

    ax.add_collection(_polyline_collection(line_points(xl, eb0, y_bot), np.array([t_small, 0.5*(t_small+t_large)]), tmin, tmax, lw=overlay_lw, z=4))
    ax.add_collection(_polyline_collection(pts_bot, np.linspace(t_small, t_large, len(pts_bot)), tmin, tmax, lw=overlay_lw, z=4))
    ax.add_collection(_polyline_collection(line_points(eb1, xr, y_bot), np.array([0.5*(t_small+t_large), t_large]), tmin, tmax, lw=overlay_lw, z=4))


def draw_cylinder_U(ax, cyl, piston_x, gas_color, lw_u=7.0, lw_piston=9.5):
    side = cyl["side"]
    w = cyl["width"]
    x_head = cyl["x_inner"]
    x_open = cyl["x_outer"]

    x0 = min(x_head, piston_x)
    x1 = max(x_head, piston_x)
    ax.add_patch(Rectangle((x0, -w/2.0), x1 - x0, w, facecolor=gas_color, edgecolor="none", alpha=0.85, zorder=1))

    if side == "left":
        ax.plot([x_head, x_head], [-w/2.0, w/2.0], color="black", lw=lw_u, solid_capstyle="round")
        ax.plot([x_open, x_head], [w/2.0, w/2.0], color="black", lw=lw_u, solid_capstyle="round")
        ax.plot([x_open, x_head], [-w/2.0, -w/2.0], color="black", lw=lw_u, solid_capstyle="round")
    else:
        ax.plot([x_head, x_head], [-w/2.0, w/2.0], color="black", lw=lw_u, solid_capstyle="round")
        ax.plot([x_head, x_open], [w/2.0, w/2.0], color="black", lw=lw_u, solid_capstyle="round")
        ax.plot([x_head, x_open], [-w/2.0, -w/2.0], color="black", lw=lw_u, solid_capstyle="round")

    ax.plot([piston_x, piston_x], [-w/2.0, w/2.0], color=COLORS["piston"], lw=lw_piston, solid_capstyle="butt", zorder=5)


def draw_mechanism(ax, s, joint_r=0.080, lw=2.6, ground_scale=0.80):
    A, B, C, D, E = s["A"], s["B"], s["C"], s["D"], s["E"]
    G, F, H, P = s["G"], s["F"], s["H"], s["P"]

    ax.plot([A[0], B[0]], [A[1], B[1]], color=COLORS["crank"], lw=lw, solid_capstyle="round")
    ax.plot([B[0], C[0]], [B[1], C[1]], color=COLORS["coupler"], lw=lw, solid_capstyle="round")
    ax.plot([C[0], D[0]], [C[1], D[1]], color=COLORS["rocker"], lw=lw, solid_capstyle="round")
    ax.plot([B[0], E[0]], [B[1], E[1]], color=COLORS["coupler"], lw=lw*0.9, solid_capstyle="round")
    ax.plot([C[0], E[0]], [C[1], E[1]], color=COLORS["coupler"], lw=lw*0.9, solid_capstyle="round")

    ax.plot([G[0], F[0]], [G[1], F[1]], color=COLORS["secondary_input"], lw=lw, solid_capstyle="round")
    ax.plot([E[0], F[0]], [E[1], F[1]], color=COLORS["secondary_plate"], lw=lw, solid_capstyle="round")
    ax.plot([E[0], H[0]], [E[1], H[1]], color=COLORS["secondary_plate"], lw=lw*0.95, solid_capstyle="round")
    ax.plot([F[0], H[0]], [F[1], H[1]], color=COLORS["secondary_plate"], lw=lw*0.95, solid_capstyle="round")

    ax.plot([H[0], P[0]], [H[1], P[1]], color=COLORS["rod"], lw=lw, solid_capstyle="round")

    draw_ground_symbol(ax, A, ground_scale, lw=lw*0.75)
    draw_ground_symbol(ax, D, ground_scale, lw=lw*0.75)
    draw_ground_symbol(ax, G, ground_scale, lw=lw*0.75)
    for pt in (A, B, C, D, E, F, G, H):
        draw_joint(ax, pt, r=joint_r, lw=lw*0.40)


def collect_cloud(samples):
    pts = []
    for s in samples:
        for key in ("A", "B", "C", "D", "E", "F", "G", "H", "P"):
            pts.append(s[key])
    return np.array(pts)


def make_animation(
    cycle_csv: Path,
    cycle_metadata_json: Path | None,
    small_json: Path,
    large_json: Path,
    output_gif: Path,
    small_restart: int | None,
    small_branch: int | None,
    large_restart: int | None,
    large_branch: int | None,
    theta_column: str | None,
    small_temp_column: str | None,
    large_temp_column: str | None,
    print_cycle_columns: bool,
    frames: int,
    fps: int,
    dpi: int,
    inter_cylinder_gap: float,
    large_cylinder_width: float,
    small_width_ratio: float,
):
    cycle = load_cycle_series(
        cycle_csv,
        metadata_path=cycle_metadata_json,
        theta_column=theta_column,
        small_temp_column=small_temp_column,
        large_temp_column=large_temp_column,
        print_columns=print_cycle_columns,
    )
    if cycle is None:
        return

    small_best = load_candidate(small_json, small_restart, small_branch)
    large_best = load_candidate(large_json, large_restart, large_branch)

    small_samples_local = build_aligned_samples(small_best, frames, "left")
    large_samples_local = build_aligned_samples(large_best, frames, "right")

    frame_thetas = np.linspace(0.0, 2.0*math.pi, frames, endpoint=False)
    small_temps = resample_periodic(cycle["theta"], cycle["small_temp"], frame_thetas)
    large_temps = resample_periodic(cycle["theta"], cycle["large_temp"], frame_thetas)

    temp_min = float(min(np.min(small_temps), np.min(large_temps)))
    temp_max = float(max(np.max(small_temps), np.max(large_temps)))

    small_width = small_width_ratio * large_cylinder_width
    x_left_head = -inter_cylinder_gap / 2.0
    x_right_head = inter_cylinder_gap / 2.0

    small_samples, small_cyl = place_samples(small_samples_local, "left", x_left_head, small_width)
    large_samples, large_cyl = place_samples(large_samples_local, "right", x_right_head, large_cylinder_width)

    cloud = np.vstack((collect_cloud(small_samples), collect_cloud(large_samples)))
    xmin = min(float(np.min(cloud[:, 0])), small_cyl["x_outer"], small_cyl["x_inner"], large_cyl["x_outer"], large_cyl["x_inner"])
    xmax = max(float(np.max(cloud[:, 0])), small_cyl["x_outer"], small_cyl["x_inner"], large_cyl["x_outer"], large_cyl["x_inner"])
    ymin = min(float(np.min(cloud[:, 1])), -1.45*large_cylinder_width)
    ymax = max(float(np.max(cloud[:, 1])), 1.45*large_cylinder_width)

    spanx = xmax - xmin
    spany = ymax - ymin
    span = max(spanx, spany)
    margin_x = 0.012 * span + 0.05
    margin_y = 0.020 * span + 0.10

    fig, ax = plt.subplots(figsize=(12.0, 5.8))  # 1200 px wide at dpi=100
    ax.set_aspect("equal")

    def update(i):
        ax.clear()
        ax.set_aspect("equal")
        ax.set_xlim(xmin - margin_x, xmax + margin_x)
        ax.set_ylim(ymin - margin_y, ymax + margin_y)
        ax.axis("off")

        sS = small_samples[i]
        sL = large_samples[i]
        cS = temp_to_rgb(float(small_temps[i]), temp_min, temp_max)
        cL = temp_to_rgb(float(large_temps[i]), temp_min, temp_max)

        draw_transfer_block(
            ax, small_cyl["x_inner"], large_cyl["x_inner"], large_cylinder_width,
            t_small=float(small_temps[i]), t_large=float(large_temps[i]),
            tmin=temp_min, tmax=temp_max, lw_pipe=2.5,
        )
        draw_cylinder_U(ax, small_cyl, sS["P"][0], gas_color=cS, lw_u=7.0, lw_piston=9.5)
        draw_cylinder_U(ax, large_cyl, sL["P"][0], gas_color=cL, lw_u=7.2, lw_piston=9.8)
        draw_mechanism(ax, sS, joint_r=0.070, lw=2.45, ground_scale=0.44*small_width)
        draw_mechanism(ax, sL, joint_r=0.082, lw=2.55, ground_scale=0.44*large_cylinder_width)
        return ()

    anim = FuncAnimation(fig, update, frames=frames, interval=1000/fps, blit=False, repeat=True)
    output_gif.parent.mkdir(parents=True, exist_ok=True)
    anim.save(output_gif, writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cycle_csv", type=Path)
    ap.add_argument("cycle_metadata_json", type=Path, nargs="?", default=None)
    ap.add_argument("small_json", type=Path)
    ap.add_argument("large_json", type=Path)

    ap.add_argument("--small-restart", type=int, default=None)
    ap.add_argument("--small-branch", type=int, choices=[-1, 1], default=None)
    ap.add_argument("--large-restart", type=int, default=None)
    ap.add_argument("--large-branch", type=int, choices=[-1, 1], default=None)

    ap.add_argument("--theta-column", type=str, default=None)
    ap.add_argument("--small-temp-column", type=str, default=None)
    ap.add_argument("--large-temp-column", type=str, default=None)
    ap.add_argument("--print-cycle-columns", action="store_true")

    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--frames", type=int, default=144)
    ap.add_argument("--fps", type=int, default=20)
    ap.add_argument("--dpi", type=int, default=100)

    ap.add_argument("--inter-cylinder-gap", type=float, default=9.0)
    ap.add_argument("--large-cylinder-width", type=float, default=8.0)
    ap.add_argument("--small-width-ratio", type=float, default=0.8)

    args = ap.parse_args()

    output = args.output if args.output else Path("outputs/sixbar_pair_v3.gif")
    make_animation(
        cycle_csv=args.cycle_csv,
        cycle_metadata_json=args.cycle_metadata_json,
        small_json=args.small_json,
        large_json=args.large_json,
        output_gif=output,
        small_restart=args.small_restart,
        small_branch=args.small_branch,
        large_restart=args.large_restart,
        large_branch=args.large_branch,
        theta_column=args.theta_column,
        small_temp_column=args.small_temp_column,
        large_temp_column=args.large_temp_column,
        print_cycle_columns=args.print_cycle_columns,
        frames=args.frames,
        fps=args.fps,
        dpi=args.dpi,
        inter_cylinder_gap=args.inter_cylinder_gap,
        large_cylinder_width=args.large_cylinder_width,
        small_width_ratio=args.small_width_ratio,
    )
    if not args.print_cycle_columns:
        print(f"Wrote {output}")


if __name__ == "__main__":
    main()
