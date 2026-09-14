#!/usr/bin/env python3
"""DADA six-bar pair animation, V3.8.

Changes vs V3.6
---------------
- valve positions are kept fixed;
- each exchanger keeps the same dimensions as V3.6 (same length) and the same
  height as V3.4, and is centered between the fixed valve and the opposite
  cylinder;
- conduit / exchanger thickness is reduced by 3 px;
- ports are extended farther into both cylinders, by one extra cylinder wall
  thickness, so the coloured conduit clearly overlaps the cylinder stroke line.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle, Polygon, Rectangle
import numpy as np


ROOT = Path.cwd()

DEFAULT_MODEL = ROOT / "outputs" / "motor_champion_sixbar_k2.json"
DEFAULT_CYCLE = ROOT / "outputs" / "champion-6bar_thermodynamic_cycle.csv"
DEFAULT_METADATA = ROOT / "outputs" / "champion-6bar_thermodynamic_cycle_metadata.json"
DEFAULT_OUTPUT = ROOT / "outputs" / "sixbar_pair_v3_8.gif"


COLORS = {
    "crank": "#1f77b4",
    "coupler": "#ff7f0e",
    "rocker": "#2ca02c",
    "secondary_dyad": "#9467bd",
    "secondary_plate": "#d62728",
    "rod": "#444444",
    "ground": "#111111",
    "joint_fill": "#ffffff",
    "piston": "#555555",
    "cylinder": "#111111",
    "valve_open": "#00ff00",
    "valve_closed": "#000000",
}


def load_mechanisms(model_path: Path) -> tuple[dict, dict]:
    data = json.loads(model_path.read_text(encoding="utf-8"))
    return data["kinematics"]["small"], data["kinematics"]["large"]


def primary_state(theta: float, p: dict) -> dict[str, np.ndarray]:
    g = float(p["primary_ground"])
    c = float(p["primary_coupler"])
    r = float(p["primary_rocker"])
    branch = int(p["primary_branch"])
    phase = float(p["primary_phase"])

    a = theta + phase
    A = np.array((0.0, 0.0))
    B = np.array((math.cos(a), math.sin(a)))
    D = np.array((g, 0.0))

    delta = D - B
    dist = float(np.linalg.norm(delta))
    u = delta / dist
    along = (c*c - r*r + dist*dist) / (2.0*dist)
    h2 = c*c - along*along
    n = np.array((-u[1], u[0]))
    C = B + along*u + branch*math.sqrt(max(1e-12, h2))*n

    ubc = (C - B) / c
    nbc = np.array((-ubc[1], ubc[0]))
    E = B + float(p["primary_e_along"]) * ubc + float(p["primary_e_normal"]) * nbc
    return {"A": A, "B": B, "C": C, "D": D, "E": E}


def sixbar_state(theta: float, p: dict) -> dict[str, np.ndarray]:
    s = primary_state(theta, p)
    E = s["E"]

    G = np.array((float(p["second_pivot_x"]), float(p["second_pivot_y"])))
    lef = float(p["link_ef"])
    lgf = float(p["link_gf"])
    branch = int(p["second_branch"])

    delta = G - E
    dist = float(np.linalg.norm(delta))
    u = delta / dist
    along = (lef*lef - lgf*lgf + dist*dist) / (2.0*dist)
    h2 = lef*lef - along*along
    n = np.array((-u[1], u[0]))
    F = E + along*u + branch*math.sqrt(max(1e-12, h2))*n

    uef = (F - E) / lef
    nef = np.array((-uef[1], uef[0]))
    H = (
        E
        + float(p["h_along_over_ef"]) * lef * uef
        + float(p["h_normal_over_ef"]) * lef * nef
    )

    rod = float(p["piston_rod"])
    axis_angle = float(p["slider_axis_angle"])
    axis_offset = float(p["slider_axis_offset"])

    axis = np.array((math.cos(axis_angle), math.sin(axis_angle)))
    normal = np.array((-axis[1], axis[0]))
    origin = axis_offset * normal

    rel = H - origin
    longitudinal = float(rel @ axis)
    transverse = float(rel @ normal)

    slider = longitudinal + math.sqrt(max(1e-12, rod*rod - transverse*transverse))
    P = origin + slider*axis

    s.update({
        "F": F,
        "G": G,
        "H": H,
        "P": P,
        "axis_origin": origin,
        "axis_angle": axis_angle,
    })
    return s


POINT_KEYS = ("A", "B", "C", "D", "E", "F", "G", "H", "P", "axis_origin")


def rotate_point(p: np.ndarray, angle: float) -> np.ndarray:
    c = math.cos(angle)
    s = math.sin(angle)
    return np.array((c*p[0] - s*p[1], s*p[0] + c*p[1]))


def orient_state(state: dict, side: str) -> dict:
    target_axis = 0.0 if side == "left" else math.pi
    rot = target_axis - float(state["axis_angle"])

    out = {}
    for key in POINT_KEYS:
        out[key] = rotate_point(state[key], rot)

    dy = -out["axis_origin"][1]
    for key in POINT_KEYS:
        out[key] = out[key] + np.array((0.0, dy))
    return out


def mechanism_samples(mech: dict, study_angles_rad: np.ndarray, side: str) -> list[dict]:
    return [orient_state(sixbar_state(float(theta), mech), side) for theta in study_angles_rad]


def place_samples(samples: list[dict], side: str, inner_head_x: float, cylinder_width: float):
    piston_x = np.asarray([s["P"][0] for s in samples])
    pmin = float(np.min(piston_x))
    pmax = float(np.max(piston_x))

    head_clear = 0.065 * cylinder_width

    if side == "left":
        shift_x = inner_head_x - head_clear - pmax
    else:
        shift_x = inner_head_x + head_clear - pmin

    placed = []
    for s in samples:
        q = {}
        for k, v in s.items():
            q[k] = v + np.array((shift_x, 0.0)) if isinstance(v, np.ndarray) else v
        placed.append(q)

    px = np.asarray([s["P"][0] for s in placed])

    rod_side_extra = 0.50
    if side == "left":
        outer_x = float(np.min(px)) - rod_side_extra
    else:
        outer_x = float(np.max(px)) + rod_side_extra

    return placed, {
        "side": side,
        "x_inner": inner_head_x,
        "x_outer": outer_x,
        "width": cylinder_width,
    }


REQUIRED_COLUMNS = (
    "motor_angle_deg",
    "study_angle_deg",
    "S_temperature_K",
    "L_temperature_K",
    "Hi_temperature_K",
    "Ho_temperature_K",
    "Hi_to_L_valve_open",
    "Ho_to_S_valve_open",
)


def load_cycle(path: Path) -> dict[str, np.ndarray]:
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    out = {c: np.asarray([float(row[c]) for row in rows], dtype=float) for c in REQUIRED_COLUMNS}
    if len(out["motor_angle_deg"]) >= 2 and abs(out["motor_angle_deg"][-1] - out["motor_angle_deg"][0] - 360.0) < 1e-9:
        for k in out:
            out[k] = out[k][:-1]
    return out


def periodic_interp(motor_deg: np.ndarray, values: np.ndarray, query_motor_deg: np.ndarray) -> np.ndarray:
    x = np.asarray(motor_deg, dtype=float)
    y = np.asarray(values, dtype=float)
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    xext = np.concatenate([x, [x[0] + 360.0]])
    yext = np.concatenate([y, [y[0]]])
    q = np.mod(query_motor_deg, 360.0)
    return np.interp(q, xext, yext)


def load_metadata(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def temperature_rgb(T: float, Tmin: float, Tmax: float) -> tuple[float, float, float]:
    a = 0.5 if Tmax <= Tmin else (float(T) - Tmin) / (Tmax - Tmin)
    a = min(1.0, max(0.0, a))
    return (a, 0.0, 1.0 - a)


def colored_polyline(ax, pts: np.ndarray, temps: np.ndarray, Tmin: float, Tmax: float, linewidth=3.0, zorder=3, capstyle="butt", joinstyle="miter"):
    segs = np.stack([pts[:-1], pts[1:]], axis=1)
    mids = 0.5 * (temps[:-1] + temps[1:])
    lc = LineCollection(
        segs,
        colors=[temperature_rgb(t, Tmin, Tmax) for t in mids],
        linewidths=linewidth,
        capstyle=capstyle,
        joinstyle=joinstyle,
        zorder=zorder,
    )
    ax.add_collection(lc)


def spring_points(x0, x1, y, amplitude, turns=4.5, n=180):
    s = np.linspace(0.0, 1.0, n)
    return np.column_stack((x0 + (x1 - x0) * s, y + amplitude * np.sin(2 * math.pi * turns * s)))


def data_dx_for_points(ax, points: float) -> float:
    pixels = points * ax.figure.dpi / 72.0
    x0_disp, y0_disp = ax.transData.transform((0.0, 0.0))
    x1_data = ax.transData.inverted().transform((x0_disp + pixels, y0_disp))[0]
    return float(x1_data)


def draw_ground_symbol(ax, p, scale, lw=2.0):
    x, y = p
    ax.plot([x, x], [y - 0.08*scale, y - 0.22*scale], color=COLORS["ground"], lw=lw, solid_capstyle="round")
    base_y = y - 0.24*scale
    half = 0.18*scale
    ax.plot([x-half, x+half], [base_y, base_y], color=COLORS["ground"], lw=lw, solid_capstyle="round")
    for k in range(5):
        xi = x - half + k * (2*half/4)
        ax.plot([xi - 0.04*scale, xi + 0.04*scale], [base_y - 0.07*scale, base_y], color=COLORS["ground"], lw=lw*0.7, solid_capstyle="round")


def draw_joint(ax, p, r=0.075, lw=1.1):
    ax.add_patch(Circle((p[0], p[1]), r, facecolor=COLORS["joint_fill"], edgecolor="black", lw=lw, zorder=8))


def draw_mechanism(ax, s, width, lw=2.5):
    A, B, C, D, E = s["A"], s["B"], s["C"], s["D"], s["E"]
    F, G, H, P = s["F"], s["G"], s["H"], s["P"]

    crank_radius = float(np.linalg.norm(B - A))
    ax.add_patch(Circle((A[0], A[1]), crank_radius, fill=False, edgecolor="black", lw=0.85, zorder=0))

    ax.plot(*zip(A, B), color=COLORS["crank"], lw=lw, solid_capstyle="round")
    ax.plot(*zip(B, C), color=COLORS["coupler"], lw=lw, solid_capstyle="round")
    ax.plot(*zip(C, D), color=COLORS["rocker"], lw=lw, solid_capstyle="round")
    ax.plot(*zip(B, E), color=COLORS["coupler"], lw=lw*0.92, solid_capstyle="round")
    ax.plot(*zip(C, E), color=COLORS["coupler"], lw=lw*0.92, solid_capstyle="round")

    ax.plot(*zip(G, F), color=COLORS["secondary_dyad"], lw=lw, solid_capstyle="round")
    ax.plot(*zip(E, F), color=COLORS["secondary_plate"], lw=lw, solid_capstyle="round")
    ax.plot(*zip(E, H), color=COLORS["secondary_plate"], lw=lw, solid_capstyle="round")
    ax.plot(*zip(F, H), color=COLORS["secondary_plate"], lw=lw, solid_capstyle="round")
    ax.plot(*zip(H, P), color=COLORS["rod"], lw=lw, solid_capstyle="round")

    gs = 0.40 * width
    draw_ground_symbol(ax, A, gs, lw=lw*0.72)
    draw_ground_symbol(ax, D, gs, lw=lw*0.72)
    draw_ground_symbol(ax, G, gs, lw=lw*0.72)

    jr = 0.075 if width < 7 else 0.085
    for pt in (A, B, C, D, E, F, G, H):
        draw_joint(ax, pt, r=jr, lw=lw*0.42)


def draw_cylinder(ax, cyl, piston_x, gas_T, Tmin, Tmax, wall_lw=7.0, piston_lw=10.0):
    w = float(cyl["width"])
    head = float(cyl["x_inner"])
    outer = float(cyl["x_outer"])

    x0, x1 = sorted((head, float(piston_x)))
    if x1 > x0:
        ax.add_patch(Rectangle((x0, -w/2), x1-x0, w, facecolor=temperature_rgb(gas_T, Tmin, Tmax), edgecolor="none", alpha=0.92, zorder=1))

    ax.plot([piston_x, piston_x], [-w/2, w/2], color=COLORS["piston"], lw=piston_lw, solid_capstyle="butt", zorder=4)
    ax.plot([head, head], [-w/2, w/2], color=COLORS["cylinder"], lw=wall_lw, solid_capstyle="round", zorder=6)
    ax.plot([min(head, outer), max(head, outer)], [w/2, w/2], color=COLORS["cylinder"], lw=wall_lw, solid_capstyle="round", zorder=6)
    ax.plot([min(head, outer), max(head, outer)], [-w/2, -w/2], color=COLORS["cylinder"], lw=wall_lw, solid_capstyle="round", zorder=6)


def draw_check_valve(ax, center, size, direction, is_open):
    x, y = map(float, center)
    h = 0.18 * size
    length = 0.44 * size
    seat_h = 0.54 * size
    gap = (0.10 if is_open else 0.03) * size
    fill = COLORS["valve_open"] if is_open else COLORS["valve_closed"]

    tri_lw = 1.7
    bar_lw = 6.8

    if direction == "right":
        seat_x = x + length/2
        tip_x = seat_x - gap
        base_x = x - length/2
    else:
        seat_x = x - length/2
        tip_x = seat_x + gap
        base_x = x + length/2

    tri = Polygon(
        [(base_x, y-h), (base_x, y+h), (tip_x, y)],
        closed=True,
        facecolor=fill,
        edgecolor="black",
        lw=tri_lw,
        joinstyle="round",
        zorder=9,
    )
    ax.add_patch(tri)
    ax.plot([seat_x, seat_x], [y-seat_h/2, y+seat_h/2], color="black", lw=bar_lw, solid_capstyle="round", zorder=9)


def draw_transfer_block(
    ax,
    left_head,
    right_head,
    small_width,
    large_width,
    Ts,
    Tl,
    Thi,
    Tho,
    hi_valve_open,
    ho_valve_open,
    Tmin,
    Tmax,
    small_wall_lw,
    large_wall_lw,
):
    xL = float(left_head)
    xR = float(right_head)
    total_inside = xR - xL

    y_top = 0.31 * large_width
    y_bot = -0.31 * large_width

    path_lw = 3.0              # reduced by 3 px
    coil_amp = 0.18 * large_width  # same as V3.4
    coil_len = 0.36 * total_inside  # same as V3.5 / V3.6

    # Extend farther into cylinders on the CORRECT side of each head wall.
    # Small cylinder head is on the right side of its chamber -> go left into gas.
    # Large cylinder head is on the left side of its chamber -> go right into gas.
    left_inside = xL - data_dx_for_points(ax, 1.5 * small_wall_lw)
    right_inside = xR + data_dx_for_points(ax, 1.5 * large_wall_lw)

    # Keep valves fixed.
    top_valve_x = xL + 0.22 * total_inside
    bot_valve_x = xL + 0.78 * total_inside

    # Center exchanger between fixed valve and opposite cylinder.
    top_center = 0.5 * (top_valve_x + right_inside)
    top_coil0_x = top_center - 0.5 * coil_len
    top_coil1_x = top_center + 0.5 * coil_len
    top_spring = spring_points(top_coil0_x, top_coil1_x, y_top, coil_amp, turns=4.5)

    bot_center = 0.5 * (left_inside + bot_valve_x)
    bot_coil0_x = bot_center - 0.5 * coil_len
    bot_coil1_x = bot_center + 0.5 * coil_len
    bot_spring = spring_points(bot_coil0_x, bot_coil1_x, y_bot, coil_amp, turns=4.5)

    # Top / Ho : L -> exchanger -> valve -> S
    pts_top_r = np.array([[right_inside, y_top], [top_coil1_x, y_top]], dtype=float)
    colored_polyline(ax, pts_top_r, np.array([Tl, Tho]), Tmin, Tmax, linewidth=path_lw, zorder=7, capstyle="butt", joinstyle="miter")
    colored_polyline(ax, top_spring[::-1], np.full(len(top_spring), Tho), Tmin, Tmax, linewidth=path_lw, zorder=7, capstyle="round", joinstyle="round")
    pts_top_l = np.array([[top_coil0_x, y_top], [top_valve_x, y_top], [left_inside, y_top]], dtype=float)
    colored_polyline(ax, pts_top_l, np.array([Tho, Tho, Ts]), Tmin, Tmax, linewidth=path_lw, zorder=7, capstyle="butt", joinstyle="miter")
    draw_check_valve(ax, (top_valve_x, y_top), (0.90 * large_width) * (2.0/3.0), direction="left", is_open=bool(round(ho_valve_open)))

    # Bottom / Hi : S -> exchanger -> valve -> L
    pts_bot_l = np.array([[left_inside, y_bot], [bot_coil0_x, y_bot]], dtype=float)
    colored_polyline(ax, pts_bot_l, np.array([Ts, Thi]), Tmin, Tmax, linewidth=path_lw, zorder=7, capstyle="butt", joinstyle="miter")
    colored_polyline(ax, bot_spring, np.full(len(bot_spring), Thi), Tmin, Tmax, linewidth=path_lw, zorder=7, capstyle="round", joinstyle="round")
    pts_bot_r = np.array([[bot_coil1_x, y_bot], [bot_valve_x, y_bot], [right_inside, y_bot]], dtype=float)
    colored_polyline(ax, pts_bot_r, np.array([Thi, Thi, Tl]), Tmin, Tmax, linewidth=path_lw, zorder=7, capstyle="butt", joinstyle="miter")
    draw_check_valve(ax, (bot_valve_x, y_bot), (0.90 * large_width) * (2.0/3.0), direction="right", is_open=bool(round(hi_valve_open)))


def cloud(samples):
    pts = []
    for s in samples:
        for key in ("A", "B", "C", "D", "E", "F", "G", "H", "P"):
            pts.append(s[key])
    return np.asarray(pts)


def make_animation(model_json: Path, cycle_csv: Path, metadata_json: Path, output: Path, frames: int, fps: int, dpi: int, large_width: float, small_width_ratio: float, exchanger_gap: float):
    small_mech, large_mech = load_mechanisms(model_json)
    cycle = load_cycle(cycle_csv)
    _ = load_metadata(metadata_json)

    motor_query = np.linspace(0.0, 360.0, frames, endpoint=False)
    study_deg = periodic_interp(cycle["motor_angle_deg"], cycle["study_angle_deg"], motor_query)
    study_rad = np.deg2rad(study_deg)

    Ts = periodic_interp(cycle["motor_angle_deg"], cycle["S_temperature_K"], motor_query)
    Tl = periodic_interp(cycle["motor_angle_deg"], cycle["L_temperature_K"], motor_query)
    Thi = periodic_interp(cycle["motor_angle_deg"], cycle["Hi_temperature_K"], motor_query)
    Tho = periodic_interp(cycle["motor_angle_deg"], cycle["Ho_temperature_K"], motor_query)
    hi_open = periodic_interp(cycle["motor_angle_deg"], cycle["Hi_to_L_valve_open"], motor_query)
    ho_open = periodic_interp(cycle["motor_angle_deg"], cycle["Ho_to_S_valve_open"], motor_query)

    all_internal_T = np.concatenate([cycle["S_temperature_K"], cycle["L_temperature_K"], cycle["Hi_temperature_K"], cycle["Ho_temperature_K"]])
    Tmin = float(np.min(all_internal_T))
    Tmax = float(np.max(all_internal_T))

    small_width = small_width_ratio * large_width
    small_local = mechanism_samples(small_mech, study_rad, side="left")
    large_local = mechanism_samples(large_mech, study_rad, side="right")

    left_head = -exchanger_gap / 2
    right_head = exchanger_gap / 2

    small_samples, small_cyl = place_samples(small_local, "left", left_head, small_width)
    large_samples, large_cyl = place_samples(large_local, "right", right_head, large_width)

    pts = np.vstack([cloud(small_samples), cloud(large_samples)])

    xmin = min(float(np.min(pts[:, 0])), small_cyl["x_outer"], small_cyl["x_inner"], large_cyl["x_outer"], large_cyl["x_inner"])
    xmax = max(float(np.max(pts[:, 0])), small_cyl["x_outer"], small_cyl["x_inner"], large_cyl["x_outer"], large_cyl["x_inner"])
    ymin = min(float(np.min(pts[:, 1])), -1.05 * large_width)
    ymax = max(float(np.max(pts[:, 1])), +1.05 * large_width)

    xspan = xmax - xmin
    yspan = ymax - ymin
    mx = 0.003 * xspan + 0.07
    my = 0.020 * yspan + 0.04

    fig = plt.figure(figsize=(12.0, 6.0), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_aspect("equal", adjustable="box")

    small_wall_lw = 6.5
    small_piston_lw = 9.5
    large_wall_lw = 7.0
    large_piston_lw = 10.0

    def update(i):
        ax.clear()
        ax.set_xlim(xmin - mx, xmax + mx)
        ax.set_ylim(ymin - my, ymax + my)
        ax.set_aspect("equal", adjustable="box")
        ax.axis("off")

        sS = small_samples[i]
        sL = large_samples[i]

        draw_cylinder(ax, small_cyl, sS["P"][0], gas_T=float(Ts[i]), Tmin=Tmin, Tmax=Tmax, wall_lw=small_wall_lw, piston_lw=small_piston_lw)
        draw_cylinder(ax, large_cyl, sL["P"][0], gas_T=float(Tl[i]), Tmin=Tmin, Tmax=Tmax, wall_lw=large_wall_lw, piston_lw=large_piston_lw)

        draw_transfer_block(
            ax,
            small_cyl["x_inner"],
            large_cyl["x_inner"],
            small_cyl["width"],
            large_cyl["width"],
            Ts=float(Ts[i]),
            Tl=float(Tl[i]),
            Thi=float(Thi[i]),
            Tho=float(Tho[i]),
            hi_valve_open=float(hi_open[i]),
            ho_valve_open=float(ho_open[i]),
            Tmin=Tmin,
            Tmax=Tmax,
            small_wall_lw=small_wall_lw,
            large_wall_lw=large_wall_lw,
        )

        draw_mechanism(ax, sS, small_width, lw=2.55)
        draw_mechanism(ax, sL, large_width, lw=2.65)
        return ()

    anim = FuncAnimation(fig, update, frames=frames, interval=1000 / fps, blit=False, repeat=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    anim.save(output, writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    print(f"Wrote {output}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    ap.add_argument("--cycle", type=Path, default=DEFAULT_CYCLE)
    ap.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)

    ap.add_argument("--frames", type=int, default=144)
    ap.add_argument("--fps", type=int, default=20)
    ap.add_argument("--dpi", type=int, default=100)

    ap.add_argument("--large-cylinder-width", type=float, default=8.0)
    ap.add_argument("--small-width-ratio", type=float, default=0.8)
    ap.add_argument("--exchanger-gap", type=float, default=9.6)

    args = ap.parse_args()

    make_animation(
        model_json=args.model,
        cycle_csv=args.cycle,
        metadata_json=args.metadata,
        output=args.output,
        frames=args.frames,
        fps=args.fps,
        dpi=args.dpi,
        large_width=args.large_cylinder_width,
        small_width_ratio=args.small_width_ratio,
        exchanger_gap=args.exchanger_gap,
    )


if __name__ == "__main__":
    main()
