#!/usr/bin/env python3
"""DADA six-bar pair animation, V3.1.

Main refinements over V3:
- piston is drawn first, cylinder afterwards, so the piston visually passes under
  the cylinder head;
- a small head clearance is kept so the piston line does not visually cross the
  closed cylinder head;
- no T-manifolds in the exchanger block: just two direct colored transfer lines;
- transfer lines/exchangers are drawn without black outlines, directly in gas
  temperature colors;
- thicker conduits, narrower spring exchangers;
- cylinders moved a bit closer by default;
- tighter margins to make the mechanisms appear larger;
- thin black circles show the crank-pin trajectories.

Run:
    PYTHONPATH=src python3 examples/animate_sixbar_pair_v3_1.py

or:
    PYTHONPATH=src python3 examples/animate_sixbar_pair_v3_1.py \
        --output outputs/sixbar_pair_v3_1.gif
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
DEFAULT_OUTPUT = ROOT / "outputs" / "sixbar_pair_v3_1.gif"


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
}


# ---------------------------------------------------------------------------
# Mechanical model
# ---------------------------------------------------------------------------

def load_mechanisms(model_path: Path) -> tuple[dict, dict]:
    data = json.loads(model_path.read_text(encoding="utf-8"))
    try:
        return data["kinematics"]["small"], data["kinematics"]["large"]
    except KeyError as exc:
        raise ValueError(
            f"{model_path} does not contain kinematics.small / kinematics.large"
        ) from exc


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
    if dist <= 1e-12:
        raise ValueError("Primary four-bar circle centres coincide.")

    u = delta / dist
    along = (c*c - r*r + dist*dist) / (2.0*dist)
    h2 = c*c - along*along
    if h2 <= 0:
        raise ValueError("Primary four-bar does not close.")
    n = np.array((-u[1], u[0]))
    C = B + along*u + branch*math.sqrt(h2)*n

    ubc = (C - B) / c
    nbc = np.array((-ubc[1], ubc[0]))
    E = (
        B
        + float(p["primary_e_along"]) * ubc
        + float(p["primary_e_normal"]) * nbc
    )
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
    if dist <= 1e-12:
        raise ValueError("Secondary dyad circle centres coincide.")

    u = delta / dist
    along = (lef*lef - lgf*lgf + dist*dist) / (2.0*dist)
    h2 = lef*lef - along*along
    if h2 <= 0:
        raise ValueError("Secondary dyad does not close.")
    n = np.array((-u[1], u[0]))
    F = E + along*u + branch*math.sqrt(h2)*n

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

    m2 = rod*rod - transverse*transverse
    if m2 <= 0:
        raise ValueError("Finite piston rod cannot reach slider.")

    slider = longitudinal + math.sqrt(m2)
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


def mechanism_samples(
    mech: dict,
    study_angles_rad: np.ndarray,
    side: str,
) -> list[dict]:
    return [
        orient_state(sixbar_state(float(theta), mech), side)
        for theta in study_angles_rad
    ]


def place_samples(
    samples: list[dict],
    side: str,
    inner_head_x: float,
    cylinder_width: float,
) -> tuple[list[dict], dict]:
    piston_x = np.asarray([s["P"][0] for s in samples])
    pmin = float(np.min(piston_x))
    pmax = float(np.max(piston_x))

    # Visual clearance at the closed head so the piston line does not appear
    # to cut through the cylinder head. Approximate screen-thickness-based gap.
    head_clear = 0.035 * cylinder_width

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


# ---------------------------------------------------------------------------
# Thermodynamic cycle
# ---------------------------------------------------------------------------

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
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no CSV header.")
        missing = [c for c in REQUIRED_COLUMNS if c not in reader.fieldnames]
        if missing:
            raise ValueError(f"{path} is missing columns: {missing}")
        rows = list(reader)

    out = {
        c: np.asarray([float(row[c]) for row in rows], dtype=float)
        for c in REQUIRED_COLUMNS
    }

    if len(out["motor_angle_deg"]) >= 2:
        if abs(out["motor_angle_deg"][-1] - out["motor_angle_deg"][0] - 360.0) < 1e-9:
            for k in out:
                out[k] = out[k][:-1]

    return out


def periodic_interp(
    motor_deg: np.ndarray,
    values: np.ndarray,
    query_motor_deg: np.ndarray,
) -> np.ndarray:
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


# ---------------------------------------------------------------------------
# Temperature colors
# ---------------------------------------------------------------------------

def temperature_rgb(T: float, Tmin: float, Tmax: float) -> tuple[float, float, float]:
    if Tmax <= Tmin:
        a = 0.5
    else:
        a = (float(T) - Tmin) / (Tmax - Tmin)
    a = min(1.0, max(0.0, a))
    return (a, 0.0, 1.0 - a)


def colored_polyline(
    ax,
    pts: np.ndarray,
    temps: np.ndarray,
    Tmin: float,
    Tmax: float,
    linewidth=5.5,
    zorder=3,
):
    segs = np.stack([pts[:-1], pts[1:]], axis=1)
    mids = 0.5 * (temps[:-1] + temps[1:])
    lc = LineCollection(
        segs,
        colors=[temperature_rgb(t, Tmin, Tmax) for t in mids],
        linewidths=linewidth,
        capstyle="round",
        joinstyle="round",
        zorder=zorder,
    )
    ax.add_collection(lc)


def spring_points(x0, x1, y, amplitude, turns=5, n=180):
    s = np.linspace(0.0, 1.0, n)
    return np.column_stack((
        x0 + (x1-x0)*s,
        y + amplitude*np.sin(2*math.pi*turns*s),
    ))


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------

def draw_ground_symbol(ax, p, scale, lw=2.0):
    x, y = p
    ax.plot(
        [x, x], [y - 0.08*scale, y - 0.22*scale],
        color=COLORS["ground"], lw=lw, solid_capstyle="round",
    )
    base_y = y - 0.24*scale
    half = 0.18*scale
    ax.plot(
        [x-half, x+half], [base_y, base_y],
        color=COLORS["ground"], lw=lw, solid_capstyle="round",
    )
    for k in range(5):
        xi = x - half + k * (2*half/4)
        ax.plot(
            [xi - 0.04*scale, xi + 0.04*scale],
            [base_y - 0.07*scale, base_y],
            color=COLORS["ground"], lw=lw*0.7, solid_capstyle="round",
        )


def draw_joint(ax, p, r=0.075, lw=1.1):
    ax.add_patch(
        Circle(
            (p[0], p[1]), r,
            facecolor=COLORS["joint_fill"],
            edgecolor="black",
            lw=lw,
            zorder=8,
        )
    )


def draw_mechanism(ax, s, width, lw=2.5):
    A, B, C, D, E = s["A"], s["B"], s["C"], s["D"], s["E"]
    F, G, H, P = s["F"], s["G"], s["H"], s["P"]

    # Thin crank-pin trajectory.
    crank_radius = float(np.linalg.norm(B - A))
    ax.add_patch(
        Circle(
            (A[0], A[1]),
            crank_radius,
            fill=False,
            edgecolor="black",
            lw=0.85,
            zorder=0,
        )
    )

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


def draw_cylinder(
    ax,
    cyl,
    piston_x,
    gas_T,
    Tmin,
    Tmax,
    wall_lw=7.0,
    piston_lw=10.0,
):
    w = float(cyl["width"])
    head = float(cyl["x_inner"])
    outer = float(cyl["x_outer"])

    x0, x1 = sorted((head, float(piston_x)))
    if x1 > x0:
        ax.add_patch(
            Rectangle(
                (x0, -w/2),
                x1-x0,
                w,
                facecolor=temperature_rgb(gas_T, Tmin, Tmax),
                edgecolor="none",
                alpha=0.92,
                zorder=1,
            )
        )

    # Draw piston first so it stays visually underneath the cylinder outline.
    ax.plot(
        [piston_x, piston_x], [-w/2, w/2],
        color=COLORS["piston"], lw=piston_lw,
        solid_capstyle="butt", zorder=4,
    )

    ax.plot(
        [head, head], [-w/2, w/2],
        color=COLORS["cylinder"], lw=wall_lw,
        solid_capstyle="round", zorder=6,
    )
    ax.plot(
        [min(head, outer), max(head, outer)], [w/2, w/2],
        color=COLORS["cylinder"], lw=wall_lw,
        solid_capstyle="round", zorder=6,
    )
    ax.plot(
        [min(head, outer), max(head, outer)], [-w/2, -w/2],
        color=COLORS["cylinder"], lw=wall_lw,
        solid_capstyle="round", zorder=6,
    )


def draw_check_valve(
    ax,
    center,
    size,
    direction,
    is_open,
    color="black",
    lw=2.2,
):
    x, y = map(float, center)
    h = 0.14 * size
    length = 0.28 * size
    seat_h = 0.34 * size
    gap = (0.085 if is_open else 0.025) * size

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
        facecolor="white",
        edgecolor=color,
        lw=lw,
        joinstyle="round",
        zorder=9,
    )
    ax.add_patch(tri)
    ax.plot(
        [seat_x, seat_x], [y-seat_h/2, y+seat_h/2],
        color=color, lw=lw, solid_capstyle="round", zorder=9,
    )


def draw_transfer_block(
    ax,
    left_head,
    right_head,
    large_width,
    Ts,
    Tl,
    Thi,
    Tho,
    hi_valve_open,
    ho_valve_open,
    Tmin,
    Tmax,
):
    """Two direct exchange lines, no T manifolds.

    Top line represents Ho path: L -> Ho -> valve -> S
    Bottom line represents Hi path: S -> Hi -> valve -> L
    """
    xL = float(left_head)
    xR = float(right_head)
    span = xR - xL

    y_top = 0.88 * large_width
    y_bot = -0.88 * large_width
    path_lw = 6.0
    coil_amp = 0.16 * large_width

    # Top branch geometry.
    top_up_l = np.array([xL, 0.0])
    top_plateau_l = np.array([xL + 0.10*span, y_top])
    top_coil0 = np.array([xL + 0.36*span, y_top])
    top_coil1 = np.array([xL + 0.63*span, y_top])
    top_valve = np.array([xL + 0.80*span, y_top])
    top_plateau_r = np.array([xR - 0.10*span, y_top])
    top_down_r = np.array([xR, 0.0])

    top_spring = spring_points(top_coil0[0], top_coil1[0], y_top, coil_amp, turns=5, n=180)

    colored_polyline(
        ax,
        np.vstack([top_down_r, top_plateau_r, [top_coil1[0], y_top]]),
        np.array([Tl, Tho, Tho]),
        Tmin, Tmax, linewidth=path_lw, zorder=3,
    )
    colored_polyline(
        ax,
        top_spring[::-1],
        np.linspace(Tho, Tho, len(top_spring)),
        Tmin, Tmax, linewidth=path_lw, zorder=3,
    )
    colored_polyline(
        ax,
        np.vstack([[top_coil0[0], y_top], top_plateau_l, top_up_l]),
        np.array([Tho, Ts, Ts]),
        Tmin, Tmax, linewidth=path_lw, zorder=3,
    )
    draw_check_valve(
        ax,
        top_valve,
        0.56*large_width,
        direction="left",
        is_open=bool(round(ho_valve_open)),
        lw=2.2,
    )

    # Bottom branch geometry.
    bot_down_l = np.array([xL, 0.0])
    bot_plateau_l = np.array([xL + 0.10*span, y_bot])
    bot_coil0 = np.array([xL + 0.36*span, y_bot])
    bot_coil1 = np.array([xL + 0.63*span, y_bot])
    bot_valve = np.array([xL + 0.80*span, y_bot])
    bot_plateau_r = np.array([xR - 0.10*span, y_bot])
    bot_up_r = np.array([xR, 0.0])

    bot_spring = spring_points(bot_coil0[0], bot_coil1[0], y_bot, coil_amp, turns=5, n=180)

    colored_polyline(
        ax,
        np.vstack([bot_down_l, bot_plateau_l, [bot_coil0[0], y_bot]]),
        np.array([Ts, Thi, Thi]),
        Tmin, Tmax, linewidth=path_lw, zorder=3,
    )
    colored_polyline(
        ax,
        bot_spring,
        np.linspace(Thi, Thi, len(bot_spring)),
        Tmin, Tmax, linewidth=path_lw, zorder=3,
    )
    colored_polyline(
        ax,
        np.vstack([[bot_coil1[0], y_bot], bot_plateau_r, bot_up_r]),
        np.array([Thi, Tl, Tl]),
        Tmin, Tmax, linewidth=path_lw, zorder=3,
    )
    draw_check_valve(
        ax,
        bot_valve,
        0.56*large_width,
        direction="right",
        is_open=bool(round(hi_valve_open)),
        lw=2.2,
    )


def cloud(samples):
    pts = []
    for s in samples:
        for key in ("A", "B", "C", "D", "E", "F", "G", "H", "P"):
            pts.append(s[key])
    return np.asarray(pts)


# ---------------------------------------------------------------------------
# Animation
# ---------------------------------------------------------------------------

def make_animation(
    model_json: Path,
    cycle_csv: Path,
    metadata_json: Path,
    output: Path,
    frames: int,
    fps: int,
    dpi: int,
    large_width: float,
    small_width_ratio: float,
    exchanger_gap: float,
):
    small_mech, large_mech = load_mechanisms(model_json)
    cycle = load_cycle(cycle_csv)
    metadata = load_metadata(metadata_json)
    _ = metadata  # kept for compatibility / future display

    motor_query = np.linspace(0.0, 360.0, frames, endpoint=False)

    study_deg = periodic_interp(
        cycle["motor_angle_deg"],
        cycle["study_angle_deg"],
        motor_query,
    )
    study_rad = np.deg2rad(study_deg)

    Ts = periodic_interp(cycle["motor_angle_deg"], cycle["S_temperature_K"], motor_query)
    Tl = periodic_interp(cycle["motor_angle_deg"], cycle["L_temperature_K"], motor_query)
    Thi = periodic_interp(cycle["motor_angle_deg"], cycle["Hi_temperature_K"], motor_query)
    Tho = periodic_interp(cycle["motor_angle_deg"], cycle["Ho_temperature_K"], motor_query)
    hi_open = periodic_interp(cycle["motor_angle_deg"], cycle["Hi_to_L_valve_open"], motor_query)
    ho_open = periodic_interp(cycle["motor_angle_deg"], cycle["Ho_to_S_valve_open"], motor_query)

    all_internal_T = np.concatenate([
        cycle["S_temperature_K"],
        cycle["L_temperature_K"],
        cycle["Hi_temperature_K"],
        cycle["Ho_temperature_K"],
    ])
    Tmin = float(np.min(all_internal_T))
    Tmax = float(np.max(all_internal_T))

    small_width = small_width_ratio * large_width

    small_local = mechanism_samples(small_mech, study_rad, side="left")
    large_local = mechanism_samples(large_mech, study_rad, side="right")

    left_head = -exchanger_gap/2
    right_head = exchanger_gap/2

    small_samples, small_cyl = place_samples(
        small_local, "left", left_head, small_width
    )
    large_samples, large_cyl = place_samples(
        large_local, "right", right_head, large_width
    )

    pts = np.vstack([cloud(small_samples), cloud(large_samples)])

    xmin = min(
        float(np.min(pts[:,0])),
        small_cyl["x_outer"],
        small_cyl["x_inner"],
        large_cyl["x_outer"],
        large_cyl["x_inner"],
    )
    xmax = max(
        float(np.max(pts[:,0])),
        small_cyl["x_outer"],
        small_cyl["x_inner"],
        large_cyl["x_outer"],
        large_cyl["x_inner"],
    )

    ymin = min(float(np.min(pts[:,1])), -1.10*large_width)
    ymax = max(float(np.max(pts[:,1])), +1.10*large_width)

    xspan = xmax-xmin
    yspan = ymax-ymin
    mx = 0.003*xspan + 0.02
    my = 0.020*yspan + 0.04

    fig = plt.figure(figsize=(12.0, 6.0), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_aspect("equal", adjustable="box")

    def update(i):
        ax.clear()
        ax.set_xlim(xmin-mx, xmax+mx)
        ax.set_ylim(ymin-my, ymax+my)
        ax.set_aspect("equal", adjustable="box")
        ax.axis("off")

        sS = small_samples[i]
        sL = large_samples[i]

        draw_transfer_block(
            ax,
            small_cyl["x_inner"],
            large_cyl["x_inner"],
            large_width,
            Ts=float(Ts[i]),
            Tl=float(Tl[i]),
            Thi=float(Thi[i]),
            Tho=float(Tho[i]),
            hi_valve_open=float(hi_open[i]),
            ho_valve_open=float(ho_open[i]),
            Tmin=Tmin,
            Tmax=Tmax,
        )

        draw_cylinder(
            ax, small_cyl, sS["P"][0],
            gas_T=float(Ts[i]), Tmin=Tmin, Tmax=Tmax,
            wall_lw=6.5, piston_lw=9.5,
        )
        draw_cylinder(
            ax, large_cyl, sL["P"][0],
            gas_T=float(Tl[i]), Tmin=Tmin, Tmax=Tmax,
            wall_lw=7.0, piston_lw=10.0,
        )

        draw_mechanism(ax, sS, small_width, lw=2.55)
        draw_mechanism(ax, sL, large_width, lw=2.65)
        return ()

    anim = FuncAnimation(
        fig,
        update,
        frames=frames,
        interval=1000/fps,
        blit=False,
        repeat=True,
    )

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
    ap.add_argument("--exchanger-gap", type=float, default=12.0)

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
