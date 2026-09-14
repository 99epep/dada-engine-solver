#!/usr/bin/env python3
"""DADA six-bar pair animation, V4.4.

Top panel:
- six-bar pair with thermally coloured gas circuit.

Bottom panels:
- volumes
- superposed P-V diagrams
- combined P/T vs angle
- heat exchanged between working gas and exchangers

Recent visual tweaks:
- large-cylinder ports extended by +1 px; small-cylinder ports unchanged;
- valve vertical bar uses square ends;
- 1200 px width, with 20 px margins left/right/top.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Circle, Polygon, Rectangle
import numpy as np


ROOT = Path.cwd()

DEFAULT_MODEL = ROOT / "outputs" / "motor_champion_sixbar_k2.json"
DEFAULT_CYCLE = ROOT / "outputs" / "champion-6bar_thermodynamic_cycle.csv"
DEFAULT_METADATA = ROOT / "outputs" / "champion-6bar_thermodynamic_cycle_metadata.json"
DEFAULT_OUTPUT = ROOT / "outputs" / "sixbar_pair_v4_5.gif"

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
    "panel_bg": "#f0f0f0",
    "S": "#1f77b4",
    "L": "#ff7f0e",
    "Hi": "#d62728",
    "Ho": "#2ca02c",
}

ENTITY_ALIASES = {
    "small": ("small", "s"),
    "large": ("large", "l", "big"),
    "hi": ("hi",),
    "ho": ("ho",),
}


@dataclass
class SeriesInfo:
    key: str
    label: str
    values: np.ndarray
    scale: float
    unit: str

    @property
    def display_values(self) -> np.ndarray:
        return self.values * self.scale


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
    along = (c * c - r * r + dist * dist) / (2.0 * dist)
    h2 = c * c - along * along
    n = np.array((-u[1], u[0]))
    C = B + along * u + branch * math.sqrt(max(1e-12, h2)) * n

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
    along = (lef * lef - lgf * lgf + dist * dist) / (2.0 * dist)
    h2 = lef * lef - along * along
    n = np.array((-u[1], u[0]))
    F = E + along * u + branch * math.sqrt(max(1e-12, h2)) * n

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

    slider = longitudinal + math.sqrt(max(1e-12, rod * rod - transverse * transverse))
    P = origin + slider * axis

    s.update(
        {
            "F": F,
            "G": G,
            "H": H,
            "P": P,
            "axis_origin": origin,
            "axis_angle": axis_angle,
        }
    )
    return s


POINT_KEYS = ("A", "B", "C", "D", "E", "F", "G", "H", "P", "axis_origin")


def rotate_point(p: np.ndarray, angle: float) -> np.ndarray:
    c = math.cos(angle)
    s = math.sin(angle)
    return np.array((c * p[0] - s * p[1], s * p[0] + c * p[1]))


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


# ---------- cycle data discovery ----------

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
        headers = reader.fieldnames or []

    data: dict[str, np.ndarray] = {}
    for h in headers:
        vals = []
        ok = True
        for row in rows:
            try:
                vals.append(float(row[h]))
            except Exception:
                ok = False
                break
        if ok:
            data[h] = np.asarray(vals, dtype=float)

    if "motor_angle_deg" in data and len(data["motor_angle_deg"]) >= 2:
        if abs(data["motor_angle_deg"][-1] - data["motor_angle_deg"][0] - 360.0) < 1e-9:
            for k in list(data.keys()):
                data[k] = data[k][:-1]

    for req in REQUIRED_COLUMNS:
        if req not in data:
            raise KeyError(f"Missing required column in cycle CSV: {req}")
    return data


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
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def tokenize(name: str) -> list[str]:
    return [t for t in re.split(r"[^a-z0-9]+", name.lower()) if t]


def has_any_token(tokens: list[str], choices: tuple[str, ...] | list[str] | set[str]) -> bool:
    return any(c in tokens for c in choices)


def discover_entity_measure_key(data: dict[str, np.ndarray], entity: str, measure: str) -> str | None:
    aliases = ENTITY_ALIASES[entity]
    best_key = None
    best_score = -10**9
    for key in data:
        low = key.lower()
        toks = tokenize(key)
        if measure not in low:
            continue
        if not has_any_token(toks, aliases):
            continue
        score = 0
        if low.startswith(f"{entity}_") or low.startswith(f"{aliases[0]}_"):
            score += 10
        if measure in toks:
            score += 6
        if "cylinder" in toks:
            score += 2
        if low.endswith("_" + measure) or low.endswith("_" + measure + "_m3") or low.endswith("_" + measure + "_pa"):
            score += 5
        score -= len(key) * 0.01
        if score > best_score:
            best_score = score
            best_key = key
    return best_key


def discover_heat_keys(data: dict[str, np.ndarray]) -> list[str]:
    candidates: list[tuple[float, str]] = []
    for key in data:
        low = key.lower()
        toks = tokenize(key)
        if "heat" not in low:
            continue
        if "temperature" in low:
            continue
        if not ("rate" in low or "power" in low):
            continue
        score = 0.0
        if "external" in toks or "reservoir" in toks:
            score += 10
        if "heat_in" in low or ("heat" in toks and "in" in toks):
            score += 4
        if "heat_out" in low or ("heat" in toks and "out" in toks):
            score += 4
        if "source" in toks or "sink" in toks:
            score += 2
        if "total" in toks:
            score -= 1
        score -= 0.01 * len(key)
        candidates.append((score, key))
    candidates.sort(reverse=True)
    out = [k for _, k in candidates[:4]]
    # stable de-dup
    seen = set()
    final = []
    for k in out:
        if k not in seen:
            seen.add(k)
            final.append(k)
    return final


def discover_mass_keys(data: dict[str, np.ndarray]) -> list[str]:
    found = []
    order = ("small", "hi", "ho", "large")
    for entity in order:
        key = discover_entity_measure_key(data, entity, "mass")
        if key is not None:
            found.append(key)
    # fallback: any other mass columns
    for key in data:
        if key not in found and "mass" in key.lower():
            found.append(key)
    return found[:6]


def make_short_label(key: str) -> str:
    toks = tokenize(key)
    if toks[:2] == ["small", "cylinder"]:
        base = "S"
    elif toks[:2] == ["large", "cylinder"]:
        base = "L"
    elif "small" in toks or toks[:1] == ["s"]:
        base = "S"
    elif "large" in toks or toks[:1] == ["l"]:
        base = "L"
    elif "hi" in toks:
        base = "Hi"
    elif "ho" in toks:
        base = "Ho"
    elif "heat" in toks and "in" in toks:
        base = "Q̇in"
    elif "heat" in toks and "out" in toks:
        base = "Q̇out"
    else:
        base = key

    if "mass" in toks and base in {"S", "L", "Hi", "Ho"}:
        return base
    if "volume" in toks and base in {"S", "L"}:
        return base
    if "pressure" in toks and base in {"S", "L"}:
        return base
    return base


def choose_scale_and_unit(values: np.ndarray, quantity: str) -> tuple[float, str]:
    vmax = float(np.max(np.abs(values))) if len(values) else 1.0
    if quantity == "volume":
        if vmax < 1e-4:
            return 1e9, "mm³"
        if vmax < 1e-2:
            return 1e6, "cm³"
        if vmax < 1.0:
            return 1e3, "L"
        return 1.0, "m³"
    if quantity == "pressure":
        if vmax > 2e4:
            return 1e-5, "bar"
        if vmax > 50:
            return 1e-3, "kPa"
        return 1.0, "Pa"
    if quantity == "mass":
        if vmax < 1.0:
            return 1e3, "g"
        return 1.0, "kg"
    if quantity == "heat":
        if vmax >= 1000:
            return 1e-3, "kW"
        return 1.0, "W"
    return 1.0, ""


def make_series_info(key: str, label: str, values: np.ndarray, quantity: str) -> SeriesInfo:
    scale, unit = choose_scale_and_unit(values, quantity)
    return SeriesInfo(key=key, label=label, values=np.asarray(values, dtype=float), scale=scale, unit=unit)


# ---------- drawing helpers ----------

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


def data_dx_for_pixels(ax, pixels: float) -> float:
    x0_disp, y0_disp = ax.transData.transform((0.0, 0.0))
    x1_data = ax.transData.inverted().transform((x0_disp + pixels, y0_disp))[0]
    return float(x1_data)


def draw_ground_symbol(ax, p, scale, lw=2.0):
    x, y = p
    ax.plot([x, x], [y - 0.08 * scale, y - 0.22 * scale], color=COLORS["ground"], lw=lw, solid_capstyle="round")
    base_y = y - 0.24 * scale
    half = 0.18 * scale
    ax.plot([x - half, x + half], [base_y, base_y], color=COLORS["ground"], lw=lw, solid_capstyle="round")
    for k in range(5):
        xi = x - half + k * (2 * half / 4)
        ax.plot([xi - 0.04 * scale, xi + 0.04 * scale], [base_y - 0.07 * scale, base_y], color=COLORS["ground"], lw=lw * 0.7, solid_capstyle="round")


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
    ax.plot(*zip(B, E), color=COLORS["coupler"], lw=lw * 0.92, solid_capstyle="round")
    ax.plot(*zip(C, E), color=COLORS["coupler"], lw=lw * 0.92, solid_capstyle="round")

    ax.plot(*zip(G, F), color=COLORS["secondary_dyad"], lw=lw, solid_capstyle="round")
    ax.plot(*zip(E, F), color=COLORS["secondary_plate"], lw=lw, solid_capstyle="round")
    ax.plot(*zip(E, H), color=COLORS["secondary_plate"], lw=lw, solid_capstyle="round")
    ax.plot(*zip(F, H), color=COLORS["secondary_plate"], lw=lw, solid_capstyle="round")
    ax.plot(*zip(H, P), color=COLORS["rod"], lw=lw, solid_capstyle="round")

    gs = 0.40 * width
    draw_ground_symbol(ax, A, gs, lw=lw * 0.72)
    draw_ground_symbol(ax, D, gs, lw=lw * 0.72)
    draw_ground_symbol(ax, G, gs, lw=lw * 0.72)

    jr = 0.075 if width < 7 else 0.085
    for pt in (A, B, C, D, E, F, G, H):
        draw_joint(ax, pt, r=jr, lw=lw * 0.42)


def draw_cylinder(ax, cyl, piston_x, gas_T, Tmin, Tmax, wall_color, wall_lw=7.0, piston_lw=10.0):
    w = float(cyl["width"])
    head = float(cyl["x_inner"])
    outer = float(cyl["x_outer"])

    x0, x1 = sorted((head, float(piston_x)))
    if x1 > x0:
        ax.add_patch(Rectangle((x0, -w / 2), x1 - x0, w, facecolor=temperature_rgb(gas_T, Tmin, Tmax), edgecolor="none", alpha=0.92, zorder=1))

    ax.plot([piston_x, piston_x], [-w / 2, w / 2], color=COLORS["piston"], lw=piston_lw, solid_capstyle="butt", zorder=4)
    ax.plot([head, head], [-w / 2, w / 2], color=wall_color, lw=wall_lw, solid_capstyle="round", zorder=6)
    ax.plot([min(head, outer), max(head, outer)], [w / 2, w / 2], color=wall_color, lw=wall_lw, solid_capstyle="round", zorder=6)
    ax.plot([min(head, outer), max(head, outer)], [-w / 2, -w / 2], color=wall_color, lw=wall_lw, solid_capstyle="round", zorder=6)


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
        seat_x = x + length / 2
        tip_x = seat_x - gap
        base_x = x - length / 2
    else:
        seat_x = x - length / 2
        tip_x = seat_x + gap
        base_x = x + length / 2

    tri = Polygon(
        [(base_x, y - h), (base_x, y + h), (tip_x, y)],
        closed=True,
        facecolor=fill,
        edgecolor="black",
        lw=tri_lw,
        joinstyle="round",
        zorder=9,
    )
    ax.add_patch(tri)
    ax.plot([seat_x, seat_x], [y - seat_h / 2, y + seat_h / 2], color="black", lw=bar_lw, solid_capstyle="butt", zorder=9)


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

    path_lw = 3.0
    coil_amp = 0.18 * large_width
    coil_len = 0.36 * total_inside

    # Small cylinder unchanged from V3.9. Large cylinder: +1 px longer.
    pen_small = max(0.02, data_dx_for_points(ax, 1.5 * small_wall_lw) - data_dx_for_pixels(ax, 10.0))
    pen_large = max(0.02, data_dx_for_points(ax, 1.5 * large_wall_lw) - data_dx_for_pixels(ax, 9.0))

    left_inside = xL - pen_small
    right_inside = xR + pen_large

    top_valve_x = xL + 0.22 * total_inside
    bot_valve_x = xL + 0.78 * total_inside

    top_center = 0.5 * (top_valve_x + right_inside)
    top_coil0_x = top_center - 0.5 * coil_len
    top_coil1_x = top_center + 0.5 * coil_len
    top_spring = spring_points(top_coil0_x, top_coil1_x, y_top, coil_amp, turns=4.5)

    bot_center = 0.5 * (left_inside + bot_valve_x)
    bot_coil0_x = bot_center - 0.5 * coil_len
    bot_coil1_x = bot_center + 0.5 * coil_len
    bot_spring = spring_points(bot_coil0_x, bot_coil1_x, y_bot, coil_amp, turns=4.5)

    # Top / Ho: right cylinder -> exchanger -> valve -> small cylinder
    pts_top_r = np.array([[right_inside, y_top], [top_coil1_x, y_top]], dtype=float)
    colored_polyline(ax, pts_top_r, np.array([Tl, Tho]), Tmin, Tmax, linewidth=path_lw, zorder=7, capstyle="butt", joinstyle="miter")
    colored_polyline(ax, top_spring[::-1], np.full(len(top_spring), Tho), Tmin, Tmax, linewidth=path_lw, zorder=7, capstyle="round", joinstyle="round")
    pts_top_mid = np.array([[top_coil0_x, y_top], [top_valve_x, y_top]], dtype=float)
    colored_polyline(ax, pts_top_mid, np.array([Tho, Tho]), Tmin, Tmax, linewidth=path_lw, zorder=7, capstyle="butt", joinstyle="miter")
    pts_top_l = np.array([[top_valve_x, y_top], [left_inside, y_top]], dtype=float)
    colored_polyline(ax, pts_top_l, np.array([Ts, Ts]), Tmin, Tmax, linewidth=path_lw, zorder=7, capstyle="butt", joinstyle="miter")
    draw_check_valve(ax, (top_valve_x, y_top), (0.90 * large_width) * (2.0 / 3.0), direction="left", is_open=bool(round(ho_valve_open)))

    # Bottom / Hi: small cylinder -> exchanger -> valve -> right cylinder
    pts_bot_l = np.array([[left_inside, y_bot], [bot_coil0_x, y_bot]], dtype=float)
    colored_polyline(ax, pts_bot_l, np.array([Ts, Thi]), Tmin, Tmax, linewidth=path_lw, zorder=7, capstyle="butt", joinstyle="miter")
    colored_polyline(ax, bot_spring, np.full(len(bot_spring), Thi), Tmin, Tmax, linewidth=path_lw, zorder=7, capstyle="round", joinstyle="round")
    pts_bot_mid = np.array([[bot_coil1_x, y_bot], [bot_valve_x, y_bot]], dtype=float)
    colored_polyline(ax, pts_bot_mid, np.array([Thi, Thi]), Tmin, Tmax, linewidth=path_lw, zorder=7, capstyle="butt", joinstyle="miter")
    pts_bot_r = np.array([[bot_valve_x, y_bot], [right_inside, y_bot]], dtype=float)
    colored_polyline(ax, pts_bot_r, np.array([Tl, Tl]), Tmin, Tmax, linewidth=path_lw, zorder=7, capstyle="butt", joinstyle="miter")
    draw_check_valve(ax, (bot_valve_x, y_bot), (0.90 * large_width) * (2.0 / 3.0), direction="right", is_open=bool(round(hi_valve_open)))


def cloud(samples):
    pts = []
    for s in samples:
        for key in ("A", "B", "C", "D", "E", "F", "G", "H", "P"):
            pts.append(s[key])
    return np.asarray(pts)


# ---------- chart helpers ----------

def init_chart_axis(ax, title: str):
    ax.set_facecolor("white")
    ax.set_title(title, fontsize=9, pad=2)
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.tick_params(labelsize=7, length=2)


def add_static_angle_cursor(ax, angle0: float):
    return ax.axvline(angle0, color="black", lw=0.9, alpha=0.75)


def build_time_series_panel(
    ax,
    xdeg: np.ndarray,
    series_list: list[SeriesInfo],
    title: str,
    y_label: str,
    color_map: dict[str, str] | None = None,
):
    init_chart_axis(ax, title)
    lines = []
    markers = []
    for s in series_list:
        color = None if color_map is None else color_map.get(s.label)
        line, = ax.plot(xdeg, s.display_values, lw=1.4, label=s.label, color=color)
        marker, = ax.plot(
            [xdeg[0]], [s.display_values[0]],
            marker="o", ms=3.2, color=line.get_color(), linestyle="None",
            zorder=12,
        )
        lines.append(line)
        markers.append(marker)
    ax.set_xlim(0.0, 360.0)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_ylabel(y_label, fontsize=8)
    legend = None
    if len(series_list) <= 6:
        legend = ax.legend(
            loc="center right",
            fontsize=7,
            frameon=False,
            handlelength=1.8,
            borderpad=0.2,
        )
        legend.set_zorder(4)
    cursor = add_static_angle_cursor(ax, xdeg[0])
    cursor.set_zorder(20)
    return {"lines": lines, "markers": markers, "cursor": cursor, "legend": legend}


def build_pv_panel(ax, vol_small: SeriesInfo | None, vol_large: SeriesInfo | None, p_small: SeriesInfo | None, p_large: SeriesInfo | None):
    init_chart_axis(ax, "P-V")
    artists = {"markers": []}
    if vol_small is not None and p_small is not None:
        line_s, = ax.plot(
            vol_small.display_values, p_small.display_values,
            lw=1.4, label="S", color=COLORS["S"],
        )
        mark_s, = ax.plot(
            [vol_small.display_values[0]], [p_small.display_values[0]],
            marker="o", ms=3.5, color=COLORS["S"], linestyle="None", zorder=12,
        )
        artists["markers"].append((mark_s, vol_small.display_values, p_small.display_values))
    if vol_large is not None and p_large is not None:
        line_l, = ax.plot(
            vol_large.display_values, p_large.display_values,
            lw=1.4, label="L", color=COLORS["L"],
        )
        mark_l, = ax.plot(
            [vol_large.display_values[0]], [p_large.display_values[0]],
            marker="o", ms=3.5, color=COLORS["L"], linestyle="None", zorder=12,
        )
        artists["markers"].append((mark_l, vol_large.display_values, p_large.display_values))
    if ax.lines:
        leg = ax.legend(loc="center right", fontsize=7, frameon=False, handlelength=1.8, borderpad=0.2)
        leg.set_zorder(4)
    ax.set_xlabel(vol_small.unit if vol_small is not None else (vol_large.unit if vol_large is not None else ""), fontsize=8)
    ax.set_ylabel(p_small.unit if p_small is not None else (p_large.unit if p_large is not None else ""), fontsize=8)
    return artists


def update_time_series_panel(panel, idx: int, xcur: float):
    panel["cursor"].set_xdata([xcur, xcur])
    for line, marker in zip(panel["lines"], panel["markers"]):
        xd = line.get_xdata()
        yd = line.get_ydata()
        marker.set_data([xd[idx]], [yd[idx]])


def update_pv_panel(panel, idx: int):
    for marker, xs, ys in panel["markers"]:
        marker.set_data([xs[idx]], [ys[idx]])


def build_pt_panel(
    ax,
    xdeg: np.ndarray,
    p_small: SeriesInfo,
    p_large: SeriesInfo,
    t_small_c: np.ndarray,
    t_large_c: np.ndarray,
):
    init_chart_axis(ax, "P — / T ···")
    ax_t = ax.twinx()
    ax_t.tick_params(labelsize=7, length=2)

    lp_s, = ax.plot(xdeg, p_small.display_values, color=COLORS["S"], lw=1.35, label="S")
    lp_l, = ax.plot(xdeg, p_large.display_values, color=COLORS["L"], lw=1.35, label="L")
    lt_s, = ax_t.plot(xdeg, t_small_c, color=COLORS["S"], lw=1.15, ls=":", alpha=0.95)
    lt_l, = ax_t.plot(xdeg, t_large_c, color=COLORS["L"], lw=1.15, ls=":", alpha=0.95)

    mp_s, = ax.plot([xdeg[0]], [p_small.display_values[0]], "o", ms=3.0, color=COLORS["S"], zorder=12)
    mp_l, = ax.plot([xdeg[0]], [p_large.display_values[0]], "o", ms=3.0, color=COLORS["L"], zorder=12)
    mt_s, = ax_t.plot([xdeg[0]], [t_small_c[0]], "o", ms=2.7, color=COLORS["S"], zorder=12)
    mt_l, = ax_t.plot([xdeg[0]], [t_large_c[0]], "o", ms=2.7, color=COLORS["L"], zorder=12)

    ax.set_xlim(0.0, 360.0)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_ylabel(p_small.unit, fontsize=8)
    ax_t.set_ylabel("°C", fontsize=8)

    leg = ax.legend(loc="center right", fontsize=7, frameon=False, handlelength=1.8, borderpad=0.2)
    leg.set_zorder(4)
    cursor = add_static_angle_cursor(ax, xdeg[0])
    cursor.set_zorder(20)

    return {
        "ax_t": ax_t,
        "p_lines": (lp_s, lp_l),
        "t_lines": (lt_s, lt_l),
        "p_markers": (mp_s, mp_l),
        "t_markers": (mt_s, mt_l),
        "cursor": cursor,
    }


def update_pt_panel(panel, idx: int, xcur: float):
    panel["cursor"].set_xdata([xcur, xcur])
    for line, marker in zip(panel["p_lines"], panel["p_markers"]):
        marker.set_data([line.get_xdata()[idx]], [line.get_ydata()[idx]])
    for line, marker in zip(panel["t_lines"], panel["t_markers"]):
        marker.set_data([line.get_xdata()[idx]], [line.get_ydata()[idx]])


def choose_ylim(series_list: list[SeriesInfo], pad=0.06):
    vals = np.concatenate([s.display_values for s in series_list]) if series_list else np.array([0.0, 1.0])
    lo = float(np.min(vals))
    hi = float(np.max(vals))
    if hi <= lo:
        hi = lo + 1.0
    margin = pad * (hi - lo)
    return lo - margin, hi + margin


# ---------- main ----------

def make_animation(
    model_json: Path,
    cycle_csv: Path,
    metadata_json: Path,
    output: Path,
    frames: int,
    fps: int,
    dpi: int,
    width_px: int,
    height_px: int,
    large_width: float,
    small_width_ratio: float,
    exchanger_gap: float,
):
    small_mech, large_mech = load_mechanisms(model_json)
    cycle = load_cycle(cycle_csv)
    _meta = load_metadata(metadata_json)

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
    mx = 0.0015 * xspan + 0.025
    my = 0.008 * yspan + 0.02

    # Exact columns from outputs/champion-6bar_thermodynamic_cycle.csv.
    vol_small_key = "S_volume_m3"
    vol_large_key = "L_volume_m3"
    p_small_key = "S_pressure_Pa"
    p_large_key = "L_pressure_Pa"
    heat_keys = [
        "Hi_wall_to_gas_heat_W",
        "Ho_wall_to_gas_heat_W",
    ]
    mass_keys = [
        "S_mass_kg",
        "Hi_mass_kg",
        "Ho_mass_kg",
        "L_mass_kg",
    ]

    def resampled(key: str) -> np.ndarray:
        return periodic_interp(cycle["motor_angle_deg"], cycle[key], motor_query)

    vol_series = []
    if vol_small_key:
        vol_series.append(make_series_info(vol_small_key, "S", resampled(vol_small_key), "volume"))
    if vol_large_key:
        vol_series.append(make_series_info(vol_large_key, "L", resampled(vol_large_key), "volume"))

    p_small = make_series_info(p_small_key, "S", resampled(p_small_key), "pressure") if p_small_key else None
    p_large = make_series_info(p_large_key, "L", resampled(p_large_key), "pressure") if p_large_key else None
    temp_small_c = resampled("S_temperature_K") - 273.15
    temp_large_c = resampled("L_temperature_K") - 273.15
    vol_small = next((s for s in vol_series if s.label == "S"), None)
    vol_large = next((s for s in vol_series if s.label == "L"), None)

    heat_labels = {
        "Hi_wall_to_gas_heat_W": "Hi",
        "Ho_wall_to_gas_heat_W": "Ho",
    }
    heat_series = [
        make_series_info(k, heat_labels[k], resampled(k), "heat")
        for k in heat_keys
    ]
    mass_series = [
        make_series_info(k, make_short_label(k), resampled(k), "mass")
        for k in mass_keys
    ]

    fig = plt.figure(figsize=(width_px / dpi, height_px / dpi), dpi=dpi, facecolor=COLORS["panel_bg"])
    fig.subplots_adjust(0, 0, 1, 1)

    left = 28 / width_px
    right = 1 - 22 / width_px
    top = 1 - 20 / height_px
    bottom = 30 / height_px

    gs = GridSpec(
        3,
        2,
        figure=fig,
        left=left,
        right=right,
        top=top,
        bottom=bottom,
        hspace=0.11,
        wspace=0.18,
        height_ratios=[2.02, 0.42, 0.42],
    )

    ax_mech = fig.add_subplot(gs[0, :])
    ax_vol = fig.add_subplot(gs[1, 0])
    ax_pv = fig.add_subplot(gs[1, 1])
    ax_pt = fig.add_subplot(gs[2, 0])
    ax_heat = fig.add_subplot(gs[2, 1])

    for ax in (ax_mech, ax_vol, ax_pv, ax_pt, ax_heat):
        ax.set_facecolor(COLORS["panel_bg"] if ax is ax_mech else "white")

    # Static charts.
    state_colors = {
        "S": COLORS["S"],
        "L": COLORS["L"],
        "Hi": COLORS["Hi"],
        "Ho": COLORS["Ho"],
    }

    vol_panel = build_time_series_panel(
        ax_vol, motor_query, vol_series, "Volumes",
        vol_series[0].unit if vol_series else "",
        color_map=state_colors,
    )
    if vol_series:
        ax_vol.set_ylim(*choose_ylim(vol_series))

    pv_panel = build_pv_panel(ax_pv, vol_small, vol_large, p_small, p_large)

    pt_panel = build_pt_panel(
        ax_pt, motor_query, p_small, p_large, temp_small_c, temp_large_c
    )

    heat_panel = build_time_series_panel(
        ax_heat, motor_query, heat_series, "Heat gas↔HX",
        heat_series[0].unit if heat_series else "W",
        color_map=state_colors,
    )
    if heat_series:
        ax_heat.set_ylim(*choose_ylim(heat_series, pad=0.08))
    ax_heat.axhline(0.0, color="black", lw=1.2, alpha=0.85, zorder=2)
    ax_heat.set_xlabel("θ (deg)", fontsize=8)



    small_wall_lw = 6.5
    small_piston_lw = 9.5
    large_wall_lw = 7.0
    large_piston_lw = 10.0

    def update(i):
        ax_mech.clear()
        ax_mech.set_facecolor(COLORS["panel_bg"])
        ax_mech.set_xlim(xmin - mx, xmax + mx)
        ax_mech.set_ylim(ymin - my, ymax + my)
        ax_mech.set_aspect("equal", adjustable="box")
        ax_mech.axis("off")

        sS = small_samples[i]
        sL = large_samples[i]

        draw_cylinder(
            ax_mech, small_cyl, sS["P"][0],
            gas_T=float(Ts[i]), Tmin=Tmin, Tmax=Tmax,
            wall_color=COLORS["S"],
            wall_lw=small_wall_lw, piston_lw=small_piston_lw,
        )
        draw_cylinder(
            ax_mech, large_cyl, sL["P"][0],
            gas_T=float(Tl[i]), Tmin=Tmin, Tmax=Tmax,
            wall_color=COLORS["L"],
            wall_lw=large_wall_lw, piston_lw=large_piston_lw,
        )

        draw_transfer_block(
            ax_mech,
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

        draw_mechanism(ax_mech, sS, small_width, lw=2.55)
        draw_mechanism(ax_mech, sL, large_width, lw=2.65)

        update_time_series_panel(vol_panel, i, motor_query[i])
        update_pv_panel(pv_panel, i)
        update_pt_panel(pt_panel, i, motor_query[i])
        update_time_series_panel(heat_panel, i, motor_query[i])
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

    ap.add_argument("--width-px", type=int, default=1200)
    ap.add_argument("--height-px", type=int, default=768)
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
        width_px=args.width_px,
        height_px=args.height_px,
        large_width=args.large_cylinder_width,
        small_width_ratio=args.small_width_ratio,
        exchanger_gap=args.exchanger_gap,
    )


if __name__ == "__main__":
    main()
