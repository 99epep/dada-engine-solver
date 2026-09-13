#!/usr/bin/env python3
"""Diagnose Stage 2C H trajectory and crankshaft clearance."""

from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np


def load_best(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("best") is not None:
        return data["best"]
    candidates = [r.get("best_dense") for r in data.get("runs", []) if r.get("best_dense") is not None]
    if not candidates:
        raise ValueError("No best candidate found.")
    return min(candidates, key=lambda x: x["score"])


def point_segment_distance(p, a, b):
    ab = b - a
    den = float(ab @ ab)
    if den <= 1e-20:
        return float(np.linalg.norm(p - a))
    t = float((p - a) @ ab / den)
    t = min(1.0, max(0.0, t))
    q = a + t * ab
    return float(np.linalg.norm(p - q))


def point_in_triangle(p, a, b, c):
    def cross2(u, v):
        return float(u[0]*v[1] - u[1]*v[0])
    s1 = cross2(b-a, p-a)
    s2 = cross2(c-b, p-b)
    s3 = cross2(a-c, p-c)
    eps = 1e-12
    has_neg = (s1 < -eps) or (s2 < -eps) or (s3 < -eps)
    has_pos = (s1 > eps) or (s2 > eps) or (s3 > eps)
    return not (has_neg and has_pos)


def triangle_signed_clearance_from_origin(e, f, h):
    o = np.zeros(2)
    d = min(
        point_segment_distance(o, e, f),
        point_segment_distance(o, f, h),
        point_segment_distance(o, h, e),
    )
    return -d if point_in_triangle(o, e, f, h) else d


def primary_state(theta, primary):
    g = float(primary["ground"])
    c = float(primary["coupler"])
    r = float(primary["rocker"])
    branch = int(primary.get("assembly_branch", -1))
    e_along = float(primary["E_along"])
    e_normal = float(primary["E_normal"])
    phase = float(primary["phase"])

    angle = theta + phase
    B = np.column_stack((np.cos(angle), np.sin(angle)))
    D = np.tile(np.array((g, 0.0)), (len(theta), 1))

    delta = D - B
    dist = np.linalg.norm(delta, axis=1)
    direction = delta / dist[:, None]
    along = (c*c - r*r + dist*dist) / (2.0*dist)
    h2 = c*c - along*along
    if np.any(h2 <= 0):
        raise ValueError("Primary does not close.")

    h = np.sqrt(h2)
    left = np.column_stack((-direction[:, 1], direction[:, 0]))
    C = B + along[:, None]*direction + branch*h[:, None]*left

    u = (C - B) / c
    n = np.column_stack((-u[:, 1], u[:, 0]))
    E = B + e_along*u + e_normal*n
    return B, C, D, E


def mechanism_states(best, frames):
    theta = np.linspace(0.0, 2.0*math.pi, frames, endpoint=False)
    B, C, D, E = primary_state(theta, best["primary"])

    p = best["parameters"]
    branch = int(best["second_branch"])
    G = np.array((float(p["second_pivot_x"]), float(p["second_pivot_y"])))
    lef = float(p["link_EF"])
    lgf = float(p["link_GF"])
    ha = float(p["H_along_over_EF"])
    hn = float(p["H_normal_over_EF"])

    delta = G - E
    dist = np.linalg.norm(delta, axis=1)
    direction = delta / dist[:, None]
    along = (lef*lef - lgf*lgf + dist*dist) / (2.0*dist)
    h2 = lef*lef - along*along
    if np.any(h2 <= 0):
        raise ValueError("Second dyad does not close.")
    h = np.sqrt(h2)
    left = np.column_stack((-direction[:, 1], direction[:, 0]))
    F = E + along[:, None]*direction + branch*h[:, None]*left

    ef_u = (F - E) / lef
    ef_n = np.column_stack((-ef_u[:, 1], ef_u[:, 0]))
    H = E + ha*lef*ef_u + hn*lef*ef_n

    rod = float(p["piston_rod"])
    axis_angle = float(p["slider_axis_angle"])
    axis_offset = float(p["slider_axis_offset"])
    axis = np.array((math.cos(axis_angle), math.sin(axis_angle)))
    normal = np.array((-axis[1], axis[0]))
    axis_origin = axis_offset * normal

    rel = H - axis_origin
    longitudinal = rel @ axis
    transverse = rel @ normal
    margin2 = rod*rod - transverse*transverse
    if np.any(margin2 <= 0):
        raise ValueError("Piston rod does not close.")
    margin = np.sqrt(margin2)
    slider = longitudinal + margin
    P = axis_origin + slider[:, None]*axis

    return theta, E, F, G, H, P, axis, rod, slider


def acute_angle_deg(u, v):
    u = u / np.linalg.norm(u)
    v = v / np.linalg.norm(v)
    c = abs(float(u @ v))
    c = min(1.0, max(-1.0, c))
    return math.degrees(math.acos(c))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("result_json", type=Path)
    ap.add_argument("--frames", type=int, default=1440)
    ap.add_argument("--shaft-radius", type=float, default=0.0,
                    help="Crankshaft radius in crank-radius units.")
    args = ap.parse_args()

    best = load_best(args.result_json)
    theta, E, F, G, H, P, axis, rod, slider = mechanism_states(best, args.frames)
    stroke = float(np.ptp(slider))

    center = np.mean(H, axis=0)
    X = H - center
    _, _, vh = np.linalg.svd(X, full_matrices=False)
    line_dir = vh[0].copy()
    if float(line_dir @ axis) < 0:
        line_dir *= -1
    line_normal = np.array((-line_dir[1], line_dir[0]))

    along = X @ line_dir
    perp = X @ line_normal
    line_span = float(np.ptp(along))
    rms_perp = float(np.sqrt(np.mean(perp**2)))
    max_perp = float(np.max(np.abs(perp)))
    perp_span = float(np.ptp(perp))
    line_axis_angle = acute_angle_deg(line_dir, axis)

    g_to_line = abs(float((G - center) @ line_normal))
    gh = np.linalg.norm(H - G, axis=1)

    rod_vec = P - H
    rod_cos = np.clip(np.abs((rod_vec @ axis) / rod), 0.0, 1.0)
    rod_angles = np.degrees(np.arccos(rod_cos))

    clearances = np.array([
        triangle_signed_clearance_from_origin(e, f, h)
        for e, f, h in zip(E, F, H)
    ])
    inside_count = int(np.sum(clearances < 0.0))
    min_clearance = float(np.min(clearances))
    closest_boundary = float(np.min(np.abs(clearances)))

    print(f"Candidate: {args.result_json}\n")

    print("H NATURAL BEST-FIT LINE")
    print(f"  line angle                     = {math.degrees(math.atan2(line_dir[1], line_dir[0])): .3f} deg")
    print(f"  cylinder axis angle            = {math.degrees(math.atan2(axis[1], axis[0])): .3f} deg")
    print(f"  acute line/cylinder angle      = {line_axis_angle: .3f} deg")
    print(f"  H span along line / stroke     = {line_span/stroke: .6f}")
    print(f"  RMS departure / stroke         = {rms_perp/stroke: .6f}")
    print(f"  max departure / stroke         = {max_perp/stroke: .6f}")
    print(f"  full transverse span / stroke  = {perp_span/stroke: .6f}\n")

    print("RELATION TO FIXED PIVOT G")
    print(f"  distance G -> best H line / r  = {g_to_line: .6f}")
    print(f"  min |GH| / r                   = {float(np.min(gh)): .6f}")
    print(f"  max |GH| / r                   = {float(np.max(gh)): .6f}\n")

    print("PISTON ROD")
    print(f"  piston stroke / crank          = {stroke: .6f}")
    print(f"  rod / stroke                   = {rod/stroke: .6f}")
    print(f"  maximum rod angle              = {float(np.max(rod_angles)): .3f} deg")
    print(f"  RMS rod angle                  = {float(np.sqrt(np.mean(rod_angles**2))): .3f} deg\n")

    print("CRANKSHAFT / SECONDARY PLATE EFH")
    print(f"  frames with A inside EFH       = {inside_count} / {args.frames}")
    print(f"  signed min clearance A -> EFH  = {min_clearance: .6f} r")
    print(f"  closest plate boundary to A    = {closest_boundary: .6f} r")
    print(f"  min |AE|                       = {float(np.min(np.linalg.norm(E,axis=1))): .6f} r")
    print(f"  min |AF|                       = {float(np.min(np.linalg.norm(F,axis=1))): .6f} r")
    print(f"  min |AH|                       = {float(np.min(np.linalg.norm(H,axis=1))): .6f} r")

    if args.shaft_radius > 0:
        margin = min_clearance - args.shaft_radius
        print(f"  requested shaft radius         = {args.shaft_radius: .6f} r")
        print(f"  shaft radial margin            = {margin: .6f} r")
        print(f"  2-D shaft clearance OK?        = {'YES' if margin > 0 else 'NO'}")
    else:
        print(f"  zero-radius shaft axis clear?  = {'YES' if inside_count == 0 and min_clearance > 0 else 'NO'}")

    print("\nJSON METRICS")
    for key in (
        "position_rms",
        "derivative_rms",
        "stroke_over_crank",
        "H_lateral_span_over_piston_stroke",
        "minimum_primary_transmission_sine",
        "minimum_secondary_transmission_sine",
        "minimum_rod_axis_cosine",
    ):
        if key in best:
            print(f"  {key:40s} = {best[key]}")


if __name__ == "__main__":
    main()
