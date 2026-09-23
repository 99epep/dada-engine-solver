#!/usr/bin/env python3
"""Stage 2F: local refinement around one selected Stage 2E candidate.

Default target:
- restart 2
- secondary branch +1

The local search keeps the attractive architecture while giving the optimizer
enough room to improve the final piston motion.

Recommended first run:

PYTHONPATH=src python3 examples/search_small_sixbar_stage2f.py \
  outputs/small_sixbar_stage2e.json \
  --seed-restart 2 \
  --seed-branch +1 \
  --iterations 300 \
  --restarts 3 \
  --stroke-floor 2.5 \
  --h-line-rms-max 0.07 \
  --h-axis-angle-max 5 \
  --crank-clearance-floor 0.50 \
  --output outputs/small_sixbar_stage2f_r2_bp1.json
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution


ROOT = Path.cwd()
TARGET = ROOT / "outputs" / "motor_champion_motion_target.csv"
PRIMARY_BRANCH = -1

NAMES = [
    "primary_ground",
    "primary_coupler",
    "primary_rocker",
    "primary_E_along",
    "primary_E_normal",
    "primary_phase",
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


def wrap_pi(a: float) -> float:
    return (a + math.pi) % (2.0 * math.pi) - math.pi


def load_seed(path: Path, restart: int, branch: int) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))

    for run in data.get("runs", []):
        if (
            int(run.get("restart", -999)) == restart
            and int(run.get("second_branch", 0)) == branch
            and run.get("best_dense") is not None
        ):
            return run["best_dense"]

    raise ValueError(
        f"No candidate restart={restart}, branch={branch:+d} found in {path}"
    )


def seed_vector(best: dict) -> np.ndarray:
    p = best["parameters"]

    if all(name in p for name in NAMES):
        return np.array([float(p[name]) for name in NAMES], dtype=float)

    primary = best["primary"]
    return np.array([
        float(primary["ground"]),
        float(primary["coupler"]),
        float(primary["rocker"]),
        float(primary["E_along"]),
        float(primary["E_normal"]),
        wrap_pi(float(primary["phase"])),
        float(p["second_pivot_x"]),
        float(p["second_pivot_y"]),
        float(p["link_EF"]),
        float(p["link_GF"]),
        float(p["H_along_over_EF"]),
        float(p["H_normal_over_EF"]),
        float(p["piston_rod"]),
        float(p["slider_axis_offset"]),
        float(p["slider_axis_angle"]),
    ], dtype=float)


def positive_relative(v: float, frac: float, lo_floor: float = 1e-3):
    lo = max(lo_floor, v * (1.0 - frac))
    hi = max(lo + 1e-6, v * (1.0 + frac))
    return (lo, hi)


def signed_local(v: float, frac: float, min_span: float):
    span = max(abs(v) * frac, min_span)
    return (v - span, v + span)


def angle_local(v: float, span_deg: float):
    span = math.radians(span_deg)
    lo = max(-math.pi, v - span)
    hi = min(math.pi, v + span)
    if hi - lo < math.radians(2.0):
        # Rare boundary case: allow a useful interval on the available side.
        if v > 0:
            lo = max(-math.pi, hi - 2.0 * span)
        else:
            hi = min(math.pi, lo + 2.0 * span)
    return (lo, hi)


def unwrapped_angle_local(v: float, span_deg: float):
    span = math.radians(span_deg)
    return (v - span, v + span)


def local_bounds(x: np.ndarray, args) -> list[tuple[float, float]]:
    b = []

    # Primary lengths: modest local freedom.
    b.append(positive_relative(x[0], args.primary_length_fraction, 0.5))
    b.append(positive_relative(x[1], args.primary_length_fraction, 0.5))
    b.append(positive_relative(x[2], args.primary_length_fraction, 0.5))

    # Primary output point: deliberately somewhat freer.
    b.append(signed_local(x[3], args.primary_point_fraction, 0.50))
    b.append(signed_local(x[4], args.primary_point_fraction, 0.50))

    b.append(unwrapped_angle_local(x[5], args.primary_phase_span_deg))

    # Fixed secondary pivot: local XY box.
    b.append((x[6] - args.pivot_span, x[6] + args.pivot_span))
    b.append((x[7] - args.pivot_span, x[7] + args.pivot_span))

    # Compact second dyad lengths.
    ef = positive_relative(x[8], args.secondary_length_fraction, 0.25)
    gf = positive_relative(x[9], args.secondary_length_fraction, 0.25)
    b.append((ef[0], min(args.max_EF, ef[1])))
    b.append((gf[0], min(args.max_GF, gf[1])))

    # H position on EF plate.
    b.append(signed_local(x[10], args.h_point_fraction, 0.25))
    b.append(signed_local(x[11], args.h_point_fraction, 0.25))

    # Rod may vary more; it is not the main compactness metric.
    b.append(positive_relative(x[12], args.rod_fraction, 1.0))

    b.append((x[13] - args.axis_offset_span, x[13] + args.axis_offset_span))
    b.append(angle_local(x[14], args.axis_angle_span_deg))

    return b


def primary(theta: np.ndarray, v: np.ndarray):
    g, c, r, ea, en, phase = map(float, v[:6])

    if abs(g - 1.0) <= abs(c - r) + 1e-10:
        raise ValueError("primary inner closure")
    if c + r <= g + 1.0 + 1e-10:
        raise ValueError("primary outer closure")

    a = theta + phase
    B = np.column_stack((np.cos(a), np.sin(a)))
    Bd = np.column_stack((-np.sin(a), np.cos(a)))
    D = np.array((g, 0.0))

    delta = D - B
    dist = np.linalg.norm(delta, axis=1)
    if np.any(dist <= 1e-12):
        raise ValueError("primary coincident centers")

    u0 = delta / dist[:, None]
    along = (c*c - r*r + dist*dist) / (2.0 * dist)
    h2 = c*c - along*along
    if np.any(h2 <= 1e-10):
        raise ValueError("primary toggle")

    h = np.sqrt(h2)
    n0 = np.column_stack((-u0[:, 1], u0[:, 0]))
    C = B + along[:, None]*u0 + PRIMARY_BRANCH*h[:, None]*n0

    bc = C - B
    dc = C - D
    det = bc[:, 0]*dc[:, 1] - bc[:, 1]*dc[:, 0]
    if np.any(np.abs(det) <= 1e-10):
        raise ValueError("primary singular")

    rhs = np.sum(bc * Bd, axis=1)
    Cd = np.column_stack((
        rhs * dc[:, 1] / det,
        -rhs * dc[:, 0] / det,
    ))

    u = bc / c
    ud = (Cd - Bd) / c
    n = np.column_stack((-u[:, 1], u[:, 0]))
    nd = np.column_stack((-ud[:, 1], ud[:, 0]))

    E = B + ea*u + en*n
    Ed = Bd + ea*ud + en*nd

    sine = float(np.min(np.abs(det) / (c*r)))
    pdata = {
        "ground": g,
        "coupler": c,
        "rocker": r,
        "assembly_branch": PRIMARY_BRANCH,
        "E_along": ea,
        "E_normal": en,
        "phase": phase,
    }
    return E, Ed, sine, pdata


def best_line_metrics(H: np.ndarray, stroke: float, cylinder_axis: np.ndarray):
    center = np.mean(H, axis=0)
    X = H - center
    cov = X.T @ X
    values, vectors = np.linalg.eigh(cov)
    line_dir = vectors[:, int(np.argmax(values))]
    line_dir /= np.linalg.norm(line_dir)

    if float(line_dir @ cylinder_axis) < 0:
        line_dir *= -1.0

    normal = np.array((-line_dir[1], line_dir[0]))
    along = X @ line_dir
    perp = X @ normal

    rms = float(np.sqrt(np.mean(perp*perp))) / stroke
    max_abs = float(np.max(np.abs(perp))) / stroke
    full_span = float(np.ptp(perp)) / stroke
    along_span = float(np.ptp(along)) / stroke

    c = min(1.0, max(0.0, abs(float(line_dir @ cylinder_axis))))
    angle_deg = math.degrees(math.acos(c))

    return {
        "H_best_line_rms_over_stroke": rms,
        "H_best_line_max_over_stroke": max_abs,
        "H_best_line_full_span_over_stroke": full_span,
        "H_best_line_along_span_over_stroke": along_span,
        "H_best_line_axis_angle_deg": angle_deg,
        "H_best_line_angle_deg": math.degrees(math.atan2(line_dir[1], line_dir[0])),
    }


def segment_distance_origin_batch(A: np.ndarray, B: np.ndarray):
    AB = B - A
    den = np.sum(AB*AB, axis=1)
    t = -np.sum(A*AB, axis=1) / np.maximum(den, 1e-30)
    t = np.clip(t, 0.0, 1.0)
    Q = A + t[:, None]*AB
    return np.linalg.norm(Q, axis=1)


def signed_triangle_clearance_origin(E: np.ndarray, F: np.ndarray, H: np.ndarray):
    d = np.minimum.reduce([
        segment_distance_origin_batch(E, F),
        segment_distance_origin_batch(F, H),
        segment_distance_origin_batch(H, E),
    ])

    def cross_to_origin(A, B):
        AB = B - A
        AO = -A
        return AB[:, 0]*AO[:, 1] - AB[:, 1]*AO[:, 0]

    s1 = cross_to_origin(E, F)
    s2 = cross_to_origin(F, H)
    s3 = cross_to_origin(H, E)
    eps = 1e-12

    inside = (
        ((s1 >= -eps) & (s2 >= -eps) & (s3 >= -eps))
        | ((s1 <= eps) & (s2 <= eps) & (s3 <= eps))
    )
    signed = np.where(inside, -d, d)
    return float(np.min(signed)), int(np.sum(inside))


def evaluate(
    v: np.ndarray,
    second_branch: int,
    theta: np.ndarray,
    target_q: np.ndarray,
    target_dq: np.ndarray,
    maximum_eh: float,
):
    E, Ed, primary_sine, pdata = primary(theta, v)

    gx, gy, lef, lgf, ha, hn, rod, axis_offset, axis_angle = map(float, v[6:])

    eh_ratio = math.hypot(ha, hn)
    eh = lef * eh_ratio
    if eh > maximum_eh:
        raise ValueError("EH too large")

    G = np.array((gx, gy))
    delta = G - E
    dist = np.linalg.norm(delta, axis=1)
    if np.any(dist <= 1e-10):
        raise ValueError("secondary coincident centers")

    u0 = delta / dist[:, None]
    along = (lef*lef - lgf*lgf + dist*dist) / (2.0 * dist)
    h2 = lef*lef - along*along
    if np.any(h2 <= 1e-10):
        raise ValueError("secondary closure")

    h = np.sqrt(h2)
    n0 = np.column_stack((-u0[:, 1], u0[:, 0]))
    F = E + along[:, None]*u0 + second_branch*h[:, None]*n0

    ef = F - E
    gf = F - G
    det = ef[:, 0]*gf[:, 1] - ef[:, 1]*gf[:, 0]
    if np.any(np.abs(det) <= 1e-10):
        raise ValueError("secondary singular")

    rhs = np.sum(ef * Ed, axis=1)
    Fd = np.column_stack((
        rhs * gf[:, 1] / det,
        -rhs * gf[:, 0] / det,
    ))
    secondary_sine = float(np.min(np.abs(det) / (lef*lgf)))

    u = ef / lef
    ud = (Fd - Ed) / lef
    n = np.column_stack((-u[:, 1], u[:, 0]))
    nd = np.column_stack((-ud[:, 1], ud[:, 0]))

    H = E + ha*lef*u + hn*lef*n
    Hd = Ed + ha*lef*ud + hn*lef*nd

    axis = np.array((math.cos(axis_angle), math.sin(axis_angle)))
    axis_normal = np.array((-axis[1], axis[0]))
    origin = axis_offset * axis_normal

    rel = H - origin
    longitudinal = rel @ axis
    transverse = rel @ axis_normal
    longitudinal_d = Hd @ axis
    transverse_d = Hd @ axis_normal

    margin2 = rod*rod - transverse*transverse
    if np.any(margin2 <= 1e-10):
        raise ValueError("rod closure")

    margin = np.sqrt(margin2)
    slider = longitudinal + margin
    slider_d = longitudinal_d - transverse*transverse_d/margin

    stroke = float(np.ptp(slider))
    if stroke <= 1e-10:
        raise ValueError("zero stroke")

    q = 1.0 - (slider - float(np.min(slider))) / stroke
    dq = -slider_d / stroke

    pos = float(np.sqrt(np.mean((q - target_q)**2)))
    der = float(np.sqrt(np.mean((dq - target_dq)**2)))
    motion = math.sqrt(pos*pos + 0.1*der*der)

    rod_cos = float(np.min(margin / rod))
    line = best_line_metrics(H, stroke, axis)
    clearance, inside_count = signed_triangle_clearance_origin(E, F, H)

    out = {
        "score": motion,
        "motion_score": motion,
        "position_rms": pos,
        "derivative_rms": der,
        "stroke_over_crank": stroke,
        "minimum_primary_transmission_sine": primary_sine,
        "minimum_secondary_transmission_sine": secondary_sine,
        "minimum_rod_axis_cosine": rod_cos,
        "EF_over_crank": lef,
        "GF_over_crank": lgf,
        "EH_over_EF": eh_ratio,
        "EH_over_crank": eh,
        "piston_rod_over_crank": rod,
        "piston_rod_over_stroke": rod / stroke,
        "crank_axis_to_EFH_clearance_over_crank": clearance,
        "crank_axis_inside_EFH_frames": inside_count,
        "primary": pdata,
        "parameters": dict(zip(NAMES, map(float, v))),
        "second_branch": second_branch,
    }
    out.update(line)
    return out


def feasible(r: dict, args) -> bool:
    return (
        r["stroke_over_crank"] >= args.stroke_floor
        and r["minimum_primary_transmission_sine"] >= args.primary_sine_floor
        and r["minimum_secondary_transmission_sine"] >= args.secondary_sine_floor
        and r["minimum_rod_axis_cosine"] >= args.rod_cos_floor
        and r["H_best_line_rms_over_stroke"] <= args.h_line_rms_max
        and r["H_best_line_axis_angle_deg"] <= args.h_axis_angle_max
        and r["crank_axis_to_EFH_clearance_over_crank"] >= args.crank_clearance_floor
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source_json", type=Path)
    ap.add_argument("--seed-restart", type=int, default=2)
    ap.add_argument("--seed-branch", type=int, choices=[-1, 1], default=1)

    ap.add_argument("--iterations", type=int, default=300)
    ap.add_argument("--population-size", type=int, default=14)
    ap.add_argument("--restarts", type=int, default=3)
    ap.add_argument("--seed", type=int, default=2600)

    ap.add_argument("--stroke-floor", type=float, default=2.5)
    ap.add_argument("--primary-sine-floor", type=float, default=0.30)
    ap.add_argument("--secondary-sine-floor", type=float, default=0.30)
    ap.add_argument("--rod-cos-floor", type=float, default=0.95)
    ap.add_argument("--maximum-EH", type=float, default=4.5)

    ap.add_argument("--h-line-rms-max", type=float, default=0.07)
    ap.add_argument("--h-axis-angle-max", type=float, default=5.0)
    ap.add_argument("--crank-clearance-floor", type=float, default=0.50)

    # Local neighborhood controls.
    ap.add_argument("--primary-length-fraction", type=float, default=0.20)
    ap.add_argument("--primary-point-fraction", type=float, default=0.35)
    ap.add_argument("--primary-phase-span-deg", type=float, default=20.0)
    ap.add_argument("--pivot-span", type=float, default=1.25)
    ap.add_argument("--secondary-length-fraction", type=float, default=0.25)
    ap.add_argument("--max-EF", type=float, default=3.0)
    ap.add_argument("--max-GF", type=float, default=5.5)
    ap.add_argument("--h-point-fraction", type=float, default=0.30)
    ap.add_argument("--rod-fraction", type=float, default=0.40)
    ap.add_argument("--axis-offset-span", type=float, default=1.5)
    ap.add_argument("--axis-angle-span-deg", type=float, default=12.0)

    ap.add_argument(
        "--output",
        type=Path,
        default=ROOT / "outputs" / "small_sixbar_stage2f.json",
    )
    args = ap.parse_args()

    seed_best = load_seed(args.source_json, args.seed_restart, args.seed_branch)
    x0 = seed_vector(seed_best)
    bounds = local_bounds(x0, args)

    raw = np.genfromtxt(TARGET, delimiter=",", names=True)
    if raw["theta_rad"][-1] >= 2.0*math.pi - 1e-9:
        raw = raw[:-1]

    coarse = raw[::12]
    theta = coarse["theta_rad"]
    tq = coarse["small_fraction_0_1"]
    tdq = coarse["small_dq_dtheta_per_rad"]

    dense_theta = raw["theta_rad"]
    dense_q = raw["small_fraction_0_1"]
    dense_dq = raw["small_dq_dtheta_per_rad"]

    # Validate the seed under the requested Stage 2F constraints.
    try:
        seed_dense = evaluate(
            x0,
            args.seed_branch,
            dense_theta,
            dense_q,
            dense_dq,
            args.maximum_EH,
        )
    except Exception:
        seed_dense = None

    report = {
        "description": (
            "Stage 2F local refinement around a selected Stage 2E candidate."
        ),
        "source_json": str(args.source_json),
        "seed_restart": args.seed_restart,
        "second_branch": args.seed_branch,
        "length_unit": "crank_radius",
        "bounds": dict(zip(NAMES, bounds)),
        "constraints": {
            "stroke_over_crank_minimum": args.stroke_floor,
            "primary_sine_minimum": args.primary_sine_floor,
            "secondary_sine_minimum": args.secondary_sine_floor,
            "rod_axis_cosine_minimum": args.rod_cos_floor,
            "EH_over_crank_maximum": args.maximum_EH,
            "H_best_line_rms_over_stroke_maximum": args.h_line_rms_max,
            "H_best_line_axis_angle_deg_maximum": args.h_axis_angle_max,
            "crank_axis_to_EFH_clearance_over_crank_minimum": args.crank_clearance_floor,
        },
        "seed_candidate": seed_dense,
        "runs": [],
    }

    start = time.monotonic()

    for restart in range(args.restarts):
        archive = []
        run_seed = args.seed + restart

        def objective(v):
            try:
                r = evaluate(
                    v,
                    args.seed_branch,
                    theta,
                    tq,
                    tdq,
                    args.maximum_EH,
                )
            except (ValueError, FloatingPointError):
                return 1e4

            violation = 0.0
            violation += max(
                0.0,
                args.stroke_floor - r["stroke_over_crank"],
            ) / max(args.stroke_floor, 1e-9)

            violation += max(
                0.0,
                args.primary_sine_floor - r["minimum_primary_transmission_sine"],
            ) / max(args.primary_sine_floor, 1e-9)

            violation += max(
                0.0,
                args.secondary_sine_floor - r["minimum_secondary_transmission_sine"],
            ) / max(args.secondary_sine_floor, 1e-9)

            violation += max(
                0.0,
                args.rod_cos_floor - r["minimum_rod_axis_cosine"],
            ) / max(1.0 - args.rod_cos_floor, 1e-3)

            violation += max(
                0.0,
                r["H_best_line_rms_over_stroke"] - args.h_line_rms_max,
            ) / max(args.h_line_rms_max, 1e-9)

            violation += max(
                0.0,
                r["H_best_line_axis_angle_deg"] - args.h_axis_angle_max,
            ) / max(args.h_axis_angle_max, 1e-9)

            violation += max(
                0.0,
                args.crank_clearance_floor
                - r["crank_axis_to_EFH_clearance_over_crank"],
            ) / max(args.crank_clearance_floor, 0.1)

            if violation > 0.0:
                return 10.0 + 100.0*violation

            if not archive or r["score"] < archive[-1]["score"]:
                archive.append(r)

            return float(r["score"])

        opt = differential_evolution(
            objective,
            bounds,
            seed=run_seed,
            popsize=args.population_size,
            maxiter=args.iterations,
            polish=True,
            tol=1e-8,
            updating="immediate",
            workers=1,
            x0=x0,
        )

        dense_candidates = []

        # Always include optimizer x and best archived improvements.
        candidate_vectors = [np.array(opt.x, dtype=float)]
        for c in sorted(archive, key=lambda x: x["score"])[:60]:
            candidate_vectors.append(
                np.array([c["parameters"][name] for name in NAMES], dtype=float)
            )

        for v in candidate_vectors:
            try:
                d = evaluate(
                    v,
                    args.seed_branch,
                    dense_theta,
                    dense_q,
                    dense_dq,
                    args.maximum_EH,
                )
            except (ValueError, FloatingPointError):
                continue

            if feasible(d, args):
                d["feasible"] = True
                dense_candidates.append(d)

        best = (
            min(dense_candidates, key=lambda x: x["score"])
            if dense_candidates else None
        )

        report["runs"].append({
            "restart": restart,
            "seed": run_seed,
            "second_branch": args.seed_branch,
            "optimizer_success": bool(opt.success),
            "optimizer_message": str(opt.message),
            "function_evaluations": int(opt.nfev),
            "best_dense": best,
        })

        candidates = [
            r["best_dense"]
            for r in report["runs"]
            if r["best_dense"] is not None
        ]
        if seed_dense is not None and feasible(seed_dense, args):
            candidates.append(seed_dense)

        report["best"] = (
            min(candidates, key=lambda x: x["score"])
            if candidates else None
        )
        report["elapsed_seconds"] = time.monotonic() - start

        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")

        if best is None:
            print(
                f"restart {restart}: no dense feasible candidate",
                flush=True,
            )
        else:
            print(
                f"restart {restart}: "
                f"pos={100*best['position_rms']:.4f}% "
                f"dRMS={best['derivative_rms']:.6f} "
                f"stroke/r={best['stroke_over_crank']:.3f} "
                f"HlineRMS={best['H_best_line_rms_over_stroke']:.3f} "
                f"Haxis={best['H_best_line_axis_angle_deg']:.2f}deg "
                f"clearA={best['crank_axis_to_EFH_clearance_over_crank']:.3f} "
                f"EF={best['EF_over_crank']:.3f} "
                f"GF={best['GF_over_crank']:.3f} "
                f"EH={best['EH_over_crank']:.3f}",
                flush=True,
            )

    print()
    b = report.get("best")
    if b:
        print("BEST DENSE CANDIDATE")
        print(f"position RMS = {100*b['position_rms']:.6f}%")
        print(f"derivative RMS = {b['derivative_rms']:.9f} rad^-1")
        print(f"stroke / crank = {b['stroke_over_crank']:.6f}")
        print(f"H best-line RMS / stroke = {b['H_best_line_rms_over_stroke']:.6f}")
        print(f"H best-line max / stroke = {b['H_best_line_max_over_stroke']:.6f}")
        print(f"H best-line / cylinder angle = {b['H_best_line_axis_angle_deg']:.6f} deg")
        print(f"crank-axis clearance / crank = {b['crank_axis_to_EFH_clearance_over_crank']:.6f}")
        print(f"EF / crank = {b['EF_over_crank']:.6f}")
        print(f"GF / crank = {b['GF_over_crank']:.6f}")
        print(f"EH / crank = {b['EH_over_crank']:.6f}")
        print(f"piston rod / stroke = {b['piston_rod_over_stroke']:.6f}")
        print(f"primary transmission sine = {b['minimum_primary_transmission_sine']:.6f}")
        print(f"secondary transmission sine = {b['minimum_secondary_transmission_sine']:.6f}")
        print(f"rod-axis cosine = {b['minimum_rod_axis_cosine']:.6f}")
        print("primary =", json.dumps(b["primary"], indent=2))
    else:
        print("No dense feasible candidate found.")

    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
