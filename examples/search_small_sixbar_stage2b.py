#!/usr/bin/env python3
"""Stage 2B: fixed primary four-bar + second RR dyad + free point H + finite piston rod."""

from __future__ import annotations
import argparse, json, math, time
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution

ROOT = Path.cwd()
TARGET = ROOT / "outputs" / "motor_champion_motion_target.csv"

PRIMARY = {
    "ground": 2.045193376592419,
    "coupler": 2.1152136168766678,
    "rocker": 1.8202134329709354,
    "branch": -1,
    "e_along": 1.873315626010214,
    "e_normal": -1.5362965718510164,
    "phase": -3.0863675667472057,
}

PARAMETER_NAMES = [
    "second_pivot_x", "second_pivot_y",
    "link_EF", "link_GF",
    "H_along_over_EF", "H_normal_over_EF",
    "piston_rod", "slider_axis_offset", "slider_axis_angle",
]

BOUNDS = [
    (-8.0, 8.0), (-8.0, 8.0),
    (0.4, 5.0), (0.4, 5.0),
    (-1.5, 1.5), (-1.5, 1.5),
    (2.0, 30.0), (-10.0, 10.0), (-math.pi, math.pi),
]


def primary_motion(theta):
    g, c, r = PRIMARY["ground"], PRIMARY["coupler"], PRIMARY["rocker"]
    branch = PRIMARY["branch"]
    angle = theta + PRIMARY["phase"]

    B = np.column_stack((np.cos(angle), np.sin(angle)))
    Bd = np.column_stack((-np.sin(angle), np.cos(angle)))
    D = np.array((g, 0.0))

    delta = D - B
    distance = np.linalg.norm(delta, axis=1)
    if np.any(distance <= 1e-12):
        raise ValueError("Primary centers coincide")

    direction = delta / distance[:, None]
    along = (c*c - r*r + distance*distance) / (2.0*distance)
    height2 = c*c - along*along
    if np.any(height2 <= 1e-12):
        raise ValueError("Primary does not close")

    height = np.sqrt(height2)
    left = np.column_stack((-direction[:, 1], direction[:, 0]))
    C = B + along[:, None]*direction + branch*height[:, None]*left

    coupler = C - B
    rocker = C - D
    det = coupler[:, 0]*rocker[:, 1] - coupler[:, 1]*rocker[:, 0]
    if np.any(np.abs(det) <= 1e-10):
        raise ValueError("Primary singular")

    rhs = np.sum(coupler * Bd, axis=1)
    Cd = np.column_stack((rhs*rocker[:, 1]/det, -rhs*rocker[:, 0]/det))

    u = coupler / c
    ud = (Cd - Bd) / c
    n = np.column_stack((-u[:, 1], u[:, 0]))
    nd = np.column_stack((-ud[:, 1], ud[:, 0]))

    E = B + PRIMARY["e_along"]*u + PRIMARY["e_normal"]*n
    Ed = Bd + PRIMARY["e_along"]*ud + PRIMARY["e_normal"]*nd

    min_sine = float(np.min(np.abs(det)/(c*r)))
    return E, Ed, min_sine


def evaluate(vector, second_branch, E, Ed, target_q, target_dq, lateral_weight):
    gx, gy, lef, lgf, ha, hn, rod, axis_offset, axis_angle = map(float, vector)
    G = np.array((gx, gy))

    # Second RR dyad E-F-G.
    delta = G - E
    distance = np.linalg.norm(delta, axis=1)
    if np.any(distance <= 1e-10):
        raise ValueError("Second centers coincide")

    along = (lef*lef - lgf*lgf + distance*distance) / (2.0*distance)
    height2 = lef*lef - along*along
    if np.any(height2 <= 1e-10):
        raise ValueError("Second dyad does not close")

    height = np.sqrt(height2)
    direction = delta / distance[:, None]
    left = np.column_stack((-direction[:, 1], direction[:, 0]))
    F = E + along[:, None]*direction + second_branch*height[:, None]*left

    ef = F - E
    gf = F - G
    det = ef[:, 0]*gf[:, 1] - ef[:, 1]*gf[:, 0]
    if np.any(np.abs(det) <= 1e-10):
        raise ValueError("Second dyad singular")

    rhs = np.sum(ef * Ed, axis=1)
    Fd = np.column_stack((rhs*gf[:, 1]/det, -rhs*gf[:, 0]/det))
    secondary_sine = float(np.min(np.abs(det)/(lef*lgf)))

    # Free rigid point H on link EF.
    u = ef / lef
    ud = (Fd - Ed) / lef
    n = np.column_stack((-u[:, 1], u[:, 0]))
    nd = np.column_stack((-ud[:, 1], ud[:, 0]))

    H = E + (ha*lef)*u + (hn*lef)*n
    Hd = Ed + (ha*lef)*ud + (hn*lef)*nd

    # Finite piston rod.
    axis = np.array((math.cos(axis_angle), math.sin(axis_angle)))
    normal = np.array((-axis[1], axis[0]))
    origin = axis_offset * normal

    rel = H - origin
    longitudinal = rel @ axis
    transverse = rel @ normal
    longitudinal_d = Hd @ axis
    transverse_d = Hd @ normal

    margin2 = rod*rod - transverse*transverse
    if np.any(margin2 <= 1e-10):
        raise ValueError("Piston rod cannot close")
    margin = np.sqrt(margin2)

    # Physical convention: piston on +axis side, working chamber beyond piston.
    slider = longitudinal + margin
    slider_d = longitudinal_d - transverse*transverse_d/margin

    stroke = float(np.ptp(slider))
    if stroke <= 1e-10:
        raise ValueError("Zero stroke")

    q = 1.0 - (slider - float(np.min(slider))) / stroke
    dq = -slider_d / stroke

    pos = float(np.sqrt(np.mean((q-target_q)**2)))
    der = float(np.sqrt(np.mean((dq-target_dq)**2)))
    motion_score = math.sqrt(pos*pos + 0.1*der*der)

    rod_cos = float(np.min(margin/rod))
    h_lat = float(np.ptp(H @ normal))
    h_lat_ratio = h_lat / stroke
    score = motion_score + lateral_weight*h_lat_ratio

    return {
        "score": score,
        "motion_score": motion_score,
        "position_rms": pos,
        "derivative_rms": der,
        "stroke_over_crank": stroke,
        "H_transverse_span_over_crank": h_lat,
        "H_lateral_span_over_piston_stroke": h_lat_ratio,
        "minimum_secondary_transmission_sine": secondary_sine,
        "minimum_rod_axis_cosine": rod_cos,
        "parameters": dict(zip(PARAMETER_NAMES, map(float, vector))),
        "second_branch": second_branch,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iterations", type=int, default=200)
    ap.add_argument("--population-size", type=int, default=10)
    ap.add_argument("--restarts", type=int, default=1)
    ap.add_argument("--seed", type=int, default=2200)
    ap.add_argument("--stroke-floor", type=float, default=1.5)
    ap.add_argument("--secondary-sine-floor", type=float, default=0.30)
    ap.add_argument("--rod-cos-floor", type=float, default=0.95)
    ap.add_argument("--lateral-weight", type=float, default=0.002)
    ap.add_argument("--lateral-ratio-max", type=float, default=None)
    ap.add_argument("--output", type=Path, default=ROOT/"outputs"/"small_sixbar_stage2b.json")
    args = ap.parse_args()

    raw = np.genfromtxt(TARGET, delimiter=",", names=True)
    if raw["theta_rad"][-1] >= 2*math.pi - 1e-9:
        raw = raw[:-1]

    dense_theta = raw["theta_rad"]
    dense_q = raw["small_fraction_0_1"]
    dense_dq = raw["small_dq_dtheta_per_rad"]

    coarse = raw[::12]
    theta = coarse["theta_rad"]
    tq = coarse["small_fraction_0_1"]
    tdq = coarse["small_dq_dtheta_per_rad"]

    E, Ed, _ = primary_motion(theta)
    dense_E, dense_Ed, primary_sine = primary_motion(dense_theta)

    report = {
        "description": "Stage 2B: frozen primary four-bar + second RR dyad + free point H + finite piston rod",
        "length_unit": "crank_radius",
        "primary": PRIMARY,
        "primary_minimum_transmission_sine": primary_sine,
        "objective": {
            "motion": "sqrt(position_RMS^2 + 0.1 * derivative_RMS^2)",
            "straightening": "motion_score + lateral_weight * (H lateral span / piston stroke)",
            "lateral_weight": args.lateral_weight,
        },
        "constraints": {
            "stroke_over_crank_minimum": args.stroke_floor,
            "minimum_secondary_transmission_sine": args.secondary_sine_floor,
            "minimum_rod_axis_cosine": args.rod_cos_floor,
            "H_lateral_span_over_piston_stroke_maximum": args.lateral_ratio_max,
        },
        "bounds": dict(zip(PARAMETER_NAMES, BOUNDS)),
        "runs": [],
    }

    start = time.monotonic()
    run_index = 0

    for restart in range(args.restarts):
        for second_branch in (+1, -1):
            archive = []
            seed = args.seed + run_index
            run_index += 1

            def objective(v):
                try:
                    r = evaluate(v, second_branch, E, Ed, tq, tdq, args.lateral_weight)
                except (ValueError, FloatingPointError):
                    return 1e4

                violation = 0.0
                violation += max(0.0, args.stroke_floor - r["stroke_over_crank"])
                violation += max(0.0, args.secondary_sine_floor - r["minimum_secondary_transmission_sine"])
                violation += max(0.0, args.rod_cos_floor - r["minimum_rod_axis_cosine"])
                if args.lateral_ratio_max is not None:
                    violation += max(0.0, r["H_lateral_span_over_piston_stroke"] - args.lateral_ratio_max)

                if violation:
                    return 10.0 + 100.0*violation

                if not archive or r["score"] < archive[-1]["score"]:
                    archive.append(r)
                return r["score"]

            opt = differential_evolution(
                objective, BOUNDS, seed=seed,
                popsize=args.population_size,
                maxiter=args.iterations,
                polish=True, tol=1e-8,
                updating="immediate", workers=1,
            )

            dense_candidates = []
            for c in sorted(archive, key=lambda x: x["score"])[:30]:
                v = np.array([c["parameters"][name] for name in PARAMETER_NAMES])
                try:
                    d = evaluate(v, second_branch, dense_E, dense_Ed, dense_q, dense_dq, args.lateral_weight)
                except (ValueError, FloatingPointError):
                    continue
                feasible = (
                    d["stroke_over_crank"] >= args.stroke_floor
                    and d["minimum_secondary_transmission_sine"] >= args.secondary_sine_floor
                    and d["minimum_rod_axis_cosine"] >= args.rod_cos_floor
                )
                if args.lateral_ratio_max is not None:
                    feasible &= d["H_lateral_span_over_piston_stroke"] <= args.lateral_ratio_max
                d["feasible"] = bool(feasible)
                if feasible:
                    dense_candidates.append(d)

            best = min(dense_candidates, key=lambda x: x["score"]) if dense_candidates else None

            report["runs"].append({
                "restart": restart,
                "seed": seed,
                "second_branch": second_branch,
                "optimizer_success": bool(opt.success),
                "optimizer_message": str(opt.message),
                "function_evaluations": int(opt.nfev),
                "best_dense": best,
            })

            valid = [r["best_dense"] for r in report["runs"] if r["best_dense"] is not None]
            report["best"] = min(valid, key=lambda x: x["score"]) if valid else None
            report["elapsed_seconds"] = time.monotonic() - start

            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(report, indent=2) + "\n")

            if best is None:
                print(f"restart {restart}, branch {second_branch:+d}: no dense feasible candidate", flush=True)
            else:
                print(
                    f"restart {restart}, branch {second_branch:+d}: "
                    f"pos={100*best['position_rms']:.4f}% "
                    f"dRMS={best['derivative_rms']:.6f} "
                    f"stroke/r={best['stroke_over_crank']:.3f} "
                    f"Hlat/stroke={best['H_lateral_span_over_piston_stroke']:.3f} "
                    f"rodcos={best['minimum_rod_axis_cosine']:.4f}",
                    flush=True,
                )

    print()
    if report.get("best"):
        b = report["best"]
        print("BEST DENSE CANDIDATE")
        print(f"position RMS = {100*b['position_rms']:.6f}%")
        print(f"derivative RMS = {b['derivative_rms']:.9f} rad^-1")
        print(f"stroke / crank = {b['stroke_over_crank']:.6f}")
        print(f"H lateral span / stroke = {b['H_lateral_span_over_piston_stroke']:.6f}")
        print(f"secondary transmission sine = {b['minimum_secondary_transmission_sine']:.6f}")
        print(f"rod-axis cosine = {b['minimum_rod_axis_cosine']:.6f}")
    else:
        print("No dense feasible candidate found.")

    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
