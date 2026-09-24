#!/usr/bin/env python3
"""Local 15-D family polish for selected V4.1 LARGE six-bar children.

Purpose
-------
Start from already synthesized six-bar children and allow SMALL continuous
changes to both the primary four-bar and downstream dyad, while preserving the
discrete primary and secondary assembly branches.

This is not a new global 15-D search. It is a trust-region family polish.

Default families:
    primary ranks 1,4,12,50

Default local freedom:
Primary:
    AD, BC, CD       +/- 5 %
    E_along, E_normal +/- max(0.20 crank, 8 % of |value|)
    primary phase    +/- 8 deg

Downstream:
    Gx, Gy           +/- 1.0 crank
    EF, GF           +/- 15 %
    H along/normal   +/- 0.25
    piston rod       +/- 20 %
    slider offset    +/- 1.0 crank
    slider angle     +/- 10 deg

The exact current child is injected as x0. Constraints and the full-cycle V4.1
motion objective are inherited from search_large_second_dyad_v41_panel_3952.py.

Typical:
    PYTHONPATH=src python3 \
      examples/polish_large_sixbar_v41_families_3952.py

Output:
    outputs/large_sixbar_v41_family_polish_3952.json
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import differential_evolution

import search_large_second_dyad_v41_panel_3952 as v41


ROOT = Path.cwd()
SOURCE = ROOT / "outputs" / "large_second_dyad_v41_panel_3952.json"
TARGET = (
    ROOT
    / "outputs"
    / "motor_hybrid_c2_15p_260k"
    / "candidate_3952_motion_target.csv"
)
OUTPUT = ROOT / "outputs" / "large_sixbar_v41_family_polish_3952.json"

NAMES = v41.base.stage.NAMES


def load_json_relaxed(path: Path):
    return json.loads(path.read_text(encoding="utf-8").replace("NaN", "null"))


def wrap_pi(a: float) -> float:
    return (a + math.pi) % (2.0 * math.pi) - math.pi


def child_vector(best: dict) -> np.ndarray:
    p = best["parameters"]
    return np.array([float(p[name]) for name in NAMES], dtype=float)


def signed_span(v: float, rel: float, minimum: float):
    span = max(abs(v) * rel, minimum)
    return (v - span, v + span)


def positive_span(v: float, rel: float, floor: float = 1e-3):
    lo = max(floor, v * (1.0 - rel))
    hi = max(lo + 1e-6, v * (1.0 + rel))
    return (lo, hi)


def local_bounds(x0: np.ndarray, args):
    b = []

    # Primary 0..5
    b.append(positive_span(x0[0], args.primary_length_fraction, 0.50))
    b.append(positive_span(x0[1], args.primary_length_fraction, 0.50))
    b.append(positive_span(x0[2], args.primary_length_fraction, 0.35))

    b.append(signed_span(
        x0[3], args.primary_point_fraction, args.primary_point_min_span
    ))
    b.append(signed_span(
        x0[4], args.primary_point_fraction, args.primary_point_min_span
    ))

    phase_span = math.radians(args.primary_phase_span_deg)
    b.append((x0[5] - phase_span, x0[5] + phase_span))

    # Downstream 6..14
    b.append((x0[6] - args.pivot_span, x0[6] + args.pivot_span))
    b.append((x0[7] - args.pivot_span, x0[7] + args.pivot_span))

    b.append(positive_span(x0[8], args.secondary_length_fraction, 0.25))
    b.append(positive_span(x0[9], args.secondary_length_fraction, 0.25))

    b.append((x0[10] - args.h_ratio_span, x0[10] + args.h_ratio_span))
    b.append((x0[11] - args.h_ratio_span, x0[11] + args.h_ratio_span))

    b.append(positive_span(x0[12], args.rod_fraction, 0.75))

    b.append((x0[13] - args.axis_offset_span, x0[13] + args.axis_offset_span))

    angle_span = math.radians(args.axis_angle_span_deg)
    b.append((x0[14] - angle_span, x0[14] + angle_span))

    return b


def primary_closure_penalty(v: np.ndarray) -> float:
    g, c, r = map(float, v[:3])
    inner_margin = abs(g - 1.0) - abs(c - r)
    outer_margin = c + r - (g + 1.0)

    violation = 0.0
    if inner_margin <= 1e-10:
        violation += -inner_margin + 1e-6
    if outer_margin <= 1e-10:
        violation += -outer_margin + 1e-6
    return violation


def evaluate_full(
    v: np.ndarray,
    *,
    primary_branch: int,
    second_branch: int,
    theta: np.ndarray,
    tq: np.ndarray,
    tdq: np.ndarray,
    mask: np.ndarray,
    velocity_weight: float,
):
    vv = np.asarray(v, dtype=float).copy()
    vv[5] = wrap_pi(float(vv[5]))
    vv[14] = wrap_pi(float(vv[14]))

    pv = vv[:6]
    downstream = vv[6:]

    E, Ed, psine, pdata = v41.base.stage.primary(
        theta,
        pv,
        primary_branch,
    )

    result = v41.evaluate_v41(
        downstream,
        E=E,
        Ed=Ed,
        psine=psine,
        pdata=pdata,
        pv=pv,
        second_branch=second_branch,
        tq=tq,
        tdq=tdq,
        mask=mask,
        velocity_weight=velocity_weight,
    )
    return result


def violation_with_primary(result: dict, args) -> float:
    violation = v41.violation_v41(result, args)
    violation += max(
        0.0,
        args.primary_sine_floor
        - result["minimum_primary_transmission_sine"],
    ) / args.primary_sine_floor
    return float(violation)


def save(path: Path, obj: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, default=SOURCE)
    ap.add_argument("--target", type=Path, default=TARGET)
    ap.add_argument("--output", type=Path, default=OUTPUT)
    ap.add_argument("--ranks", default="1,4,12,50")

    ap.add_argument("--restarts", type=int, default=2)
    ap.add_argument("--generations", type=int, default=260)
    ap.add_argument("--population-size", type=int, default=10)
    ap.add_argument("--seed", type=int, default=3956200)
    ap.add_argument("--coarse-stride", type=int, default=4)
    ap.add_argument("--archive-per-run", type=int, default=40)

    # Trust region: intentionally modest.
    ap.add_argument("--primary-length-fraction", type=float, default=0.05)
    ap.add_argument("--primary-point-fraction", type=float, default=0.08)
    ap.add_argument("--primary-point-min-span", type=float, default=0.20)
    ap.add_argument("--primary-phase-span-deg", type=float, default=8.0)

    ap.add_argument("--pivot-span", type=float, default=1.0)
    ap.add_argument("--secondary-length-fraction", type=float, default=0.15)
    ap.add_argument("--h-ratio-span", type=float, default=0.25)
    ap.add_argument("--rod-fraction", type=float, default=0.20)
    ap.add_argument("--axis-offset-span", type=float, default=1.0)
    ap.add_argument("--axis-angle-span-deg", type=float, default=10.0)

    # Same motion/feasibility settings as V4.1 panel.
    ap.add_argument("--noise-fronts", type=int, default=6)
    ap.add_argument("--velocity-weight", type=float, default=0.0005)
    ap.add_argument("--stroke-floor", type=float, default=1.0)
    ap.add_argument("--stroke-ceiling", type=float, default=3.0)
    ap.add_argument("--primary-sine-floor", type=float, default=0.30)
    ap.add_argument("--secondary-sine-floor", type=float, default=0.30)
    ap.add_argument("--rod-cos-floor", type=float, default=0.95)
    ap.add_argument("--maximum-EH", type=float, default=7.0)
    ap.add_argument("--crank-clearance-floor", type=float, default=0.50)
    ap.add_argument("--h-lateral-rms-max", type=float, default=0.20)
    ap.add_argument("--h-lateral-span-max", type=float, default=0.65)

    args = ap.parse_args()

    source = load_json_relaxed(args.source)
    by_rank = {
        int(f["primary_rank"]): f
        for f in source.get("families", [])
    }

    ranks = [int(x) for x in args.ranks.split(",") if x.strip()]
    missing = [
        rank for rank in ranks
        if rank not in by_rank or not by_rank[rank].get("best")
    ]
    if missing:
        raise SystemExit(f"Missing feasible source children for ranks: {missing}")

    raw = np.genfromtxt(args.target, delimiter=",", names=True)
    if raw["theta_rad"][-1] >= 2.0 * math.pi - 1e-9:
        raw = raw[:-1]

    band = v41.base.motor3952._detect_noise_bands(
        args.target,
        args.noise_fronts,
    )["large"]

    dense_theta = np.asarray(raw["theta_rad"], float)
    dense_q = np.asarray(raw["large_fraction_0_1"], float)
    dense_dq = np.asarray(raw["large_dq_dtheta_per_rad"], float)
    dense_mask = ~v41.base.motor3952._theta_band_mask(
        dense_theta,
        band,
    )

    coarse = raw[::args.coarse_stride]
    theta = np.asarray(coarse["theta_rad"], float)
    tq = np.asarray(coarse["large_fraction_0_1"], float)
    tdq = np.asarray(coarse["large_dq_dtheta_per_rad"], float)
    mask = ~v41.base.motor3952._theta_band_mask(theta, band)

    report = {
        "description": (
            "Local 15-D family polish of selected V4.1 LARGE six-bar children. "
            "Discrete branches fixed; small trust region around each existing child."
        ),
        "source": str(args.source),
        "target": str(args.target),
        "ranks": ranks,
        "trust_region": {
            "primary_length_fraction": args.primary_length_fraction,
            "primary_point_fraction": args.primary_point_fraction,
            "primary_point_min_span": args.primary_point_min_span,
            "primary_phase_span_deg": args.primary_phase_span_deg,
            "pivot_span": args.pivot_span,
            "secondary_length_fraction": args.secondary_length_fraction,
            "h_ratio_span": args.h_ratio_span,
            "rod_fraction": args.rod_fraction,
            "axis_offset_span": args.axis_offset_span,
            "axis_angle_span_deg": args.axis_angle_span_deg,
        },
        "families": [],
        "best": None,
    }

    wall0 = time.time()
    all_best = []

    for family_index, rank in enumerate(ranks):
        source_family = by_rank[rank]
        seed_best = source_family["best"]
        x0 = child_vector(seed_best)
        primary_branch = int(seed_best["primary"]["assembly_branch"])
        second_branch = int(seed_best["second_branch"])
        bounds = local_bounds(x0, args)

        family = {
            "primary_rank": rank,
            "primary_branch": primary_branch,
            "second_branch": second_branch,
            "seed": seed_best,
            "bounds": dict(zip(NAMES, bounds)),
            "runs": [],
            "best": seed_best,
        }

        print(
            f"\nFAMILY rank={rank} "
            f"primary={primary_branch:+d} second={second_branch:+d} "
            f"seed pos={100*seed_best['position_rms_full_cycle']:.5f}%",
            flush=True,
        )

        for restart in range(args.restarts):
            run_seed = args.seed + 1000 * family_index + restart
            archive = []

            def objective(v):
                v = np.asarray(v, dtype=float)

                cp = primary_closure_penalty(v)
                if cp > 0.0:
                    return 10.0 + 100.0 * cp

                try:
                    result = evaluate_full(
                        v,
                        primary_branch=primary_branch,
                        second_branch=second_branch,
                        theta=theta,
                        tq=tq,
                        tdq=tdq,
                        mask=mask,
                        velocity_weight=args.velocity_weight,
                    )
                except (
                    ValueError,
                    FloatingPointError,
                    OverflowError,
                    np.linalg.LinAlgError,
                ):
                    return 1000.0

                violation = violation_with_primary(result, args)
                if violation > 0.0:
                    return 10.0 + 100.0 * violation

                archive.append((float(result["score"]), v.copy()))
                return float(result["score"])

            t0 = time.time()
            p0 = time.process_time()

            opt = differential_evolution(
                objective,
                bounds,
                seed=run_seed,
                popsize=args.population_size,
                maxiter=args.generations,
                polish=True,
                tol=1e-8,
                updating="immediate",
                workers=1,
                init="latinhypercube",
                x0=x0,
            )

            vectors = [np.asarray(opt.x, float), x0.copy()]

            if getattr(opt, "population", None) is not None:
                order = np.argsort(np.asarray(opt.population_energies, float))
                vectors += [
                    np.asarray(opt.population[int(i)], float)
                    for i in order[:args.archive_per_run]
                ]

            vectors += [
                v
                for _, v in sorted(archive, key=lambda z: z[0])[
                    :args.archive_per_run
                ]
            ]

            dense_results = []
            seen = set()

            for v in vectors:
                key = tuple(round(float(x), 10) for x in v)
                if key in seen:
                    continue
                seen.add(key)

                try:
                    result = evaluate_full(
                        v,
                        primary_branch=primary_branch,
                        second_branch=second_branch,
                        theta=dense_theta,
                        tq=dense_q,
                        tdq=dense_dq,
                        mask=dense_mask,
                        velocity_weight=args.velocity_weight,
                    )
                except (
                    ValueError,
                    FloatingPointError,
                    OverflowError,
                    np.linalg.LinAlgError,
                ):
                    continue

                if violation_with_primary(result, args) <= 0.0:
                    result["feasible"] = True
                    result["source_primary_rank"] = rank
                    dense_results.append(result)

            best = (
                min(dense_results, key=lambda r: r["score"])
                if dense_results else None
            )

            family["runs"].append({
                "restart": restart,
                "seed": run_seed,
                "function_evaluations": int(opt.nfev),
                "optimizer_success": bool(opt.success),
                "optimizer_message": str(opt.message),
                "wall_seconds": time.time() - t0,
                "cpu_seconds": time.process_time() - p0,
                "best_dense": best,
            })

            candidates = [seed_best]
            candidates += [
                run["best_dense"]
                for run in family["runs"]
                if run.get("best_dense") is not None
            ]
            family["best"] = min(candidates, key=lambda r: r["score"])

            if best:
                print(
                    f"  restart={restart}: "
                    f"pos={100*best['position_rms_full_cycle']:.5f}% "
                    f"HP={100*best['position_rms_hp_band']:.3f}% "
                    f"nonHP={100*best['position_rms_non_hp']:.5f}% "
                    f"p-sine={best['minimum_primary_transmission_sine']:.3f} "
                    f"s-sine={best['minimum_secondary_transmission_sine']:.3f} "
                    f"rod/stroke={best['piston_rod_over_stroke']:.3f}",
                    flush=True,
                )
            else:
                print(f"  restart={restart}: no dense feasible candidate", flush=True)

            temp = report["families"] + [family]
            temp_best = [
                f["best"] for f in temp
                if f.get("best") is not None
            ]
            report["best"] = (
                min(temp_best, key=lambda r: r["score"])
                if temp_best else None
            )
            report["elapsed_wall_hours"] = (time.time() - wall0) / 3600.0

            old = report["families"]
            report["families"] = temp
            save(args.output, report)
            report["families"] = old

        report["families"].append(family)
        all_best.append(family["best"])
        report["best"] = min(all_best, key=lambda r: r["score"])
        report["elapsed_wall_hours"] = (time.time() - wall0) / 3600.0
        save(args.output, report)

    print("\nFAMILY POLISH SUMMARY")
    for family in report["families"]:
        seed = family["seed"]
        best = family["best"]
        improvement_pp = 100.0 * (
            seed["position_rms_full_cycle"]
            - best["position_rms_full_cycle"]
        )
        print(
            f"rank {family['primary_rank']:2d}: "
            f"{100*seed['position_rms_full_cycle']:.5f}% -> "
            f"{100*best['position_rms_full_cycle']:.5f}% "
            f"(gain {improvement_pp:+.5f} pp)"
        )

    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
