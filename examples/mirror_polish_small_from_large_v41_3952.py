#!/usr/bin/env python3
"""Mirror polished LARGE six-bar families and locally adapt them to SMALL, candidate 3952.

For each selected LARGE family:
1. take the polished LARGE six-bar;
2. apply the exact physical mirror used historically:
   - primary branch flips;
   - E normal coordinate flips;
   - primary phase changes sign;
   - second-pivot Y flips;
   - secondary branch flips;
   - H normal coordinate flips;
   - slider signed offset and axis angle reflect;
3. fit only the mirrored primary phase over a complete revolution against SMALL;
4. perform a modest 15-D trust-region polish while keeping both discrete
   assembly branches fixed;
5. export loader-compatible SMALL and LARGE mechanism JSON files as a pair.

The mirror is used as the canonical transform.  For the scalar piston law it is
equivalent to traversing the original mechanism with reversed crank direction,
but preserves the same visible crank rotation direction for S and L.

Default source:
  outputs/large_sixbar_v41_family_polish_3952.json

Default families:
  1,4,12,50

Outputs:
  outputs/sixbar_pairs_3952/mirror_family_summary.json
  outputs/sixbar_pairs_3952/rank_01_large.json
  outputs/sixbar_pairs_3952/rank_01_small.json
  ...

Typical run:
  PYTHONPATH=src python3 \
    examples/mirror_polish_small_from_large_v41_3952.py
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import differential_evolution, minimize_scalar

import search_large_sixbar_stageL1 as stage
import search_large_second_dyad_v41_panel_3952 as v41


ROOT = Path.cwd()
SOURCE = ROOT / "outputs" / "large_sixbar_v41_family_polish_3952.json"
TARGET = (
    ROOT
    / "outputs"
    / "motor_hybrid_c2_15p_260k"
    / "candidate_3952_motion_target.csv"
)
OUTDIR = ROOT / "outputs" / "sixbar_pairs_3952"
SUMMARY = OUTDIR / "mirror_family_summary.json"

NAMES = stage.NAMES


def load_json_relaxed(path: Path):
    return json.loads(path.read_text(encoding="utf-8").replace("NaN", "null"))


def wrap_pi(a: float) -> float:
    return (a + math.pi) % (2.0 * math.pi) - math.pi


def candidate_vector(best: dict) -> np.ndarray:
    p = best["parameters"]
    return np.array([float(p[name]) for name in NAMES], dtype=float)


def mirror_candidate(best: dict) -> tuple[np.ndarray, int, int]:
    """Exact physical mirror; stage.mirrored_vector is an involution."""
    source = candidate_vector(best)
    mirrored = stage.mirrored_vector(source)
    primary_branch = -int(best["primary"]["assembly_branch"])
    second_branch = -int(best["second_branch"])
    return mirrored, primary_branch, second_branch


def signed_span(v: float, rel: float, minimum: float):
    span = max(abs(v) * rel, minimum)
    return (v - span, v + span)


def positive_span(v: float, rel: float, floor: float):
    lo = max(floor, v * (1.0 - rel))
    hi = max(lo + 1e-6, v * (1.0 + rel))
    return (lo, hi)


def local_bounds(x0: np.ndarray, args):
    b = []

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
    inner = abs(g - 1.0) - abs(c - r)
    outer = c + r - (g + 1.0)
    violation = 0.0
    if inner <= 1e-10:
        violation += -inner + 1e-6
    if outer <= 1e-10:
        violation += -outer + 1e-6
    return float(violation)


def evaluate_full(
    v: np.ndarray,
    *,
    primary_branch: int,
    second_branch: int,
    theta: np.ndarray,
    target_q: np.ndarray,
    target_dq: np.ndarray,
    valid_mask: np.ndarray,
    velocity_weight: float,
):
    vv = np.asarray(v, dtype=float).copy()
    vv[5] = wrap_pi(float(vv[5]))
    vv[14] = wrap_pi(float(vv[14]))

    pv = vv[:6]
    downstream = vv[6:]

    E, Ed, primary_sine, pdata = stage.primary(
        theta, pv, primary_branch
    )

    result = v41.evaluate_v41(
        downstream,
        E=E,
        Ed=Ed,
        psine=primary_sine,
        pdata=pdata,
        pv=pv,
        second_branch=second_branch,
        tq=target_q,
        tdq=target_dq,
        mask=valid_mask,
        velocity_weight=velocity_weight,
    )
    return result


def total_violation(result: dict, args) -> float:
    value = v41.violation_v41(result, args)
    value += max(
        0.0,
        args.primary_sine_floor
        - result["minimum_primary_transmission_sine"],
    ) / args.primary_sine_floor
    return float(value)


def constrained_score(
    v,
    *,
    args,
    primary_branch,
    second_branch,
    theta,
    target_q,
    target_dq,
    valid_mask,
):
    cp = primary_closure_penalty(np.asarray(v, float))
    if cp > 0.0:
        return 10.0 + 100.0 * cp

    try:
        result = evaluate_full(
            v,
            primary_branch=primary_branch,
            second_branch=second_branch,
            theta=theta,
            target_q=target_q,
            target_dq=target_dq,
            valid_mask=valid_mask,
            velocity_weight=args.velocity_weight,
        )
    except (
        ValueError,
        FloatingPointError,
        OverflowError,
        np.linalg.LinAlgError,
    ):
        return 1000.0

    violation = total_violation(result, args)
    if violation > 0.0:
        return 10.0 + 100.0 * violation

    return float(result["score"])


def fit_phase(
    x: np.ndarray,
    *,
    args,
    primary_branch: int,
    second_branch: int,
    theta: np.ndarray,
    target_q: np.ndarray,
    target_dq: np.ndarray,
    valid_mask: np.ndarray,
):
    """Global phase scan, then bounded scalar refinement."""
    grid = np.linspace(-math.pi, math.pi, args.phase_scan_points, endpoint=False)
    scores = np.full(grid.shape, np.inf)

    for i, phase in enumerate(grid):
        v = np.array(x, dtype=float, copy=True)
        v[5] = phase
        scores[i] = constrained_score(
            v,
            args=args,
            primary_branch=primary_branch,
            second_branch=second_branch,
            theta=theta,
            target_q=target_q,
            target_dq=target_dq,
            valid_mask=valid_mask,
        )

    if not np.isfinite(scores).any():
        raise ValueError("Mirrored mechanism could not be evaluated at any phase.")

    i = int(np.argmin(scores))
    p0 = float(grid[i])
    step = 2.0 * math.pi / args.phase_scan_points

    def objective(p):
        v = np.array(x, dtype=float, copy=True)
        v[5] = wrap_pi(float(p))
        return constrained_score(
            v,
            args=args,
            primary_branch=primary_branch,
            second_branch=second_branch,
            theta=theta,
            target_q=target_q,
            target_dq=target_dq,
            valid_mask=valid_mask,
        )

    opt = minimize_scalar(
        objective,
        bounds=(p0 - 2.0 * step, p0 + 2.0 * step),
        method="bounded",
        options={"xatol": 1e-11},
    )

    y = np.array(x, dtype=float, copy=True)
    y[5] = wrap_pi(float(opt.x))
    return y


def loader_file(candidate: dict, description: str):
    return {
        "description": description,
        "length_unit": "crank_radius",
        "best": candidate,
    }


def save(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(obj, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, default=SOURCE)
    ap.add_argument("--target", type=Path, default=TARGET)
    ap.add_argument("--outdir", type=Path, default=OUTDIR)
    ap.add_argument("--summary", type=Path, default=SUMMARY)
    ap.add_argument("--ranks", default="1,4,12,50")

    ap.add_argument("--restarts", type=int, default=2)
    ap.add_argument("--generations", type=int, default=260)
    ap.add_argument("--population-size", type=int, default=10)
    ap.add_argument("--seed", type=int, default=3958200)
    ap.add_argument("--coarse-stride", type=int, default=4)
    ap.add_argument("--archive-per-run", type=int, default=40)
    ap.add_argument("--phase-scan-points", type=int, default=720)

    # Same modest trust region that proved useful for LARGE family polish.
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

    # Motion objective / mechanical envelope.
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
        if f.get("best") is not None
    }

    ranks = [int(x) for x in args.ranks.split(",") if x.strip()]
    missing = [rank for rank in ranks if rank not in by_rank]
    if missing:
        raise SystemExit(f"Missing polished LARGE families: {missing}")

    raw = np.genfromtxt(args.target, delimiter=",", names=True)
    if raw["theta_rad"][-1] >= 2.0 * math.pi - 1e-9:
        raw = raw[:-1]

    dense_theta = np.asarray(raw["theta_rad"], float)
    dense_q = np.asarray(raw["small_fraction_0_1"], float)
    dense_dq = np.asarray(raw["small_dq_dtheta_per_rad"], float)

    noise_band = v41.base.motor3952._detect_noise_bands(
        args.target, args.noise_fronts
    )["small"]

    dense_valid = ~v41.base.motor3952._theta_band_mask(
        dense_theta, noise_band
    )

    coarse = raw[::args.coarse_stride]
    theta = np.asarray(coarse["theta_rad"], float)
    target_q = np.asarray(coarse["small_fraction_0_1"], float)
    target_dq = np.asarray(coarse["small_dq_dtheta_per_rad"], float)
    valid = ~v41.base.motor3952._theta_band_mask(theta, noise_band)

    report = {
        "description": (
            "Polished LARGE V4.1 six-bar families physically mirrored and locally "
            "adapted to candidate-3952 SMALL motion."
        ),
        "mirror_policy": (
            "physical mirror is canonical; equivalent scalar piston law to "
            "reverse crank traversal, while keeping visible crank rotation direction "
            "the same for S and L"
        ),
        "source": str(args.source),
        "target": str(args.target),
        "small_noise_band": noise_band,
        "ranks": ranks,
        "families": [],
    }

    wall0 = time.time()

    for family_index, rank in enumerate(ranks):
        large = by_rank[rank]["best"]
        mirror_x, primary_branch, second_branch = mirror_candidate(large)

        phase_x = fit_phase(
            mirror_x,
            args=args,
            primary_branch=primary_branch,
            second_branch=second_branch,
            theta=dense_theta,
            target_q=dense_q,
            target_dq=dense_dq,
            valid_mask=dense_valid,
        )

        phase_seed = evaluate_full(
            phase_x,
            primary_branch=primary_branch,
            second_branch=second_branch,
            theta=dense_theta,
            target_q=dense_q,
            target_dq=dense_dq,
            valid_mask=dense_valid,
            velocity_weight=args.velocity_weight,
        )
        phase_seed["feasible"] = total_violation(phase_seed, args) <= 0.0
        phase_seed["source_large_rank"] = rank

        bounds = local_bounds(phase_x, args)

        family = {
            "primary_rank": rank,
            "large": large,
            "mirrored_primary_branch": primary_branch,
            "mirrored_second_branch": second_branch,
            "phase_fitted_mirror_seed": phase_seed,
            "bounds": dict(zip(NAMES, bounds)),
            "runs": [],
            "best_small": phase_seed if phase_seed["feasible"] else None,
        }

        print(
            f"\nPAIR rank={rank}: "
            f"L={100*large['position_rms_full_cycle']:.5f}% "
            f"mirror+phase S={100*phase_seed['position_rms_full_cycle']:.5f}% "
            f"branches={primary_branch:+d}/{second_branch:+d}",
            flush=True,
        )

        for restart in range(args.restarts):
            run_seed = args.seed + 1000 * family_index + restart
            archive = []

            def objective(v):
                score = constrained_score(
                    v,
                    args=args,
                    primary_branch=primary_branch,
                    second_branch=second_branch,
                    theta=theta,
                    target_q=target_q,
                    target_dq=target_dq,
                    valid_mask=valid,
                )
                if score < 10.0:
                    archive.append((float(score), np.asarray(v, float).copy()))
                return score

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
                x0=phase_x,
            )

            vectors = [
                np.asarray(opt.x, float),
                phase_x.copy(),
            ]

            if getattr(opt, "population", None) is not None:
                order = np.argsort(np.asarray(opt.population_energies, float))
                vectors += [
                    np.asarray(opt.population[int(i)], float)
                    for i in order[:args.archive_per_run]
                ]

            vectors += [
                vector
                for _, vector in sorted(archive, key=lambda item: item[0])[
                    :args.archive_per_run
                ]
            ]

            dense_results = []
            seen = set()
            for vector in vectors:
                key = tuple(round(float(x), 10) for x in vector)
                if key in seen:
                    continue
                seen.add(key)

                try:
                    candidate = evaluate_full(
                        vector,
                        primary_branch=primary_branch,
                        second_branch=second_branch,
                        theta=dense_theta,
                        target_q=dense_q,
                        target_dq=dense_dq,
                        valid_mask=dense_valid,
                        velocity_weight=args.velocity_weight,
                    )
                except (
                    ValueError,
                    FloatingPointError,
                    OverflowError,
                    np.linalg.LinAlgError,
                ):
                    continue

                if total_violation(candidate, args) <= 0.0:
                    candidate["feasible"] = True
                    candidate["source_large_rank"] = rank
                    dense_results.append(candidate)

            best = (
                min(dense_results, key=lambda item: item["score"])
                if dense_results else None
            )

            family["runs"].append({
                "restart": restart,
                "seed": run_seed,
                "optimizer_success": bool(opt.success),
                "optimizer_message": str(opt.message),
                "function_evaluations": int(opt.nfev),
                "wall_seconds": time.time() - t0,
                "cpu_seconds": time.process_time() - p0,
                "best_dense": best,
            })

            candidates = [
                run["best_dense"]
                for run in family["runs"]
                if run.get("best_dense") is not None
            ]
            if phase_seed["feasible"]:
                candidates.append(phase_seed)

            family["best_small"] = (
                min(candidates, key=lambda item: item["score"])
                if candidates else None
            )

            if best:
                print(
                    f"  restart={restart}: "
                    f"S pos={100*best['position_rms_full_cycle']:.5f}% "
                    f"HP={100*best['position_rms_hp_band']:.3f}% "
                    f"nonHP={100*best['position_rms_non_hp']:.5f}% "
                    f"p-sine={best['minimum_primary_transmission_sine']:.3f} "
                    f"s-sine={best['minimum_secondary_transmission_sine']:.3f} "
                    f"rod/stroke={best['piston_rod_over_stroke']:.3f}",
                    flush=True,
                )
            else:
                print(f"  restart={restart}: no dense feasible SMALL candidate", flush=True)

            report["elapsed_wall_hours"] = (time.time() - wall0) / 3600.0
            temp = report["families"] + [family]
            old = report["families"]
            report["families"] = temp
            save(args.summary, report)
            report["families"] = old

        report["families"].append(family)

        small = family["best_small"]
        if small is not None:
            large_path = args.outdir / f"rank_{rank:02d}_large.json"
            small_path = args.outdir / f"rank_{rank:02d}_small.json"

            save(
                large_path,
                loader_file(
                    large,
                    f"Candidate 3952 LARGE polished family rank {rank}.",
                ),
            )
            save(
                small_path,
                loader_file(
                    small,
                    f"Candidate 3952 SMALL physical-mirror adaptation of LARGE family rank {rank}.",
                ),
            )

            family["large_mechanism_file"] = str(large_path)
            family["small_mechanism_file"] = str(small_path)

        report["elapsed_wall_hours"] = (time.time() - wall0) / 3600.0
        save(args.summary, report)

    print("\nMIRRORED PAIR SUMMARY")
    for family in report["families"]:
        large = family["large"]
        small = family["best_small"]
        if small is None:
            print(f"rank {family['primary_rank']:2d}: no feasible SMALL mirror adaptation")
            continue
        print(
            f"rank {family['primary_rank']:2d}: "
            f"L={100*large['position_rms_full_cycle']:.5f}% "
            f"S={100*small['position_rms_full_cycle']:.5f}% "
            f"combined-RMS="
            f"{100*math.sqrt(0.5*(large['position_rms_full_cycle']**2 + small['position_rms_full_cycle']**2)):.5f}%"
        )

    print(f"\nWrote {args.summary}")


if __name__ == "__main__":
    main()
