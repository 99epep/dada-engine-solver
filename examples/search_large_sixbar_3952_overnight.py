#!/usr/bin/env python3
"""Overnight global multi-basin search for the candidate-3952 LARGE six-bar.

This deliberately does NOT synthesize the small cylinder and does NOT use the
mirrored small mechanism to define the search box.

Strategy
--------
- large cylinder only;
- candidate-3952 position target;
- automatically detected six-front HP acceleration artefact excluded from fit;
- position-dominant score, tiny velocity tie-break;
- broad ABSOLUTE geometry bounds in crank-radius units;
- all four (primary_branch, second_branch) combinations explored;
- independent Latin-hypercube populations;
- two search modes:
    * free: current mechanical constraints;
    * comfortable: additionally primary transmission sine >= 0.40 and
      piston rod / stroke <= 2.5;
- current best is retained only as a reference/incumbent, never as the centre
  of the global bounds;
- JSON checkpoint after every island and a lightweight live checkpoint during
  each island.

Example:
    PYTHONPATH=src python3 examples/search_large_sixbar_3952_overnight.py \
        --hours 7
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import differential_evolution

import search_large_sixbar_stageL1 as stage
from search_motor_hybrid_3952_sixbar import (
    _detect_noise_bands,
    _theta_band_mask,
)


ROOT = Path.cwd()
TARGET = ROOT / "outputs/motor_hybrid_c2_15p_260k/candidate_3952_motion_target.csv"
CURRENT = ROOT / "outputs/large_sixbar_3952.json"
DEFAULT_OUTPUT = ROOT / "outputs/large_sixbar_3952_overnight.json"
DEFAULT_LIVE = ROOT / "outputs/large_sixbar_3952_overnight_live.json"

# Absolute bounds, all lengths in primary crank radii (AB = 1).
# These are intentionally much broader than any local seed neighbourhood.
BOUNDS = [
    (0.50, 12.0),               # primary_ground AD
    (0.50, 12.0),               # primary_coupler BC
    (0.35, 6.0),                # primary_rocker CD
    (-12.0, 18.0),              # primary_E_along
    (-15.0, 15.0),              # primary_E_normal
    (-math.pi, math.pi),        # primary_phase
    (-12.0, 20.0),              # second_pivot_x
    (-15.0, 15.0),              # second_pivot_y
    (0.35, 8.0),                # link_EF
    (0.50, 20.0),               # link_GF
    (-4.0, 5.0),                # H_along_over_EF
    (-4.0, 4.0),                # H_normal_over_EF
    (0.75, 12.0),               # piston_rod
    (-12.0, 12.0),              # slider_axis_offset
    (-math.pi, math.pi),        # slider_axis_angle
]

BRANCHES = [
    (+1, +1),
    (+1, -1),
    (-1, +1),
    (-1, -1),
]


def rms_score(r: dict, velocity_weight: float) -> float:
    return math.sqrt(
        float(r["position_rms"]) ** 2
        + velocity_weight * float(r["derivative_rms"]) ** 2
    )


def constraint_violation(
    r: dict,
    *,
    mode: str,
    stroke_min: float,
    stroke_max: float,
    primary_sine_min: float,
    secondary_sine_min: float,
    rod_cos_min: float,
    h_line_rms_max: float,
    h_axis_deg_max: float,
    clearance_min: float,
    comfortable_primary_sine_min: float,
    comfortable_rod_stroke_max: float,
) -> float:
    v = 0.0

    v += max(0.0, stroke_min - r["stroke_over_crank"]) / max(stroke_min, 1e-9)
    v += max(0.0, r["stroke_over_crank"] - stroke_max) / max(stroke_max, 1e-9)
    v += max(
        0.0, primary_sine_min - r["minimum_primary_transmission_sine"]
    ) / max(primary_sine_min, 1e-9)
    v += max(
        0.0, secondary_sine_min - r["minimum_secondary_transmission_sine"]
    ) / max(secondary_sine_min, 1e-9)
    v += max(
        0.0, rod_cos_min - r["minimum_rod_axis_cosine"]
    ) / max(1.0 - rod_cos_min, 1e-3)
    v += max(
        0.0, r["H_best_line_rms_over_stroke"] - h_line_rms_max
    ) / max(h_line_rms_max, 1e-9)
    v += max(
        0.0, r["H_best_line_axis_angle_deg"] - h_axis_deg_max
    ) / max(h_axis_deg_max, 1e-9)
    v += max(
        0.0, clearance_min - r["crank_axis_to_EFH_clearance_over_crank"]
    ) / max(clearance_min, 0.1)

    if mode == "comfortable":
        v += max(
            0.0,
            comfortable_primary_sine_min
            - r["minimum_primary_transmission_sine"],
        ) / max(comfortable_primary_sine_min, 1e-9)
        v += max(
            0.0,
            r["piston_rod_over_stroke"] - comfortable_rod_stroke_max,
        ) / max(comfortable_rod_stroke_max, 1e-9)

    return float(v)


def evaluate_candidate(
    vector: np.ndarray,
    *,
    primary_branch: int,
    second_branch: int,
    theta: np.ndarray,
    q: np.ndarray,
    dq: np.ndarray,
    fit_mask: np.ndarray,
    maximum_eh: float,
    velocity_weight: float,
    mode: str,
    settings: dict,
) -> tuple[float, dict | None]:
    try:
        r = stage.evaluate(
            np.asarray(vector, dtype=float),
            primary_branch,
            second_branch,
            theta,
            q,
            dq,
            maximum_eh,
            score_mask=fit_mask,
        )
    except (ValueError, FloatingPointError, OverflowError):
        return 1e6, None

    violation = constraint_violation(r, mode=mode, **settings)
    if violation > 0.0:
        return 10.0 + 100.0 * violation, r

    return rms_score(r, velocity_weight), r


def dense_result(
    vector: np.ndarray,
    *,
    primary_branch: int,
    second_branch: int,
    dense_theta: np.ndarray,
    dense_q: np.ndarray,
    dense_dq: np.ndarray,
    dense_fit_mask: np.ndarray,
    maximum_eh: float,
    velocity_weight: float,
    mode: str,
    settings: dict,
) -> dict | None:
    try:
        raw = stage.evaluate(
            vector,
            primary_branch,
            second_branch,
            dense_theta,
            dense_q,
            dense_dq,
            maximum_eh,
        )
        fit = stage.evaluate(
            vector,
            primary_branch,
            second_branch,
            dense_theta,
            dense_q,
            dense_dq,
            maximum_eh,
            score_mask=dense_fit_mask,
        )
    except (ValueError, FloatingPointError, OverflowError):
        return None

    if constraint_violation(fit, mode=mode, **settings) > 0.0:
        return None

    out = dict(raw)
    out["position_rms_raw"] = float(raw["position_rms"])
    out["derivative_rms_raw"] = float(raw["derivative_rms"])
    out["position_rms"] = float(fit["position_rms"])
    out["derivative_rms"] = float(fit["derivative_rms"])
    out["fit_position_rms"] = float(fit["position_rms"])
    out["fit_derivative_rms"] = float(fit["derivative_rms"])
    out["score"] = rms_score(fit, velocity_weight)
    out["motion_score"] = out["score"]
    out["primary_branch"] = int(primary_branch)
    out["second_branch"] = int(second_branch)
    out["search_mode"] = mode
    out["feasible"] = True
    return out


def vector_from_candidate(candidate: dict) -> np.ndarray:
    p = candidate["parameters"]
    return np.array([float(p[name]) for name in stage.NAMES], dtype=float)


def inside_bounds(x: np.ndarray) -> bool:
    return all(lo <= float(v) <= hi for v, (lo, hi) in zip(x, BOUNDS))


def best_or_none(items: list[dict]) -> dict | None:
    items = [x for x in items if x is not None]
    return min(items, key=lambda x: x["score"]) if items else None


def json_dump(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", type=Path, default=TARGET)
    ap.add_argument("--current", type=Path, default=CURRENT)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--live-output", type=Path, default=DEFAULT_LIVE)

    ap.add_argument("--hours", type=float, default=7.0)
    ap.add_argument("--seed", type=int, default=3952000)

    # Independent-island size. One island is deliberately bounded so all four
    # branch combinations get repeated chances during a long run.
    ap.add_argument("--generations-per-island", type=int, default=220)
    ap.add_argument("--population-size", type=int, default=14)
    ap.add_argument("--max-islands", type=int, default=1000)
    ap.add_argument("--dense-stride", type=int, default=1)
    ap.add_argument("--coarse-stride", type=int, default=12)
    ap.add_argument("--progress-generations", type=int, default=20)

    ap.add_argument("--velocity-weight", type=float, default=0.001)
    ap.add_argument("--maximum-EH", type=float, default=7.0)

    ap.add_argument("--stroke-min", type=float, default=1.0)
    ap.add_argument("--stroke-max", type=float, default=3.0)
    ap.add_argument("--primary-sine-min", type=float, default=0.30)
    ap.add_argument("--secondary-sine-min", type=float, default=0.30)
    ap.add_argument("--rod-cos-min", type=float, default=0.95)
    ap.add_argument("--h-line-rms-max", type=float, default=0.25)
    ap.add_argument("--h-axis-deg-max", type=float, default=10.0)
    ap.add_argument("--clearance-min", type=float, default=0.50)

    ap.add_argument("--comfortable-primary-sine-min", type=float, default=0.40)
    ap.add_argument("--comfortable-rod-stroke-max", type=float, default=2.50)
    ap.add_argument("--noise-fronts", type=int, default=6)

    args = ap.parse_args()
    if args.hours <= 0:
        ap.error("--hours must be > 0")
    if args.population_size < 4:
        ap.error("--population-size must be >= 4")

    raw = np.genfromtxt(args.target, delimiter=",", names=True)
    if raw["theta_rad"][-1] >= 2.0 * math.pi - 1e-9:
        raw = raw[:-1]

    noise = _detect_noise_bands(args.target, args.noise_fronts)["large"]

    dense = raw[::args.dense_stride]
    dense_theta = np.asarray(dense["theta_rad"], dtype=float)
    dense_q = np.asarray(dense["large_fraction_0_1"], dtype=float)
    dense_dq = np.asarray(dense["large_dq_dtheta_per_rad"], dtype=float)
    dense_fit_mask = ~_theta_band_mask(dense_theta, noise)

    coarse = raw[::args.coarse_stride]
    theta = np.asarray(coarse["theta_rad"], dtype=float)
    q = np.asarray(coarse["large_fraction_0_1"], dtype=float)
    dq = np.asarray(coarse["large_dq_dtheta_per_rad"], dtype=float)
    fit_mask = ~_theta_band_mask(theta, noise)

    settings = {
        "stroke_min": args.stroke_min,
        "stroke_max": args.stroke_max,
        "primary_sine_min": args.primary_sine_min,
        "secondary_sine_min": args.secondary_sine_min,
        "rod_cos_min": args.rod_cos_min,
        "h_line_rms_max": args.h_line_rms_max,
        "h_axis_deg_max": args.h_axis_deg_max,
        "clearance_min": args.clearance_min,
        "comfortable_primary_sine_min": args.comfortable_primary_sine_min,
        "comfortable_rod_stroke_max": args.comfortable_rod_stroke_max,
    }

    started = time.monotonic()
    deadline = started + args.hours * 3600.0

    report = {
        "description": "Candidate 3952 large-cylinder overnight global multi-basin search.",
        "target": str(args.target),
        "noise_band": noise,
        "length_unit": "primary_crank_radius",
        "bounds": dict(zip(stage.NAMES, BOUNDS)),
        "settings": {
            **vars(args),
            "target": str(args.target),
            "current": str(args.current),
            "output": str(args.output),
            "live_output": str(args.live_output),
        },
        "branch_combinations": [list(x) for x in BRANCHES],
        "modes": {
            "free": "baseline mechanical constraints",
            "comfortable": (
                f"primary sine >= {args.comfortable_primary_sine_min} and "
                f"rod/stroke <= {args.comfortable_rod_stroke_max}"
            ),
        },
        "reference_current_best": None,
        "islands": [],
        "best": None,
        "best_by_mode": {},
        "best_by_branch": {},
    }

    # Current candidate is a reference only. It does not define BOUNDS.
    current_vector = None
    current_branch = None
    if args.current.exists():
        current_data = json.loads(args.current.read_text(encoding="utf-8"))
        current = current_data.get("best")
        if current is not None:
            current_vector = vector_from_candidate(current)
            current_branch = (
                int(current["primary"].get("assembly_branch", 1)),
                int(current["second_branch"]),
            )
            if inside_bounds(current_vector):
                ref = dense_result(
                    current_vector,
                    primary_branch=current_branch[0],
                    second_branch=current_branch[1],
                    dense_theta=dense_theta,
                    dense_q=dense_q,
                    dense_dq=dense_dq,
                    dense_fit_mask=dense_fit_mask,
                    maximum_eh=args.maximum_EH,
                    velocity_weight=args.velocity_weight,
                    mode="free",
                    settings=settings,
                )
                report["reference_current_best"] = ref

    schedule = []
    # First eight islands guarantee one pass through every branch in both modes.
    for mode in ("free", "comfortable"):
        for branches in BRANCHES:
            schedule.append((mode, branches))

    island_index = 0
    while island_index < args.max_islands and time.monotonic() < deadline:
        mode, (primary_branch, second_branch) = schedule[
            island_index % len(schedule)
        ]
        island_seed = args.seed + island_index

        generation = 0
        best_live_score = math.inf
        best_live_vector = None

        def objective(v):
            score, _ = evaluate_candidate(
                v,
                primary_branch=primary_branch,
                second_branch=second_branch,
                theta=theta,
                q=q,
                dq=dq,
                fit_mask=fit_mask,
                maximum_eh=args.maximum_EH,
                velocity_weight=args.velocity_weight,
                mode=mode,
                settings=settings,
            )
            return score

        def callback(xk, convergence):
            nonlocal generation, best_live_score, best_live_vector
            generation += 1
            score = float(objective(xk))
            if score < best_live_score:
                best_live_score = score
                best_live_vector = np.asarray(xk, dtype=float).copy()

            if (
                generation == 1
                or generation % args.progress_generations == 0
                or time.monotonic() >= deadline
            ):
                live = {
                    "island": island_index,
                    "mode": mode,
                    "primary_branch": primary_branch,
                    "second_branch": second_branch,
                    "generation": generation,
                    "coarse_score": best_live_score,
                    "convergence": float(convergence),
                    "elapsed_hours": (time.monotonic() - started) / 3600.0,
                    "vector": (
                        best_live_vector.tolist()
                        if best_live_vector is not None
                        else None
                    ),
                }
                json_dump(args.live_output, live)
                print(
                    f"island {island_index:03d} {mode:11s} "
                    f"branches={primary_branch:+d}/{second_branch:+d} "
                    f"gen={generation:4d} coarse={best_live_score:.6f} "
                    f"conv={convergence:.3g}",
                    flush=True,
                )

            return time.monotonic() >= deadline

        # Only one reference island is allowed to contain the current incumbent.
        # Every other island is a fresh Latin-hypercube population over BOUNDS.
        x0 = None
        if (
            island_index == 0
            and current_vector is not None
            and current_branch == (primary_branch, second_branch)
            and inside_bounds(current_vector)
        ):
            x0 = current_vector

        t0 = time.monotonic()
        opt = differential_evolution(
            objective,
            BOUNDS,
            seed=island_seed,
            popsize=args.population_size,
            maxiter=args.generations_per_island,
            polish=True,
            tol=1e-8,
            updating="immediate",
            workers=1,
            init="latinhypercube",
            x0=x0,
            callback=callback,
        )
        elapsed = time.monotonic() - t0

        candidate = dense_result(
            np.asarray(opt.x, dtype=float),
            primary_branch=primary_branch,
            second_branch=second_branch,
            dense_theta=dense_theta,
            dense_q=dense_q,
            dense_dq=dense_dq,
            dense_fit_mask=dense_fit_mask,
            maximum_eh=args.maximum_EH,
            velocity_weight=args.velocity_weight,
            mode=mode,
            settings=settings,
        )

        island = {
            "index": island_index,
            "seed": island_seed,
            "mode": mode,
            "primary_branch": primary_branch,
            "second_branch": second_branch,
            "generations_completed": generation,
            "optimizer_success": bool(opt.success),
            "optimizer_message": str(opt.message),
            "function_evaluations": int(opt.nfev),
            "elapsed_seconds": elapsed,
            "candidate": candidate,
        }
        report["islands"].append(island)

        all_candidates = [
            item["candidate"]
            for item in report["islands"]
            if item.get("candidate") is not None
        ]
        if report["reference_current_best"] is not None:
            all_candidates.append(report["reference_current_best"])

        report["best"] = best_or_none(all_candidates)

        for m in ("free", "comfortable"):
            c = [
                item["candidate"]
                for item in report["islands"]
                if item.get("candidate") is not None
                and item["candidate"]["search_mode"] == m
            ]
            report["best_by_mode"][m] = best_or_none(c)

        for pb, sb in BRANCHES:
            key = f"{pb:+d}/{sb:+d}"
            c = [
                item["candidate"]
                for item in report["islands"]
                if item.get("candidate") is not None
                and int(item["candidate"]["primary_branch"]) == pb
                and int(item["candidate"]["second_branch"]) == sb
            ]
            report["best_by_branch"][key] = best_or_none(c)

        report["elapsed_hours"] = (time.monotonic() - started) / 3600.0
        json_dump(args.output, report)

        if candidate is None:
            print(
                f"completed island {island_index:03d}: no dense feasible candidate",
                flush=True,
            )
        else:
            print(
                f"completed island {island_index:03d}: "
                f"fit={100*candidate['position_rms']:.4f}% "
                f"raw={100*candidate['position_rms_raw']:.4f}% "
                f"stroke/r={candidate['stroke_over_crank']:.3f} "
                f"rod/stroke={candidate['piston_rod_over_stroke']:.3f} "
                f"primary-sine={candidate['minimum_primary_transmission_sine']:.3f}",
                flush=True,
            )

        island_index += 1

    report["elapsed_hours"] = (time.monotonic() - started) / 3600.0
    report["stopped_by_time_budget"] = time.monotonic() >= deadline
    json_dump(args.output, report)

    print("\nOVERNIGHT SEARCH SUMMARY")
    b = report.get("best")
    if b is None:
        print("No dense feasible candidate found.")
    else:
        print(
            f"best fit position RMS = {100*b['position_rms']:.6f}%\n"
            f"full-cycle position RMS = {100*b['position_rms_raw']:.6f}%\n"
            f"branches = {b['primary_branch']:+d}/{b['second_branch']:+d}\n"
            f"mode = {b['search_mode']}\n"
            f"stroke/r = {b['stroke_over_crank']:.6f}\n"
            f"rod/stroke = {b['piston_rod_over_stroke']:.6f}\n"
            f"primary sine = {b['minimum_primary_transmission_sine']:.6f}"
        )
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
