#!/usr/bin/env python3
"""Fresh-island saturation study for candidate-3952 LARGE V4.1 primary search.

This experiment estimates algorithmic basin-capture saturation, not literal
geometric coverage of a continuous six-dimensional parameter space.

Every island is fresh: no historical x0 seed, Latin-hypercube initialization
across the full current bounds, deterministic independent seeds, and strict
alternation of the two primary assembly branches.

Default overnight run:
  PYTHONPATH=src python3 examples/search_large_primary_cadence_3952_v4_1_saturation.py

Default: 512 islands, cumulative checkpoints every 64 islands. The JSON is
saved after every island and the run resumes automatically if interrupted.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import differential_evolution

import search_large_primary_cadence_3952_v4_1 as v41

ROOT = Path.cwd()
OUTPUT = ROOT / "outputs" / "large_primary_cadence_3952_v4_1_saturation.json"
NAMES = v41.v4.NAMES
BOUNDS = v41.v4.BOUNDS


def geometry_distance(a: dict, b: dict) -> float:
    return v41.v4.geometry_distance(a, b)


def connected_components(records: list[dict], threshold: float) -> list[list[int]]:
    """Order-independent operational families under a distance threshold."""
    n = len(records)
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if geometry_distance(records[i], records[j]) < threshold:
                adj[i].append(j)
                adj[j].append(i)

    seen = [False] * n
    comps = []
    for start in range(n):
        if seen[start]:
            continue
        stack = [start]
        seen[start] = True
        comp = []
        while stack:
            i = stack.pop()
            comp.append(i)
            for j in adj[i]:
                if not seen[j]:
                    seen[j] = True
                    stack.append(j)
        comps.append(comp)
    return comps


def chao1(s_obs: int, f1: int, f2: int) -> float:
    if f2 > 0:
        return float(s_obs + f1 * f1 / (2.0 * f2))
    return float(s_obs + f1 * max(f1 - 1, 0) / 2.0)


def pareto_count(records: list[dict]) -> int:
    """Non-dominated count in score(min), sine(max), BE/BC(min)."""
    count = 0
    for i, a in enumerate(records):
        dominated = False
        for j, b in enumerate(records):
            if i == j:
                continue
            weak = (
                b["score"] <= a["score"]
                and b["minimum_primary_transmission_sine"]
                >= a["minimum_primary_transmission_sine"]
                and b["BE_over_BC"] <= a["BE_over_BC"]
            )
            strict = (
                b["score"] < a["score"]
                or b["minimum_primary_transmission_sine"]
                > a["minimum_primary_transmission_sine"]
                or b["BE_over_BC"] < a["BE_over_BC"]
            )
            if weak and strict:
                dominated = True
                break
        if not dominated:
            count += 1
    return count


def family_summary(records: list[dict], threshold: float) -> dict:
    comps = connected_components(records, threshold)
    counts = sorted((len(c) for c in comps), reverse=True)
    n = len(records)
    f1 = sum(c == 1 for c in counts)
    f2 = sum(c == 2 for c in counts)

    families = []
    for comp in comps:
        members = [records[i] for i in comp]
        rep = min(members, key=lambda x: x["score"])
        families.append({
            "count": len(comp),
            "capture_share": len(comp) / max(n, 1),
            "best_score": rep["score"],
            "best_branch": rep["branch"],
            "best_sine": rep["minimum_primary_transmission_sine"],
            "best_BE_over_BC": rep["BE_over_BC"],
            "best_parameters": rep["parameters"],
        })
    families.sort(key=lambda x: (-x["count"], x["best_score"]))

    return {
        "distance_threshold": float(threshold),
        "observed_families": len(comps),
        "singletons_f1": f1,
        "doubletons_f2": f2,
        "good_turing_unseen_capture_mass_estimate": f1 / n if n else None,
        "chao1_family_count_estimate": chao1(len(comps), f1, f2) if n else None,
        "top1_capture_share": counts[0] / n if counts else None,
        "top4_capture_share": sum(counts[:4]) / n if counts else None,
        "top10_capture_share": sum(counts[:10]) / n if counts else None,
        "families_by_capture": families[:25],
    }


def summarize(records: list[dict], thresholds: list[float]) -> dict:
    best = min(records, key=lambda x: x["score"]) if records else None
    return {
        "completed_winners": len(records),
        "best_score": best["score"] if best else None,
        "best_branch": best["branch"] if best else None,
        "best_sine": best["minimum_primary_transmission_sine"] if best else None,
        "best_BE_over_BC": best["BE_over_BC"] if best else None,
        "pareto_count_score_sine_BEBC": pareto_count(records) if records else 0,
        "families": [family_summary(records, t) for t in thresholds],
    }


def dense_best(opt, branch, target, dense_theta, args, top_population):
    vectors = [np.asarray(opt.x, dtype=float)]
    if getattr(opt, "population", None) is not None:
        order = np.argsort(np.asarray(opt.population_energies, dtype=float))
        vectors += [
            np.asarray(opt.population[int(i)], dtype=float)
            for i in order[:top_population]
        ]

    records, seen = [], set()
    for vector in vectors:
        key = tuple(round(float(x), 10) for x in vector)
        if key in seen:
            continue
        seen.add(key)
        record = v41.v4.make_record(vector, branch, target, dense_theta, args)
        if record is not None:
            records.append(record)
    return min(records, key=lambda x: x["score"]) if records else None


def save(path: Path, obj: dict):
    v41.v4.save(path, obj)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--total-islands", type=int, default=512)
    ap.add_argument("--checkpoint-every", type=int, default=64)
    ap.add_argument("--generations", type=int, default=260)
    ap.add_argument("--population-size", type=int, default=16)
    ap.add_argument("--seed", type=int, default=3954100)
    ap.add_argument("--output", type=Path, default=OUTPUT)
    ap.add_argument("--family-thresholds", nargs="+", type=float, default=[0.03, 0.04, 0.05])
    ap.add_argument("--dense-top-population", type=int, default=12)
    ap.add_argument("--restart", action="store_true")
    cli = ap.parse_args()

    if cli.total_islands < 2:
        ap.error("--total-islands must be >= 2")
    if cli.checkpoint_every < 1:
        ap.error("--checkpoint-every must be >= 1")

    args = v41.default_args()
    args.generations = cli.generations
    args.population_size = cli.population_size
    args.historical_seed_count = 0
    v41._install_v41(args)
    target = v41.v4.load_target(args.target, args)

    nc = int(round(360.0 / args.coarse_step_deg))
    nd = int(round(360.0 / args.dense_step_deg))
    coarse_theta = np.linspace(0.0, 2.0 * math.pi, nc, endpoint=False)
    dense_theta = np.linspace(0.0, 2.0 * math.pi, nd, endpoint=False)

    islands = []
    if cli.output.exists() and not cli.restart:
        old = json.loads(cli.output.read_text(encoding="utf-8"))
        settings = old["experiment_settings"]
        required = {
            "base_seed": cli.seed,
            "generations": cli.generations,
            "population_size": cli.population_size,
            "family_thresholds": cli.family_thresholds,
        }
        for key, value in required.items():
            if settings.get(key) != value:
                raise SystemExit(
                    f"Cannot resume: {key} changed from {settings.get(key)!r} to {value!r}. "
                    "Use --restart or another output file."
                )
        islands = old.get("islands", [])
        print(f"Resuming with {len(islands)} completed islands", flush=True)

    report = {
        "description": (
            "Fresh-island V4.1 saturation experiment. Family frequencies estimate "
            "algorithmic basin-capture probability, not literal geometric volume."
        ),
        "experiment_settings": {
            "total_islands_requested": cli.total_islands,
            "checkpoint_every": cli.checkpoint_every,
            "base_seed": cli.seed,
            "generations": cli.generations,
            "population_size": cli.population_size,
            "initial_population_per_island": len(BOUNDS) * cli.population_size,
            "family_thresholds": cli.family_thresholds,
            "dense_top_population": cli.dense_top_population,
            "branch_policy": "strict alternation +1,-1",
            "historical_seeds": 0,
        },
        "search_bounds": dict(zip(NAMES, BOUNDS)),
        "target_structure": {
            "turnarounds_deg": target["turnarounds_deg"],
            "short_span_deg": target["short_span_deg"],
            "long_span_deg": target["long_span_deg"],
            "short_split": target["short_split"],
            "fast_slow_speed_ratio": target["target_fast_slow_speed_ratio"],
            "fast_displacement_fraction": target["target_fast_displacement_fraction"],
        },
        "islands": islands,
        "milestones": [],
    }

    session_start = time.time()
    for island_index in range(len(islands), cli.total_islands):
        branch = +1 if island_index % 2 == 0 else -1
        run_seed = cli.seed + island_index

        def objective(vector):
            vector = np.asarray(vector, dtype=float)
            closure = v41.v4.primary_closure_penalty(vector)
            if closure > 0.0:
                return 20.0 + 100.0 * closure
            try:
                features = v41.primary_features_v41(vector, branch, coarse_theta, args)
                score, _ = v41.v4.cadence_score(features, target, args)
                return float(score)
            except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
                return 1e3

        t0, p0 = time.time(), time.process_time()
        opt = differential_evolution(
            objective,
            BOUNDS,
            seed=run_seed,
            popsize=cli.population_size,
            maxiter=cli.generations,
            polish=True,
            tol=1e-8,
            updating="immediate",
            workers=1,
            init="latinhypercube",
            x0=None,
        )
        best = dense_best(opt, branch, target, dense_theta, args, cli.dense_top_population)
        islands.append({
            "index": island_index,
            "branch": branch,
            "optimizer_seed": run_seed,
            "optimizer_success": bool(opt.success),
            "optimizer_message": str(opt.message),
            "function_evaluations": int(opt.nfev),
            "wall_seconds": time.time() - t0,
            "cpu_seconds": time.process_time() - p0,
            "best": best,
        })

        winners = [x["best"] for x in islands if x.get("best") is not None]
        report["completed_islands"] = len(islands)
        report["successful_dense_winners"] = len(winners)
        report["cumulative_island_wall_hours"] = sum(x["wall_seconds"] for x in islands) / 3600.0
        report["current_session_wall_hours"] = (time.time() - session_start) / 3600.0
        report["best"] = min(winners, key=lambda x: x["score"]) if winners else None

        is_checkpoint = (
            len(islands) % cli.checkpoint_every == 0
            or len(islands) == cli.total_islands
        )
        if is_checkpoint:
            current = summarize(winners, cli.family_thresholds) if winners else None
            report["current_family_analysis"] = current
            # Replace or append this exact cumulative checkpoint only.
            report["milestones"] = [
                m for m in report.get("milestones", [])
                if m.get("completed_winners") != len(winners)
            ]
            if current is not None:
                report["milestones"].append(current)
                report["milestones"].sort(key=lambda x: x["completed_winners"])
        save(cli.output, report)

        if best is None:
            print(f"island {island_index:04d} branch={branch:+d}: no dense winner", flush=True)
        else:
            m = best["cadence_match"]
            print(
                f"island {island_index:04d} branch={branch:+d}: "
                f"score={best['score']:.5f} turn={m['turning_rms_deg']:.2f}deg "
                f"ratio={m['fast_slow_speed_ratio']:.2f} "
                f"symBP={m['long_mirror_asymmetry_rms']:.3f} "
                f"sine={best['minimum_primary_transmission_sine']:.3f} "
                f"BE/BC={best['BE_over_BC']:.2f}",
                flush=True,
            )

        if is_checkpoint:
            s = report["current_family_analysis"]
            fam = min(s["families"], key=lambda x: abs(x["distance_threshold"] - 0.04))
            print(
                f"\nSATURATION {len(islands)}/{cli.total_islands}: "
                f"best={s['best_score']:.6f}, families@{fam['distance_threshold']:.2f}="
                f"{fam['observed_families']}, f1={fam['singletons_f1']}, f2={fam['doubletons_f2']}, "
                f"unseen_mass~{fam['good_turing_unseen_capture_mass_estimate']:.3f}, "
                f"top4_share={fam['top4_capture_share']:.3f}\n",
                flush=True,
            )

    print(f"Completed {len(islands)} fresh islands. Wrote {cli.output}", flush=True)


if __name__ == "__main__":
    main()
