#!/usr/bin/env python3
"""Search LARGE primary four-bars whose coupler point E has the target cadence.

The primary is not asked to reproduce piston position.  It is scored mainly on
the timing of three dominant tangential-acceleration blocks of E.

Target blocks are extracted automatically from candidate 3952.  The artificial
six-front HP acceleration-noise band is bridged in target VELOCITY before the
reference acceleration lobes are detected.

Typical:
  PYTHONPATH=src python3 examples/search_large_primary_cadence_3952.py

Long:
  PYTHONPATH=src python3 examples/search_large_primary_cadence_3952.py \
      --hours 6 --islands 10000
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import differential_evolution, linear_sum_assignment

import search_large_sixbar_stageL1 as stage
import search_motor_hybrid_3952_sixbar as motor3952


ROOT = Path.cwd()
TARGET = ROOT / "outputs/motor_hybrid_c2_15p_260k/candidate_3952_motion_target.csv"
OUTPUT = ROOT / "outputs/large_primary_cadence_3952.json"
TRACE = ROOT / "outputs/large_primary_cadence_3952_best_trace.csv"

NAMES = [
    "primary_ground", "primary_coupler", "primary_rocker",
    "primary_E_along", "primary_E_normal", "primary_phase",
]
BOUNDS = [
    (0.50, 12.0), (0.50, 12.0), (0.35, 6.0),
    (-12.0, 18.0), (-15.0, 15.0), (-math.pi, math.pi),
]


def cyclic_distance_deg(a, b):
    return abs((float(a) - float(b) + 180.0) % 360.0 - 180.0)


def smooth_periodic(y, step_deg, width_deg):
    y = np.asarray(y, float)
    if width_deg <= 0:
        return y.copy()
    half = max(1, int(round(0.5 * width_deg / step_deg)))
    out = np.zeros_like(y)
    for k in range(-half, half + 1):
        out += np.roll(y, k)
    return out / (2 * half + 1)


def derivative_periodic(y, step_rad):
    y = np.asarray(y, float)
    return (np.roll(y, -1) - np.roll(y, 1)) / (2.0 * step_rad)


def bridge_mask(y, mask):
    """Linearly bridge one contiguous cyclic True interval."""
    y = np.asarray(y, float)
    mask = np.asarray(mask, bool)
    n = len(y)
    if not np.any(mask):
        return y.copy()
    starts = [i for i in range(n) if mask[i] and not mask[(i - 1) % n]]
    ends = [i for i in range(n) if mask[i] and not mask[(i + 1) % n]]
    if len(starts) != 1 or len(ends) != 1:
        raise ValueError("expected one contiguous cyclic band")
    left = (starts[0] - 1) % n
    right = (ends[0] + 1) % n
    total = (right - left) % n
    out = y.copy()
    for k in range(1, total):
        i = (left + k) % n
        f = k / total
        out[i] = (1 - f) * y[left] + f * y[right]
    return out


def extract_lobes(theta_deg, signal, epsilon_fraction):
    """Circular signed acceleration lobes, strongest first."""
    theta_deg = np.asarray(theta_deg, float)
    signal = np.asarray(signal, float)
    n = len(signal)
    step_deg = 360.0 / n
    step_rad = 2.0 * math.pi / n
    peak = float(np.max(np.abs(signal)))
    if peak <= 1e-14:
        return []

    eps = epsilon_fraction * peak
    signs = np.zeros(n, np.int8)
    signs[signal > eps] = 1
    signs[signal < -eps] = -1

    # Close only tiny zero holes.
    for _ in range(4):
        old = signs.copy()
        for i in range(n):
            if old[i] == 0:
                l, r = old[(i - 1) % n], old[(i + 1) % n]
                if l != 0 and l == r:
                    signs[i] = l

    segs = []
    start, current = 0, int(signs[0])
    for i in range(1, n):
        s = int(signs[i])
        if s != current:
            segs.append([start, i - 1, current])
            start, current = i, s
    segs.append([start, n - 1, current])

    if len(segs) > 1 and segs[0][2] != 0 and segs[0][2] == segs[-1][2]:
        first, last = segs.pop(0), segs.pop(-1)
        segs.insert(0, [last[0], first[1] + n, first[2]])

    out = []
    for start, end, sign in segs:
        if sign == 0:
            continue
        idx = np.array([i % n for i in range(start, end + 1)], int)
        values = signal[idx]
        delta = float(np.sum(values) * step_rad)
        strength = abs(delta)
        if strength <= 1e-12:
            continue

        w = np.abs(values)
        ang = np.radians(theta_deg[idx])
        center = math.degrees(math.atan2(
            float(np.sum(w * np.sin(ang))),
            float(np.sum(w * np.cos(ang))),
        )) % 360.0

        out.append({
            "sign": int(sign),
            "start_deg": float(theta_deg[start % n] % 360.0),
            "end_deg": float(theta_deg[end % n] % 360.0),
            "center_deg": center,
            "span_deg": float(len(idx) * step_deg),
            "delta_speed": delta,
            "strength": strength,
            "peak_abs_acceleration": float(np.max(w)),
        })

    out.sort(key=lambda x: x["strength"], reverse=True)
    return out


def load_target(path, smoothing_deg, epsilon_fraction, noise_fronts):
    raw = np.genfromtxt(path, delimiter=",", names=True)
    if raw["theta_rad"][-1] >= 2 * math.pi - 1e-9:
        raw = raw[:-1]
    theta = np.asarray(raw["theta_rad"], float)
    theta_deg = np.mod(np.degrees(theta), 360.0)
    dq = np.asarray(raw["large_dq_dtheta_per_rad"], float)

    band = motor3952._detect_noise_bands(path, noise_fronts)["large"]
    noise = motor3952._theta_band_mask(theta, band)
    dq_bridge = bridge_mask(dq, noise)

    step_deg = 360.0 / len(theta)
    step_rad = 2 * math.pi / len(theta)
    dq_smooth = smooth_periodic(dq_bridge, step_deg, smoothing_deg)
    acc = derivative_periodic(dq_smooth, step_rad)
    lobes = extract_lobes(theta_deg, acc, epsilon_fraction)
    if len(lobes) < 4:
        raise RuntimeError("target lobe extraction failed")

    blocks = sorted(lobes[:3], key=lambda x: x["center_deg"])
    return {
        "theta": theta, "theta_deg": theta_deg,
        "acceleration": acc,
        "noise_band": band,
        "all_lobes": lobes,
        "blocks": blocks,
    }


def primary_cadence(v, branch, theta, smoothing_deg, epsilon_fraction):
    E, Ed, sine, pdata = stage.primary(theta, v, branch)
    n = len(theta)
    step_deg = 360.0 / n
    step_rad = 2 * math.pi / n

    speed = np.linalg.norm(Ed, axis=1)
    if np.min(speed) <= 1e-8:
        raise ValueError("E speed nearly vanishes")
    speed_smooth = smooth_periodic(speed, step_deg, smoothing_deg)
    at = derivative_periodic(speed_smooth, step_rad)

    Edd = np.column_stack([
        derivative_periodic(Ed[:, 0], step_rad),
        derivative_periodic(Ed[:, 1], step_rad),
    ])
    an = (Ed[:, 0] * Edd[:, 1] - Ed[:, 1] * Edd[:, 0]) / speed

    direction = np.unwrap(np.arctan2(Ed[:, 1], Ed[:, 0]))
    direction_tv = float(np.sum(np.abs(np.diff(direction, append=direction[0] + 2*math.pi*round((direction[-1]-direction[0])/(2*math.pi))))))

    lobes = extract_lobes(np.degrees(theta) % 360.0, at, epsilon_fraction)
    span = max(float(np.ptp(E[:, 0])), float(np.ptp(E[:, 1])))
    path_length = float(np.sum(np.linalg.norm(np.roll(E, -1, axis=0) - E, axis=1)))

    return {
        "E": E, "speed": speed, "speed_smooth": speed_smooth,
        "tangential_acceleration": at, "normal_acceleration": an,
        "lobes": lobes,
        "minimum_primary_transmission_sine": float(sine),
        "primary": pdata,
        "E_span_over_crank": span,
        "E_path_length_over_crank": path_length,
        "direction_total_variation_rad": direction_tv,
        "normal_acceleration_rms": float(np.sqrt(np.mean(an * an))),
        "tangential_acceleration_rms": float(np.sqrt(np.mean(at * at))),
    }


def score_cadence(c, target_blocks, args):
    sine = c["minimum_primary_transmission_sine"]
    if sine < args.primary_sine_floor:
        return 20 + 100 * (args.primary_sine_floor - sine) / args.primary_sine_floor, None
    if c["E_span_over_crank"] < args.E_span_floor:
        return 20 + 100 * (args.E_span_floor - c["E_span_over_crank"]) / args.E_span_floor, None

    lobes = c["lobes"]
    if len(lobes) < 3:
        return 50.0, None

    strongest = lobes[0]["strength"]
    pool = [x for x in lobes if x["strength"] >= 0.02 * strongest][:10]
    if len(pool) < 3:
        pool = lobes[:3]

    t_strength = np.array([x["strength"] for x in target_blocks], float)
    t_strength /= np.sum(t_strength)
    best = None

    # Relative speed-up/slow-down structure matters weakly; allow global sign flip.
    for flip in (1, -1):
        cost = np.empty((3, len(pool)))
        for i, t in enumerate(target_blocks):
            for j, p in enumerate(pool):
                center = cyclic_distance_deg(p["center_deg"], t["center_deg"]) / args.center_scale_deg
                width = math.log(max(p["span_deg"], 1.0) / max(t["span_deg"], 1.0))
                sign = 0.0 if p["sign"] == flip * t["sign"] else 1.0
                cost[i, j] = center**2 + args.width_weight*width**2 + args.sign_weight*sign

        rows, cols = linear_sum_assignment(cost)
        assigned = [pool[int(j)] for j in cols]
        center_errors = np.array([
            cyclic_distance_deg(assigned[i]["center_deg"], target_blocks[i]["center_deg"])
            for i in range(3)
        ])
        widths = np.array([
            math.log(max(assigned[i]["span_deg"], 1.0) / max(target_blocks[i]["span_deg"], 1.0))
            for i in range(3)
        ])
        c_strength = np.array([x["strength"] for x in assigned], float)
        c_strength /= np.sum(c_strength)
        strength_error = float(np.sqrt(np.mean((c_strength - t_strength)**2)))

        matched = float(sum(x["strength"] for x in assigned))
        total = float(sum(x["strength"] for x in lobes))
        extra = max(0.0, (total - matched) / max(matched, 1e-12))

        base = math.sqrt(float(np.mean(cost[rows, cols])))
        score = base + args.strength_weight*strength_error + args.extra_activity_weight*extra
        meta = {
            "score": float(score),
            "sign_multiplier": flip,
            "assigned_blocks": assigned,
            "center_rms_deg": float(np.sqrt(np.mean(center_errors**2))),
            "center_max_abs_deg": float(np.max(center_errors)),
            "width_log_rms": float(np.sqrt(np.mean(widths**2))),
            "relative_strength_rms": strength_error,
            "extra_activity_ratio": extra,
            "all_lobe_count": len(lobes),
            "significant_lobe_count": len(pool),
        }
        if best is None or score < best[0]:
            best = (score, meta)

    return best



def primary_closure_penalty(v):
    """Continuous full-revolution closure penalty for DE exploration."""
    g, c, r = map(float, v[:3])
    outer_margin = c + r - (g + 1.0)
    inner_margin = abs(g - 1.0) - abs(c - r)
    scale = max(g + c + r + 1.0, 1.0)
    violation = max(0.0, -outer_margin) + max(0.0, -inner_margin)
    return violation / scale

def make_record(v, branch, target, theta, args):
    try:
        c = primary_cadence(v, branch, theta, args.primary_smoothing_deg, args.lobe_epsilon_fraction)
        score, match = score_cadence(c, target["blocks"], args)
    except (ValueError, FloatingPointError, OverflowError):
        return None
    if match is None or score >= 20:
        return None

    return {
        "score": float(score),
        "branch": int(branch),
        "parameters": {k: float(x) for k, x in zip(NAMES, v)},
        "primary": c["primary"],
        "cadence_match": match,
        "minimum_primary_transmission_sine": c["minimum_primary_transmission_sine"],
        "E_span_over_crank": c["E_span_over_crank"],
        "E_path_length_over_crank": c["E_path_length_over_crank"],
        "direction_total_variation_rad": c["direction_total_variation_rad"],
        "normal_acceleration_rms": c["normal_acceleration_rms"],
        "tangential_acceleration_rms": c["tangential_acceleration_rms"],
        "lobes": c["lobes"],
    }


def geom_distance(a, b):
    if a["branch"] != b["branch"]:
        return math.inf
    xa = np.array([a["parameters"][k] for k in NAMES])
    xb = np.array([b["parameters"][k] for k in NAMES])
    d = []
    for i, (lo, hi) in enumerate(BOUNDS):
        if i == 5:
            z = abs((xa[i] - xb[i] + math.pi) % (2*math.pi) - math.pi) / math.pi
        else:
            z = abs(xa[i] - xb[i]) / (hi - lo)
        d.append(z)
    return float(np.sqrt(np.mean(np.square(d))))


def diverse(records, count, min_distance):
    kept = []
    for x in sorted(records, key=lambda z: z["score"]):
        if all(geom_distance(x, y) >= min_distance for y in kept):
            kept.append(x)
            if len(kept) >= count:
                break
    for i, x in enumerate(kept):
        x["library_rank"] = i
    return kept


def save(path, obj):
    def conv(x):
        if isinstance(x, Path): return str(x)
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, dict): return {k: conv(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)): return [conv(v) for v in x]
        return x
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(conv(obj), indent=2) + "\n", encoding="utf-8")


def write_trace(path, best, target, theta, args):
    v = np.array([best["parameters"][k] for k in NAMES])
    c = primary_cadence(v, best["branch"], theta, args.primary_smoothing_deg, args.lobe_epsilon_fraction)
    target_acc = np.interp(np.degrees(theta) % 360.0, target["theta_deg"], target["acceleration"], period=360.0)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        fields = ["theta_deg","E_x","E_y","E_speed","E_speed_smooth",
                  "E_tangential_acceleration","E_normal_acceleration",
                  "target_large_acceleration_smoothed"]
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        for i in range(len(theta)):
            w.writerow({
                "theta_deg": float(np.degrees(theta[i]) % 360.0),
                "E_x": float(c["E"][i,0]), "E_y": float(c["E"][i,1]),
                "E_speed": float(c["speed"][i]), "E_speed_smooth": float(c["speed_smooth"][i]),
                "E_tangential_acceleration": float(c["tangential_acceleration"][i]),
                "E_normal_acceleration": float(c["normal_acceleration"][i]),
                "target_large_acceleration_smoothed": float(target_acc[i]),
            })


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=Path, default=TARGET)
    ap.add_argument("--output", type=Path, default=OUTPUT)
    ap.add_argument("--trace-output", type=Path, default=TRACE)
    ap.add_argument("--islands", type=int, default=24)
    ap.add_argument("--hours", type=float, default=0.0)
    ap.add_argument("--generations", type=int, default=220)
    ap.add_argument("--population-size", type=int, default=16)
    ap.add_argument("--seed", type=int, default=3952100)
    ap.add_argument("--coarse-step-deg", type=float, default=1.0)
    ap.add_argument("--dense-step-deg", type=float, default=0.25)
    ap.add_argument("--target-smoothing-deg", type=float, default=12.0)
    ap.add_argument("--primary-smoothing-deg", type=float, default=8.0)
    ap.add_argument("--lobe-epsilon-fraction", type=float, default=0.01)
    ap.add_argument("--noise-fronts", type=int, default=6)

    # Temporal concordance dominates.
    ap.add_argument("--center-scale-deg", type=float, default=20.0)
    ap.add_argument("--width-weight", type=float, default=0.20)
    ap.add_argument("--strength-weight", type=float, default=0.10)
    ap.add_argument("--sign-weight", type=float, default=0.10)
    ap.add_argument("--extra-activity-weight", type=float, default=0.40)

    ap.add_argument("--primary-sine-floor", type=float, default=0.30)
    ap.add_argument("--E-span-floor", type=float, default=0.75)
    ap.add_argument("--library-size", type=int, default=200)
    ap.add_argument("--diversity-distance", type=float, default=0.08)
    ap.add_argument("--population-candidates-per-island", type=int, default=24)
    args = ap.parse_args()

    if args.islands < 1 or args.hours < 0:
        ap.error("invalid islands/hours")
    for v, name in [(args.coarse_step_deg,"coarse"),(args.dense_step_deg,"dense")]:
        n = round(360.0/v)
        if v <= 0 or not math.isclose(n*v, 360.0, abs_tol=1e-9):
            ap.error(f"--{name}-step-deg must divide 360")

    target = load_target(args.target, args.target_smoothing_deg,
                         args.lobe_epsilon_fraction, args.noise_fronts)

    print("Target LARGE cadence blocks:")
    for i, b in enumerate(target["blocks"], 1):
        print(f"  {i}: {b['start_deg']:.2f}->{b['end_deg']:.2f} deg "
              f"center={b['center_deg']:.2f} span={b['span_deg']:.2f} "
              f"sign={b['sign']:+d} |Dv|={b['strength']:.6f}")
    ratio = target["blocks"][-1]["strength"] / target["all_lobes"][3]["strength"]
    print(f"  3rd/4th lobe strength ratio = {ratio:.3f}")

    nc = round(360.0/args.coarse_step_deg)
    nd = round(360.0/args.dense_step_deg)
    coarse = np.linspace(0, 2*math.pi, nc, endpoint=False)
    dense_theta = np.linspace(0, 2*math.pi, nd, endpoint=False)

    wall0, cpu0 = time.time(), time.process_time()
    deadline = wall0 + args.hours*3600 if args.hours > 0 else None
    report = {
        "description": "L primary cadence library for candidate 3952",
        "length_unit": "primary_crank_radius",
        "target": str(args.target),
        "target_noise_band": target["noise_band"],
        "target_blocks": target["blocks"],
        "target_all_lobes": target["all_lobes"],
        "bounds": dict(zip(NAMES, BOUNDS)),
        "settings": vars(args),
        "islands": [], "candidates": [], "best": None,
    }
    raw_records = []

    island = 0
    while island < args.islands and (deadline is None or time.time() < deadline):
        branch = 1 if island % 2 == 0 else -1
        seed = args.seed + island
        archive = []

        def objective(v):
            v = np.asarray(v, dtype=float)
            closure = primary_closure_penalty(v)
            if closure > 0.0:
                # Unlike the old 15-D global search, infeasible populations
                # still see a gradient pointing toward full-cycle closure.
                return 20.0 + 100.0 * closure
            rec = make_record(v, branch, target, coarse, args)
            if rec is None:
                return 1e3
            archive.append((rec["score"], v.copy()))
            return rec["score"]

        t0, p0 = time.time(), time.process_time()
        opt = differential_evolution(
            objective, BOUNDS, seed=seed, popsize=args.population_size,
            maxiter=args.generations, polish=True, tol=1e-8,
            updating="immediate", workers=1, init="latinhypercube",
            callback=lambda xk, conv: bool(deadline is not None and time.time() >= deadline),
        )

        vectors = [np.asarray(opt.x)]
        if getattr(opt, "population", None) is not None:
            order = np.argsort(np.asarray(opt.population_energies))
            vectors += [np.asarray(opt.population[int(j)]) for j in order[:args.population_candidates_per_island]]
        vectors += [x for _, x in sorted(archive, key=lambda z: z[0])[:24]]

        dense_records, seen = [], set()
        for v in vectors:
            key = tuple(round(float(x), 10) for x in v)
            if key in seen: continue
            seen.add(key)
            rec = make_record(v, branch, target, dense_theta, args)
            if rec is not None:
                dense_records.append(rec); raw_records.append(rec)
        dense_records.sort(key=lambda x: x["score"])
        best = dense_records[0] if dense_records else None

        report["islands"].append({
            "index": island, "branch": branch, "seed": seed,
            "optimizer_success": bool(opt.success),
            "optimizer_message": str(opt.message),
            "function_evaluations": int(opt.nfev),
            "wall_seconds": time.time()-t0,
            "cpu_seconds": time.process_time()-p0,
            "best": best,
        })

        lib = diverse(raw_records, args.library_size, args.diversity_distance)
        report["candidates"], report["best"] = lib, (lib[0] if lib else None)
        report["elapsed_wall_hours"] = (time.time()-wall0)/3600
        report["elapsed_cpu_hours"] = (time.process_time()-cpu0)/3600
        save(args.output, report)

        if best:
            m = best["cadence_match"]
            print(f"island {island:03d} br={branch:+d} score={best['score']:.5f} "
                  f"centerRMS={m['center_rms_deg']:.2f}deg "
                  f"max={m['center_max_abs_deg']:.2f}deg "
                  f"extra={m['extra_activity_ratio']:.3f} "
                  f"sine={best['minimum_primary_transmission_sine']:.3f} "
                  f"Espan={best['E_span_over_crank']:.3f}", flush=True)
        else:
            print(f"island {island:03d} br={branch:+d}: no dense candidate", flush=True)
        island += 1

    lib = diverse(raw_records, args.library_size, args.diversity_distance)
    report["candidates"], report["best"] = lib, (lib[0] if lib else None)
    report["elapsed_wall_hours"] = (time.time()-wall0)/3600
    report["elapsed_cpu_hours"] = (time.process_time()-cpu0)/3600
    report["completed_islands"] = island
    report["stopped_by_wall_time"] = bool(deadline is not None and time.time() >= deadline)
    save(args.output, report)

    print("\nPRIMARY CADENCE SEARCH SUMMARY")
    print(f"wall={report['elapsed_wall_hours']:.3f} h "
          f"cpu={report['elapsed_cpu_hours']:.3f} h "
          f"islands={island} diverse={len(lib)}")
    if not lib:
        print("No feasible primary found."); return

    best = lib[0]; m = best["cadence_match"]
    print(f"best score={best['score']:.6f}\n"
          f"branch={best['branch']:+d}\n"
          f"center RMS={m['center_rms_deg']:.3f} deg\n"
          f"center max={m['center_max_abs_deg']:.3f} deg\n"
          f"extra activity={m['extra_activity_ratio']:.4f}\n"
          f"primary sine={best['minimum_primary_transmission_sine']:.4f}\n"
          f"E span/r={best['E_span_over_crank']:.4f}\n"
          f"AD={best['primary']['ground']:.6f} "
          f"BC={best['primary']['coupler']:.6f} "
          f"CD={best['primary']['rocker']:.6f}")
    print("matched blocks:")
    for t, c in zip(target["blocks"], m["assigned_blocks"]):
        print(f"  target {t['center_deg']:7.2f} <- E {c['center_deg']:7.2f} deg "
              f"span={c['span_deg']:6.2f} sign={c['sign']:+d}")

    write_trace(args.trace_output, best, target, dense_theta, args)
    print(f"Wrote {args.output}")
    print(f"Wrote {args.trace_output}")


if __name__ == "__main__":
    main()
