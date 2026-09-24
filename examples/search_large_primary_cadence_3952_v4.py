#!/usr/bin/env python3
"""Primary cadence search V4 for candidate 3952, LARGE cylinder.

Physical model:
  * two real turnarounds / two monotonic branches;
  * short high->low branch = FAST cadence then SLOW cadence;
  * long low->high return branch = smooth/regular, but NOT a plateau;
  * HP transition itself gets a guard band: fit the cadence on either side,
    not the exact sharpness of the target kink.

V2/V3 candidates and the known old 2.08%-class primary are rescored under V4
and may seed some of 128 broad-bound DE islands. Remaining islands are fresh.

Typical run:
  PYTHONPATH=src python3 examples/search_large_primary_cadence_3952_v4.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import differential_evolution

import search_large_sixbar_stageL1 as stage
import search_motor_hybrid_3952_sixbar as motor3952

ROOT = Path.cwd()
TARGET = ROOT / "outputs" / "motor_hybrid_c2_15p_260k" / "candidate_3952_motion_target.csv"
V2_LIBRARY = ROOT / "outputs" / "large_primary_cadence_3952_v2.json"
V3_LIBRARY = ROOT / "outputs" / "large_primary_cadence_3952_v3.json"
KNOWN_SIXBAR = ROOT / "outputs" / "large_sixbar_3952.json"
OUTPUT = ROOT / "outputs" / "large_primary_cadence_3952_v4.json"
TRACE = ROOT / "outputs" / "large_primary_cadence_3952_v4_best_trace.csv"

NAMES = [
    "primary_ground", "primary_coupler", "primary_rocker",
    "primary_E_along", "primary_E_normal", "primary_phase",
]
BOUNDS = [
    (0.50, 12.0), (0.50, 12.0), (0.35, 6.0),
    (-12.0, 18.0), (-15.0, 15.0), (-math.pi, math.pi),
]


def cyclic_distance_deg(a: float, b: float) -> float:
    return abs((float(a) - float(b) + 180.0) % 360.0 - 180.0)


def smooth_periodic(y: np.ndarray, step_deg: float, width_deg: float) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    if width_deg <= 0.0:
        return y.copy()
    half = max(1, int(round(0.5 * width_deg / step_deg)))
    out = np.zeros_like(y)
    for shift in range(-half, half + 1):
        out += np.roll(y, shift)
    return out / float(2 * half + 1)


def bridge_mask(y: np.ndarray, mask: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    if not np.any(mask):
        return y.copy()
    if np.all(mask):
        raise ValueError("noise mask covers complete cycle")
    n = len(y)
    starts = [i for i in range(n) if mask[i] and not mask[(i - 1) % n]]
    ends = [i for i in range(n) if mask[i] and not mask[(i + 1) % n]]
    if len(starts) != 1 or len(ends) != 1:
        raise ValueError("expected one contiguous cyclic noise band")
    start, end = starts[0], ends[0]
    left, right = (start - 1) % n, (end + 1) % n
    total = (right - left) % n
    out = y.copy()
    for k in range(1, total):
        i = (left + k) % n
        f = k / total
        out[i] = (1.0 - f) * y[left] + f * y[right]
    return out


def principal_axis(E: np.ndarray) -> tuple[np.ndarray, float]:
    X = np.asarray(E, dtype=float) - np.mean(E, axis=0)
    cov = X.T @ X
    values, vectors = np.linalg.eigh(cov)
    i = int(np.argmax(values))
    axis = vectors[:, i] / np.linalg.norm(vectors[:, i])
    if axis[0] < 0.0 or (abs(axis[0]) < 1e-12 and axis[1] < 0.0):
        axis = -axis
    return axis, float(values[i]) / max(float(np.sum(values)), 1e-30)


def interp_periodic(theta_deg: np.ndarray, values: np.ndarray, query_deg):
    return np.interp(
        np.mod(np.asarray(query_deg, dtype=float), 360.0),
        np.asarray(theta_deg, dtype=float),
        np.asarray(values, dtype=float),
        period=360.0,
    )


def zero_crossings_periodic(theta_deg: np.ndarray, velocity: np.ndarray) -> list[float]:
    theta_deg = np.asarray(theta_deg, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    n = len(velocity)
    crossings = []
    for i in range(n):
        j = (i + 1) % n
        vi, vj = float(velocity[i]), float(velocity[j])
        if vi == 0.0:
            crossings.append(float(theta_deg[i] % 360.0))
        elif vi * vj < 0.0:
            a = float(theta_deg[i])
            b = 360.0 if i == n - 1 else float(theta_deg[j])
            f = abs(vi) / (abs(vi) + abs(vj))
            crossings.append(float((a + f * (b - a)) % 360.0))
    merged = []
    for x in sorted(crossings):
        if not merged or cyclic_distance_deg(x, merged[-1]) > 0.25:
            merged.append(x)
    if len(merged) > 1 and cyclic_distance_deg(merged[0], merged[-1]) <= 0.25:
        merged.pop()
    return merged


def cyclic_u(theta_deg: np.ndarray, start_deg: float) -> np.ndarray:
    return (np.asarray(theta_deg, dtype=float) - float(start_deg)) % 360.0


def periodic_integral(theta_deg, values, start_deg, u0, u1, samples_per_degree=4.0):
    span = float(u1 - u0)
    if span <= 0.0:
        return 0.0
    n = max(8, int(math.ceil(span * samples_per_degree)) + 1)
    u = np.linspace(u0, u1, n)
    deg = (float(start_deg) + u) % 360.0
    y = interp_periodic(theta_deg, values, deg)
    return float(np.trapz(y, np.radians(u)))


def discover_short_split(theta_deg, velocity, short_start_deg, short_span_deg, args):
    """Best ordered two-level speed split on the target short branch.

    This identifies two cadences without requiring either cadence to be flat.
    Turnaround margins are excluded from the split fit.
    """
    step = 0.25
    u = np.arange(
        args.turnaround_guard_deg,
        short_span_deg - args.turnaround_guard_deg + 1e-12,
        step,
    )
    speed = np.abs(interp_periodic(theta_deg, velocity, short_start_deg + u))
    best = None
    for split in np.arange(
        args.minimum_short_subphase_deg,
        short_span_deg - args.minimum_short_subphase_deg + 1e-12,
        step,
    ):
        left, right = u <= split, u > split
        if np.count_nonzero(left) < 8 or np.count_nonzero(right) < 8:
            continue
        m1, m2 = float(np.mean(speed[left])), float(np.mean(speed[right]))
        sse = float(np.sum((speed[left] - m1) ** 2) + np.sum((speed[right] - m2) ** 2))
        if best is None or sse < best["sse"]:
            best = {
                "sse": sse,
                "split_offset_deg": float(split),
                "split_deg": float((short_start_deg + split) % 360.0),
            }
    if best is None:
        raise RuntimeError("could not discover target fast/slow split")
    return best


def load_target(path: Path, args) -> dict:
    raw = np.genfromtxt(path, delimiter=",", names=True)
    if raw["theta_rad"][-1] >= 2.0 * math.pi - 1e-9:
        raw = raw[:-1]
    theta = np.asarray(raw["theta_rad"], dtype=float)
    theta_deg = np.mod(np.degrees(theta), 360.0)
    q = np.asarray(raw["large_fraction_0_1"], dtype=float)
    dq = np.asarray(raw["large_dq_dtheta_per_rad"], dtype=float)

    band = motor3952._detect_noise_bands(path, args.noise_fronts)["large"]
    noise_mask = motor3952._theta_band_mask(theta, band)
    dq_bridge = bridge_mask(dq, noise_mask)
    velocity = smooth_periodic(dq_bridge, 360.0 / len(theta), args.target_smoothing_deg)

    crossings = zero_crossings_periodic(theta_deg, velocity)
    if len(crossings) != 2:
        raise RuntimeError(f"target should have exactly two turnarounds, got {crossings}")
    q_at = [float(interp_periodic(theta_deg, q, x)) for x in crossings]
    i_high = int(np.argmax(q_at))
    high_turn = float(crossings[i_high])
    low_turn = float(crossings[1 - i_high])
    short_span = (low_turn - high_turn) % 360.0
    long_span = 360.0 - short_span
    if short_span > long_span:
        raise RuntimeError("high->low target branch unexpectedly longer than return branch")

    split = discover_short_split(theta_deg, velocity, high_turn, short_span, args)
    split_u = split["split_offset_deg"]
    fast_lo = args.turnaround_guard_deg
    fast_hi = split_u - args.transition_guard_deg
    slow_lo = split_u + args.transition_guard_deg
    slow_hi = short_span - args.turnaround_guard_deg

    u = cyclic_u(theta_deg, high_turn)
    fast_mask = (u >= fast_lo) & (u <= fast_hi)
    slow_mask = (u >= slow_lo) & (u <= slow_hi)
    if np.count_nonzero(fast_mask) < 5 or np.count_nonzero(slow_mask) < 5:
        raise RuntimeError("target fast/slow core masks too small")

    fast_mean = float(np.mean(np.abs(velocity[fast_mask])))
    slow_mean = float(np.mean(np.abs(velocity[slow_mask])))
    speed_ratio = fast_mean / max(slow_mean, 1e-12)

    fast_disp = abs(periodic_integral(theta_deg, velocity, high_turn, 0.0, split_u))
    slow_disp = abs(periodic_integral(theta_deg, velocity, high_turn, split_u, short_span))
    fast_disp_fraction = fast_disp / max(fast_disp + slow_disp, 1e-12)

    profile_u = np.linspace(
        short_span + args.turnaround_guard_deg,
        360.0 - args.turnaround_guard_deg,
        args.long_profile_points,
    )
    profile_deg = (high_turn + profile_u) % 360.0
    long_profile = interp_periodic(theta_deg, velocity, profile_deg)
    if float(np.mean(long_profile)) <= 0.0:
        raise RuntimeError("target long return branch is not positive")
    long_profile_norm = long_profile / max(float(np.mean(long_profile)), 1e-12)

    return {
        "theta": theta,
        "theta_deg": theta_deg,
        "q": q,
        "velocity_smooth": velocity,
        "noise_band": band,
        "turnarounds_deg": {"high_position": high_turn, "low_position": low_turn},
        "short_span_deg": float(short_span),
        "long_span_deg": float(long_span),
        "short_split": split,
        "short_fast_core": [float(fast_lo), float(fast_hi)],
        "short_slow_core": [float(slow_lo), float(slow_hi)],
        "target_fast_mean_speed": fast_mean,
        "target_slow_mean_speed": slow_mean,
        "target_fast_slow_speed_ratio": float(speed_ratio),
        "target_fast_displacement_fraction": float(fast_disp_fraction),
        "long_profile_u": profile_u,
        "long_profile_norm": long_profile_norm,
    }


def primary_features(v: np.ndarray, branch: int, theta: np.ndarray, args) -> dict:
    E, Ed, sine, pdata = stage.primary(theta, v, branch)
    theta_deg = np.mod(np.degrees(theta), 360.0)
    axis, axis_fraction = principal_axis(E)
    projected_velocity = smooth_periodic(
        Ed @ axis,
        360.0 / len(theta),
        args.primary_smoothing_deg,
    )
    X = E - np.mean(E, axis=0)
    normal = np.array((-axis[1], axis[0]))
    along, lateral = X @ axis, X @ normal
    return {
        "theta_deg": theta_deg,
        "E": E,
        "axis_angle_deg": math.degrees(math.atan2(axis[1], axis[0])),
        "axis_position_variance_fraction": axis_fraction,
        "projected_velocity": projected_velocity,
        "minimum_primary_transmission_sine": float(sine),
        "primary": pdata,
        "E_span_over_crank": max(float(np.ptp(E[:, 0])), float(np.ptp(E[:, 1]))),
        "E_principal_span_over_crank": float(np.ptp(along)),
        "E_lateral_span_over_crank": float(np.ptp(lateral)),
        "E_path_length_over_crank": float(np.sum(np.linalg.norm(np.roll(E, -1, axis=0) - E, axis=1))),
        "BE_over_BC": float(math.hypot(float(pdata["E_along"]), float(pdata["E_normal"])) / max(float(pdata["coupler"]), 1e-12)),
    }


def turning_metrics(theta_deg, velocity, target, sign_multiplier, args):
    v = sign_multiplier * np.asarray(velocity, dtype=float)
    scale = max(float(np.sqrt(np.mean(v * v))), 1e-10)
    high = target["turnarounds_deg"]["high_position"]
    low = target["turnarounds_deg"]["low_position"]
    short_span = target["short_span_deg"]
    u = cyclic_u(theta_deg, high)
    g = args.turnaround_sign_guard_deg
    short_core = (u >= g) & (u <= short_span - g)
    long_core = (u >= short_span + g) & (u <= 360.0 - g)

    wrong_short = np.maximum(v[short_core], 0.0)   # should be negative
    wrong_long = np.maximum(-v[long_core], 0.0)    # should be positive
    monotonicity = math.sqrt(
        (float(np.sum(wrong_short ** 2)) + float(np.sum(wrong_long ** 2)))
        / max(np.count_nonzero(short_core) + np.count_nonzero(long_core), 1)
    ) / scale

    endpoint_v = interp_periodic(theta_deg, v, np.array([high, low]))
    endpoint_zero = float(np.sqrt(np.mean(endpoint_v ** 2)) / scale)

    crossings = zero_crossings_periodic(theta_deg, v)
    if len(crossings) >= 2:
        best = None
        for i in range(len(crossings)):
            for j in range(len(crossings)):
                if i == j:
                    continue
                err = np.array([
                    cyclic_distance_deg(crossings[i], high),
                    cyclic_distance_deg(crossings[j], low),
                ])
                rms = float(np.sqrt(np.mean(err ** 2)))
                if best is None or rms < best[0]:
                    best = (rms, i, j)
        turning_rms = float(best[0])
        matched = [float(crossings[best[1]]), float(crossings[best[2]])]
    else:
        turning_rms, matched = 180.0, []

    return {
        "monotonicity_wrong_sign_rms": float(monotonicity),
        "endpoint_zero_rms": endpoint_zero,
        "zero_crossings_deg": [float(x) for x in crossings],
        "zero_crossing_count": len(crossings),
        "matched_turnarounds_deg": matched,
        "turning_rms_deg": turning_rms,
        "extra_crossings": max(0, len(crossings) - 2),
    }


def cadence_metrics(features, target, sign_multiplier):
    theta_deg = features["theta_deg"]
    v = sign_multiplier * features["projected_velocity"]
    high = target["turnarounds_deg"]["high_position"]
    short_span = target["short_span_deg"]
    split_u = target["short_split"]["split_offset_deg"]
    u = cyclic_u(theta_deg, high)
    fast_lo, fast_hi = target["short_fast_core"]
    slow_lo, slow_hi = target["short_slow_core"]
    fast_mask = (u >= fast_lo) & (u <= fast_hi)
    slow_mask = (u >= slow_lo) & (u <= slow_hi)

    fast_mean = float(np.mean(np.abs(v[fast_mask])))
    slow_mean = float(np.mean(np.abs(v[slow_mask])))
    ratio = fast_mean / max(slow_mean, 1e-12)
    ratio_error = abs(math.log(max(ratio, 1e-12) / target["target_fast_slow_speed_ratio"]))

    fast_disp = abs(periodic_integral(theta_deg, v, high, 0.0, split_u))
    slow_disp = abs(periodic_integral(theta_deg, v, high, split_u, short_span))
    fraction = fast_disp / max(fast_disp + slow_disp, 1e-12)
    fraction_error = abs(fraction - target["target_fast_displacement_fraction"])

    profile_deg = (high + target["long_profile_u"]) % 360.0
    cand_profile = interp_periodic(theta_deg, v, profile_deg)
    cand_mean = float(np.mean(cand_profile))
    if cand_mean <= 1e-12:
        long_rms = 2.0
    else:
        cand_norm = cand_profile / cand_mean
        long_rms = float(np.sqrt(np.mean((cand_norm - target["long_profile_norm"]) ** 2)))

    return {
        "fast_mean_speed": fast_mean,
        "slow_mean_speed": slow_mean,
        "fast_slow_speed_ratio": float(ratio),
        "speed_ratio_log_error": float(ratio_error),
        "fast_displacement_fraction": float(fraction),
        "displacement_fraction_error": float(fraction_error),
        "long_profile_normalized_rms": float(long_rms),
    }


def cadence_score(features, target, args):
    sine = features["minimum_primary_transmission_sine"]
    if sine < args.primary_sine_floor:
        return 20.0 + 100.0 * (args.primary_sine_floor - sine) / max(args.primary_sine_floor, 1e-9), None
    if features["E_span_over_crank"] < args.E_span_floor:
        return 20.0 + 100.0 * (args.E_span_floor - features["E_span_over_crank"]) / max(args.E_span_floor, 1e-9), None

    best = None
    for sign_multiplier in (+1, -1):
        turn = turning_metrics(features["theta_deg"], features["projected_velocity"], target, sign_multiplier, args)
        cadence = cadence_metrics(features, target, sign_multiplier)
        turning_term = turn["turning_rms_deg"] / max(args.turning_scale_deg, 1e-9)
        extra_penalty = args.extra_crossing_weight * turn["extra_crossings"]

        sine_penalty = 0.0 if sine >= args.preferred_primary_sine else (
            args.sine_preference_weight * (args.preferred_primary_sine - sine)
            / max(args.preferred_primary_sine - args.primary_sine_floor, 1e-9)
        )
        axis_fraction = features["axis_position_variance_fraction"]
        axis_penalty = 0.0 if axis_fraction >= args.preferred_axis_variance_fraction else (
            args.axis_preference_weight
            * (args.preferred_axis_variance_fraction - axis_fraction)
            / max(args.preferred_axis_variance_fraction, 1e-9)
        )

        score = (
            args.monotonicity_weight * turn["monotonicity_wrong_sign_rms"]
            + args.endpoint_zero_weight * turn["endpoint_zero_rms"]
            + args.turning_weight * turning_term
            + extra_penalty
            + args.speed_ratio_weight * cadence["speed_ratio_log_error"]
            + args.displacement_fraction_weight * cadence["displacement_fraction_error"]
            + args.long_shape_weight * cadence["long_profile_normalized_rms"]
            + sine_penalty + axis_penalty
        )
        meta = {
            "score": float(score),
            "sign_multiplier": int(sign_multiplier),
            **turn, **cadence,
            "sine_preference_penalty": float(sine_penalty),
            "axis_preference_penalty": float(axis_penalty),
        }
        if best is None or score < best[0]:
            best = (float(score), meta)
    return best


def primary_closure_penalty(v: np.ndarray) -> float:
    g, c, r = map(float, v[:3])
    outer_margin = c + r - (g + 1.0)
    inner_margin = abs(g - 1.0) - abs(c - r)
    scale = max(g + c + r + 1.0, 1.0)
    return (max(0.0, -outer_margin) + max(0.0, -inner_margin)) / scale


def make_record(v, branch, target, theta, args):
    try:
        features = primary_features(v, branch, theta, args)
        score, match = cadence_score(features, target, args)
    except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
        return None
    if match is None or score >= 20.0:
        return None
    return {
        "score": float(score),
        "branch": int(branch),
        "parameters": {name: float(x) for name, x in zip(NAMES, v)},
        "primary": features["primary"],
        "cadence_match": match,
        "minimum_primary_transmission_sine": features["minimum_primary_transmission_sine"],
        "projection_axis_angle_deg": features["axis_angle_deg"],
        "axis_position_variance_fraction": features["axis_position_variance_fraction"],
        "E_span_over_crank": features["E_span_over_crank"],
        "E_principal_span_over_crank": features["E_principal_span_over_crank"],
        "E_lateral_span_over_crank": features["E_lateral_span_over_crank"],
        "E_path_length_over_crank": features["E_path_length_over_crank"],
        "BE_over_BC": features["BE_over_BC"],
    }


def vector_from_record(item):
    if "parameters" in item and all(name in item["parameters"] for name in NAMES):
        v = np.array([float(item["parameters"][name]) for name in NAMES], dtype=float)
        branch = int(item.get("branch", item.get("primary", {}).get("assembly_branch", 0)))
        if branch not in (-1, 1):
            raise ValueError("invalid primary branch")
        return v, branch
    p = item["primary"]
    return np.array([
        float(p["ground"]), float(p["coupler"]), float(p["rocker"]),
        float(p["E_along"]), float(p["E_normal"]), float(p["phase"]),
    ]), int(p["assembly_branch"])


def geometry_distance(a, b):
    if int(a["branch"]) != int(b["branch"]):
        return math.inf
    xa = np.array([a["parameters"][k] for k in NAMES], dtype=float)
    xb = np.array([b["parameters"][k] for k in NAMES], dtype=float)
    d = []
    for i, (lo, hi) in enumerate(BOUNDS):
        if i == 5:
            z = abs((xa[i] - xb[i] + math.pi) % (2.0 * math.pi) - math.pi) / math.pi
        else:
            z = abs(xa[i] - xb[i]) / (hi - lo)
        d.append(z)
    return float(np.sqrt(np.mean(np.square(d))))


def diverse(records, count, minimum_distance):
    kept = []
    for item in sorted(records, key=lambda x: x["score"]):
        if all(geometry_distance(item, other) >= minimum_distance for other in kept):
            kept.append(item)
            if len(kept) >= count:
                break
    for rank, item in enumerate(kept):
        item["library_rank"] = rank
    return kept


def load_json_relaxed(path: Path):
    return json.loads(path.read_text(encoding="utf-8").replace("NaN", "null"))


def source_records(path: Path, prefix: str):
    if not path.exists():
        return []
    data = load_json_relaxed(path)
    out = []
    for item in data.get("candidates", []):
        if item is not None:
            out.append((f"{prefix}_library", item))
    for island in data.get("islands", []):
        item = island.get("best")
        if item is not None:
            out.append((f"{prefix}_island_{int(island.get('index', -1))}", item))
    return out


def rescore_seed_sources(paths, known_sixbar_path, target, dense_theta, args):
    items = []
    for path, prefix in paths:
        items.extend(source_records(path, prefix))
    if known_sixbar_path.exists():
        data = load_json_relaxed(known_sixbar_path)
        if data.get("best"):
            items.append(("known_old_2p08_primary", data["best"]))

    rescored, seen = [], set()
    for source, item in items:
        try:
            v, branch = vector_from_record(item)
        except (KeyError, TypeError, ValueError):
            continue
        key = (branch, *[round(float(x), 10) for x in v])
        if key in seen:
            continue
        seen.add(key)
        record = make_record(v, branch, target, dense_theta, args)
        if record is not None:
            record["seed_source"] = source
            rescored.append(record)
    return sorted(rescored, key=lambda x: x["score"])


def select_seed_panel(rescored, count, minimum_distance):
    kept = []
    for item in rescored:
        if all(geometry_distance(item, other) >= minimum_distance for other in kept):
            kept.append(item)
            if len(kept) >= count:
                break
    return kept


def build_jobs(islands, seeds):
    jobs = []
    for item in seeds[:islands]:
        v, branch = vector_from_record(item)
        jobs.append({
            "branch": branch, "x0": v,
            "seed_source": item.get("seed_source", "historical"),
            "seed_v4_score": float(item["score"]),
        })
    target_plus = islands // 2
    target_minus = islands - target_plus
    plus = sum(j["branch"] == 1 for j in jobs)
    minus = sum(j["branch"] == -1 for j in jobs)
    while len(jobs) < islands:
        if plus < target_plus and (plus <= minus or minus >= target_minus):
            branch, plus = 1, plus + 1
        elif minus < target_minus:
            branch, minus = -1, minus + 1
        else:
            branch = 1 if plus < target_plus else -1
            plus += int(branch == 1)
            minus += int(branch == -1)
        jobs.append({"branch": branch, "x0": None, "seed_source": None, "seed_v4_score": None})
    return jobs


def safe_json_value(x):
    if isinstance(x, Path):
        return str(x)
    if isinstance(x, np.generic):
        x = x.item()
    if isinstance(x, float):
        return x if math.isfinite(x) else None
    if isinstance(x, np.ndarray):
        return [safe_json_value(v) for v in x.tolist()]
    if isinstance(x, dict):
        return {k: safe_json_value(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [safe_json_value(v) for v in x]
    return x


def save(path: Path, obj: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(safe_json_value(obj), indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_trace(path, best, target, dense_theta, args):
    v = np.array([best["parameters"][k] for k in NAMES], dtype=float)
    f = primary_features(v, int(best["branch"]), dense_theta, args)
    sign = int(best["cadence_match"]["sign_multiplier"])
    cv = sign * f["projected_velocity"]
    cand_scale = max(float(np.mean(np.abs(cv))), 1e-12)
    target_scale = max(float(np.mean(np.abs(target["velocity_smooth"]))), 1e-12)
    high = target["turnarounds_deg"]["high_position"]
    u = cyclic_u(f["theta_deg"], high)
    split_u = target["short_split"]["split_offset_deg"]
    short_span = target["short_span_deg"]
    tv = interp_periodic(target["theta_deg"], target["velocity_smooth"], f["theta_deg"])

    def phase_name(ui):
        if ui <= split_u:
            return "short_fast"
        if ui <= short_span:
            return "short_slow"
        return "long_regular"

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        fields = [
            "theta_deg", "phase", "u_from_high_turn_deg", "E_x", "E_y",
            "candidate_projected_velocity", "candidate_velocity_normalized",
            "target_velocity", "target_velocity_normalized",
        ]
        w = csv.DictWriter(stream, fieldnames=fields)
        w.writeheader()
        for i in range(len(dense_theta)):
            w.writerow({
                "theta_deg": float(f["theta_deg"][i]),
                "phase": phase_name(float(u[i])),
                "u_from_high_turn_deg": float(u[i]),
                "E_x": float(f["E"][i, 0]), "E_y": float(f["E"][i, 1]),
                "candidate_projected_velocity": float(cv[i]),
                "candidate_velocity_normalized": float(cv[i] / cand_scale),
                "target_velocity": float(tv[i]),
                "target_velocity_normalized": float(tv[i] / target_scale),
            })


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=Path, default=TARGET)
    ap.add_argument("--v2-library", type=Path, default=V2_LIBRARY)
    ap.add_argument("--v3-library", type=Path, default=V3_LIBRARY)
    ap.add_argument("--known-sixbar", type=Path, default=KNOWN_SIXBAR)
    ap.add_argument("--output", type=Path, default=OUTPUT)
    ap.add_argument("--trace-output", type=Path, default=TRACE)

    ap.add_argument("--islands", type=int, default=128)
    ap.add_argument("--generations", type=int, default=260)
    ap.add_argument("--population-size", type=int, default=16)
    ap.add_argument("--seed", type=int, default=3952900)
    ap.add_argument("--coarse-step-deg", type=float, default=1.0)
    ap.add_argument("--dense-step-deg", type=float, default=0.25)
    ap.add_argument("--target-smoothing-deg", type=float, default=12.0)
    ap.add_argument("--primary-smoothing-deg", type=float, default=8.0)
    ap.add_argument("--noise-fronts", type=int, default=6)

    ap.add_argument("--turnaround-guard-deg", type=float, default=8.0)
    ap.add_argument("--turnaround-sign-guard-deg", type=float, default=3.0)
    ap.add_argument("--transition-guard-deg", type=float, default=6.0)
    ap.add_argument("--minimum-short-subphase-deg", type=float, default=25.0)
    ap.add_argument("--long-profile-points", type=int, default=64)

    ap.add_argument("--monotonicity-weight", type=float, default=2.0)
    ap.add_argument("--endpoint-zero-weight", type=float, default=0.40)
    ap.add_argument("--turning-weight", type=float, default=0.65)
    ap.add_argument("--turning-scale-deg", type=float, default=15.0)
    ap.add_argument("--extra-crossing-weight", type=float, default=0.50)
    ap.add_argument("--speed-ratio-weight", type=float, default=0.55)
    ap.add_argument("--displacement-fraction-weight", type=float, default=0.80)
    ap.add_argument("--long-shape-weight", type=float, default=0.45)

    ap.add_argument("--primary-sine-floor", type=float, default=0.30)
    ap.add_argument("--preferred-primary-sine", type=float, default=0.35)
    ap.add_argument("--sine-preference-weight", type=float, default=0.12)
    ap.add_argument("--preferred-axis-variance-fraction", type=float, default=0.70)
    ap.add_argument("--axis-preference-weight", type=float, default=0.03)
    ap.add_argument("--E-span-floor", type=float, default=0.75)

    ap.add_argument("--historical-seed-count", type=int, default=32)
    ap.add_argument("--seed-diversity-distance", type=float, default=0.025)
    ap.add_argument("--library-size", type=int, default=400)
    ap.add_argument("--diversity-distance", type=float, default=0.04)
    ap.add_argument("--population-candidates-per-island", type=int, default=20)
    args = ap.parse_args()

    if args.islands < 1:
        ap.error("--islands must be >= 1")
    if not (0.0 < args.primary_sine_floor <= args.preferred_primary_sine <= 1.0):
        ap.error("invalid sine thresholds")

    target = load_target(args.target, args)
    print("TARGET LARGE V4 PHYSICAL STRUCTURE", flush=True)
    print(f"  high-position turn = {target['turnarounds_deg']['high_position']:.3f} deg", flush=True)
    print(f"  low-position turn  = {target['turnarounds_deg']['low_position']:.3f} deg", flush=True)
    print(f"  short branch       = {target['short_span_deg']:.3f} deg", flush=True)
    print(f"  long branch        = {target['long_span_deg']:.3f} deg", flush=True)
    print(
        f"  fast->slow split   = {target['short_split']['split_deg']:.3f} deg "
        f"(+{target['short_split']['split_offset_deg']:.3f} deg from high turn)",
        flush=True,
    )
    print(f"  fast/slow core speed ratio = {target['target_fast_slow_speed_ratio']:.4f}", flush=True)
    print(f"  fast displacement fraction = {target['target_fast_displacement_fraction']:.4f}", flush=True)
    print("  long branch = normalized smooth target profile (NOT a plateau)", flush=True)

    nc = int(round(360.0 / args.coarse_step_deg))
    nd = int(round(360.0 / args.dense_step_deg))
    coarse_theta = np.linspace(0.0, 2.0 * math.pi, nc, endpoint=False)
    dense_theta = np.linspace(0.0, 2.0 * math.pi, nd, endpoint=False)

    rescored = rescore_seed_sources(
        [(args.v3_library, "v3"), (args.v2_library, "v2")],
        args.known_sixbar,
        target,
        dense_theta,
        args,
    )
    seed_panel = select_seed_panel(
        rescored,
        min(args.historical_seed_count, args.islands),
        args.seed_diversity_distance,
    )
    print(
        f"HISTORICAL candidates rescored={len(rescored)}, selected as V4 seeds={len(seed_panel)}",
        flush=True,
    )
    for i, item in enumerate(seed_panel[:16]):
        m = item["cadence_match"]
        print(
            f"  seed {i:02d}: branch={item['branch']:+d} V4={item['score']:.5f} "
            f"source={item.get('seed_source')} turnRMS={m['turning_rms_deg']:.2f}deg "
            f"ratio={m['fast_slow_speed_ratio']:.2f} sine={item['minimum_primary_transmission_sine']:.3f} "
            f"BE/BC={item['BE_over_BC']:.2f}",
            flush=True,
        )
    if len(seed_panel) > 16:
        print("  ...", flush=True)

    jobs = build_jobs(args.islands, seed_panel)
    wall0, cpu0 = time.time(), time.process_time()
    report = {
        "description": (
            "Candidate 3952 LARGE primary cadence V4: two physical monotonic branches; "
            "short branch split into fast/slow cadences; long branch regular but not plateau-like."
        ),
        "length_unit": "primary_crank_radius",
        "target": str(args.target),
        "target_noise_band": target["noise_band"],
        "target_structure": {
            "turnarounds_deg": target["turnarounds_deg"],
            "short_span_deg": target["short_span_deg"],
            "long_span_deg": target["long_span_deg"],
            "short_split": target["short_split"],
            "short_fast_core": target["short_fast_core"],
            "short_slow_core": target["short_slow_core"],
            "fast_slow_speed_ratio": target["target_fast_slow_speed_ratio"],
            "fast_displacement_fraction": target["target_fast_displacement_fraction"],
        },
        "bounds": dict(zip(NAMES, BOUNDS)),
        "settings": vars(args),
        "historical_rescored": rescored,
        "historical_seed_panel": seed_panel,
        "islands": [], "candidates": [], "best": None,
    }
    raw_records = []

    for island, job in enumerate(jobs):
        branch = int(job["branch"])
        x0 = job["x0"]
        run_seed = args.seed + island
        archive = []

        def objective(v):
            v = np.asarray(v, dtype=float)
            closure = primary_closure_penalty(v)
            if closure > 0.0:
                return 20.0 + 100.0 * closure
            try:
                features = primary_features(v, branch, coarse_theta, args)
                score, match = cadence_score(features, target, args)
            except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
                return 1e3
            if match is None:
                return float(score)
            archive.append((float(score), v.copy()))
            return float(score)

        t0, p0 = time.time(), time.process_time()
        opt = differential_evolution(
            objective, BOUNDS, seed=run_seed,
            popsize=args.population_size, maxiter=args.generations,
            polish=True, tol=1e-8, updating="immediate", workers=1,
            init="latinhypercube", x0=x0,
        )

        vectors = [np.asarray(opt.x, dtype=float)]
        if getattr(opt, "population", None) is not None:
            order = np.argsort(np.asarray(opt.population_energies, dtype=float))
            vectors += [
                np.asarray(opt.population[int(j)], dtype=float)
                for j in order[: args.population_candidates_per_island]
            ]
        vectors += [
            x for _, x in sorted(archive, key=lambda z: z[0])[: args.population_candidates_per_island]
        ]

        dense, seen = [], set()
        for v in vectors:
            key = tuple(round(float(x), 10) for x in v)
            if key in seen:
                continue
            seen.add(key)
            record = make_record(v, branch, target, dense_theta, args)
            if record is not None:
                record["island_index"] = island
                record["island_seed_source"] = job["seed_source"]
                dense.append(record)
                raw_records.append(record)
        dense.sort(key=lambda x: x["score"])
        best = dense[0] if dense else None

        report["islands"].append({
            "index": island, "branch": branch, "optimizer_seed": run_seed,
            "seed_source": job["seed_source"], "seed_v4_score": job["seed_v4_score"],
            "optimizer_success": bool(opt.success), "optimizer_message": str(opt.message),
            "function_evaluations": int(opt.nfev),
            "wall_seconds": time.time() - t0, "cpu_seconds": time.process_time() - p0,
            "best": best,
        })
        library = diverse(raw_records, args.library_size, args.diversity_distance)
        report["candidates"] = library
        report["best"] = library[0] if library else None
        report["elapsed_wall_hours"] = (time.time() - wall0) / 3600.0
        report["elapsed_cpu_hours"] = (time.process_time() - cpu0) / 3600.0
        save(args.output, report)

        if best is None:
            print(f"island {island:03d} branch={branch:+d}: no dense candidate", flush=True)
        else:
            m = best["cadence_match"]
            print(
                f"island {island:03d} branch={branch:+d} "
                f"{'seeded' if job['seed_source'] else 'fresh ':6s}: "
                f"score={best['score']:.5f} turn={m['turning_rms_deg']:.2f}deg "
                f"mono={m['monotonicity_wrong_sign_rms']:.4f} ratio={m['fast_slow_speed_ratio']:.2f} "
                f"dispFast={m['fast_displacement_fraction']:.3f} long={m['long_profile_normalized_rms']:.3f} "
                f"zeros={m['zero_crossing_count']} sine={best['minimum_primary_transmission_sine']:.3f} "
                f"BE/BC={best['BE_over_BC']:.2f}",
                flush=True,
            )

    library = diverse(raw_records, args.library_size, args.diversity_distance)
    report["candidates"] = library
    report["best"] = library[0] if library else None
    report["completed_islands"] = len(jobs)
    report["elapsed_wall_hours"] = (time.time() - wall0) / 3600.0
    report["elapsed_cpu_hours"] = (time.process_time() - cpu0) / 3600.0
    save(args.output, report)

    print("\nPRIMARY CADENCE V4 SUMMARY", flush=True)
    print(
        f"wall={report['elapsed_wall_hours']:.3f} h cpu={report['elapsed_cpu_hours']:.3f} h "
        f"islands={len(jobs)} diverse={len(library)}",
        flush=True,
    )
    if not library:
        print("No feasible V4 primary found.", flush=True)
        return

    best = library[0]
    m = best["cadence_match"]
    print(
        f"best score={best['score']:.6f}\n"
        f"branch={best['branch']:+d}\n"
        f"turning RMS={m['turning_rms_deg']:.3f} deg\n"
        f"zero crossings={m['zero_crossing_count']}\n"
        f"monotonicity wrong-sign RMS={m['monotonicity_wrong_sign_rms']:.6f}\n"
        f"fast/slow ratio={m['fast_slow_speed_ratio']:.6f} "
        f"(target {target['target_fast_slow_speed_ratio']:.6f})\n"
        f"fast displacement fraction={m['fast_displacement_fraction']:.6f} "
        f"(target {target['target_fast_displacement_fraction']:.6f})\n"
        f"long profile nRMS={m['long_profile_normalized_rms']:.6f}\n"
        f"primary sine={best['minimum_primary_transmission_sine']:.6f}\n"
        f"axis variance fraction={best['axis_position_variance_fraction']:.6f}\n"
        f"BE/BC={best['BE_over_BC']:.6f}\n"
        f"AD={best['primary']['ground']:.6f} BC={best['primary']['coupler']:.6f} "
        f"CD={best['primary']['rocker']:.6f}",
        flush=True,
    )
    write_trace(args.trace_output, best, target, dense_theta, args)
    print(f"Wrote {args.output}", flush=True)
    print(f"Wrote {args.trace_output}", flush=True)


if __name__ == "__main__":
    main()
