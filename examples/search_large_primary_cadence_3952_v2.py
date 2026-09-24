#!/usr/bin/env python3
"""Primary-cadence search V2 for candidate 3952, LARGE cylinder.

V1 compared the tangential acceleration d|E'|/dtheta of a 2-D point E with
the signed piston acceleration d²q/dtheta².  Those are not equivalent
quantities.  V2 instead extracts a signed scalar motion from E by projecting
its velocity onto the principal axis of E's trajectory.

The primary is rewarded for:
  1. a few dominant projected-acceleration blocks with roughly the target timing;
  2. simple motion BETWEEN those blocks: projected velocity should be close to
     a straight line versus crank angle, which covers both nearly constant
     velocity and nearly constant acceleration;
  3. little extra acceleration activity;
  4. acceptable four-bar transmission.

Timing is deliberately softer than in V1.  The downstream dyad is expected to
transform the primary motion substantially.

Typical run:
    PYTHONPATH=src python3 examples/search_large_primary_cadence_3952_v2.py

Longer exploration:
    PYTHONPATH=src python3 examples/search_large_primary_cadence_3952_v2.py \
        --islands 48 --generations 300 --population-size 16

Outputs:
    outputs/large_primary_cadence_3952_v2.json
    outputs/large_primary_cadence_3952_v2_best_trace.csv
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
TARGET = (
    ROOT
    / "outputs"
    / "motor_hybrid_c2_15p_260k"
    / "candidate_3952_motion_target.csv"
)
OUTPUT = ROOT / "outputs" / "large_primary_cadence_3952_v2.json"
TRACE = ROOT / "outputs" / "large_primary_cadence_3952_v2_best_trace.csv"
KNOWN_SIXBAR = ROOT / "outputs" / "large_sixbar_3952.json"

NAMES = [
    "primary_ground",
    "primary_coupler",
    "primary_rocker",
    "primary_E_along",
    "primary_E_normal",
    "primary_phase",
]

BOUNDS = [
    (0.50, 12.0),
    (0.50, 12.0),
    (0.35, 6.0),
    (-12.0, 18.0),
    (-15.0, 15.0),
    (-math.pi, math.pi),
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


def derivative_periodic(y: np.ndarray, step_rad: float) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    return (np.roll(y, -1) - np.roll(y, 1)) / (2.0 * step_rad)


def bridge_mask(y: np.ndarray, mask: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    if not np.any(mask):
        return y.copy()
    if np.all(mask):
        raise ValueError("noise mask covers full cycle")

    n = len(y)
    starts = [i for i in range(n) if mask[i] and not mask[(i - 1) % n]]
    ends = [i for i in range(n) if mask[i] and not mask[(i + 1) % n]]
    if len(starts) != 1 or len(ends) != 1:
        raise ValueError("expected one cyclic noise interval")

    start, end = starts[0], ends[0]
    left, right = (start - 1) % n, (end + 1) % n
    total = (right - left) % n
    out = y.copy()
    for k in range(1, total):
        i = (left + k) % n
        f = k / total
        out[i] = (1.0 - f) * y[left] + f * y[right]
    return out


def extract_lobes(
    theta_deg: np.ndarray,
    acceleration: np.ndarray,
    epsilon_fraction: float,
) -> list[dict]:
    theta_deg = np.asarray(theta_deg, dtype=float)
    acceleration = np.asarray(acceleration, dtype=float)
    n = len(acceleration)
    step_deg = 360.0 / n
    step_rad = 2.0 * math.pi / n

    peak = float(np.max(np.abs(acceleration)))
    if not math.isfinite(peak) or peak <= 1e-14:
        return []

    eps = epsilon_fraction * peak
    signs = np.zeros(n, dtype=np.int8)
    signs[acceleration > eps] = 1
    signs[acceleration < -eps] = -1

    # Close only tiny holes inside one signed block.
    for _ in range(4):
        old = signs.copy()
        for i in range(n):
            if old[i] != 0:
                continue
            left = old[(i - 1) % n]
            right = old[(i + 1) % n]
            if left != 0 and left == right:
                signs[i] = left

    segments = []
    start = 0
    current = int(signs[0])
    for i in range(1, n):
        s = int(signs[i])
        if s != current:
            segments.append([start, i - 1, current])
            start, current = i, s
    segments.append([start, n - 1, current])

    if (
        len(segments) > 1
        and segments[0][2] != 0
        and segments[0][2] == segments[-1][2]
    ):
        first = segments.pop(0)
        last = segments.pop(-1)
        segments.insert(0, [last[0], first[1] + n, first[2]])

    out = []
    for start, end, sign in segments:
        if sign == 0:
            continue
        idx = np.array([i % n for i in range(start, end + 1)], dtype=int)
        values = acceleration[idx]
        delta = float(np.sum(values) * step_rad)
        strength = abs(delta)
        if strength <= 1e-12:
            continue

        weights = np.abs(values)
        angles = np.radians(theta_deg[idx])
        sx = float(np.sum(weights * np.cos(angles)))
        sy = float(np.sum(weights * np.sin(angles)))
        center = math.degrees(math.atan2(sy, sx)) % 360.0

        out.append({
            "sign": int(sign),
            "start_deg": float(theta_deg[start % n] % 360.0),
            "end_deg": float(theta_deg[end % n] % 360.0),
            "center_deg": center,
            "span_deg": float(len(idx) * step_deg),
            "delta_velocity": delta,
            "strength": strength,
            "peak_abs_acceleration": float(np.max(weights)),
        })

    out.sort(key=lambda x: x["strength"], reverse=True)
    return out


def principal_axis(E: np.ndarray) -> tuple[np.ndarray, float]:
    """Principal motion axis and fraction of E-position variance it captures."""
    X = E - np.mean(E, axis=0)
    cov = X.T @ X
    values, vectors = np.linalg.eigh(cov)
    order = np.argsort(values)
    major = float(values[order[-1]])
    total = float(np.sum(values))
    axis = vectors[:, order[-1]]
    axis = axis / np.linalg.norm(axis)

    # Deterministic sign only; cadence matching may still use global sign flip.
    if axis[0] < 0.0 or (abs(axis[0]) < 1e-12 and axis[1] < 0.0):
        axis = -axis

    fraction = major / max(total, 1e-30)
    return axis, fraction


def block_bounds_unwrapped(block: dict, center_unwrapped: float) -> tuple[float, float]:
    c = float(block["center_deg"])
    ds = (float(block["start_deg"]) - c + 180.0) % 360.0 - 180.0
    de = (float(block["end_deg"]) - c + 180.0) % 360.0 - 180.0
    if ds > 0.0:
        ds -= 360.0
    if de < 0.0:
        de += 360.0
    return center_unwrapped + ds, center_unwrapped + de


def phase_linearity(
    theta_deg: np.ndarray,
    velocity: np.ndarray,
    assigned_blocks: list[dict],
    *,
    minimum_phase_span_deg: float,
) -> dict:
    """Fit v(theta)=a*theta+b in each calm interval between selected blocks."""
    if len(assigned_blocks) != 3:
        raise ValueError("phase_linearity expects exactly three blocks")

    # Keep the target-order assignment, unwrap its centers in positive-angle order.
    centers = [float(x["center_deg"]) % 360.0 for x in assigned_blocks]
    unwrapped_centers = [centers[0]]
    for raw in centers[1:]:
        prev_raw = centers[len(unwrapped_centers) - 1]
        delta = (raw - prev_raw) % 360.0
        unwrapped_centers.append(unwrapped_centers[-1] + delta)

    intervals = []
    block_ranges = [
        block_bounds_unwrapped(block, center)
        for block, center in zip(assigned_blocks, unwrapped_centers)
    ]

    scale = max(
        float(np.ptp(velocity)),
        0.25 * float(np.mean(np.abs(velocity))),
        1e-8,
    )

    theta0 = np.asarray(theta_deg, dtype=float)
    v0 = np.asarray(velocity, dtype=float)
    theta_ext = np.concatenate((theta0 - 360.0, theta0, theta0 + 360.0, theta0 + 720.0))
    v_ext = np.tile(v0, 4)

    weighted_sq = 0.0
    weight_sum = 0.0

    for i in range(3):
        current_end = block_ranges[i][1]
        if i < 2:
            next_start = block_ranges[i + 1][0]
        else:
            next_start = block_ranges[0][0] + 360.0

        span = float(next_start - current_end)
        if span < minimum_phase_span_deg:
            intervals.append({
                "span_deg": span,
                "used": False,
                "slope_per_rad": None,
                "normalized_rms": None,
            })
            continue

        mask = (theta_ext >= current_end) & (theta_ext <= next_start)
        xdeg = theta_ext[mask]
        y = v_ext[mask]
        if len(y) < 4:
            intervals.append({
                "span_deg": span,
                "used": False,
                "slope_per_rad": None,
                "normalized_rms": None,
            })
            continue

        x = np.radians(xdeg - xdeg[0])
        A = np.column_stack((x, np.ones_like(x)))
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        fitted = A @ coef
        rms = float(np.sqrt(np.mean((y - fitted) ** 2)))
        nrms = rms / scale

        intervals.append({
            "span_deg": span,
            "used": True,
            "slope_per_rad": float(coef[0]),
            "normalized_rms": nrms,
        })
        weighted_sq += span * nrms * nrms
        weight_sum += span

    overall = math.sqrt(weighted_sq / weight_sum) if weight_sum > 0.0 else 1.0
    return {
        "normalized_rms": float(overall),
        "intervals": intervals,
        "used_span_deg": float(weight_sum),
    }


def load_target(path: Path, args) -> dict:
    raw = np.genfromtxt(path, delimiter=",", names=True)
    if raw["theta_rad"][-1] >= 2.0 * math.pi - 1e-9:
        raw = raw[:-1]

    theta = np.asarray(raw["theta_rad"], dtype=float)
    theta_deg = np.mod(np.degrees(theta), 360.0)
    dq = np.asarray(raw["large_dq_dtheta_per_rad"], dtype=float)

    band = motor3952._detect_noise_bands(path, args.noise_fronts)["large"]
    noise_mask = motor3952._theta_band_mask(theta, band)

    # Preserve the broad HP transition but remove the artificial six-front ripple.
    dq_bridge = bridge_mask(dq, noise_mask)
    step_deg = 360.0 / len(theta)
    step_rad = 2.0 * math.pi / len(theta)
    velocity = smooth_periodic(dq_bridge, step_deg, args.target_smoothing_deg)
    acceleration = derivative_periodic(velocity, step_rad)

    lobes = extract_lobes(theta_deg, acceleration, args.lobe_epsilon_fraction)
    if len(lobes) < 3:
        raise RuntimeError("target lobe extraction found fewer than three blocks")

    blocks = sorted(lobes[:3], key=lambda x: x["center_deg"])
    linearity = phase_linearity(
        theta_deg,
        velocity,
        blocks,
        minimum_phase_span_deg=args.minimum_phase_span_deg,
    )

    return {
        "theta": theta,
        "theta_deg": theta_deg,
        "q": np.asarray(raw["large_fraction_0_1"], dtype=float),
        "velocity_raw": dq,
        "velocity_bridged": dq_bridge,
        "velocity_smooth": velocity,
        "acceleration": acceleration,
        "noise_band": band,
        "all_lobes": lobes,
        "blocks": blocks,
        "phase_linearity": linearity,
    }


def primary_features(
    v: np.ndarray,
    branch: int,
    theta: np.ndarray,
    args,
) -> dict:
    E, Ed, sine, pdata = stage.primary(theta, v, branch)
    theta_deg = np.mod(np.degrees(theta), 360.0)
    step_deg = 360.0 / len(theta)
    step_rad = 2.0 * math.pi / len(theta)

    axis, axis_fraction = principal_axis(E)
    projected_velocity = Ed @ axis
    projected_velocity = smooth_periodic(
        projected_velocity, step_deg, args.primary_smoothing_deg
    )
    projected_acceleration = derivative_periodic(projected_velocity, step_rad)

    lobes = extract_lobes(
        theta_deg,
        projected_acceleration,
        args.lobe_epsilon_fraction,
    )

    speed = np.linalg.norm(Ed, axis=1)
    center = np.mean(E, axis=0)
    X = E - center
    normal = np.array((-axis[1], axis[0]))
    along = X @ axis
    lateral = X @ normal

    return {
        "theta_deg": theta_deg,
        "E": E,
        "Ed": Ed,
        "speed": speed,
        "axis": axis,
        "axis_angle_deg": math.degrees(math.atan2(axis[1], axis[0])),
        "axis_position_variance_fraction": axis_fraction,
        "projected_velocity": projected_velocity,
        "projected_acceleration": projected_acceleration,
        "lobes": lobes,
        "minimum_primary_transmission_sine": float(sine),
        "primary": pdata,
        "E_span_over_crank": max(float(np.ptp(E[:, 0])), float(np.ptp(E[:, 1]))),
        "E_path_length_over_crank": float(
            np.sum(np.linalg.norm(np.roll(E, -1, axis=0) - E, axis=1))
        ),
        "E_principal_span_over_crank": float(np.ptp(along)),
        "E_lateral_span_over_crank": float(np.ptp(lateral)),
    }


def cadence_score(features: dict, target: dict, args) -> tuple[float, dict | None]:
    sine = features["minimum_primary_transmission_sine"]
    if sine < args.primary_sine_floor:
        return (
            20.0
            + 100.0
            * (args.primary_sine_floor - sine)
            / max(args.primary_sine_floor, 1e-9),
            None,
        )

    if features["E_span_over_crank"] < args.E_span_floor:
        return (
            20.0
            + 100.0
            * (args.E_span_floor - features["E_span_over_crank"])
            / max(args.E_span_floor, 1e-9),
            None,
        )

    lobes = features["lobes"]
    if len(lobes) < 3:
        return 50.0, None

    strongest = lobes[0]["strength"]
    pool = [x for x in lobes if x["strength"] >= 0.02 * strongest][:10]
    if len(pool) < 3:
        pool = lobes[:3]

    target_blocks = target["blocks"]
    target_strength = np.array([x["strength"] for x in target_blocks], dtype=float)
    target_strength /= np.sum(target_strength)

    best = None

    # Principal-axis sign is arbitrary, hence one global sign inversion is free.
    for sign_multiplier in (+1, -1):
        cost = np.empty((3, len(pool)), dtype=float)

        for i, target_block in enumerate(target_blocks):
            for j, candidate_block in enumerate(pool):
                center_error = cyclic_distance_deg(
                    candidate_block["center_deg"],
                    target_block["center_deg"],
                ) / args.center_scale_deg

                width_error = math.log(
                    max(candidate_block["span_deg"], 1.0)
                    / max(target_block["span_deg"], 1.0)
                )

                sign_error = (
                    0.0
                    if candidate_block["sign"]
                    == sign_multiplier * target_block["sign"]
                    else 1.0
                )

                cost[i, j] = (
                    center_error**2
                    + args.width_weight * width_error**2
                    + args.sign_weight * sign_error
                )

        rows, cols = linear_sum_assignment(cost)
        assigned = [pool[int(j)] for j in cols]

        phase = phase_linearity(
            features["theta_deg"],
            features["projected_velocity"],
            assigned,
            minimum_phase_span_deg=args.minimum_phase_span_deg,
        )

        center_errors = np.array([
            cyclic_distance_deg(
                assigned[i]["center_deg"],
                target_blocks[i]["center_deg"],
            )
            for i in range(3)
        ])

        width_errors = np.array([
            math.log(
                max(assigned[i]["span_deg"], 1.0)
                / max(target_blocks[i]["span_deg"], 1.0)
            )
            for i in range(3)
        ])

        candidate_strength = np.array(
            [assigned[i]["strength"] for i in range(3)],
            dtype=float,
        )
        candidate_strength /= np.sum(candidate_strength)
        strength_error = float(
            np.sqrt(np.mean((candidate_strength - target_strength) ** 2))
        )

        matched_strength = float(sum(x["strength"] for x in assigned))
        total_strength = float(sum(x["strength"] for x in lobes))
        extra_activity = max(
            0.0,
            (total_strength - matched_strength)
            / max(matched_strength, 1e-12),
        )

        timing_part = math.sqrt(float(np.mean(cost[rows, cols])))

        score = (
            timing_part
            + args.strength_weight * strength_error
            + args.extra_activity_weight * extra_activity
            + args.phase_linearity_weight * phase["normalized_rms"]
        )

        meta = {
            "score": float(score),
            "timing_part": timing_part,
            "sign_multiplier": int(sign_multiplier),
            "assigned_blocks": assigned,
            "center_rms_deg": float(np.sqrt(np.mean(center_errors**2))),
            "center_max_abs_deg": float(np.max(center_errors)),
            "width_log_rms": float(np.sqrt(np.mean(width_errors**2))),
            "relative_strength_rms": strength_error,
            "extra_activity_ratio": extra_activity,
            "phase_linearity": phase,
            "all_lobe_count": len(lobes),
            "significant_lobe_count": len(pool),
        }

        if best is None or score < best[0]:
            best = (score, meta)

    assert best is not None
    return best


def primary_closure_penalty(v: np.ndarray) -> float:
    g, c, r = map(float, v[:3])
    outer_margin = c + r - (g + 1.0)
    inner_margin = abs(g - 1.0) - abs(c - r)
    scale = max(g + c + r + 1.0, 1.0)
    violation = max(0.0, -outer_margin) + max(0.0, -inner_margin)
    return violation / scale


def make_record(
    v: np.ndarray,
    branch: int,
    target: dict,
    theta: np.ndarray,
    args,
) -> dict | None:
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
        "minimum_primary_transmission_sine": features[
            "minimum_primary_transmission_sine"
        ],
        "projection_axis_angle_deg": features["axis_angle_deg"],
        "axis_position_variance_fraction": features[
            "axis_position_variance_fraction"
        ],
        "E_span_over_crank": features["E_span_over_crank"],
        "E_principal_span_over_crank": features["E_principal_span_over_crank"],
        "E_lateral_span_over_crank": features["E_lateral_span_over_crank"],
        "E_path_length_over_crank": features["E_path_length_over_crank"],
        "lobes": features["lobes"],
    }


def geometry_distance(a: dict, b: dict) -> float:
    if a["branch"] != b["branch"]:
        return math.inf
    xa = np.array([a["parameters"][k] for k in NAMES], dtype=float)
    xb = np.array([b["parameters"][k] for k in NAMES], dtype=float)

    values = []
    for i, (lo, hi) in enumerate(BOUNDS):
        if i == 5:
            delta = abs((xa[i] - xb[i] + math.pi) % (2.0 * math.pi) - math.pi)
            values.append(delta / math.pi)
        else:
            values.append(abs(xa[i] - xb[i]) / (hi - lo))
    return float(np.sqrt(np.mean(np.square(values))))


def diverse(records: list[dict], count: int, minimum_distance: float) -> list[dict]:
    kept = []
    for item in sorted(records, key=lambda x: x["score"]):
        if all(
            geometry_distance(item, other) >= minimum_distance
            for other in kept
        ):
            kept.append(item)
            if len(kept) >= count:
                break
    for rank, item in enumerate(kept):
        item["library_rank"] = rank
    return kept


def known_sixbar_reference(path: Path, target: dict, theta: np.ndarray, args) -> dict | None:
    if not path.exists():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        best = data.get("best")
        if not best:
            return None
        p = best["primary"]
        v = np.array([
            float(p["ground"]),
            float(p["coupler"]),
            float(p["rocker"]),
            float(p["E_along"]),
            float(p["E_normal"]),
            float(p["phase"]),
        ])
        branch = int(p["assembly_branch"])
        record = make_record(v, branch, target, theta, args)
        if record is None:
            return {
                "available": True,
                "scored": False,
                "reason": "reference primary failed V2 feasibility/score evaluation",
            }
        return {
            "available": True,
            "scored": True,
            "known_sixbar_position_rms": best.get("position_rms"),
            "primary_record": record,
        }
    except Exception as exc:
        return {
            "available": True,
            "scored": False,
            "reason": f"{type(exc).__name__}: {exc}",
        }


def save(path: Path, obj: dict) -> None:
    def convert(x):
        if isinstance(x, Path):
            return str(x)
        if isinstance(x, np.generic):
            return x.item()
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, dict):
            return {k: convert(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [convert(v) for v in x]
        return x

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(convert(obj), indent=2) + "\n", encoding="utf-8")


def write_trace(path: Path, best: dict, target: dict, theta: np.ndarray, args) -> None:
    v = np.array([best["parameters"][k] for k in NAMES], dtype=float)
    features = primary_features(v, int(best["branch"]), theta, args)

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        fields = [
            "theta_deg",
            "E_x",
            "E_y",
            "E_speed",
            "E_projected_velocity",
            "E_projected_acceleration",
            "target_velocity",
            "target_acceleration",
        ]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for i in range(len(theta)):
            writer.writerow({
                "theta_deg": float(np.degrees(theta[i]) % 360.0),
                "E_x": float(features["E"][i, 0]),
                "E_y": float(features["E"][i, 1]),
                "E_speed": float(features["speed"][i]),
                "E_projected_velocity": float(features["projected_velocity"][i]),
                "E_projected_acceleration": float(
                    features["projected_acceleration"][i]
                ),
                "target_velocity": float(target["velocity_smooth"][i]),
                "target_acceleration": float(target["acceleration"][i]),
            })


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=Path, default=TARGET)
    ap.add_argument("--output", type=Path, default=OUTPUT)
    ap.add_argument("--trace-output", type=Path, default=TRACE)
    ap.add_argument("--known-sixbar", type=Path, default=KNOWN_SIXBAR)

    ap.add_argument("--islands", type=int, default=32)
    ap.add_argument("--generations", type=int, default=240)
    ap.add_argument("--population-size", type=int, default=16)
    ap.add_argument("--seed", type=int, default=3952500)

    ap.add_argument("--coarse-step-deg", type=float, default=1.0)
    ap.add_argument("--dense-step-deg", type=float, default=0.25)
    ap.add_argument("--target-smoothing-deg", type=float, default=12.0)
    ap.add_argument("--primary-smoothing-deg", type=float, default=8.0)
    ap.add_argument("--lobe-epsilon-fraction", type=float, default=0.01)
    ap.add_argument("--noise-fronts", type=int, default=6)

    # Timing is now deliberately soft; phase simplicity matters strongly.
    ap.add_argument("--center-scale-deg", type=float, default=45.0)
    ap.add_argument("--width-weight", type=float, default=0.10)
    ap.add_argument("--strength-weight", type=float, default=0.05)
    ap.add_argument("--sign-weight", type=float, default=0.05)
    ap.add_argument("--extra-activity-weight", type=float, default=0.20)
    ap.add_argument("--phase-linearity-weight", type=float, default=1.50)
    ap.add_argument("--minimum-phase-span-deg", type=float, default=8.0)

    ap.add_argument("--primary-sine-floor", type=float, default=0.30)
    ap.add_argument("--E-span-floor", type=float, default=0.75)

    ap.add_argument("--library-size", type=int, default=200)
    ap.add_argument("--diversity-distance", type=float, default=0.08)
    ap.add_argument("--population-candidates-per-island", type=int, default=24)
    args = ap.parse_args()

    target = load_target(args.target, args)

    print("TARGET LARGE V2 BLOCKS (signed piston velocity):", flush=True)
    for i, block in enumerate(target["blocks"], 1):
        print(
            f"  {i}: {block['start_deg']:.2f}->{block['end_deg']:.2f} deg "
            f"center={block['center_deg']:.2f} "
            f"span={block['span_deg']:.2f} "
            f"sign={block['sign']:+d} "
            f"|Dv|={block['strength']:.6f}",
            flush=True,
        )
    print(
        "target inter-block linearity RMS = "
        f"{target['phase_linearity']['normalized_rms']:.6f}",
        flush=True,
    )

    nc = int(round(360.0 / args.coarse_step_deg))
    nd = int(round(360.0 / args.dense_step_deg))
    coarse_theta = np.linspace(0.0, 2.0 * math.pi, nc, endpoint=False)
    dense_theta = np.linspace(0.0, 2.0 * math.pi, nd, endpoint=False)

    reference = known_sixbar_reference(
        args.known_sixbar,
        target,
        dense_theta,
        args,
    )
    if reference and reference.get("scored"):
        rr = reference["primary_record"]
        print(
            "KNOWN 2.08%-CLASS SIX-BAR PRIMARY REFERENCE: "
            f"V2 score={rr['score']:.5f} "
            f"centerRMS={rr['cadence_match']['center_rms_deg']:.2f}deg "
            f"phaseLinear={rr['cadence_match']['phase_linearity']['normalized_rms']:.4f} "
            f"sine={rr['minimum_primary_transmission_sine']:.3f}",
            flush=True,
        )

    wall0 = time.time()
    cpu0 = time.process_time()

    report = {
        "description": (
            "Candidate 3952 LARGE primary cadence V2: signed principal-axis "
            "velocity, soft timing match, and near-linear velocity between "
            "dominant acceleration blocks."
        ),
        "length_unit": "primary_crank_radius",
        "target": str(args.target),
        "target_noise_band": target["noise_band"],
        "target_blocks": target["blocks"],
        "target_all_lobes": target["all_lobes"],
        "target_phase_linearity": target["phase_linearity"],
        "known_sixbar_reference": reference,
        "bounds": dict(zip(NAMES, BOUNDS)),
        "settings": vars(args),
        "islands": [],
        "candidates": [],
        "best": None,
    }

    raw_records = []

    for island in range(args.islands):
        branch = 1 if island % 2 == 0 else -1
        seed = args.seed + island
        archive = []

        def objective(v):
            v = np.asarray(v, dtype=float)
            closure = primary_closure_penalty(v)
            if closure > 0.0:
                return 20.0 + 100.0 * closure

            try:
                features = primary_features(v, branch, coarse_theta, args)
                score, match = cadence_score(features, target, args)
            except (
                ValueError,
                FloatingPointError,
                OverflowError,
                np.linalg.LinAlgError,
            ):
                return 1e3

            if match is None:
                return float(score)

            archive.append((float(score), v.copy()))
            return float(score)

        t0 = time.time()
        p0 = time.process_time()

        opt = differential_evolution(
            objective,
            BOUNDS,
            seed=seed,
            popsize=args.population_size,
            maxiter=args.generations,
            polish=True,
            tol=1e-8,
            updating="immediate",
            workers=1,
            init="latinhypercube",
        )

        vectors = [np.asarray(opt.x, dtype=float)]
        if getattr(opt, "population", None) is not None:
            order = np.argsort(np.asarray(opt.population_energies))
            vectors += [
                np.asarray(opt.population[int(j)], dtype=float)
                for j in order[: args.population_candidates_per_island]
            ]
        vectors += [
            x for _, x in sorted(archive, key=lambda z: z[0])[
                : args.population_candidates_per_island
            ]
        ]

        dense = []
        seen = set()
        for v in vectors:
            key = tuple(round(float(x), 10) for x in v)
            if key in seen:
                continue
            seen.add(key)
            record = make_record(v, branch, target, dense_theta, args)
            if record is not None:
                dense.append(record)
                raw_records.append(record)

        dense.sort(key=lambda x: x["score"])
        best = dense[0] if dense else None

        report["islands"].append({
            "index": island,
            "branch": branch,
            "seed": seed,
            "optimizer_success": bool(opt.success),
            "optimizer_message": str(opt.message),
            "function_evaluations": int(opt.nfev),
            "wall_seconds": time.time() - t0,
            "cpu_seconds": time.process_time() - p0,
            "best": best,
        })

        library = diverse(
            raw_records,
            args.library_size,
            args.diversity_distance,
        )
        report["candidates"] = library
        report["best"] = library[0] if library else None
        report["elapsed_wall_hours"] = (time.time() - wall0) / 3600.0
        report["elapsed_cpu_hours"] = (time.process_time() - cpu0) / 3600.0
        save(args.output, report)

        if best is None:
            print(
                f"island {island:03d} branch={branch:+d}: no dense candidate",
                flush=True,
            )
        else:
            m = best["cadence_match"]
            print(
                f"island {island:03d} branch={branch:+d}: "
                f"score={best['score']:.5f} "
                f"centerRMS={m['center_rms_deg']:.2f}deg "
                f"phaseLinear={m['phase_linearity']['normalized_rms']:.4f} "
                f"extra={m['extra_activity_ratio']:.3f} "
                f"sine={best['minimum_primary_transmission_sine']:.3f} "
                f"axisVar={best['axis_position_variance_fraction']:.3f}",
                flush=True,
            )

    library = diverse(
        raw_records,
        args.library_size,
        args.diversity_distance,
    )
    report["candidates"] = library
    report["best"] = library[0] if library else None
    report["completed_islands"] = args.islands
    report["elapsed_wall_hours"] = (time.time() - wall0) / 3600.0
    report["elapsed_cpu_hours"] = (time.process_time() - cpu0) / 3600.0
    save(args.output, report)

    print("\nPRIMARY CADENCE V2 SUMMARY", flush=True)
    print(
        f"wall={report['elapsed_wall_hours']:.3f}h "
        f"cpu={report['elapsed_cpu_hours']:.3f}h "
        f"diverse={len(library)}",
        flush=True,
    )

    if not library:
        print("No feasible V2 primary found.", flush=True)
        return

    best = library[0]
    m = best["cadence_match"]
    print(
        f"best score={best['score']:.6f}\n"
        f"branch={best['branch']:+d}\n"
        f"center RMS={m['center_rms_deg']:.3f} deg\n"
        f"phase linearity={m['phase_linearity']['normalized_rms']:.6f}\n"
        f"extra activity={m['extra_activity_ratio']:.4f}\n"
        f"primary sine={best['minimum_primary_transmission_sine']:.4f}\n"
        f"projection axis={best['projection_axis_angle_deg']:.3f} deg\n"
        f"axis position variance fraction={best['axis_position_variance_fraction']:.4f}\n"
        f"AD={best['primary']['ground']:.6f} "
        f"BC={best['primary']['coupler']:.6f} "
        f"CD={best['primary']['rocker']:.6f}",
        flush=True,
    )

    print("matched blocks:", flush=True)
    for t, c in zip(target["blocks"], m["assigned_blocks"]):
        print(
            f"  target {t['center_deg']:7.2f}deg <- "
            f"Eproj {c['center_deg']:7.2f}deg "
            f"span={c['span_deg']:6.2f}deg sign={c['sign']:+d}",
            flush=True,
        )

    print("inter-block projected-velocity fits:", flush=True)
    for i, phase in enumerate(m["phase_linearity"]["intervals"], 1):
        if not phase["used"]:
            print(
                f"  phase {i}: span={phase['span_deg']:.2f}deg (too short)",
                flush=True,
            )
        else:
            print(
                f"  phase {i}: span={phase['span_deg']:.2f}deg "
                f"slope={phase['slope_per_rad']:.6f} "
                f"nRMS={phase['normalized_rms']:.5f}",
                flush=True,
            )

    write_trace(args.trace_output, best, target, dense_theta, args)
    print(f"Wrote {args.output}", flush=True)
    print(f"Wrote {args.trace_output}", flush=True)


if __name__ == "__main__":
    main()
