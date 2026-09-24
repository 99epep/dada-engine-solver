#!/usr/bin/env python3
"""Primary cadence search V3 for candidate 3952, LARGE cylinder.

V3 abandons "three acceleration peaks" as the primary representation.
Instead it automatically approximates the target piston velocity by THREE
cyclic linear regimes:

    v(theta) ~= a_i * theta + b_i

A linear velocity regime means constant acceleration; a near-zero slope means
nearly constant velocity.  Candidate E motion is projected onto its principal
axis, then judged on the same three temporal regimes.

The target phase boundaries are discovered automatically from the bridged and
smoothed LARGE target velocity.  For candidate 3952 they are expected near:

    31 deg / 181 deg / 349 deg
    spans ~150 deg / 168 deg / 42 deg

Candidate scoring:
  - coarse regime-shape agreement after the best affine velocity rescaling;
  - projected-velocity linearity inside each regime;
  - soft preference for primary transmission sine >= 0.35;
  - weak preference for a reasonably directional E trajectory;
  - hard minimum primary sine 0.30 and E amplitude floor.

The exact transition neighborhoods are excluded from the regime fits so a
physical four-bar may round the transitions naturally.

V2 candidates are NOT discarded.  V3 rescored candidates from the V2 library
and V2 island winners, keeps a diverse best subset, and injects up to 32 of
them as x0 seeds into broad-bound V3 islands.  The remaining islands are fresh
Latin-hypercube searches.

Typical run:
    PYTHONPATH=src python3 examples/search_large_primary_cadence_3952_v3.py

Outputs:
    outputs/large_primary_cadence_3952_v3.json
    outputs/large_primary_cadence_3952_v3_best_trace.csv
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
TARGET = (
    ROOT
    / "outputs"
    / "motor_hybrid_c2_15p_260k"
    / "candidate_3952_motion_target.csv"
)
V2_LIBRARY = ROOT / "outputs" / "large_primary_cadence_3952_v2.json"
KNOWN_SIXBAR = ROOT / "outputs" / "large_sixbar_3952.json"
OUTPUT = ROOT / "outputs" / "large_primary_cadence_3952_v3.json"
TRACE = ROOT / "outputs" / "large_primary_cadence_3952_v3_best_trace.csv"

NAMES = [
    "primary_ground",
    "primary_coupler",
    "primary_rocker",
    "primary_E_along",
    "primary_E_normal",
    "primary_phase",
]

# Primary crank AB = 1.
BOUNDS = [
    (0.50, 12.0),               # AD
    (0.50, 12.0),               # BC
    (0.35, 6.0),                # CD
    (-12.0, 18.0),              # E along BC
    (-15.0, 15.0),              # E normal to BC
    (-math.pi, math.pi),        # cycle phase
]


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
    """Linearly bridge one contiguous cyclic True interval."""
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
    """Principal position axis and fraction of position variance it captures."""
    X = np.asarray(E, dtype=float) - np.mean(E, axis=0)
    cov = X.T @ X
    values, vectors = np.linalg.eigh(cov)
    i = int(np.argmax(values))
    axis = vectors[:, i]
    axis = axis / np.linalg.norm(axis)

    # Deterministic orientation.  Global velocity sign remains free through the
    # affine comparison with the target.
    if axis[0] < 0.0 or (abs(axis[0]) < 1e-12 and axis[1] < 0.0):
        axis = -axis

    major = float(values[i])
    total = float(np.sum(values))
    return axis, major / max(total, 1e-30)


def line_fit(x_rad: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Return slope, intercept, RMS."""
    x_rad = np.asarray(x_rad, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(y) < 3:
        raise ValueError("too few points for line fit")
    A = np.column_stack((x_rad, np.ones_like(x_rad)))
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    fitted = A @ coef
    rms = float(np.sqrt(np.mean((y - fitted) ** 2)))
    return float(coef[0]), float(coef[1]), rms


def interval_sse_from_prefix(
    a,
    b,
    px,
    py,
    pxx,
    pxy,
    pyy,
):
    """Vectorized line-fit SSE on half-open intervals [a,b)."""
    a = np.asarray(a)
    b = np.asarray(b)
    n = b - a

    sx = px[b] - px[a]
    sy = py[b] - py[a]
    sxx = pxx[b] - pxx[a]
    sxy = pxy[b] - pxy[a]
    syy = pyy[b] - pyy[a]

    den = n * sxx - sx * sx
    slope = np.divide(
        n * sxy - sx * sy,
        den,
        out=np.zeros_like(np.asarray(den, dtype=float)),
        where=np.abs(den) > 1e-15,
    )
    intercept = (sy - slope * sx) / n

    sse = (
        syy
        + slope * slope * sxx
        + n * intercept * intercept
        - 2.0 * slope * sxy
        - 2.0 * intercept * sy
        + 2.0 * slope * intercept * sx
    )
    return np.maximum(sse, 0.0)


def discover_target_regimes(
    theta_deg: np.ndarray,
    velocity: np.ndarray,
    *,
    minimum_span_deg: int,
) -> dict:
    """Best cyclic 3-line partition of the target velocity on a 1-degree grid."""
    theta_deg = np.asarray(theta_deg, dtype=float)
    velocity = np.asarray(velocity, dtype=float)

    sample_deg = np.arange(360, dtype=float)
    sample_v = np.interp(
        sample_deg,
        theta_deg,
        velocity,
        period=360.0,
    )
    doubled = np.concatenate((sample_v, sample_v))
    x = np.arange(720, dtype=float)

    def prefix(z):
        return np.concatenate(([0.0], np.cumsum(z, dtype=float)))

    px = prefix(x)
    py = prefix(doubled)
    pxx = prefix(x * x)
    pxy = prefix(x * doubled)
    pyy = prefix(doubled * doubled)

    best = None
    m = int(minimum_span_deg)

    # About 100k Python iterations, with the expensive b2 search vectorized.
    for start in range(360):
        end = start + 360
        for b1 in range(start + m, end - 2 * m + 1):
            b2 = np.arange(b1 + m, end - m + 1, dtype=int)
            if len(b2) == 0:
                continue

            sse1 = float(
                interval_sse_from_prefix(
                    np.array(start),
                    np.array(b1),
                    px, py, pxx, pxy, pyy,
                )
            )
            sse23 = (
                interval_sse_from_prefix(
                    np.full_like(b2, b1),
                    b2,
                    px, py, pxx, pxy, pyy,
                )
                + interval_sse_from_prefix(
                    b2,
                    np.full_like(b2, end),
                    px, py, pxx, pxy, pyy,
                )
            )
            j = int(np.argmin(sse23))
            total = sse1 + float(sse23[j])
            if best is None or total < best["sse"]:
                best = {
                    "sse": total,
                    "start": start,
                    "b1": b1,
                    "b2": int(b2[j]),
                    "end": end,
                }

    if best is None:
        raise RuntimeError("target 3-regime segmentation failed")

    boundaries = [
        float(best["start"] % 360),
        float(best["b1"] % 360),
        float(best["b2"] % 360),
    ]
    spans = [
        float(best["b1"] - best["start"]),
        float(best["b2"] - best["b1"]),
        float(best["end"] - best["b2"]),
    ]

    target_fit = fit_fixed_regimes(
        theta_deg,
        velocity,
        boundaries,
        transition_half_width_deg=0.0,
    )

    return {
        "boundaries_deg": boundaries,
        "spans_deg": spans,
        "full_fit_rms": math.sqrt(best["sse"] / 360.0),
        "full_fit_normalized_rms": (
            math.sqrt(best["sse"] / 360.0)
            / max(float(np.std(sample_v)), 1e-12)
        ),
        "slopes_per_rad": [
            x["slope_per_rad"] for x in target_fit["regimes"]
        ],
    }


def phase_coordinates(
    theta_deg: np.ndarray,
    boundaries_deg: list[float],
) -> tuple[np.ndarray, list[float]]:
    """Unwrap theta from first boundary; return u in [0,360) and cuts."""
    start = float(boundaries_deg[0]) % 360.0
    u = (np.asarray(theta_deg, dtype=float) - start) % 360.0

    d1 = (float(boundaries_deg[1]) - start) % 360.0
    d2 = (float(boundaries_deg[2]) - start) % 360.0
    cuts = [0.0, d1, d2, 360.0]

    if not (0.0 < d1 < d2 < 360.0):
        raise ValueError("regime boundaries are not in positive cyclic order")
    return u, cuts


def fit_fixed_regimes(
    theta_deg: np.ndarray,
    velocity: np.ndarray,
    boundaries_deg: list[float],
    *,
    transition_half_width_deg: float,
) -> dict:
    """Fit one velocity line inside the core of each fixed temporal regime."""
    theta_deg = np.asarray(theta_deg, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    if theta_deg.shape != velocity.shape:
        raise ValueError("theta/velocity shape mismatch")

    u, cuts = phase_coordinates(theta_deg, boundaries_deg)

    model = np.full_like(velocity, np.nan, dtype=float)
    core_mask = np.zeros_like(velocity, dtype=bool)
    regimes = []
    weighted_sq = 0.0
    weight_sum = 0.0

    scale = max(float(np.std(velocity)), 1e-12)

    for i in range(3):
        lo = cuts[i]
        hi = cuts[i + 1]
        span = hi - lo

        margin = min(
            float(transition_half_width_deg),
            max(0.0, 0.20 * span),
        )
        core_lo = lo + margin
        core_hi = hi - margin

        mask = (u >= core_lo) & (u < core_hi)
        if np.count_nonzero(mask) < 5:
            raise ValueError("regime core too short")

        x = np.radians(u[mask])
        y = velocity[mask]
        slope, intercept, rms = line_fit(x, y)
        fitted = slope * x + intercept

        model[mask] = fitted
        core_mask |= mask
        core_span = core_hi - core_lo
        nrms = rms / scale

        regimes.append({
            "index": i,
            "start_deg": float(
                (boundaries_deg[0] + lo) % 360.0
            ),
            "end_deg": float(
                (boundaries_deg[0] + hi) % 360.0
            ),
            "span_deg": float(span),
            "core_span_deg": float(core_span),
            "slope_per_rad": float(slope),
            "rms": float(rms),
            "normalized_rms": float(nrms),
        })

        weighted_sq += core_span * nrms * nrms
        weight_sum += core_span

    overall = math.sqrt(weighted_sq / max(weight_sum, 1e-12))
    return {
        "model": model,
        "core_mask": core_mask,
        "regimes": regimes,
        "normalized_linearity_rms": float(overall),
        "used_span_deg": float(weight_sum),
    }


def affine_shape_error(
    candidate_model: np.ndarray,
    target_model: np.ndarray,
    mask: np.ndarray,
) -> dict:
    """Best target ~= scale*candidate + offset on common regime cores."""
    c = np.asarray(candidate_model, dtype=float)[mask]
    t = np.asarray(target_model, dtype=float)[mask]
    if len(c) < 5:
        raise ValueError("too few points for shape comparison")

    A = np.column_stack((c, np.ones_like(c)))
    coef, *_ = np.linalg.lstsq(A, t, rcond=None)
    fitted = A @ coef
    rms = float(np.sqrt(np.mean((t - fitted) ** 2)))
    normalized = rms / max(float(np.std(t)), 1e-12)
    return {
        "normalized_rms": float(normalized),
        "rms": rms,
        "velocity_scale": float(coef[0]),
        "velocity_offset": float(coef[1]),
    }


def load_target(path: Path, args) -> dict:
    raw = np.genfromtxt(path, delimiter=",", names=True)
    if raw["theta_rad"][-1] >= 2.0 * math.pi - 1e-9:
        raw = raw[:-1]

    theta = np.asarray(raw["theta_rad"], dtype=float)
    theta_deg = np.mod(np.degrees(theta), 360.0)
    dq = np.asarray(raw["large_dq_dtheta_per_rad"], dtype=float)

    band = motor3952._detect_noise_bands(
        path,
        args.noise_fronts,
    )["large"]
    noise_mask = motor3952._theta_band_mask(theta, band)

    dq_bridge = bridge_mask(dq, noise_mask)
    step_deg = 360.0 / len(theta)
    velocity = smooth_periodic(
        dq_bridge,
        step_deg,
        args.target_smoothing_deg,
    )

    segmentation = discover_target_regimes(
        theta_deg,
        velocity,
        minimum_span_deg=args.minimum_regime_span_deg,
    )

    target_regimes = fit_fixed_regimes(
        theta_deg,
        velocity,
        segmentation["boundaries_deg"],
        transition_half_width_deg=args.transition_half_width_deg,
    )

    return {
        "theta": theta,
        "theta_deg": theta_deg,
        "q": np.asarray(raw["large_fraction_0_1"], dtype=float),
        "velocity_raw": dq,
        "velocity_bridged": dq_bridge,
        "velocity_smooth": velocity,
        "noise_band": band,
        "segmentation": segmentation,
        "regime_fit": target_regimes,
    }


def target_model_on_grid(
    target: dict,
    theta_deg: np.ndarray,
    args,
) -> dict:
    """Refit target velocity on an arbitrary uniform candidate grid."""
    target_velocity = np.interp(
        np.asarray(theta_deg, dtype=float),
        target["theta_deg"],
        target["velocity_smooth"],
        period=360.0,
    )
    return fit_fixed_regimes(
        theta_deg,
        target_velocity,
        target["segmentation"]["boundaries_deg"],
        transition_half_width_deg=args.transition_half_width_deg,
    )


def primary_features(
    v: np.ndarray,
    branch: int,
    theta: np.ndarray,
    args,
) -> dict:
    E, Ed, sine, pdata = stage.primary(theta, v, branch)
    theta_deg = np.mod(np.degrees(theta), 360.0)
    step_deg = 360.0 / len(theta)

    axis, axis_fraction = principal_axis(E)
    projected_velocity = Ed @ axis
    projected_velocity = smooth_periodic(
        projected_velocity,
        step_deg,
        args.primary_smoothing_deg,
    )

    center = np.mean(E, axis=0)
    X = E - center
    normal = np.array((-axis[1], axis[0]))
    along = X @ axis
    lateral = X @ normal

    return {
        "theta_deg": theta_deg,
        "E": E,
        "Ed": Ed,
        "axis": axis,
        "axis_angle_deg": math.degrees(math.atan2(axis[1], axis[0])),
        "axis_position_variance_fraction": axis_fraction,
        "projected_velocity": projected_velocity,
        "minimum_primary_transmission_sine": float(sine),
        "primary": pdata,
        "E_span_over_crank": max(
            float(np.ptp(E[:, 0])),
            float(np.ptp(E[:, 1])),
        ),
        "E_principal_span_over_crank": float(np.ptp(along)),
        "E_lateral_span_over_crank": float(np.ptp(lateral)),
        "E_path_length_over_crank": float(
            np.sum(np.linalg.norm(np.roll(E, -1, axis=0) - E, axis=1))
        ),
    }


def cadence_score(
    features: dict,
    target: dict,
    target_fit_on_grid: dict,
    args,
) -> tuple[float, dict | None]:
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

    try:
        candidate_fit = fit_fixed_regimes(
            features["theta_deg"],
            features["projected_velocity"],
            target["segmentation"]["boundaries_deg"],
            transition_half_width_deg=args.transition_half_width_deg,
        )
    except (ValueError, np.linalg.LinAlgError):
        return 50.0, None

    common_mask = (
        candidate_fit["core_mask"]
        & target_fit_on_grid["core_mask"]
        & np.isfinite(candidate_fit["model"])
        & np.isfinite(target_fit_on_grid["model"])
    )

    shape = affine_shape_error(
        candidate_fit["model"],
        target_fit_on_grid["model"],
        common_mask,
    )

    linearity = candidate_fit["normalized_linearity_rms"]

    # Soft mechanical preference, not a hard rejection.  0.30 remains legal.
    if sine >= args.preferred_primary_sine:
        sine_penalty = 0.0
    else:
        sine_penalty = args.sine_preference_weight * (
            args.preferred_primary_sine - sine
        ) / max(
            args.preferred_primary_sine - args.primary_sine_floor,
            1e-9,
        )

    axis_fraction = features["axis_position_variance_fraction"]
    if axis_fraction >= args.preferred_axis_variance_fraction:
        axis_penalty = 0.0
    else:
        axis_penalty = args.axis_preference_weight * (
            args.preferred_axis_variance_fraction - axis_fraction
        ) / max(args.preferred_axis_variance_fraction, 1e-9)

    score = (
        args.shape_weight * shape["normalized_rms"]
        + args.linearity_weight * linearity
        + sine_penalty
        + axis_penalty
    )

    candidate_slopes = np.array([
        x["slope_per_rad"] for x in candidate_fit["regimes"]
    ])
    target_slopes = np.array([
        x["slope_per_rad"] for x in target_fit_on_grid["regimes"]
    ])

    # Apply the best affine velocity scale to candidate slopes so the printed
    # slope comparison is dimensionless with respect to mechanism size.
    scaled_candidate_slopes = shape["velocity_scale"] * candidate_slopes
    slope_rms = float(
        np.sqrt(np.mean((scaled_candidate_slopes - target_slopes) ** 2))
    )
    target_slope_scale = max(
        float(np.sqrt(np.mean(target_slopes**2))),
        1e-12,
    )

    meta = {
        "score": float(score),
        "shape_normalized_rms": shape["normalized_rms"],
        "linearity_normalized_rms": float(linearity),
        "sine_preference_penalty": float(sine_penalty),
        "axis_preference_penalty": float(axis_penalty),
        "velocity_scale_to_target": shape["velocity_scale"],
        "velocity_offset_to_target": shape["velocity_offset"],
        "candidate_regimes": candidate_fit["regimes"],
        "target_regimes_on_same_grid": target_fit_on_grid["regimes"],
        "scaled_candidate_slopes_per_rad": scaled_candidate_slopes.tolist(),
        "target_slopes_per_rad": target_slopes.tolist(),
        "normalized_slope_rms": slope_rms / target_slope_scale,
    }
    return float(score), meta


def primary_closure_penalty(v: np.ndarray) -> float:
    """Continuous full-revolution closure guidance."""
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
    target_fit_on_grid: dict,
    theta: np.ndarray,
    args,
) -> dict | None:
    try:
        features = primary_features(v, branch, theta, args)
        score, match = cadence_score(
            features,
            target,
            target_fit_on_grid,
            args,
        )
    except (
        ValueError,
        FloatingPointError,
        OverflowError,
        np.linalg.LinAlgError,
    ):
        return None

    if match is None or score >= 20.0:
        return None

    return {
        "score": float(score),
        "branch": int(branch),
        "parameters": {
            name: float(x) for name, x in zip(NAMES, v)
        },
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
        "E_principal_span_over_crank": features[
            "E_principal_span_over_crank"
        ],
        "E_lateral_span_over_crank": features[
            "E_lateral_span_over_crank"
        ],
        "E_path_length_over_crank": features[
            "E_path_length_over_crank"
        ],
    }


def geometry_distance(a: dict, b: dict) -> float:
    if int(a["branch"]) != int(b["branch"]):
        return math.inf

    xa = np.array([a["parameters"][k] for k in NAMES], dtype=float)
    xb = np.array([b["parameters"][k] for k in NAMES], dtype=float)
    d = []

    for i, (lo, hi) in enumerate(BOUNDS):
        if i == 5:
            z = abs(
                (xa[i] - xb[i] + math.pi) % (2.0 * math.pi) - math.pi
            ) / math.pi
        else:
            z = abs(xa[i] - xb[i]) / (hi - lo)
        d.append(z)

    return float(np.sqrt(np.mean(np.square(d))))


def diverse(
    records: list[dict],
    count: int,
    minimum_distance: float,
) -> list[dict]:
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


def vector_from_record(item: dict) -> tuple[np.ndarray, int]:
    if "parameters" in item and all(
        name in item["parameters"] for name in NAMES
    ):
        v = np.array(
            [float(item["parameters"][name]) for name in NAMES],
            dtype=float,
        )
        branch = int(item.get(
            "branch",
            item.get("primary", {}).get("assembly_branch", 0),
        ))
        if branch not in (-1, 1):
            raise ValueError("record has no valid primary branch")
        return v, branch

    p = item["primary"]
    return (
        np.array([
            float(p["ground"]),
            float(p["coupler"]),
            float(p["rocker"]),
            float(p["E_along"]),
            float(p["E_normal"]),
            float(p["phase"]),
        ], dtype=float),
        int(p["assembly_branch"]),
    )


def rescore_v2_seeds(
    path: Path,
    target: dict,
    target_fit_dense: dict,
    dense_theta: np.ndarray,
    args,
) -> list[dict]:
    """Rescore both V2 library entries and V2 island winners under V3."""
    if not path.exists():
        return []

    data = json.loads(path.read_text(encoding="utf-8"))
    source_items = []

    for item in data.get("candidates", []):
        source_items.append(("v2_library", item))

    for island in data.get("islands", []):
        item = island.get("best")
        if item is not None:
            source_items.append(
                (f"v2_island_{int(island.get('index', -1))}", item)
            )

    rescored = []
    seen = set()

    for source, item in source_items:
        try:
            v, branch = vector_from_record(item)
        except (KeyError, ValueError, TypeError):
            continue

        key = (
            branch,
            *[round(float(x), 10) for x in v],
        )
        if key in seen:
            continue
        seen.add(key)

        record = make_record(
            v,
            branch,
            target,
            target_fit_dense,
            dense_theta,
            args,
        )
        if record is None:
            continue

        record["seed_source"] = source
        rescored.append(record)

    return sorted(rescored, key=lambda x: x["score"])


def select_v2_seed_panel(
    rescored: list[dict],
    *,
    count: int,
    minimum_distance: float,
) -> list[dict]:
    if count <= 0:
        return []

    # Aim for a balanced seed panel so one V2 assembly branch cannot occupy
    # every seeded V3 island.
    per_branch = count // 2
    chosen = []

    for branch in (+1, -1):
        branch_items = [
            x for x in rescored if int(x["branch"]) == branch
        ]
        for item in branch_items:
            if all(
                geometry_distance(item, other) >= minimum_distance
                for other in chosen
                if int(other["branch"]) == branch
            ):
                chosen.append(item)
                if sum(
                    int(x["branch"]) == branch for x in chosen
                ) >= per_branch:
                    break

    # Fill unused slots globally if one branch has too few candidates.
    for item in rescored:
        if len(chosen) >= count:
            break
        if item in chosen:
            continue
        if all(
            geometry_distance(item, other) >= minimum_distance
            for other in chosen
            if int(other["branch"]) == int(item["branch"])
        ):
            chosen.append(item)

    return sorted(chosen, key=lambda x: x["score"])[:count]


def known_sixbar_reference(
    path: Path,
    target: dict,
    target_fit_dense: dict,
    dense_theta: np.ndarray,
    args,
) -> dict | None:
    if not path.exists():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        best = data.get("best")
        if not best:
            return None

        v, branch = vector_from_record(best)
        record = make_record(
            v,
            branch,
            target,
            target_fit_dense,
            dense_theta,
            args,
        )
        return {
            "known_sixbar_position_rms": best.get("position_rms"),
            "primary_record": record,
        }
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}"}


def build_jobs(
    islands: int,
    seeds: list[dict],
) -> list[dict]:
    jobs = []

    for item in seeds[:islands]:
        v, branch = vector_from_record(item)
        jobs.append({
            "branch": branch,
            "x0": v,
            "seed_source": item.get("seed_source", "v2"),
            "seed_v3_score": float(item["score"]),
        })

    target_plus = islands // 2
    target_minus = islands - target_plus

    plus = sum(job["branch"] == 1 for job in jobs)
    minus = sum(job["branch"] == -1 for job in jobs)

    while len(jobs) < islands:
        if plus < target_plus and (plus <= minus or minus >= target_minus):
            branch = 1
            plus += 1
        elif minus < target_minus:
            branch = -1
            minus += 1
        else:
            branch = 1 if plus < target_plus else -1
            if branch == 1:
                plus += 1
            else:
                minus += 1

        jobs.append({
            "branch": branch,
            "x0": None,
            "seed_source": None,
            "seed_v3_score": None,
        })

    return jobs


def save(path: Path, obj: dict) -> None:
    def conv(x):
        if isinstance(x, Path):
            return str(x)
        if isinstance(x, np.generic):
            return x.item()
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, dict):
            return {k: conv(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [conv(v) for v in x]
        return x

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(conv(obj), indent=2) + "\n",
        encoding="utf-8",
    )


def write_trace(
    path: Path,
    best: dict,
    target: dict,
    dense_theta: np.ndarray,
    target_fit_dense: dict,
    args,
) -> None:
    v = np.array(
        [best["parameters"][k] for k in NAMES],
        dtype=float,
    )
    features = primary_features(
        v,
        int(best["branch"]),
        dense_theta,
        args,
    )
    candidate_fit = fit_fixed_regimes(
        features["theta_deg"],
        features["projected_velocity"],
        target["segmentation"]["boundaries_deg"],
        transition_half_width_deg=args.transition_half_width_deg,
    )

    shape = affine_shape_error(
        candidate_fit["model"],
        target_fit_dense["model"],
        candidate_fit["core_mask"] & target_fit_dense["core_mask"],
    )

    candidate_scaled = (
        shape["velocity_scale"] * features["projected_velocity"]
        + shape["velocity_offset"]
    )
    model_scaled = (
        shape["velocity_scale"] * candidate_fit["model"]
        + shape["velocity_offset"]
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        fields = [
            "theta_deg",
            "E_x",
            "E_y",
            "E_projected_velocity",
            "E_projected_velocity_scaled_to_target",
            "E_piecewise_model_scaled_to_target",
            "target_velocity",
            "target_piecewise_model",
            "regime_core",
        ]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()

        for i in range(len(dense_theta)):
            writer.writerow({
                "theta_deg": float(features["theta_deg"][i]),
                "E_x": float(features["E"][i, 0]),
                "E_y": float(features["E"][i, 1]),
                "E_projected_velocity": float(
                    features["projected_velocity"][i]
                ),
                "E_projected_velocity_scaled_to_target": float(
                    candidate_scaled[i]
                ),
                "E_piecewise_model_scaled_to_target": (
                    float(model_scaled[i])
                    if math.isfinite(float(model_scaled[i]))
                    else ""
                ),
                "target_velocity": float(
                    np.interp(
                        features["theta_deg"][i],
                        target["theta_deg"],
                        target["velocity_smooth"],
                        period=360.0,
                    )
                ),
                "target_piecewise_model": (
                    float(target_fit_dense["model"][i])
                    if math.isfinite(float(target_fit_dense["model"][i]))
                    else ""
                ),
                "regime_core": int(
                    candidate_fit["core_mask"][i]
                    and target_fit_dense["core_mask"][i]
                ),
            })


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=Path, default=TARGET)
    ap.add_argument("--v2-library", type=Path, default=V2_LIBRARY)
    ap.add_argument("--known-sixbar", type=Path, default=KNOWN_SIXBAR)
    ap.add_argument("--output", type=Path, default=OUTPUT)
    ap.add_argument("--trace-output", type=Path, default=TRACE)

    # Deliberately much broader than V2.
    ap.add_argument("--islands", type=int, default=128)
    ap.add_argument("--generations", type=int, default=240)
    ap.add_argument("--population-size", type=int, default=16)
    ap.add_argument("--seed", type=int, default=3952700)

    ap.add_argument("--coarse-step-deg", type=float, default=1.0)
    ap.add_argument("--dense-step-deg", type=float, default=0.25)
    ap.add_argument("--target-smoothing-deg", type=float, default=12.0)
    ap.add_argument("--primary-smoothing-deg", type=float, default=8.0)
    ap.add_argument("--noise-fronts", type=int, default=6)

    ap.add_argument("--minimum-regime-span-deg", type=int, default=25)
    ap.add_argument("--transition-half-width-deg", type=float, default=8.0)

    ap.add_argument("--shape-weight", type=float, default=1.0)
    ap.add_argument("--linearity-weight", type=float, default=1.25)

    ap.add_argument("--primary-sine-floor", type=float, default=0.30)
    ap.add_argument("--preferred-primary-sine", type=float, default=0.35)
    ap.add_argument("--sine-preference-weight", type=float, default=0.12)

    ap.add_argument("--preferred-axis-variance-fraction", type=float, default=0.70)
    ap.add_argument("--axis-preference-weight", type=float, default=0.05)

    ap.add_argument("--E-span-floor", type=float, default=0.75)

    ap.add_argument("--v2-seed-count", type=int, default=32)
    ap.add_argument("--v2-seed-diversity-distance", type=float, default=0.025)

    ap.add_argument("--library-size", type=int, default=400)
    ap.add_argument("--diversity-distance", type=float, default=0.04)
    ap.add_argument("--population-candidates-per-island", type=int, default=20)
    args = ap.parse_args()

    if args.islands < 1:
        ap.error("--islands must be >= 1")
    if not (0.0 < args.primary_sine_floor <= args.preferred_primary_sine <= 1.0):
        ap.error("invalid primary-sine thresholds")
    if args.transition_half_width_deg < 0.0:
        ap.error("--transition-half-width-deg must be >= 0")

    target = load_target(args.target, args)

    print("TARGET LARGE V3 REGIMES", flush=True)
    print(
        "  boundaries = "
        + " / ".join(
            f"{x:.2f}deg"
            for x in target["segmentation"]["boundaries_deg"]
        ),
        flush=True,
    )
    print(
        "  spans      = "
        + " / ".join(
            f"{x:.2f}deg"
            for x in target["segmentation"]["spans_deg"]
        ),
        flush=True,
    )
    print(
        "  slopes     = "
        + " / ".join(
            f"{x:+.6f}/rad"
            for x in target["segmentation"]["slopes_per_rad"]
        ),
        flush=True,
    )
    print(
        f"  3-line target RMS = "
        f"{target['segmentation']['full_fit_rms']:.6f} "
        f"(normalized {target['segmentation']['full_fit_normalized_rms']:.4f})",
        flush=True,
    )

    nc = int(round(360.0 / args.coarse_step_deg))
    nd = int(round(360.0 / args.dense_step_deg))
    coarse_theta = np.linspace(
        0.0, 2.0 * math.pi, nc, endpoint=False
    )
    dense_theta = np.linspace(
        0.0, 2.0 * math.pi, nd, endpoint=False
    )
    coarse_deg = np.mod(np.degrees(coarse_theta), 360.0)
    dense_deg = np.mod(np.degrees(dense_theta), 360.0)

    target_fit_coarse = target_model_on_grid(
        target,
        coarse_deg,
        args,
    )
    target_fit_dense = target_model_on_grid(
        target,
        dense_deg,
        args,
    )

    known = known_sixbar_reference(
        args.known_sixbar,
        target,
        target_fit_dense,
        dense_theta,
        args,
    )
    if known and known.get("primary_record"):
        k = known["primary_record"]
        print(
            "KNOWN 2.08%-CLASS PRIMARY UNDER V3: "
            f"score={k['score']:.5f} "
            f"shape={k['cadence_match']['shape_normalized_rms']:.4f} "
            f"linear={k['cadence_match']['linearity_normalized_rms']:.4f} "
            f"sine={k['minimum_primary_transmission_sine']:.3f}",
            flush=True,
        )

    rescored_v2 = rescore_v2_seeds(
        args.v2_library,
        target,
        target_fit_dense,
        dense_theta,
        args,
    )
    seed_panel = select_v2_seed_panel(
        rescored_v2,
        count=min(args.v2_seed_count, args.islands),
        minimum_distance=args.v2_seed_diversity_distance,
    )

    print(
        f"V2 candidates rescored={len(rescored_v2)}, "
        f"selected as V3 seeds={len(seed_panel)}",
        flush=True,
    )
    for i, item in enumerate(seed_panel[:12]):
        print(
            f"  seed {i:02d}: branch={item['branch']:+d} "
            f"V3={item['score']:.5f} "
            f"V2source={item.get('seed_source')} "
            f"sine={item['minimum_primary_transmission_sine']:.3f}",
            flush=True,
        )
    if len(seed_panel) > 12:
        print("  ...", flush=True)

    jobs = build_jobs(args.islands, seed_panel)

    wall0 = time.time()
    cpu0 = time.process_time()

    report = {
        "description": (
            "Candidate 3952 LARGE primary cadence V3: automatic cyclic "
            "three-regime linear velocity decomposition, with V2 seed reuse "
            "and broad global exploration."
        ),
        "length_unit": "primary_crank_radius",
        "target": str(args.target),
        "target_noise_band": target["noise_band"],
        "target_segmentation": target["segmentation"],
        "target_regime_fit": target["regime_fit"],
        "known_sixbar_reference": known,
        "v2_library": str(args.v2_library),
        "v2_rescored": rescored_v2,
        "v2_seed_panel": seed_panel,
        "bounds": dict(zip(NAMES, BOUNDS)),
        "settings": vars(args),
        "islands": [],
        "candidates": [],
        "best": None,
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
                features = primary_features(
                    v,
                    branch,
                    coarse_theta,
                    args,
                )
                score, match = cadence_score(
                    features,
                    target,
                    target_fit_coarse,
                    args,
                )
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

        vectors = [np.asarray(opt.x, dtype=float)]
        if getattr(opt, "population", None) is not None:
            order = np.argsort(
                np.asarray(opt.population_energies, dtype=float)
            )
            vectors += [
                np.asarray(opt.population[int(j)], dtype=float)
                for j in order[
                    : args.population_candidates_per_island
                ]
            ]
        vectors += [
            x for _, x in sorted(
                archive,
                key=lambda z: z[0],
            )[: args.population_candidates_per_island]
        ]

        dense = []
        seen = set()

        for v in vectors:
            key = tuple(round(float(x), 10) for x in v)
            if key in seen:
                continue
            seen.add(key)

            record = make_record(
                v,
                branch,
                target,
                target_fit_dense,
                dense_theta,
                args,
            )
            if record is not None:
                record["island_index"] = island
                record["island_seed_source"] = job["seed_source"]
                dense.append(record)
                raw_records.append(record)

        dense.sort(key=lambda x: x["score"])
        best = dense[0] if dense else None

        report["islands"].append({
            "index": island,
            "branch": branch,
            "optimizer_seed": run_seed,
            "seed_source": job["seed_source"],
            "seed_v3_score": job["seed_v3_score"],
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
        report["elapsed_wall_hours"] = (
            time.time() - wall0
        ) / 3600.0
        report["elapsed_cpu_hours"] = (
            time.process_time() - cpu0
        ) / 3600.0
        save(args.output, report)

        if best is None:
            print(
                f"island {island:03d} branch={branch:+d}: "
                "no dense candidate",
                flush=True,
            )
        else:
            m = best["cadence_match"]
            print(
                f"island {island:03d} branch={branch:+d} "
                f"{'seeded' if job['seed_source'] else 'fresh ':6s}: "
                f"score={best['score']:.5f} "
                f"shape={m['shape_normalized_rms']:.4f} "
                f"linear={m['linearity_normalized_rms']:.4f} "
                f"slopeRMS={m['normalized_slope_rms']:.4f} "
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
    report["completed_islands"] = len(jobs)
    report["elapsed_wall_hours"] = (
        time.time() - wall0
    ) / 3600.0
    report["elapsed_cpu_hours"] = (
        time.process_time() - cpu0
    ) / 3600.0
    save(args.output, report)

    print("\nPRIMARY CADENCE V3 SUMMARY", flush=True)
    print(
        f"wall={report['elapsed_wall_hours']:.3f} h "
        f"cpu={report['elapsed_cpu_hours']:.3f} h "
        f"islands={len(jobs)} "
        f"diverse={len(library)}",
        flush=True,
    )

    if not library:
        print("No feasible V3 primary found.", flush=True)
        return

    best = library[0]
    m = best["cadence_match"]
    print(
        f"best score={best['score']:.6f}\n"
        f"branch={best['branch']:+d}\n"
        f"shape RMS={m['shape_normalized_rms']:.6f}\n"
        f"linearity RMS={m['linearity_normalized_rms']:.6f}\n"
        f"normalized slope RMS={m['normalized_slope_rms']:.6f}\n"
        f"primary sine={best['minimum_primary_transmission_sine']:.6f}\n"
        f"axis variance fraction={best['axis_position_variance_fraction']:.6f}\n"
        f"projection axis={best['projection_axis_angle_deg']:.3f} deg\n"
        f"AD={best['primary']['ground']:.6f} "
        f"BC={best['primary']['coupler']:.6f} "
        f"CD={best['primary']['rocker']:.6f}",
        flush=True,
    )

    print("regimes:", flush=True)
    for c, t in zip(
        m["candidate_regimes"],
        m["target_regimes_on_same_grid"],
    ):
        print(
            f"  {c['start_deg']:7.2f}->{c['end_deg']:7.2f} deg "
            f"candidate slope={c['slope_per_rad']:+.6f}/rad "
            f"scaled={m['velocity_scale_to_target']*c['slope_per_rad']:+.6f}/rad "
            f"target={t['slope_per_rad']:+.6f}/rad "
            f"nRMS={c['normalized_rms']:.5f}",
            flush=True,
        )

    write_trace(
        args.trace_output,
        best,
        target,
        dense_theta,
        target_fit_dense,
        args,
    )
    print(f"Wrote {args.output}", flush=True)
    print(f"Wrote {args.trace_output}", flush=True)


if __name__ == "__main__":
    main()
