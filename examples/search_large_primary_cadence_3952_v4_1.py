#!/usr/bin/env python3
"""V4.1 primary-cadence search for candidate 3952, LARGE cylinder.

V4.1 keeps the V4 physical decomposition:
  * exactly two meaningful turnarounds / two monotonic branches;
  * short high->low branch: FAST cadence then SLOW cadence;
  * long low->high branch: BP exchange branch.

The only conceptual change from V4 concerns the long BP branch.

V4 compared the long-branch velocity profile to candidate 3952.  V4.1 drops
that requirement.  The long branch is allowed to be strongly non-uniform or
oscillatory (for example FAST -> SLOW -> FAST), provided that:
  * position remains monotonic (already enforced by the V4 monotonicity term);
  * the branch is approximately mirror-symmetric about its own midpoint.

NO smoothness term, NO curvature term, NO penalty on the number of local
velocity extrema is used on the long branch.

Mirror asymmetry is measured on the candidate's own long branch, between its
own two matched zero crossings, so turnaround-timing error and BP symmetry are
independent criteria.

This script reuses the tested V4 machinery and changes only:
  * primary_features(): also retains exact unsmoothed projected velocity;
  * cadence_metrics(): replaces target-profile matching by mirror asymmetry.

The old key ``long_profile_normalized_rms`` is retained as a compatibility
alias for ``long_mirror_asymmetry_rms`` because the V4 driver uses that key in
its score/console output.  In V4.1, "long" therefore means MIRROR ASYMMETRY,
not resemblance to 3952.

Recommended first run (fast, no DE):
    PYTHONPATH=src python3 examples/search_large_primary_cadence_3952_v4_1.py \
        --rescore-only

This rescales/re-ranks V4, V3, V2 and the known old 2.08%-class primary under
the V4.1 criterion.

Full global search afterwards, if useful:
    PYTHONPATH=src python3 examples/search_large_primary_cadence_3952_v4_1.py

Full-search outputs:
    outputs/large_primary_cadence_3952_v4_1.json
    outputs/large_primary_cadence_3952_v4_1_best_trace.csv

Rescore-only output:
    outputs/large_primary_cadence_3952_v4_1_rescore.json
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np

import search_large_primary_cadence_3952_v4 as v4


ROOT = Path.cwd()
V4_LIBRARY = ROOT / "outputs" / "large_primary_cadence_3952_v4.json"
V3_LIBRARY = ROOT / "outputs" / "large_primary_cadence_3952_v3.json"
V2_LIBRARY = ROOT / "outputs" / "large_primary_cadence_3952_v2.json"
KNOWN_SIXBAR = ROOT / "outputs" / "large_sixbar_3952.json"

OUTPUT = ROOT / "outputs" / "large_primary_cadence_3952_v4_1.json"
TRACE = ROOT / "outputs" / "large_primary_cadence_3952_v4_1_best_trace.csv"
RESCORE_OUTPUT = ROOT / "outputs" / "large_primary_cadence_3952_v4_1_rescore.json"


_ORIGINAL_CADENCE_METRICS = v4.cadence_metrics


def primary_features_v41(
    values: np.ndarray,
    branch: int,
    theta: np.ndarray,
    args,
) -> dict:
    """V4 primary features plus exact, unsmoothed projected velocity."""
    E, Ed, sine, pdata = v4.stage.primary(theta, values, branch)
    theta_deg = np.mod(np.degrees(theta), 360.0)

    axis, axis_fraction = v4.principal_axis(E)
    projected_velocity_raw = Ed @ axis
    projected_velocity = v4.smooth_periodic(
        projected_velocity_raw,
        360.0 / len(theta),
        args.primary_smoothing_deg,
    )

    X = E - np.mean(E, axis=0)
    normal = np.array((-axis[1], axis[0]))
    along = X @ axis
    lateral = X @ normal

    return {
        "theta_deg": theta_deg,
        "E": E,
        "axis_angle_deg": math.degrees(math.atan2(axis[1], axis[0])),
        "axis_position_variance_fraction": axis_fraction,
        "projected_velocity_raw": projected_velocity_raw,
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
            np.sum(
                np.linalg.norm(
                    np.roll(E, -1, axis=0) - E,
                    axis=1,
                )
            )
        ),
        "BE_over_BC": float(
            math.hypot(
                float(pdata["E_along"]),
                float(pdata["E_normal"]),
            )
            / max(float(pdata["coupler"]), 1e-12)
        ),
    }


def _match_candidate_turnarounds(
    theta_deg: np.ndarray,
    signed_smoothed_velocity: np.ndarray,
    target: dict,
):
    """Match candidate zero crossings to target high/low turns."""
    crossings = v4.zero_crossings_periodic(
        theta_deg,
        signed_smoothed_velocity,
    )
    if len(crossings) < 2:
        return None

    target_high = target["turnarounds_deg"]["high_position"]
    target_low = target["turnarounds_deg"]["low_position"]

    best = None
    for i in range(len(crossings)):
        for j in range(len(crossings)):
            if i == j:
                continue
            err_high = v4.cyclic_distance_deg(
                crossings[i],
                target_high,
            )
            err_low = v4.cyclic_distance_deg(
                crossings[j],
                target_low,
            )
            rms = math.sqrt(
                0.5 * (err_high * err_high + err_low * err_low)
            )
            if best is None or rms < best[0]:
                best = (
                    rms,
                    float(crossings[i]),
                    float(crossings[j]),
                )

    return best


def long_mirror_symmetry_metrics(
    features: dict,
    target: dict,
    sign_multiplier: int,
    args,
) -> dict:
    """Mirror symmetry of the candidate's own BP/long branch.

    No smoothness or local-extrema penalty is present.  The raw analytic
    projected velocity is deliberately used so symmetric oscillations are
    preserved rather than filtered away.
    """
    theta_deg = features["theta_deg"]

    signed_smooth = (
        sign_multiplier * features["projected_velocity"]
    )
    signed_raw = (
        sign_multiplier * features["projected_velocity_raw"]
    )

    match = _match_candidate_turnarounds(
        theta_deg,
        signed_smooth,
        target,
    )
    if match is None:
        return {
            "long_mirror_asymmetry_rms": 2.0,
            "candidate_long_span_deg": None,
            "candidate_long_midpoint_deg": None,
        }

    _, high_turn, low_turn = match

    # With target sign convention: high->low is the short branch, low->high
    # is the long BP return branch.
    short_span = (low_turn - high_turn) % 360.0
    long_span = 360.0 - short_span

    guard = min(
        float(args.turnaround_guard_deg),
        0.20 * long_span,
    )
    half_usable = 0.5 * long_span - guard
    if half_usable <= 1.0:
        return {
            "long_mirror_asymmetry_rms": 2.0,
            "candidate_long_span_deg": float(long_span),
            "candidate_long_midpoint_deg": float(
                (low_turn + 0.5 * long_span) % 360.0
            ),
        }

    # Points equally distant from the two ends of the candidate's own long
    # branch.  A FAST->SLOW->FAST law can score exactly zero here.
    n = max(32, int(getattr(args, "long_profile_points", 64)))
    distance = np.linspace(
        guard,
        0.5 * long_span,
        n,
        endpoint=True,
    )

    left_deg = (low_turn + distance) % 360.0
    right_deg = (high_turn - distance) % 360.0

    left_v = v4.interp_periodic(
        theta_deg,
        signed_raw,
        left_deg,
    )
    right_v = v4.interp_periodic(
        theta_deg,
        signed_raw,
        right_deg,
    )

    numerator = float(
        np.sqrt(np.mean((left_v - right_v) ** 2))
    )
    denominator = float(
        np.sqrt(
            np.mean(
                0.5 * (left_v * left_v + right_v * right_v)
            )
        )
    )
    asymmetry = numerator / max(denominator, 1e-12)

    return {
        "long_mirror_asymmetry_rms": float(asymmetry),
        "candidate_long_span_deg": float(long_span),
        "candidate_long_midpoint_deg": float(
            (low_turn + 0.5 * long_span) % 360.0
        ),
    }


def cadence_metrics_v41(
    features: dict,
    target: dict,
    sign_multiplier: int,
):
    """V4 short-branch metrics + V4.1 free-form BP symmetry."""
    base = _ORIGINAL_CADENCE_METRICS(
        features,
        target,
        sign_multiplier,
    )

    symmetry = long_mirror_symmetry_metrics(
        features,
        target,
        sign_multiplier,
        _ACTIVE_ARGS,
    )

    # Remove conceptual reliance on the target's long-branch velocity shape.
    old_target_profile_error = base.pop(
        "long_profile_normalized_rms",
        None,
    )
    base["v4_target_long_profile_rms_diagnostic_only"] = (
        old_target_profile_error
    )

    base.update(symmetry)

    # Compatibility alias used by V4 cadence_score/main.  In V4.1 this value
    # is mirror asymmetry, NOT target-profile RMS.
    base["long_profile_normalized_rms"] = (
        symmetry["long_mirror_asymmetry_rms"]
    )

    return base


# The V4 driver calls cadence_metrics(features,target,sign) without args.
_ACTIVE_ARGS = None


def _install_v41(args):
    global _ACTIVE_ARGS
    _ACTIVE_ARGS = args
    v4.primary_features = primary_features_v41
    v4.cadence_metrics = cadence_metrics_v41


def default_args():
    """Settings matching the pushed V4 defaults."""
    return SimpleNamespace(
        target=v4.TARGET,
        v2_library=V2_LIBRARY,
        v3_library=V3_LIBRARY,
        known_sixbar=KNOWN_SIXBAR,
        output=OUTPUT,
        trace_output=TRACE,
        islands=128,
        generations=260,
        population_size=16,
        seed=3953100,
        coarse_step_deg=1.0,
        dense_step_deg=0.25,
        target_smoothing_deg=12.0,
        primary_smoothing_deg=8.0,
        noise_fronts=6,
        turnaround_guard_deg=8.0,
        turnaround_sign_guard_deg=3.0,
        transition_guard_deg=6.0,
        minimum_short_subphase_deg=25.0,
        long_profile_points=64,
        monotonicity_weight=2.0,
        endpoint_zero_weight=0.40,
        turning_weight=0.65,
        turning_scale_deg=15.0,
        extra_crossing_weight=0.50,
        speed_ratio_weight=0.55,
        displacement_fraction_weight=0.80,
        # Reuse V4's score weight, but it now multiplies MIRROR ASYMMETRY.
        long_shape_weight=0.45,
        primary_sine_floor=0.30,
        preferred_primary_sine=0.35,
        sine_preference_weight=0.12,
        preferred_axis_variance_fraction=0.70,
        axis_preference_weight=0.03,
        E_span_floor=0.75,
        historical_seed_count=32,
        seed_diversity_distance=0.025,
        library_size=400,
        diversity_distance=0.04,
        population_candidates_per_island=20,
    )


def rescore_only():
    args = default_args()
    _install_v41(args)

    target = v4.load_target(args.target, args)

    nd = int(round(360.0 / args.dense_step_deg))
    dense_theta = np.linspace(
        0.0,
        2.0 * math.pi,
        nd,
        endpoint=False,
    )

    rescored = v4.rescore_seed_sources(
        [
            (V4_LIBRARY, "v4"),
            (V3_LIBRARY, "v3"),
            (V2_LIBRARY, "v2"),
        ],
        KNOWN_SIXBAR,
        target,
        dense_theta,
        args,
    )

    library = v4.diverse(
        rescored,
        args.library_size,
        args.diversity_distance,
    )

    report = {
        "description": (
            "Candidate 3952 LARGE V4.1 RESCORE ONLY. "
            "Short branch retains V4 fast/slow cadence criteria. "
            "Long BP branch is scored only for mirror symmetry plus the "
            "existing monotonicity constraint; no smoothness/extrema penalty."
        ),
        "mode": "rescore_only",
        "source_v4": str(V4_LIBRARY),
        "source_v3": str(V3_LIBRARY),
        "source_v2": str(V2_LIBRARY),
        "known_sixbar": str(KNOWN_SIXBAR),
        "target": str(args.target),
        "target_structure": {
            "turnarounds_deg": target["turnarounds_deg"],
            "short_span_deg": target["short_span_deg"],
            "long_span_deg": target["long_span_deg"],
            "short_split": target["short_split"],
            "fast_slow_speed_ratio": target[
                "target_fast_slow_speed_ratio"
            ],
            "fast_displacement_fraction": target[
                "target_fast_displacement_fraction"
            ],
        },
        "settings": vars(args),
        "rescored_count": len(rescored),
        "candidates": library,
        "best": library[0] if library else None,
    }

    v4.save(RESCORE_OUTPUT, report)

    print("PRIMARY CADENCE V4.1 — RESCORE ONLY")
    print(
        "Long BP score = mirror asymmetry only "
        "(no smoothness/local-extrema penalty)"
    )
    print(
        f"rescored={len(rescored)} diverse={len(library)} "
        f"output={RESCORE_OUTPUT}"
    )

    print("\nTOP V4.1 EXISTING CANDIDATES")
    for rank, item in enumerate(library[:20]):
        m = item["cadence_match"]
        print(
            f"{rank:02d} "
            f"score={item['score']:.5f} "
            f"branch={item['branch']:+d} "
            f"turn={m['turning_rms_deg']:.2f}deg "
            f"mono={m['monotonicity_wrong_sign_rms']:.4f} "
            f"ratio={m['fast_slow_speed_ratio']:.2f} "
            f"dispFast={m['fast_displacement_fraction']:.3f} "
            f"symBP={m['long_mirror_asymmetry_rms']:.4f} "
            f"zeros={m['zero_crossing_count']} "
            f"sine={item['minimum_primary_transmission_sine']:.3f} "
            f"BE/BC={item['BE_over_BC']:.2f} "
            f"source={item.get('seed_source')}"
        )


def full_search():
    args = default_args()
    _install_v41(args)

    # Reuse V4's tested full-search driver, but seed it primarily from V4
    # instead of V3 and write distinct V4.1 files.
    v4.OUTPUT = OUTPUT
    v4.TRACE = TRACE
    v4.V3_LIBRARY = V4_LIBRARY
    v4.V2_LIBRARY = V2_LIBRARY
    v4.KNOWN_SIXBAR = KNOWN_SIXBAR

    # V4's argparse is still used.  Its option --long-shape-weight now means
    # "long BP mirror-symmetry weight" in V4.1.
    v4.main()


if __name__ == "__main__":
    if "--rescore-only" in sys.argv:
        sys.argv.remove("--rescore-only")
        if len(sys.argv) != 1:
            raise SystemExit(
                "--rescore-only currently uses V4.1 defaults and accepts "
                "no additional command-line options."
            )
        rescore_only()
    else:
        # For a normal run the V4 argparse remains available, including
        # --islands, --generations, etc.
        full_search()
