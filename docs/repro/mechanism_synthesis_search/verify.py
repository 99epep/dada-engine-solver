#!/usr/bin/env python3
"""Verify the documentary mechanism-synthesis search policy and saturation result."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import tomllib


HERE = Path(__file__).resolve().parent
POLICY_PATH = HERE / "policy.toml"
SATURATION_PATH = HERE / "saturation_reference.toml"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def close(a: float, b: float, *, atol: float = 1e-12) -> bool:
    return math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=atol)


def load_toml(path: Path) -> dict:
    with path.open("rb") as stream:
        return tomllib.load(stream)


def verify_policy(policy: dict) -> dict:
    require(policy.get("schema_version") == 1, "Unsupported policy schema.")

    geometry = policy["primary_geometry"]
    names = geometry["parameter_names"]
    bounds = geometry["bounds"]

    expected_names = [
        "ground",
        "coupler",
        "rocker",
        "e_along",
        "e_normal",
        "phase_rad",
    ]
    require(names == expected_names, "Unexpected primary coordinate list.")
    require(geometry["branches"] == [-1, 1], "Both primary branches must be present.")
    require(close(geometry["crank"], 1.0), "Crank normalization must be one.")

    for name in names:
        interval = bounds[name]
        require(len(interval) == 2, f"{name}: expected a two-value bound.")
        lo, hi = map(float, interval)
        require(math.isfinite(lo) and math.isfinite(hi), f"{name}: non-finite bound.")
        require(lo < hi, f"{name}: invalid bound order.")

    target = policy["reference_target"]
    require(
        close(
            target["short_branch_span_deg"] + target["long_branch_span_deg"],
            360.0,
            atol=1e-9,
        ),
        "Reference branch spans must sum to one revolution.",
    )
    require(
        target["fast_slow_mean_speed_ratio"] > 1.0,
        "Reference short branch must contain distinct fast and slow cadences.",
    )
    require(
        0.0 < target["fast_displacement_fraction"] < 1.0,
        "Fast displacement fraction must lie inside (0, 1).",
    )

    hard = policy["hard_constraints"]
    soft = policy["soft_preferences"]

    require(hard["full_revolution_closure"] is True, "Closure must be mandatory.")
    require(
        0.0 < hard["minimum_primary_transmission_sine"] < 1.0,
        "Invalid transmission-sine floor.",
    )
    require(
        soft["preferred_primary_transmission_sine"]
        >= hard["minimum_primary_transmission_sine"],
        "Soft transmission preference cannot be below the hard floor.",
    )
    require(
        hard["minimum_e_span_over_crank"] > 0.0,
        "E-span floor must be positive.",
    )
    require(
        0.0 < soft["preferred_axis_variance_fraction"] <= 1.0,
        "Axis-variance preference must lie inside (0, 1].",
    )

    long_branch = policy["long_branch"]
    require(
        long_branch["criterion"] == "mirror_symmetry",
        "The retained long-branch proxy must be mirror symmetry.",
    )
    require(
        long_branch["uses_raw_projected_velocity"] is True,
        "Long-branch symmetry must use raw projected velocity.",
    )
    require(
        long_branch["smoothness_penalty"] is False
        and long_branch["curvature_penalty"] is False
        and long_branch["local_extrema_penalty"] is False
        and long_branch["target_profile_matching"] is False,
        "The long branch must remain free of rejected shape constraints.",
    )

    saturation = policy["fresh_island_saturation"]
    expected_population = len(names) * int(saturation["population_size"])
    require(
        expected_population == saturation["expected_initial_population_per_island"],
        "Initial DE population does not match dimensions × population size.",
    )
    require(
        saturation["historical_seeds"] == 0,
        "Fresh-island saturation must not use historical seeds.",
    )
    require(
        saturation["branch_policy"] == "strict_alternation",
        "Fresh-island saturation must balance the two branches.",
    )

    clustering = policy["family_clustering"]
    require(clustering["branches_are_separate"] is True, "Branches must cluster separately.")
    require(clustering["phase_is_cyclic"] is True, "Phase distance must be cyclic.")
    require(clustering["connected_components"] is True, "Family clustering must be order independent.")
    require(
        clustering["edge_rule"] == "distance_strictly_less_than_threshold",
        "Unexpected clustering edge rule.",
    )
    require(
        clustering["reference_threshold"] in clustering["reported_thresholds"],
        "Reference clustering threshold must be among the reported thresholds.",
    )

    return {
        "dimensions": len(names),
        "initial_population_per_island": expected_population,
    }


def verify_saturation(reference: dict) -> dict:
    require(reference.get("schema_version") == 1, "Unsupported saturation schema.")

    counts = [int(value) for value in reference["family_counts_descending"]]
    require(counts == sorted(counts, reverse=True), "Family counts must be descending.")
    require(all(value > 0 for value in counts), "Family counts must be positive.")

    total = sum(counts)
    require(total == int(reference["total_captures"]), "Family counts do not sum to total captures.")

    observed = len(counts)
    singletons = sum(value == 1 for value in counts)
    doubletons = sum(value == 2 for value in counts)

    unseen = singletons / total
    top1 = sum(counts[:1]) / total
    top4 = sum(counts[:4]) / total
    top10 = sum(counts[:10]) / total

    require(observed == reference["expected_observed_families"], "Observed-family count mismatch.")
    require(singletons == reference["expected_singletons"], "Singleton count mismatch.")
    require(doubletons == reference["expected_doubletons"], "Doubleton count mismatch.")

    require(
        close(unseen, reference["expected_good_turing_unseen_capture_mass"]),
        "Good-Turing unseen capture-mass mismatch.",
    )
    require(
        close(top1, reference["expected_top1_capture_share"]),
        "Top-1 capture share mismatch.",
    )
    require(
        close(top4, reference["expected_top4_capture_share"]),
        "Top-4 capture share mismatch.",
    )
    require(
        close(top10, reference["expected_top10_capture_share"]),
        "Top-10 capture share mismatch.",
    )

    return {
        "total": total,
        "observed": observed,
        "singletons": singletons,
        "doubletons": doubletons,
        "unseen": unseen,
        "top1": top1,
        "top4": top4,
        "top10": top10,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check",
        action="store_true",
        help="Run regression checks and print only a compact success summary.",
    )
    args = parser.parse_args()

    policy = load_toml(POLICY_PATH)
    reference = load_toml(SATURATION_PATH)

    p = verify_policy(policy)
    s = verify_saturation(reference)

    if args.check:
        print(
            "OK mechanism_synthesis_search "
            f"dimensions={p['dimensions']} "
            f"population={p['initial_population_per_island']} "
            f"captures={s['total']} "
            f"families={s['observed']}"
        )
        return

    print("Mechanism-synthesis documentary reproduction")
    print(f"reference motion: {policy['reference_motion_id']}")
    print(
        "primary search: "
        f"{p['dimensions']} dimensions, branches -1/+1, "
        f"{p['initial_population_per_island']} initial DE individuals/island"
    )
    print(
        "saturation @ d=0.04: "
        f"N={s['total']}, families={s['observed']}, "
        f"singletons={s['singletons']}, doubletons={s['doubletons']}"
    )
    print(f"Good-Turing unseen capture mass: {100.0 * s['unseen']:.4f}%")
    print(f"top-1 capture share: {100.0 * s['top1']:.4f}%")
    print(f"top-4 capture share: {100.0 * s['top4']:.4f}%")
    print(f"top-10 capture share: {100.0 * s['top10']:.4f}%")


if __name__ == "__main__":
    main()
