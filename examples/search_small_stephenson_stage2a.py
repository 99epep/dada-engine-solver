#!/usr/bin/env python3
"""Stage 2A pilot: fixed primary four-bar + second RR dyad + finite piston rod.

Run from the dada-engine-solver repository root, for example:

    PYTHONPATH=src python3 search_small_stephenson_stage2a.py \
        --iterations 200 \
        --stroke-floor 1.5 \
        --output outputs/small_stephenson_stage2a.json

The primary four-bar is frozen to the current extended small-cylinder F winner.
Its former coupler output point becomes joint E of a second RR dyad E-F-G.
The piston rod attaches at F.

Physical slider convention:
- +axis points from the linkage toward the cylinder head;
- piston is the +sqrt intersection along that axis;
- working chamber lies beyond the piston in +axis direction;
- therefore normalized volume decreases when slider coordinate increases.

Crank radius is the unit length.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import differential_evolution


ROOT = Path.cwd()
TARGET = ROOT / "outputs" / "motor_champion_motion_target.csv"

# Current small-cylinder extended F winner, frozen as primary stage.
PRIMARY = {
    "ground": 2.045193376592419,
    "coupler": 2.1152136168766678,
    "rocker": 1.8202134329709354,
    "branch": -1,
    "e_along": 1.873315626010214,
    "e_normal": -1.5362965718510164,
    "phase": -3.0863675667472057,
}

PARAMETER_NAMES = [
    "second_pivot_x",
    "second_pivot_y",
    "link_EF",
    "link_GF",
    "piston_rod",
    "slider_axis_offset",
    "slider_axis_angle",
]

BOUNDS = [
    (-8.0, 8.0),                # Gx
    (-8.0, 8.0),                # Gy
    (0.4, 12.0),                # EF
    (0.4, 12.0),                # GF
    (2.0, 30.0),                # piston rod
    (-10.0, 10.0),              # signed perpendicular axis offset
    (-math.pi, math.pi),        # slider axis angle
]


def primary_motion(theta: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Return E, dE/dtheta and minimum primary transmission sine."""
    g = PRIMARY["ground"]
    c = PRIMARY["coupler"]
    r = PRIMARY["rocker"]
    branch = PRIMARY["branch"]

    angle = theta + PRIMARY["phase"]
    B = np.column_stack((np.cos(angle), np.sin(angle)))
    Bd = np.column_stack((-np.sin(angle), np.cos(angle)))
    D = np.array((g, 0.0))

    delta = D - B
    distance = np.linalg.norm(delta, axis=1)
    direction = delta / distance[:, None]
    along = (c * c - r * r + distance * distance) / (2.0 * distance)
    height2 = c * c - along * along
    if np.any(height2 <= 0.0):
        raise ValueError("Frozen primary four-bar does not close.")
    height = np.sqrt(height2)
    left = np.column_stack((-direction[:, 1], direction[:, 0]))

    C = B + along[:, None] * direction + branch * height[:, None] * left

    coupler = C - B
    rocker = C - D
    determinant = (
        coupler[:, 0] * rocker[:, 1]
        - coupler[:, 1] * rocker[:, 0]
    )
    if np.any(np.abs(determinant) <= 1e-10):
        raise ValueError("Frozen primary four-bar reaches a velocity singularity.")

    rhs = np.sum(coupler * Bd, axis=1)
    Cd = np.column_stack(
        (
            rhs * rocker[:, 1] / determinant,
            -rhs * rocker[:, 0] / determinant,
        )
    )

    u = coupler / c
    ud = (Cd - Bd) / c
    n = np.column_stack((-u[:, 1], u[:, 0]))
    nd = np.column_stack((-ud[:, 1], ud[:, 0]))

    E = B + PRIMARY["e_along"] * u + PRIMARY["e_normal"] * n
    Ed = Bd + PRIMARY["e_along"] * ud + PRIMARY["e_normal"] * nd

    minimum_sine = float(np.min(np.abs(determinant) / (c * r)))
    return E, Ed, minimum_sine


def evaluate(
    vector: np.ndarray,
    second_branch: int,
    E: np.ndarray,
    Ed: np.ndarray,
    target_q: np.ndarray,
    target_dq: np.ndarray,
) -> dict[str, object]:
    gx, gy, link_ef, link_gf, rod, axis_offset, axis_angle = map(float, vector)
    G = np.array((gx, gy))

    # Second RR dyad: E-F-G.
    delta = G - E
    distance = np.linalg.norm(delta, axis=1)
    if np.any(distance <= 1e-9):
        raise ValueError("Second dyad centers coincide.")

    along = (
        link_ef * link_ef
        - link_gf * link_gf
        + distance * distance
    ) / (2.0 * distance)
    height2 = link_ef * link_ef - along * along
    if np.any(height2 <= 1e-10):
        raise ValueError("Second dyad cannot close throughout the sampled cycle.")

    height = np.sqrt(height2)
    direction = delta / distance[:, None]
    left = np.column_stack((-direction[:, 1], direction[:, 0]))
    F = (
        E
        + along[:, None] * direction
        + second_branch * height[:, None] * left
    )

    ef = F - E
    gf = F - G
    determinant = ef[:, 0] * gf[:, 1] - ef[:, 1] * gf[:, 0]
    if np.any(np.abs(determinant) <= 1e-10):
        raise ValueError("Second dyad reaches a velocity singularity.")

    rhs = np.sum(ef * Ed, axis=1)
    Fd = np.column_stack(
        (
            rhs * gf[:, 1] / determinant,
            -rhs * gf[:, 0] / determinant,
        )
    )

    secondary_sine = float(
        np.min(np.abs(determinant) / (link_ef * link_gf))
    )

    # Slider axis. Longitudinal origin is a gauge; only perpendicular offset matters.
    axis = np.array((math.cos(axis_angle), math.sin(axis_angle)))
    normal = np.array((-axis[1], axis[0]))
    origin = axis_offset * normal

    relative = F - origin
    longitudinal = relative @ axis
    transverse = relative @ normal
    longitudinal_d = Fd @ axis
    transverse_d = Fd @ normal

    margin2 = rod * rod - transverse * transverse
    if np.any(margin2 <= 1e-10):
        raise ValueError("Finite piston rod cannot reach slider over full cycle.")
    margin = np.sqrt(margin2)

    # Physical convention: piston beyond linkage toward head.
    slider = longitudinal + margin
    slider_d = longitudinal_d - transverse * transverse_d / margin

    stroke = float(np.ptp(slider))
    if stroke <= 1e-10:
        raise ValueError("Zero piston stroke.")

    # Working chamber lies beyond piston toward +axis, hence volume decreases
    # as slider coordinate increases.
    q = 1.0 - (slider - float(np.min(slider))) / stroke
    dq = -slider_d / stroke

    position_rms = float(np.sqrt(np.mean((q - target_q) ** 2)))
    derivative_rms = float(np.sqrt(np.mean((dq - target_dq) ** 2)))
    score = math.sqrt(position_rms**2 + 0.1 * derivative_rms**2)

    rod_axis_cosine = float(np.min(margin / rod))

    return {
        "score": score,
        "position_rms": position_rms,
        "derivative_rms": derivative_rms,
        "stroke_over_crank": stroke,
        "minimum_secondary_transmission_sine": secondary_sine,
        "minimum_rod_axis_cosine": rod_axis_cosine,
        "parameters": dict(zip(PARAMETER_NAMES, map(float, vector))),
        "second_branch": second_branch,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--iterations", type=int, default=160)
    parser.add_argument("--population-size", type=int, default=10)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--stroke-floor", type=float, default=1.5)
    parser.add_argument("--secondary-sine-floor", type=float, default=0.30)
    parser.add_argument("--rod-cos-floor", type=float, default=0.95)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "outputs" / "small_stephenson_stage2a.json",
    )
    args = parser.parse_args()

    if not TARGET.exists():
        raise SystemExit(f"Target not found: {TARGET}")

    raw = np.genfromtxt(TARGET, delimiter=",", names=True)
    if raw["theta_rad"][-1] >= 2.0 * math.pi - 1e-9:
        raw = raw[:-1]

    dense_theta = raw["theta_rad"]
    dense_q = raw["small_fraction_0_1"]
    dense_dq = raw["small_dq_dtheta_per_rad"]

    # Same coarse 3-degree search spacing as current finite-rod search.
    coarse = raw[::12]
    theta = coarse["theta_rad"]
    target_q = coarse["small_fraction_0_1"]
    target_dq = coarse["small_dq_dtheta_per_rad"]

    E, Ed, primary_sine = primary_motion(theta)
    dense_E, dense_Ed, dense_primary_sine = primary_motion(dense_theta)

    report: dict[str, object] = {
        "description": (
            "Stage 2A: frozen current small-cylinder primary four-bar, "
            "second RR dyad, finite piston rod, physical chamber-side convention."
        ),
        "length_unit": "crank_radius",
        "primary": PRIMARY,
        "physical_slider_convention": {
            "slider_branch": +1,
            "axis_direction": "linkage_toward_cylinder_head",
            "working_chamber": "beyond_piston_along_positive_axis",
            "volume_increases_with_slider_coordinate": False,
        },
        "constraints": {
            "stroke_over_crank_minimum": args.stroke_floor,
            "minimum_secondary_transmission_sine": args.secondary_sine_floor,
            "minimum_rod_axis_cosine": args.rod_cos_floor,
        },
        "primary_minimum_transmission_sine": dense_primary_sine,
        "bounds": dict(zip(PARAMETER_NAMES, BOUNDS)),
        "runs": [],
    }

    start = time.monotonic()

    for run_index, second_branch in enumerate((+1, -1)):
        archive: list[dict[str, object]] = []

        def objective(v: np.ndarray) -> float:
            try:
                result = evaluate(v, second_branch, E, Ed, target_q, target_dq)
            except (ValueError, FloatingPointError):
                return 1.0e4

            violations = (
                max(0.0, args.stroke_floor - result["stroke_over_crank"])
                + max(
                    0.0,
                    args.secondary_sine_floor
                    - result["minimum_secondary_transmission_sine"],
                )
                + max(
                    0.0,
                    args.rod_cos_floor
                    - result["minimum_rod_axis_cosine"],
                )
            )
            if violations > 0.0:
                return 10.0 + 100.0 * violations

            if not archive or result["score"] < archive[-1]["score"]:
                archive.append(result)
            return float(result["score"])

        optimization = differential_evolution(
            objective,
            BOUNDS,
            seed=args.seed + run_index,
            popsize=args.population_size,
            maxiter=args.iterations,
            polish=True,
            tol=1e-8,
            updating="immediate",
            workers=1,
        )

        dense_candidates: list[dict[str, object]] = []
        for candidate in sorted(archive, key=lambda x: x["score"])[:20]:
            vector = np.array(
                [candidate["parameters"][name] for name in PARAMETER_NAMES],
                dtype=float,
            )
            try:
                dense = evaluate(
                    vector,
                    second_branch,
                    dense_E,
                    dense_Ed,
                    dense_q,
                    dense_dq,
                )
            except (ValueError, FloatingPointError):
                continue

            dense["feasible"] = (
                dense["stroke_over_crank"] >= args.stroke_floor
                and dense["minimum_secondary_transmission_sine"]
                >= args.secondary_sine_floor
                and dense["minimum_rod_axis_cosine"] >= args.rod_cos_floor
            )
            if dense["feasible"]:
                dense_candidates.append(dense)

        best = (
            min(dense_candidates, key=lambda x: x["score"])
            if dense_candidates
            else None
        )

        report["runs"].append(
            {
                "second_branch": second_branch,
                "optimizer_success": bool(optimization.success),
                "optimizer_message": str(optimization.message),
                "function_evaluations": int(optimization.nfev),
                "best_dense": best,
            }
        )

        valid = [
            run["best_dense"]
            for run in report["runs"]
            if run["best_dense"] is not None
        ]
        report["best"] = min(valid, key=lambda x: x["score"]) if valid else None
        report["elapsed_seconds"] = time.monotonic() - start

        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")

        print(
            f"branch {second_branch:+d}:",
            json.dumps(best, indent=2) if best else "no dense feasible candidate",
            flush=True,
        )

    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
