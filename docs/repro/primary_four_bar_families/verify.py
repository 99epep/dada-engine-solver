#!/usr/bin/env python3
"""Reproduce the geometric properties documented in PRIMARY_FOUR_BAR_FAMILIES.md."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import tomllib

import numpy as np
from scipy.optimize import minimize_scalar


DATA = Path(__file__).with_name("families.toml")
COARSE_SAMPLES = 1440
TRANSMISSION_TOLERANCE = 1.0e-10
RATIO_TOLERANCE = 1.0e-12


def load_families() -> list[dict]:
    data = tomllib.loads(DATA.read_text(encoding="utf-8"))
    if data.get("schema_version") != 1:
        raise ValueError("Unsupported primary-family catalogue schema.")
    if data.get("length_unit") != "crank_radius" or data.get("angle_unit") != "rad":
        raise ValueError("Expected crank-radius lengths and radian angles.")
    families = data.get("family")
    if not isinstance(families, list) or not families:
        raise ValueError("The catalogue must contain at least one family.")
    ids = [item.get("id") for item in families]
    if any(not isinstance(item, str) or not item for item in ids) or len(ids) != len(
        set(ids)
    ):
        raise ValueError("Family IDs must be unique non-empty strings.")
    return families


def primary_state(
    family: dict,
    theta: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return B, C and E for one study angle, with AB normalized to 1."""
    g = float(family["ground"])
    c = float(family["coupler"])
    r = float(family["rocker"])
    branch = int(family["branch"])
    if branch not in (-1, 1):
        raise ValueError(f"{family['id']}: branch must be -1 or +1.")
    if min(g, c, r) <= 0.0:
        raise ValueError(f"{family['id']}: primary link lengths must be positive.")

    angle = float(theta) + float(family["phase_rad"])
    b = np.array((math.cos(angle), math.sin(angle)), dtype=float)
    d = np.array((g, 0.0), dtype=float)
    delta = d - b
    distance = float(np.linalg.norm(delta))
    if distance <= 0.0:
        raise ValueError(f"{family['id']}: coincident circle centers.")

    along = (c * c - r * r + distance * distance) / (2.0 * distance)
    height2 = c * c - along * along
    if height2 <= 0.0:
        raise ValueError(f"{family['id']}: primary loop cannot close over this angle.")

    unit = delta / distance
    left = np.array((-unit[1], unit[0]), dtype=float)
    c_point = b + along * unit + branch * math.sqrt(height2) * left

    bc = c_point - b
    e_point = (
        b
        + float(family["e_along"]) / c * bc
        + float(family["e_normal"])
        / c
        * np.array((-bc[1], bc[0]), dtype=float)
    )
    return b, c_point, e_point


def transmission_sine(family: dict, theta: float) -> float:
    b, c_point, _ = primary_state(family, theta)
    d = np.array((float(family["ground"]), 0.0), dtype=float)
    bc = c_point - b
    dc = c_point - d
    denominator = float(family["coupler"]) * float(family["rocker"])
    cross = bc[0] * dc[1] - bc[1] * dc[0]
    return abs(float(cross)) / denominator


def minimum_transmission_sine(family: dict) -> float:
    """Find the full-cycle minimum after bracketing it on a periodic scan."""
    step = 2.0 * math.pi / COARSE_SAMPLES
    angles = np.arange(COARSE_SAMPLES, dtype=float) * step
    values = np.array(
        [transmission_sine(family, float(theta)) for theta in angles]
    )
    index = int(np.argmin(values))
    center = float(angles[index])

    result = minimize_scalar(
        lambda theta: transmission_sine(family, float(theta)),
        bounds=(center - step, center + step),
        method="bounded",
        options={"xatol": 1.0e-14},
    )
    if not result.success:
        raise RuntimeError(f"{family['id']}: transmission minimum search failed.")
    return float(result.fun)


def be_over_bc(family: dict) -> float:
    return math.hypot(
        float(family["e_along"]),
        float(family["e_normal"]),
    ) / float(family["coupler"])


def verify_family(family: dict) -> dict:
    # This scan also verifies full-revolution closure for the selected branch.
    for theta in np.linspace(
        0.0,
        2.0 * math.pi,
        COARSE_SAMPLES,
        endpoint=False,
    ):
        primary_state(family, float(theta))

    minimum_sine = minimum_transmission_sine(family)
    ratio = be_over_bc(family)

    return {
        "id": family["id"],
        "branch": int(family["branch"]),
        "phase_deg": math.degrees(float(family["phase_rad"])),
        "minimum_transmission_sine": minimum_sine,
        "be_over_bc": ratio,
        "expected_minimum_transmission_sine": float(
            family["expected_minimum_transmission_sine"]
        ),
        "expected_be_over_bc": float(family["expected_be_over_bc"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail if reproduced values differ from the catalogue references",
    )
    args = parser.parse_args()

    rows = [verify_family(family) for family in load_families()]

    print(
        f"{'family':<28} {'branch':>6} {'phase_deg':>12} "
        f"{'min|sin(mu)|':>14} {'BE/BC':>12}"
    )
    for row in rows:
        print(
            f"{row['id']:<28} "
            f"{row['branch']:>6d} "
            f"{row['phase_deg']:>12.6f} "
            f"{row['minimum_transmission_sine']:>14.9f} "
            f"{row['be_over_bc']:>12.9f}"
        )

    if args.check:
        failures = []
        for row in rows:
            if not math.isclose(
                row["minimum_transmission_sine"],
                row["expected_minimum_transmission_sine"],
                rel_tol=0.0,
                abs_tol=TRANSMISSION_TOLERANCE,
            ):
                failures.append(
                    f"{row['id']}: minimum transmission sine "
                    f"{row['minimum_transmission_sine']:.12g} != "
                    f"{row['expected_minimum_transmission_sine']:.12g}"
                )

            if not math.isclose(
                row["be_over_bc"],
                row["expected_be_over_bc"],
                rel_tol=0.0,
                abs_tol=RATIO_TOLERANCE,
            ):
                failures.append(
                    f"{row['id']}: BE/BC "
                    f"{row['be_over_bc']:.12g} != "
                    f"{row['expected_be_over_bc']:.12g}"
                )

        if failures:
            raise SystemExit("\n".join(failures))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
