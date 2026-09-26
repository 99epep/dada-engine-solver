#!/usr/bin/env python3
"""Third-stage compact-rod search with selected mechanical limits relaxed.

Starts from:
  outputs/sixbar_rank01_compact_rods_wide/best_pair.json

Keeps fixed:
- final rank_01 thermo-5D hardware and gas inventory;
- original rank_01 rod-shortening requirement (rod = original - 2..3 R);
- maximum instantaneous rod/slider-axis angle = 10 deg;
- primary and secondary transmission sine floors = 0.30;
- crank/EFH clearance floor = 0.50;
- H lateral RMS/stroke limit = 0.25;
- all other mechanical screens from the previous campaigns.

Relaxes only:
- stroke/R ceiling: 3.0 -> 3.3;
- H lateral span/stroke ceiling: 0.65 -> 0.75;
- |mean rod/axis angle|: 3 deg -> 4 deg.

The same widened 30-D trust region as the previous "wide" campaign is retained.

Run from repository root:

    PYTHONPATH=src python3 examples/optimize_rank01_compact_piston_rods_relaxed.py
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

import optimize_rank01_compact_piston_rods as base


ROOT = Path.cwd()

ORIGINAL_PAIR = (
    ROOT
    / "outputs"
    / "sixbar_thermo_coupled_3952_fine_hlat25"
    / "rank_01"
    / "best_pair.json"
)

SOURCE_PAIR = (
    ROOT
    / "outputs"
    / "sixbar_rank01_compact_rods_wide"
    / "best_pair.json"
)

OUTPUT = (
    ROOT
    / "outputs"
    / "sixbar_rank01_compact_rods_relaxed"
)


def original_rank01_rod_bounds(_source_vector):
    """Keep rods 2..3 crank radii shorter than the original rank_01 pair."""
    original = base.load_json(ORIGINAL_PAIR)
    x = base.source_vector(original)
    return {
        "small": (float(x[12]) - 3.0, float(x[12]) - 2.0),
        "large": (float(x[27]) - 3.0, float(x[27]) - 2.0),
    }


def keep_source_seed(v, _primary_branch, _second_branch, _theta):
    """Start exactly from the current wide-campaign champion."""
    return np.asarray(v, dtype=float).copy()


def add_default_arg(flag: str, value) -> None:
    if flag not in sys.argv:
        sys.argv.extend([flag, str(value)])


def main():
    if not SOURCE_PAIR.exists():
        raise FileNotFoundError(
            f"Missing wide compact-rod champion: {SOURCE_PAIR}"
        )
    if not ORIGINAL_PAIR.exists():
        raise FileNotFoundError(
            f"Missing original rank_01 pair: {ORIGINAL_PAIR}"
        )

    # Start from the current wide champion.
    base.PAIR_PATH = SOURCE_PAIR
    base.DEFAULT_OUTPUT = OUTPUT
    base.rod_bounds = original_rank01_rod_bounds
    base.aligned_seed_side = keep_source_seed

    # Retain the broader trust region used in stage 2.
    base.fine.BASE_PRIMARY_LENGTH_FRAC = 0.05
    base.fine.BASE_PRIMARY_POINT_FRAC = 0.07
    base.fine.BASE_PRIMARY_POINT_MIN = 0.15
    base.fine.BASE_PRIMARY_PHASE_DEG = 5.0
    base.fine.BASE_PIVOT_SPAN = 0.75
    base.fine.BASE_SECONDARY_LENGTH_FRAC = 0.10
    base.fine.BASE_H_RATIO_SPAN = 0.20

    base.BASE_ROD_SPAN = 0.75
    base.BASE_AXIS_OFFSET_SPAN = 1.00
    base.BASE_AXIS_ANGLE_DEG = 8.0
    base.DEFAULT_RADII = (1.0, 0.70, 0.45, 0.30)

    # Relax ONLY the two retained fine-mechanical limits selected after the
    # wide campaign. side_mechanical_metrics() reads this namespace directly.
    base.fine.STROKE_CEILING = 3.3
    base.fine.H_LATERAL_SPAN_MAX = 0.75
    base.fine.MECH_ARGS.stroke_ceiling = 3.3
    base.fine.MECH_ARGS.h_lateral_span_max = 0.75

    # Keep the serious mechanical safeguards unchanged.
    assert base.fine.PRIMARY_SINE_FLOOR == 0.30
    assert base.fine.SECONDARY_SINE_FLOOR == 0.30
    assert base.fine.CRANK_CLEARANCE_FLOOR == 0.50
    assert base.fine.H_LATERAL_RMS_MAX == 0.25

    add_default_arg("--output-directory", OUTPUT)
    add_default_arg("--evaluations", 1024)
    add_default_arg("--evaluations-per-radius", 256)
    add_default_arg("--radius-multipliers", "1.0,0.70,0.45,0.30")
    add_default_arg("--budget-seconds", 3600)
    add_default_arg("--max-rod-angle-deg", 10.0)
    add_default_arg("--max-mean-rod-angle-deg", 4.0)
    add_default_arg("--seed", 3973001)

    print(
        "Relaxed compact-rod campaign:\n"
        f"  source: {SOURCE_PAIR}\n"
        f"  output: {OUTPUT}\n"
        "  fixed thermo: final rank_01 thermo-5D\n"
        "  rods: original rank_01 minus 2..3 R\n"
        "  stroke/R ceiling: 3.3\n"
        "  H lateral span/stroke ceiling: 0.75\n"
        "  |mean rod/axis angle| ceiling: 4 deg\n"
        "  max instantaneous rod/axis angle: 10 deg\n"
        "  transmission sine floors: 0.30 / 0.30\n"
        "  crank/EFH clearance floor: 0.50 R\n"
        "  H lateral RMS/stroke ceiling: 0.25\n",
        flush=True,
    )

    base.main()


if __name__ == "__main__":
    main()
