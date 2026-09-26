#!/usr/bin/env python3
"""Wide second-stage search around the compact rank_01 piston-rod champion.

This is deliberately a thin driver around
examples/optimize_rank01_compact_piston_rods.py so that the tested thermodynamic
evaluation, mechanical screening, warm starts and result format remain exactly
the same.

Differences from the first compact-rod campaign:
- start from outputs/sixbar_rank01_compact_rods/best_pair.json;
- keep rod-length bounds referenced to the ORIGINAL rank_01 pair
  (original rod minus 2..3 crank radii);
- widen the local 30-D mechanical search substantially;
- relax rod-axis limits from 8/2 deg to 10/3 deg by default;
- write to a new output directory.

Run from repository root:

    PYTHONPATH=src python3 examples/optimize_rank01_compact_piston_rods_wide.py
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

COMPACT_PAIR = (
    ROOT
    / "outputs"
    / "sixbar_rank01_compact_rods"
    / "best_pair.json"
)

WIDE_OUTPUT = (
    ROOT
    / "outputs"
    / "sixbar_rank01_compact_rods_wide"
)


def original_rank01_rod_bounds(_source_vector):
    """Preserve the original requirement: rods are original rank_01 minus 2..3 R."""
    original = base.load_json(ORIGINAL_PAIR)
    x = base.source_vector(original)
    return {
        "small": (float(x[12]) - 3.0, float(x[12]) - 2.0),
        "large": (float(x[27]) - 3.0, float(x[27]) - 2.0),
    }


def keep_current_compact_seed(v, _primary_branch, _second_branch, _theta):
    """The first campaign already found the compact/aligned seed; start there."""
    return np.asarray(v, dtype=float).copy()


def add_default_arg(flag: str, value) -> None:
    """Insert a CLI default only when the caller did not explicitly provide it."""
    if flag not in sys.argv:
        sys.argv.extend([flag, str(value)])


def main():
    if not COMPACT_PAIR.exists():
        raise FileNotFoundError(
            f"Missing first-stage compact champion: {COMPACT_PAIR}\n"
            "Run optimize_rank01_compact_piston_rods.py first."
        )
    if not ORIGINAL_PAIR.exists():
        raise FileNotFoundError(
            f"Missing original rank_01 pair: {ORIGINAL_PAIR}"
        )

    # Start from the current compact champion, while preserving the original
    # 2..3-R shortening requirement through the custom rod-bounds function.
    base.PAIR_PATH = COMPACT_PAIR
    base.DEFAULT_OUTPUT = WIDE_OUTPUT
    base.rod_bounds = original_rank01_rod_bounds
    base.aligned_seed_side = keep_current_compact_seed

    # Wider 30-D geometry trust region.
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

    # Broader defaults; explicit CLI options still win.
    add_default_arg("--output-directory", WIDE_OUTPUT)
    add_default_arg("--evaluations", 160)
    add_default_arg("--evaluations-per-radius", 40)
    add_default_arg("--radius-multipliers", "1.0,0.70,0.45,0.30")
    add_default_arg("--budget-seconds", 3600)
    add_default_arg("--max-rod-angle-deg", 10.0)
    add_default_arg("--max-mean-rod-angle-deg", 3.0)
    add_default_arg("--seed", 3972001)

    print(
        "Wide compact-rod campaign:\n"
        f"  source: {COMPACT_PAIR}\n"
        f"  output: {WIDE_OUTPUT}\n"
        "  primary lengths: +/-5%\n"
        "  E point: +/-7% (minimum +/-0.15R)\n"
        "  primary phase: +/-5 deg\n"
        "  G pivot: +/-0.75R\n"
        "  EF/GF: +/-10%\n"
        "  H ratios: +/-0.20\n"
        "  slider offset: +/-1.0R\n"
        "  slider angle: +/-8 deg\n"
        "  rod alignment limits: max 10 deg, |mean| 3 deg\n"
        "  rod lengths: still original rank_01 minus 2..3R\n",
        flush=True,
    )

    base.main()


if __name__ == "__main__":
    main()
