#!/usr/bin/env python3
"""Deep thermo-5D retuning of the final compact rank_01 mechanism.

Fixed mechanism:
    outputs/sixbar_rank01_compact_rods_hlat90/best_pair.json

Thermo starting point:
    existing final rank_01 thermo-5D champion from
    outputs/sixbar_thermo5d_3952/rank_01/report.json

Free variables only:
    swept_ratio
    n_i
    length_i_m
    n_o
    length_o_m

Everything mechanical is strictly frozen.

The campaign deliberately starts from the already-polished thermo point that
was held fixed during all compact-mechanism searches, then performs a broad-to-
fine 5D incumbent-centred Sobol search.

Default:
    1024 non-seed evaluations
    256 evaluations per radius
    7200 s wall-clock budget

Run from repository root:

    PYTHONPATH=src python3 examples/optimize_rank01_compact_thermo5d_deep.py
"""

from __future__ import annotations

from pathlib import Path
import sys

import optimize_sixbar_pairs_thermo5d_3952_hlat25 as base


ROOT = Path.cwd()

PAIR_PATH = (
    ROOT / "outputs" / "sixbar_rank01_compact_rods_hlat90" / "best_pair.json"
)
PAIR_REPORT = (
    ROOT / "outputs" / "sixbar_rank01_compact_rods_hlat90" / "report.json"
)

THERMO_SOURCE_REPORT = (
    ROOT / "outputs" / "sixbar_thermo5d_3952" / "rank_01" / "report.json"
)

OUTPUT = ROOT / "outputs" / "sixbar_rank01_compact_thermo5d_deep"


def add_default_arg(flag: str, value) -> None:
    if flag not in sys.argv:
        sys.argv.extend([flag, str(value)])


_original_make_basis = base.make_3952_basis


def make_basis_from_current_thermo():
    """Rebuild the 260 K basis and move its thermo seed to the known final point."""
    (
        seed_machine,
        definition,
        total_mass,
        geometry_info,
        _old_thermo_seed,
        basis_identity,
    ) = _original_make_basis()

    thermo_report = base.load_json(THERMO_SOURCE_REPORT)
    thermo_best = thermo_report.get("best_feasible")
    if thermo_best is None:
        raise RuntimeError(
            f"No best_feasible in thermo source report: {THERMO_SOURCE_REPORT}"
        )

    current = {
        "swept_ratio": float(thermo_best["parameters"]["swept_ratio"]),
        "n_i": int(thermo_best["parameters"]["n_i"]),
        "length_i_m": float(thermo_best["parameters"]["length_i_m"]),
        "n_o": int(thermo_best["parameters"]["n_o"]),
        "length_o_m": float(thermo_best["parameters"]["length_o_m"]),
    }

    pair = base.load_json(PAIR_PATH)
    small_mech = base.mechanism_from_candidate(pair["small"])
    large_mech = base.mechanism_from_candidate(pair["large"])

    # Build an exact machine at the current known thermo point so the generic
    # optimizer's seed hardware and wall-capacity warm-start scaling are both
    # referenced to the correct hardware.
    tuned_seed_machine = base.build_design(
        seed_machine,
        geometry_info,
        total_mass,
        current,
        small_mech,
        large_mech,
    )

    identity = dict(basis_identity)
    identity.update({
        "local_thermo_seed_source": str(THERMO_SOURCE_REPORT),
        "local_thermo_seed_candidate_id": thermo_best.get("candidate_id"),
        "local_thermo_seed": current,
        "fixed_pair_source": str(PAIR_PATH),
    })

    return (
        tuned_seed_machine,
        definition,
        total_mass,
        geometry_info,
        current,
        identity,
    )


def main():
    for path in (PAIR_PATH, PAIR_REPORT, THERMO_SOURCE_REPORT):
        if not path.exists():
            raise FileNotFoundError(path)

    pair = base.load_json(PAIR_PATH)
    mech_report = base.load_json(PAIR_REPORT)
    mech_best = mech_report.get("best_feasible")
    if mech_best is None:
        raise RuntimeError(f"No best_feasible in {PAIR_REPORT}")
    if pair.get("candidate_id") != mech_best.get("candidate_id"):
        raise RuntimeError(
            "Compact mechanism best_pair.json and report.json disagree."
        )
    if mech_best.get("last_complete_state") is None:
        raise RuntimeError(
            "Compact mechanism champion has no periodic state for warm start."
        )

    # Restrict the generic thermo-5D engine to this single mechanism.
    base.PAIR_SOURCES = {1: (PAIR_PATH, PAIR_REPORT)}
    base.OUTPUT = OUTPUT
    base.make_3952_basis = make_basis_from_current_thermo

    # Broad enough to let the new piston laws choose noticeably different
    # hardware, but centred on the already-good final thermo champion.
    #
    # Tuple = swept-ratio absolute radius : tube-count radius :
    #         tube-length relative radius
    local_radii = (
        (0.0800, 240, 0.10000),
        (0.0400, 120, 0.05000),
        (0.0200,  60, 0.02500),
        (0.0100,  30, 0.01250),
    )

    add_default_arg("--families", "1")
    add_default_arg("--output-directory", OUTPUT)
    add_default_arg("--evaluations-per-family", 1024)
    add_default_arg("--evaluations-per-radius", 256)
    add_default_arg(
        "--radii",
        ",".join(f"{r}:{n}:{l}" for r, n, l in local_radii),
    )
    add_default_arg("--budget-seconds", 7200)
    add_default_arg("--candidate-seconds", 180)
    add_default_arg("--seed", 3975001)

    print(
        "Deep compact rank_01 thermo-5D retuning:\n"
        f"  fixed mechanism: {PAIR_PATH}\n"
        f"  thermo centre: {THERMO_SOURCE_REPORT}\n"
        f"  output: {OUTPUT}\n"
        "  free coordinates: swept_ratio, n_i, length_i_m, n_o, length_o_m\n"
        "  mechanism: strictly frozen\n"
        "  evaluations: 1024 + exact seed\n"
        "  radii:\n"
        "    0.080 : 240 : 10.0%\n"
        "    0.040 : 120 :  5.0%\n"
        "    0.020 :  60 :  2.5%\n"
        "    0.010 :  30 :  1.25%\n",
        flush=True,
    )

    base.main()


if __name__ == "__main__":
    main()
