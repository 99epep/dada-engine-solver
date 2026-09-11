"""Stage F6: local seven-variable search around the F5 crossed-flatness ridge.

This stage keeps the F4/F5 crossed before/after symmetry and searches the same
seven quantities:

    small_min_deg
    small_max_deg
    large_max_deg
    h_s_low_before
    h_s_low_after
    h_s_high_before
    h_s_high_after

The L-side endpoint flatnesses remain imposed by the mechanical mirror symmetry:

    h_L_low_before  = h_S_low_after
    h_L_low_after   = h_S_low_before
    h_L_high_before = h_S_high_after
    h_L_high_after  = h_S_high_before

F6 is a local refinement of F5, but deliberately keeps both low-side h ranges
crossing zero.  In particular it retains the F5 rank-9 region where
h_s_low_before is positive.

No thermodynamic equation is changed here.  The additional V_S + V_L diagnostic
is only a kinematic indicator; it must not be interpreted as the actual flow
through either exchanger.  Actual exchanger branch flow must be checked from
the dynamic cycle solution separately.

Default search: 2**8 = 256 Sobol points, with no control evaluations.

Outputs:
    outputs/motor_two_extrema_stageF6.csv
    outputs/motor_two_extrema_stageF6.json
    outputs/motor_two_extrema_stageF6_motion_overlay.png
    outputs/motor_two_extrema_stageF6_total_volume.png

The CSV is checkpointed after every evaluation.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np
from scipy.stats import qmc

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.kinematics import HarmonicVolumeKinematics

from compare_motor_motion_laws_stage7A5 import (
    ROOT,
    A5_CAMPAIGN,
    _candidate_design_and_mass,
    _evaluate,
    _same_inventory_design,
)
from optimize_motor_two_extrema_stageF5 import (
    TAU,
    LARGE_MIN_DEG,
    PARAMETERS,
    _make_kinematics,
    _row,
    _write_csv,
    _load_checkpoint,
    _verify_exact_harmonic,
)


F5_CSV = ROOT / "outputs" / "motor_two_extrema_stageF5.csv"

OUTPUT_CSV = ROOT / "outputs" / "motor_two_extrema_stageF6.csv"
OUTPUT_JSON = ROOT / "outputs" / "motor_two_extrema_stageF6.json"
OUTPUT_PLOT = ROOT / "outputs" / "motor_two_extrema_stageF6_motion_overlay.png"
OUTPUT_TOTAL_PLOT = ROOT / "outputs" / "motor_two_extrema_stageF6_total_volume.png"


# Local F6 box.
#
# Compared with F5:
# - small_max remains open downward because F5 #9 sat near the lower F5 bound;
# - large_max remains open upward enough to include the high-performing ridge;
# - both low-side h parameters still cross zero;
# - h_s_low_before reaches +0.35 to retain and extend the F5 #9 direction;
# - high-side ranges retain the asymmetric hard-before / softer-after family.
BOUNDS = {
    "small_min_deg": (67.0, 79.0),
    "small_max_deg": (233.0, 245.0),
    "large_max_deg": (361.0, 371.0),

    "h_s_low_before": (-0.60, 0.35),
    "h_s_low_after":  (-0.45, 0.30),

    "h_s_high_before": (-1.80, -1.00),
    "h_s_high_after":  (-1.20, -0.55),
}


def _scale(point: np.ndarray) -> dict[str, float]:
    result: dict[str, float] = {}
    for x, name in zip(point, PARAMETERS, strict=True):
        lo, hi = BOUNDS[name]
        result[name] = lo + float(x) * (hi - lo)
    return result


def _verify_monotone_flatness_domain() -> None:
    """Verify monotonicity for the actual four F6 branch endpoint pairings."""
    y = np.linspace(0.0, 1.0, 10001)

    # S min -> max : low_after  -> high_before
    # S max -> min : high_after -> low_before
    # L min -> max : low_before -> high_after
    # L max -> min : high_before -> low_after
    branch_pairs = (
        ("h_s_low_after", "h_s_high_before"),
        ("h_s_high_after", "h_s_low_before"),
        ("h_s_low_before", "h_s_high_after"),
        ("h_s_high_before", "h_s_low_after"),
    )

    for start_name, end_name in branch_pairs:
        for h0 in BOUNDS[start_name]:
            for h1 in BOUNDS[end_name]:
                derivative = (
                    1.0
                    - h0 * (1.0 - 4.0 * y + 3.0 * y * y)
                    + h1 * (2.0 * y - 3.0 * y * y)
                )
                minimum = float(np.min(derivative))
                if minimum < -1e-12:
                    raise RuntimeError(
                        "F6 bounds can create a non-monotone branch: "
                        f"{start_name}={h0}, {end_name}={h1}, "
                        f"minimum dE/dy={minimum}"
                    )


def _load_f5_rows() -> list[dict]:
    if not F5_CSV.exists():
        raise FileNotFoundError(
            f"F5 reference CSV not found: {F5_CSV.relative_to(ROOT)}"
        )

    with F5_CSV.open(newline="") as stream:
        raw_rows = list(csv.DictReader(stream))

    rows: list[dict] = []
    for raw in raw_rows:
        row = dict(raw)
        row["index"] = int(float(row["index"]))
        row["feasible"] = row["feasible"].strip().lower() == "true"
        for key in row:
            if key in ("label", "feasible"):
                continue
            try:
                row[key] = float(row[key])
            except (TypeError, ValueError):
                pass
        rows.append(row)

    return rows


def _f5_references() -> tuple[dict, dict]:
    rows = _load_f5_rows()

    feasible = [
        row
        for row in rows
        if row["feasible"]
        and isinstance(row.get("indicated_thermal_efficiency"), float)
        and math.isfinite(row["indicated_thermal_efficiency"])
    ]
    if not feasible:
        raise RuntimeError("F5 contains no feasible reference row.")

    feasible.sort(
        key=lambda row: row["indicated_thermal_efficiency"],
        reverse=True,
    )
    champion = feasible[0]

    # The deliberately retained soft-low candidate discussed after F5.
    soft = next((row for row in rows if row["index"] == 40), None)
    if soft is None:
        raise RuntimeError("F5 index 40 was not found.")
    if not soft["feasible"]:
        raise RuntimeError("F5 index 40 is not feasible in the saved CSV.")

    return champion, soft


def _kinematic_arrays(kinematics, angle_rad: np.ndarray):
    small = np.array([
        kinematics.small_cylinder_volume(float(x))
        for x in angle_rad
    ])
    large = np.array([
        kinematics.large_cylinder_volume(float(x))
        for x in angle_rad
    ])
    dsmall = np.array([
        kinematics.small_cylinder_volume_derivative(float(x))
        for x in angle_rad
    ])
    dlarge = np.array([
        kinematics.large_cylinder_volume_derivative(float(x))
        for x in angle_rad
    ])
    return small, large, dsmall, dlarge


def _row_values(row: dict) -> dict[str, float]:
    return {name: float(row[name]) for name in PARAMETERS}


def _plot_results(
    small_limits,
    large_limits,
    rows: list[dict],
    f5_champion: dict,
    f5_soft: dict,
) -> None:
    feasible = [
        row
        for row in rows
        if bool(row["feasible"])
        and row.get("indicated_thermal_efficiency") not in (None, "")
        and math.isfinite(float(row["indicated_thermal_efficiency"]))
    ]
    if not feasible:
        return

    feasible.sort(
        key=lambda row: float(row["indicated_thermal_efficiency"]),
        reverse=True,
    )
    top = feasible[:3]

    cache = Path(tempfile.gettempdir()) / "dada_solver_matplotlib"
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib.pyplot as plt

    angle_deg = np.linspace(0.0, 360.0, 1441)
    angle_rad = np.radians(angle_deg)

    harmonic = HarmonicVolumeKinematics(
        small_limits,
        large_limits,
        math.radians(250.0),
    )
    f5_best_kin = _make_kinematics(
        small_limits,
        large_limits,
        _row_values(f5_champion),
    )
    f5_soft_kin = _make_kinematics(
        small_limits,
        large_limits,
        _row_values(f5_soft),
    )

    # ------------------------------------------------------------------
    # Physical cylinder volumes
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(2, 1, figsize=(12, 11), sharex=True)

    references = [
        ("Exact harmonic 250°", harmonic, "--", 1.8),
        (
            f"F5 champion #{1}, index {int(f5_champion['index'])}, "
            f"eta={100*float(f5_champion['indicated_thermal_efficiency']):.4f}%",
            f5_best_kin,
            ":",
            2.0,
        ),
        (
            f"F5 soft-low #9, index {int(f5_soft['index'])}, "
            f"eta={100*float(f5_soft['indicated_thermal_efficiency']):.4f}%",
            f5_soft_kin,
            "-.",
            2.0,
        ),
    ]

    for label, kin, linestyle, linewidth in references:
        small, large, _, _ = _kinematic_arrays(kin, angle_rad)
        axes[0].plot(
            angle_deg,
            small * 1000.0,
            linestyle=linestyle,
            linewidth=linewidth,
            label=label,
        )
        axes[1].plot(
            angle_deg,
            large * 1000.0,
            linestyle=linestyle,
            linewidth=linewidth,
            label=label,
        )

    for rank, row in enumerate(top, 1):
        kin = _make_kinematics(
            small_limits,
            large_limits,
            _row_values(row),
        )
        small, large, _, _ = _kinematic_arrays(kin, angle_rad)
        label = (
            f"F6 #{rank}, index {int(row['index'])}, "
            f"eta={100*float(row['indicated_thermal_efficiency']):.4f}%"
        )
        axes[0].plot(angle_deg, small * 1000.0, label=label)
        axes[1].plot(angle_deg, large * 1000.0, label=label)

    axes[0].set_ylabel("Small-cylinder volume (L)")
    axes[1].set_ylabel("Large-cylinder volume (L)")
    axes[1].set_xlabel("Study angle (deg)")

    for axis in axes:
        axis.set_xlim(0.0, 360.0)
        axis.set_xticks(np.arange(0.0, 361.0, 45.0))
        axis.grid(True, alpha=0.3)
        axis.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.15),
            ncol=2,
            frameon=False,
        )

    fig.suptitle(
        "Stage F6 — cylinder motion overlay",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.97), h_pad=3.0)
    OUTPUT_PLOT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PLOT, dpi=160)
    plt.close(fig)

    # ------------------------------------------------------------------
    # Cylinder-volume sum diagnostic
    # ------------------------------------------------------------------
    best = top[0]
    best_kin = _make_kinematics(
        small_limits,
        large_limits,
        _row_values(best),
    )

    total_curves = [
        ("Exact harmonic 250°", harmonic, "--", 1.8),
        (
            f"F5 champion, index {int(f5_champion['index'])}",
            f5_best_kin,
            ":",
            2.0,
        ),
        (
            f"F5 soft-low #9, index {int(f5_soft['index'])}",
            f5_soft_kin,
            "-.",
            2.0,
        ),
        (
            f"F6 #1, index {int(best['index'])}",
            best_kin,
            "-",
            2.2,
        ),
    ]

    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

    for label, kin, linestyle, linewidth in total_curves:
        small, large, dsmall, dlarge = _kinematic_arrays(
            kin,
            angle_rad,
        )
        axes[0].plot(
            angle_deg,
            (small + large) * 1000.0,
            linestyle=linestyle,
            linewidth=linewidth,
            label=label,
        )
        axes[1].plot(
            angle_deg,
            (dsmall + dlarge) * 1000.0,
            linestyle=linestyle,
            linewidth=linewidth,
            label=label,
        )

    axes[0].set_ylabel(r"$V_S + V_L$ (L)")
    axes[1].set_ylabel(r"$d(V_S+V_L)/d\theta$ (L/rad)")
    axes[1].set_xlabel("Study angle (deg)")
    axes[1].axhline(0.0, linewidth=0.8)

    best_values = _row_values(best)
    event_angles = (
        best_values["small_min_deg"],
        LARGE_MIN_DEG,
        best_values["small_max_deg"],
        best_values["large_max_deg"],
    )
    for axis in axes:
        for angle in event_angles:
            axis.axvline(angle, linestyle=":", linewidth=0.8, alpha=0.5)
        axis.set_xlim(0.0, 360.0)
        axis.set_xticks(np.arange(0.0, 361.0, 45.0))
        axis.grid(True, alpha=0.3)
        axis.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.15),
            ncol=2,
            frameon=False,
        )

    fig.suptitle(
        "Stage F6 — cylinder-volume sum diagnostic\n"
        "Kinematic indicator only: not the actual exchanger mass flow",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.95), h_pad=3.0)
    fig.savefig(OUTPUT_TOTAL_PLOT, dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sobol-power",
        type=int,
        default=8,
        help="2**m Sobol points; default m=8 gives 256.",
    )
    parser.add_argument("--seed", type=int, default=2718)
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse matching rows from a compatible F6 CSV.",
    )
    parser.add_argument("--csv", type=Path, default=OUTPUT_CSV)
    parser.add_argument("--output", type=Path, default=OUTPUT_JSON)
    args = parser.parse_args()

    if args.sobol_power < 1:
        raise ValueError("--sobol-power must be positive.")

    definition = CampaignDefinition(A5_CAMPAIGN)
    base_design, total_mass, _saved, _a5_best, candidate_id = (
        _candidate_design_and_mass(definition)
    )
    small_limits = base_design.configuration.machine_volumes.small_cylinder
    large_limits = base_design.configuration.machine_volumes.large_cylinder

    _verify_exact_harmonic(small_limits, large_limits)
    _verify_monotone_flatness_domain()

    f5_champion, f5_soft = _f5_references()

    # No controls in F6: Sobol starts directly at index 0.
    sampler = qmc.Sobol(
        d=len(PARAMETERS),
        scramble=True,
        seed=args.seed,
    )
    candidates = [
        ("sobol", _scale(point))
        for point in sampler.random_base2(args.sobol_power)
    ]

    print("Stage F6: local crossed before/after flatness search")
    print(f"Fixed gas mass: {total_mass:.12g} kg")
    print(f"L_min fixed at {LARGE_MIN_DEG:.1f} deg")
    print(
        "Symmetry: L low/high before = S low/high after; "
        "L low/high after = S low/high before"
    )
    print(
        f"{len(candidates)} evaluations "
        f"({2**args.sobol_power} Sobol, no controls)"
    )
    print(
        "\nF5 champion benchmark: "
        f"eta={100*float(f5_champion['indicated_thermal_efficiency']):.6f}% "
        f"P={float(f5_champion['indicated_power_w']):.3f} W "
        f"(index {int(f5_champion['index'])})"
    )
    print(
        "F5 soft-low #9 reference: "
        f"eta={100*float(f5_soft['indicated_thermal_efficiency']):.6f}% "
        f"P={float(f5_soft['indicated_power_w']):.3f} W "
        f"(index {int(f5_soft['index'])})"
    )
    print()

    checkpoint = _load_checkpoint(args.csv) if args.resume else {}
    rows: list[dict] = []

    for index, (label, values) in enumerate(candidates):
        cached = checkpoint.get(index)
        if cached is not None:
            matches = (
                cached.get("label") == label
                and all(
                    abs(float(cached[name]) - values[name]) < 1e-12
                    for name in PARAMETERS
                )
            )
            if matches:
                row = dict(cached)
                row["index"] = int(float(row["index"]))
                row["feasible"] = (
                    str(row["feasible"]).strip().lower() == "true"
                )
                for key in row:
                    if key not in ("label", "feasible"):
                        try:
                            row[key] = float(row[key])
                        except (TypeError, ValueError):
                            pass
                rows.append(row)

                eta = row.get("indicated_thermal_efficiency")
                eta_text = (
                    "n/a"
                    if eta in (None, "")
                    else f"{100*float(eta):8.5f}%"
                )
                print(f"{index:3d} CK {eta_text:>10s} {label}")
                continue

        kinematics = _make_kinematics(
            small_limits,
            large_limits,
            values,
        )
        design = _same_inventory_design(
            base_design,
            total_mass,
            kinematics,
        )
        result, _ = _evaluate(
            f"two_extrema_f6_{index}",
            design,
            definition,
        )
        row = _row(
            index,
            label,
            values,
            result,
            kinematics,
        )
        rows.append(row)
        _write_csv(args.csv, rows)

        eta = row["indicated_thermal_efficiency"]
        eta_text = "n/a" if eta is None else f"{100*eta:8.5f}%"
        power = row["indicated_power_w"]
        power_text = "n/a" if power is None else f"{power:7.3f} W"
        flag = "OK" if row["feasible"] else "X"

        smin = values["small_min_deg"]
        smax = values["small_max_deg"]
        lmax = values["large_max_deg"]
        sectors = (
            LARGE_MIN_DEG - smin,
            smax - LARGE_MIN_DEG,
            lmax - smax,
            smin + 360.0 - lmax,
        )

        print(
            f"{index:3d} {flag:2s} {eta_text:>10s} "
            f"P={power_text:>10s}  "
            f"Smin={smin:6.2f} "
            f"Smax={smax:6.2f} "
            f"Lmax={lmax:6.2f} "
            f"hSLb={values['h_s_low_before']:6.3f} "
            f"hSLa={values['h_s_low_after']:6.3f} "
            f"hSHb={values['h_s_high_before']:6.3f} "
            f"hSHa={values['h_s_high_after']:6.3f}  "
            f"sectors="
            f"{sectors[0]:.1f}/"
            f"{sectors[1]:.1f}/"
            f"{sectors[2]:.1f}/"
            f"{sectors[3]:.1f}"
        )

    rows.sort(key=lambda row: int(row["index"]))
    _write_csv(args.csv, rows)

    feasible = [
        row
        for row in rows
        if bool(row["feasible"])
        and row.get("indicated_thermal_efficiency") not in (None, "")
        and math.isfinite(float(row["indicated_thermal_efficiency"]))
    ]
    feasible.sort(
        key=lambda row: float(row["indicated_thermal_efficiency"]),
        reverse=True,
    )

    output = {
        "experiment": (
            "Stage F6 local two-extrema study with crossed "
            "before/after endpoint flatness"
        ),
        "comparison_basis": {
            "same_total_working_gas_mass": True,
            "total_mass_kg": total_mass,
            "same_cylinder_volume_limits": True,
            "same_operation": True,
            "same_microtube_hardware": True,
            "mechanical_losses_modeled": False,
            "stage7A5_candidate_id": candidate_id,
        },
        "kinematics": {
            "large_min_deg_fixed": LARGE_MIN_DEG,
            "flatness_symmetry": {
                "h_L_low_before": "h_S_low_after",
                "h_L_low_after": "h_S_low_before",
                "h_L_high_before": "h_S_high_after",
                "h_L_high_after": "h_S_high_before",
                "before_after_direction": "increasing study angle",
            },
            "search_bounds": BOUNDS,
        },
        "diagnostic_note": (
            "V_S + V_L and its angular derivative are kinematic cylinder-volume "
            "diagnostics only. They do not identify the actual hydraulic path or "
            "the actual exchanger mass flow."
        ),
        "f5_references": {
            "champion": f5_champion,
            "soft_low_rank_9_index_40": f5_soft,
        },
        "sobol": {
            "power": args.sobol_power,
            "sample_count": 2 ** args.sobol_power,
            "seed": args.seed,
            "scramble": True,
            "control_count": 0,
        },
        "best_feasible": feasible[0] if feasible else None,
        "ranking": feasible,
        "all_rows": rows,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n")

    _plot_results(
        small_limits,
        large_limits,
        rows,
        f5_champion,
        f5_soft,
    )

    print("\nTOP FEASIBLE")
    for rank, row in enumerate(feasible[:12], 1):
        print(
            f"{rank:2d} "
            f"eta={100*float(row['indicated_thermal_efficiency']):.6f}% "
            f"P={float(row['indicated_power_w']):.3f} W  "
            f"Smin={float(row['small_min_deg']):.3f} "
            f"Smax={float(row['small_max_deg']):.3f} "
            f"Lmax={float(row['large_max_deg']):.3f} "
            f"hSLb={float(row['h_s_low_before']):.4f} "
            f"hSLa={float(row['h_s_low_after']):.4f} "
            f"hSHb={float(row['h_s_high_before']):.4f} "
            f"hSHa={float(row['h_s_high_after']):.4f}"
        )

    print(f"\nWrote {args.csv.relative_to(ROOT)}")
    print(f"Wrote {args.output.relative_to(ROOT)}")
    if feasible:
        print(f"Wrote {OUTPUT_PLOT.relative_to(ROOT)}")
        print(f"Wrote {OUTPUT_TOTAL_PLOT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
