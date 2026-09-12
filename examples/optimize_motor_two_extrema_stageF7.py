"""Stage F7: dense local 11D refinement around the K2/F6 basin.

The broad F6-11D pass used 512 Sobol points in 11 dimensions and did not beat
the K2/F6 reference, even though the 11D family contains that exact motion.
F7 therefore keeps the same 11-variable family but narrows the search box and
samples it much more densely.

Frozen:
- Stage K2 champion exchanger asymmetry;
- Stage F6 cylinder-volume limits, total gas mass and 2 Hz operation;
- large-cylinder minimum at 180 deg.

Free (11 variables):
- small_min_deg
- small_max_deg
- large_max_deg
- four one-sided h values for S
- four one-sided h values for L

The exact K2/F6 motion and the top three broad F6-11D points are evaluated as
explicit controls. Default local search: 1024 scrambled Sobol points.

Outputs:
    outputs/motor_two_extrema_stageF7.csv
    outputs/motor_two_extrema_stageF7.json

Run:
    PYTHONPATH=src python3 examples/optimize_motor_two_extrema_stageF7.py

Resume:
    PYTHONPATH=src python3 examples/optimize_motor_two_extrema_stageF7.py --resume
"""

from __future__ import annotations

from dataclasses import replace
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import qmc

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.exchangers.microtube import MicrotubeExchanger

from compare_motor_motion_laws_stage7A5 import (
    ROOT,
    A5_CAMPAIGN,
    _candidate_design_and_mass,
    _evaluate,
    _same_inventory_design,
)
from optimize_motor_exchanger_asymmetry_stageK1 import _scaled_exchanger
from optimize_motor_two_extrema_stageF5 import LARGE_MIN_DEG
from optimize_motor_two_extrema_stageF6_11d import (
    PARAMETERS,
    _independent_reference_values,
    _load_f6_champion,
    _load_k2_champion,
    _make_kinematics,
    _row as _base_row,
)

F6_11D_JSON = ROOT / "outputs" / "motor_two_extrema_stageF6_11d.json"
OUTPUT_CSV = ROOT / "outputs" / "motor_two_extrema_stageF7.csv"
OUTPUT_JSON = ROOT / "outputs" / "motor_two_extrema_stageF7.json"

# Local box containing:
# - the exact crossed F6/K2 reference;
# - F6-11D ranks 1, 2 and 3;
# - margin around those points.
BOUNDS = {
    "small_min_deg": (67.0, 77.5),
    "small_max_deg": (230.0, 239.0),
    "large_max_deg": (362.0, 371.0),

    "h_s_low_before": (-0.35, 0.35),
    "h_s_low_after": (-0.55, 0.05),
    "h_s_high_before": (-1.80, -1.20),
    "h_s_high_after": (-1.10, -0.35),

    "h_l_low_before": (-0.40, 0.20),
    "h_l_low_after": (-0.60, 0.35),
    "h_l_high_before": (-1.05, -0.30),
    "h_l_high_after": (-1.70, -0.80),
}


def _load_first_pass():
    data = json.loads(F6_11D_JSON.read_text())
    ranking = data.get("ranking", [])
    if len(ranking) < 3:
        raise RuntimeError("F6-11D must contain at least three feasible ranks.")
    return data


def _values(row):
    return {name: float(row[name]) for name in PARAMETERS}


def _scale(point):
    out = {}
    for x, name in zip(point, PARAMETERS, strict=True):
        lo, hi = BOUNDS[name]
        out[name] = lo + float(x) * (hi - lo)
    return out


def _verify_controls(controls):
    errors = []
    for label, values in controls:
        for name in PARAMETERS:
            lo, hi = BOUNDS[name]
            if not lo <= values[name] <= hi:
                errors.append(
                    f"{label}: {name}={values[name]} outside [{lo}, {hi}]"
                )
    if errors:
        raise RuntimeError(
            "F7 bounds do not contain mandatory controls:\n"
            + "\n".join(errors)
        )


def _verify_monotonicity():
    y = np.linspace(0.0, 1.0, 20001)
    branch_pairs = (
        ("h_s_low_after", "h_s_high_before"),
        ("h_s_high_after", "h_s_low_before"),
        ("h_l_low_after", "h_l_high_before"),
        ("h_l_high_after", "h_l_low_before"),
    )
    worst = math.inf
    worst_case = None
    for start_name, end_name in branch_pairs:
        for h0 in BOUNDS[start_name]:
            for h1 in BOUNDS[end_name]:
                derivative = (
                    1.0
                    - h0 * (1.0 - 4.0*y + 3.0*y*y)
                    + h1 * (2.0*y - 3.0*y*y)
                )
                minimum = float(np.min(derivative))
                if minimum < worst:
                    worst = minimum
                    worst_case = (start_name, h0, end_name, h1)
    if worst < -1e-12:
        raise RuntimeError(
            f"F7 bounds can create a non-monotone branch: "
            f"{worst_case}, min dE/dy={worst}"
        )
    print(f"Monotonicity box check passed: min dE/dy={worst:.6f}")


def _edge_distance(values):
    distances = []
    for name in PARAMETERS:
        lo, hi = BOUNDS[name]
        u = (values[name] - lo) / (hi - lo)
        distances.append(min(u, 1.0-u))
    return float(min(distances))


def _row(index, label, values, result, kinematics):
    row = _base_row(index, values, result, kinematics)
    row["label"] = label
    row["minimum_normalized_bound_distance"] = _edge_distance(values)
    return row


def _write_csv(path, rows):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _load_checkpoint(path):
    if not path.exists():
        return {}
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        return {}
    required = set(PARAMETERS) | {"index", "label"}
    if not required.issubset(rows[0]):
        return {}
    return {int(float(row["index"])): row for row in rows}


def _restore(row):
    out = dict(row)
    out["index"] = int(float(out["index"]))
    for key in list(out):
        if key == "label":
            continue
        if key == "feasible" or key.startswith("constraint_"):
            out[key] = str(out[key]).strip().lower() == "true"
            continue
        if out[key] in ("", None):
            out[key] = None
            continue
        try:
            out[key] = float(out[key])
        except (TypeError, ValueError):
            pass
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sobol-power",
        type=int,
        default=10,
        help="2**m local Sobol points; default 10 = 1024.",
    )
    parser.add_argument("--seed", type=int, default=7711)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--csv", type=Path, default=OUTPUT_CSV)
    parser.add_argument("--output", type=Path, default=OUTPUT_JSON)
    args = parser.parse_args()

    if args.sobol_power < 1:
        raise ValueError("--sobol-power must be positive.")

    definition = CampaignDefinition(A5_CAMPAIGN)
    base_design, total_mass, _saved, _a5_best, candidate_id = (
        _candidate_design_and_mass(definition)
    )
    if not isinstance(base_design.heat_in, MicrotubeExchanger):
        raise TypeError("Expected microtube H_i design.")
    if not isinstance(base_design.heat_out, MicrotubeExchanger):
        raise TypeError("Expected microtube H_o design.")

    f6_best, f6_values = _load_f6_champion()
    k2_best = _load_k2_champion()
    first_pass = _load_first_pass()
    first_ranking = first_pass["ranking"]

    crossed_reference = _independent_reference_values(f6_values)

    controls = [
        ("k2_f6_crossed_reference", crossed_reference),
        ("f6_11d_rank1", _values(first_ranking[0])),
        ("f6_11d_rank2", _values(first_ranking[1])),
        ("f6_11d_rank3", _values(first_ranking[2])),
    ]
    _verify_controls(controls)
    _verify_monotonicity()

    small_limits = base_design.configuration.machine_volumes.small_cylinder
    large_limits = base_design.configuration.machine_volumes.large_cylinder

    reference_design = _same_inventory_design(
        base_design,
        total_mass,
        _make_kinematics(
            small_limits,
            large_limits,
            crossed_reference,
        ),
    )
    heat_in = _scaled_exchanger(
        reference_design.heat_in,
        float(k2_best["k_i"]),
    )
    heat_out = _scaled_exchanger(
        reference_design.heat_out,
        float(k2_best["k_o"]),
    )

    sampler = qmc.Sobol(
        d=len(PARAMETERS),
        scramble=True,
        seed=args.seed,
    )
    sobol = [
        _scale(point)
        for point in sampler.random_base2(args.sobol_power)
    ]

    candidates = list(controls)
    candidates.extend(
        (f"sobol_{i}", values)
        for i, values in enumerate(sobol)
    )

    print("Stage F7 — dense local 11D refinement")
    print(
        f"K2/F6 reference: "
        f"eta={100*float(k2_best['indicated_thermal_efficiency']):.6f}% "
        f"P={float(k2_best['indicated_power_w']):.3f} W"
    )
    print(
        f"Frozen exchanger asymmetry: "
        f"k_i={float(k2_best['k_i']):.6f}, "
        f"k_o={float(k2_best['k_o']):.6f}"
    )
    print(
        f"First-pass 11D best: "
        f"eta={100*float(first_ranking[0]['indicated_thermal_efficiency']):.6f}% "
        f"P={float(first_ranking[0]['indicated_power_w']):.3f} W"
    )
    print(f"Fixed total gas mass: {total_mass:.12g} kg")
    print(
        f"{len(controls)} controls + {len(sobol)} local Sobol points"
    )
    print()

    checkpoint = _load_checkpoint(args.csv) if args.resume else {}
    rows = []

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
                row = _restore(cached)
                rows.append(row)
                eta = row.get("indicated_thermal_efficiency")
                eta_text = "n/a" if eta is None else f"{100*eta:8.5f}%"
                print(f"{index:4d} CK {eta_text:>10s} {label}")
                continue

        kinematics = _make_kinematics(
            small_limits,
            large_limits,
            values,
        )
        design = replace(
            _same_inventory_design(
                base_design,
                total_mass,
                kinematics,
            ),
            heat_in=heat_in,
            heat_out=heat_out,
        )

        result, _last_state = _evaluate(
            f"two_extrema_f7_{index}",
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
        eta_text = "n/a" if eta is None else f"{100*float(eta):8.5f}%"
        power = row["indicated_power_w"]
        power_text = "n/a" if power is None else f"{float(power):7.3f} W"
        flag = "OK" if row["feasible"] else "X"

        print(
            f"{index:4d} {flag:2s} {eta_text:>10s} "
            f"P={power_text:>10s} "
            f"Smin={values['small_min_deg']:6.2f} "
            f"Smax={values['small_max_deg']:6.2f} "
            f"Lmax={values['large_max_deg']:6.2f} "
            f"symRMS={float(row['cross_symmetry_rms_h']):.3f} "
            f"{label}"
        )

        if index == 0 and eta is not None:
            delta_pp = 100.0 * (
                float(eta)
                - float(k2_best["indicated_thermal_efficiency"])
            )
            print(
                f"     K2 reference re-evaluation delta: "
                f"{delta_pp:+.6f} percentage point"
            )

    rows.sort(key=lambda row: int(row["index"]))
    _write_csv(args.csv, rows)

    feasible = [
        row for row in rows
        if bool(row["feasible"])
        and row.get("indicated_thermal_efficiency") is not None
        and math.isfinite(float(row["indicated_thermal_efficiency"]))
    ]
    feasible.sort(
        key=lambda row: float(row["indicated_thermal_efficiency"]),
        reverse=True,
    )

    reference_row = next(
        (
            row for row in rows
            if row["label"] == "k2_f6_crossed_reference"
        ),
        None,
    )

    reference_eta = (
        float(reference_row["indicated_thermal_efficiency"])
        if reference_row is not None
        and reference_row.get("indicated_thermal_efficiency") is not None
        else float(k2_best["indicated_thermal_efficiency"])
    )

    best = feasible[0] if feasible else None
    improvement_pp = (
        100.0 * (
            float(best["indicated_thermal_efficiency"])
            - reference_eta
        )
        if best is not None
        else None
    )

    report = {
        "experiment": "Stage F7 dense local 11D refinement",
        "comparison_basis": {
            "same_total_working_gas_mass": True,
            "total_mass_kg": total_mass,
            "same_cylinder_volume_limits": True,
            "frequency_hz": 2.0,
            "K2_exchanger_geometry_frozen": True,
            "K2_k_i": float(k2_best["k_i"]),
            "K2_k_o": float(k2_best["k_o"]),
            "large_min_deg_fixed": LARGE_MIN_DEG,
            "crossed_flatness_symmetry_imposed": False,
            "mechanical_losses_modeled": False,
            "external_aerodynamic_losses_in_objective": False,
            "stage7A5_candidate_id": candidate_id,
        },
        "search": {
            "parameter_count": len(PARAMETERS),
            "parameters": PARAMETERS,
            "bounds": BOUNDS,
            "sobol_power": args.sobol_power,
            "sobol_count": len(sobol),
            "explicit_control_count": len(controls),
            "seed": args.seed,
            "scramble": True,
        },
        "references": {
            "F6_crossed_champion": f6_best,
            "K2_champion_with_F6_motion": k2_best,
            "first_pass_11D_top3": first_ranking[:3],
            "K2_F6_reference_11D_coordinates": crossed_reference,
        },
        "reference_reevaluation": reference_row,
        "best_feasible": best,
        "improvement_over_reference_percentage_points": improvement_pp,
        "ranking": feasible,
        "all_rows": rows,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")

    print("\nTOP FEASIBLE")
    for rank, row in enumerate(feasible[:20], 1):
        print(
            f"{rank:2d} "
            f"eta={100*float(row['indicated_thermal_efficiency']):.6f}% "
            f"P={float(row['indicated_power_w']):.3f} W "
            f"Smin={float(row['small_min_deg']):.3f} "
            f"Smax={float(row['small_max_deg']):.3f} "
            f"Lmax={float(row['large_max_deg']):.3f} "
            f"symRMS={float(row['cross_symmetry_rms_h']):.4f} "
            f"edge={float(row['minimum_normalized_bound_distance']):.3f} "
            f"{row['label']}"
        )

    if best is not None:
        print("\nBEST VS K2/F6 REFERENCE")
        print(f"  delta eta = {improvement_pp:+.6f} percentage point")
        print(
            f"  delta P   = "
            f"{float(best['indicated_power_w']) - float(k2_best['indicated_power_w']):+.6f} W"
        )
        print(
            f"  timings   = "
            f"Smin {float(best['small_min_deg']):.4f}, "
            f"Smax {float(best['small_max_deg']):.4f}, "
            f"Lmax {float(best['large_max_deg']):.4f}"
        )
        print(
            f"  S low     = "
            f"{float(best['h_s_low_before']):+.5f} / "
            f"{float(best['h_s_low_after']):+.5f}"
        )
        print(
            f"  S high    = "
            f"{float(best['h_s_high_before']):+.5f} / "
            f"{float(best['h_s_high_after']):+.5f}"
        )
        print(
            f"  L low     = "
            f"{float(best['h_l_low_before']):+.5f} / "
            f"{float(best['h_l_low_after']):+.5f}"
        )
        print(
            f"  L high    = "
            f"{float(best['h_l_high_before']):+.5f} / "
            f"{float(best['h_l_high_after']):+.5f}"
        )

    print(f"\nWrote {args.csv.relative_to(ROOT)}")
    print(f"Wrote {args.output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
