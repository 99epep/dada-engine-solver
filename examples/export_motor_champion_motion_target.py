"""Export the current thermodynamic champion motion as a 4-bar synthesis target.

The current champion is the K2 exchanger optimum evaluated with the crossed-F6
motion law. This script exports that motion independently of exchanger details
so a mechanical synthesis tool can fit rocker-point or coupler-point four-bar
outputs to it.

Outputs
-------
    outputs/motor_champion_motion_target.csv
    outputs/motor_champion_motion_target.json

Run
---
    PYTHONPATH=src python3 examples/export_motor_champion_motion_target.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from dada_solver.campaign.definition import CampaignDefinition

from compare_motor_motion_laws_stage7A5 import (
    ROOT,
    A5_CAMPAIGN,
    _candidate_design_and_mass,
)
from optimize_motor_two_extrema_stageF6_11d import (
    _independent_reference_values,
    _load_f6_champion,
    _load_k2_champion,
    _make_kinematics,
)

OUTPUT_CSV = ROOT / "outputs" / "motor_champion_motion_target.csv"
OUTPUT_JSON = ROOT / "outputs" / "motor_champion_motion_target.json"


def _fourier_coefficients(values, theta, maximum_harmonic):
    n_samples = len(values)
    result = {
        "a0": float(np.mean(values)),
        "harmonics": [],
    }

    for harmonic in range(1, maximum_harmonic + 1):
        a_n = (
            2.0 / n_samples
            * float(np.sum(values * np.cos(harmonic * theta)))
        )
        b_n = (
            2.0 / n_samples
            * float(np.sum(values * np.sin(harmonic * theta)))
        )
        result["harmonics"].append(
            {
                "n": harmonic,
                "a_cos": a_n,
                "b_sin": b_n,
                "amplitude": math.hypot(a_n, b_n),
                "phase_deg_cosine_convention": math.degrees(
                    math.atan2(-b_n, a_n)
                ),
            }
        )

    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--step-deg", type=float, default=0.25)
    parser.add_argument("--fourier-harmonics", type=int, default=12)
    parser.add_argument("--csv", type=Path, default=OUTPUT_CSV)
    parser.add_argument("--output", type=Path, default=OUTPUT_JSON)
    args = parser.parse_args()

    if not math.isfinite(args.step_deg) or args.step_deg <= 0.0:
        raise ValueError("--step-deg must be finite and positive.")
    if args.fourier_harmonics < 1:
        raise ValueError("--fourier-harmonics must be positive.")

    definition = CampaignDefinition(A5_CAMPAIGN)
    base_design, total_mass, _saved, _a5_best, candidate_id = (
        _candidate_design_and_mass(definition)
    )

    f6_best, f6_values = _load_f6_champion()
    k2_best = _load_k2_champion()
    motion_values = _independent_reference_values(f6_values)

    small_limits = base_design.configuration.machine_volumes.small_cylinder
    large_limits = base_design.configuration.machine_volumes.large_cylinder

    kinematics = _make_kinematics(
        small_limits,
        large_limits,
        motion_values,
    )

    sample_count = int(round(360.0 / args.step_deg))
    if not math.isclose(
        sample_count * args.step_deg,
        360.0,
        rel_tol=0.0,
        abs_tol=1e-10,
    ):
        raise ValueError("--step-deg must divide 360 deg exactly.")

    angle_deg = np.linspace(0.0, 360.0, sample_count + 1)
    angle_rad = np.radians(angle_deg)

    rows = []
    for deg, theta in zip(angle_deg, angle_rad, strict=True):
        t = float(theta)

        v_s = float(kinematics.small_cylinder_volume(t))
        v_l = float(kinematics.large_cylinder_volume(t))
        dv_s = float(kinematics.small_cylinder_volume_derivative(t))
        dv_l = float(kinematics.large_cylinder_volume_derivative(t))
        d2v_s = float(
            kinematics.small_cylinder_volume_second_derivative(t)
        )
        d2v_l = float(
            kinematics.large_cylinder_volume_second_derivative(t)
        )

        q_s = (v_s - small_limits.minimum) / small_limits.swept
        q_l = (v_l - large_limits.minimum) / large_limits.swept
        dq_s = dv_s / small_limits.swept
        dq_l = dv_l / large_limits.swept
        d2q_s = d2v_s / small_limits.swept
        d2q_l = d2v_l / large_limits.swept

        rows.append(
            {
                "theta_deg": float(deg),
                "theta_rad": t,
                "small_volume_m3": v_s,
                "large_volume_m3": v_l,
                "small_volume_l": 1000.0 * v_s,
                "large_volume_l": 1000.0 * v_l,
                "small_fraction_0_1": q_s,
                "large_fraction_0_1": q_l,
                "small_centered_minus1_plus1": 2.0 * q_s - 1.0,
                "large_centered_minus1_plus1": 2.0 * q_l - 1.0,
                "small_dq_dtheta_per_rad": dq_s,
                "large_dq_dtheta_per_rad": dq_l,
                "small_d2q_dtheta2_per_rad2": d2q_s,
                "large_d2q_dtheta2_per_rad2": d2q_l,
                "small_centered_dq_dtheta_per_rad": 2.0 * dq_s,
                "large_centered_dq_dtheta_per_rad": 2.0 * dq_l,
                "small_centered_d2q_dtheta2_per_rad2": 2.0 * d2q_s,
                "large_centered_d2q_dtheta2_per_rad2": 2.0 * d2q_l,
                "total_cylinder_volume_l": 1000.0 * (v_s + v_l),
                "dtotal_volume_dtheta_l_per_rad": 1000.0 * (dv_s + dv_l),
            }
        )

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fourier_count = max(
        4096,
        2 ** int(math.ceil(math.log2(sample_count))),
    )
    theta_f = np.linspace(
        0.0,
        2.0 * math.pi,
        fourier_count,
        endpoint=False,
    )

    q_s_f = np.asarray(
        [
            (
                kinematics.small_cylinder_volume(float(t))
                - small_limits.minimum
            )
            / small_limits.swept
            for t in theta_f
        ],
        dtype=float,
    )
    q_l_f = np.asarray(
        [
            (
                kinematics.large_cylinder_volume(float(t))
                - large_limits.minimum
            )
            / large_limits.swept
            for t in theta_f
        ],
        dtype=float,
    )

    x_s_f = 2.0 * q_s_f - 1.0
    x_l_f = 2.0 * q_l_f - 1.0

    report = {
        "purpose": (
            "Mechanism-synthesis target for the current thermodynamic champion."
        ),
        "angle_convention": (
            "Study angle before motor-operation direction transformation; "
            "physical crank direction/sign may be chosen later."
        ),
        "source": {
            "stage7A5_candidate_id": candidate_id,
            "total_working_gas_mass_kg": total_mass,
            "F6_crossed_motion": f6_best,
            "K2_exchanger_champion": {
                "k_i": float(k2_best["k_i"]),
                "k_o": float(k2_best["k_o"]),
                "indicated_thermal_efficiency": float(
                    k2_best["indicated_thermal_efficiency"]
                ),
                "indicated_power_w": float(k2_best["indicated_power_w"]),
            },
            "motion_parameters_11D_coordinates": motion_values,
        },
        "cylinder_limits": {
            "small": {
                "minimum_m3": small_limits.minimum,
                "maximum_m3": small_limits.maximum,
                "swept_m3": small_limits.swept,
            },
            "large": {
                "minimum_m3": large_limits.minimum,
                "maximum_m3": large_limits.maximum,
                "swept_m3": large_limits.swept,
            },
        },
        "csv_sampling": {
            "step_deg": args.step_deg,
            "rows_including_360_deg": len(rows),
        },
        "recommended_mechanism_fit_coordinates": {
            "small": "small_centered_minus1_plus1",
            "large": "large_centered_minus1_plus1",
            "derivative_small": "small_centered_dq_dtheta_per_rad",
            "derivative_large": "large_centered_dq_dtheta_per_rad",
            "second_derivative_small": "small_centered_d2q_dtheta2_per_rad2",
            "second_derivative_large": "large_centered_d2q_dtheta2_per_rad2",
            "note": (
                "Fit normalized shape first; absolute piston stroke and area "
                "can be chosen during physical mechanism sizing."
            ),
        },
        "fourier": {
            "convention": (
                "f(theta)=a0+sum[a_n*cos(n*theta)+b_n*sin(n*theta)]"
            ),
            "sample_count": fourier_count,
            "maximum_harmonic": args.fourier_harmonics,
            "small_centered_motion": _fourier_coefficients(
                x_s_f,
                theta_f,
                args.fourier_harmonics,
            ),
            "large_centered_motion": _fourier_coefficients(
                x_l_f,
                theta_f,
                args.fourier_harmonics,
            ),
        },
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")

    print(
        "Current champion motion target\n"
        f"  eta = {100*float(k2_best['indicated_thermal_efficiency']):.6f}%\n"
        f"  P   = {float(k2_best['indicated_power_w']):.6f} W\n"
        f"  F6 timings: "
        f"Smin={motion_values['small_min_deg']:.6f} deg, "
        f"Smax={motion_values['small_max_deg']:.6f} deg, "
        f"Lmin=180.000000 deg, "
        f"Lmax={motion_values['large_max_deg']:.6f} deg"
    )
    print(f"Wrote {args.csv.relative_to(ROOT)}")
    print(f"Wrote {args.output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
