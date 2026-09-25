#!/usr/bin/env python3
"""Export the final rank_01 thermo-5D six-bar champion cycle for animation.

Outputs:
  outputs/rank_01_thermodynamic_cycle.csv
  outputs/rank_01_thermodynamic_cycle_raw.csv
  outputs/rank_01_thermodynamic_cycle_metadata.json
  outputs/rank_01_animation_model.json

Run from the repository root:

  PYTHONPATH=src python3 examples/export_rank01_thermodynamic_cycle.py
"""

from __future__ import annotations

from dataclasses import asdict, replace
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math

import numpy as np

from dada_solver.dynamics import ValveTopology
from dada_solver.exchangers.wall_cycle import (
    solve_periodic_wall_motor,
    wall_cycle_performance,
)
from dada_solver.heat_transfer import PrescribedHeatRate
from dada_solver.state import ThermodynamicState
from dada_solver.valves import ValveState

import optimize_sixbar_pairs_thermo5d_3952_hlat25 as thermo5d


ROOT = Path.cwd()

DEFAULT_REPORT = (
    ROOT
    / "outputs"
    / "sixbar_thermo5d_3952"
    / "rank_01"
    / "report.json"
)

DEFAULT_PREFIX = ROOT / "outputs" / "rank_01_thermodynamic_cycle"
DEFAULT_MODEL = ROOT / "outputs" / "rank_01_animation_model.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def mechanism_for_animation(mech) -> dict:
    """Return the field names expected by animate_sixbar_graphs.py."""
    return {
        "primary_ground": float(mech.primary_ground),
        "primary_coupler": float(mech.primary_coupler),
        "primary_rocker": float(mech.primary_rocker),
        "primary_e_along": float(mech.primary_e_along),
        "primary_e_normal": float(mech.primary_e_normal),
        "primary_phase": float(mech.primary_phase),
        "second_pivot_x": float(mech.second_pivot_x),
        "second_pivot_y": float(mech.second_pivot_y),
        "link_ef": float(mech.link_ef),
        "link_gf": float(mech.link_gf),
        "h_along_over_ef": float(mech.h_along_over_ef),
        "h_normal_over_ef": float(mech.h_normal_over_ef),
        "piston_rod": float(mech.piston_rod),
        "slider_axis_offset": float(mech.slider_axis_offset),
        "slider_axis_angle": float(mech.slider_axis_angle),
        "primary_branch": int(mech.primary_branch),
        "second_branch": int(mech.second_branch),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    ap.add_argument(
        "--output-prefix",
        type=Path,
        default=DEFAULT_PREFIX,
        help="Prefix without extension for CSV/metadata output.",
    )
    ap.add_argument(
        "--model-output",
        type=Path,
        default=DEFAULT_MODEL,
        help="Animation-compatible mechanism JSON.",
    )
    ap.add_argument(
        "--samples",
        type=int,
        default=1441,
        help="Number of samples in the regular-grid CSV, including both endpoints.",
    )
    args = ap.parse_args()

    if args.samples < 3:
        raise ValueError("--samples must be >= 3")

    report_path = args.report.resolve()
    report = load_json(report_path)

    if int(report.get("family_rank", -1)) != 1:
        raise ValueError(
            f"Expected rank_01 report, got family_rank={report.get('family_rank')!r}"
        )

    best = report.get("best_feasible")
    if not isinstance(best, dict):
        raise RuntimeError("rank_01 report has no best_feasible candidate.")

    if best.get("result", {}).get("status") != "converged":
        raise RuntimeError("rank_01 best_feasible is not converged.")

    params = best.get("parameters")
    if not isinstance(params, dict):
        raise RuntimeError("rank_01 best_feasible has no thermo-5D parameters.")

    initial_state = best.get("last_complete_state")
    if initial_state is None:
        raise RuntimeError("rank_01 best_feasible has no last_complete_state.")

    # Rebuild exactly the same 3952 basis and the exact fixed family-1 pair
    # used by the thermo-5D campaign.
    (
        seed_machine,
        definition,
        total_mass,
        geometry_info,
        _thermo_seed,
        basis_identity,
    ) = thermo5d.make_3952_basis()

    source = thermo5d.load_pair(1)

    expected_pair_id = report.get("fixed_pair", {}).get("candidate_id")
    if expected_pair_id and source["pair_id"] != expected_pair_id:
        raise RuntimeError(
            "Family-1 pair mismatch between thermo5D report and source pair:\n"
            f"  report: {expected_pair_id}\n"
            f"  source: {source['pair_id']}"
        )

    design = thermo5d.build_design(
        seed_machine,
        geometry_info,
        total_mass,
        params,
        source["small"],
        source["large"],
    )

    wrapper = design.build()

    periodic = solve_periodic_wall_motor(
        wrapper,
        np.asarray(initial_state, dtype=float),
        maximum_cycles=design.configuration.numerical.maximum_cycles,
        settings=definition.wall_numerical_settings,
        backend=definition.wall_backend,
        exact_kinematics_cache=definition.exact_kinematics_cache,
    )

    if not periodic.converged:
        raise RuntimeError(
            "Export requires a converged physical cycle; "
            f"solver returned {periodic.status}: {periodic.message}"
        )

    if periodic.angles is None or periodic.trajectory is None:
        raise RuntimeError("Periodic solve returned no trajectory.")

    model = wrapper.model
    stored = ValveTopology(ValveState.CLOSED, ValveState.CLOSED)

    def row(angle: float, values: np.ndarray) -> dict:
        gas = ThermodynamicState.from_array(values[:8])

        volumes = model.volumes(angle)
        temperatures, pressures = gas.temperatures_and_pressures(
            model.gas,
            volumes,
        )

        # Use the current contextual exchanger path rather than calling
        # exchanger.rates() directly. This keeps the export consistent with the
        # microtube gas-side model used by the thermo-5D optimization.
        hi, ho = wrapper.thermal_rates(angle, values)

        instantaneous = replace(
            model,
            cold_heat_transfer=PrescribedHeatRate(hi["gas_heat_w"]),
            hot_heat_transfer=PrescribedHeatRate(ho["gas_heat_w"]),
        )

        rates = instantaneous.evaluate(angle, gas, stored)
        topology = model.effective_topology(angle, gas, stored)

        # For these final IndependentSixBarVolumeKinematics mechanisms, the
        # mechanism phases are already expressed in the solver's study-angle
        # convention and are driven directly by the solver angle. Therefore:
        # study_angle_deg == motor_angle_deg.
        out = {
            "time_s": angle / model.angular_speed,
            "cycle_fraction": angle / (2 * math.pi),
            "motor_angle_deg": math.degrees(angle),
            "study_angle_deg": -math.degrees(angle),
        }

        for i, region in enumerate(("S", "L", "Hi", "Ho")):
            out.update(
                {
                    f"{region}_volume_m3": volumes.as_tuple()[i],
                    f"{region}_pressure_Pa": pressures[i],
                    f"{region}_temperature_K": temperatures[i],
                    f"{region}_mass_kg": values[2 * i],
                    f"{region}_internal_energy_J": values[2 * i + 1],
                }
            )

        for region, index, exchanger, thermal in (
            ("Hi", 8, wrapper.heat_in, hi),
            ("Ho", 9, wrapper.heat_out, ho),
        ):
            out.update(
                {
                    f"{region}_wall_temperature_K":
                        values[index] / exchanger.wall_capacity_j_k,
                    f"{region}_wall_energy_J": values[index],
                    f"{region}_wall_to_gas_heat_W": thermal["gas_heat_w"],
                    f"{region}_external_to_wall_heat_W":
                        thermal["air_heat_w"],
                    f"{region}_air_inlet_temperature_K":
                        exchanger.air_inlet_temperature_k,
                    f"{region}_air_outlet_temperature_K":
                        thermal["air_outlet_temperature_k"],
                }
            )

        for name, field in (
            ("S_to_Hi", "small_to_cold"),
            ("Hi_to_L", "cold_to_large"),
            ("L_to_Ho", "large_to_hot"),
            ("Ho_to_S", "hot_to_small"),
        ):
            out[f"{name}_mass_flow_kg_s"] = getattr(rates.flows, field)

        out["Hi_to_L_valve_open"] = int(
            topology.cold_to_large == ValveState.OPEN
        )
        out["Ho_to_S_valve_open"] = int(
            topology.hot_to_small == ValveState.OPEN
        )

        vs, vl = model.cylinder_volume_rates(angle)
        out.update(
            {
                "S_volume_rate_m3_s": vs,
                "L_volume_rate_m3_s": vl,
                "S_indicated_power_W": pressures[0] * vs,
                "L_indicated_power_W": pressures[1] * vl,
                "total_indicated_power_W": rates.gas_work_rate,
                "cumulative_external_Hi_heat_J": values[10],
                "cumulative_external_Ho_heat_J": values[11],
                "cumulative_Hi_gas_heat_J": values[12],
                "cumulative_Ho_gas_heat_J": values[13],
                "cumulative_indicated_work_J": values[14],
            }
        )

        return out

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)

    raw_path = prefix.with_name(prefix.name + "_raw").with_suffix(".csv")
    main_path = prefix.with_suffix(".csv")
    metadata_path = prefix.with_name(prefix.name + "_metadata").with_suffix(".json")

    files = []

    regular_angles = np.linspace(
        0.0,
        2 * math.pi,
        args.samples,
        endpoint=True,
    )
    regular_trajectory = np.array(
        [
            np.interp(regular_angles, periodic.angles, values)
            for values in periodic.trajectory
        ]
    )

    for path, angles, trajectory, kind in (
        (
            raw_path,
            periodic.angles,
            periodic.trajectory,
            "raw adaptive solver grid",
        ),
        (
            main_path,
            regular_angles,
            regular_trajectory,
            f"regular {args.samples}-point grid",
        ),
    ):
        rows = [
            row(float(angle), values)
            for angle, values in zip(angles, trajectory.T)
        ]

        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

        files.append(
            {
                "file": path.name,
                "rows": len(rows),
                "columns": len(rows[0]),
                "sampling": kind,
            }
        )

    performance = wall_cycle_performance(wrapper, periodic.trajectory)

    metadata = {
        "study": "final rank_01 six-bar thermo-5D champion cycle export",
        "family_rank": 1,
        "thermo5d_report": str(report_path),
        "thermo5d_report_sha256":
            hashlib.sha256(report_path.read_bytes()).hexdigest(),
        "thermo5d_candidate_id": best.get("candidate_id"),
        "pair_candidate_id": source["pair_id"],
        "thermo5d_parameters": {
            "swept_ratio": float(params["swept_ratio"]),
            "n_i": int(params["n_i"]),
            "length_i_m": float(params["length_i_m"]),
            "n_o": int(params["n_o"]),
            "length_o_m": float(params["length_o_m"]),
        },
        "fixed_total_gas_mass_kg": total_mass,
        "basis": basis_identity,
        "files": files,
        "indicated_thermal_efficiency":
            performance.thermal_efficiency,
        "indicated_power_W": performance.gas_power,
        "reported_champion_efficiency":
            best["result"]["indicated_thermal_efficiency"],
        "reported_champion_power_W":
            best["result"]["indicated_power_w"],
        "convergence_history": periodic.history,
        "period_s": 2 * math.pi / model.angular_speed,
        "sampling":
            "Main CSV: regular angular grid, conservative states and "
            "quadratures linearly interpolated from the converged adaptive "
            "cycle; instantaneous rates recomputed. Raw CSV: adaptive solver "
            "samples.",
        "signs":
            "Heat positive into gas/wall as named; indicated work positive "
            "on expansion; signed flows positive in column-name direction.",
        "angle_convention":
            "For final IndependentSixBarVolumeKinematics, study angle equals "
            "solver/motor angle. No extra sign reversal is applied.",
        "units":
            "SI units encoded in column names. Valve open: 1, closed: 0.",
        "endpoint":
            "Both cycle endpoints retained, without forcing the final state "
            "equal to the initial state.",
        "mechanical_losses": "unknown",
    }

    metadata_path.write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )

    model_output = args.model_output
    model_output.parent.mkdir(parents=True, exist_ok=True)
    animation_model = {
        "experiment":
            "Final family-1 six-bar pair used by the rank_01 thermo-5D champion",
        "family_rank": 1,
        "thermo5d_candidate_id": best.get("candidate_id"),
        "pair_candidate_id": source["pair_id"],
        "source_pair": str(source["pair_path"]),
        "kinematics": {
            "small": mechanism_for_animation(source["small"]),
            "large": mechanism_for_animation(source["large"]),
        },
    }
    model_output.write_text(
        json.dumps(animation_model, indent=2) + "\n",
        encoding="utf-8",
    )

    print(json.dumps(
        {
            "csv": str(main_path),
            "raw_csv": str(raw_path),
            "metadata": str(metadata_path),
            "animation_model": str(model_output),
            "efficiency": performance.thermal_efficiency,
            "power_W": performance.gas_power,
            "reported_efficiency":
                best["result"]["indicated_thermal_efficiency"],
            "reported_power_W":
                best["result"]["indicated_power_w"],
        },
        indent=2,
    ))


if __name__ == "__main__":
    main()
