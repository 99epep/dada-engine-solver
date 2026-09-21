"""Controlled motion-law comparison around the Stage 7A5 motor champion.

Compare, with identical gas inventory, cylinder volume limits, operation and
microtube hardware:

1. the Stage 7A5 four-bar champion;
2. a purely harmonic two-cylinder law, optimizing only the relative small-
   cylinder phase;
3. the repository's historical ideal_piecewise_linear reference law.

The comparison fixes total working-gas mass to the Stage 7A5 periodic inventory.
This avoids changing gas inventory merely because another motion law has
different cylinder volumes at the arbitrary theta=0 filling origin.

No mechanical loss, force, stress or inertia model is introduced here.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import json
import math
from pathlib import Path

import numpy as np

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.campaign.evaluator import (
    MachineEvaluator,
    finalize_microtube_validity,
    json_values,
)
from dada_solver.configuration import (
    ChargeConfiguration,
    load_simulation_configuration,
)
from dada_solver.exchangers.air_wall import AirWallMotor
from dada_solver.exchangers.wall_cycle import (
    WallDiagnosticCycle,
    convergence_summary,
    solve_periodic_wall_motor,
    wall_cycle_performance,
)
from dada_solver.factory import build_initial_state
from dada_solver.kinematics import (
    HarmonicVolumeKinematics,
    IdealPiecewiseLinearVolumeKinematics,
)
from dada_solver.machine import MachineDesign
from dada_solver.results import extract_cycle_diagnostics
from dada_solver.state import ThermodynamicState, UniformCharge
from dada_solver.validity import assess_cycle_validity


ROOT = Path(__file__).resolve().parents[1]
A5_CAMPAIGN = ROOT / "examples" / "motor_mechanics_stage7A5_edge.toml"
A5_REPORT = ROOT / "outputs" / "motor_mechanics_stage7A5_edge" / "report.json"
PIECEWISE_REFERENCE = ROOT / "examples" / "motor_demonstrator_piecewise.toml"
OUTPUT_JSON = ROOT / "outputs" / "motor_motion_law_comparison_stage7A5.json"
PHASE_CSV = ROOT / "outputs" / "motor_motion_law_harmonic_phase_scan_stage7A5.csv"


def _uniform_wall_initial(config, wrapper: AirWallMotor) -> np.ndarray:
    base_initial = build_initial_state(config, wrapper.model)
    volumes = wrapper.model.volumes(0.0)
    pressure = (
        base_initial.total_mass
        * config.gas.gas_constant
        * config.charge.temperature
        / volumes.total
    )
    gas = UniformCharge(pressure, config.charge.temperature).create_state(
        config.gas, volumes
    ).as_array()
    return np.r_[
        gas,
        wrapper.heat_in.wall_capacity_j_k * wrapper.heat_in.air_inlet_temperature_k,
        wrapper.heat_out.wall_capacity_j_k * wrapper.heat_out.air_inlet_temperature_k,
    ]


def _evaluate(
    label: str,
    design: MachineDesign,
    definition: CampaignDefinition,
    *,
    initial_state: np.ndarray | None = None,
    progress_callback=None,
    statistics_callback=None, periodic_observer=None, measure_rhs_time=False,
) -> tuple[dict, np.ndarray | None]:
    import time
    def measured(phase,function,*args,**kwargs):
        before=time.perf_counter()
        try: return function(*args,**kwargs)
        finally:
            if statistics_callback is not None:
                statistics_callback(dict(phase=phase,elapsed_seconds=time.perf_counter()-before))
    wrapper = measured('build',design.build)
    if not isinstance(wrapper, AirWallMotor):
        raise TypeError("This comparison requires the dynamic-wall microtube motor.")

    if initial_state is None:
        initial_state = _uniform_wall_initial(design.configuration, wrapper)

    periodic = measured(
        "periodic_integration", solve_periodic_wall_motor, wrapper,
        initial_state,
        maximum_cycles=design.configuration.numerical.maximum_cycles,
        settings=definition.wall_numerical_settings,
        **({"backend": definition.wall_backend} if hasattr(definition,"wall_backend") else {}),
        exact_kinematics_cache=getattr(definition, "exact_kinematics_cache", True),
        adaptive_acceleration=getattr(definition, "adaptive_wall_acceleration", None),
        progress_callback=progress_callback,
        statistics_callback=statistics_callback, measure_rhs_time=measure_rhs_time,
    )
    if periodic_observer is not None: periodic_observer(periodic)

    result = {
        "label": label,
        "status": periodic.status,
        "message": periodic.message,
        "convergence": convergence_summary(periodic.history),
    }
    if not periodic.converged:
        return result, periodic.last_complete_state

    trajectory = periodic.trajectory
    angles = periodic.angles
    assert trajectory is not None
    assert angles is not None

    performance = wall_cycle_performance(wrapper, trajectory)
    cycle = WallDiagnosticCycle(angles, trajectory)
    from dada_solver.diagnostic_replay import replay_wall_trajectory
    replay = measured('shared_replay', replay_wall_trajectory, wrapper, cycle, trajectory) if getattr(definition,'shared_replay',True) else None

    def wall_heat_rates(index: int, _angle: float, gas_state: ThermodynamicState):
        incoming, outgoing = wrapper.thermal_rates(_angle, trajectory[:, index])
        return incoming['gas_heat_w'], outgoing['gas_heat_w']

    diagnostics = measured(
        "generic_diagnostics", extract_cycle_diagnostics, cycle,
        wrapper.model,
        heat_rate_provider=wall_heat_rates if replay is None else None, replay=replay,
    )
    validity = measured(
        "generic_validity", assess_cycle_validity, cycle,
        wrapper.model,
        design.configuration.validity, replay=replay,
    )
    helper = MachineEvaluator(definition)
    from dada_solver.exchangers.gas_diagnostics import cycle_microtube_diagnostics
    gas_domains = measured("microtube_diagnostics",cycle_microtube_diagnostics,wrapper, angles, trajectory,replay=replay)
    maximum_reynolds, maximum_mach = helper._tube_validity(
        wrapper, angles, trajectory, gas_domains=gas_domains,replay=replay
    )
    validity = finalize_microtube_validity(
        validity,
        maximum_reynolds,
        maximum_mach,
        design.configuration.validity.maximum_mach_number,
        requires_laminar=gas_domains is None,
        domain_failures=gas_domains['failed_criteria'] if gas_domains else (),
    )

    max_abs_flow = max(
        max(abs(ext.minimum), abs(ext.maximum))
        for ext in diagnostics.mass_flow_extrema.values()
    )

    result.update(
        {
            "indicated_power_w": performance.gas_power,
            "heat_input_w": performance.heat_in_power,
            "heat_out_w": performance.heat_out_power,
            "indicated_thermal_efficiency": performance.thermal_efficiency,
            "maximum_tube_reynolds": maximum_reynolds,
            "maximum_tube_mach_number": maximum_mach,
            "maximum_absolute_mass_flow_kg_s": max_abs_flow,
            "validity": json_values(asdict(validity)),
            "diagnostics": json_values(asdict(diagnostics)),
            "microtube_gas_domains": gas_domains,
            "total_mass_kg": float(trajectory[:8:2, -1].sum()),
        }
    )
    return result, periodic.last_complete_state


def _normalized_motion(kinematics, angles: np.ndarray) -> np.ndarray:
    limits_s = kinematics.small_volume_limits
    limits_l = kinematics.large_volume_limits
    small = np.array(
        [kinematics.small_cylinder_volume(float(theta)) for theta in angles]
    )
    large = np.array(
        [kinematics.large_cylinder_volume(float(theta)) for theta in angles]
    )
    return np.vstack(
        (
            (small - limits_s.minimum) / limits_s.swept,
            (large - limits_l.minimum) / limits_l.swept,
        )
    )


def _shape_rms(reference, candidate, sample_count: int = 1440) -> dict:
    angles = np.linspace(0.0, 2.0 * math.pi, sample_count, endpoint=False)
    ref = _normalized_motion(reference, angles)
    test = _normalized_motion(candidate, angles)
    rms = np.sqrt(np.mean((test - ref) ** 2, axis=1))
    return {
        "small_normalized_volume_rms": float(rms[0]),
        "large_normalized_volume_rms": float(rms[1]),
        "combined_normalized_volume_rms": float(
            np.sqrt(np.mean((test - ref) ** 2))
        ),
    }


def _candidate_design_and_mass(definition: CampaignDefinition):
    report = json.loads(A5_REPORT.read_text())
    best = report["best_at_phase_end"]
    candidate_id = best["candidate_id"]
    candidate_path = (
        ROOT
        / "outputs"
        / "motor_mechanics_stage7A5_edge"
        / "candidates"
        / f"{candidate_id}.json"
    )
    record = json.loads(candidate_path.read_text())

    physical = dict(definition.fixed_parameters)
    physical.update(best["physical"])
    design = definition.adapter.build(physical)

    saved_state = record["final_periodic_state"]
    total_mass = float(saved_state["total_mass_kg"])
    state = np.asarray(saved_state["values"], dtype=float)
    return design, total_mass, state, best, candidate_id


def _same_inventory_design(
    base: MachineDesign,
    total_mass: float,
    kinematics,
) -> MachineDesign:
    charge = ChargeConfiguration(
        temperature=base.configuration.charge.temperature,
        total_mass=total_mass,
    )
    config = replace(base.configuration, charge=charge)
    return replace(base, configuration=config, kinematics=kinematics)


def _phase_key(degrees: float) -> float:
    return round(degrees % 360.0, 9)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse-step-deg", type=float, default=30.0)
    parser.add_argument("--refine-step-deg", type=float, default=5.0)
    parser.add_argument("--refine-half-width-deg", type=float, default=20.0)
    parser.add_argument("--output", type=Path, default=OUTPUT_JSON)
    parser.add_argument("--phase-csv", type=Path, default=PHASE_CSV)
    args = parser.parse_args()

    for name, value in (
        ("coarse step", args.coarse_step_deg),
        ("refine step", args.refine_step_deg),
        ("refine half width", args.refine_half_width_deg),
    ):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive.")

    definition = CampaignDefinition(A5_CAMPAIGN)
    base_design, total_mass, saved_a5_state, best, candidate_id = (
        _candidate_design_and_mass(definition)
    )
    a5_kinematics = base_design.kinematics
    assert a5_kinematics is not None

    # Re-evaluate A5 with exactly the same inventory, using its saved periodic
    # state as an initial guess.
    a5_design = _same_inventory_design(
        base_design, total_mass, a5_kinematics
    )
    a5_result, _ = _evaluate(
        "stage7A5_four_bar",
        a5_design,
        definition,
        initial_state=saved_a5_state,
    )

    small_limits = base_design.configuration.machine_volumes.small_cylinder
    large_limits = base_design.configuration.machine_volumes.large_cylinder

    phase_results: dict[float, dict] = {}

    def evaluate_phase(degrees: float) -> dict:
        key = _phase_key(degrees)
        if key in phase_results:
            return phase_results[key]
        kin = HarmonicVolumeKinematics(
            small_limits,
            large_limits,
            math.radians(key),
        )
        design = _same_inventory_design(base_design, total_mass, kin)
        result, _ = _evaluate(
            f"harmonic_{key:.3f}_deg",
            design,
            definition,
        )
        result["small_phase_offset_degrees"] = key
        result["shape_rms_vs_stage7A5"] = _shape_rms(
            a5_kinematics, kin
        )
        phase_results[key] = result

        eta = result.get("indicated_thermal_efficiency")
        eta_text = (
            "unavailable" if eta is None else f"{100.0 * eta:.6f}%"
        )
        power = result.get("indicated_power_w")
        power_text = (
            "unavailable" if power is None else f"{power:.3f} W"
        )
        print(
            f"harmonic {key:8.3f} deg  "
            f"eta={eta_text:>12s}  P={power_text}"
        )
        return result

    print("\nCoarse harmonic phase scan")
    coarse = np.arange(0.0, 360.0, args.coarse_step_deg)
    for phase in coarse:
        evaluate_phase(float(phase))

    valid_coarse = [
        r
        for r in phase_results.values()
        if r.get("indicated_thermal_efficiency") is not None
        and math.isfinite(r["indicated_thermal_efficiency"])
    ]
    if not valid_coarse:
        raise RuntimeError(
            "No harmonic phase produced a usable motor result."
        )
    coarse_best = max(
        valid_coarse,
        key=lambda r: r["indicated_thermal_efficiency"],
    )
    center = float(coarse_best["small_phase_offset_degrees"])

    print("\nRefined harmonic phase scan")
    offsets = np.arange(
        -args.refine_half_width_deg,
        args.refine_half_width_deg + 0.5 * args.refine_step_deg,
        args.refine_step_deg,
    )
    for offset in offsets:
        evaluate_phase(center + float(offset))

    valid_harmonic = [
        r
        for r in phase_results.values()
        if r.get("indicated_thermal_efficiency") is not None
        and math.isfinite(r["indicated_thermal_efficiency"])
    ]
    harmonic_best = max(
        valid_harmonic,
        key=lambda r: r["indicated_thermal_efficiency"],
    )
    best_phase = float(
        harmonic_best["small_phase_offset_degrees"]
    )
    best_harmonic_kinematics = HarmonicVolumeKinematics(
        small_limits,
        large_limits,
        math.radians(best_phase),
    )

    # Historical repository reference, not claimed to be an optimum.
    piecewise_config = load_simulation_configuration(
        PIECEWISE_REFERENCE
    )
    piecewise_kinematics = IdealPiecewiseLinearVolumeKinematics(
        small_limits,
        large_limits,
        piecewise_config.small_lambda_target,
        piecewise_config.large_lambda_target,
        piecewise_config.adiabatic_sector_fraction,
    )
    piecewise_design = _same_inventory_design(
        base_design,
        total_mass,
        piecewise_kinematics,
    )
    print("\nHistorical ideal piecewise-linear reference")
    piecewise_result, _ = _evaluate(
        "historical_ideal_piecewise_linear",
        piecewise_design,
        definition,
    )
    piecewise_result["parameters"] = {
        "small_lambda_target": piecewise_config.small_lambda_target,
        "large_lambda_target": piecewise_config.large_lambda_target,
        "adiabatic_sector_fraction": (
            piecewise_config.adiabatic_sector_fraction
        ),
    }
    piecewise_result["shape_rms_vs_stage7A5"] = _shape_rms(
        a5_kinematics,
        piecewise_kinematics,
    )

    harmonic_best["shape_rms_vs_stage7A5"] = _shape_rms(
        a5_kinematics,
        best_harmonic_kinematics,
    )

    output = {
        "experiment": (
            "Stage 7A5 controlled motion-law comparison"
        ),
        "comparison_basis": {
            "same_total_working_gas_mass": True,
            "total_mass_kg": total_mass,
            "same_cylinder_volume_limits": True,
            "same_operation": True,
            "same_microtube_hardware": True,
            "mechanical_losses_modeled": False,
            "stage7A5_candidate_id": candidate_id,
            "stage7A5_physical": best["physical"],
            "piecewise_reference_source": str(
                PIECEWISE_REFERENCE.relative_to(ROOT)
            ),
        },
        "stage7A5_four_bar": a5_result,
        "best_harmonic": harmonic_best,
        "historical_ideal_piecewise_linear": piecewise_result,
        "harmonic_phase_scan": sorted(
            phase_results.values(),
            key=lambda r: r["small_phase_offset_degrees"],
        ),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.phase_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(json_values(output), indent=2) + "\n"
    )

    rows = [
        "phase_deg,eta,indicated_power_w,heat_input_w,"
        "max_mach,combined_shape_rms"
    ]
    for result in sorted(
        phase_results.values(),
        key=lambda r: r["small_phase_offset_degrees"],
    ):
        rows.append(
            ",".join(
                (
                    f'{result["small_phase_offset_degrees"]:.9g}',
                    (
                        ""
                        if result.get(
                            "indicated_thermal_efficiency"
                        ) is None
                        else f'{result["indicated_thermal_efficiency"]:.12g}'
                    ),
                    (
                        ""
                        if result.get("indicated_power_w") is None
                        else f'{result["indicated_power_w"]:.12g}'
                    ),
                    (
                        ""
                        if result.get("heat_input_w") is None
                        else f'{result["heat_input_w"]:.12g}'
                    ),
                    (
                        ""
                        if result.get(
                            "maximum_tube_mach_number"
                        ) is None
                        else f'{result["maximum_tube_mach_number"]:.12g}'
                    ),
                    f'{result["shape_rms_vs_stage7A5"]["combined_normalized_volume_rms"]:.12g}',
                )
            )
        )
    args.phase_csv.write_text("\n".join(rows) + "\n")

    def summary(label: str, result: dict) -> None:
        eta = result.get("indicated_thermal_efficiency")
        power = result.get("indicated_power_w")
        qin = result.get("heat_input_w")
        eta_text = "n/a" if eta is None else f"{100 * eta:.6f}%"
        power_text = "n/a" if power is None else f"{power:.3f} W"
        qin_text = "n/a" if qin is None else f"{qin:.3f} W"
        print(
            f"{label:36s} "
            f"eta={eta_text:>11s}  "
            f"P={power_text:>10s}  "
            f"Qin={qin_text:>11s}"
        )

    print("\nSUMMARY")
    summary("Stage 7A5 four-bar", a5_result)
    summary(
        f"Best harmonic ({best_phase:.3f} deg)",
        harmonic_best,
    )
    summary(
        "Historical piecewise-linear",
        piecewise_result,
    )
    print(f"\nWrote {args.output.relative_to(ROOT)}")
    print(f"Wrote {args.phase_csv.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
