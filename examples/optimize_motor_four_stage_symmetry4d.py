"""Targeted 4D symmetry test around the retained four-stage motor champion.

Coordinates are absolute event separations in forward motor-cycle degrees:
    delta0_deg = S(T0) - L(T0)
    delta1_deg = S(T1) - L(T1)
    delta2_deg = S(T2) - L(T2)
    delta3_deg = S(T3) - L(T3)

T0 of the large cylinder remains the phase gauge. For T1/T2/T3, the pair mean
is frozen at the retained shared-7D timing. Motion levels a/b and K2 hardware
are frozen. Each candidate is refilled at 100 kPa at global t=0.

Search order: exact symmetric seed; axial probes at +/-0.5, +/-1, +/-2,
+/-4 degrees; then local scrambled 4D Sobol points.
"""
from __future__ import annotations

from dataclasses import replace
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.stats import qmc

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.configuration import ChargeConfiguration
from dada_solver.independent_four_stage_kinematics import IndependentFourStageVolumeKinematics
from dada_solver.integration import IntegrationInterrupted
from compare_motor_motion_laws_stage7A5 import A5_CAMPAIGN, ROOT, _candidate_design_and_mass, _evaluate
from optimize_motor_exchanger_asymmetry_stageK1 import _scaled_exchanger
from optimize_motor_piecewise_stageP3 import _feasibility

THERMO3D_REPORT = ROOT / "outputs" / "motor_four_stage_thermo3d" / "report.json"
DEFAULT_DIRECTORY = ROOT / "outputs" / "motor_four_stage_symmetry4d"
NAMES = ("delta0_deg", "delta1_deg", "delta2_deg", "delta3_deg")
AXIAL_MAGNITUDES_DEG = (0.5, 1.0, 2.0, 4.0)
SOBOL_SEED = 260919
CHARGE_PRESSURE_PA = 100_000.0
MINIMUM_STAGE_FRACTION = 0.02


def _load_history(path):
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def _efficiency(record):
    return float(record["result"]["indicated_thermal_efficiency"])


def _eligible(record):
    return bool(record.get("feasible") and record.get("result", {}).get("status") == "converged"
                and record.get("result", {}).get("indicated_thermal_efficiency") is not None)


def _candidate_id(mapping):
    return hashlib.sha256(json.dumps(mapping, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _target_filling_mass(design):
    wrapper = design.build()
    model = getattr(wrapper, "model", wrapper)
    volumes = model.volumes(0.0)
    config = design.configuration
    mass = CHARGE_PRESSURE_PA * volumes.total / (config.gas.gas_constant * config.charge.temperature)
    return float(mass), float(volumes.total)


def _rescale_gas_inventory(state, target_mass):
    values = np.asarray(state, dtype=float).copy()
    source_mass = float(np.sum(values[:8:2]))
    if values.ndim != 1 or len(values) < 8 or not math.isfinite(source_mass) or source_mass <= 0:
        raise ValueError("Invalid warm-start gas inventory.")
    values[:8] *= target_mass / source_mass
    return values


def _kinematic_parameters(delta_deg, motion):
    d0, d1, d2, d3 = np.asarray(delta_deg, dtype=float) / 360.0
    base = np.asarray([motion["t1"], motion["t2"], motion["t3"]], dtype=float)
    event_delta = np.asarray([d1, d2, d3])
    large_abs = base - 0.5 * event_delta
    small_abs = base + 0.5 * event_delta
    small_local = small_abs - d0
    return {
        "t0_s": float(d0),
        "t1_l": float(large_abs[0]), "t2_l": float(large_abs[1]), "t3_l": float(large_abs[2]),
        "t1_s": float(small_local[0]), "t2_s": float(small_local[1]), "t3_s": float(small_local[2]),
        "a_l": float(motion["a_l"]), "b_l": float(motion["b_l"]),
        "a_s": float(motion["a_s"]), "b_s": float(motion["b_s"]),
    }


def _valid(parameters):
    large = np.asarray([0.0, parameters["t1_l"], parameters["t2_l"], parameters["t3_l"], 1.0])
    small = np.asarray([0.0, parameters["t1_s"], parameters["t2_s"], parameters["t3_s"], 1.0])
    return (-0.5 <= parameters["t0_s"] < 0.5
            and np.all(np.diff(large) >= MINIMUM_STAGE_FRACTION)
            and np.all(np.diff(small) >= MINIMUM_STAGE_FRACTION))


def _diagnostics(delta_deg, parameters, motion):
    d = np.asarray(delta_deg, dtype=float)
    base = 360.0 * np.asarray([motion["t1"], motion["t2"], motion["t3"]], dtype=float)
    large_abs = base - 0.5 * d[1:]
    small_abs = base + 0.5 * d[1:]
    large_stages = 360.0 * np.diff([0, parameters["t1_l"], parameters["t2_l"], parameters["t3_l"], 1])
    small_stages = 360.0 * np.diff([0, parameters["t1_s"], parameters["t2_s"], parameters["t3_s"], 1])
    return {
        "event_separation_deg_S_minus_L": d.tolist(),
        "large_absolute_event_degrees": [0.0, *large_abs.tolist()],
        "small_absolute_event_degrees": [float(d[0]), *small_abs.tolist()],
        "pair_mean_T1_T2_T3_degrees": base.tolist(),
        "large_stage_degrees": large_stages.tolist(),
        "small_stage_degrees": small_stages.tolist(),
        "stage_duration_delta_S_minus_L_degrees": (small_stages - large_stages).tolist(),
        "event_separation_rms_deg": float(np.sqrt(np.mean(d*d))),
    }


def _proposals(total, half_width):
    yield {"kind": "symmetric_seed", "parameters": {name: 0.0 for name in NAMES}}
    for magnitude in AXIAL_MAGNITUDES_DEG:
        for coordinate, name in enumerate(NAMES):
            for sign in (-1.0, 1.0):
                v = np.zeros(4); v[coordinate] = sign * magnitude
                yield {"kind": "axial_probe", "axis": name, "axis_magnitude_deg": magnitude,
                       "parameters": {key: float(value) for key, value in zip(NAMES, v, strict=True)}}
    fixed = 1 + 2 * len(NAMES) * len(AXIAL_MAGNITUDES_DEG)
    remaining = max(0, total - fixed)
    if remaining:
        sampler = qmc.Sobol(d=4, scramble=True, seed=SOBOL_SEED)
        exponent = max(1, math.ceil(math.log2(remaining)))
        for point in sampler.random_base2(exponent)[:remaining]:
            v = (2.0 * point - 1.0) * half_width
            yield {"kind": "local_sobol",
                   "parameters": {key: float(value) for key, value in zip(NAMES, v, strict=True)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget-seconds", type=float, default=3600.0)
    parser.add_argument("--evaluations", type=int, default=128)
    parser.add_argument("--candidate-seconds", type=float, default=120.0)
    parser.add_argument("--sobol-half-width-deg", type=float, default=2.0)
    parser.add_argument("--output-directory", type=Path, default=DEFAULT_DIRECTORY)
    args = parser.parse_args()

    thermo = json.loads(THERMO3D_REPORT.read_text())
    control = thermo["control_1bar"]
    hardware = control["parameters"]
    motion = control["motion_parameters"]

    definition = CampaignDefinition(A5_CAMPAIGN)
    base, _legacy_mass, _saved, _a5, _candidate = _candidate_design_and_mass(definition)
    charge = ChargeConfiguration(temperature=base.configuration.charge.temperature, pressure=CHARGE_PRESSURE_PA)
    base = replace(base, configuration=replace(base.configuration, charge=charge),
                   heat_in=_scaled_exchanger(base.heat_in, float(hardware["k_i"])),
                   heat_out=_scaled_exchanger(base.heat_out, float(hardware["k_o"])))
    limits = base.configuration.machine_volumes

    directory = args.output_directory; directory.mkdir(parents=True, exist_ok=True)
    history_path = directory / "history.jsonl"; definition_path = directory / "definition.json"
    identity = {
        "names": list(NAMES), "axial_magnitudes_deg": list(AXIAL_MAGNITUDES_DEG),
        "sobol_half_width_deg": args.sobol_half_width_deg, "sobol_seed": SOBOL_SEED,
        "charge_pressure_pa": CHARGE_PRESSURE_PA,
        "frozen_pair_mean_times": {k: float(motion[k]) for k in ("t1", "t2", "t3")},
        "frozen_levels": {k: float(motion[k]) for k in ("a_l", "b_l", "a_s", "b_s")},
        "K2": {"k_i": float(hardware["k_i"]), "k_o": float(hardware["k_o"])},
        "strategy": "Exact symmetric seed, signed axial probes, then local 4D Sobol; only relative S-L event timing varies.",
    }
    if definition_path.exists() and json.loads(definition_path.read_text()) != identity:
        raise ValueError("Symmetry-4D definition changed; use a new output directory.")
    definition_path.write_text(json.dumps(identity, indent=2) + "\n")

    history = _load_history(history_path)
    existing = {r["candidate_id"] for r in history if r.get("result", {}).get("status") != "interrupted"}
    eligible = [r for r in history if _eligible(r)]
    best = max(eligible, key=_efficiency) if eligible else None
    start = time.monotonic(); deadline = start + args.budget_seconds

    def progress(_):
        if time.monotonic() >= deadline:
            raise IntegrationInterrupted("4D symmetry-test wall-clock budget exhausted.")
        if time.monotonic() >= candidate_deadline:
            raise IntegrationInterrupted("Candidate wall-clock budget exhausted.")

    for proposal_index, proposal in enumerate(_proposals(args.evaluations, args.sobol_half_width_deg)):
        if time.monotonic() >= deadline: break
        delta = proposal["parameters"]; cid = _candidate_id(delta)
        if cid in existing: continue
        vector = np.asarray([delta[name] for name in NAMES], dtype=float)
        kp = _kinematic_parameters(vector, motion)
        if not _valid(kp): continue
        kin = IndependentFourStageVolumeKinematics(limits.small_cylinder, limits.large_cylinder, **kp)
        design = replace(base, kinematics=kin)
        target_mass, filling_volume = _target_filling_mass(design)

        reusable = [r for r in history if r.get("result", {}).get("status") == "converged" and r.get("last_complete_state")]
        if reusable:
            source = min(reusable, key=lambda r: np.linalg.norm(vector - np.asarray([r["parameters"][n] for n in NAMES])))
            initial = _rescale_gas_inventory(source["last_complete_state"], target_mass); warm = source["candidate_id"]
        else:
            initial = _rescale_gas_inventory(control["last_complete_state"], target_mass); warm = "thermo3d_K2_1bar_control"

        before = time.monotonic(); candidate_deadline = min(deadline, before + args.candidate_seconds)
        try:
            result, state = _evaluate(f"four_stage_symmetry4d_{len(history)}", design, definition,
                                      initial_state=initial, progress_callback=progress)
        except IntegrationInterrupted as exc:
            result, state = {"status": "interrupted", "message": str(exc)}, None
        except (ValueError, RuntimeError) as exc:
            result, state = {"status": "integration_failure", "message": str(exc)}, None
        feasible, reasons = _feasibility(result) if result.get("status") == "converged" else (False, [])
        record = {
            "index": len(history), "proposal_index": proposal_index, "candidate_id": cid,
            "kind": proposal["kind"], "axis": proposal.get("axis"),
            "axis_magnitude_deg": proposal.get("axis_magnitude_deg"), "parameters": delta,
            "independent_kinematics_parameters": kp, "timing_diagnostics": _diagnostics(vector, kp, motion),
            "filling_total_volume_m3": filling_volume, "target_filling_mass_kg": target_mass,
            "result": result, "feasible": feasible, "physical_constraint_failures": reasons,
            "elapsed_seconds": time.monotonic() - before, "warm_start_source": warm,
            "last_complete_state": state.tolist() if state is not None else None,
        }
        with history_path.open("a") as stream: stream.write(json.dumps(record) + "\n")
        history.append(record)
        if result.get("status") != "interrupted": existing.add(cid)
        if feasible and (best is None or _efficiency(record) > _efficiency(best)): best = record

        seed = next((r for r in history if r.get("kind") == "symmetric_seed" and _eligible(r)), None)
        axis_summary = {}
        for name in NAMES:
            rr = [r for r in history if r.get("kind") == "axial_probe" and r.get("axis") == name and _eligible(r)]
            axis_summary[name] = sorted([{
                "value_deg": r["parameters"][name], "efficiency": _efficiency(r),
                "delta_efficiency_vs_seed": _efficiency(r)-_efficiency(seed) if seed else None,
                "power_W": r["result"]["indicated_power_w"]} for r in rr], key=lambda x: x["value_deg"])
        report = {"best_feasible": best, "symmetric_seed": seed,
                  "external_7d_control": {"efficiency": control["result"]["indicated_thermal_efficiency"],
                                          "power_W": control["result"]["indicated_power_w"]},
                  "axis_summary": axis_summary, "attempted_total": len(history), "definition": identity}
        (directory / "report.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps({"index": record["index"], "kind": record["kind"], "delta_deg": delta,
                          "status": result.get("status"), "feasible": feasible,
                          "efficiency": result.get("indicated_thermal_efficiency"),
                          "best_efficiency": _efficiency(best) if best else None,
                          "elapsed_seconds": record["elapsed_seconds"]}), flush=True)

    print(f"Saved {directory / 'report.json'}", flush=True)


if __name__ == "__main__":
    main()
