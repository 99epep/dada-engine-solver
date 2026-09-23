#!/usr/bin/env python3
"""Synthesize an independent six-bar pair for structured-C2 candidate 3952.

This driver reuses the validated Stage-2F / Stage-L1 six-bar synthesis code
without changing the historical crossed-F6 target or the old synthesis outputs.

Pipeline
--------
1. Extract one structured 15-parameter candidate from
   outputs/motor_hybrid_c2_15p_260k/history.jsonl.
2. Export its normalized q, dq/dtheta and d2q/dtheta2 motion as a synthesis target.
3. Locally refine the existing small-cylinder six-bar around the previously
   selected mechanism.
4. Mirror/refine that result for the large-cylinder law with Stage L1.
5. Compare position, velocity and acceleration, including dedicated HP-side
   diagnostics around the two structured kink centers.
6. Replay the actual six-bar pair in the same 260 K thermodynamic configuration.

The geometry search deliberately keeps the historical Stage-2F/L1 objective:
position is dominant and velocity has the existing moderate weight. Acceleration
is *diagnostic only*. In particular, the sharp HP-side acceleration feature of
the structured target is not forced onto the mechanism. The thermodynamic replay
decides whether any natural mechanical smoothing is acceptable.

Run
---
PYTHONPATH=src python3 examples/search_motor_hybrid_3952_sixbar.py

A lighter first pass can be run with e.g.
PYTHONPATH=src python3 examples/search_motor_hybrid_3952_sixbar.py \
  --small-iterations 180 --large-iterations 240 --restarts 2
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
import argparse
import csv
import json
import math
import sys
import time

import numpy as np

from dada_solver.campaign.evaluator import json_values
from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
from dada_solver.six_bar import (
    IndependentSixBarVolumeKinematics,
    load_six_bar_mechanism,
)
from dada_solver.wall_backend import WallBackendSettings

import optimize_motor_hybrid_c2_15p_260k as hybrid15
import search_small_sixbar_stage2f as small_search
import search_large_sixbar_stageL1 as large_search
from fit_motor_hybrid_c2_15p_to_source16 import StructuredMotion15


ROOT = Path.cwd()

DEFAULT_HISTORY = ROOT / "outputs/motor_hybrid_c2_15p_260k/history.jsonl"
DEFAULT_TARGET = ROOT / "outputs/motor_hybrid_c2_15p_260k/candidate_3952_motion_target.csv"
DEFAULT_TARGET_META = ROOT / "outputs/motor_hybrid_c2_15p_260k/candidate_3952_motion_target.json"

OLD_SMALL = ROOT / "outputs/small_sixbar_stage2f_r4_freeH_tightaxis.json"
OLD_LARGE = ROOT / "outputs/large_sixbar_stageL1.json"

DEFAULT_SMALL = ROOT / "outputs/small_sixbar_3952.json"
DEFAULT_SMALL_SEED = ROOT / "outputs/small_sixbar_3952_best_seed.json"
DEFAULT_LARGE = ROOT / "outputs/large_sixbar_3952.json"
DEFAULT_COMPARE = ROOT / "outputs/sixbar_3952_motion_comparison.csv"
DEFAULT_REPORT = ROOT / "outputs/sixbar_3952_synthesis_report.json"
DEFAULT_THERMO = ROOT / "outputs/motor_hybrid_3952_sixbar_thermo.json"


def _load_candidate(path: Path, requested: int) -> tuple[dict, str]:
    """Prefer exact history index; fall back to proposal_index."""
    by_index = None
    by_proposal = None

    with path.open("r", encoding="utf-8") as stream:
        for line in stream:
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)
            if record.get("index") == requested:
                if by_index is not None:
                    raise RuntimeError(f"Duplicate index={requested} in {path}")
                by_index = record
            if record.get("proposal_index") == requested:
                if by_proposal is not None:
                    raise RuntimeError(f"Duplicate proposal_index={requested} in {path}")
                by_proposal = record

    if by_index is not None:
        record, field = by_index, "index"
    elif by_proposal is not None:
        record, field = by_proposal, "proposal_index"
    else:
        raise RuntimeError(
            f"Candidate {requested} was not found as index or proposal_index in {path}"
        )

    params = record.get("parameters")
    if not isinstance(params, dict):
        raise RuntimeError("Selected history record has no structured15 parameters.")

    return record, field


def _hp_centers(params: dict) -> dict[str, float]:
    """Return the two HP-side structured feature centers in normalized t degrees."""
    small = (
        float(params["small_max_deg"])
        + float(params["small_down_kink_u"])
        * float(params["small_down_duration_deg"])
    ) % 360.0

    ld = float(params["large_down_duration_deg"])
    large_up_duration = 360.0 - ld
    large = (
        ld
        + float(params["large_up_kink_u"]) * large_up_duration
    ) % 360.0

    return {"small_t_deg": small, "large_t_deg": large}


def _export_target(
    record: dict,
    selected_field: str,
    requested: int,
    csv_path: Path,
    json_path: Path,
    step_deg: float,
) -> dict:
    params = {k: float(v) for k, v in record["parameters"].items()}
    motion = StructuredMotion15(params)

    samples = int(round(360.0 / step_deg))
    if samples < 360 or not math.isclose(samples * step_deg, 360.0, abs_tol=1e-10):
        raise ValueError("--target-step-deg must divide 360 exactly and be <= 1 degree.")

    theta_deg = np.linspace(0.0, 360.0, samples + 1)
    theta = np.radians(theta_deg)

    # StructuredMotion15 uses normalized time t in the opposite sense to the
    # solver's study-angle convention.
    t = np.mod(-theta / (2.0 * math.pi), 1.0)
    q, dq_dt, ddq_dt2 = motion.normalized(t)
    dq_dtheta = -dq_dt / (2.0 * math.pi)
    ddq_dtheta2 = ddq_dt2 / (2.0 * math.pi) ** 2

    fields = [
        "theta_deg", "theta_rad", "normalized_t_deg",
        "small_fraction_0_1", "large_fraction_0_1",
        "small_dq_dtheta_per_rad", "large_dq_dtheta_per_rad",
        "small_d2q_dtheta2_per_rad2", "large_d2q_dtheta2_per_rad2",
    ]

    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for i in range(len(theta)):
            writer.writerow({
                "theta_deg": float(theta_deg[i]),
                "theta_rad": float(theta[i]),
                "normalized_t_deg": float((360.0 * t[i]) % 360.0),
                "small_fraction_0_1": float(q[0, i]),
                "large_fraction_0_1": float(q[1, i]),
                "small_dq_dtheta_per_rad": float(dq_dtheta[0, i]),
                "large_dq_dtheta_per_rad": float(dq_dtheta[1, i]),
                "small_d2q_dtheta2_per_rad2": float(ddq_dtheta2[0, i]),
                "large_d2q_dtheta2_per_rad2": float(ddq_dtheta2[1, i]),
            })

    hp = _hp_centers(params)
    result = record.get("result", {})
    meta = {
        "requested_candidate_number": requested,
        "selected_history_field": selected_field,
        "index": record.get("index"),
        "proposal_index": record.get("proposal_index"),
        "candidate_id": record.get("candidate_id"),
        "kind": record.get("kind"),
        "parameters": params,
        "shape_diagnostics": record.get("shape_diagnostics"),
        "baseline": {
            "status": result.get("status"),
            "feasible": record.get("feasible"),
            "indicated_thermal_efficiency": result.get("indicated_thermal_efficiency"),
            "indicated_power_w": result.get("indicated_power_w"),
            "heat_input_w": result.get("heat_input_w"),
        },
        "hp_feature_centers_normalized_t_deg": hp,
        "target_step_deg": step_deg,
        "target_csv": str(csv_path),
    }
    json_path.write_text(json.dumps(json_values(meta), indent=2) + "\n", encoding="utf-8")
    return meta


def _call_search_main(module, argv: list[str], target: Path) -> None:
    old_argv = sys.argv
    old_target = module.TARGET
    module.TARGET = target
    sys.argv = [str(Path(module.__file__).name), *argv]
    try:
        module.main()
    finally:
        sys.argv = old_argv
        module.TARGET = old_target




def _seed_link_caps(path: Path, restart: int, branch: int, secondary_fraction: float = 0.25) -> tuple[float, float]:
    """Return EF/GF caps that include the complete intended local seed neighborhood."""
    data = json.loads(path.read_text(encoding="utf-8"))
    candidate = None
    for run in data.get("runs", []):
        if (
            int(run.get("restart", -999)) == int(restart)
            and int(run.get("second_branch", 0)) == int(branch)
            and run.get("best_dense") is not None
        ):
            candidate = run["best_dense"]
            break
    if candidate is None:
        raise ValueError(
            f"No seed restart={restart}, branch={branch:+d} found in {path}"
        )
    params = candidate["parameters"]
    ef = float(params["link_EF"])
    gf = float(params["link_GF"])
    # Stage 2F local_bounds() first forms +/- secondary_fraction bounds and then
    # clips their upper end to --max-EF/--max-GF.  Make those caps large enough
    # that the warm start and its complete requested local neighborhood survive.
    scale = 1.0 + float(secondary_fraction)
    return max(3.0, ef * scale + 1e-9), max(5.5, gf * scale + 1e-9)

def _normalize_small_best(source: Path, output: Path) -> tuple[int, int]:
    data = json.loads(source.read_text(encoding="utf-8"))
    best = data.get("best")
    if not best:
        raise RuntimeError(f"{source} has no global best candidate.")
    # Stage 2F marks dense optimizer candidates explicitly.  Its retained
    # seed is also admitted only after the same feasibility screen, but older
    # reports may omit the literal "feasible" field on that seed.
    if best.get("feasible") is False:
        raise RuntimeError(f"{source} global best is explicitly infeasible.")
    best["feasible"] = True

    branch = int(best["second_branch"])
    normalized = {
        "description": "Global best small six-bar repackaged as a Stage-L1 seed.",
        "source_json": str(source),
        "length_unit": "crank_radius",
        "best": best,
        "runs": [{
            "restart": 0,
            "second_branch": branch,
            "best_dense": best,
        }],
    }
    output.write_text(json.dumps(normalized, indent=2) + "\n", encoding="utf-8")
    return 0, branch


def _periodic_second_derivative(dq: np.ndarray, dtheta: float) -> np.ndarray:
    return (np.roll(dq, -1) - np.roll(dq, 1)) / (2.0 * dtheta)


def _cyclic_distance_deg(x: np.ndarray, center: float) -> np.ndarray:
    return np.abs((x - center + 180.0) % 360.0 - 180.0)


def _side_metrics(
    target_q: np.ndarray,
    target_dq: np.ndarray,
    target_ddq: np.ndarray,
    mech_q: np.ndarray,
    mech_dq: np.ndarray,
    mech_ddq: np.ndarray,
    hp_mask: np.ndarray,
) -> dict:
    def rms(x):
        return float(np.sqrt(np.mean(np.asarray(x, dtype=float) ** 2)))

    def metrics(mask):
        return {
            "position_rms": rms((mech_q - target_q)[mask]),
            "velocity_rms_per_rad": rms((mech_dq - target_dq)[mask]),
            "acceleration_rms_per_rad2": rms((mech_ddq - target_ddq)[mask]),
            "target_peak_abs_velocity_per_rad": float(np.max(np.abs(target_dq[mask]))),
            "mechanism_peak_abs_velocity_per_rad": float(np.max(np.abs(mech_dq[mask]))),
            "target_peak_abs_acceleration_per_rad2": float(np.max(np.abs(target_ddq[mask]))),
            "mechanism_peak_abs_acceleration_per_rad2": float(np.max(np.abs(mech_ddq[mask]))),
        }

    return {
        "whole_cycle": metrics(np.ones_like(hp_mask, dtype=bool)),
        "hp_window": metrics(hp_mask),
    }


def _mechanism_motion(mechanism, theta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    q = np.empty_like(theta)
    dq = np.empty_like(theta)
    for i, angle in enumerate(theta):
        q[i], dq[i] = mechanism.normalized_motion(float(angle))
    return q, dq


def _compare_pair(
    target_path: Path,
    target_meta: dict,
    small_path: Path,
    small_restart: int,
    small_branch: int,
    large_path: Path,
    compare_csv: Path,
    hp_half_width_deg: float,
) -> dict:
    raw = np.genfromtxt(target_path, delimiter=",", names=True)
    if raw["theta_rad"][-1] >= 2.0 * math.pi - 1e-9:
        raw = raw[:-1]

    theta = np.asarray(raw["theta_rad"], dtype=float)
    t_deg = np.asarray(raw["normalized_t_deg"], dtype=float)
    dtheta = float(theta[1] - theta[0])

    small = load_six_bar_mechanism(
        small_path, restart=small_restart, second_branch=small_branch
    )
    large = load_six_bar_mechanism(large_path)

    sq, sdq = _mechanism_motion(small, theta)
    lq, ldq = _mechanism_motion(large, theta)
    sddq = _periodic_second_derivative(sdq, dtheta)
    lddq = _periodic_second_derivative(ldq, dtheta)

    stq = np.asarray(raw["small_fraction_0_1"], dtype=float)
    stdq = np.asarray(raw["small_dq_dtheta_per_rad"], dtype=float)
    stddq = np.asarray(raw["small_d2q_dtheta2_per_rad2"], dtype=float)
    ltq = np.asarray(raw["large_fraction_0_1"], dtype=float)
    ltdq = np.asarray(raw["large_dq_dtheta_per_rad"], dtype=float)
    ltddq = np.asarray(raw["large_d2q_dtheta2_per_rad2"], dtype=float)

    hp = target_meta["hp_feature_centers_normalized_t_deg"]
    smask = _cyclic_distance_deg(t_deg, float(hp["small_t_deg"])) <= hp_half_width_deg
    lmask = _cyclic_distance_deg(t_deg, float(hp["large_t_deg"])) <= hp_half_width_deg

    compare_csv.parent.mkdir(parents=True, exist_ok=True)
    with compare_csv.open("w", newline="", encoding="utf-8") as stream:
        fields = [
            "theta_deg", "theta_rad", "normalized_t_deg",
            "small_target_q", "small_mechanism_q",
            "small_target_dq", "small_mechanism_dq",
            "small_target_ddq", "small_mechanism_ddq",
            "small_hp_window",
            "large_target_q", "large_mechanism_q",
            "large_target_dq", "large_mechanism_dq",
            "large_target_ddq", "large_mechanism_ddq",
            "large_hp_window",
        ]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for i in range(len(theta)):
            writer.writerow({
                "theta_deg": math.degrees(float(theta[i])),
                "theta_rad": float(theta[i]),
                "normalized_t_deg": float(t_deg[i]),
                "small_target_q": float(stq[i]),
                "small_mechanism_q": float(sq[i]),
                "small_target_dq": float(stdq[i]),
                "small_mechanism_dq": float(sdq[i]),
                "small_target_ddq": float(stddq[i]),
                "small_mechanism_ddq": float(sddq[i]),
                "small_hp_window": int(smask[i]),
                "large_target_q": float(ltq[i]),
                "large_mechanism_q": float(lq[i]),
                "large_target_dq": float(ltdq[i]),
                "large_mechanism_dq": float(ldq[i]),
                "large_target_ddq": float(ltddq[i]),
                "large_mechanism_ddq": float(lddq[i]),
                "large_hp_window": int(lmask[i]),
            })

    return {
        "hp_half_width_deg_in_normalized_t": hp_half_width_deg,
        "hp_feature_centers_normalized_t_deg": hp,
        "small": _side_metrics(stq, stdq, stddq, sq, sdq, sddq, smask),
        "large": _side_metrics(ltq, ltdq, ltddq, lq, ldq, lddq, lmask),
        "comparison_csv": str(compare_csv),
    }


def _thermo_replay(
    record: dict,
    small_path: Path,
    small_restart: int,
    small_branch: int,
    large_path: Path,
    output: Path,
) -> dict:
    sp = hybrid15.resolve(hybrid15.DEFAULT_SOURCE)
    fp = hybrid15.resolve(hybrid15.DEFAULT_FOURIER)

    source_report = json.loads(sp.read_text(encoding="utf-8"))
    sb = source_report.get("best_feasible")
    if sb is None:
        raise RuntimeError("Source report has no best_feasible.")

    fourier_report = json.loads(fp.read_text(encoding="utf-8"))
    fb = hybrid15._source_best(fourier_report)
    harmonics = hybrid15._source_harmonics(fb)

    cd, base = hybrid15._load_basis()
    geometry = hybrid15.base_geometry(base)
    machine = hybrid15.fourier_build_design(
        base,
        geometry,
        fb["thermo"],
        np.asarray(fb["small_coefficients"]),
        np.asarray(fb["large_coefficients"]),
        harmonics,
    )
    total_mass = float(sb["result"]["total_mass_kg"])
    machine = hybrid15._fixed_inventory_design(machine, total_mass)
    limits = machine.configuration.machine_volumes

    small = load_six_bar_mechanism(
        small_path, restart=small_restart, second_branch=small_branch
    )
    large = load_six_bar_mechanism(large_path)
    kin = IndependentSixBarVolumeKinematics(
        small, large, limits.small_cylinder, limits.large_cylinder
    )
    design = replace(machine, kinematics=kin)

    definition = SimpleNamespace(
        wall_numerical_settings=cd.wall_numerical_settings,
        wall_backend=WallBackendSettings("numba"),
        exact_kinematics_cache=True,
        shared_replay=True,
    )

    saved_state = record.get("last_complete_state")
    if saved_state is not None:
        initial = np.asarray(saved_state, dtype=float)
        warm_start = "candidate_3952_periodic_state"
    else:
        wall_source = np.asarray(sb["last_complete_state"], dtype=float)
        initial = hybrid15._safe_uniform_initial(design, wall_source=wall_source)
        warm_start = "safe_uniform_from_source_wall_state"

    retry = False
    started = time.monotonic()
    try:
        try:
            result, state = hybrid15._evaluate(
                "motor_hybrid_3952_sixbar",
                design,
                definition,
                initial_state=initial,
            )
        except MicrotubeDomainError:
            retry = True
            safe = hybrid15._safe_uniform_initial(design, wall_source=initial)
            result, state = hybrid15._evaluate(
                "motor_hybrid_3952_sixbar_safe",
                design,
                definition,
                initial_state=safe,
            )
            warm_start = "safe_uniform_retry_after_domain_error"
    except Exception as exc:
        result = {"status": "failed", "message": f"{type(exc).__name__}: {exc}"}
        state = None

    baseline = record.get("result", {})
    eta0 = baseline.get("indicated_thermal_efficiency")
    p0 = baseline.get("indicated_power_w")
    eta = result.get("indicated_thermal_efficiency")
    power = result.get("indicated_power_w")

    report = {
        "description": "Thermodynamic replay of the actual six-bar pair against candidate 3952 hardware.",
        "baseline_candidate": {
            "index": record.get("index"),
            "proposal_index": record.get("proposal_index"),
            "candidate_id": record.get("candidate_id"),
            "indicated_thermal_efficiency": eta0,
            "indicated_power_w": p0,
        },
        "sixbar_result": result,
        "delta": {
            "efficiency_percentage_points": (
                100.0 * (float(eta) - float(eta0))
                if eta is not None and eta0 is not None else None
            ),
            "power_w": (
                float(power) - float(p0)
                if power is not None and p0 is not None else None
            ),
        },
        "warm_start": warm_start,
        "safe_retry_used": retry,
        "elapsed_seconds": time.monotonic() - started,
        "last_complete_state": state.tolist() if state is not None else None,
        "small_mechanism": str(small_path),
        "large_mechanism": str(large_path),
    }
    output.write_text(json.dumps(json_values(report), indent=2) + "\n", encoding="utf-8")
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidate", type=int, default=3952)
    ap.add_argument("--history", type=Path, default=DEFAULT_HISTORY)
    ap.add_argument("--target-step-deg", type=float, default=0.25)
    ap.add_argument("--target", type=Path, default=DEFAULT_TARGET)
    ap.add_argument("--target-metadata", type=Path, default=DEFAULT_TARGET_META)

    ap.add_argument("--old-small", type=Path, default=OLD_SMALL)
    ap.add_argument("--old-large", type=Path, default=OLD_LARGE)
    ap.add_argument("--small-output", type=Path, default=DEFAULT_SMALL)
    ap.add_argument("--small-seed-output", type=Path, default=DEFAULT_SMALL_SEED)
    ap.add_argument("--large-output", type=Path, default=DEFAULT_LARGE)

    ap.add_argument("--small-iterations", type=int, default=450)
    ap.add_argument("--large-iterations", type=int, default=600)
    ap.add_argument("--restarts", type=int, default=4)
    ap.add_argument("--population-size", type=int, default=14)
    ap.add_argument("--small-seed", type=int, default=395200)
    ap.add_argument("--large-seed", type=int, default=395300)

    ap.add_argument(
        "--hp-half-width-deg",
        type=float,
        default=18.0,
        help="Diagnostic HP half-window in normalized structured-motion degrees.",
    )
    ap.add_argument("--comparison-csv", type=Path, default=DEFAULT_COMPARE)
    ap.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    ap.add_argument("--thermo-output", type=Path, default=DEFAULT_THERMO)

    ap.add_argument(
        "--skip-search",
        action="store_true",
        help="Reuse existing small/large synthesis files and only analyze/replay them.",
    )
    ap.add_argument(
        "--skip-thermo",
        action="store_true",
        help="Stop after geometric/kinematic comparison.",
    )
    args = ap.parse_args()

    record, field = _load_candidate(args.history, args.candidate)
    target_meta = _export_target(
        record,
        field,
        args.candidate,
        args.target,
        args.target_metadata,
        args.target_step_deg,
    )

    baseline = target_meta["baseline"]
    print(
        f"Selected candidate via {field}: "
        f"index={target_meta['index']} proposal={target_meta['proposal_index']} "
        f"eta={baseline.get('indicated_thermal_efficiency')} "
        f"P={baseline.get('indicated_power_w')} W",
        flush=True,
    )
    hp = target_meta["hp_feature_centers_normalized_t_deg"]
    print(
        "HP feature centers in structured normalized angle: "
        f"S={hp['small_t_deg']:.3f} deg, L={hp['large_t_deg']:.3f} deg",
        flush=True,
    )

    small_restart = 0
    small_branch = 1

    if not args.skip_search:
        if not args.old_small.exists():
            raise FileNotFoundError(args.old_small)

        small_max_ef, small_max_gf = _seed_link_caps(
            args.old_small, 0, +1, secondary_fraction=0.25
        )
        print(
            f"Small warm-start link caps: max-EF={small_max_ef:.6f}, "
            f"max-GF={small_max_gf:.6f}",
            flush=True,
        )

        _call_search_main(
            small_search,
            [
                str(args.old_small),
                "--seed-restart", "0",
                "--seed-branch", "+1",
                "--iterations", str(args.small_iterations),
                "--population-size", str(args.population_size),
                "--restarts", str(args.restarts),
                "--seed", str(args.small_seed),
                "--stroke-floor", "2.5",
                "--h-line-rms-max", "0.07",
                "--h-axis-angle-max", "5",
                "--crank-clearance-floor", "0.50",
                "--max-EF", f"{small_max_ef:.12g}",
                "--max-GF", f"{small_max_gf:.12g}",
                "--output", str(args.small_output),
            ],
            args.target,
        )

        small_restart, small_branch = _normalize_small_best(
            args.small_output, args.small_seed_output
        )

        large_max_ef, large_max_gf = _seed_link_caps(
            args.small_seed_output,
            small_restart,
            small_branch,
            secondary_fraction=0.30,
        )
        large_max_ef = max(3.5, large_max_ef)
        large_max_gf = max(8.0, large_max_gf)
        print(
            f"Large mirrored-seed link caps: max-EF={large_max_ef:.6f}, "
            f"max-GF={large_max_gf:.6f}",
            flush=True,
        )

        _call_search_main(
            large_search,
            [
                str(args.small_seed_output),
                "--seed-restart", str(small_restart),
                "--seed-branch", f"{small_branch:+d}",
                "--iterations", str(args.large_iterations),
                "--population-size", str(args.population_size),
                "--restarts", str(args.restarts),
                "--seed", str(args.large_seed),
                "--stroke-floor", "2.5",
                "--h-line-rms-max", "0.15",
                "--h-axis-angle-max", "3",
                "--crank-clearance-floor", "0.50",
                "--max-EF", f"{large_max_ef:.12g}",
                "--max-GF", f"{large_max_gf:.12g}",
                "--output", str(args.large_output),
            ],
            args.target,
        )
    else:
        # The normalized seed is preferred because it always contains the
        # selected global small best at restart 0.
        if args.small_seed_output.exists():
            data = json.loads(args.small_seed_output.read_text(encoding="utf-8"))
            best = data["best"]
            small_branch = int(best["second_branch"])
        elif args.small_output.exists():
            small_restart, small_branch = _normalize_small_best(
                args.small_output, args.small_seed_output
            )
        else:
            raise FileNotFoundError(
                "No synthesized small mechanism found for --skip-search."
            )

    if not args.small_seed_output.exists():
        small_restart, small_branch = _normalize_small_best(
            args.small_output, args.small_seed_output
        )
    if not args.large_output.exists():
        raise FileNotFoundError(args.large_output)

    motion_report = _compare_pair(
        args.target,
        target_meta,
        args.small_seed_output,
        small_restart,
        small_branch,
        args.large_output,
        args.comparison_csv,
        args.hp_half_width_deg,
    )

    report = {
        "candidate": target_meta,
        "search_policy": {
            "small_warm_start": str(args.old_small),
            "historical_large_reference": str(args.old_large),
            "small_iterations": args.small_iterations,
            "large_iterations": args.large_iterations,
            "restarts": args.restarts,
            "population_size": args.population_size,
            "objective_note": (
                "Historical six-bar objective retained: position plus moderate "
                "velocity mismatch. Acceleration, including the sharp HP feature, "
                "is diagnostic only and is not forced onto the mechanism."
            ),
        },
        "motion_comparison": motion_report,
        "small_output": str(args.small_output),
        "small_seed_output": str(args.small_seed_output),
        "large_output": str(args.large_output),
    }

    if not args.skip_thermo:
        thermo = _thermo_replay(
            record,
            args.small_seed_output,
            small_restart,
            small_branch,
            args.large_output,
            args.thermo_output,
        )
        report["thermodynamic_replay"] = thermo

    args.report.write_text(
        json.dumps(json_values(report), indent=2) + "\n",
        encoding="utf-8",
    )

    print()
    print("SYNTHESIS SUMMARY")
    for side in ("small", "large"):
        m = motion_report[side]
        full = m["whole_cycle"]
        hp_m = m["hp_window"]
        print(
            f"{side}: pos RMS={100*full['position_rms']:.4f}% "
            f"vel RMS={full['velocity_rms_per_rad']:.5f} "
            f"HP accel peak target/mech="
            f"{hp_m['target_peak_abs_acceleration_per_rad2']:.3f}/"
            f"{hp_m['mechanism_peak_abs_acceleration_per_rad2']:.3f}",
            flush=True,
        )

    thermo = report.get("thermodynamic_replay")
    if thermo:
        delta = thermo["delta"]
        print(
            "thermo: "
            f"status={thermo['sixbar_result'].get('status')} "
            f"dEta={delta.get('efficiency_percentage_points')} pp "
            f"dP={delta.get('power_w')} W",
            flush=True,
        )

    print(f"Wrote {args.report}", flush=True)


if __name__ == "__main__":
    main()
