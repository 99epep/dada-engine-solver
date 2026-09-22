#!/usr/bin/env python3
"""Focused local refinement of the 16-control free-spline motor law at 260 K.

This campaign starts from the exact best spline candidate of a previous report.
Hardware, total gas inventory, source temperatures, topology and all non-motion
parameters remain frozen to the underlying Fourier-8H machine.

Typical use:

    PYTHONPATH=src python3 examples/refine_motor_free_spline_260k.py \
        --source-report outputs/motor_free_spline_260k_v3/report.json \
        --output-directory outputs/motor_free_spline_260k_refine \
        --evaluations 4096 \
        --budget-seconds 21600

The run is append-only and resumable.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
import argparse
import hashlib
import json
import math
import time

import numpy as np
from scipy.stats import qmc

from dada_solver.configuration import ChargeConfiguration
from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
from dada_solver.integration import IntegrationInterrupted
from dada_solver.wall_backend import WallBackendSettings

from compare_motor_motion_laws_stage7A5 import ROOT, _evaluate
from optimize_motor_temperature_point import base_geometry, feasibility
from refine_motor_four_stage_hx9d_variable_gas import _load_basis
from refine_motor_fourier_c2_260k import build_design as fourier_build_design

from optimize_motor_free_spline_260k_v3 import (
    _admissible,
    _build_candidate_design,
    _candidate_id,
    _eta,
    _eligible,
    _fixed_inventory_design,
    _load_history,
    _make_kinematics,
    _normalized_motion,
    _record_distance,
    _safe_uniform_initial,
    _save_motion,
    _shape_diagnostics,
    _shape_rms,
    _source_best,
    _source_harmonics,
    _canonical,
)

POWER_FLOOR_W = 25.0
DEFAULT_SOURCE = ROOT / "outputs" / "motor_free_spline_260k_v3" / "report.json"
DEFAULT_OUTPUT = ROOT / "outputs" / "motor_free_spline_260k_refine"
DEFAULT_FOURIER = ROOT / "outputs" / "motor_fourier_c2_8h_refine" / "report.json"
DEFAULT_SEED = 26092252

# Fine local refinement around the current spline champion.
SHAPE_RADII = (0.0025, 0.0015, 0.0008, 0.0004, 0.0002)
PHASE_RADII_DEG = (0.18, 0.10, 0.05, 0.025, 0.01)
EVALS_PER_RADIUS = 192

# Rare wider probes preserve a small chance of escaping the local basin.
WIDE_PROBE_EVERY = 32
WIDE_SHAPE_RADIUS = 0.007
WIDE_PHASE_RADIUS_DEG = 0.5

DIRECT_PERIODIC_WARM_DISTANCE = 0.040


def _direction(u, center):
    d = 2.0 * np.asarray(u, dtype=float) - 1.0
    d -= np.mean(d)
    d -= center * float(d @ center)
    n = float(np.linalg.norm(d))
    if n < 1e-12:
        d = np.roll(center, 1) - center
        d -= np.mean(d)
        d -= center * float(d @ center)
        n = float(np.linalg.norm(d))
    return d / n


def _phase_difference_deg(a, b):
    return ((float(a) - float(b) + 180.0) % 360.0) - 180.0


def _propose(center, u, radius, phase_radius_deg):
    small0 = np.asarray(center["small_controls"], dtype=float)
    large0 = np.asarray(center["large_controls"], dtype=float)
    n = len(small0)

    ds = _direction(u[:n], small0)
    dl = _direction(u[n:2*n], large0)

    # Continuous radial spread avoids sampling only one shell.
    rs = radius * (0.20 + 0.80 * float(u[2*n]))
    rl = radius * (0.20 + 0.80 * float(u[2*n+1]))
    small = _canonical(small0 + rs * ds)
    large = _canonical(large0 + rl * dl)

    us = (float(u[0]) + float(u[n])) % 1.0
    ul = (float(u[n-1]) + float(u[2*n-1])) % 1.0
    small_phase = (
        float(center["small_phase_deg"])
        + phase_radius_deg * (2.0 * us - 1.0)
    ) % 360.0
    large_phase = (
        float(center["large_phase_deg"])
        + phase_radius_deg * (2.0 * ul - 1.0)
    ) % 360.0
    return small, large, small_phase, large_phase


def _resolve_path(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-report", type=Path, default=DEFAULT_SOURCE)
    ap.add_argument("--fourier-report", type=Path, default=DEFAULT_FOURIER)
    ap.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--evaluations", type=int, default=4096)
    ap.add_argument("--budget-seconds", type=float, default=21600.0)
    ap.add_argument("--candidate-seconds", type=float, default=120.0)
    ap.add_argument("--evaluations-per-radius", type=int, default=EVALS_PER_RADIUS)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = ap.parse_args()

    source_path = _resolve_path(args.source_report)
    fourier_path = _resolve_path(args.fourier_report)
    out = _resolve_path(args.output_directory)
    out.mkdir(parents=True, exist_ok=True)

    source_report = json.loads(source_path.read_text())
    spline_source = source_report.get("best_feasible")
    if spline_source is None:
        raise RuntimeError("Spline source report has no best_feasible.")
    controls = int(spline_source["controls_per_piston"])
    if len(spline_source["small_controls"]) != controls or len(spline_source["large_controls"]) != controls:
        raise RuntimeError("Invalid spline control count.")

    fourier_report = json.loads(fourier_path.read_text())
    fourier_best = _source_best(fourier_report)
    harmonics = _source_harmonics(fourier_best)

    campaign_definition, base = _load_basis()
    geometry = base_geometry(base)
    source_design = fourier_build_design(
        base,
        geometry,
        fourier_best["thermo"],
        np.asarray(fourier_best["small_coefficients"], dtype=float),
        np.asarray(fourier_best["large_coefficients"], dtype=float),
        harmonics,
    )

    total_mass_kg = float(
        source_report.get("source_fourier", {}).get(
            "total_mass_kg",
            spline_source["result"]["total_mass_kg"],
        )
    )
    source_design = _fixed_inventory_design(source_design, total_mass_kg)
    fourier_kinematics = source_design.kinematics

    seed_small = np.asarray(spline_source["small_controls"], dtype=float)
    seed_large = np.asarray(spline_source["large_controls"], dtype=float)
    seed_small_phase = float(spline_source["small_phase_deg"])
    seed_large_phase = float(spline_source["large_phase_deg"])

    seed_kinematics = _make_kinematics(
        source_design.configuration.machine_volumes,
        seed_small,
        seed_large,
        seed_small_phase,
        seed_large_phase,
    )
    if not _admissible(seed_kinematics):
        raise RuntimeError("Source spline champion is not structurally admissible.")

    definition = SimpleNamespace(
        wall_numerical_settings=campaign_definition.wall_numerical_settings,
        wall_backend=WallBackendSettings("numba"),
        exact_kinematics_cache=True,
        shared_replay=True,
    )

    identity = {
        "study": "260 K focused free-spline refinement v1",
        "source_report": str(source_path),
        "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "fourier_report": str(fourier_path),
        "fourier_sha256": hashlib.sha256(fourier_path.read_bytes()).hexdigest(),
        "source_spline_candidate_id": spline_source["candidate_id"],
        "source_spline_efficiency": spline_source["result"]["indicated_thermal_efficiency"],
        "source_spline_power_W": spline_source["result"]["indicated_power_w"],
        "controls_per_piston": controls,
        "fixed_total_mass_kg": total_mass_kg,
        "fixed_thermo": fourier_best["thermo"],
        "shape_radii": list(SHAPE_RADII),
        "phase_radii_deg": list(PHASE_RADII_DEG),
        "evaluations_per_radius": args.evaluations_per_radius,
        "wide_probe_every": WIDE_PROBE_EVERY,
        "wide_shape_radius": WIDE_SHAPE_RADIUS,
        "wide_phase_radius_deg": WIDE_PHASE_RADIUS_DEG,
        "sobol_seed": args.seed,
        "constraints": {
            "exactly_one_minimum_and_one_maximum_per_piston": True,
            "C2_periodic_position": True,
            "hardware_frozen": True,
            "gas_inventory_frozen": True,
            "four_stage_boundaries_imposed": False,
            "velocity_limit_imposed": False,
            "acceleration_limit_imposed": False,
        },
    }

    definition_path = out / "definition.json"
    history_path = out / "history.jsonl"
    report_path = out / "report.json"

    if definition_path.exists():
        if json.loads(definition_path.read_text()) != identity:
            raise RuntimeError("Campaign definition changed. Use a new output directory.")
    else:
        definition_path.write_text(json.dumps(identity, indent=2) + "\n")

    history = _load_history(history_path)
    valid = [r for r in history if _eligible(r)]
    best = max(valid, key=_eta) if valid else None
    finished = {
        r["candidate_id"]
        for r in history
        if r.get("result", {}).get("status") != "interrupted"
    }

    dimension = 2 * controls + 2
    sobol = qmc.Sobol(d=dimension, scramble=True, seed=args.seed)
    points = sobol.random_base2(16)

    seed_record = {
        "small_controls": list(map(float, seed_small)),
        "large_controls": list(map(float, seed_large)),
        "small_phase_deg": seed_small_phase,
        "large_phase_deg": seed_large_phase,
    }

    start = time.monotonic()
    deadline = start + args.budget_seconds
    completed = 0
    proposal_index = max((r.get("proposal_index", -1) for r in history), default=-1) + 1

    while completed < args.evaluations and time.monotonic() < deadline:
        this_proposal = proposal_index
        proposal_index += 1

        if not history and this_proposal == 0:
            small = seed_small.copy()
            large = seed_large.copy()
            small_phase = seed_small_phase
            large_phase = seed_large_phase
            kind = "source_spline_seed"
            radius_index = -1
            radius = 0.0
            phase_radius = 0.0
        else:
            searched = sum(1 for r in history if r.get("kind") != "source_spline_seed")
            radius_index = min(
                searched // args.evaluations_per_radius,
                len(SHAPE_RADII) - 1,
            )

            if searched > 0 and searched % WIDE_PROBE_EVERY == WIDE_PROBE_EVERY - 1:
                center = seed_record
                radius = WIDE_SHAPE_RADIUS
                phase_radius = WIDE_PHASE_RADIUS_DEG
                kind = "rare_wide_probe"
            else:
                center = best if best is not None else seed_record
                radius = SHAPE_RADII[radius_index]
                phase_radius = PHASE_RADII_DEG[radius_index]
                kind = "incumbent_spline_refinement"

            u = points[this_proposal % len(points)]
            small, large, small_phase, large_phase = _propose(
                center, u, radius, phase_radius
            )

        cid = _candidate_id(small, large, small_phase, large_phase)
        if cid in finished:
            continue

        candidate_design = _build_candidate_design(
            source_design, small, large, small_phase, large_phase
        )
        candidate_kinematics = candidate_design.kinematics
        if not _admissible(candidate_kinematics):
            continue

        before = time.monotonic()
        candidate_deadline = min(deadline, before + args.candidate_seconds)

        reusable = [
            r for r in history
            if r.get("result", {}).get("status") == "converged"
            and r.get("last_complete_state") is not None
        ]

        warm_distance = None
        if reusable:
            warm = min(
                reusable,
                key=lambda r: _record_distance(
                    small, large, small_phase, large_phase, r
                ),
            )
            warm_distance = _record_distance(
                small, large, small_phase, large_phase, warm
            )
            warm_source = warm["candidate_id"]
            if warm_distance <= DIRECT_PERIODIC_WARM_DISTANCE:
                initial = np.asarray(warm["last_complete_state"], dtype=float)
                warm_mode = "nearest_periodic_state"
            else:
                initial = _safe_uniform_initial(
                    candidate_design,
                    wall_source=np.asarray(warm["last_complete_state"], dtype=float),
                )
                warm_mode = "direct_safe_uniform_for_distant_shape"
        else:
            initial = np.asarray(spline_source["last_complete_state"], dtype=float)
            warm_source = spline_source["candidate_id"]
            warm_mode = "source_spline_periodic_state"

        def progress(_):
            now = time.monotonic()
            if now >= deadline:
                raise IntegrationInterrupted("Spline refinement wall-clock budget exhausted.")
            if now >= candidate_deadline:
                raise IntegrationInterrupted("Candidate wall-clock budget exhausted.")

        backend_box = {}
        def observe_periodic(periodic):
            backend_box["statistics"] = periodic.backend_statistics

        safe_retry = False
        try:
            try:
                result, state = _evaluate(
                    f"free_spline_refine_260k_{this_proposal}",
                    candidate_design,
                    definition,
                    initial_state=initial,
                    progress_callback=progress,
                    periodic_observer=observe_periodic,
                )
            except MicrotubeDomainError:
                safe_retry = True
                safe = _safe_uniform_initial(candidate_design, wall_source=initial)
                result, state = _evaluate(
                    f"free_spline_refine_260k_{this_proposal}_safe",
                    candidate_design,
                    definition,
                    initial_state=safe,
                    progress_callback=progress,
                    periodic_observer=observe_periodic,
                )
                warm_mode = "safe_uniform_retry_after_domain_error"
        except MicrotubeDomainError as exc:
            result = {"status": "invalid_exchanger", "message": str(exc)}
            state = None
        except IntegrationInterrupted as exc:
            result = {"status": "interrupted", "message": str(exc)}
            state = None
        except (ValueError, RuntimeError, ArithmeticError) as exc:
            result = {"status": "integration_failure", "message": str(exc)}
            state = None

        feasible, reasons = (
            feasibility(result, POWER_FLOOR_W)
            if result.get("status") == "converged"
            else (False, [])
        )

        shape = _shape_diagnostics(candidate_kinematics)
        distance_vs_fourier = _shape_rms(
            fourier_kinematics, candidate_kinematics, sample_count=2048
        )
        distance_vs_seed = _shape_rms(
            seed_kinematics, candidate_kinematics, sample_count=2048
        )

        record = {
            "index": len(history),
            "proposal_index": this_proposal,
            "candidate_id": cid,
            "kind": kind,
            "controls_per_piston": controls,
            "radius_index": radius_index,
            "shape_radius": radius,
            "phase_radius_deg": phase_radius,
            "small_controls": list(map(float, small)),
            "large_controls": list(map(float, large)),
            "small_phase_deg": float(small_phase) % 360.0,
            "large_phase_deg": float(large_phase) % 360.0,
            "shape_diagnostics": shape,
            "shape_distance_vs_source_fourier": distance_vs_fourier,
            "shape_distance_vs_source_spline": distance_vs_seed,
            "result": result,
            "feasible": feasible,
            "physical_constraint_failures": reasons,
            "elapsed_seconds": time.monotonic() - before,
            "warm_start_source": warm_source,
            "warm_start_mode": warm_mode,
            "warm_start_distance": warm_distance,
            "backend_statistics": backend_box.get("statistics"),
            "safe_retry_used": safe_retry,
            "last_complete_state": state.tolist() if state is not None else None,
        }

        with history_path.open("a") as stream:
            stream.write(json.dumps(record) + "\n")
            stream.flush()
        history.append(record)
        completed += 1

        moved = False
        if feasible and (best is None or _eta(record) > _eta(best)):
            best = record
            moved = True
            _save_motion(
                out / "best_motion.csv",
                fourier_kinematics,
                candidate_kinematics,
            )

        report = {
            "best_feasible": best,
            "source_spline": {
                "candidate_id": spline_source["candidate_id"],
                "efficiency": spline_source["result"]["indicated_thermal_efficiency"],
                "power_W": spline_source["result"]["indicated_power_w"],
                "shape_diagnostics": spline_source["shape_diagnostics"],
            },
            "source_fourier": {
                "candidate_id": fourier_best.get("candidate_id"),
                "harmonics": harmonics,
                "efficiency": fourier_best["result"]["indicated_thermal_efficiency"],
                "power_W": fourier_best["result"]["indicated_power_w"],
                "total_mass_kg": total_mass_kg,
                "thermo": fourier_best["thermo"],
            },
            "attempted_total": len(history),
            "completed_this_run": completed,
            "requested_duration_seconds": args.budget_seconds,
            "actual_duration_seconds": time.monotonic() - start,
            "definition": identity,
        }
        report_path.write_text(json.dumps(report, indent=2) + "\n")

        print(json.dumps({
            "index": record["index"],
            "proposal_index": this_proposal,
            "kind": kind,
            "radius": radius,
            "phase_radius_deg": phase_radius,
            "status": result.get("status"),
            "feasible": feasible,
            "efficiency": result.get("indicated_thermal_efficiency"),
            "power_W": result.get("indicated_power_w"),
            "rms_vs_source_spline": distance_vs_seed["combined_position_rms"],
            "center_moved": moved,
            "best_efficiency": _eta(best) if best else None,
            "actual_backend": (
                (backend_box.get("statistics") or {}).get("actual_backend")
            ),
            "elapsed_seconds": record["elapsed_seconds"],
        }), flush=True)

    print(f"Saved {report_path}", flush=True)


if __name__ == "__main__":
    main()
