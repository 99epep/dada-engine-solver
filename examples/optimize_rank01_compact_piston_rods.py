#!/usr/bin/env python3
"""Local 30-D re-optimization of final rank_01 with shorter, aligned piston rods.

Fixed: final rank_01 thermo-5D hardware, gas inventory, source temperatures,
volume split, exchanger geometry and assembly branches.

Free: all 30 continuous six-bar coordinates (15 SMALL + 15 LARGE).

Extra rod constraints:
- each rod is 2..3 crank radii shorter than the current rank_01 rod;
- max instantaneous rod/slider-axis angle <= 8 deg by default;
- absolute mean signed rod/slider-axis angle <= 2 deg by default.

The first seed aligns each slider axis with the PCA best-fit line of H(theta),
recentres the slider through mean(H), and shortens the rod by 2.5 crank radii.

Run from repository root:
  PYTHONPATH=src python3 examples/optimize_rank01_compact_piston_rods.py
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import argparse
import hashlib
import json
import math
import time

import numpy as np
from scipy.stats import qmc

from dada_solver.six_bar import IndependentSixBarVolumeKinematics

import optimize_sixbar_pairs_thermo_3952_fine as fine
import optimize_sixbar_pairs_thermo5d_3952_hlat25 as thermo5d

ROOT = Path.cwd()
PAIR_PATH = ROOT / "outputs/sixbar_thermo_coupled_3952_fine_hlat25/rank_01/best_pair.json"
THERMO_REPORT = ROOT / "outputs/sixbar_thermo5d_3952/rank_01/report.json"
DEFAULT_OUTPUT = ROOT / "outputs/sixbar_rank01_compact_rods"
RANK = 1

DEFAULT_RADII = (1.0, 0.65, 0.40, 0.25)
BASE_ROD_SPAN = 0.50
BASE_AXIS_OFFSET_SPAN = 0.40
BASE_AXIS_ANGLE_DEG = 4.0


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def append_jsonl(path: Path, record: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(record, separators=(",", ":")) + "\n")
        f.flush()


def load_history(path: Path):
    if not path.exists():
        return []
    return [json.loads(x) for x in path.read_text(encoding="utf-8").splitlines() if x.strip()]


def branches_from_pair(pair):
    return {
        "small_primary": int(pair["small"]["primary"]["assembly_branch"]),
        "small_second": int(pair["small"]["second_branch"]),
        "large_primary": int(pair["large"]["primary"]["assembly_branch"]),
        "large_second": int(pair["large"]["second_branch"]),
    }


def source_vector(pair):
    return np.concatenate((fine.vector_from_candidate(pair["small"]), fine.vector_from_candidate(pair["large"])))


def h_trajectory(v, primary_branch, second_branch, theta):
    x = np.asarray(v, float)
    E, _, _, _ = fine.stage.primary(theta, x[:6], primary_branch)
    G = np.array([x[6], x[7]], float)
    ef, gf = float(x[8]), float(x[9])
    d = G[None, :] - E
    dist = np.linalg.norm(d, axis=1)
    if np.any(dist <= 1e-12):
        raise ValueError("secondary centres coincide")
    u = d / dist[:, None]
    along = (ef*ef - gf*gf + dist*dist) / (2*dist)
    h2 = ef*ef - along*along
    if np.any(h2 <= 1e-10):
        raise ValueError("secondary closure/toggle")
    n = np.column_stack((-u[:, 1], u[:, 0]))
    F = E + along[:, None]*u + second_branch*np.sqrt(h2)[:, None]*n
    EF = F - E
    nEF = np.column_stack((-EF[:, 1], EF[:, 0]))
    return E + float(x[10])*EF + float(x[11])*nEF


def aligned_seed_side(v, primary_branch, second_branch, theta):
    y = np.asarray(v, float).copy()
    H = h_trajectory(y, primary_branch, second_branch, theta)
    center = H.mean(axis=0)
    X = H - center
    values, vectors = np.linalg.eigh(X.T @ X)
    line = vectors[:, int(np.argmax(values))]
    line /= np.linalg.norm(line)
    old_axis = np.array([math.cos(y[14]), math.sin(y[14])])
    if float(line @ old_axis) < 0:
        line *= -1
    angle = math.atan2(float(line[1]), float(line[0]))
    normal = np.array([-line[1], line[0]])
    y[12] = max(1e-6, y[12] - 2.5)
    y[13] = float(center @ normal)
    y[14] = fine.wrap_pi(angle)
    return y


def rod_alignment(v, primary_branch, second_branch, theta):
    H = h_trajectory(v, primary_branch, second_branch, theta)
    rod = float(v[12])
    angle = float(v[14])
    normal = np.array([-math.sin(angle), math.cos(angle)])
    transverse = H @ normal - float(v[13])
    if np.any(np.abs(transverse) >= rod):
        raise ValueError("rod cannot close")
    axial = np.sqrt(np.maximum(1e-30, rod*rod - transverse*transverse))
    a = np.degrees(np.arctan2(-transverse, axial))
    return {
        "rod_axis_mean_angle_deg": float(np.mean(a)),
        "rod_axis_rms_angle_deg": float(np.sqrt(np.mean(a*a))),
        "rod_axis_max_abs_angle_deg": float(np.max(np.abs(a))),
        "rod_axis_min_cosine_direct": float(np.min(np.cos(np.radians(a)))),
    }


def rod_bounds(src):
    return {
        "small": (float(src[12]) - 3.0, float(src[12]) - 2.0),
        "large": (float(src[27]) - 3.0, float(src[27]) - 2.0),
    }


def screen_pair(x, branches, motion, bounds, max_angle, max_mean):
    x = np.asarray(x, float)
    sv, lv = x[:15].copy(), x[15:].copy()
    for v in (sv, lv):
        v[5] = fine.wrap_pi(v[5]); v[14] = fine.wrap_pi(v[14])

    sm = fine.side_mechanical_metrics(
        sv, primary_branch=branches["small_primary"], second_branch=branches["small_second"],
        side="small", motion=motion)
    lm = fine.side_mechanical_metrics(
        lv, primary_branch=branches["large_primary"], second_branch=branches["large_second"],
        side="large", motion=motion)
    sm, lm = dict(sm), dict(lm)

    for side, v, m, pb, sb in (
        ("small", sv, sm, branches["small_primary"], branches["small_second"]),
        ("large", lv, lm, branches["large_primary"], branches["large_second"]),
    ):
        lo, hi = bounds[side]
        if not (lo <= float(v[12]) <= hi):
            raise ValueError(f"{side} rod outside [{lo},{hi}]")
        a = rod_alignment(v, pb, sb, motion["theta"])
        m.update(a)
        if a["rod_axis_max_abs_angle_deg"] > max_angle:
            raise ValueError(f"{side} max rod angle")
        if abs(a["rod_axis_mean_angle_deg"]) > max_mean:
            raise ValueError(f"{side} mean rod angle")

    smech = fine.mechanism_from_vector(sv, branches["small_primary"], branches["small_second"])
    lmech = fine.mechanism_from_vector(lv, branches["large_primary"], branches["large_second"])
    return sv, lv, sm, lm, smech, lmech


def perturb_side(c, u, radius, bounds):
    c = np.asarray(c, float); z = 2*np.asarray(u, float)-1; y = c.copy()
    for i in (0,1,2):
        y[i] = max(1e-4, c[i]*(1 + radius*fine.BASE_PRIMARY_LENGTH_FRAC*z[i]))
    for i in (3,4):
        span = max(abs(c[i])*fine.BASE_PRIMARY_POINT_FRAC, fine.BASE_PRIMARY_POINT_MIN)
        y[i] = c[i] + radius*span*z[i]
    y[5] = fine.wrap_pi(c[5] + math.radians(fine.BASE_PRIMARY_PHASE_DEG)*radius*z[5])
    y[6] = c[6] + fine.BASE_PIVOT_SPAN*radius*z[6]
    y[7] = c[7] + fine.BASE_PIVOT_SPAN*radius*z[7]
    for i in (8,9):
        y[i] = max(1e-4, c[i]*(1 + radius*fine.BASE_SECONDARY_LENGTH_FRAC*z[i]))
    y[10] = c[10] + fine.BASE_H_RATIO_SPAN*radius*z[10]
    y[11] = c[11] + fine.BASE_H_RATIO_SPAN*radius*z[11]
    lo, hi = bounds
    y[12] = float(np.clip(c[12] + radius*BASE_ROD_SPAN*z[12], lo, hi))
    y[13] = c[13] + radius*BASE_AXIS_OFFSET_SPAN*z[13]
    y[14] = fine.wrap_pi(c[14] + math.radians(BASE_AXIS_ANGLE_DEG)*radius*z[14])
    return y


def propose(center, u, radius, bounds):
    return np.concatenate((
        perturb_side(center[:15], u[:15], radius, bounds["small"]),
        perturb_side(center[15:], u[15:], radius, bounds["large"]),
    ))


def build_fixed_final_thermo():
    report = load_json(THERMO_REPORT); best = report["best_feasible"]
    seed_machine, definition, mass, geom, _, basis = thermo5d.make_3952_basis()
    source = thermo5d.load_pair(RANK)
    machine = thermo5d.build_design(seed_machine, geom, mass, best["parameters"], source["small"], source["large"])
    identity = {
        "thermo5d_report": str(THERMO_REPORT),
        "thermo5d_report_sha256": hashlib.sha256(THERMO_REPORT.read_bytes()).hexdigest(),
        "thermo5d_candidate_id": best["candidate_id"],
        "pair_candidate_id": source["pair_id"],
        "fixed_thermo5d_parameters": best["parameters"],
        "reference_efficiency": best["result"]["indicated_thermal_efficiency"],
        "reference_power_W": best["result"]["indicated_power_w"],
        "fixed_total_gas_mass_kg": mass,
        "basis": basis,
    }
    return machine, definition, np.asarray(best["last_complete_state"], float), identity


def build_design(machine, smech, lmech):
    limits = machine.configuration.machine_volumes
    kin = IndependentSixBarVolumeKinematics(smech, lmech, limits.small_cylinder, limits.large_cylinder)
    return replace(machine, kinematics=kin)


def vector_from_record(r):
    return np.asarray([r["parameters"][n] for n in fine.PAIR_NAMES], float)


def eligible(r): return fine.eligible(r)
def eta(r): return fine.eta(r)


def parse_radii(s):
    out = tuple(float(x) for x in s.split(",") if x.strip())
    if not out or min(out) <= 0: raise argparse.ArgumentTypeError("positive radii required")
    return out


def save_best(path, rec):
    fine.save_best_pair(path, rec)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--evaluations", type=int, default=96)
    ap.add_argument("--evaluations-per-radius", type=int, default=24)
    ap.add_argument("--radius-multipliers", type=parse_radii, default=DEFAULT_RADII)
    ap.add_argument("--budget-seconds", type=float, default=1800)
    ap.add_argument("--candidate-seconds", type=float, default=120)
    ap.add_argument("--seed", type=int, default=3971001)
    ap.add_argument("--max-rod-angle-deg", type=float, default=8.0)
    ap.add_argument("--max-mean-rod-angle-deg", type=float, default=2.0)
    ap.add_argument("--max-structural-attempts", type=int, default=500)
    args = ap.parse_args()

    pair = load_json(PAIR_PATH); branches = branches_from_pair(pair); src = source_vector(pair)
    bounds = rod_bounds(src); motion = fine.load_motion_target()
    machine, definition, reference_state, fixed_identity = build_fixed_final_thermo()

    s0 = aligned_seed_side(src[:15], branches["small_primary"], branches["small_second"], motion["theta"])
    l0 = aligned_seed_side(src[15:], branches["large_primary"], branches["large_second"], motion["theta"])
    seed_vec = np.concatenate((s0, l0))

    # If the exactly aligned seed is too strict, cheaply repair only rod/offset/angle.
    try:
        screen_pair(seed_vec, branches, motion, bounds, args.max_rod_angle_deg, args.max_mean_rod_angle_deg)
        seed_kind = "direct_PCA_aligned_seed"
    except Exception:
        pts6 = qmc.Sobol(d=6, scramble=True, seed=args.seed+17).random_base2(12)
        best_seed = None; best_score = math.inf
        for u in pts6:
            x = seed_vec.copy()
            for j, off in enumerate((0,15)):
                side = "small" if j == 0 else "large"; lo, hi = bounds[side]
                x[off+12] = lo + u[3*j]*(hi-lo)
                x[off+13] += (2*u[3*j+1]-1)*1.0
                x[off+14] = fine.wrap_pi(x[off+14] + math.radians(12)*(2*u[3*j+2]-1))
            try:
                _, _, sm, lm, _, _ = screen_pair(x, branches, motion, bounds, args.max_rod_angle_deg, args.max_mean_rod_angle_deg)
            except Exception:
                continue
            score = sm["position_rms_full_cycle"]**2 + lm["position_rms_full_cycle"]**2
            if score < best_score: best_score, best_seed = score, x.copy()
        if best_seed is None:
            raise RuntimeError("No admissible compact aligned seed. Relax rod-angle limits.")
        seed_vec = best_seed; seed_kind = "repaired_axis_rod_seed"

    out = args.output_directory; out.mkdir(parents=True, exist_ok=True)
    hist_path = out/"history.jsonl"; report_path = out/"report.json"; best_path = out/"best_pair.json"; def_path = out/"definition.json"
    definition_payload = {
        "study": "final rank_01 thermo fixed; local 30-D six-bar search with compact aligned piston rods",
        "fixed_thermo": fixed_identity,
        "source_pair": str(PAIR_PATH),
        "branches": branches,
        "rod_constraints": {
            "small_rod_over_crank": list(bounds["small"]),
            "large_rod_over_crank": list(bounds["large"]),
            "max_abs_rod_axis_angle_deg": args.max_rod_angle_deg,
            "max_abs_mean_rod_axis_angle_deg": args.max_mean_rod_angle_deg,
            "source_small_rod": float(src[12]), "source_large_rod": float(src[27]),
        },
        "seed_kind": seed_kind,
        "search": {"evaluations": args.evaluations, "evaluations_per_radius": args.evaluations_per_radius,
                   "radius_multipliers": list(args.radius_multipliers), "sobol_seed": args.seed,
                   "free_dimensions": 30, "objective": "maximize indicated thermal efficiency"},
    }
    if def_path.exists() and load_json(def_path) != definition_payload:
        raise ValueError("Definition changed; use a new output directory.")
    if not def_path.exists(): def_path.write_text(json.dumps(definition_payload, indent=2)+"\n", encoding="utf-8")

    history = load_history(hist_path); good = [r for r in history if eligible(r)]; best = max(good, key=eta) if good else None
    points = qmc.Sobol(d=30, scramble=True, seed=args.seed).random_base2(16)
    skipped = 0; started = time.monotonic(); deadline = started + args.budget_seconds; done = 0

    while time.monotonic() < deadline and sum(r.get("kind") != "compact_aligned_seed" for r in history) < args.evaluations:
        proposal_index = max((r.get("proposal_index", -1) for r in history), default=-1)+1
        if not history:
            x = seed_vec.copy(); kind = "compact_aligned_seed"; ridx = -1; radius = 0.0
        else:
            n = sum(r.get("kind") != "compact_aligned_seed" for r in history)
            ridx = min(n//args.evaluations_per_radius, len(args.radius_multipliers)-1)
            wanted = float(args.radius_multipliers[ridx]); center = vector_from_record(best) if best is not None else seed_vec
            found = None
            for attempt in range(args.max_structural_attempts):
                u = points[(proposal_index + skipped + attempt) % len(points)]
                for backoff in fine.DIRECTION_BACKOFFS:
                    rr = wanted*backoff; trial = propose(center, u, rr, bounds)
                    try: screen_pair(trial, branches, motion, bounds, args.max_rod_angle_deg, args.max_mean_rod_angle_deg)
                    except Exception: continue
                    found = trial; radius = rr; skipped += attempt; break
                if found is not None: break
            if found is None: raise RuntimeError("No feasible compact 30-D direction after backoff.")
            x = found; kind = "compact_rods_30d_local_sobol"

        cid = fine.candidate_id(RANK, branches, x)
        if any(r.get("candidate_id") == cid and r.get("result", {}).get("status") != "interrupted" for r in history):
            skipped += 1; continue

        sv, lv, sm, lm, smech, lmech = screen_pair(x, branches, motion, bounds, args.max_rod_angle_deg, args.max_mean_rod_angle_deg)
        design = build_design(machine, smech, lmech)
        reusable = [r for r in history if r.get("result", {}).get("status") == "converged" and r.get("last_complete_state")]
        if reusable:
            warm = min(reusable, key=lambda r: fine.normalized_pair_distance(x, vector_from_record(r)))
            initial = np.asarray(warm["last_complete_state"], float); warm_source = warm["candidate_id"]
            warm_distance = fine.normalized_pair_distance(x, vector_from_record(warm))
        else:
            initial = reference_state; warm_source = fixed_identity["thermo5d_candidate_id"]; warm_distance = None

        thermo = fine.evaluate_thermo(label=f"rank01_compact_{proposal_index}", design=design, definition=definition,
                                      initial=initial, deadline=deadline, candidate_seconds=args.candidate_seconds)
        params = {n: float(v) for n, v in zip(fine.PAIR_NAMES, x)}
        sdict = fine.mechanism_candidate_dict(sv, branches["small_primary"], branches["small_second"], sm)
        ldict = fine.mechanism_candidate_dict(lv, branches["large_primary"], branches["large_second"], lm)
        combined = float(math.sqrt(0.5*(sm["position_rms_full_cycle"]**2 + lm["position_rms_full_cycle"]**2)))
        rec = {
            "index": len(history), "proposal_index": proposal_index, "family_rank": RANK, "candidate_id": cid,
            "kind": kind, "radius_index": ridx, "radius_multiplier": radius, "parameters": params, "branches": branches,
            "mechanical": {"small": sdict, "large": ldict, "combined_position_rms_vs_3952": combined},
            "result": thermo["compact"], "feasible": thermo["feasible"], "physical_constraint_failures": thermo["reasons"],
            "elapsed_seconds": thermo["elapsed_seconds"], "warm_start_source": warm_source, "warm_start_distance": warm_distance,
            "warm_start_mode": thermo["warm_mode"], "safe_retry_used": thermo["safe_retry"], "backend": thermo["backend"],
            "last_complete_state": thermo["state"].tolist() if thermo["state"] is not None else None,
        }
        append_jsonl(hist_path, rec); history.append(rec); done += 1
        moved = False
        if eligible(rec) and (best is None or eta(rec) > eta(best)):
            best = rec; moved = True; save_best(best_path, rec)
        e = rec["result"].get("indicated_thermal_efficiency")
        print(json.dumps({
            "index": rec["index"], "status": rec["result"].get("status"), "feasible": rec["feasible"], "efficiency": e,
            "loss_pp_vs_rank01": 100*(fixed_identity["reference_efficiency"]-e) if e is not None else None,
            "power_W": rec["result"].get("indicated_power_w"), "S_rod_R": float(sv[12]), "L_rod_R": float(lv[12]),
            "S_mean_deg": sm["rod_axis_mean_angle_deg"], "S_max_deg": sm["rod_axis_max_abs_angle_deg"],
            "L_mean_deg": lm["rod_axis_mean_angle_deg"], "L_max_deg": lm["rod_axis_max_abs_angle_deg"],
            "motion_rms_pct": 100*combined, "center_moved": moved,
            "best_efficiency": eta(best) if best is not None else None,
        }, separators=(",", ":")), flush=True)

    best = max((r for r in history if eligible(r)), key=eta, default=None)
    if best is not None: save_best(best_path, best)
    summary = {"description": definition_payload["study"], "definition": definition_payload, "attempted_total": len(history),
               "completed_this_run": done, "skipped_structural_total": skipped,
               "actual_duration_seconds": time.monotonic()-started, "best_feasible": best,
               "history": str(hist_path), "best_pair": str(best_path)}
    report_path.write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(f"Saved {report_path}")


if __name__ == "__main__":
    main()
