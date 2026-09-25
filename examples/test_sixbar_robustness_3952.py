#!/usr/bin/env python3
"""Manufacturing-tolerance robustness test for the four final 3952 six-bar pairs.

This is deliberately a screening robustness study, not an industrial reliability
certification.

Question
--------
If the final paired six-bar mechanisms are manufactured with small independent
dimensional errors, do they remain mechanically valid and retain their
thermodynamic performance without any re-optimization?

Families
--------
1, 4, 12 and 50 use the final geometry from:
    outputs/sixbar_thermo_coupled_3952_fine_hlat25/rank_XX/best_pair.json

Each family keeps its own final thermo-5D setting from:
    outputs/sixbar_thermo5d_3952/report.json

No perturbed sample is repaired or optimized.

Two-stage workflow
------------------
1. Cheap mechanical screen on many deterministic low-discrepancy samples.
2. Thermodynamic replay on a reproducible random subset of mechanically valid
   samples.

Default tolerance levels
------------------------
The scalar linear tolerance tau is interpreted as a bounded manufacturing
envelope.

    tau = 0.001  -> +/-0.1% link-length error
    tau = 0.005  -> +/-0.5%
    tau = 0.010  -> +/-1.0%

Angular indexing/axis error scales independently but proportionally:
    +/-0.05 deg, +/-0.25 deg, +/-0.50 deg respectively.

Perturbation model
------------------
For each SMALL and LARGE mechanism independently:

- AD, BC, CD, EF, GF and piston rod:
    multiplicative uniform error +/- tau.
- E attachment coordinates:
    additive error +/- tau * max(BC, 1 crank radius) on each local coordinate.
- second fixed pivot Gx, Gy:
    additive error +/- tau crank radius on each coordinate.
- H attachment ratios:
    additive error +/- tau on H_along/EF and H_normal/EF.
- slider-axis offset:
    additive error +/- tau crank radius.
- primary phase and slider-axis angle:
    bounded angular error; default 0.5 deg at tau=1%.

AB=1 remains the dimensional datum and assembly branches remain fixed.
All four families receive the same normalized 30-D Sobol perturbations.

Mechanical acceptance
---------------------
A perturbed pair must satisfy, on both S and L sides:

- full-cycle closure through the production SixBarCylinderMechanism;
- stroke/crank in [1, 3];
- primary transmission sine >= 0.30;
- secondary transmission sine >= 0.30;
- rod/slider-axis cosine >= 0.95;
- EH/crank <= 7;
- crank-axis clearance >= 0.50 crank radius;
- H lateral RMS/stroke <= 0.25;
- H lateral span/stroke <= 0.65;
- exactly two piston velocity zero crossings.

Thermodynamic acceptance
------------------------
The mechanism pair is inserted into the family's own final thermo-5D machine.
No thermal parameter is re-tuned. The normal solver validity criteria and 25 W
indicated-power floor apply.

Typical first pass
------------------
    PYTHONPATH=src python3 examples/test_sixbar_robustness_3952.py

Screen only:
    PYTHONPATH=src python3 examples/test_sixbar_robustness_3952.py --stage mechanical

A deeper thermo replay can resume the existing screening:
    PYTHONPATH=src python3 examples/test_sixbar_robustness_3952.py \
      --stage thermo --thermo-samples 128

Outputs
-------
    outputs/sixbar_robustness_3952/
        definition.json
        mechanical_samples.jsonl
        thermo_samples.jsonl
        report.json
        summary.csv
"""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
import argparse
import csv
import json
import math
import time

import numpy as np
from scipy.stats import qmc

import optimize_sixbar_pairs_thermo_3952_fine as fine
import optimize_sixbar_pairs_thermo5d_3952_hlat25 as thermo5d

ROOT = Path.cwd()
DEFAULT_OUTPUT = ROOT / "outputs" / "sixbar_robustness_3952"
FINAL_PAIR_ROOT = ROOT / "outputs" / "sixbar_thermo_coupled_3952_fine_hlat25"
THERMO_REPORT = ROOT / "outputs" / "sixbar_thermo5d_3952" / "report.json"
FAMILIES = (1, 4, 12, 50)

STROKE_FLOOR = 1.0
STROKE_CEILING = 3.0
PRIMARY_SINE_FLOOR = 0.30
SECONDARY_SINE_FLOOR = 0.30
ROD_COS_FLOOR = 0.95
MAXIMUM_EH = 7.0
CRANK_CLEARANCE_FLOOR = 0.50
H_LATERAL_RMS_MAX = 0.25
H_LATERAL_SPAN_MAX = 0.65
POWER_FLOOR_W = 25.0

PAIR_NAMES = tuple(f"S_{x}" for x in fine.NAMES) + tuple(f"L_{x}" for x in fine.NAMES)


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def load_jsonl(path: Path):
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def append_jsonl(path: Path, row: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(row, separators=(",", ":")) + "\n")
        f.flush()


def wrap_pi(a: float) -> float:
    return (float(a) + math.pi) % (2.0 * math.pi) - math.pi


def parse_float_list(text: str):
    vals = tuple(float(x) for x in text.split(",") if x.strip())
    if not vals or any(x <= 0 for x in vals):
        raise argparse.ArgumentTypeError("Expected positive comma-separated floats.")
    return vals


def sobol_points(n: int, seed: int) -> np.ndarray:
    if n <= 0:
        raise ValueError("mechanical sample count must be positive")
    m = int(math.ceil(math.log2(n)))
    pts = qmc.Sobol(d=30, scramble=True, seed=seed).random_base2(m)
    return np.asarray(pts[:n], float)


def candidate_from_pair_side(side: dict) -> np.ndarray:
    p = side["parameters"]
    return np.asarray([float(p[name]) for name in fine.NAMES], dtype=float)


def load_family(rank: int, thermo_report: dict):
    pair_path = FINAL_PAIR_ROOT / f"rank_{rank:02d}" / "best_pair.json"
    report_path = FINAL_PAIR_ROOT / f"rank_{rank:02d}" / "report.json"
    pair = load_json(pair_path)
    mech_report = load_json(report_path)
    mech_best = mech_report.get("best_feasible")
    if mech_best is None:
        raise RuntimeError(f"Family {rank}: missing mechanical champion.")
    if pair.get("candidate_id") != mech_best.get("candidate_id"):
        raise RuntimeError(f"Family {rank}: best_pair/report candidate IDs disagree.")

    t = thermo_report["families"][str(rank)]["best_feasible"]
    if t["pair_candidate_id"] != pair["candidate_id"]:
        raise RuntimeError(f"Family {rank}: thermo-5D champion refers to another pair.")
    if t.get("last_complete_state") is None:
        raise RuntimeError(f"Family {rank}: no final thermo periodic state.")

    center = np.concatenate((candidate_from_pair_side(pair["small"]), candidate_from_pair_side(pair["large"])))
    branches = {
        "small_primary": int(pair["small"]["primary"]["assembly_branch"]),
        "small_second": int(pair["small"]["second_branch"]),
        "large_primary": int(pair["large"]["primary"]["assembly_branch"]),
        "large_second": int(pair["large"]["second_branch"]),
    }
    return {
        "rank": rank,
        "pair_path": pair_path,
        "pair_candidate_id": pair["candidate_id"],
        "center": center,
        "branches": branches,
        "thermo_parameters": dict(t["parameters"]),
        "thermo_state": np.asarray(t["last_complete_state"], float),
        "nominal_efficiency": float(t["result"]["indicated_thermal_efficiency"]),
        "nominal_power_W": float(t["result"]["indicated_power_w"]),
    }


def perturb_side(center, u, *, tau, angle_deg_at_1pct):
    c = np.asarray(center, float)
    z = 2.0 * np.asarray(u, float) - 1.0
    y = c.copy()

    for i in (0, 1, 2):
        y[i] = max(1e-6, c[i] * (1.0 + tau * z[i]))

    e_scale = max(abs(c[1]), 1.0)
    y[3] = c[3] + tau * e_scale * z[3]
    y[4] = c[4] + tau * e_scale * z[4]

    angle_span = math.radians(angle_deg_at_1pct * tau / 0.01)
    y[5] = wrap_pi(c[5] + angle_span * z[5])

    y[6] = c[6] + tau * z[6]
    y[7] = c[7] + tau * z[7]

    for i in (8, 9):
        y[i] = max(1e-6, c[i] * (1.0 + tau * z[i]))

    y[10] = c[10] + tau * z[10]
    y[11] = c[11] + tau * z[11]
    y[12] = max(1e-6, c[12] * (1.0 + tau * z[12]))
    y[13] = c[13] + tau * z[13]
    y[14] = wrap_pi(c[14] + angle_span * z[14])
    return y


def perturb_pair(center, u, *, tau, angle_deg_at_1pct):
    return np.concatenate((
        perturb_side(center[:15], u[:15], tau=tau, angle_deg_at_1pct=angle_deg_at_1pct),
        perturb_side(center[15:], u[15:], tau=tau, angle_deg_at_1pct=angle_deg_at_1pct),
    ))


def raw_side_metrics(v, *, primary_branch, second_branch, side, motion):
    if not fine.primary_closure_ok(v):
        raise ValueError("primary_closure")
    pv = np.asarray(v[:6], float)
    downstream = np.asarray(v[6:], float)
    theta = motion["theta"]
    try:
        E, Ed, psine, pdata = fine.stage.primary(theta, pv, primary_branch)
    except ValueError as exc:
        raise ValueError(f"primary_geometry:{exc}") from exc
    try:
        return fine.v41.evaluate_v41(
            downstream,
            E=E,
            Ed=Ed,
            psine=psine,
            pdata=pdata,
            pv=pv,
            second_branch=second_branch,
            tq=motion[f"{side}_q"],
            tdq=motion[f"{side}_dq"],
            mask=motion["masks"][side],
            velocity_weight=0.0005,
        )
    except ValueError as exc:
        raise ValueError(f"downstream_geometry:{exc}") from exc


def constraint_failures(m):
    f = []
    if m["stroke_over_crank"] < STROKE_FLOOR: f.append("stroke_low")
    if m["stroke_over_crank"] > STROKE_CEILING: f.append("stroke_high")
    if m["minimum_primary_transmission_sine"] < PRIMARY_SINE_FLOOR: f.append("primary_transmission")
    if m["minimum_secondary_transmission_sine"] < SECONDARY_SINE_FLOOR: f.append("secondary_transmission")
    if m["minimum_rod_axis_cosine"] < ROD_COS_FLOOR: f.append("rod_angle")
    if m["EH_over_crank"] > MAXIMUM_EH: f.append("EH")
    if m["crank_axis_to_EFH_clearance_over_crank"] < CRANK_CLEARANCE_FLOOR: f.append("crank_clearance")
    if m["H_axis_lateral_rms_over_stroke"] > H_LATERAL_RMS_MAX: f.append("H_lateral_rms")
    if m["H_axis_lateral_span_over_stroke"] > H_LATERAL_SPAN_MAX: f.append("H_lateral_span")
    if int(m["zero_crossing_count"]) != 2: f.append("motion_topology")
    return f


def compact_mechanical(m):
    keys = (
        "position_rms_full_cycle", "stroke_over_crank",
        "minimum_primary_transmission_sine", "minimum_secondary_transmission_sine",
        "minimum_rod_axis_cosine", "EH_over_crank",
        "H_axis_lateral_rms_over_stroke", "H_axis_lateral_span_over_stroke",
        "crank_axis_to_EFH_clearance_over_crank", "zero_crossing_count", "extra_reversals",
    )
    return {k: m[k] for k in keys}


def screen_pair(x, *, branches, motion):
    sv, lv = np.asarray(x[:15], float).copy(), np.asarray(x[15:], float).copy()
    sv[5], sv[14] = wrap_pi(sv[5]), wrap_pi(sv[14])
    lv[5], lv[14] = wrap_pi(lv[5]), wrap_pi(lv[14])
    try:
        sm = raw_side_metrics(sv, primary_branch=branches["small_primary"], second_branch=branches["small_second"], side="small", motion=motion)
    except ValueError as exc:
        return False, [f"small:{exc}"], None, None
    try:
        lm = raw_side_metrics(lv, primary_branch=branches["large_primary"], second_branch=branches["large_second"], side="large", motion=motion)
    except ValueError as exc:
        return False, [f"large:{exc}"], compact_mechanical(sm), None

    failures = [f"small:{q}" for q in constraint_failures(sm)] + [f"large:{q}" for q in constraint_failures(lm)]
    if not failures:
        try:
            fine.mechanism_from_vector(sv, branches["small_primary"], branches["small_second"])
        except ValueError as exc:
            failures.append(f"small:production_geometry:{exc}")
        try:
            fine.mechanism_from_vector(lv, branches["large_primary"], branches["large_second"])
        except ValueError as exc:
            failures.append(f"large:production_geometry:{exc}")
    return not failures, failures, compact_mechanical(sm), compact_mechanical(lm)


def sample_key(rank, tau, sample_index):
    return f"{rank}:{tau:.12g}:{sample_index}"


def percentile_dict(values):
    a = np.asarray(list(values), float)
    if a.size == 0:
        return None
    return {
        "n": int(a.size), "min": float(np.min(a)), "p05": float(np.percentile(a, 5)),
        "p25": float(np.percentile(a, 25)), "median": float(np.percentile(a, 50)),
        "p75": float(np.percentile(a, 75)), "p95": float(np.percentile(a, 95)),
        "max": float(np.max(a)), "mean": float(np.mean(a)), "std": float(np.std(a)),
    }


def stable_subset(indices, n, *, seed):
    indices = np.asarray(sorted(indices), dtype=int)
    if len(indices) <= n:
        return indices.tolist()
    rng = np.random.default_rng(seed)
    return sorted(map(int, rng.choice(indices, size=n, replace=False)))


def write_report(output, definition, mechanical_rows, thermo_rows, families):
    groups = {}
    for rank in definition["families"]:
        fam = families[rank]
        for tau in definition["tolerances"]:
            mk = [r for r in mechanical_rows if r["family_rank"] == rank and math.isclose(r["tolerance"], tau, rel_tol=0, abs_tol=1e-15)]
            tk = [r for r in thermo_rows if r["family_rank"] == rank and math.isclose(r["tolerance"], tau, rel_tol=0, abs_tol=1e-15)]
            valid = [r for r in mk if r["mechanically_valid"]]
            failures = Counter(q for r in mk for q in r.get("failures", []))
            mech = {
                "total": len(mk), "valid": len(valid),
                "valid_fraction": len(valid)/len(mk) if mk else None,
                "failure_counts": dict(failures.most_common()),
            }
            if valid:
                for side in ("small", "large"):
                    mech[side] = {
                        "stroke_over_crank": percentile_dict(r[side]["stroke_over_crank"] for r in valid),
                        "primary_sine": percentile_dict(r[side]["minimum_primary_transmission_sine"] for r in valid),
                        "secondary_sine": percentile_dict(r[side]["minimum_secondary_transmission_sine"] for r in valid),
                        "rod_cosine": percentile_dict(r[side]["minimum_rod_axis_cosine"] for r in valid),
                        "H_lateral_rms": percentile_dict(r[side]["H_axis_lateral_rms_over_stroke"] for r in valid),
                        "H_lateral_span": percentile_dict(r[side]["H_axis_lateral_span_over_stroke"] for r in valid),
                    }
            converged = [r for r in tk if r["result"].get("status") == "converged"]
            feasible = [r for r in converged if r.get("thermo_feasible")]
            etas = [r["result"]["indicated_thermal_efficiency"] for r in feasible if r["result"].get("indicated_thermal_efficiency") is not None]
            powers = [r["result"]["indicated_power_w"] for r in feasible if r["result"].get("indicated_power_w") is not None]
            losses = [100*(fam["nominal_efficiency"]-e) for e in etas]
            thermo = {
                "selected": len(tk), "converged": len(converged), "feasible": len(feasible),
                "converged_fraction": len(converged)/len(tk) if tk else None,
                "feasible_fraction": len(feasible)/len(tk) if tk else None,
                "efficiency": percentile_dict(etas), "power_W": percentile_dict(powers),
                "efficiency_loss_percentage_points": percentile_dict(losses),
                "within_0p1pp_of_nominal_fraction": (sum(x <= 0.1 for x in losses)/len(losses) if losses else None),
                "nominal_efficiency": fam["nominal_efficiency"], "nominal_power_W": fam["nominal_power_W"],
            }
            groups[f"{rank}:{tau:.12g}"] = {"family_rank": rank, "tolerance": tau, "mechanical": mech, "thermodynamic": thermo}

    report = {"description": "Bounded manufacturing-tolerance robustness screen for final candidate-3952 paired six-bar families.", "definition": definition, "groups": groups}
    (output/"report.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    with (output/"summary.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["family_rank","tolerance","mechanical_samples","mechanical_valid","mechanical_valid_fraction","thermo_selected","thermo_converged","thermo_feasible","eta_nominal","eta_median","eta_p05","eta_min","eta_loss_pp_median","eta_loss_pp_p95","power_nominal_W","power_median_W","power_p05_W"])
        for g in groups.values():
            m,t = g["mechanical"],g["thermodynamic"]
            e,l,p = t["efficiency"] or {}, t["efficiency_loss_percentage_points"] or {}, t["power_W"] or {}
            w.writerow([g["family_rank"],g["tolerance"],m["total"],m["valid"],m["valid_fraction"],t["selected"],t["converged"],t["feasible"],t["nominal_efficiency"],e.get("median"),e.get("p05"),e.get("min"),l.get("median"),l.get("p95"),t["nominal_power_W"],p.get("median"),p.get("p05")])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", choices=("mechanical","thermo","all"), default="all")
    ap.add_argument("--families", default="1,4,12,50")
    ap.add_argument("--tolerances", type=parse_float_list, default=(0.001,0.005,0.01))
    ap.add_argument("--angle-deg-at-1pct", type=float, default=0.5)
    ap.add_argument("--mechanical-samples", type=int, default=2048)
    ap.add_argument("--thermo-samples", type=int, default=64)
    ap.add_argument("--seed", type=int, default=3952701)
    ap.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--budget-seconds", type=float, default=21600.0)
    ap.add_argument("--candidate-seconds", type=float, default=180.0)
    args = ap.parse_args()

    ranks = tuple(int(x) for x in args.families.split(",") if x.strip())
    if any(r not in FAMILIES for r in ranks): raise ValueError(f"Allowed families are {FAMILIES}.")
    if args.mechanical_samples <= 0 or args.thermo_samples <= 0: raise ValueError("sample counts must be positive")
    if args.angle_deg_at_1pct <= 0: raise ValueError("angle-deg-at-1pct must be positive")

    output = args.output_directory
    output.mkdir(parents=True, exist_ok=True)
    mechanical_path, thermo_path, definition_path = output/"mechanical_samples.jsonl", output/"thermo_samples.jsonl", output/"definition.json"
    thermo_report = load_json(THERMO_REPORT)
    families = {rank: load_family(rank, thermo_report) for rank in ranks}

    definition = {
        "study": "sixbar manufacturing-tolerance robustness around final 3952 pairs",
        "families": list(ranks), "tolerances": list(args.tolerances),
        "angle_deg_at_1pct": args.angle_deg_at_1pct,
        "mechanical_samples_per_family_tolerance": args.mechanical_samples,
        "sampling": "30-D scrambled Sobol; identical normalized perturbations for all families/tolerances",
        "distribution": "independent bounded uniform coordinates in [-1,+1]",
        "perturbation_model": {
            "AD_BC_CD_EF_GF_piston_rod": "multiplicative +/- tau",
            "E_coordinates": "additive +/- tau*max(BC,1)",
            "G_coordinates": "additive +/- tau crank radius",
            "H_ratios": "additive +/- tau",
            "slider_axis_offset": "additive +/- tau crank radius",
            "phase_and_axis_angle": f"+/- {args.angle_deg_at_1pct} deg at tau=0.01, linear with tau",
            "crank_AB": "fixed dimensional datum = 1", "assembly_branches": "fixed",
        },
        "mechanical_constraints": {
            "stroke_over_crank": [STROKE_FLOOR,STROKE_CEILING],
            "primary_sine_minimum": PRIMARY_SINE_FLOOR, "secondary_sine_minimum": SECONDARY_SINE_FLOOR,
            "rod_axis_cosine_minimum": ROD_COS_FLOOR, "EH_over_crank_maximum": MAXIMUM_EH,
            "crank_axis_clearance_minimum": CRANK_CLEARANCE_FLOOR,
            "H_lateral_rms_over_stroke_maximum": H_LATERAL_RMS_MAX,
            "H_lateral_span_over_stroke_maximum": H_LATERAL_SPAN_MAX,
            "piston_zero_crossings_required": 2,
        },
        "thermo_policy": {"family_specific_final_thermo5d_machine": True, "no_reoptimization_after_perturbation": True, "power_floor_W": POWER_FLOOR_W, "source_report": str(THERMO_REPORT)},
        "seed": args.seed,
        "family_sources": {str(rank): {"pair_path": str(families[rank]["pair_path"]), "pair_candidate_id": families[rank]["pair_candidate_id"], "nominal_efficiency": families[rank]["nominal_efficiency"], "nominal_power_W": families[rank]["nominal_power_W"], "thermo_parameters": families[rank]["thermo_parameters"]} for rank in ranks},
    }
    if definition_path.exists():
        if load_json(definition_path) != definition: raise ValueError("Robustness definition changed; use a new output directory.")
    else:
        definition_path.write_text(json.dumps(definition, indent=2)+"\n", encoding="utf-8")

    motion = fine.load_motion_target()
    common_u = sobol_points(args.mechanical_samples, args.seed)
    deadline = time.monotonic() + args.budget_seconds
    mechanical_rows = load_jsonl(mechanical_path)
    mechanical_done = {sample_key(int(r["family_rank"]),float(r["tolerance"]),int(r["sample_index"])) for r in mechanical_rows}

    if args.stage in ("mechanical","all"):
        for tau in args.tolerances:
            for sample_index,u in enumerate(common_u):
                for rank in ranks:
                    if time.monotonic() >= deadline: break
                    key = sample_key(rank,tau,sample_index)
                    if key in mechanical_done: continue
                    fam = families[rank]
                    x = perturb_pair(fam["center"],u,tau=tau,angle_deg_at_1pct=args.angle_deg_at_1pct)
                    valid,failures,sm,lm = screen_pair(x,branches=fam["branches"],motion=motion)
                    row = {"family_rank":rank,"tolerance":tau,"sample_index":sample_index,"sample_key":key,"mechanically_valid":bool(valid),"failures":failures,"small":sm,"large":lm,"parameters":{name:float(value) for name,value in zip(PAIR_NAMES,x)}}
                    append_jsonl(mechanical_path,row); mechanical_done.add(key)
                if time.monotonic() >= deadline: break
            if time.monotonic() >= deadline: break

    mechanical_rows = load_jsonl(mechanical_path)

    if args.stage in ("thermo","all") and time.monotonic() < deadline:
        seed_machine, solver_definition, total_mass, geometry_info, _ts, _id = thermo5d.make_3952_basis()
        thermo_rows = load_jsonl(thermo_path)
        thermo_done = {sample_key(int(r["family_rank"]),float(r["tolerance"]),int(r["sample_index"])) for r in thermo_rows}
        valid_by_group = defaultdict(list)
        for r in mechanical_rows:
            if r["mechanically_valid"]: valid_by_group[(int(r["family_rank"]),float(r["tolerance"]))].append(r)

        for tau_i,tau in enumerate(args.tolerances):
            for rank in ranks:
                if time.monotonic() >= deadline: break
                candidates = valid_by_group.get((rank,float(tau)),[])
                chosen = stable_subset([int(r["sample_index"]) for r in candidates],args.thermo_samples,seed=args.seed+rank*1009+tau_i*100003)
                by_index = {int(r["sample_index"]):r for r in candidates}
                fam = families[rank]
                for sample_index in chosen:
                    if time.monotonic() >= deadline: break
                    key = sample_key(rank,tau,sample_index)
                    if key in thermo_done: continue
                    mrow = by_index[sample_index]
                    x = np.asarray([mrow["parameters"][name] for name in PAIR_NAMES],float)
                    small_mech = fine.mechanism_from_vector(x[:15],fam["branches"]["small_primary"],fam["branches"]["small_second"])
                    large_mech = fine.mechanism_from_vector(x[15:],fam["branches"]["large_primary"],fam["branches"]["large_second"])
                    design = thermo5d.build_design(seed_machine,geometry_info,total_mass,fam["thermo_parameters"],small_mech,large_mech)
                    ev = thermo5d.evaluate(label=f"robust3952_r{rank}_tau{tau:g}_i{sample_index}",design=design,definition=solver_definition,initial=fam["thermo_state"],deadline=deadline,candidate_seconds=args.candidate_seconds)
                    row = {"family_rank":rank,"tolerance":tau,"sample_index":sample_index,"sample_key":key,"result":ev["compact"],"thermo_feasible":bool(ev["feasible"]),"physical_constraint_failures":ev["reasons"],"safe_retry_used":ev["safe_retry"],"warm_start_mode":"family_nominal_final_thermo_state","elapsed_seconds":ev["elapsed_seconds"],"backend":ev["backend"]}
                    append_jsonl(thermo_path,row); thermo_done.add(key)
                    print(json.dumps({"family":rank,"tau":tau,"sample":sample_index,"status":row["result"].get("status"),"feasible":row["thermo_feasible"],"eta":row["result"].get("indicated_thermal_efficiency"),"power_W":row["result"].get("indicated_power_w")},separators=(",",":")),flush=True)
            if time.monotonic() >= deadline: break

    write_report(output,definition,load_jsonl(mechanical_path),load_jsonl(thermo_path),families)
    print(f"Saved {output/'report.json'}")
    print(f"Saved {output/'summary.csv'}")


if __name__ == "__main__":
    main()
