#!/usr/bin/env python3
"""Free-spline search seeded independently from the retained linear UU motion.

Hardware, gas inventory and thermal conditions are frozen to the final 8H
machine.  The known refined spline is diagnostic-only and never guides search.

Example:
PYTHONPATH=src python3 examples/optimize_motor_spline_from_linear_260k.py \
  --fourier-report outputs/motor_fourier_c2_8h_refine/report.json \
  --reference-spline-report outputs/motor_free_spline_260k_refine/report.json \
  --output-directory outputs/motor_spline_from_linear_260k \
  --controls 16 --evaluations 4096 --budget-seconds 21600
"""
from __future__ import annotations
from pathlib import Path
from types import SimpleNamespace
import argparse, hashlib, json, math, time
import numpy as np
from scipy.stats import qmc

from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
from dada_solver.four_stage_kinematics import FourStageVolumeKinematics
from dada_solver.integration import IntegrationInterrupted
from dada_solver.wall_backend import WallBackendSettings
from compare_motor_motion_laws_stage7A5 import ROOT, _evaluate
from optimize_motor_temperature_point import base_geometry, feasibility
from refine_motor_four_stage_hx9d_variable_gas import _load_basis
from refine_motor_fourier_c2_260k import build_design as fourier_build_design
from optimize_motor_free_spline_260k_v3 import (
    _admissible, _build_candidate_design, _candidate_id, _eta, _eligible,
    _fixed_inventory_design, _load_history, _make_kinematics,
    _normalized_motion, _record_distance, _safe_uniform_initial,
    _shape_diagnostics, _shape_rms, _source_best, _source_harmonics,
    _canonical,
)

POWER_FLOOR_W = 25.0
DEFAULT_FOURIER = ROOT/"outputs/motor_fourier_c2_8h_refine/report.json"
DEFAULT_REFERENCE = ROOT/"outputs/motor_free_spline_260k_refine/report.json"
DEFAULT_OUTPUT = ROOT/"outputs/motor_spline_from_linear_260k"
SHAPE_RADII = (0.10,0.070,0.048,0.032,0.021,0.013,0.008,0.0045,0.0025,0.0012)
PHASE_RADII = (7.0,5.0,3.5,2.5,1.7,1.1,0.7,0.4,0.22,0.10)
EVALS_PER_RADIUS = 160
WIDE_EVERY = 24
WIDE_SHAPE = 0.14
WIDE_PHASE = 9.0
DIRECT_WARM = 0.050

def resolve(path):
    return path if path.is_absolute() else ROOT/path

def linear_kin(limits,p):
    return FourStageVolumeKinematics(
        limits.small_cylinder, limits.large_cylinder,
        float(p["t1"]),float(p["t2"]),float(p["t3"]),
        float(p["a_l"]),float(p["b_l"]),float(p["a_s"]),float(p["b_s"]))

def sample_controls(kin,n):
    th=2*math.pi*np.arange(n,dtype=float)/n
    sl,ll=kin.small_volume_limits,kin.large_volume_limits
    s=np.array([(kin.small_cylinder_volume(float(a))-sl.minimum)/sl.swept for a in th])
    l=np.array([(kin.large_cylinder_volume(float(a))-ll.minimum)/ll.swept for a in th])
    return _canonical(s),_canonical(l)

def circular_smooth(v):
    v=np.asarray(v,float)
    return _canonical((np.roll(v,2)+4*np.roll(v,1)+6*v+4*np.roll(v,-1)+np.roll(v,-2))/16)

def make_seed(limits,s,l,max_passes=8):
    for passes in range(max_passes+1):
        kin=_make_kinematics(limits,s,l,0.0,0.0)
        if _admissible(kin): return s,l,kin,passes
        s,l=circular_smooth(s),circular_smooth(l)
    raise RuntimeError("Could not create admissible spline seed from linear law.")

def direction(u,center):
    d=2*np.asarray(u,float)-1
    d-=np.mean(d); d-=center*float(d@center)
    n=float(np.linalg.norm(d))
    if n<1e-12:
        d=np.roll(center,1)-center
        d-=np.mean(d); d-=center*float(d@center)
        n=float(np.linalg.norm(d))
    return d/n

def propose(center,u,r,pr):
    s0=np.asarray(center["small_controls"],float)
    l0=np.asarray(center["large_controls"],float); n=len(s0)
    ds,dl=direction(u[:n],s0),direction(u[n:2*n],l0)
    rs=r*(0.18+0.82*float(u[2*n])); rl=r*(0.18+0.82*float(u[2*n+1]))
    s,l=_canonical(s0+rs*ds),_canonical(l0+rl*dl)
    us=(float(u[0])+float(u[n]))%1
    ul=(float(u[n-1])+float(u[2*n-1]))%1
    ps=(float(center["small_phase_deg"])+pr*(2*us-1))%360
    pl=(float(center["large_phase_deg"])+pr*(2*ul-1))%360
    return s,l,ps,pl

def reference_spline(path,limits):
    if not path.exists(): return None,None
    r=json.loads(path.read_text()); b=r.get("best_feasible")
    if b is None:return None,None
    k=_make_kinematics(limits,np.array(b["small_controls"]),np.array(b["large_controls"]),
                       float(b["small_phase_deg"]),float(b["large_phase_deg"]))
    return b,k

def save_csv(path,linear,reference,candidate):
    t=np.linspace(0,1,2881)
    ql,dl=_normalized_motion(linear,t); qc,dc=_normalized_motion(candidate,t)
    cols=[t,360*t,ql[0],ql[1],qc[0],qc[1],dl[0],dl[1],dc[0],dc[1]]
    names=["motor_time_fraction","motor_angle_deg","linear_small","linear_large",
           "candidate_small","candidate_large","linear_dsmall_dt","linear_dlarge_dt",
           "candidate_dsmall_dt","candidate_dlarge_dt"]
    if reference is not None:
        qr,dr=_normalized_motion(reference,t)
        cols += [qr[0],qr[1],dr[0],dr[1]]
        names += ["reference_spline_small","reference_spline_large",
                  "reference_spline_dsmall_dt","reference_spline_dlarge_dt"]
    np.savetxt(path,np.column_stack(cols),delimiter=",",header=",".join(names),comments="")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--fourier-report",type=Path,default=DEFAULT_FOURIER)
    ap.add_argument("--reference-spline-report",type=Path,default=DEFAULT_REFERENCE)
    ap.add_argument("--output-directory",type=Path,default=DEFAULT_OUTPUT)
    ap.add_argument("--controls",type=int,default=16)
    ap.add_argument("--evaluations",type=int,default=4096)
    ap.add_argument("--budget-seconds",type=float,default=21600)
    ap.add_argument("--candidate-seconds",type=float,default=150)
    ap.add_argument("--evaluations-per-radius",type=int,default=EVALS_PER_RADIUS)
    ap.add_argument("--seed",type=int,default=26092273)
    args=ap.parse_args()
    if args.controls<8: raise ValueError("--controls must be >= 8")

    fp,rp=resolve(args.fourier_report),resolve(args.reference_spline_report)
    out=resolve(args.output_directory); out.mkdir(parents=True,exist_ok=True)
    fr=json.loads(fp.read_text()); fb=_source_best(fr); h=_source_harmonics(fb)
    lr=fr.get("reference_linear_UU")
    if not lr or not lr.get("parameters"):
        raise RuntimeError("Fourier report lacks reference_linear_UU.parameters")
    lp=dict(lr["parameters"])

    cd,base=_load_basis(); geometry=base_geometry(base)
    design=fourier_build_design(base,geometry,fb["thermo"],
        np.array(fb["small_coefficients"]),np.array(fb["large_coefficients"]),h)
    fstate=np.array(fb["last_complete_state"],float)
    total_mass=float(fb.get("result",{}).get("total_mass_kg",np.sum(fstate[:8:2])))
    design=_fixed_inventory_design(design,total_mass)
    limits=design.configuration.machine_volumes

    linear=linear_kin(limits,lp)
    ss,ls=sample_controls(linear,args.controls)
    ss,ls,seedkin,smooth_passes=make_seed(limits,ss,ls)
    seed_dist=_shape_rms(linear,seedkin,4096)
    refbest,refkin=reference_spline(rp,limits)

    definition=SimpleNamespace(
        wall_numerical_settings=cd.wall_numerical_settings,
        wall_backend=WallBackendSettings("numba"),
        exact_kinematics_cache=True,shared_replay=True)

    identity={
      "study":"260 K spline independently seeded from retained linear UU v1",
      "fourier_report":str(fp),"fourier_sha256":hashlib.sha256(fp.read_bytes()).hexdigest(),
      "reference_spline_report":str(rp) if rp.exists() else None,
      "reference_spline_sha256":hashlib.sha256(rp.read_bytes()).hexdigest() if rp.exists() else None,
      "reference_spline_is_diagnostic_only":True,
      "linear_parameters":lp,"controls_per_piston":args.controls,
      "linear_seed_smoothing_passes":smooth_passes,
      "linear_seed_distance_vs_exact_linear":seed_dist,
      "fixed_total_mass_kg":total_mass,"fixed_thermo":fb["thermo"],
      "shape_radii":list(SHAPE_RADII),"phase_radii_deg":list(PHASE_RADII),
      "evaluations_per_radius":args.evaluations_per_radius,
      "wide_probe_every":WIDE_EVERY,"wide_shape_radius":WIDE_SHAPE,
      "wide_phase_radius_deg":WIDE_PHASE,"sobol_seed":args.seed,
      "constraints":{"exactly_one_minimum_and_one_maximum_per_piston":True,
        "C2_periodic_position":True,"hardware_frozen_to_final_8H":True,
        "gas_inventory_frozen_to_final_8H":True,"four_stage_boundaries_imposed_after_seed":False},
    }
    dp=out/"definition.json"; hp=out/"history.jsonl"; reportp=out/"report.json"
    if dp.exists() and json.loads(dp.read_text())!=identity:
        raise RuntimeError("Campaign definition changed; use a new output directory.")
    if not dp.exists(): dp.write_text(json.dumps(identity,indent=2)+"\n")

    history=_load_history(hp)
    valid=[r for r in history if _eligible(r)]
    best=max(valid,key=_eta) if valid else None
    finished={r["candidate_id"] for r in history if r.get("result",{}).get("status")!="interrupted"}
    sobol=qmc.Sobol(d=2*args.controls+2,scramble=True,seed=args.seed).random_base2(16)
    seedrec={"small_controls":ss.tolist(),"large_controls":ls.tolist(),
             "small_phase_deg":0.0,"large_phase_deg":0.0}
    start=time.monotonic(); deadline=start+args.budget_seconds; completed=0
    pi=max((r.get("proposal_index",-1) for r in history),default=-1)+1

    while completed<args.evaluations and time.monotonic()<deadline:
        idx=pi; pi+=1
        if not history and idx==0:
            s,l,ps,pl=ss.copy(),ls.copy(),0.0,0.0
            kind="linear_smoothed_spline_seed"; ri=-1; radius=phase_radius=0.0
        else:
            searched=sum(1 for r in history if r.get("kind")!="linear_smoothed_spline_seed")
            ri=min(searched//args.evaluations_per_radius,len(SHAPE_RADII)-1)
            if searched>0 and searched%WIDE_EVERY==WIDE_EVERY-1:
                center=seedrec; radius=WIDE_SHAPE; phase_radius=WIDE_PHASE
                kind="linear_seed_wide_probe"
            else:
                center=best if best is not None else seedrec
                radius=SHAPE_RADII[ri]; phase_radius=PHASE_RADII[ri]
                kind="incumbent_spline_refinement"
            s,l,ps,pl=propose(center,sobol[idx%len(sobol)],radius,phase_radius)

        cid=_candidate_id(s,l,ps,pl)
        if cid in finished: continue
        cdesign=_build_candidate_design(design,s,l,ps,pl); ckin=cdesign.kinematics
        if not _admissible(ckin): continue

        before=time.monotonic(); cand_deadline=min(deadline,before+args.candidate_seconds)
        reusable=[r for r in history if r.get("result",{}).get("status")=="converged"
                  and r.get("last_complete_state") is not None]
        warm_distance=None
        if reusable:
            warm=min(reusable,key=lambda r:_record_distance(s,l,ps,pl,r))
            warm_distance=_record_distance(s,l,ps,pl,warm)
            warm_source=warm["candidate_id"]
            if warm_distance<=DIRECT_WARM:
                initial=np.array(warm["last_complete_state"],float); warm_mode="nearest_periodic_state"
            else:
                initial=_safe_uniform_initial(cdesign,np.array(warm["last_complete_state"],float))
                warm_mode="direct_safe_uniform_for_distant_shape"
        else:
            initial=_safe_uniform_initial(cdesign,fstate)
            warm_source="uniform_gas_with_8H_wall_energy"; warm_mode="safe_uniform_linear_seed"

        def progress(_):
            now=time.monotonic()
            if now>=deadline: raise IntegrationInterrupted("Campaign budget exhausted.")
            if now>=cand_deadline: raise IntegrationInterrupted("Candidate budget exhausted.")

        box={}
        def observe(periodic): box["statistics"]=periodic.backend_statistics
        safe_retry=False
        try:
            try:
                result,state=_evaluate(f"spline_from_linear_260k_{idx}",cdesign,definition,
                    initial_state=initial,progress_callback=progress,periodic_observer=observe)
            except MicrotubeDomainError:
                safe_retry=True
                safe=_safe_uniform_initial(cdesign,initial)
                result,state=_evaluate(f"spline_from_linear_260k_{idx}_safe",cdesign,definition,
                    initial_state=safe,progress_callback=progress,periodic_observer=observe)
                warm_mode="safe_uniform_retry_after_domain_error"
        except MicrotubeDomainError as exc:
            result={"status":"invalid_exchanger","message":str(exc)}; state=None
        except IntegrationInterrupted as exc:
            result={"status":"interrupted","message":str(exc)}; state=None
        except (ValueError,RuntimeError,ArithmeticError) as exc:
            result={"status":"integration_failure","message":str(exc)}; state=None

        feasible,reasons=(feasibility(result,POWER_FLOOR_W)
            if result.get("status")=="converged" else (False,[]))
        dlin=_shape_rms(linear,ckin,2048); dseed=_shape_rms(seedkin,ckin,2048)
        dref=_shape_rms(refkin,ckin,2048) if refkin is not None else None
        rec={
          "index":len(history),"proposal_index":idx,"candidate_id":cid,"kind":kind,
          "controls_per_piston":args.controls,"radius_index":ri,"shape_radius":radius,
          "phase_radius_deg":phase_radius,"small_controls":s.tolist(),"large_controls":l.tolist(),
          "small_phase_deg":float(ps)%360,"large_phase_deg":float(pl)%360,
          "shape_diagnostics":_shape_diagnostics(ckin),
          "shape_distance_vs_exact_linear":dlin,
          "shape_distance_vs_linear_spline_seed":dseed,
          "shape_distance_vs_known_spline_diagnostic_only":dref,
          "result":result,"feasible":feasible,"physical_constraint_failures":reasons,
          "elapsed_seconds":time.monotonic()-before,"warm_start_source":warm_source,
          "warm_start_mode":warm_mode,"warm_start_distance":warm_distance,
          "backend_statistics":box.get("statistics"),"safe_retry_used":safe_retry,
          "last_complete_state":state.tolist() if state is not None else None}
        with hp.open("a") as f:f.write(json.dumps(rec)+"\n");f.flush()
        history.append(rec);completed+=1
        if result.get("status")!="interrupted":finished.add(cid)
        moved=False
        if feasible and (best is None or _eta(rec)>_eta(best)):
            best=rec;moved=True;save_csv(out/"best_motion.csv",linear,refkin,ckin)

        report={
          "best_feasible":best,
          "linear_seed":{"parameters":lp,"smoothing_passes":smooth_passes,
                         "shape_diagnostics":_shape_diagnostics(seedkin),
                         "shape_distance_vs_exact_linear":seed_dist},
          "reference_linear_UU_original_machine":lr,
          "fixed_machine_fourier_8H":{"candidate_id":fb.get("candidate_id"),
            "efficiency":fb["result"]["indicated_thermal_efficiency"],
            "power_W":fb["result"]["indicated_power_w"],"total_mass_kg":total_mass,
            "thermo":fb["thermo"]},
          "known_refined_spline_diagnostic_only":(
            {"candidate_id":refbest["candidate_id"],
             "efficiency":refbest["result"]["indicated_thermal_efficiency"],
             "power_W":refbest["result"]["indicated_power_w"],
             "shape_diagnostics":refbest["shape_diagnostics"]} if refbest else None),
          "attempted_total":len(history),"completed_this_run":completed,
          "requested_duration_seconds":args.budget_seconds,
          "actual_duration_seconds":time.monotonic()-start,"definition":identity}
        reportp.write_text(json.dumps(report,indent=2)+"\n")
        print(json.dumps({"index":rec["index"],"proposal_index":idx,"kind":kind,
          "radius":radius,"phase_radius_deg":phase_radius,"status":result.get("status"),
          "feasible":feasible,"efficiency":result.get("indicated_thermal_efficiency"),
          "power_W":result.get("indicated_power_w"),
          "rms_vs_linear":dlin["combined_position_rms"],
          "rms_vs_known_spline":dref["combined_position_rms"] if dref else None,
          "center_moved":moved,"best_efficiency":_eta(best) if best else None,
          "actual_backend":(box.get("statistics") or {}).get("actual_backend"),
          "elapsed_seconds":rec["elapsed_seconds"]}),flush=True)
    print(f"Saved {reportp}",flush=True)

if __name__=="__main__":
    main()
