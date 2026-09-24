#!/usr/bin/env python3
"""Optimize downstream dyad + H + piston slider for a panel of fixed L primaries.

Candidate 3952 LARGE side.  The primary A-B-C-D+E is frozen.  For each selected
primary, optimize the 9 downstream variables and both secondary branches.

The objective is position-first outside the automatically detected HP-noise
band.  H transverse motion is allowed, but directly bounded relative to piston
stroke so the finite connecting rod may help without excessive "gigotement".

Default primary panel: ranks 0,1,5,10,13 from large_primary_cadence_3952.json.
"""

from __future__ import annotations
import argparse, json, math, time
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution

import search_large_sixbar_stageL1 as stage
import search_motor_hybrid_3952_sixbar as motor3952

ROOT = Path.cwd()
PRIMARY_LIBRARY = ROOT/"outputs/large_primary_cadence_3952.json"
TARGET = ROOT/"outputs/motor_hybrid_c2_15p_260k/candidate_3952_motion_target.csv"
OUTPUT = ROOT/"outputs/large_downstream_primary_panel_3952.json"

DNAMES = [
    "second_pivot_x","second_pivot_y","link_EF","link_GF",
    "H_along_over_EF","H_normal_over_EF","piston_rod",
    "slider_axis_offset","slider_axis_angle",
]

def load_panel(path, ranks):
    d=json.loads(path.read_text())
    by={int(x["library_rank"]):x for x in d["candidates"]}
    missing=[r for r in ranks if r not in by]
    if missing: raise ValueError(f"Missing primary ranks: {missing}")
    return [by[r] for r in ranks]

def pvector(item):
    p=item["primary"]
    return np.array([p["ground"],p["coupler"],p["rocker"],
                     p["E_along"],p["E_normal"],p["phase"]],float)

def bounds_for(E):
    c=np.mean(E,axis=0)
    return [
        (float(c[0])-10,float(c[0])+10),
        (float(c[1])-10,float(c[1])+10),
        (.35,8.0),(.50,20.0),(-2.0,3.0),(-2.5,2.5),
        (.75,12.0),(-16.0,16.0),(-math.pi,math.pi),
    ]

def closure_penalty(E,gx,gy,lef,lgf):
    d=np.linalg.norm(np.array((gx,gy))-E,axis=1)
    low=abs(lef-lgf); high=lef+lgf; scale=max(high,1.)
    return (max(0.,low-float(np.min(d)))+
            max(0.,float(np.max(d))-high))/scale

def evaluate(x, *, E, Ed, psine, pdata, pv, second_branch,
             tq, tdq, mask, velocity_weight):
    gx,gy,lef,lgf,ha,hn,rod,axis_offset,axis_angle=map(float,x)
    G=np.array((gx,gy))
    delta=G-E
    dist=np.linalg.norm(delta,axis=1)
    if np.any(dist<=1e-10): raise ValueError("secondary centers")
    u0=delta/dist[:,None]
    along=(lef*lef-lgf*lgf+dist*dist)/(2*dist)
    h2=lef*lef-along*along
    if np.any(h2<=1e-10): raise ValueError("secondary closure")
    h=np.sqrt(h2)
    n0=np.column_stack((-u0[:,1],u0[:,0]))
    F=E+along[:,None]*u0+second_branch*h[:,None]*n0

    ef=F-E; gf=F-G
    det=ef[:,0]*gf[:,1]-ef[:,1]*gf[:,0]
    if np.any(np.abs(det)<=1e-10): raise ValueError("secondary singular")
    rhs=np.sum(ef*Ed,axis=1)
    Fd=np.column_stack((rhs*gf[:,1]/det,-rhs*gf[:,0]/det))
    ssine=float(np.min(np.abs(det)/(lef*lgf)))

    u=ef/lef; ud=(Fd-Ed)/lef
    n=np.column_stack((-u[:,1],u[:,0]))
    nd=np.column_stack((-ud[:,1],ud[:,0]))
    H=E+ha*lef*u+hn*lef*n
    Hd=Ed+ha*lef*ud+hn*lef*nd

    axis=np.array((math.cos(axis_angle),math.sin(axis_angle)))
    normal=np.array((-axis[1],axis[0]))
    origin=axis_offset*normal
    rel=H-origin
    longitudinal=rel@axis; transverse=rel@normal
    longitudinal_d=Hd@axis; transverse_d=Hd@normal
    margin2=rod*rod-transverse*transverse
    if np.any(margin2<=1e-10): raise ValueError("rod closure")
    margin=np.sqrt(margin2)
    slider=longitudinal+margin
    slider_d=longitudinal_d-transverse*transverse_d/margin
    stroke=float(np.ptp(slider))
    if stroke<=1e-10: raise ValueError("zero stroke")

    q=1-(slider-float(np.min(slider)))/stroke
    dq=-slider_d/stroke
    m=np.asarray(mask,bool)
    pos=float(np.sqrt(np.mean((q[m]-tq[m])**2)))
    der=float(np.sqrt(np.mean((dq[m]-tdq[m])**2)))
    raw_pos=float(np.sqrt(np.mean((q-tq)**2)))
    raw_der=float(np.sqrt(np.mean((dq-tdq)**2)))
    score=math.sqrt(pos*pos+velocity_weight*der*der)

    rod_cos=float(np.min(margin/rod))
    clearance,inside=stage.signed_triangle_clearance_origin(E,F,H)
    line=stage.best_line_metrics(H,stroke,axis)

    tc=transverse-float(np.mean(transverse))
    lat_rms=float(np.sqrt(np.mean(tc*tc)))/stroke
    lat_span=float(np.ptp(transverse))/stroke
    eh=lef*math.hypot(ha,hn)

    full=np.concatenate((pv,np.asarray(x,float)))
    out={
        "score":score,"motion_score":score,
        "position_rms":pos,"derivative_rms":der,
        "position_rms_raw":raw_pos,"derivative_rms_raw":raw_der,
        "stroke_over_crank":stroke,
        "minimum_primary_transmission_sine":float(psine),
        "minimum_secondary_transmission_sine":ssine,
        "minimum_rod_axis_cosine":rod_cos,
        "EF_over_crank":lef,"GF_over_crank":lgf,
        "EH_over_EF":math.hypot(ha,hn),"EH_over_crank":eh,
        "piston_rod_over_crank":rod,"piston_rod_over_stroke":rod/stroke,
        "H_axis_lateral_rms_over_stroke":lat_rms,
        "H_axis_lateral_span_over_stroke":lat_span,
        "crank_axis_to_EFH_clearance_over_crank":clearance,
        "crank_axis_inside_EFH_frames":inside,
        "primary":dict(pdata),
        "parameters":dict(zip(stage.NAMES,map(float,full))),
        "second_branch":int(second_branch),
    }
    out.update(line)
    return out

def violation(r,args):
    v=0.
    v+=max(0.,args.stroke_floor-r["stroke_over_crank"])/args.stroke_floor
    v+=max(0.,r["stroke_over_crank"]-args.stroke_ceiling)/args.stroke_ceiling
    v+=max(0.,args.secondary_sine_floor-r["minimum_secondary_transmission_sine"])/args.secondary_sine_floor
    v+=max(0.,args.rod_cos_floor-r["minimum_rod_axis_cosine"])/max(1.-args.rod_cos_floor,1e-3)
    v+=max(0.,r["EH_over_crank"]-args.maximum_EH)/args.maximum_EH
    v+=max(0.,args.crank_clearance_floor-r["crank_axis_to_EFH_clearance_over_crank"])/max(args.crank_clearance_floor,.1)
    v+=max(0.,r["piston_rod_over_stroke"]-args.rod_stroke_ceiling)/args.rod_stroke_ceiling
    v+=max(0.,r["H_axis_lateral_rms_over_stroke"]-args.h_lateral_rms_max)/args.h_lateral_rms_max
    v+=max(0.,r["H_axis_lateral_span_over_stroke"]-args.h_lateral_span_max)/args.h_lateral_span_max
    return float(v)

def rod_reach_penalty(x,E,second_branch):
    """Continuous fallback after secondary closure, for DE guidance."""
    gx,gy,lef,lgf,ha,hn,rod,axis_offset,axis_angle=map(float,x)
    G=np.array((gx,gy)); delta=G-E; dist=np.linalg.norm(delta,axis=1)
    u0=delta/dist[:,None]
    along=(lef*lef-lgf*lgf+dist*dist)/(2*dist)
    h2=lef*lef-along*along
    if np.any(h2<=0): return 1000.
    h=np.sqrt(h2); n0=np.column_stack((-u0[:,1],u0[:,0]))
    F=E+along[:,None]*u0+second_branch*h[:,None]*n0
    u=(F-E)/lef; n=np.column_stack((-u[:,1],u[:,0]))
    H=E+ha*lef*u+hn*lef*n
    axis=np.array((math.cos(axis_angle),math.sin(axis_angle)))
    normal=np.array((-axis[1],axis[0]))
    trans=(H-axis_offset*normal)@normal
    req=float(np.max(np.abs(trans)))
    if rod<req: return 10.+100.*(req-rod)/max(rod,1.)
    return 1000.

def dense_check(x, **kw):
    args=kw.pop("args")
    try: r=evaluate(x,**kw)
    except (ValueError,FloatingPointError,OverflowError): return None
    if violation(r,args)>0: return None
    r["feasible"]=True
    return r

def save(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2)+"\n")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--primary-library",type=Path,default=PRIMARY_LIBRARY)
    ap.add_argument("--target",type=Path,default=TARGET)
    ap.add_argument("--output",type=Path,default=OUTPUT)
    ap.add_argument("--primary-ranks",default="0,1,5,10,13")
    ap.add_argument("--restarts",type=int,default=2)
    ap.add_argument("--generations",type=int,default=280)
    ap.add_argument("--population-size",type=int,default=12)
    ap.add_argument("--seed",type=int,default=3952300)
    ap.add_argument("--coarse-stride",type=int,default=4)
    ap.add_argument("--noise-fronts",type=int,default=6)
    ap.add_argument("--velocity-weight",type=float,default=.001)
    ap.add_argument("--stroke-floor",type=float,default=1.)
    ap.add_argument("--stroke-ceiling",type=float,default=3.)
    ap.add_argument("--secondary-sine-floor",type=float,default=.30)
    ap.add_argument("--rod-cos-floor",type=float,default=.95)
    ap.add_argument("--maximum-EH",type=float,default=7.)
    ap.add_argument("--crank-clearance-floor",type=float,default=.50)
    ap.add_argument("--rod-stroke-ceiling",type=float,default=3.)
    ap.add_argument("--h-lateral-rms-max",type=float,default=.20)
    ap.add_argument("--h-lateral-span-max",type=float,default=.65)
    ap.add_argument("--archive-per-run",type=int,default=30)
    args=ap.parse_args()

    ranks=[int(x) for x in args.primary_ranks.split(",") if x.strip()]
    panel=load_panel(args.primary_library,ranks)

    raw=np.genfromtxt(args.target,delimiter=",",names=True)
    if raw["theta_rad"][-1]>=2*math.pi-1e-9: raw=raw[:-1]
    band=motor3952._detect_noise_bands(args.target,args.noise_fronts)["large"]

    dth=np.asarray(raw["theta_rad"],float)
    dq=np.asarray(raw["large_fraction_0_1"],float)
    ddq=np.asarray(raw["large_dq_dtheta_per_rad"],float)
    dmask=~motor3952._theta_band_mask(dth,band)

    coarse=raw[::args.coarse_stride]
    th=np.asarray(coarse["theta_rad"],float)
    tq=np.asarray(coarse["large_fraction_0_1"],float)
    tdq=np.asarray(coarse["large_dq_dtheta_per_rad"],float)
    mask=~motor3952._theta_band_mask(th,band)

    wall0=time.time(); cpu0=time.process_time()
    report={
        "description":"3952 LARGE downstream synthesis over fixed cadence-matched primaries",
        "length_unit":"crank_radius","primary_library":str(args.primary_library),
        "primary_ranks":ranks,"target":str(args.target),"noise_band":band,
        "objective":{"position_rms":"outside automatic HP noise band",
                     "velocity_weight":args.velocity_weight},
        "constraints":{
            "stroke_over_crank":[args.stroke_floor,args.stroke_ceiling],
            "secondary_sine_minimum":args.secondary_sine_floor,
            "rod_axis_cosine_minimum":args.rod_cos_floor,
            "EH_over_crank_maximum":args.maximum_EH,
            "crank_axis_clearance_minimum":args.crank_clearance_floor,
            "piston_rod_over_stroke_maximum":args.rod_stroke_ceiling,
            "H_axis_lateral_rms_over_stroke_maximum":args.h_lateral_rms_max,
            "H_axis_lateral_span_over_stroke_maximum":args.h_lateral_span_max,
        },
        "families":[],"best_by_primary_rank":{},"best":None,
    }

    global_best=[]
    for pi,item in enumerate(panel):
        rank=int(item["library_rank"]); pbranch=int(item["branch"]); pv=pvector(item)
        E,Ed,psine,pdata=stage.primary(th,pv,pbranch)
        dE,dEd,dpsine,dpdata=stage.primary(dth,pv,pbranch)
        bounds=bounds_for(E)
        fam={"primary_rank":rank,"primary_cadence_score":float(item["score"]),
             "primary_branch":pbranch,"primary":pdata,
             "minimum_primary_transmission_sine":float(psine),
             "downstream_bounds":dict(zip(DNAMES,bounds)),"runs":[],"best":None}
        print(f"\nPRIMARY rank={rank} branch={pbranch:+d} cadence={item['score']:.5f} sine={psine:.4f}",flush=True)

        for sbranch in (+1,-1):
            for restart in range(args.restarts):
                seed=args.seed+1000*pi+100*(sbranch<0)+restart
                archive=[]
                def obj(x):
                    x=np.asarray(x,float)
                    cp=closure_penalty(E,*map(float,x[:4]))
                    if cp>0: return 10.+100.*cp
                    try:
                        r=evaluate(x,E=E,Ed=Ed,psine=psine,pdata=pdata,pv=pv,
                                   second_branch=sbranch,tq=tq,tdq=tdq,mask=mask,
                                   velocity_weight=args.velocity_weight)
                    except (ValueError,FloatingPointError,OverflowError):
                        return rod_reach_penalty(x,E,sbranch)
                    v=violation(r,args)
                    if v>0: return 10.+100.*v
                    archive.append((float(r["score"]),x.copy()))
                    return float(r["score"])

                t0=time.time(); p0=time.process_time()
                opt=differential_evolution(
                    obj,bounds,seed=seed,popsize=args.population_size,
                    maxiter=args.generations,polish=True,tol=1e-8,
                    updating="immediate",workers=1,init="latinhypercube")

                vec=[np.asarray(opt.x,float)]
                if getattr(opt,"population",None) is not None:
                    order=np.argsort(np.asarray(opt.population_energies))
                    vec += [np.asarray(opt.population[int(j)],float)
                            for j in order[:args.archive_per_run]]
                vec += [x for _,x in sorted(archive,key=lambda z:z[0])[:args.archive_per_run]]
                dense=[]; seen=set()
                for x in vec:
                    k=tuple(round(float(v),10) for v in x)
                    if k in seen: continue
                    seen.add(k)
                    r=dense_check(
                        x,E=dE,Ed=dEd,psine=dpsine,pdata=dpdata,pv=pv,
                        second_branch=sbranch,tq=dq,tdq=ddq,mask=dmask,
                        velocity_weight=args.velocity_weight,args=args)
                    if r:
                        r["primary_library_rank"]=rank
                        r["primary_cadence_score"]=float(item["score"])
                        dense.append(r)
                best=min(dense,key=lambda z:z["score"]) if dense else None
                fam["runs"].append({
                    "restart":restart,"seed":seed,"second_branch":sbranch,
                    "optimizer_success":bool(opt.success),
                    "optimizer_message":str(opt.message),
                    "function_evaluations":int(opt.nfev),
                    "wall_seconds":time.time()-t0,"cpu_seconds":time.process_time()-p0,
                    "best_dense":best})
                fc=[r["best_dense"] for r in fam["runs"] if r["best_dense"]]
                fam["best"]=min(fc,key=lambda z:z["score"]) if fc else None

                if best:
                    print(f"  second={sbranch:+d} restart={restart}: "
                          f"pos={100*best['position_rms']:.4f}% raw={100*best['position_rms_raw']:.4f}% "
                          f"stroke/r={best['stroke_over_crank']:.3f} "
                          f"rod/stroke={best['piston_rod_over_stroke']:.3f} "
                          f"HlatRMS={best['H_axis_lateral_rms_over_stroke']:.3f} "
                          f"HlatSpan={best['H_axis_lateral_span_over_stroke']:.3f} "
                          f"sec-sine={best['minimum_secondary_transmission_sine']:.3f}",flush=True)
                else:
                    print(f"  second={sbranch:+d} restart={restart}: no dense feasible candidate",flush=True)

                temp=report["families"]+[fam]
                allb=[f["best"] for f in temp if f.get("best")]
                report["best"]=min(allb,key=lambda z:z["score"]) if allb else None
                report["best_by_primary_rank"]={str(f["primary_rank"]):f.get("best") for f in temp}
                report["elapsed_wall_hours"]=(time.time()-wall0)/3600
                report["elapsed_cpu_hours"]=(time.process_time()-cpu0)/3600
                old=report["families"]; report["families"]=temp; save(args.output,report); report["families"]=old

        report["families"].append(fam)
        if fam["best"]: global_best.append(fam["best"])
        report["best_by_primary_rank"]={str(f["primary_rank"]):f.get("best") for f in report["families"]}
        report["best"]=min(global_best,key=lambda z:z["score"]) if global_best else None
        report["elapsed_wall_hours"]=(time.time()-wall0)/3600
        report["elapsed_cpu_hours"]=(time.process_time()-cpu0)/3600
        save(args.output,report)

    print("\nDOWNSTREAM PANEL SUMMARY")
    for f in report["families"]:
        b=f["best"]
        if not b:
            print(f"rank {f['primary_rank']:2d}: no feasible downstream")
        else:
            print(f"rank {f['primary_rank']:2d}: pos={100*b['position_rms']:.5f}% "
                  f"raw={100*b['position_rms_raw']:.5f}% second={b['second_branch']:+d} "
                  f"stroke/r={b['stroke_over_crank']:.3f} "
                  f"rod/stroke={b['piston_rod_over_stroke']:.3f} "
                  f"HlatRMS={b['H_axis_lateral_rms_over_stroke']:.3f}")
    b=report["best"]
    if b:
        print("\nGLOBAL BEST")
        print(f"primary rank={b['primary_library_rank']} second={b['second_branch']:+d}\n"
              f"position RMS={100*b['position_rms']:.6f}%\n"
              f"full-cycle position RMS={100*b['position_rms_raw']:.6f}%\n"
              f"stroke/r={b['stroke_over_crank']:.6f}\n"
              f"rod/stroke={b['piston_rod_over_stroke']:.6f}\n"
              f"H lateral RMS/stroke={b['H_axis_lateral_rms_over_stroke']:.6f}\n"
              f"H lateral span/stroke={b['H_axis_lateral_span_over_stroke']:.6f}\n"
              f"primary sine={b['minimum_primary_transmission_sine']:.6f}\n"
              f"secondary sine={b['minimum_secondary_transmission_sine']:.6f}")
    save(args.output,report)
    print(f"Wrote {args.output}")

if __name__=="__main__":
    main()
