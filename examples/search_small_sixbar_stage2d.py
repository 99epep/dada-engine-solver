#!/usr/bin/env python3
"""Stage 2D: optimize a compact six-bar while exploiting H's natural near-line motion.

Compared with Stage 2C:
- the FINAL piston law remains the main objective;
- H is evaluated against its own best-fit straight line (PCA);
- the cylinder axis is constrained to align with that natural H line;
- optional crank-axis clearance is enforced for the moving EFH plate.

Typical first run:

PYTHONPATH=src python3 examples/search_small_sixbar_stage2d.py \
  --iterations 350 \
  --stroke-floor 1.5 \
  --h-line-rms-max 0.05 \
  --h-axis-angle-max 20 \
  --crank-clearance-floor 0.20 \
  --warm-start outputs/small_sixbar_stage2c_lat080.json \
  --output outputs/small_sixbar_stage2d_align20.json
"""

from __future__ import annotations
import argparse, json, math, time
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution

ROOT = Path.cwd()
TARGET = ROOT / "outputs" / "motor_champion_motion_target.csv"

BASE = dict(
    ground=2.045193376592419,
    coupler=2.1152136168766678,
    rocker=1.8202134329709354,
    e_along=1.873315626010214,
    e_normal=-1.5362965718510164,
    phase=-3.0863675667472057,
)
PRIMARY_BRANCH = -1

NAMES = [
    "primary_ground","primary_coupler","primary_rocker",
    "primary_E_along","primary_E_normal","primary_phase_delta",
    "second_pivot_x","second_pivot_y","link_EF","link_GF",
    "H_along_over_EF","H_normal_over_EF",
    "piston_rod","slider_axis_offset","slider_axis_angle",
]

def around(v, f):
    a, b = v*(1-f), v*(1+f)
    return (min(a,b), max(a,b))

BOUNDS = [
    around(BASE["ground"], .25),
    around(BASE["coupler"], .25),
    around(BASE["rocker"], .25),
    around(BASE["e_along"], .30),
    around(BASE["e_normal"], .30),
    (-math.pi/4, math.pi/4),
    (-8.,8.), (-8.,8.),
    (.4,5.), (.4,5.),
    (-2.,2.), (-2.,2.),
    (2.,30.), (-10.,10.), (-math.pi,math.pi),
]

def primary(theta, v):
    g,c,r,ea,en,pd = map(float, v[:6])
    phase = BASE["phase"] + pd

    if abs(g-1.) <= abs(c-r) + 1e-10:
        raise ValueError("primary inner closure")
    if c+r <= g+1. + 1e-10:
        raise ValueError("primary outer closure")

    a = theta + phase
    B = np.column_stack((np.cos(a), np.sin(a)))
    Bd = np.column_stack((-np.sin(a), np.cos(a)))
    D = np.array((g,0.))

    delta = D-B
    dist = np.linalg.norm(delta,axis=1)
    u0 = delta/dist[:,None]
    along = (c*c-r*r+dist*dist)/(2*dist)
    h2 = c*c-along*along
    if np.any(h2 <= 1e-10):
        raise ValueError("primary toggle")
    h = np.sqrt(h2)
    n0 = np.column_stack((-u0[:,1],u0[:,0]))
    C = B + along[:,None]*u0 + PRIMARY_BRANCH*h[:,None]*n0

    bc = C-B
    dc = C-D
    det = bc[:,0]*dc[:,1]-bc[:,1]*dc[:,0]
    if np.any(np.abs(det) <= 1e-10):
        raise ValueError("primary singular")

    rhs = np.sum(bc*Bd,axis=1)
    Cd = np.column_stack((rhs*dc[:,1]/det, -rhs*dc[:,0]/det))

    u = bc/c
    ud = (Cd-Bd)/c
    n = np.column_stack((-u[:,1],u[:,0]))
    nd = np.column_stack((-ud[:,1],ud[:,0]))

    E = B + ea*u + en*n
    Ed = Bd + ea*ud + en*nd
    sine = float(np.min(np.abs(det)/(c*r)))

    pdata = dict(
        ground=g,coupler=c,rocker=r,assembly_branch=PRIMARY_BRANCH,
        E_along=ea,E_normal=en,phase=phase,phase_delta=pd,
    )
    return E,Ed,sine,pdata


def best_line_metrics(H, stroke, cylinder_axis):
    center = np.mean(H, axis=0)
    X = H-center

    # Principal direction of H trajectory.
    cov = X.T @ X
    values, vectors = np.linalg.eigh(cov)
    line_dir = vectors[:, int(np.argmax(values))]
    line_dir = line_dir / np.linalg.norm(line_dir)

    if float(line_dir @ cylinder_axis) < 0:
        line_dir = -line_dir

    line_normal = np.array((-line_dir[1], line_dir[0]))
    along = X @ line_dir
    perp = X @ line_normal

    rms = float(np.sqrt(np.mean(perp*perp))) / stroke
    max_abs = float(np.max(np.abs(perp))) / stroke
    full_span = float(np.ptp(perp)) / stroke
    along_span = float(np.ptp(along)) / stroke

    c = abs(float(line_dir @ cylinder_axis))
    c = min(1.0, max(0.0, c))
    angle_deg = math.degrees(math.acos(c))

    return dict(
        H_best_line_rms_over_stroke=rms,
        H_best_line_max_over_stroke=max_abs,
        H_best_line_full_span_over_stroke=full_span,
        H_best_line_along_span_over_stroke=along_span,
        H_best_line_axis_angle_deg=angle_deg,
        H_best_line_angle_deg=math.degrees(math.atan2(line_dir[1],line_dir[0])),
    )


def segment_distance_origin_batch(A, B):
    AB = B-A
    den = np.sum(AB*AB,axis=1)
    t = -np.sum(A*AB,axis=1) / np.maximum(den,1e-30)
    t = np.clip(t,0.0,1.0)
    Q = A + t[:,None]*AB
    return np.linalg.norm(Q,axis=1)


def signed_triangle_clearance_origin(E,F,H):
    # Distance from origin A=(0,0) to EFH plate boundary.
    d = np.minimum.reduce([
        segment_distance_origin_batch(E,F),
        segment_distance_origin_batch(F,H),
        segment_distance_origin_batch(H,E),
    ])

    # Origin inside triangle? Use consistent edge-cross signs.
    def cross_to_origin(A,B):
        AB=B-A
        AO=-A
        return AB[:,0]*AO[:,1]-AB[:,1]*AO[:,0]

    s1=cross_to_origin(E,F)
    s2=cross_to_origin(F,H)
    s3=cross_to_origin(H,E)
    eps=1e-12
    inside = ((s1>=-eps)&(s2>=-eps)&(s3>=-eps)) | ((s1<=eps)&(s2<=eps)&(s3<=eps))
    signed=np.where(inside,-d,d)
    return float(np.min(signed)), int(np.sum(inside))


def evaluate(v, second_branch, theta, tq, tdq, max_eh):
    E,Ed,primary_sine,pdata = primary(theta,v)
    gx,gy,lef,lgf,ha,hn,rod,axis_offset,axis_angle = map(float,v[6:])

    eh_ratio = math.hypot(ha,hn)
    eh = lef*eh_ratio
    if eh > max_eh:
        raise ValueError("EH too large")

    G = np.array((gx,gy))
    delta = G-E
    dist = np.linalg.norm(delta,axis=1)
    if np.any(dist <= 1e-10):
        raise ValueError("secondary coincident")

    u0 = delta/dist[:,None]
    along = (lef*lef-lgf*lgf+dist*dist)/(2*dist)
    h2 = lef*lef-along*along
    if np.any(h2 <= 1e-10):
        raise ValueError("secondary closure")
    h = np.sqrt(h2)
    n0 = np.column_stack((-u0[:,1],u0[:,0]))
    F = E + along[:,None]*u0 + second_branch*h[:,None]*n0

    ef=F-E
    gf=F-G
    det=ef[:,0]*gf[:,1]-ef[:,1]*gf[:,0]
    if np.any(np.abs(det)<=1e-10):
        raise ValueError("secondary singular")

    rhs=np.sum(ef*Ed,axis=1)
    Fd=np.column_stack((rhs*gf[:,1]/det,-rhs*gf[:,0]/det))
    secondary_sine=float(np.min(np.abs(det)/(lef*lgf)))

    u=ef/lef
    ud=(Fd-Ed)/lef
    n=np.column_stack((-u[:,1],u[:,0]))
    nd=np.column_stack((-ud[:,1],ud[:,0]))
    H=E+ha*lef*u+hn*lef*n
    Hd=Ed+ha*lef*ud+hn*lef*nd

    axis=np.array((math.cos(axis_angle),math.sin(axis_angle)))
    normal=np.array((-axis[1],axis[0]))
    origin=axis_offset*normal
    rel=H-origin
    longitudinal=rel@axis
    transverse=rel@normal
    longitudinal_d=Hd@axis
    transverse_d=Hd@normal

    margin2=rod*rod-transverse*transverse
    if np.any(margin2<=1e-10):
        raise ValueError("rod closure")
    margin=np.sqrt(margin2)

    slider=longitudinal+margin
    slider_d=longitudinal_d-transverse*transverse_d/margin
    stroke=float(np.ptp(slider))
    if stroke<=1e-10:
        raise ValueError("zero stroke")

    q=1.-(slider-float(np.min(slider)))/stroke
    dq=-slider_d/stroke

    pos=float(np.sqrt(np.mean((q-tq)**2)))
    der=float(np.sqrt(np.mean((dq-tdq)**2)))
    motion=math.sqrt(pos*pos+.1*der*der)

    rod_cos=float(np.min(margin/rod))
    line=best_line_metrics(H,stroke,axis)
    crank_clearance,inside_count=signed_triangle_clearance_origin(E,F,H)

    out=dict(
        score=motion,
        motion_score=motion,
        position_rms=pos,derivative_rms=der,
        stroke_over_crank=stroke,
        minimum_primary_transmission_sine=primary_sine,
        minimum_secondary_transmission_sine=secondary_sine,
        minimum_rod_axis_cosine=rod_cos,
        EF_over_crank=lef,GF_over_crank=lgf,
        EH_over_EF=eh_ratio,EH_over_crank=eh,
        piston_rod_over_crank=rod,
        piston_rod_over_stroke=rod/stroke,
        crank_axis_to_EFH_clearance_over_crank=crank_clearance,
        crank_axis_inside_EFH_frames=inside_count,
        primary=pdata,
        parameters=dict(zip(NAMES,map(float,v))),
        second_branch=second_branch,
    )
    out.update(line)
    return out


def extract_best(data):
    if data.get("best") is not None:
        return data["best"]
    vals=[r.get("best_dense") for r in data.get("runs",[]) if r.get("best_dense")]
    return min(vals,key=lambda x:x["score"]) if vals else None


def warm_start(path):
    if path is None:
        return None
    data=json.loads(path.read_text())
    b=extract_best(data)
    if b is None:
        raise ValueError("warm-start has no candidate")
    p=b["parameters"]

    if "primary_ground" in p:
        x=np.array([p[n] for n in NAMES],float)
    else:
        x=np.array([
            BASE["ground"],BASE["coupler"],BASE["rocker"],
            BASE["e_along"],BASE["e_normal"],0.,
            p["second_pivot_x"],p["second_pivot_y"],
            p["link_EF"],p["link_GF"],
            p["H_along_over_EF"],p["H_normal_over_EF"],
            p["piston_rod"],p["slider_axis_offset"],p["slider_axis_angle"],
        ],float)

    lo=np.array([b[0] for b in BOUNDS])
    hi=np.array([b[1] for b in BOUNDS])
    return np.clip(x,lo,hi)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--iterations",type=int,default=350)
    ap.add_argument("--population-size",type=int,default=12)
    ap.add_argument("--restarts",type=int,default=1)
    ap.add_argument("--seed",type=int,default=2400)

    ap.add_argument("--stroke-floor",type=float,default=1.5)
    ap.add_argument("--primary-sine-floor",type=float,default=.30)
    ap.add_argument("--secondary-sine-floor",type=float,default=.30)
    ap.add_argument("--rod-cos-floor",type=float,default=.95)
    ap.add_argument("--maximum-EH",type=float,default=6.)

    ap.add_argument("--h-line-rms-max",type=float,default=.05)
    ap.add_argument("--h-axis-angle-max",type=float,default=20.)
    ap.add_argument("--crank-clearance-floor",type=float,default=.20)

    ap.add_argument("--warm-start",type=Path)
    ap.add_argument("--output",type=Path,default=ROOT/"outputs/small_sixbar_stage2d.json")
    args=ap.parse_args()

    raw=np.genfromtxt(TARGET,delimiter=",",names=True)
    if raw["theta_rad"][-1] >= 2*math.pi-1e-9:
        raw=raw[:-1]

    coarse=raw[::12]
    theta=coarse["theta_rad"]
    tq=coarse["small_fraction_0_1"]
    tdq=coarse["small_dq_dtheta_per_rad"]
    dense_theta=raw["theta_rad"]
    dense_q=raw["small_fraction_0_1"]
    dense_dq=raw["small_dq_dtheta_per_rad"]

    x0=warm_start(args.warm_start)

    report=dict(
        description="Stage 2D: final piston fit + intrinsic H straightness + H-line/cylinder alignment + crank clearance",
        length_unit="crank_radius",
        primary_branch=PRIMARY_BRANCH,
        base_primary=BASE,
        bounds=dict(zip(NAMES,BOUNDS)),
        constraints=dict(
            stroke_over_crank_minimum=args.stroke_floor,
            primary_sine_minimum=args.primary_sine_floor,
            secondary_sine_minimum=args.secondary_sine_floor,
            rod_axis_cosine_minimum=args.rod_cos_floor,
            EH_over_crank_maximum=args.maximum_EH,
            H_best_line_rms_over_stroke_maximum=args.h_line_rms_max,
            H_best_line_axis_angle_deg_maximum=args.h_axis_angle_max,
            crank_axis_to_EFH_clearance_over_crank_minimum=args.crank_clearance_floor,
        ),
        warm_start=str(args.warm_start) if args.warm_start else None,
        runs=[],
    )

    start=time.monotonic()
    run_id=0

    for restart in range(args.restarts):
        for branch in (+1,-1):
            archive=[]
            seed=args.seed+run_id
            run_id+=1

            def obj(v):
                try:
                    r=evaluate(v,branch,theta,tq,tdq,args.maximum_EH)
                except (ValueError,FloatingPointError):
                    return 1e4

                # Scale each violation by its own allowed magnitude so the
                # optimizer sees comparable pressure from all constraints.
                violation=0.
                violation += max(0., args.stroke_floor-r["stroke_over_crank"]) / max(args.stroke_floor,1e-9)
                violation += max(0., args.primary_sine_floor-r["minimum_primary_transmission_sine"]) / max(args.primary_sine_floor,1e-9)
                violation += max(0., args.secondary_sine_floor-r["minimum_secondary_transmission_sine"]) / max(args.secondary_sine_floor,1e-9)
                violation += max(0., args.rod_cos_floor-r["minimum_rod_axis_cosine"]) / max(1.-args.rod_cos_floor,1e-3)
                violation += max(0., r["H_best_line_rms_over_stroke"]-args.h_line_rms_max) / max(args.h_line_rms_max,1e-9)
                violation += max(0., r["H_best_line_axis_angle_deg"]-args.h_axis_angle_max) / max(args.h_axis_angle_max,1e-9)
                violation += max(0., args.crank_clearance_floor-r["crank_axis_to_EFH_clearance_over_crank"]) / max(args.crank_clearance_floor,0.1)

                if violation>0:
                    return 10.+100.*violation

                if not archive or r["score"]<archive[-1]["score"]:
                    archive.append(r)
                return r["score"]

            kw={}
            if x0 is not None:
                kw["x0"]=x0

            opt=differential_evolution(
                obj,BOUNDS,seed=seed,popsize=args.population_size,
                maxiter=args.iterations,polish=True,tol=1e-8,
                updating="immediate",workers=1,**kw
            )

            dense=[]
            for c in sorted(archive,key=lambda x:x["score"])[:50]:
                v=np.array([c["parameters"][n] for n in NAMES])
                try:
                    d=evaluate(v,branch,dense_theta,dense_q,dense_dq,args.maximum_EH)
                except (ValueError,FloatingPointError):
                    continue

                ok=(
                    d["stroke_over_crank"]>=args.stroke_floor and
                    d["minimum_primary_transmission_sine"]>=args.primary_sine_floor and
                    d["minimum_secondary_transmission_sine"]>=args.secondary_sine_floor and
                    d["minimum_rod_axis_cosine"]>=args.rod_cos_floor and
                    d["H_best_line_rms_over_stroke"]<=args.h_line_rms_max and
                    d["H_best_line_axis_angle_deg"]<=args.h_axis_angle_max and
                    d["crank_axis_to_EFH_clearance_over_crank"]>=args.crank_clearance_floor
                )
                if ok:
                    d["feasible"]=True
                    dense.append(d)

            best=min(dense,key=lambda x:x["score"]) if dense else None
            report["runs"].append(dict(
                restart=restart,seed=seed,second_branch=branch,
                optimizer_success=bool(opt.success),
                optimizer_message=str(opt.message),
                function_evaluations=int(opt.nfev),
                best_dense=best,
            ))

            vals=[r["best_dense"] for r in report["runs"] if r["best_dense"]]
            report["best"]=min(vals,key=lambda x:x["score"]) if vals else None
            report["elapsed_seconds"]=time.monotonic()-start

            args.output.parent.mkdir(parents=True,exist_ok=True)
            args.output.write_text(json.dumps(report,indent=2)+"\n")

            if best:
                print(
                    f"restart {restart}, branch {branch:+d}: "
                    f"pos={100*best['position_rms']:.4f}% "
                    f"dRMS={best['derivative_rms']:.6f} "
                    f"stroke/r={best['stroke_over_crank']:.3f} "
                    f"HlineRMS={best['H_best_line_rms_over_stroke']:.3f} "
                    f"Haxis={best['H_best_line_axis_angle_deg']:.2f}deg "
                    f"clearA={best['crank_axis_to_EFH_clearance_over_crank']:.3f} "
                    f"EF={best['EF_over_crank']:.3f} "
                    f"GF={best['GF_over_crank']:.3f} "
                    f"EH={best['EH_over_crank']:.3f}",
                    flush=True,
                )
            else:
                print(f"restart {restart}, branch {branch:+d}: no dense feasible candidate",flush=True)

    print()
    b=report.get("best")
    if b:
        print("BEST DENSE CANDIDATE")
        print(f"position RMS = {100*b['position_rms']:.6f}%")
        print(f"derivative RMS = {b['derivative_rms']:.9f} rad^-1")
        print(f"stroke / crank = {b['stroke_over_crank']:.6f}")
        print(f"H best-line RMS / stroke = {b['H_best_line_rms_over_stroke']:.6f}")
        print(f"H best-line max / stroke = {b['H_best_line_max_over_stroke']:.6f}")
        print(f"H best-line full span / stroke = {b['H_best_line_full_span_over_stroke']:.6f}")
        print(f"H best-line / cylinder angle = {b['H_best_line_axis_angle_deg']:.6f} deg")
        print(f"crank-axis clearance / crank = {b['crank_axis_to_EFH_clearance_over_crank']:.6f}")
        print(f"EF / crank = {b['EF_over_crank']:.6f}")
        print(f"GF / crank = {b['GF_over_crank']:.6f}")
        print(f"EH / crank = {b['EH_over_crank']:.6f}")
        print(f"piston rod / stroke = {b['piston_rod_over_stroke']:.6f}")
        print(f"primary transmission sine = {b['minimum_primary_transmission_sine']:.6f}")
        print(f"secondary transmission sine = {b['minimum_secondary_transmission_sine']:.6f}")
        print(f"rod-axis cosine = {b['minimum_rod_axis_cosine']:.6f}")
        print("primary =",json.dumps(b["primary"],indent=2))
    else:
        print("No dense feasible candidate found.")

    print(f"\nWrote {args.output}")

if __name__=="__main__":
    main()
