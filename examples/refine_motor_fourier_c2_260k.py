#!/usr/bin/env python3
"""Refine the 260 K smooth UU champion, optionally with more Fourier harmonics.

Examples:
  # finish 4 harmonics
  PYTHONPATH=src python3 examples/refine_motor_fourier_c2_260k.py \
    --harmonics 4 \
    --source-report outputs/motor_fourier_c2_260k/report.json \
    --output-directory outputs/motor_fourier_c2_4h_refine \
    --evaluations 512 --budget-seconds 7200

  # expand the refined champion to 6 harmonics
  PYTHONPATH=src python3 examples/refine_motor_fourier_c2_260k.py \
    --harmonics 6 \
    --source-report outputs/motor_fourier_c2_4h_refine/report.json \
    --output-directory outputs/motor_fourier_c2_6h_refine \
    --evaluations 1024 --budget-seconds 10800

The seed is preserved exactly when the Fourier order is increased: new
coefficients start at zero. Thermodynamic/hardware variables remain hard-bound
to +/-3% around the original linear 260 K UU champion.
"""
from __future__ import annotations
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
import argparse, hashlib, json, math, time
import numpy as np
from scipy.stats import qmc

from dada_solver.integration import IntegrationInterrupted
from dada_solver.wall_backend import WallBackendSettings
from compare_motor_motion_laws_stage7A5 import ROOT, _evaluate
from optimize_motor_fourier_c2_260k import FourierVolumeKinematics
from optimize_motor_temperature_point import (
    COLD_K, base_geometry, build_design as linear_build_design,
    feasibility, uniform_metadata, warm_start,
)
from refine_motor_four_stage_hx9d_variable_gas import _load_basis

DEFAULT_SOURCE = ROOT/'outputs'/'motor_fourier_c2_260k'/'report.json'
DEFAULT_OUTPUT = ROOT/'outputs'/'motor_fourier_c2_4h_refine'
HOT_K = COLD_K + 260.0
POWER_FLOOR_W = 25.0
THERMO_LIMIT = .03
THERMO = ('swept_ratio','n_i','length_i_m','n_o','length_o_m')
RADII = (.018,.012,.008,.005,.003,.0015)
THERMO_RADII = (.004,.003,.002,.0012,.0007,.0004)
EVALS_PER_RADIUS = 96
SEED = 26092141

def load_history(path):
    return [] if not path.exists() else [
        json.loads(x) for x in path.read_text().splitlines() if x.strip()
    ]

def eligible(r):
    z=r.get('result',{})
    return bool(r.get('feasible') and z.get('status')=='converged'
                and z.get('indicated_thermal_efficiency') is not None)

def eta(r): return float(r['result']['indicated_thermal_efficiency'])

def source_best(report):
    b=report.get('best_feasible')
    if b is None or b.get('last_complete_state') is None:
        raise RuntimeError('Source report has no usable best_feasible.')
    return b

def source_harmonics(record):
    n=len(record['small_coefficients'])
    if n%2 or len(record['large_coefficients'])!=n:
        raise RuntimeError('Invalid source Fourier coefficient arrays.')
    return n//2

def pad(values, h0, h):
    if h<h0: raise ValueError('Target harmonic count cannot be lower than source.')
    out=np.zeros(2*h,float); out[:2*h0]=np.asarray(values,float); return out

def normalize(v):
    v=np.asarray(v,float); n=float(np.linalg.norm(v))
    if not math.isfinite(n) or n<1e-12: raise ValueError('Degenerate coefficient vector.')
    return v/n

def raw_curve(coeffs,h,n=4096):
    t=np.linspace(0,1,n,endpoint=False); q=np.zeros(n); d=np.zeros(n)
    for k in range(1,h+1):
        a,b=coeffs[2*k-2:2*k]; w=2*math.pi*k
        q += a*np.cos(w*t)+b*np.sin(w*t)
        d += -a*w*np.sin(w*t)+b*w*np.cos(w*t)
    span=float(q.max()-q.min())
    if span<=1e-9: raise ValueError('Degenerate Fourier motion.')
    return t,(q-q.min())/span,d/span

def extrema_count(coeffs,h):
    _,_,d=raw_curve(coeffs,h)
    z=max(1e-8,1e-4*float(np.max(np.abs(d))))
    s=np.sign(d); s[np.abs(d)<z]=0; s=s[s!=0]
    if not len(s): return 0
    return int(np.sum(s[1:]!=s[:-1])+(s[0]!=s[-1]))

def admissible(s,l,h):
    return extrema_count(s,h)==2 and extrema_count(l,h)==2

def bounds(anchor):
    b={}
    for n in ('swept_ratio','length_i_m','length_o_m'):
        x=float(anchor[n]); b[n]=(x*(1-THERMO_LIMIT),x*(1+THERMO_LIMIT))
    for n in ('n_i','n_o'):
        x=int(anchor[n]); b[n]=(math.floor(x*(1-THERMO_LIMIT)),
                               math.ceil(x*(1+THERMO_LIMIT)))
    return b

def propose_thermo(center,anchor,b,u,ri):
    rel=THERMO_RADII[min(ri,len(THERMO_RADII)-1)]; p=dict(center)
    for i,n in enumerate(('swept_ratio','length_i_m','length_o_m')):
        x=float(center[n])*(1+(2*u[i]-1)*rel); p[n]=float(np.clip(x,*b[n]))
    for j,n in enumerate(('n_i','n_o'),3):
        rad=max(2,int(round(float(anchor[n])*rel)))
        x=int(round(float(center[n])+(2*u[j]-1)*rad))
        p[n]=int(np.clip(x,*b[n]))
    return p

def build_design(base,g,thermo,s,l,h):
    placeholder=dict(thermo,t1=.25,t2=.5,t3=.75,a_l=.5,b_l=.5,a_s=.5,b_s=.5)
    d=linear_build_design(base,g,placeholder,HOT_K)
    lim=d.configuration.machine_volumes
    kin=FourierVolumeKinematics(
        lim.small_cylinder,lim.large_cylinder,
        tuple(map(float,s)),tuple(map(float,l)),harmonics=h)
    return replace(d,
        configuration=replace(d.configuration,
            heat_in_valve_placement='upstream',heat_out_valve_placement='upstream'),
        kinematics=kin)

def candidate_id(t,s,l,h):
    x={'harmonics':h,'thermo':{
        n:(int(t[n]) if n in ('n_i','n_o') else float(t[n])) for n in THERMO},
       'small':list(map(float,s)),'large':list(map(float,l))}
    return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def distance(t,s,l,r,anchor):
    x=0.
    for n in THERMO:
        scale=max(abs(float(anchor[n])),1e-12)
        x += ((float(t[n])-float(r['thermo'][n]))/scale)**2
    ds=np.asarray(s)-np.asarray(r['small_coefficients'])
    dl=np.asarray(l)-np.asarray(r['large_coefficients'])
    return math.sqrt(x+float(ds@ds)+float(dl@dl))

def shape_diag(s,l,h):
    ts,qs,ds=raw_curve(s,h,8192); tl,ql,dl=raw_curve(l,h,8192)
    return {
        'small_max_deg':360*float(ts[np.argmax(qs)]),
        'small_min_deg':360*float(ts[np.argmin(qs)]),
        'large_max_deg':360*float(tl[np.argmax(ql)]),
        'large_min_deg':360*float(tl[np.argmin(ql)]),
        'small_peak_abs_dq_dt':float(np.max(np.abs(ds))),
        'large_peak_abs_dq_dt':float(np.max(np.abs(dl))),
        'small_extrema_count':extrema_count(s,h),
        'large_extrema_count':extrema_count(l,h)}

def save_motion(path,s,l,h):
    t,qs,ds=raw_curve(s,h,2880); _,ql,dl=raw_curve(l,h,2880)
    a=np.column_stack((np.r_[t,1.],360*np.r_[t,1.],
        np.r_[qs,qs[0]],np.r_[ql,ql[0]],np.r_[ds,ds[0]],np.r_[dl,dl[0]]))
    np.savetxt(path,a,delimiter=',',
        header='motor_time_fraction,motor_angle_deg,small_fraction_0_1,large_fraction_0_1,dsmall_dt,dlarge_dt',
        comments='')

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-report',type=Path,default=DEFAULT_SOURCE)
    ap.add_argument('--output-directory',type=Path,default=DEFAULT_OUTPUT)
    ap.add_argument('--harmonics',type=int,default=4)
    ap.add_argument('--budget-seconds',type=float,default=7200)
    ap.add_argument('--evaluations',type=int,default=512)
    ap.add_argument('--candidate-seconds',type=float,default=180)
    ap.add_argument('--evaluations-per-radius',type=int,default=EVALS_PER_RADIUS)
    ap.add_argument('--seed',type=int,default=SEED)
    args=ap.parse_args()
    if args.harmonics<1 or min(args.budget_seconds,args.evaluations,
        args.candidate_seconds,args.evaluations_per_radius)<=0: raise ValueError('Invalid arguments.')

    sp=args.source_report if args.source_report.is_absolute() else ROOT/args.source_report
    report=json.loads(sp.read_text()); sb=source_best(report)
    ref=report.get('reference_linear_UU')
    if not ref or not ref.get('parameters'): raise RuntimeError('Missing reference_linear_UU.')
    anchor=ref['parameters']; h0=source_harmonics(sb); h=args.harmonics
    s0=pad(sb['small_coefficients'],h0,h); l0=pad(sb['large_coefficients'],h0,h)
    if not admissible(s0,l0,h): raise RuntimeError('Expanded source seed is not admissible.')
    t0={n:(int(sb['thermo'][n]) if n in ('n_i','n_o') else float(sb['thermo'][n])) for n in THERMO}

    cd,base=_load_basis()
    definition=SimpleNamespace(wall_numerical_settings=cd.wall_numerical_settings,
        wall_backend=WallBackendSettings('numba'),exact_kinematics_cache=True,shared_replay=True)
    g=base_geometry(base); bd=bounds(anchor)

    out=args.output_directory if args.output_directory.is_absolute() else ROOT/args.output_directory
    out.mkdir(parents=True,exist_ok=True); hp=out/'history.jsonl'; rp=out/'report.json'; dp=out/'definition.json'
    identity=dict(study='260 K UU local smooth-motion refinement',source_report=str(sp),
        source_sha256=hashlib.sha256(sp.read_bytes()).hexdigest(),source_harmonics=h0,target_harmonics=h,
        shape_radii=list(RADII),evaluations_per_radius=args.evaluations_per_radius,sobol_seed=args.seed,
        thermo_hard_relative_limit_vs_original_linear_UU=THERMO_LIMIT,
        thermo_bounds={k:list(v) for k,v in bd.items()},source_efficiency=sb['result']['indicated_thermal_efficiency'],
        reference_linear_UU_efficiency=ref['efficiency'],
        strategy='Local Sobol refinement; progressively smaller shape and thermo radii; every eighth proposal probes only newly added harmonics.')
    if dp.exists() and json.loads(dp.read_text())!=identity:
        raise RuntimeError('Definition changed; use a new output directory.')
    if not dp.exists(): dp.write_text(json.dumps(identity,indent=2)+'\n')

    history=load_history(hp); valid=[r for r in history if eligible(r)]
    best=max(valid,key=eta) if valid else None
    finished={r['candidate_id'] for r in history if r.get('result',{}).get('status')!='interrupted'}
    pts=qmc.Sobol(d=5+4*h,scramble=True,seed=args.seed).random_base2(16)
    start=time.monotonic(); deadline=start+args.budget_seconds; completed=0
    pi=max((r.get('proposal_index',-1) for r in history),default=-1)+1

    while completed<args.evaluations and time.monotonic()<deadline:
        this=pi; pi+=1
        if not history and this==0:
            thermo=dict(t0); s=s0.copy(); l=l0.copy(); kind='expanded_source_seed' if h>h0 else 'source_seed'
            ri=-1; radius=0.
        else:
            nn=sum(1 for r in history if r.get('kind') not in ('source_seed','expanded_source_seed'))
            ri=min(nn//args.evaluations_per_radius,len(RADII)-1); radius=RADII[ri]
            u=pts[this%len(pts)]
            center=best if best is not None else {'thermo':t0,'small_coefficients':s0,'large_coefficients':l0}
            thermo=propose_thermo(center['thermo'],anchor,bd,u[:5],ri)
            cs=np.asarray(center['small_coefficients']); cl=np.asarray(center['large_coefficients'])
            ds=2*u[5:5+2*h]-1; dl=2*u[5+2*h:]-1
            weights=np.repeat(1/np.sqrt(np.arange(1,h+1)),2)
            high_only=h>h0 and this%8==0
            if high_only:
                mask=np.zeros(2*h); mask[2*h0:]=1
                s=normalize(cs+1.5*radius*ds*weights*mask)
                l=normalize(cl+1.5*radius*dl*weights*mask); kind='new_harmonics_probe'
            else:
                s=normalize(cs+radius*ds*weights); l=normalize(cl+radius*dl*weights); kind='local_shape_refinement'
            if not admissible(s,l,h): continue

        cid=candidate_id(thermo,s,l,h)
        if cid in finished: continue
        before=time.monotonic(); design=build_design(base,g,thermo,s,l,h)
        try:
            uniform,hardware,geo=uniform_metadata(design)
            reusable=[r for r in history if r.get('result',{}).get('status')=='converged'
                      and r.get('last_complete_state') is not None and r.get('hardware') is not None]
            warm=None
            if reusable:
                warm=min(reusable,key=lambda r:distance(thermo,s,l,r,anchor))
                initial=warm_start(warm['last_complete_state'],warm['hardware'],uniform,hardware)
            else:
                initial=warm_start(sb['last_complete_state'],sb['hardware'],uniform,hardware)
            candidate_deadline=min(deadline,before+args.candidate_seconds)
            def progress(_):
                now=time.monotonic()
                if now>=deadline: raise IntegrationInterrupted('Campaign wall-clock budget exhausted.')
                if now>=candidate_deadline: raise IntegrationInterrupted('Candidate wall-clock budget exhausted.')
            try:
                result,state=_evaluate(f'fourier_c2_refine_{h}h_{this}',design,definition,
                    initial_state=np.asarray(initial),progress_callback=progress)
            except IntegrationInterrupted as exc: result={'status':'interrupted','message':str(exc)}; state=None
            except (ValueError,RuntimeError) as exc: result={'status':'integration_failure','message':str(exc)}; state=None
        except (ValueError,RuntimeError) as exc:
            result={'status':'setup_failure','message':str(exc)}; state=None; hardware=None; geo=None

        feasible,reasons=feasibility(result,POWER_FLOOR_W) if result.get('status')=='converged' else (False,[])
        rec=dict(index=len(history),proposal_index=this,candidate_id=cid,kind=kind,harmonics=h,
            radius_index=ri,shape_radius=radius,thermo={n:(int(thermo[n]) if n in ('n_i','n_o') else float(thermo[n])) for n in THERMO},
            small_coefficients=list(map(float,s)),large_coefficients=list(map(float,l)),
            shape_diagnostics=shape_diag(s,l,h),hardware=hardware,derived_geometry=geo,result=result,
            feasible=feasible,physical_constraint_failures=reasons,elapsed_seconds=time.monotonic()-before,
            warm_start_source=warm['candidate_id'] if warm is not None else 'source_best',
            last_complete_state=state.tolist() if state is not None else None)
        with hp.open('a') as f: f.write(json.dumps(rec)+'\n'); f.flush()
        history.append(rec); completed+=1
        if result.get('status')!='interrupted': finished.add(cid)
        moved=False
        if feasible and (best is None or eta(rec)>eta(best)):
            best=rec; moved=True; save_motion(out/'best_motion.csv',s,l,h)
        outreport=dict(best_feasible=best,source_best=dict(efficiency=sb['result']['indicated_thermal_efficiency'],
            power_W=sb['result']['indicated_power_w'],harmonics=h0,thermo=sb['thermo']),
            reference_linear_UU=ref,attempted_total=len(history),completed_this_run=completed,
            requested_duration_seconds=args.budget_seconds,actual_duration_seconds=time.monotonic()-start,
            definition=identity)
        rp.write_text(json.dumps(outreport,indent=2)+'\n')
        print(json.dumps(dict(index=rec['index'],kind=kind,harmonics=h,radius=radius,
            status=result.get('status'),feasible=feasible,efficiency=result.get('indicated_thermal_efficiency'),
            power_W=result.get('indicated_power_w'),center_moved=moved,best_efficiency=eta(best) if best else None,
            elapsed_seconds=rec['elapsed_seconds'])),flush=True)
    print(f'Saved {rp}',flush=True)

if __name__=='__main__': main()
