"""Two-pass whole-machine optimization at one source temperature.

Pass 1 freezes motion and varies:
  swept_ratio, n_i, length_i_m, n_o, length_o_m
Pass 2 freezes the best hardware and varies:
  t1, t2, t3, a_l, b_l, a_s, b_s

Every candidate is filled at 100 kPa absolute at theta=0. Total swept volume
and cylinder clearance ratios stay fixed. Histories are append-only and the
campaign resumes from them.

Example:
  PYTHONPATH=src python3 examples/optimize_motor_temperature_point.py \
    --delta-t-k 270 --power-floor-w 40
"""

from __future__ import annotations
from dataclasses import replace
import argparse, hashlib, json, math, time
from pathlib import Path

import numpy as np
from scipy.stats import qmc

from dada_solver.configuration import ChargeConfiguration
from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
from dada_solver.four_stage_kinematics import FourStageVolumeKinematics
from dada_solver.integration import IntegrationInterrupted

from compare_motor_motion_laws_stage7A5 import ROOT, _evaluate, _uniform_wall_initial
from optimize_motor_four_stage_thermo3d import _volume_limits
from optimize_motor_four_stage_k2 import _symmetry_diagnostics
from refine_motor_four_stage_hx9d_variable_gas import _hardware_metrics, _load_basis

REFERENCE_REPORT = ROOT/'outputs'/'motor_four_stage_hx9d_variable_gas'/'report.json'
OUTPUT_ROOT = ROOT/'outputs'/'motor_temperature_map'
COLD_K = 298.15
CHARGE_PA = 100000.0
MIN_STAGE = .02
MAX_P = 1_200_000.0
MAX_T = 850.0
MAX_FLOW = .08
RATIO_BOUNDS = (.50, 1.30)
COUNT_BOUNDS = (2400, 5200)
# Retain the broad physical length envelopes already explored by the historical
# thermo-3D campaign. If a temperature champion reaches a bound, expand it in a
# separately documented follow-up rather than silently changing the map domain.
HI_LENGTH_BOUNDS = (.03175, .06350)
HO_LENGTH_BOUNDS = (.016933333333333334, .035983333333333335)
KIN = ('t1','t2','t3','a_l','b_l','a_s','b_s')
THERMO = ('swept_ratio','n_i','length_i_m','n_o','length_o_m')
ALL = (*THERMO,*KIN)
DEFAULT_THERMO_RADII=((.08,240,.10),(.04,120,.05),(.02,60,.025),(.01,30,.0125))
DEFAULT_MOTION_RADII=((.01,.005),(.005,.0025),(.0025,.00125),(.001,.0005))


def load_history(path):
    return [] if not path.exists() else [json.loads(x) for x in path.read_text().splitlines() if x.strip()]

def append_record(path, record):
    with path.open('a') as f:
        f.write(json.dumps(record)+'\n'); f.flush()

def eligible(r):
    return bool(r.get('feasible') and r.get('result',{}).get('status')=='converged'
                and r.get('result',{}).get('indicated_thermal_efficiency') is not None)

def eta(r): return float(r['result']['indicated_thermal_efficiency'])

def phase_degrees(p):
    t1,t2,t3=(float(p[x]) for x in ('t1','t2','t3'))
    return dict(low_pressure_exchange_deg=360*t1,
                compression_deg=360*(t2-t1),
                high_pressure_exchange_deg=360*(t3-t2),
                expansion_deg=360*(1-t3))

def candidate_id(p):
    q={k:(int(p[k]) if k in ('n_i','n_o') else float(p[k])) for k in ALL}
    return hashlib.sha256(json.dumps(q,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def parse_thermo_radii(text):
    rows=[]
    for item in text.split(','):
        if item.strip():
            a,b,c=item.split(':'); rows.append((float(a),int(b),float(c)))
    if not rows or any(min(x)<=0 for x in rows): raise argparse.ArgumentTypeError('Invalid thermo radii.')
    return tuple(rows)

def parse_motion_radii(text):
    rows=[]
    for item in text.split(','):
        if item.strip():
            a,b=item.split(':'); rows.append((float(a),float(b)))
    if not rows or any(min(x)<=0 for x in rows): raise argparse.ArgumentTypeError('Invalid motion radii.')
    return tuple(rows)

def valid_motion(p):
    v=np.asarray([p[n] for n in KIN],float)
    return bool(np.all(np.isfinite(v)) and np.all((v[3:]>=0)&(v[3:]<=1))
                and np.min(np.diff(np.r_[0,v[:3],1]))>=MIN_STAGE)

def valid_thermo(p):
    return (RATIO_BOUNDS[0]<=float(p['swept_ratio'])<=RATIO_BOUNDS[1]
            and COUNT_BOUNDS[0]<=int(p['n_i'])<=COUNT_BOUNDS[1]
            and COUNT_BOUNDS[0]<=int(p['n_o'])<=COUNT_BOUNDS[1]
            and HI_LENGTH_BOUNDS[0]<=float(p['length_i_m'])<=HI_LENGTH_BOUNDS[1]
            and HO_LENGTH_BOUNDS[0]<=float(p['length_o_m'])<=HO_LENGTH_BOUNDS[1])

def base_geometry(base):
    s=base.configuration.machine_volumes.small_cylinder
    l=base.configuration.machine_volumes.large_cylinder
    return dict(total_swept=float(s.swept+l.swept),
                small_clear=float(s.minimum/s.swept),
                large_clear=float(l.minimum/l.swept),
                base_ratio=float(s.swept/l.swept),
                base_ni=int(base.heat_in.bank.tube_count),
                base_no=int(base.heat_out.bank.tube_count),
                base_li=float(base.heat_in.bank.tube_length_m),
                base_lo=float(base.heat_out.bank.tube_length_m),
                hi_cda=float(base.heat_in.outlet_valve_cda_m2),
                ho_cda=float(base.heat_out.outlet_valve_cda_m2),
                hot_ref_T=float(base.heat_in.inputs.air_inlet_temperature_k),
                hot_ref_rho=float(base.heat_in.inputs.air_density_kg_m3))

def load_source(path, base, g):
    data=json.loads(path.read_text())
    if 'champion' in data:
        c=data['champion']
        if c is None: raise RuntimeError(f'No champion in {path}')
        p={k:(int(c['parameters'][k]) if k in ('n_i','n_o') else float(c['parameters'][k])) for k in ALL}
        state=np.asarray(c.get('last_complete_state'),float)
        if state.shape!=(10,): raise RuntimeError('Temperature champion has no ten-state periodic state.')
        return p,state,c
    c=data.get('best_feasible')
    if c is None: raise RuntimeError('Reference 9D report has no best_feasible.')
    p=dict(swept_ratio=g['base_ratio'], n_i=int(c['parameters']['n_i']),
           length_i_m=float(c['hardware']['H_i']['tube_length_m']),
           n_o=int(c['parameters']['n_o']),
           length_o_m=float(c['hardware']['H_o']['tube_length_m']),
           **{n:float(c['parameters'][n]) for n in KIN})
    state=np.asarray(c.get('last_complete_state'),float)
    if state.shape!=(10,): raise RuntimeError('Reference champion has no ten-state periodic state.')
    return p,state,c

def build_design(base,g,p,hot_k):
    s,l=_volume_limits(g['total_swept'],float(p['swept_ratio']),g['small_clear'],g['large_clear'])
    cfg=replace(base.configuration,
        machine_volumes=replace(base.configuration.machine_volumes,small_cylinder=s,large_cylinder=l),
        charge=ChargeConfiguration(temperature=base.configuration.charge.temperature,pressure=CHARGE_PA))
    ni,no=int(p['n_i']),int(p['n_o'])
    hi=replace(base.heat_in,
        bank=replace(base.heat_in.bank,tube_count=ni,tube_length_m=float(p['length_i_m'])),
        outlet_valve_cda_m2=g['hi_cda']*ni/g['base_ni'])
    rho=g['hot_ref_rho']*g['hot_ref_T']/hot_k
    hi=replace(hi,inputs=replace(hi.inputs,air_inlet_temperature_k=hot_k,air_density_kg_m3=rho))
    ho=replace(base.heat_out,
        bank=replace(base.heat_out.bank,tube_count=no,tube_length_m=float(p['length_o_m'])),
        outlet_valve_cda_m2=g['ho_cda']*no/g['base_no'])
    ho=replace(ho,inputs=replace(ho.inputs,air_inlet_temperature_k=COLD_K))
    kin=FourStageVolumeKinematics(s,l,**{n:float(p[n]) for n in KIN})
    return replace(base,configuration=cfg,heat_in=hi,heat_out=ho,kinematics=kin)

def uniform_metadata(design):
    wrapper=design.build()
    state=np.asarray(_uniform_wall_initial(design.configuration,wrapper),float)
    s=design.configuration.machine_volumes.small_cylinder
    l=design.configuration.machine_volumes.large_cylinder
    hw={'H_i':_hardware_metrics(design.heat_in),'H_o':_hardware_metrics(design.heat_out)}
    geo=dict(small_swept_volume_m3=float(s.swept),large_swept_volume_m3=float(l.swept),
             total_swept_volume_m3=float(s.swept+l.swept),
             small_clearance_ratio=float(s.minimum/s.swept),large_clearance_ratio=float(l.minimum/l.swept),
             swept_ratio=float(s.swept/l.swept),gas_inventory_kg=float(np.sum(state[:8:2])))
    return state,hw,geo

def warm_start(source_state,source_hw,target_uniform,target_hw):
    x=np.asarray(source_state,float).copy()
    if x.shape!=(10,): raise ValueError('Expected ten-state warm start.')
    sm=float(np.sum(x[:8:2])); tm=float(np.sum(target_uniform[:8:2]))
    if sm<=0 or tm<=0: raise ValueError('Invalid warm-start inventory.')
    x[:8]*=tm/sm
    for i,side in ((8,'H_i'),(9,'H_o')):
        x[i]*=float(target_hw[side]['wall_capacity_j_k'])/float(source_hw[side]['wall_capacity_j_k'])
    return x

def feasibility(result,power_floor):
    reasons=[]
    e=result.get('indicated_thermal_efficiency'); p=result.get('indicated_power_w')
    if e is None or not math.isfinite(float(e)): reasons.append('motor_efficiency_unavailable')
    if p is None or float(p)<power_floor: reasons.append('minimum_motor_power')
    d=result.get('diagnostics',{})
    ps=d.get('pressure_extrema',{}); ts=d.get('temperature_extrema',{})
    pmax=max((float(v['maximum']) for v in ps.values()),default=None)
    tmax=max((float(v['maximum']) for v in ts.values()),default=None)
    if pmax is None or pmax>MAX_P: reasons.append('maximum_pressure')
    if tmax is None or tmax>MAX_T: reasons.append('maximum_temperature')
    flow=result.get('maximum_absolute_mass_flow_kg_s')
    if flow is None or float(flow)>MAX_FLOW: reasons.append('maximum_absolute_mass_flow')
    if result.get('validity',{}).get('verdict')!='valid': reasons.append('valid_thermodynamic_or_exchanger_model')
    result['maximum_pressure_pa']=pmax; result['maximum_temperature_k']=tmax
    return not reasons,reasons

def thermo_distance(p,r,g):
    q=r['parameters']
    a=np.asarray([p['swept_ratio'],p['n_i']/g['base_ni'],p['length_i_m']/g['base_li'],
                  p['n_o']/g['base_no'],p['length_o_m']/g['base_lo']],float)
    b=np.asarray([q['swept_ratio'],q['n_i']/g['base_ni'],q['length_i_m']/g['base_li'],
                  q['n_o']/g['base_no'],q['length_o_m']/g['base_lo']],float)
    return float(np.linalg.norm(a-b))

def motion_distance(p,r):
    return float(np.linalg.norm(np.asarray([p[n]-r['parameters'][n] for n in KIN],float)))

def eval_record(label,p,design,definition,initial,hw,geo,phase,ridx,radii,power_floor,hot_k,
                phase_deadline,candidate_seconds,warm_source):
    start=time.monotonic()
    candidate_deadline=min(phase_deadline,start+candidate_seconds)
    def progress(_):
        if time.monotonic()>=phase_deadline: raise IntegrationInterrupted(f'{phase} phase wall-clock budget exhausted.')
        if time.monotonic()>=candidate_deadline: raise IntegrationInterrupted('Candidate wall-clock budget exhausted.')
    try:
        result,state=_evaluate(label,design,definition,initial_state=np.asarray(initial),progress_callback=progress)
    except MicrotubeDomainError as exc:
        result={'status':'invalid_exchanger','message':str(exc)}; state=None
    except IntegrationInterrupted as exc:
        result={'status':'interrupted','message':str(exc)}; state=None
    except (ValueError,RuntimeError) as exc:
        result={'status':'integration_failure','message':str(exc)}; state=None
    ok,reasons=feasibility(result,power_floor) if result.get('status')=='converged' else (False,[])
    ec=1-COLD_K/hot_k; ev=result.get('indicated_thermal_efficiency')
    return dict(candidate_id=candidate_id(p),phase=phase,radius_index=ridx,**radii,
        parameters={n:(int(p[n]) if n in ('n_i','n_o') else float(p[n])) for n in ALL},
        phase_degrees=phase_degrees(p),symmetry=_symmetry_diagnostics(p),hardware=hw,derived_geometry=geo,
        result=result,carnot_efficiency=ec,fraction_of_carnot=(float(ev)/ec if ev is not None else None),
        feasible=ok,physical_constraint_failures=reasons,elapsed_seconds=time.monotonic()-start,
        warm_start_source=warm_source,last_complete_state=(state.tolist() if state is not None else None))

def write_report(path,delta,hot_k,power_floor,source_report,source_p,th,mot,th_target,mot_target,definition_path):
    th_ok=[r for r in th if eligible(r)]; mo_ok=[r for r in mot if eligible(r)]
    bt=max(th_ok,key=eta) if th_ok else None; bm=max(mo_ok,key=eta) if mo_ok else None
    nt=sum(
        r.get('kind')!='source_seed'
        and r.get('result',{}).get('status')!='interrupted'
        for r in th
    )
    nm=sum(
        r.get('kind')!='thermo_seed'
        and r.get('result',{}).get('status')!='interrupted'
        for r in mot
    )
    tc=nt>=th_target; mc=nm>=mot_target
    report=dict(study='two-pass whole-machine temperature continuation',target_delta_t_k=delta,
        cold_source_temperature_k=COLD_K,hot_source_temperature_k=hot_k,power_floor_w=power_floor,
        source_report=str(source_report),source_parameters=source_p,definition=str(definition_path),
        thermo_pass=dict(complete=tc,target_nonseed_evaluations=th_target,attempted_total=len(th),
                         completed_nonseed_evaluations=nt,best_feasible=bt),
        motion_pass=dict(complete=mc,target_nonseed_evaluations=mot_target,attempted_total=len(mot),
                         completed_nonseed_evaluations=nm,best_feasible=bm),
        campaign_complete=bool(tc and mc and bm is not None),champion=(bm if bm is not None else bt))
    path.write_text(json.dumps(report,indent=2)+'\n')
    return report

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--delta-t-k',type=float,required=True)
    ap.add_argument('--power-floor-w',type=float,required=True)
    ap.add_argument('--source-report',type=Path,default=REFERENCE_REPORT)
    ap.add_argument('--output-directory',type=Path)
    ap.add_argument('--thermo-evaluations',type=int,default=96)
    ap.add_argument('--motion-evaluations',type=int,default=96)
    ap.add_argument('--thermo-evaluations-per-radius',type=int,default=24)
    ap.add_argument('--motion-evaluations-per-radius',type=int,default=24)
    ap.add_argument('--thermo-budget-seconds',type=float,default=10800)
    ap.add_argument('--motion-budget-seconds',type=float,default=10800)
    ap.add_argument('--candidate-seconds',type=float,default=240)
    ap.add_argument('--thermo-radii',type=parse_thermo_radii,default=DEFAULT_THERMO_RADII)
    ap.add_argument('--motion-radii',type=parse_motion_radii,default=DEFAULT_MOTION_RADII)
    ap.add_argument('--seed',type=int,default=26091730)
    args=ap.parse_args()
    if args.delta_t_k<=0 or args.power_floor_w<0: raise ValueError('Invalid delta T or power floor.')
    if min(args.thermo_evaluations,args.motion_evaluations,args.thermo_evaluations_per_radius,args.motion_evaluations_per_radius)<=0: raise ValueError('Evaluation counts must be positive.')
    if min(args.thermo_budget_seconds,args.motion_budget_seconds,args.candidate_seconds)<=0: raise ValueError('Budgets must be positive.')

    definition,base=_load_basis(); g=base_geometry(base); hot_k=COLD_K+args.delta_t_k
    source=args.source_report if args.source_report.is_absolute() else ROOT/args.source_report
    if not source.exists(): raise FileNotFoundError(source)
    source_p,source_state,source_record=load_source(source,base,g)
    tag=f"dT_{args.delta_t_k:g}K".replace('.','p')
    out=args.output_directory or OUTPUT_ROOT/tag
    out=out if out.is_absolute() else ROOT/out
    out.mkdir(parents=True,exist_ok=True)
    dp=out/'definition.json'; thp=out/'thermo_history.jsonl'; mop=out/'motion_history.jsonl'; rp=out/'report.json'
    identity=dict(target_delta_t_k=args.delta_t_k,cold_source_temperature_k=COLD_K,hot_source_temperature_k=hot_k,
        power_floor_w=args.power_floor_w,charge_pressure_pa=CHARGE_PA,charge_temperature_k=base.configuration.charge.temperature,
        source_report=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),source_parameters=source_p,
        total_swept_volume_m3=g['total_swept'],clearance_ratios=[g['small_clear'],g['large_clear']],
        thermo_radii=[list(x) for x in args.thermo_radii],motion_radii=[list(x) for x in args.motion_radii],
        thermo_evaluations_per_radius=args.thermo_evaluations_per_radius,motion_evaluations_per_radius=args.motion_evaluations_per_radius,
        bounds=dict(swept_ratio=list(RATIO_BOUNDS),tube_count=list(COUNT_BOUNDS),
                    H_i_tube_length_m=list(HI_LENGTH_BOUNDS),H_o_tube_length_m=list(HO_LENGTH_BOUNDS)),
        minimum_stage_fraction=MIN_STAGE,seed=args.seed,
        strategy='Local thermo5D continuation followed by local motion7D continuation; 100 kPa theta=0 filling.')
    if dp.exists() and json.loads(dp.read_text())!=identity: raise ValueError('Definition changed; use a new output directory.')
    if not dp.exists(): dp.write_text(json.dumps(identity,indent=2)+'\n')

    # Thermo pass.
    th=load_history(thp); bt=max((r for r in th if eligible(r)),key=eta,default=None)
    finished={r['candidate_id'] for r in th if r.get('result',{}).get('status')!='interrupted'}
    n=sum(
        r.get('kind')!='source_seed'
        and r.get('result',{}).get('status')!='interrupted'
        for r in th
    )
    pts=qmc.Sobol(d=5,scramble=True,seed=args.seed+int(round(args.delta_t_k*10))).random_base2(14)
    proposal=max((r.get('proposal_index',-1) for r in th),default=-1)+1
    deadline=time.monotonic()+args.thermo_budget_seconds
    while n<args.thermo_evaluations and time.monotonic()<deadline:
        pi=proposal; proposal+=1
        if not th and pi==0:
            p=dict(source_p); kind='source_seed'; ridx=-1
            radii=dict(swept_ratio_radius=None,tube_count_radius=None,tube_length_relative_radius=None)
        else:
            ridx=min(n//args.thermo_evaluations_per_radius,len(args.thermo_radii)-1)
            rr,nr,lr=args.thermo_radii[ridx]; c=bt['parameters'] if bt is not None else source_p; u=pts[pi-1]
            p=dict(source_p)
            p.update(swept_ratio=float(c['swept_ratio']+(2*u[0]-1)*rr),
                     n_i=int(round(c['n_i']+(2*u[1]-1)*nr)),
                     length_i_m=float(c['length_i_m']*(1+(2*u[2]-1)*lr)),
                     n_o=int(round(c['n_o']+(2*u[3]-1)*nr)),
                     length_o_m=float(c['length_o_m']*(1+(2*u[4]-1)*lr)))
            for name in KIN: p[name]=float(source_p[name])
            kind='adaptive_local_thermo5d'
            radii=dict(swept_ratio_radius=rr,tube_count_radius=nr,tube_length_relative_radius=lr)
            if not valid_thermo(p): continue
        cid=candidate_id(p)
        if cid in finished: continue
        design=build_design(base,g,p,hot_k); uniform,hw,geo=uniform_metadata(design)
        reusable=[r for r in th if r.get('result',{}).get('status')=='converged' and r.get('last_complete_state') is not None]
        ws=None
        if reusable:
            ws=min(reusable,key=lambda r:thermo_distance(p,r,g))
            initial=warm_start(ws['last_complete_state'],ws['hardware'],uniform,hw); wsid=ws['candidate_id']
        else:
            initial=warm_start(source_state,source_record['hardware'],uniform,hw); wsid='source_report'
        rec=eval_record(f'temperature_{args.delta_t_k:g}_thermo_{pi}',p,design,definition,initial,hw,geo,'thermo',ridx,radii,
                        args.power_floor_w,hot_k,deadline,args.candidate_seconds,wsid)
        rec.update(index=len(th),proposal_index=pi,kind=kind)
        append_record(thp,rec); th.append(rec)
        if rec.get('result',{}).get('status')!='interrupted':
            finished.add(cid)
            if kind!='source_seed':
                n+=1
        if eligible(rec) and (bt is None or eta(rec)>eta(bt)):
            bt=rec
            print(f"[thermo] eta={100*eta(rec):.6f}% P={rec['result']['indicated_power_w']:.3f} W "
                  f"r={p['swept_ratio']:.4f} Hi={p['n_i']}x{1e3*p['length_i_m']:.2f}mm "
                  f"Ho={p['n_o']}x{1e3*p['length_o_m']:.2f}mm",flush=True)
        write_report(rp,args.delta_t_k,hot_k,args.power_floor_w,source,source_p,th,load_history(mop),
                     args.thermo_evaluations,args.motion_evaluations,dp)

    if n<args.thermo_evaluations:
        write_report(rp,args.delta_t_k,hot_k,args.power_floor_w,source,source_p,th,load_history(mop),
                     args.thermo_evaluations,args.motion_evaluations,dp)
        print(f"Thermo pass incomplete: {n}/{args.thermo_evaluations}. Rerun to resume.",flush=True); return
    if bt is None: raise RuntimeError('Thermo pass has no feasible candidate; review power floor/constraints.')

    # Motion pass.
    frozen={name:bt['parameters'][name] for name in THERMO}; seedp=dict(bt['parameters'])
    mo=load_history(mop); bm=max((r for r in mo if eligible(r)),key=eta,default=None)
    finished={r['candidate_id'] for r in mo if r.get('result',{}).get('status')!='interrupted'}
    n=sum(
        r.get('kind')!='thermo_seed'
        and r.get('result',{}).get('status')!='interrupted'
        for r in mo
    )
    pts=qmc.Sobol(d=7,scramble=True,seed=args.seed+100000+int(round(args.delta_t_k*10))).random_base2(14)
    proposal=max((r.get('proposal_index',-1) for r in mo),default=-1)+1
    deadline=time.monotonic()+args.motion_budget_seconds
    while n<args.motion_evaluations and time.monotonic()<deadline:
        pi=proposal; proposal+=1
        if not mo and pi==0:
            p=dict(seedp); kind='thermo_seed'; ridx=-1; radii=dict(timing_radius=None,level_radius=None)
        else:
            ridx=min(n//args.motion_evaluations_per_radius,len(args.motion_radii)-1)
            tr,lr=args.motion_radii[ridx]; c=bm['parameters'] if bm is not None else seedp; u=pts[pi-1]
            p=dict(seedp)
            for i,name in enumerate(('t1','t2','t3')): p[name]=float(c[name]+(2*u[i]-1)*tr)
            for i,name in enumerate(('a_l','b_l','a_s','b_s'),start=3): p[name]=float(c[name]+(2*u[i]-1)*lr)
            for name in THERMO: p[name]=frozen[name]
            kind='adaptive_local_motion7d'; radii=dict(timing_radius=tr,level_radius=lr)
            if not valid_motion(p): continue
        cid=candidate_id(p)
        if cid in finished: continue
        design=build_design(base,g,p,hot_k); uniform,hw,geo=uniform_metadata(design)
        reusable=[r for r in mo if r.get('result',{}).get('status')=='converged' and r.get('last_complete_state') is not None]
        if reusable:
            ws=min(reusable,key=lambda r:motion_distance(p,r))
            initial=warm_start(ws['last_complete_state'],ws['hardware'],uniform,hw); wsid=ws['candidate_id']
        else:
            initial=warm_start(bt['last_complete_state'],bt['hardware'],uniform,hw); wsid='thermo_best'
        rec=eval_record(f'temperature_{args.delta_t_k:g}_motion_{pi}',p,design,definition,initial,hw,geo,'motion',ridx,radii,
                        args.power_floor_w,hot_k,deadline,args.candidate_seconds,wsid)
        rec.update(index=len(mo),proposal_index=pi,kind=kind)
        append_record(mop,rec); mo.append(rec)
        if rec.get('result',{}).get('status')!='interrupted':
            finished.add(cid)
            if kind!='thermo_seed':
                n+=1
        if eligible(rec) and (bm is None or eta(rec)>eta(bm)):
            bm=rec; ph=rec['phase_degrees']
            print(f"[motion] eta={100*eta(rec):.6f}% P={rec['result']['indicated_power_w']:.3f} W "
                  f"phases={ph['low_pressure_exchange_deg']:.1f}/{ph['compression_deg']:.1f}/"
                  f"{ph['high_pressure_exchange_deg']:.1f}/{ph['expansion_deg']:.1f} deg",flush=True)
        write_report(rp,args.delta_t_k,hot_k,args.power_floor_w,source,source_p,th,mo,
                     args.thermo_evaluations,args.motion_evaluations,dp)

    report=write_report(rp,args.delta_t_k,hot_k,args.power_floor_w,source,source_p,th,mo,
                        args.thermo_evaluations,args.motion_evaluations,dp)
    if report['campaign_complete']:
        c=report['champion']; r=c['result']
        print(f"COMPLETE deltaT={args.delta_t_k:g} K eta={100*r['indicated_thermal_efficiency']:.6f}% "
              f"eta/Carnot={100*c['fraction_of_carnot']:.3f}% P={r['indicated_power_w']:.3f} W",flush=True)
    else:
        print(f"Motion pass incomplete: {n}/{args.motion_evaluations}. Rerun to resume.",flush=True)

if __name__=='__main__':
    main()
