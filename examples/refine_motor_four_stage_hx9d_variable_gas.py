"""Local 9D refinement around the balanced four-stage point using the production microtube gas model.

Nine free variables: t1,t2,t3,a_l,b_l,a_s,b_s,n_i,n_o.
The source report supplies only the starting kinematics; its legacy state/result
is never reused. K2 tube lengths and 0.33 mm ID stay fixed. Tube count changes
surface, flow area, headers/dead volume and wall capacity. Outlet-valve CdA is
scaled linearly with tube count. Every candidate is filled at 100 kPa at theta=0.

Run:
  PYTHONPATH=src python3 examples/refine_motor_four_stage_hx9d_variable_gas.py \
    --budget-seconds 21600 --evaluations 288 --candidate-seconds 180
"""
from __future__ import annotations
from dataclasses import replace
import argparse, hashlib, json, math, time
from pathlib import Path
import numpy as np
from scipy.stats import qmc

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.configuration import ChargeConfiguration
from dada_solver.exchangers.gas_correlations import MicrotubeDomainError, MicrotubeGasModel
from dada_solver.four_stage_kinematics import FourStageVolumeKinematics
from dada_solver.integration import IntegrationInterrupted
from compare_motor_motion_laws_stage7A5 import (
    A5_CAMPAIGN, ROOT, _candidate_design_and_mass, _evaluate, _uniform_wall_initial,
)
from optimize_motor_exchanger_asymmetry_stageK1 import _scaled_exchanger
from optimize_motor_four_stage_k2 import _symmetry_diagnostics

SOURCE_REPORT = ROOT/'outputs'/'motor_four_stage_from_150_30'/'report.json'
THERMO3D_REPORT = ROOT/'outputs'/'motor_four_stage_thermo3d'/'report.json'
DEFAULT_DIRECTORY = ROOT/'outputs'/'motor_four_stage_hx9d_variable_gas'
KIN_NAMES=('t1','t2','t3','a_l','b_l','a_s','b_s')
PARAMETER_NAMES=(*KIN_NAMES,'n_i','n_o')
CHARGE_PRESSURE_PA=100000.0
MINIMUM_STAGE_FRACTION=.02
MINIMUM_POWER_W=40.0
MAXIMUM_PRESSURE_PA=1_200_000.0
MAXIMUM_TEMPERATURE_K=850.0
MAXIMUM_ABSOLUTE_MASS_FLOW_KG_S=.08
COUNT_BOUNDS=(2400,5200)
SOBOL_SEED=26091623
DEFAULT_RADII=((.020,.010,500),(.010,.005,300),(.005,.0025,150),(.0025,.00125,80),(.001,.0005,40),(.0005,.00025,20))

def _load_history(path):
    return [] if not path.exists() else [json.loads(x) for x in path.read_text().splitlines() if x.strip()]

def _eligible(r):
    return bool(r.get('feasible') and r.get('result',{}).get('status')=='converged' and r.get('result',{}).get('indicated_thermal_efficiency') is not None)

def _efficiency(r): return float(r['result']['indicated_thermal_efficiency'])

def _phase_degrees(p):
    t1,t2,t3=(float(p[x]) for x in ('t1','t2','t3'))
    return dict(low_pressure_exchange_deg=360*t1,compression_deg=360*(t2-t1),high_pressure_exchange_deg=360*(t3-t2),expansion_deg=360*(1-t3))

def _valid_parameters(p):
    v=np.asarray([p[n] for n in KIN_NAMES],float)
    return bool(np.all(np.isfinite(v)) and np.all((v[3:]>=0)&(v[3:]<=1)) and np.min(np.diff(np.r_[0,v[:3],1]))>=MINIMUM_STAGE_FRACTION and all(COUNT_BOUNDS[0]<=int(p[n])<=COUNT_BOUNDS[1] for n in ('n_i','n_o')))

def _candidate_id(p):
    q={**{n:float(p[n]) for n in KIN_NAMES},'n_i':int(p['n_i']),'n_o':int(p['n_o'])}
    return hashlib.sha256(json.dumps(q,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def _load_center():
    best=json.loads(SOURCE_REPORT.read_text()).get('best_feasible')
    if best is None: raise RuntimeError('Missing best_feasible in balanced 7D source report.')
    p=best['parameters']
    return {**{n:float(p[n]) for n in KIN_NAMES},'n_i':3708,'n_o':3708}

def _load_basis():
    control=json.loads(THERMO3D_REPORT.read_text()).get('control_1bar')
    if control is None or not control.get('feasible'): raise RuntimeError('Missing feasible thermo-3D K2 control.')
    definition=CampaignDefinition(A5_CAMPAIGN)
    base,*_= _candidate_design_and_mass(definition)
    hi=_scaled_exchanger(base.heat_in,float(control['parameters']['k_i']))
    ho=_scaled_exchanger(base.heat_out,float(control['parameters']['k_o']))
    production=MicrotubeGasModel()
    hi=replace(hi,inputs=replace(hi.inputs,gas_model=production))
    ho=replace(ho,inputs=replace(ho.inputs,gas_model=production))
    if hi.bank.tube_count!=3708 or ho.bank.tube_count!=3708: raise RuntimeError('Expected K2 basis with 3708 tubes.')
    charge=ChargeConfiguration(temperature=base.configuration.charge.temperature,pressure=CHARGE_PRESSURE_PA)
    return definition,replace(base,configuration=replace(base.configuration,charge=charge),heat_in=hi,heat_out=ho)

def _build_design(base,p):
    ni,no=int(p['n_i']),int(p['n_o'])
    hi=replace(base.heat_in,bank=replace(base.heat_in.bank,tube_count=ni),outlet_valve_cda_m2=base.heat_in.outlet_valve_cda_m2*ni/base.heat_in.bank.tube_count)
    ho=replace(base.heat_out,bank=replace(base.heat_out.bank,tube_count=no),outlet_valve_cda_m2=base.heat_out.outlet_valve_cda_m2*no/base.heat_out.bank.tube_count)
    lim=base.configuration.machine_volumes
    kin=FourStageVolumeKinematics(lim.small_cylinder,lim.large_cylinder,**{n:float(p[n]) for n in KIN_NAMES})
    return replace(base,heat_in=hi,heat_out=ho,kinematics=kin)

def _hardware_metrics(ex):
    meta=dict(ex.build().metadata); d=ex.bank.dimensions()
    return dict(tube_count=int(ex.bank.tube_count),tube_length_m=float(ex.bank.tube_length_m),inner_diameter_m=float(ex.bank.inner_diameter_m),outlet_valve_cda_m2=float(ex.outlet_valve_cda_m2),tube_flow_area_m2=float(d['tube_flow_area_m2']),working_gas_volume_m3=float(d['working_gas_volume_m3']),tube_gas_volume_m3=float(d['tube_gas_volume_m3']),header_gas_volume_m3=float(d['header_gas_volume_m3']),core_width_m=float(d['core_width_m']),core_height_m=float(d['core_height_m']),wall_capacity_j_k=float(meta['wall_capacity_j_k']),static_conductance_reference_w_k=float(meta['overall_static_conductance_w_k']),static_conductance_role=meta.get('static_conductance_role'))

def _uniform_state_and_hardware(design):
    wrapper=design.build()
    state=np.asarray(_uniform_wall_initial(design.configuration,wrapper),float)
    return state,{'H_i':_hardware_metrics(design.heat_in),'H_o':_hardware_metrics(design.heat_out)}

def _warm_start(source,target_uniform,target_hw):
    state=np.asarray(source['last_complete_state'],float).copy()
    if state.shape!=(10,): raise ValueError('Expected ten-state warm start.')
    sm=float(np.sum(state[:8:2])); tm=float(np.sum(target_uniform[:8:2]))
    if sm<=0 or tm<=0: raise ValueError('Warm-start gas inventory must be positive.')
    state[:8]*=tm/sm
    for i,side in ((8,'H_i'),(9,'H_o')):
        state[i]*=float(target_hw[side]['wall_capacity_j_k'])/float(source['hardware'][side]['wall_capacity_j_k'])
    return state

def _feasibility(result):
    reasons=[]; eta=result.get('indicated_thermal_efficiency'); power=result.get('indicated_power_w')
    if eta is None or not math.isfinite(eta): reasons.append('motor_efficiency_unavailable')
    if power is None or power<MINIMUM_POWER_W: reasons.append('minimum_motor_power')
    diag=result.get('diagnostics',{}); p=diag.get('pressure_extrema',{}); t=diag.get('temperature_extrema',{})
    pmax=max((float(x['maximum']) for x in p.values()),default=None); tmax=max((float(x['maximum']) for x in t.values()),default=None)
    if pmax is None or pmax>MAXIMUM_PRESSURE_PA: reasons.append('maximum_pressure')
    if tmax is None or tmax>MAXIMUM_TEMPERATURE_K: reasons.append('maximum_temperature')
    flow=result.get('maximum_absolute_mass_flow_kg_s')
    if flow is None or flow>MAXIMUM_ABSOLUTE_MASS_FLOW_KG_S: reasons.append('maximum_absolute_mass_flow')
    if result.get('validity',{}).get('verdict')!='valid': reasons.append('valid_thermodynamic_or_exchanger_model')
    result['maximum_pressure_pa']=pmax; result['maximum_temperature_k']=tmax
    return not reasons,reasons

def _parse_radii(text):
    rows=[]
    for item in text.split(','):
        if not item.strip(): continue
        a,b,c=item.split(':'); row=(float(a),float(b),int(c))
        if min(row)<=0: raise argparse.ArgumentTypeError('All radii must be positive.')
        rows.append(row)
    if not rows: raise argparse.ArgumentTypeError('At least one radius triplet required.')
    return tuple(rows)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--budget-seconds',type=float,default=21600.0); ap.add_argument('--evaluations',type=int,default=288)
    ap.add_argument('--candidate-seconds',type=float,default=180.0); ap.add_argument('--evaluations-per-radius',type=int,default=48)
    ap.add_argument('--radii',type=_parse_radii,default=DEFAULT_RADII); ap.add_argument('--seed',type=int,default=SOBOL_SEED)
    ap.add_argument('--output-directory',type=Path,default=DEFAULT_DIRECTORY); args=ap.parse_args()
    if min(args.budget_seconds,args.evaluations,args.candidate_seconds,args.evaluations_per_radius)<=0: raise ValueError('Budgets/counts must be positive.')

    definition,base=_load_basis(); seed_params=_load_center(); directory=args.output_directory; directory.mkdir(parents=True,exist_ok=True)
    hp=directory/'history.jsonl'; dp=directory/'definition.json'
    identity=dict(parameter_names=list(PARAMETER_NAMES),source_report=str(SOURCE_REPORT.relative_to(ROOT)),seed_parameters=seed_params,
        production_gas_model=dict(species='air',domain_policy='reject',thermal_entry=True,hydraulics='compressible_pressure_squared_with_variable_viscosity',internal_film='instantaneous_flow_dependent'),
        charge_pressure_pa=CHARGE_PRESSURE_PA,minimum_stage_fraction=MINIMUM_STAGE_FRACTION,tube_count_bounds=list(COUNT_BOUNDS),
        retained_H_i_tube_length_m=base.heat_in.bank.tube_length_m,retained_H_o_tube_length_m=base.heat_out.bank.tube_length_m,tube_inner_diameter_m=base.heat_in.bank.inner_diameter_m,
        valve_policy='outlet CdA scales linearly with tube count',radii=[list(x) for x in args.radii],evaluations_per_radius=args.evaluations_per_radius,sobol_seed=args.seed,
        strategy='Pure local adaptive 9D Sobol refinement around the balanced 7D point; no global proposals and no historical state warm starts.',
        feasibility='Power/pressure/temperature/flow constraints plus current thermodynamic+microtube validity; no legacy Re<2300 blanket rule.')
    if dp.exists() and json.loads(dp.read_text())!=identity: raise ValueError('Definition changed; use a new output directory.')
    if not dp.exists(): dp.write_text(json.dumps(identity,indent=2)+'\n')

    history=_load_history(hp); eligible=[r for r in history if _eligible(r)]; best=max(eligible,key=_efficiency) if eligible else None
    finished={r['candidate_id'] for r in history if r.get('result',{}).get('status')!='interrupted'}
    points=qmc.Sobol(d=9,scramble=True,seed=args.seed).random_base2(14)
    next_prop=max((r.get('proposal_index',-1) for r in history),default=-1)+1
    start=time.monotonic(); deadline=start+args.budget_seconds; completed=0

    def progress(_):
        if time.monotonic()>=deadline: raise IntegrationInterrupted('9D campaign wall-clock budget exhausted.')
        if time.monotonic()>=candidate_deadline: raise IntegrationInterrupted('Candidate wall-clock budget exhausted.')

    while completed<args.evaluations and time.monotonic()<deadline:
        proposal_index=next_prop; next_prop+=1
        if not history and proposal_index==0:
            p=dict(seed_params); kind='balanced_variable_gas_seed'; ridx=-1; tr=lr=nr=None
        else:
            n=sum(1 for r in history if r['kind']!='balanced_variable_gas_seed'); ridx=min(n//args.evaluations_per_radius,len(args.radii)-1)
            tr,lr,nr=args.radii[ridx]; center=best['parameters'] if best is not None else seed_params; u=points[proposal_index-1]
            p={name:float(center[name]+(2*u[i]-1)*tr) for i,name in enumerate(('t1','t2','t3'))}
            for i,name in enumerate(('a_l','b_l','a_s','b_s'),start=3): p[name]=float(center[name]+(2*u[i]-1)*lr)
            p['n_i']=int(round(center['n_i']+(2*u[7]-1)*nr)); p['n_o']=int(round(center['n_o']+(2*u[8]-1)*nr)); kind='adaptive_local_9d_variable_gas'
            if not _valid_parameters(p): continue
        cid=_candidate_id(p)
        if cid in finished: continue
        before=time.monotonic(); source=None; design=_build_design(base,p)
        try:
            uniform,hardware=_uniform_state_and_hardware(design)
            reusable=[r for r in history if r.get('result',{}).get('status')=='converged' and r.get('last_complete_state') is not None and r.get('hardware')]
            initial=uniform
            if reusable:
                target=np.asarray([p[n] for n in KIN_NAMES]+[p['n_i']/3708,p['n_o']/3708],float)
                source=min(reusable,key=lambda r:np.linalg.norm(target-np.asarray([r['parameters'][n] for n in KIN_NAMES]+[r['parameters']['n_i']/3708,r['parameters']['n_o']/3708],float)))
                initial=_warm_start(source,uniform,hardware)
            candidate_deadline=min(deadline,before+args.candidate_seconds)
            try: result,state=_evaluate(f'four_stage_hx9d_variable_gas_{proposal_index}',design,definition,initial_state=np.asarray(initial),progress_callback=progress)
            except MicrotubeDomainError as exc: result={'status':'invalid_exchanger','message':str(exc)};state=None
            except IntegrationInterrupted as exc: result={'status':'interrupted','message':str(exc)};state=None
            except (ValueError,RuntimeError) as exc: result={'status':'integration_failure','message':str(exc)};state=None
        except MicrotubeDomainError as exc:
            result={'status':'invalid_exchanger','message':str(exc)};state=None;hardware={'H_i':_hardware_metrics(design.heat_in),'H_o':_hardware_metrics(design.heat_out)}
        feasible,reasons=_feasibility(result) if result.get('status')=='converged' else (False,[])
        rec=dict(index=len(history),proposal_index=proposal_index,candidate_id=cid,kind=kind,radius_index=ridx,timing_radius=tr,level_radius=lr,tube_radius=nr,
            parameters={**{n:float(p[n]) for n in KIN_NAMES},'n_i':int(p['n_i']),'n_o':int(p['n_o'])},phase_degrees=_phase_degrees(p),symmetry=_symmetry_diagnostics(p),hardware=hardware,result=result,feasible=feasible,physical_constraint_failures=reasons,
            elapsed_seconds=time.monotonic()-before,warm_start_source=source['candidate_id'] if source is not None else None,last_complete_state=state.tolist() if state is not None else None)
        with hp.open('a') as f: f.write(json.dumps(rec)+'\n');f.flush()
        history.append(rec);completed+=1
        if result.get('status')!='interrupted': finished.add(cid)
        moved=False
        if feasible and (best is None or _efficiency(rec)>_efficiency(best)): best=rec;moved=True
        report=dict(best_feasible=best,seed=next((r for r in history if r['kind']=='balanced_variable_gas_seed'),None),attempted_total=len(history),completed_this_run=completed,requested_duration_seconds=args.budget_seconds,actual_duration_seconds=time.monotonic()-start,current_radius_index=ridx,current_radii=[tr,lr,nr],definition=identity)
        (directory/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(index=rec['index'],radius=[tr,lr,nr],status=result.get('status'),feasible=feasible,efficiency=result.get('indicated_thermal_efficiency'),power_W=result.get('indicated_power_w'),n_i=int(p['n_i']),n_o=int(p['n_o']),phases_deg=rec['phase_degrees'],center_moved=moved,best_efficiency=_efficiency(best) if best else None)),flush=True)
    print(f"Saved {directory/'report.json'}",flush=True)

if __name__=='__main__': main()
