"""Stage K2: zoomed H_i/H_o exchanger-length search around the K1 ridge.

Keeps the Stage F6 champion motion and total working-gas mass fixed.
Varies only tube-length multipliers:
    k_i in [0.75, 1.55]
    k_o in [0.42, 0.78]

Default: 128 scrambled Sobol points, no controls.
Also exports the non-dominated efficiency/power Pareto front.
"""
from __future__ import annotations

from dataclasses import replace
import argparse, csv, json, math, os, tempfile
from pathlib import Path
import numpy as np
from scipy.stats import qmc

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.exchangers.microtube import MicrotubeExchanger
from compare_motor_motion_laws_stage7A5 import (
    ROOT, A5_CAMPAIGN, _candidate_design_and_mass, _evaluate,
    _same_inventory_design,
)
from optimize_motor_two_extrema_stageF6 import PARAMETERS as F6_PARAMETERS, _make_kinematics
from optimize_motor_exchanger_asymmetry_stageK1 import (
    _hardware_metrics, _scaled_exchanger, _result_row,
)

F6_JSON = ROOT / "outputs" / "motor_two_extrema_stageF6.json"
K1_JSON = ROOT / "outputs" / "motor_exchanger_asymmetry_stageK1.json"
OUTPUT_CSV = ROOT / "outputs" / "motor_exchanger_asymmetry_stageK2.csv"
OUTPUT_JSON = ROOT / "outputs" / "motor_exchanger_asymmetry_stageK2.json"
OUTPUT_PARETO_CSV = ROOT / "outputs" / "motor_exchanger_asymmetry_stageK2_pareto.csv"
OUTPUT_ETA_PLOT = ROOT / "outputs" / "motor_exchanger_asymmetry_stageK2_efficiency.png"
OUTPUT_POWER_PLOT = ROOT / "outputs" / "motor_exchanger_asymmetry_stageK2_power.png"
OUTPUT_PARETO_PLOT = ROOT / "outputs" / "motor_exchanger_asymmetry_stageK2_eta_vs_power.png"

BOUNDS = {"k_i": (0.75, 1.55), "k_o": (0.42, 0.78)}


def load_best(path):
    best = json.loads(path.read_text()).get("best_feasible")
    if best is None:
        raise RuntimeError(f"No feasible champion in {path}")
    return best


def scale(point):
    return {name: lo + float(x) * (hi - lo)
            for x, (name, (lo, hi)) in zip(point, BOUNDS.items(), strict=True)}


def write_csv(path, rows):
    if not rows: return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def load_checkpoint(path):
    if not path.exists(): return {}
    with path.open(newline="") as f: rows = list(csv.DictReader(f))
    if not rows or not {"index","label","k_i","k_o"}.issubset(rows[0]): return {}
    return {int(float(r["index"])): r for r in rows}


def restore(row):
    r = dict(row); r["index"] = int(float(r["index"]))
    for k in list(r):
        if k == "label": continue
        if k == "feasible" or k.startswith("constraint_"):
            r[k] = str(r[k]).lower() == "true"; continue
        if r[k] in ("", None): r[k] = None; continue
        try: r[k] = float(r[k])
        except (TypeError, ValueError): pass
    return r


def pareto_front(rows):
    f = [r for r in rows if r["feasible"] and r.get("indicated_thermal_efficiency") is not None and r.get("indicated_power_w") is not None]
    out = []
    for a in f:
        ea, pa = float(a["indicated_thermal_efficiency"]), float(a["indicated_power_w"])
        dominated = any(
            float(b["indicated_thermal_efficiency"]) >= ea and float(b["indicated_power_w"]) >= pa and
            (float(b["indicated_thermal_efficiency"]) > ea or float(b["indicated_power_w"]) > pa)
            for b in f if b is not a
        )
        if not dominated: out.append(a)
    return sorted(out, key=lambda r: float(r["indicated_power_w"]))


def setup_plt():
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "dada_solver_matplotlib"))
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_plane(path, rows, metric, title, cbar):
    f = [r for r in rows if r["feasible"] and r.get(metric) is not None]
    if not f: return
    plt = setup_plt(); fig, ax = plt.subplots(figsize=(8.5,7))
    z = np.array([float(r[metric]) for r in f]);
    if metric == "indicated_thermal_efficiency": z *= 100
    sc = ax.scatter([r["k_i"] for r in f],[r["k_o"] for r in f],c=z,s=50)
    fig.colorbar(sc,ax=ax).set_label(cbar)
    best=max(f,key=lambda r:float(r["indicated_thermal_efficiency"]))
    ax.scatter([best["k_i"]],[best["k_o"]],marker="x",s=140,label=f"best eta: ki={best['k_i']:.3f}, ko={best['k_o']:.3f}")
    ax.set(xlabel="k_i = H_i tube-length multiplier",ylabel="k_o = H_o tube-length multiplier",xlim=BOUNDS["k_i"],ylim=BOUNDS["k_o"],title=title)
    ax.grid(alpha=.25); ax.legend(loc="upper center",bbox_to_anchor=(.5,-.1),frameon=False)
    fig.tight_layout(rect=(0,.05,1,1)); fig.savefig(path,dpi=160); plt.close(fig)


def plot_pareto(rows, front):
    f=[r for r in rows if r["feasible"] and r.get("indicated_thermal_efficiency") is not None]
    if not f: return
    plt=setup_plt(); fig,ax=plt.subplots(figsize=(8.5,6.5))
    ax.scatter([r["indicated_power_w"] for r in f],[100*r["indicated_thermal_efficiency"] for r in f],s=35,alpha=.55,label="feasible K2 samples")
    if front:
        x=[r["indicated_power_w"] for r in front]; y=[100*r["indicated_thermal_efficiency"] for r in front]
        ax.plot(x,y,marker="o",label="non-dominated front")
        for r in front: ax.annotate(str(int(r["index"])),(r["indicated_power_w"],100*r["indicated_thermal_efficiency"]),xytext=(4,4),textcoords="offset points",fontsize=7)
    ax.set(xlabel="Indicated power (W)",ylabel="Indicated thermal efficiency (%)",title="Stage K2 — efficiency / indicated-power trade-off")
    ax.grid(alpha=.3); ax.legend(loc="upper center",bbox_to_anchor=(.5,-.12),frameon=False)
    fig.tight_layout(rect=(0,.05,1,1)); fig.savefig(OUTPUT_PARETO_PLOT,dpi=160); plt.close(fig)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--sobol-power",type=int,default=7); ap.add_argument("--seed",type=int,default=4252); ap.add_argument("--resume",action="store_true"); ap.add_argument("--csv",type=Path,default=OUTPUT_CSV); ap.add_argument("--output",type=Path,default=OUTPUT_JSON); args=ap.parse_args()
    definition=CampaignDefinition(A5_CAMPAIGN)
    base,total_mass,_saved,_a5,candidate_id=_candidate_design_and_mass(definition)
    if not isinstance(base.heat_in,MicrotubeExchanger) or not isinstance(base.heat_out,MicrotubeExchanger): raise TypeError("Expected microtube exchangers")
    f6=load_best(F6_JSON); k1=load_best(K1_JSON)
    f6v={name:float(f6[name]) for name in F6_PARAMETERS}
    kin=_make_kinematics(base.configuration.machine_volumes.small_cylinder,base.configuration.machine_volumes.large_cylinder,f6v)
    ref=_same_inventory_design(base,total_mass,kin)
    base_hi,base_ho=_hardware_metrics(ref.heat_in),_hardware_metrics(ref.heat_out)
    sampler=qmc.Sobol(d=2,scramble=True,seed=args.seed)
    candidates=[scale(p) for p in sampler.random_base2(args.sobol_power)]
    print("Stage K2: zoomed exchanger-asymmetry search")
    print(f"K1 benchmark eta={100*float(k1['indicated_thermal_efficiency']):.6f}% P={float(k1['indicated_power_w']):.3f} W ki={float(k1['k_i']):.4f} ko={float(k1['k_o']):.4f}")
    print(f"{len(candidates)} Sobol evaluations; bounds {BOUNDS}")
    checkpoint=load_checkpoint(args.csv) if args.resume else {}; rows=[]
    for index,v in enumerate(candidates):
        c=checkpoint.get(index)
        if c and c.get("label")=="sobol" and abs(float(c["k_i"])-v["k_i"])<1e-12 and abs(float(c["k_o"])-v["k_o"])<1e-12:
            row=restore(c); rows.append(row); print(f"{index:3d} CK {100*float(row['indicated_thermal_efficiency']):8.5f}% ki={v['k_i']:.4f} ko={v['k_o']:.4f}"); continue
        hi=_scaled_exchanger(ref.heat_in,v["k_i"]); ho=_scaled_exchanger(ref.heat_out,v["k_o"])
        hm_i,hm_o=_hardware_metrics(hi),_hardware_metrics(ho)
        result,_=_evaluate(f"exchanger_asymmetry_k2_{index}",replace(ref,heat_in=hi,heat_out=ho),definition)
        row=_result_row(index,v,result,hm_i,hm_o)
        row["total_tube_length_factor_vs_symmetric_baseline"]=.5*(v["k_i"]+v["k_o"])
        row["hi_to_ho_length_ratio"]=v["k_i"]/v["k_o"]
        rows.append(row); write_csv(args.csv,rows)
        eta=row["indicated_thermal_efficiency"]; p=row["indicated_power_w"]
        print(f"{index:3d} {'OK' if row['feasible'] else 'X ':2s} {100*eta if eta is not None else float('nan'):8.5f}% P={p if p is not None else float('nan'):7.3f} W ki={v['k_i']:.4f} ko={v['k_o']:.4f}")
    rows.sort(key=lambda r:int(r["index"])); write_csv(args.csv,rows)
    feasible=sorted([r for r in rows if r["feasible"] and r.get("indicated_thermal_efficiency") is not None],key=lambda r:float(r["indicated_thermal_efficiency"]),reverse=True)
    front=pareto_front(rows); write_csv(OUTPUT_PARETO_CSV,front)
    report={
      "experiment":"Stage K2 zoomed independent H_i/H_o tube-length search with Stage F6 champion motion fixed",
      "search":{"bounds":BOUNDS,"sobol_power":args.sobol_power,"sample_count":2**args.sobol_power,"seed":args.seed,"scramble":True,"control_count":0},
      "definition":{"scaled_geometry":["tube_length_m"],"same_total_working_gas_mass":True,"same_F6_motion":True,"external_air_mass_flow_unchanged":True,"external_aerodynamic_losses_in_objective":False},
      "F6_reference":{"candidate_id":candidate_id,"index":int(f6["index"]),"parameters":f6v,"indicated_thermal_efficiency":float(f6["indicated_thermal_efficiency"]),"indicated_power_w":float(f6["indicated_power_w"]),"total_mass_kg":total_mass},
      "K1_reference":k1,"baseline_hardware":{"H_i":base_hi,"H_o":base_ho},
      "best_feasible":feasible[0] if feasible else None,"pareto_front":front,"ranking":feasible,"all_rows":rows}
    args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(report,indent=2)+"\n")
    plot_plane(OUTPUT_ETA_PLOT,rows,"indicated_thermal_efficiency","Stage K2 — efficiency over exchanger-length asymmetry","Indicated thermal efficiency (%)")
    plot_plane(OUTPUT_POWER_PLOT,rows,"indicated_power_w","Stage K2 — indicated power over exchanger-length asymmetry","Indicated power (W)")
    plot_pareto(rows,front)
    print("\nTOP FEASIBLE")
    for rank,r in enumerate(feasible[:12],1): print(f"{rank:2d} eta={100*r['indicated_thermal_efficiency']:.6f}% P={r['indicated_power_w']:.3f} W ki={r['k_i']:.4f} ko={r['k_o']:.4f} ki/ko={r['hi_to_ho_length_ratio']:.3f}")
    print("\nPARETO FRONT")
    for r in front: print(f"idx={int(r['index']):3d} eta={100*r['indicated_thermal_efficiency']:.6f}% P={r['indicated_power_w']:.3f} W ki={r['k_i']:.4f} ko={r['k_o']:.4f}")

if __name__=="__main__": main()
