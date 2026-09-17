"""Run sequential whole-machine temperature branches.

This wrapper repeatedly launches optimize_motor_temperature_point.py. Each
temperature starts from the completed champion of the previous temperature in
its own branch. Different branches restart from the common 300 K-delta-T
reference champion.

Every requested delta T must have an explicit power floor.

Example after validating the 270 K point:
  PYTHONPATH=src python3 examples/run_motor_temperature_map.py \
    --branches "300,270,240,210;300,330,360" \
    --power-floors "300:40,270:40,240:25,210:10,330:40,360:40"

The 240/210 values above are examples only; choose them deliberately from the
preceding results before running the lower-temperature branch.
"""

from __future__ import annotations
import argparse, json, os, subprocess, sys
from pathlib import Path

REFERENCE_REPORT = Path("outputs/motor_four_stage_hx9d_variable_gas/report.json")
MAP_ROOT = Path("outputs/motor_temperature_map")
POINT_SCRIPT = Path("examples/optimize_motor_temperature_point.py")


def parse_branches(text: str) -> list[list[float]]:
    result=[]
    for branch in text.split(';'):
        values=[float(x.strip()) for x in branch.split(',') if x.strip()]
        if values: result.append(values)
    if not result: raise argparse.ArgumentTypeError('At least one branch is required.')
    return result

def parse_power_floors(text: str) -> dict[float,float]:
    result={}
    for item in text.split(','):
        if not item.strip(): continue
        dt,p=item.split(':')
        result[float(dt)]=float(p)
    return result

def tag(dt: float) -> str:
    return f"dT_{dt:g}K".replace('.','p')

def load_complete(path: Path):
    if not path.exists(): return None
    data=json.loads(path.read_text())
    return data if data.get('campaign_complete') and data.get('champion') is not None else None

def summary_row(report: dict) -> dict:
    c=report['champion']; r=c['result']; p=c['parameters']; ph=c['phase_degrees']
    return dict(
        delta_t_k=report['target_delta_t_k'],
        hot_source_temperature_k=report['hot_source_temperature_k'],
        power_floor_w=report['power_floor_w'],
        indicated_thermal_efficiency=r['indicated_thermal_efficiency'],
        carnot_efficiency=c['carnot_efficiency'],
        fraction_of_carnot=c['fraction_of_carnot'],
        indicated_power_w=r['indicated_power_w'],
        heat_input_w=r['heat_input_w'],
        heat_out_w=r['heat_out_w'],
        gas_inventory_kg=c['derived_geometry']['gas_inventory_kg'],
        swept_ratio=p['swept_ratio'],
        n_i=p['n_i'], length_i_m=p['length_i_m'],
        n_o=p['n_o'], length_o_m=p['length_o_m'],
        low_pressure_exchange_deg=ph['low_pressure_exchange_deg'],
        compression_deg=ph['compression_deg'],
        high_pressure_exchange_deg=ph['high_pressure_exchange_deg'],
        expansion_deg=ph['expansion_deg'],
        a_l=p['a_l'], b_l=p['b_l'], a_s=p['a_s'], b_s=p['b_s'],
    )

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--branches',type=parse_branches,required=True,
                    help='Semicolon-separated branches, e.g. 270,240,210;330,360')
    ap.add_argument('--power-floors',type=parse_power_floors,required=True,
                    help='Explicit deltaT:W mapping for every requested point.')
    ap.add_argument('--reference-report',type=Path,default=REFERENCE_REPORT)
    ap.add_argument('--map-root',type=Path,default=MAP_ROOT)
    ap.add_argument('--thermo-evaluations',type=int,default=96)
    ap.add_argument('--motion-evaluations',type=int,default=96)
    ap.add_argument('--thermo-evaluations-per-radius',type=int,default=24)
    ap.add_argument('--motion-evaluations-per-radius',type=int,default=24)
    ap.add_argument('--thermo-budget-seconds',type=float,default=10800)
    ap.add_argument('--motion-budget-seconds',type=float,default=10800)
    ap.add_argument('--candidate-seconds',type=float,default=240)
    ap.add_argument('--seed',type=int,default=26091730)
    args=ap.parse_args()

    requested=[dt for branch in args.branches for dt in branch]
    missing=[dt for dt in requested if dt not in args.power_floors]
    if missing:
        raise ValueError(f'Missing explicit power floor(s) for delta T: {missing}')
    if not args.reference_report.exists(): raise FileNotFoundError(args.reference_report)
    if not POINT_SCRIPT.exists(): raise FileNotFoundError(POINT_SCRIPT)

    args.map_root.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy()
    src=str(Path('src').resolve())
    env['PYTHONPATH']=src+(os.pathsep+env['PYTHONPATH'] if env.get('PYTHONPATH') else '')

    rows=[]
    for branch_index,branch in enumerate(args.branches):
        source=args.reference_report
        print(f"=== branch {branch_index+1}: {' -> '.join(f'{x:g}' for x in branch)} K ===",flush=True)
        for dt in branch:
            out=args.map_root/tag(dt)
            report_path=out/'report.json'
            existing=load_complete(report_path)
            if existing is not None:
                print(f"deltaT={dt:g} K already complete; reusing champion.",flush=True)
                rows.append(summary_row(existing)); source=report_path; continue

            cmd=[
                sys.executable,str(POINT_SCRIPT),
                '--delta-t-k',f'{dt:g}',
                '--power-floor-w',f'{args.power_floors[dt]:g}',
                '--source-report',str(source),
                '--output-directory',str(out),
                '--thermo-evaluations',str(args.thermo_evaluations),
                '--motion-evaluations',str(args.motion_evaluations),
                '--thermo-evaluations-per-radius',str(args.thermo_evaluations_per_radius),
                '--motion-evaluations-per-radius',str(args.motion_evaluations_per_radius),
                '--thermo-budget-seconds',str(args.thermo_budget_seconds),
                '--motion-budget-seconds',str(args.motion_budget_seconds),
                '--candidate-seconds',str(args.candidate_seconds),
                '--seed',str(args.seed),
            ]
            print(' '.join(cmd),flush=True)
            completed=subprocess.run(cmd,env=env)
            if completed.returncode!=0:
                raise SystemExit(completed.returncode)
            report=load_complete(report_path)
            if report is None:
                print(f"deltaT={dt:g} K is not complete. Rerun this wrapper to resume it.",flush=True)
                summary=args.map_root/'adapted_summary.json'
                summary.write_text(json.dumps({'results':rows},indent=2)+'\n')
                return
            rows.append(summary_row(report))
            source=report_path

    # Deduplicate points that may occur in several branches, preserving the
    # latest completed copy.
    merged={float(row['delta_t_k']):row for row in rows}
    ordered=[merged[k] for k in sorted(merged)]
    summary=args.map_root/'adapted_summary.json'
    summary.write_text(json.dumps({'results':ordered},indent=2)+'\n')
    print(f"Saved {summary}",flush=True)

if __name__=='__main__':
    main()
