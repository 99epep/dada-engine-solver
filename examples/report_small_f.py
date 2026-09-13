"""Fine-screen F candidates and compare the best with extended E."""
import json
import argparse
import os
import tempfile
from pathlib import Path
import numpy as np
from search_small_rocker_e import ROOT, assess


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',type=Path,default=ROOT/'outputs/small_coupler_f_search.json')
    parser.add_argument('--plot',type=Path,default=ROOT/'outputs/small_e_f_comparison.png')
    args=parser.parse_args()
    path=args.input
    report=json.loads(path.read_text())
    grid=np.linspace(0,2*np.pi,7200,endpoint=False)
    for run in report['runs']:
        if not run['best']:continue
        try:
            r,_,_=assess(run['best']['parameters'],run['branch'],run['volume_sign'],grid,np.zeros_like(grid),np.zeros_like(grid),'F')
            run['fine_geometry_check']=dict(step_degrees=.05,metrics=r['metrics'],margins=r['margins'],passed=all(m>=0 for m in r['margins'].values()))
        except ValueError as exc:run['fine_geometry_check']=dict(passed=False,reason=str(exc))
    valid=[r for r in report['runs'] if r.get('fine_geometry_check',{}).get('passed')]
    report['best']=min(valid,key=lambda r:r['best']['score']) if valid else None
    path.write_text(json.dumps(report,indent=2)+'\n')
    if not valid:raise ValueError('No candidate passed the fine screen')
    e=json.loads((ROOT/'outputs/small_rocker_e_extended.json').read_text())['best']
    target=np.genfromtxt(ROOT/'outputs/motor_champion_motion_target.csv',delimiter=',',names=True)[:-1]
    os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'dada_solver_matplotlib'))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
    for family,run in [('E',e),('F',report['best'])]:
        _,q,dq=assess(run['best']['parameters'],run['branch'],run['volume_sign'],target['theta_rad'],target['small_fraction_0_1'],target['small_dq_dtheta_per_rad'],family)
        axes[0].plot(target['theta_deg'],q,label='Best '+family);axes[1].plot(target['theta_deg'],dq)
    for ax,key in zip(axes,['small_fraction_0_1','small_dq_dtheta_per_rad']):
        ax.plot(target['theta_deg'],target[key],'k--',label='Free-motion target');ax.grid(alpha=.2)
    axes[0].legend();axes[0].set_ylabel('Normalized volume');axes[1].set_ylabel('Derivative (rad⁻¹)');axes[1].set_xlabel('Study angle (deg)')
    fig.suptitle('Small cylinder — finite-rod E and F synthesis')
    fig.tight_layout();fig.savefig(args.plot,dpi=140)
    print(json.dumps(dict(elapsed_seconds=report['elapsed_seconds'],fine_valid=len(valid),best=report['best']['best']),indent=2))

if __name__=='__main__':main()
