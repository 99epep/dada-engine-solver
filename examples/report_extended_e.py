"""Compare extended E fits and independently screen geometry every 0.05 degree."""
import json
import os
import tempfile
import numpy as np
from search_small_rocker_e import ROOT, assess


def main():
    path=ROOT/'outputs/small_rocker_e_extended.json'
    report=json.loads(path.read_text())
    old=json.loads((ROOT/'outputs/small_rocker_e_search.json').read_text())
    raw=np.genfromtxt(ROOT/'outputs/motor_champion_motion_target.csv',delimiter=',',names=True)[:-1]
    theta=np.linspace(0,2*np.pi,7200,endpoint=False)
    # Extra grid is used only for geometry; target comparison retains exact CSV samples.
    for run in report['runs']:
        if not run['best']:continue
        try:
            result,_,_=assess(run['best']['parameters'],run['branch'],run['volume_sign'],theta,np.zeros_like(theta),np.zeros_like(theta))
            run['fine_geometry_check']=dict(step_degrees=.05,metrics=result['metrics'],margins=result['margins'],passed=all(v>=0 for v in result['margins'].values()))
        except ValueError as exc:
            run['fine_geometry_check']=dict(passed=False,reason=str(exc))
    valid=[r for r in report['runs'] if r.get('fine_geometry_check',{}).get('passed')]
    if not valid:raise ValueError('No candidate passed the fine geometry check')
    report['best']=min(valid,key=lambda r:r['best']['score'])
    path.write_text(json.dumps(report,indent=2)+'\n')
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'outputs/.matplotlib'))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
    for label,run in [('Previous best E',old['best']),('Extended best E',report['best'])]:
        _,q,dq=assess(run['best']['parameters'],run['branch'],run['volume_sign'],raw['theta_rad'],raw['small_fraction_0_1'],raw['small_dq_dtheta_per_rad'])
        axes[0].plot(raw['theta_deg'],q,label=label)
        axes[1].plot(raw['theta_deg'],dq)
    for ax,key in zip(axes,['small_fraction_0_1','small_dq_dtheta_per_rad']):
        ax.plot(raw['theta_deg'],raw[key],'k--',label='Free-motion target',lw=2);ax.grid(alpha=.2)
    axes[0].legend();axes[0].set_ylabel('Normalized volume');axes[1].set_ylabel('Derivative (rad⁻¹)')
    axes[1].set_xlabel('Study angle (deg)');fig.suptitle('Small-cylinder E synthesis — extended search')
    fig.tight_layout();fig.savefig(ROOT/'outputs/small_rocker_e_extended_comparison.png',dpi=140)
    print(json.dumps(dict(runs=len(report['runs']),fine_valid=len(valid),elapsed_seconds=report['elapsed_seconds'],best=report['best']['best']),indent=2))

if __name__=='__main__':main()
