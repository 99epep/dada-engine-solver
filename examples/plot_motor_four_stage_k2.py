"""Plot and export the current best four-stage K2 motion in motor-time order."""
import csv
import json
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dada_solver.four_stage_kinematics import FourStageVolumeKinematics
from dada_solver.geometry import CylinderVolumeLimits
from compare_motor_motion_laws_stage7A5 import ROOT


def main():
    directory=ROOT/'outputs/motor_four_stage_k2'
    report=json.loads((directory/'report.json').read_text())
    best=report['best_feasible']
    if best is None: raise ValueError('No feasible candidate to plot.')
    ref=json.loads((ROOT/'outputs/motor_champion_sixbar_k2.json').read_text())
    limits=ref['configuration']['machine_volumes']
    kin=FourStageVolumeKinematics(CylinderVolumeLimits(**limits['small_cylinder']),
        CylinderVolumeLimits(**limits['large_cylinder']),**best['parameters'])
    target=np.genfromtxt(ROOT/'outputs/motor_champion_motion_target.csv',delimiter=',',names=True)
    origin=target['theta_rad'][target['large_fraction_0_1'].argmax()]
    t=np.linspace(0,1,1441)
    fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
    rows=[]
    for i,side in enumerate(['small','large']):
        region=['S','L'][i]
        v=np.array([getattr(kin,f'{side}_cylinder_volume')(-2*math.pi*x) for x in t])
        bounds=getattr(kin,f'{side}_volume_limits')
        q=(v-bounds.minimum)/bounds.swept
        free=np.interp((origin-2*math.pi*t)%(2*math.pi),target['theta_rad'],target[f'{side}_fraction_0_1'])
        axes[i].plot(t*360,free,'--',label='Free target (L maximum at time zero)',color='gray')
        axes[i].plot(t*360,q,label=f'Four-stage {region}',color=['tab:blue','tab:orange'][i])
        for knot in [kin.t1,kin.t2,kin.t3]: axes[i].axvline(knot*360,color='gray',lw=.6,ls=':')
        axes[i].set_ylabel(f'{region} swept-volume fraction');axes[i].legend();axes[i].grid(alpha=.2)
        rows.append((v,q))
    axes[-1].set_xlabel('Forward motor-cycle angle (degrees)')
    eta=best['result']['indicated_thermal_efficiency']*100
    fig.suptitle(f'Four-stage K2 search: best observed efficiency {eta:.3f}% (candidate {best["index"]})')
    fig.tight_layout()
    fig.savefig(directory/'best_motion.png',dpi=160);fig.savefig(directory/'best_motion.svg')
    with (directory/'best_motion.csv').open('w',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['cycle_fraction','motor_angle_deg','time_s','S_volume_m3','L_volume_m3','S_swept_fraction','L_swept_fraction'])
        for j,x in enumerate(t): writer.writerow([x,x*360,x*.5,rows[0][0][j],rows[1][0][j],rows[0][1][j],rows[1][1][j]])
    print(json.dumps(best['parameters'],indent=2))


if __name__=='__main__': main()
