"""Plot the compact finite-rod synthesis against its free-motion target."""
import json
import os
import tempfile
from pathlib import Path
import numpy as np
from search_compact_motor_motion import ROOT, candidate, evaluate
from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_model


def main():
    os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'dada_solver_matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    report=json.loads((ROOT/'outputs/compact_motor_motion_search.json').read_text())
    kin=build_model(load_simulation_configuration(ROOT/report['reference'])).kinematics
    target=np.genfromtxt(ROOT/report['target'],delimiter=',',names=True)
    coupler_path=ROOT/'outputs/compact_motor_coupler_search.json'
    coupler_report=json.loads(coupler_path.read_text()) if coupler_path.exists() else None
    angles=target['theta_rad']
    fig,axes=plt.subplots(2,2,figsize=(10,6),sharex=True)
    for col,name in enumerate(['small','large']):
        record=report['sides'][name]
        base=getattr(kin,name+'_assembly')
        fitted,phase=candidate(base,record['best']['parameters'])
        fits=[(base,0,'Reference mechanism'),(fitted,phase,'Rocker + finite rod')]
        if coupler_report and 'best' in coupler_report['sides'][name]:
            coupler,cp=candidate(base,coupler_report['sides'][name]['best']['parameters'],'coupler')
            fits.append((coupler,cp,'Coupler + finite rod'))
        for assembly,offset,label in fits:
            q,dq,_=evaluate(assembly,kin.crank_radius,kin.crank_direction,kin.crank_angle_offset+offset,angles)
            if not getattr(kin,name+'_volume_increases_with_coordinate'): q=1-q; dq=-dq
            axes[0,col].plot(np.degrees(angles),q,label=label)
            axes[1,col].plot(np.degrees(angles),dq)
        axes[0,col].plot(np.degrees(angles),target[name+'_fraction_0_1'],'k--',label='Free-motion target')
        axes[1,col].plot(np.degrees(angles),target[name+'_dq_dtheta_per_rad'],'k--')
        axes[0,col].set_title(name.capitalize()+' cylinder')
        axes[1,col].set_xlabel('Study angle (deg)')
    axes[0,0].set_ylabel('Normalized volume (0–1)')
    axes[1,0].set_ylabel('Normalized volume derivative (rad⁻¹)')
    axes[0,0].legend(fontsize=8)
    for ax in axes.flat: ax.grid(alpha=.2); ax.set_xlim(0,360)
    fig.suptitle('Compact reference family — kinematic synthesis only')
    fig.tight_layout()
    fig.savefig(ROOT/'outputs/compact_motor_motion_comparison.png',dpi=140)


if __name__=='__main__': main()
