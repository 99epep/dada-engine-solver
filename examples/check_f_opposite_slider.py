"""Move the extended F piston to the opposite rod branch without rotating links."""
from dataclasses import replace, asdict
import json
import math
import os
import tempfile
from pathlib import Path
import numpy as np
from search_small_rocker_e import ROOT, build
from search_compact_motor_motion import evaluate


def main():
    run=json.loads((ROOT/'outputs/small_coupler_f_extended.json').read_text())['best']
    upper,phase=build(run['best']['parameters'],run['branch'],'F')
    lower=replace(upper,slider=replace(upper.slider,assembly_branch=-1))
    raw=np.genfromtxt(ROOT/'outputs/motor_champion_motion_target.csv',delimiter=',',names=True)[:-1]
    theta=raw['theta_rad'];qtarget=raw['small_fraction_0_1'];dtarget=raw['small_dq_dtheta_per_rad']
    curves=[];report=dict(description='Opposite slider branch; unchanged four-bar, rod length, axis and phase',
                         volume_increases_with_coordinate=True,phase_rad=phase,candidates={})
    for name,assembly in [('upper',upper),('lower',lower)]:
        q,dq,metrics=evaluate(assembly,1.,1,phase,theta)
        metrics['stroke_over_crank']=metrics.pop('stroke_m');metrics['envelope_diagonal_over_crank']=metrics.pop('envelope_diagonal_m')
        report['candidates'][name]=dict(position_rms=float(np.sqrt(np.mean((q-qtarget)**2))),
            derivative_rms=float(np.sqrt(np.mean((dq-dtarget)**2))),metrics=metrics,assembly=asdict(assembly))
        curves.append((q,dq))
    report['upper_lower_position_rms']=float(np.sqrt(np.mean((curves[0][0]-curves[1][0])**2)))
    (ROOT/'outputs/f_opposite_slider.json').write_text(json.dumps(report,indent=2)+'\n')
    os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'dada_solver_matplotlib'))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
    for name,(q,dq) in zip(['Original upper piston','Opposite lower piston'],curves):
        axes[0].plot(raw['theta_deg'],q,label=name);axes[1].plot(raw['theta_deg'],dq)
    axes[0].plot(raw['theta_deg'],qtarget,'k--',label='Free-motion target');axes[1].plot(raw['theta_deg'],dtarget,'k--')
    axes[0].legend();axes[0].set_ylabel('Normalized volume');axes[1].set_ylabel('Derivative (rad⁻¹)');axes[1].set_xlabel('Study angle (deg)')
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('Same F mechanism — opposite piston position')
    fig.tight_layout();fig.savefig(ROOT/'outputs/f_opposite_slider_comparison.png',dpi=140);plt.close(fig)

    axis=np.array([math.cos(lower.slider.axis_angle),math.sin(lower.slider.axis_angle)])
    normal=np.array([-axis[1],axis[0]]);origin=np.array([lower.slider.axis_origin_x,lower.slider.axis_origin_y])
    pivot=np.array([lower.loop.rocker_pivot_x,lower.loop.rocker_pivot_y])
    def state(t):
        a=t+phase
        return lower.evaluate(np.array([math.cos(a),math.sin(a)]),np.array([-math.sin(a),math.cos(a)]))
    states=[state(t) for t in theta];coords=np.array([s.coordinate for s in states])
    lo,hi=coords.min(),coords.max();head=lo-.12*(hi-lo);width=.5
    fig,axes=plt.subplots(1,2,figsize=(10,9))
    for ax,i,title in zip(axes,[np.argmin(coords),np.argmax(coords)],['Minimum volume','Maximum volume']):
        s=states[i];pin=np.array(s.crank_pin);joint=np.array(s.coupler_joint);f=np.array(s.output_point);p=origin+axis*s.coordinate;hp=origin+axis*head
        ax.set_aspect('equal');ax.axis('off');ax.set_title(title)
        ax.plot([0,pin[0]],[0,pin[1]],'o-',color='#333333',lw=3)
        ax.plot(*np.array([pin,joint,pivot]).T,'o-',color='#2474b5',lw=2.5)
        ax.add_patch(Polygon(np.array([pin,joint,f]),facecolor='#2474b5',alpha=.2))
        ax.plot(*np.array([pin,f,joint]).T,color='#2474b5')
        ax.plot(*np.array([f,p]).T,'o-',color='#e87516',lw=2.5)
        ax.plot(*f,'o',color='#b31b1b');ax.text(*(f+np.array([.15,0])),'F')
        ax.plot(0,0,'ks');ax.plot(*pivot,'ks')
        for sign in [-1,1]:ax.plot(*np.array([origin+axis*c+sign*normal*width for c in [head,hi+.1*(hi-lo)]]).T,color='#777777')
        ax.add_patch(Polygon(np.array([hp-width*normal,p-width*normal,p+width*normal,hp+width*normal]),facecolor='#efb757',alpha=.4))
        ax.plot(*np.array([hp-width*normal,hp+width*normal]).T,color='#333333',lw=4)
        ax.plot(*np.array([p-width*normal,p+width*normal]).T,color='#e87516',lw=6)
        ax.annotate('Head / working gas below piston',xy=hp,xytext=(-4,-10.8),fontsize=8,arrowprops={'arrowstyle':'->'})
        ax.set_xlim(-5,5);ax.set_ylim(-12,2)
    fig.suptitle('F mechanism unchanged — piston on the lower side\nWorking chamber away from the mechanism; no rod crossing the head')
    fig.text(.5,.03,'Exact relative linkage scale; schematic cylinder width and clearance.\nBody clearances, rod sweep inside skirt and piston guidance remain unverified.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.87,bottom=.12);fig.savefig(ROOT/'outputs/f_lower_piston_layout.png',dpi=140);plt.close(fig)
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
