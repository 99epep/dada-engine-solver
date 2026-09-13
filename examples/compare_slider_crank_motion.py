"""Compare offset finite-rod slider-crank laws with the free motor target.

Independent synchronized crank axes; no cylinder layout or thermodynamics.
Lengths are normalized by crank radius. The piston coordinate is
x = cos(a) + sqrt(l*l - (e-sin(a))**2).
Volume direction is explicitly selectable; a piston-rod extension changes the
layout, not this kinematic equation.
"""
import json
import math
import os
from pathlib import Path
import tempfile
import numpy as np
from scipy.optimize import differential_evolution, brentq

ROOT=Path(__file__).resolve().parents[1]


def motion(theta, ratio, offset, phase, inverted=False):
    if ratio <= 1+abs(offset):
        raise ValueError('Rod cannot complete a nonsingular revolution')
    def raw(a):
        a=np.asarray(a)
        transverse=offset-np.sin(a)
        root=np.sqrt(ratio**2-transverse**2)
        return np.cos(a)+root, -np.sin(a)+transverse*np.cos(a)/root
    # Find physical travel extrema from analytic velocity, not plot samples.
    grid=np.linspace(0,2*np.pi,129)
    roots=[]
    for a,b in zip(grid[:-1],grid[1:]):
        da,db=raw(a)[1],raw(b)[1]
        if abs(da)<1e-14: roots.append(a)
        if da*db<0: roots.append(brentq(lambda t:raw(t)[1],a,b))
    values=[float(raw(a)[0]) for a in roots]
    lo,hi=min(values),max(values)
    x,dx=raw(np.asarray(theta)+phase)
    q=(x-lo)/(hi-lo); dq=dx/(hi-lo)
    if inverted: q=1-q; dq=-dq
    return q,dq


def main():
    os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'dada_solver_matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    target=np.genfromtxt(ROOT/'outputs/motor_champion_motion_target.csv',delimiter=',',names=True)
    theta=target['theta_rad']; sample=np.arange(0,len(theta)-1,8)
    fig,axes=plt.subplots(2,2,figsize=(12,7),sharex=True)
    report=dict(description='Independent crank axes; kinematic fit only',derivative_weight=.1,
                rod_ratio_bounds=[2.,10.],offset_over_rod_bound=.4,seed=27,sides={})
    for col,name in enumerate(['small','large']):
        tq=target[name+'_fraction_0_1']; tdq=target[name+'_dq_dtheta_per_rad']
        results={}
        for family,inverted,free_offset in [('conventional',True,False),('inverted',False,False),('inverted_offset',False,True)]:
            def decode(v):
                ratio=v[0]; offset=v[1]*ratio if free_offset else 0.
                return ratio,offset,v[-1]
            def objective(v):
                ratio,offset,phase=decode(v)
                q,dq=motion(theta[sample],ratio,offset,phase,inverted)
                return float(np.mean((q-tq[sample])**2)+.1*np.mean((dq-tdq[sample])**2))
            bounds=[(2.,10.)]+([(-.4,.4)] if free_offset else [])+[(-np.pi,np.pi)]
            fit=differential_evolution(objective,bounds,seed=27,popsize=8,maxiter=80,tol=1e-7,polish=True)
            ratio,offset,phase=decode(fit.x)
            q,dq=motion(theta,ratio,offset,phase,inverted)
            result=dict(rod_over_crank=ratio,offset_over_crank=offset,phase_rad=phase,volume_increases_with_coordinate=not inverted,
                position_rms=float(np.sqrt(np.mean((q[:-1]-tq[:-1])**2))),
                derivative_rms=float(np.sqrt(np.mean((dq[:-1]-tdq[:-1])**2))),
                maximum_rod_angle_degrees=math.degrees(math.asin((1+abs(offset))/ratio)))
            results[family]=result
            label={'conventional':'Volume decreases with x, centered','inverted':'Volume increases with x, centered','inverted_offset':'Volume increases with x, offset'}[family]
            axes[0,col].plot(np.degrees(theta),q,label=label,lw=1.6)
            axes[1,col].plot(np.degrees(theta),dq,lw=1.6)
        # Explicit retained ideal sinusoid benchmark, in the study convention.
        from dada_solver.kinematics import HarmonicVolumeKinematics
        from dada_solver.geometry import CylinderVolumeLimits
        limits=CylinderVolumeLimits(1.,2.)
        harmonic=HarmonicVolumeKinematics(limits,limits,math.radians(249))
        q=np.array([getattr(harmonic,name+'_cylinder_volume')(t)-1 for t in theta])
        dq=np.array([getattr(harmonic,name+'_cylinder_volume_derivative')(t) for t in theta])
        axes[0,col].plot(np.degrees(theta),q,':',color='gray',label='Harmonic reference, 249°')
        axes[1,col].plot(np.degrees(theta),dq,':',color='gray')
        axes[0,col].plot(np.degrees(theta),tq,'k--',label='Free-motion target',lw=2)
        axes[1,col].plot(np.degrees(theta),tdq,'k--',lw=2)
        axes[0,col].set_title(name.capitalize()+' cylinder')
        axes[1,col].set_xlabel('Study cycle angle (deg)')
        report['sides'][name]=results
    axes[0,0].set_ylabel('Normalized volume');axes[1,0].set_ylabel('Normalized volume derivative (rad⁻¹)')
    axes[0,0].legend(fontsize=8)
    for ax in axes.flat: ax.set_xlim(0,360);ax.grid(alpha=.25)
    fig.suptitle('Finite-rod slider-crank laws versus the free-motion target\nIndependent phase and rod ratio per cylinder; no thermodynamic optimization')
    fig.tight_layout();fig.savefig(ROOT/'outputs/slider_crank_target_comparison.png',dpi=150)
    (ROOT/'outputs/slider_crank_target_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
