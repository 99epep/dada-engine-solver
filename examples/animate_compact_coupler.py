"""Animate the finite-rod coupler candidates in one shared-crank frame.

Cylinder outlines are schematic guides, not sized collision envelopes.
All linkage lengths and piston travel retain their relative geometric scale.
"""
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np
from dada_solver.four_bar import FourBarLoop, CouplerOutputPoint, SliderConstraint, FourBarSliderAssembly

ROOT = Path(__file__).resolve().parents[1]


def main():
    os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'dada_solver_matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    report=json.loads((ROOT/'outputs/compact_motor_coupler_search.json').read_text())
    sides=[]
    phases=np.linspace(0,2*math.pi,120,endpoint=False)
    all_points=[np.zeros(2)]
    for name in ['small','large']:
        record=report['sides'][name]; raw=record['assembly']
        assembly=FourBarSliderAssembly(FourBarLoop(**raw['loop']),CouplerOutputPoint(**raw['output']),SliderConstraint(**raw['slider']))
        offset=record['crank_angle_offset']
        rotation=np.array([[math.cos(offset),math.sin(offset)],[-math.sin(offset),math.cos(offset)]])
        axis=np.array([math.cos(assembly.slider.axis_angle),math.sin(assembly.slider.axis_angle)])
        origin=np.array([assembly.slider.axis_origin_x,assembly.slider.axis_origin_y])
        pivot=rotation@np.array([assembly.loop.rocker_pivot_x,assembly.loop.rocker_pivot_y])
        frames=[]
        for phase in phases:
            angle=-record['crank_direction']*phase+offset
            pin=record['crank_radius']*np.array([math.cos(angle),math.sin(angle)])
            state=assembly.evaluate(pin,record['crank_radius']*np.array([-math.sin(angle),math.cos(angle)]))
            piston=origin+axis*state.coordinate
            points=np.array([rotation@np.array(p) for p in [state.crank_pin,state.coupler_joint,state.output_point,piston]])
            frames.append(points); all_points.extend(points)
        all_points.append(pivot)
        sides.append((name,record,pivot,rotation@axis,np.array(frames)))
    fig,ax=plt.subplots(figsize=(10,7))
    ax.set_aspect('equal'); ax.axis('off')
    points=np.array(all_points); low=points.min(axis=0); high=points.max(axis=0)
    pad=.12*np.max(high-low)
    ax.set_xlim(low[0]-pad,high[0]+pad); ax.set_ylim(low[1]-pad,high[1]+pad)
    fig.suptitle('DADA motor — coupler output with finite piston rods',fontsize=15)
    ax.set_title('Shared crank • actual relative linkage scale • slow playback',fontsize=10)
    crank,=ax.plot([],[],'o-',color='#333333',lw=3,ms=5)
    ax.plot(0,0,'ks',ms=7)
    moving=[]
    for (name,record,pivot,axis,frames),color in zip(sides,['#2474b5','#e87516']):
        normal=np.array([-axis[1],axis[0]])
        travel=frames[:,3,:]@axis
        first=frames[np.argmin(travel),3]; last=frames[np.argmax(travel),3]
        width=.045*np.max(high-low)
        # Guide width is a display choice; no piston bore is implied.
        for sign in [-1,1]:
            guide=np.array([first,last])+sign*normal*width
            ax.plot(*guide.T,color=color,alpha=.3,lw=1.5)
        ax.plot(*pivot,'ks',ms=6)
        links,=ax.plot([],[],'o-',color=color,lw=2.5,ms=4,label=name.capitalize()+' cylinder linkage')
        plate,=ax.plot([],[],color=color,lw=1.7)
        rod,=ax.plot([],[],'o-',color=color,lw=2,ms=4)
        piston,=ax.plot([],[],color=color,lw=5)
        moving.append((frames,pivot,normal,width,links,plate,rod,piston))
    ax.legend(loc='lower left',frameon=False)
    fig.text(.5,.035,'Schematic piston guides • cylinder dimensions and collisions not evaluated',ha='center',fontsize=10)
    fig.subplots_adjust(top=.87,bottom=.1)

    def update(i):
        common=sides[0][4][i,0]
        assert np.allclose(common,sides[1][4][i,0]), 'Crank pins must coincide'
        crank.set_data([0,common[0]],[0,common[1]])
        for frames,pivot,normal,width,links,plate,rod,piston in moving:
            pin,joint,output,slider=frames[i]
            links.set_data(*np.array([pin,joint,pivot]).T)
            plate.set_data(*np.array([pin,output,joint]).T)
            rod.set_data(*np.array([output,slider]).T)
            piston.set_data(*np.array([slider-normal*width,slider+normal*width]).T)
    animation=FuncAnimation(fig,update,frames=120,interval=1000/30)
    output=ROOT/'outputs/compact_motor_coupler.gif'
    animation.save(output,writer=PillowWriter(fps=30),dpi=110)
    update(25); fig.savefig(ROOT/'outputs/compact_motor_coupler_preview.png',dpi=110)
    plt.close(fig)
    print(output)


if __name__=='__main__': main()
