"""Animate one finite-rod F candidate at true relative linkage scale."""
import json
import math
import os
import tempfile
from pathlib import Path
import numpy as np
from search_small_rocker_e import ROOT, build


def main():
    os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'dada_solver_matplotlib'))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    source=ROOT/'outputs/small_coupler_f_search.json'
    run=json.loads(source.read_text())['best']
    assembly,phase=build(run['best']['parameters'],run['branch'],'F')
    axis=np.array([math.cos(assembly.slider.axis_angle),math.sin(assembly.slider.axis_angle)])
    normal=np.array([-axis[1],axis[0]])
    origin=np.array([assembly.slider.axis_origin_x,assembly.slider.axis_origin_y])
    pivot=np.array([assembly.loop.rocker_pivot_x,assembly.loop.rocker_pivot_y])
    frames=[]
    for theta in np.linspace(0,2*np.pi,120,endpoint=False):
        a=-theta+phase
        s=assembly.evaluate(np.array([math.cos(a),math.sin(a)]),np.array([-math.sin(a),math.cos(a)]))
        frames.append(np.array([s.crank_pin,s.coupler_joint,s.output_point,origin+axis*s.coordinate]))
    frames=np.array(frames)
    points=np.vstack([frames.reshape(-1,2),pivot,[0,0]])
    lo=points.min(axis=0);hi=points.max(axis=0);pad=.1*max(hi-lo)
    fig,ax=plt.subplots(figsize=(10,7));ax.set_aspect('equal');ax.axis('off')
    ax.set_xlim(lo[0]-pad,hi[0]+pad);ax.set_ylim(lo[1]-pad,hi[1]+pad)
    fig.suptitle('Small cylinder — coupler-point F and finite piston rod',fontsize=14)
    ax.set_title('All lengths at the same scale • slow motor-direction playback',fontsize=10)
    ax.plot(0,0,'ks',ms=6);ax.plot(*pivot,'ks',ms=6)
    ax.text(0,-.35,'A',ha='center');ax.text(pivot[0],pivot[1]-.35,'D',ha='center')
    crank,=ax.plot([],[],'o-',color='#333333',lw=3)
    links,=ax.plot([],[],'o-',color='#2474b5',lw=2.5,label='Four-bar links')
    triangle,=ax.plot([],[],color='#2474b5',lw=1.5,label='Rigid coupler extension to F')
    rod,=ax.plot([],[],'o-',color='#e87516',lw=2.5,label='Finite piston rod')
    marker,=ax.plot([],[],'o',color='#b31b1b',ms=8)
    label=ax.text(0,0,'F',color='#b31b1b',fontsize=12)
    piston,=ax.plot([],[],color='#e87516',lw=6)
    travel=frames[:,3]@axis
    start=frames[np.argmin(travel),3];end=frames[np.argmax(travel),3]
    width=.25
    for sign in [-1,1]:
        guide=np.array([start,end])+sign*width*normal
        ax.plot(*guide.T,color='#777777',lw=1)
    ax.legend(loc='lower left',fontsize=9,frameon=False)
    fig.text(.5,.025,'Schematic piston guide • stroke / envelope ≈ 12.4% • cylinder body not sized',ha='center',fontsize=10)
    fig.subplots_adjust(top=.86,bottom=.12)
    def update(i):
        pin,joint,f,p=frames[i]
        crank.set_data([0,pin[0]],[0,pin[1]])
        links.set_data(*np.array([pin,joint,pivot]).T)
        triangle.set_data(*np.array([pin,f,joint]).T)
        rod.set_data(*np.array([f,p]).T)
        marker.set_data([f[0]],[f[1]]);label.set_position(f+np.array([.15,.15]))
        piston.set_data(*np.array([p-width*normal,p+width*normal]).T)
    animation=FuncAnimation(fig,update,frames=120,interval=1000/30)
    animation.save(ROOT/'outputs/small_coupler_f.gif',writer=PillowWriter(fps=30),dpi=110)
    update(25);fig.savefig(ROOT/'outputs/small_coupler_f_preview.png',dpi=110);plt.close(fig)

if __name__=='__main__':main()
