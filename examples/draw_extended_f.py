"""Draw the extended F linkage and explicitly locate its working chamber.

Chamber width and clearance are schematic. Linkage geometry is exact.
A volume increasing with slider coordinate requires the head on the lower
coordinate side; this shows the piston-rod entry through that head.
"""
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
    from matplotlib.patches import Polygon
    run=json.loads((ROOT/'outputs/small_coupler_f_extended.json').read_text())['best']
    assembly,phase=build(run['best']['parameters'],run['branch'],'F')
    axis=np.array([math.cos(assembly.slider.axis_angle),math.sin(assembly.slider.axis_angle)])
    normal=np.array([-axis[1],axis[0]])
    origin=np.array([assembly.slider.axis_origin_x,assembly.slider.axis_origin_y])
    pivot=np.array([assembly.loop.rocker_pivot_x,assembly.loop.rocker_pivot_y])
    angles=np.linspace(0,2*np.pi,7200,endpoint=False)
    states=[assembly.evaluate(np.array([math.cos(a),math.sin(a)]),np.array([-math.sin(a),math.cos(a)])) for a in angles]
    coordinates=np.array([s.coordinate for s in states]);low=coordinates.min();high=coordinates.max()
    sign=run['volume_sign'];width=.55;clearance=.12*(high-low)
    head=low-clearance if sign>0 else high+clearance
    ends=[low-clearance,high+clearance]
    fig,axes=plt.subplots(1,2,figsize=(11,9))
    for ax,index,title in zip(axes,[np.argmin(coordinates),np.argmax(coordinates)],['Minimum volume','Maximum volume'] if sign>0 else ['Maximum volume','Minimum volume']):
        s=states[index];pin=np.array(s.crank_pin);joint=np.array(s.coupler_joint);f=np.array(s.output_point);p=origin+axis*s.coordinate
        ax.set_aspect('equal');ax.axis('off')
        ax.plot([0,pin[0]],[0,pin[1]],'o-',color='#333333',lw=3)
        ax.plot(*np.array([pin,joint,pivot]).T,'o-',color='#2474b5',lw=2.5)
        ax.add_patch(Polygon(np.array([pin,joint,f]),closed=True,facecolor='#2474b5',alpha=.13))
        ax.plot(*np.array([pin,f,joint]).T,color='#2474b5',lw=1.7)
        ax.plot(*np.array([f,p]).T,'o-',color='#e87516',lw=2.5)
        ax.plot(*f,'o',color='#b31b1b',ms=7);ax.text(*(f+np.array([.15,0])),'F',color='#b31b1b')
        ax.plot(0,0,'ks',ms=6);ax.text(0,-.3,'A');ax.plot(*pivot,'ks',ms=6);ax.text(*(pivot+np.array([.1,0])),'D')
        for side in [-1,1]:ax.plot(*np.array([origin+axis*x+side*normal*width for x in ends]).T,color='#555555',lw=1.5)
        gas=np.array([origin+axis*head-normal*width,p-normal*width,p+normal*width,origin+axis*head+normal*width])
        ax.add_patch(Polygon(gas,facecolor='#efb757',alpha=.4))
        headpoint=origin+axis*head
        ax.plot(*np.array([headpoint-normal*width,headpoint+normal*width]).T,color='#333333',lw=4)
        ax.plot(*np.array([p-normal*width,p+normal*width]).T,color='#e87516',lw=6)
        ax.annotate('Head (schematic)',xy=headpoint,xytext=(headpoint[0]+1.4,headpoint[1]-.7),arrowprops={'arrowstyle':'->'},fontsize=9)
        ax.annotate('Piston',xy=p,xytext=(p[0]-2.,p[1]+.5),arrowprops={'arrowstyle':'->'},fontsize=9)
        ax.set_title(title);ax.set_xlim(-3,7);ax.set_ylim(-4,10)
    fig.suptitle('Extended coupler F — 0.883% position RMS\nWorking chamber on the mechanism-facing side of the piston',fontsize=14)
    fig.text(.5,.03,'Exact relative linkage scale; schematic bore and clearance.\nArticulated rod cannot pass directly through a conventional head seal.\nCrosshead / axial piston-rod arrangement remains to be designed.',ha='center',fontsize=10)
    fig.subplots_adjust(top=.87,bottom=.12,wspace=.15)
    fig.savefig(ROOT/'outputs/small_coupler_f_extended_chamber.png',dpi=150)
    plt.close(fig)

if __name__=='__main__':main()
