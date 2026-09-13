"""Show the chamber-away F winner at its volume extrema, with actual linkage scale."""
import json
import os
import math
import tempfile
from pathlib import Path
import numpy as np
from search_small_rocker_e import ROOT, build


def main():
    report=json.loads((ROOT/'outputs/small_f_chamber_away_search.json').read_text())
    run=report['best'];assembly,phase=build(run['best']['parameters'],run['branch'],'F',-1)
    if run['volume_sign']!=1:raise ValueError('Expected chamber-away volume direction')
    axis=np.array([math.cos(assembly.slider.axis_angle),math.sin(assembly.slider.axis_angle)])
    normal=np.array([-axis[1],axis[0]])
    # Rotate the whole drawing only, keeping the positive guide direction up.
    rotation=np.array([[axis[1],-axis[0]],[axis[0],axis[1]]])
    origin=np.array([assembly.slider.axis_origin_x,assembly.slider.axis_origin_y])
    pivot=np.array([assembly.loop.rocker_pivot_x,assembly.loop.rocker_pivot_y])
    theta=np.linspace(0,2*np.pi,7200,endpoint=False)
    states=[assembly.evaluate(np.array([math.cos(t),math.sin(t)]),np.array([-math.sin(t),math.cos(t)])) for t in theta]
    coordinates=np.array([s.coordinate for s in states]);lo=coordinates.min();hi=coordinates.max()
    stroke=hi-lo;width=.3*stroke;head=lo-.1*stroke
    os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'dada_solver_matplotlib'))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    fig,axes=plt.subplots(1,2,figsize=(10,8))
    points=np.array([rotation@np.array(p) for s in states for p in [s.crank_pin,s.coupler_joint,s.output_point,origin+axis*s.coordinate]])
    lower=points.min(axis=0)-stroke;upper=points.max(axis=0)+stroke
    for ax,index,title in zip(axes,[np.argmin(coordinates),np.argmax(coordinates)],['Minimum volume — piston away','Maximum volume — piston toward mechanism']):
        s=states[index]
        pin,joint,f,p,hp,d=[rotation@np.array(x) for x in [s.crank_pin,s.coupler_joint,s.output_point,origin+axis*s.coordinate,origin+axis*head,pivot]]
        n=rotation@normal
        ax.set_aspect('equal');ax.axis('off');ax.set_title(title,fontsize=10)
        ax.plot([0,pin[0]],[0,pin[1]],'o-',color='#333333',lw=3)
        ax.plot(*np.array([pin,joint,d]).T,'o-',color='#2474b5',lw=2.5)
        ax.add_patch(Polygon(np.array([pin,joint,f]),facecolor='#2474b5',alpha=.15))
        ax.plot(*np.array([pin,f,joint]).T,color='#2474b5')
        ax.plot(*np.array([f,p]).T,'o-',color='#e87516',lw=2.5)
        ax.plot(*f,'o',color='#b31b1b');ax.text(*(f+np.array([.15,0])),'F')
        ax.plot(0,0,'ks');ax.plot(*d,'ks')
        for side in [-1,1]:
            ax.plot(*np.array([rotation@(origin+axis*c)+side*n*width for c in [head,hi+.1*stroke]]).T,color='#555555')
        ax.add_patch(Polygon(np.array([hp-width*n,p-width*n,p+width*n,hp+width*n]),facecolor='#efb757',alpha=.4))
        ax.plot(*np.array([hp-width*n,hp+width*n]).T,color='#333333',lw=4)
        ax.plot(*np.array([p-width*n,p+width*n]).T,color='#e87516',lw=6)
        ax.annotate('Head',xy=hp,xytext=hp+np.array([width*2,-stroke*.5]),arrowprops={'arrowstyle':'->'},fontsize=9)
        ax.set_xlim(lower[0],upper[0]);ax.set_ylim(lower[1],upper[1])
    fig.suptitle('F refitted with working chamber away from mechanism')
    fig.text(.5,.035,'Actual relative linkage scale; schematic bore and clearance.\nFull body collisions and bearing/guide loads remain unverified.',ha='center',fontsize=10)
    fig.subplots_adjust(top=.9,bottom=.12);fig.savefig(ROOT/'outputs/small_f_chamber_away_layout.png',dpi=140)

if __name__=='__main__':main()
