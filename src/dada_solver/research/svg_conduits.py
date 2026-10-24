"""Vector thermal conduits with exact exchanger ports and small seam overlaps."""
import xml.etree.ElementTree as ET
import numpy as np
from matplotlib.transforms import Affine2D
from matplotlib.colors import to_hex
from . import transfer_block as block

SVG='http://www.w3.org/2000/svg'


def prepare_vector_conduits(axes):
    """Remove the preview offset locally, then capture exact production ports.

    Only temporary export artists are changed. Read both cut-port edges from the
    approved serpentine path instead of approximating them with its bbox centre.
    """
    ports=[]
    for patch in axes.patches:
        if patch.get_path() is not block.HX_PATH: continue
        patch.set_transform(Affine2D().translate(0,-block.HX_Y_OFFSET/block.SCALE)+patch.get_transform())
        for x in (block.HX_VISIBLE_X0_SRC,block.HX_VISIBLE_X1_SRC):
            ys=block.HX_PATH.vertices[np.isclose(block.HX_PATH.vertices[:,0],x,rtol=0,atol=1e-9),1]
            if not len(ys): raise ValueError('Cannot identify the exact exchanger port.')
            ports.append(patch.get_transform().transform([[x,float(ys.min())],[x,float(ys.max())]]))
    for index,image in enumerate(axes.images): image.set_gid(f'vector_conduit_{index}')
    return ports


def conduit_polygon(image, ports, endpoints):
    """Display-space quad; overlap joins by 2.5% of the pipe thickness."""
    x0,x1,y0,y1=image.get_extent()
    polygon=image.get_transform().transform([[x0,y1],[x1,y1],[x1,y0],[x0,y0]])
    if image.get_clip_path() is None:
        for upper,lower in ((0,3),(1,2)):
            center=float(np.mean(polygon[[upper,lower],1]))
            height=abs(polygon[upper,1]-polygon[lower,1])
            matches=[port for port in ports if abs(polygon[upper,0]-port[0,0])<1e-6
                     and abs(center-float(np.mean(port[:,1])))<height/2]
            if matches:
                port=min(matches,key=lambda p:abs(center-float(np.mean(p[:,1]))))
                polygon[upper,1],polygon[lower,1]=sorted(port[:,1],reverse=True)
                polygon[[upper,lower],0]+=(1 if upper==1 else -1)*height*.025
            elif any(abs(polygon[upper,0]-end)<1e-6 for end in endpoints):
                polygon[[upper,lower],0]+=(1 if upper==1 else -1)*height*.025
    return polygon


def conduit_endpoints(axes):
    rectangular=[image for image in axes.images if image.get_clip_path() is None]
    xs=[image.get_transform().transform((image.get_extent()[i],0))[0]
        for image in rectangular for i in (0,1)]
    return (min(xs),max(xs)) if xs else ()


def draw_vector_conduits(axes):
    """Render thermal fills as vector strips in Agg and other Matplotlib backends.

    SVG animation exports use exact native linearGradient fills instead. Strips
    share edges and disable antialiasing to avoid rasterization hairline seams.
    """
    from matplotlib.collections import PolyCollection
    ports=prepare_vector_conduits(axes)
    endpoints=conduit_endpoints(axes)
    for image in list(axes.images):
        quad=axes.transData.inverted().transform(conduit_polygon(image,ports,endpoints))
        colors=np.asarray(image.get_array(),dtype=float)
        t=np.linspace(0,1,257)
        top=quad[0]+t[:,None]*(quad[1]-quad[0])
        bottom=quad[3]+t[:,None]*(quad[2]-quad[3])
        polygons=np.stack((top[:-1],top[1:],bottom[1:],bottom[:-1]),axis=1)
        rgb=colors[0,0]+((t[:-1]+t[1:])/2)[:,None]*(colors[0,-1]-colors[0,0])
        collection=PolyCollection(polygons,facecolors=rgb,edgecolors='none',
            antialiaseds=False,transform=axes.transData,zorder=image.get_zorder())
        collection.set_clip_path(image.get_clip_path())
        collection.set_clip_box(image.get_clip_box())
        axes.add_collection(collection,autolim=False)
        image.remove()


def _coordinate(value):
    # Match Matplotlib's SVG coordinate precision at component boundaries.
    return f'{float(value):.6f}'.rstrip('0').rstrip('.') or '0'


def replace_svg_conduits(frame, axes, ports):
    """Replace tagged gradient images in-place, retaining clips and draw order."""
    factor=72/axes.figure.dpi;height=axes.figure.bbox.height
    def project(points):
        result=np.asarray(points,dtype=float).copy()
        result[:,1]=height-result[:,1]
        return result*factor
    display_ports=ports
    endpoints=conduit_endpoints(axes)
    images={image.get_gid():image for image in axes.images}
    defs=ET.SubElement(frame,'{'+SVG+'}defs')
    replaced=0
    for parent in list(frame.iter()):
        for index,node in enumerate(list(parent)):
            if node.tag!='{'+SVG+'}image': continue
            image=images.get(node.get('id'))
            if image is None: raise ValueError('Unexpected untagged SVG image; vector replacement refused.')
            colors=np.asarray(image.get_array(),dtype=float)
            if colors.ndim!=3 or colors.shape[2]!=3:
                raise ValueError('Expected a horizontal RGB thermal gradient.')
            expected=np.linspace(colors[0,0],colors[0,-1],colors.shape[1])
            if not np.allclose(colors,expected[None,:,:],atol=1e-12):
                raise ValueError('Unexpected non-linear gradient; vector replacement refused.')
            polygon=project(conduit_polygon(image,display_ports,endpoints))
            identifier=node.get('id')+'_gradient'
            gradient=ET.SubElement(defs,'{'+SVG+'}linearGradient',id=identifier,
                gradientUnits='userSpaceOnUse',x1=_coordinate(polygon[0,0]),x2=_coordinate(polygon[1,0]),
                y1='0',y2='0',attrib={'color-interpolation':'sRGB'})
            ET.SubElement(gradient,'{'+SVG+'}stop',offset='0',attrib={'stop-color':to_hex(colors[0,0])})
            ET.SubElement(gradient,'{'+SVG+'}stop',offset='1',attrib={'stop-color':to_hex(colors[0,-1])})
            path='M '+' L '.join(_coordinate(x)+' '+_coordinate(y) for x,y in polygon)+' Z'
            attributes=dict(id=node.get('id'),d=path,fill=f'url(#{identifier})')
            if node.get('opacity') is not None: attributes['opacity']=node.get('opacity')
            replacement=ET.Element('{'+SVG+'}path',attributes)
            parent.remove(node);parent.insert(index,replacement);replaced+=1
    if replaced!=len(images): raise ValueError('SVG gradient replacement was incomplete.')
