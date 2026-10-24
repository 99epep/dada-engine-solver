"""Live production drawing and self-contained animated SVG frame export."""
import io
import re
import xml.etree.ElementTree as ET
from .machine_render import FRAME_COUNT, FRAME_DURATION_MS, prepare_machine_frames, draw_machine_frame

SVG='http://www.w3.org/2000/svg'
XLINK='http://www.w3.org/1999/xlink'


def prepare_animation(payload):
    from .artifacts import MechanismArtifact
    artifacts={side:MechanismArtifact.from_data(data) for side,data in payload['artifacts'].items()}
    geometry,data,low,high=prepare_machine_frames(artifacts,payload['cycle'],motor_operation=payload['motor_operation'])
    from .machine_render import visible_machine_geometry
    geometry=visible_machine_geometry(geometry,payload.get('show_mechanisms','both'))
    return dict(geometry=geometry,data=data,layout=payload['layout'],low=low,high=high,show_mechanisms=payload.get('show_mechanisms','both'))


def frame_machine_axes(ax, geometry, *, width_inches=11., max_height_inches=7.):
    """Fit the full-cycle envelope with small constant margins, without text."""
    from .machine_render import VERTICAL_MARGIN_PX
    figure=ax.figure;dpi=figure.dpi;margin=VERTICAL_MARGIN_PX
    figure.set_facecolor('#f0f0f0')
    xmin,xmax,ymin,ymax=geometry['limits']
    ratio=(xmax-xmin)/(ymax-ymin)
    inner_width=min(width_inches*dpi-2*margin,(max_height_inches*dpi-2*margin)*ratio)
    inner_height=inner_width/ratio
    width,height=inner_width+2*margin,inner_height+2*margin
    figure.set_size_inches(width/dpi,height/dpi,forward=True)
    ax.set_position((margin/width,margin/height,inner_width/width,inner_height/height))


def draw_animation_frame(ax, prepared, index, *, vector_conduits=True):
    draw_machine_frame(ax,prepared['geometry'],index,prepared['data'],prepared['layout'],prepared['low'],prepared['high'],show_mechanisms=prepared.get('show_mechanisms','both'),vector_conduits=vector_conduits)


def animated_machine_svg(prepared, *, vector_conduits=True):
    """Animate exact renderer frames with discrete SMIL; do not interpolate physics.

    Thermal conduits use native SVG linear gradients aligned to exact ports.
    WebP and live views use shared vector-strip fills of the same geometry.
    Readers without SMIL support display the first frame as a static fallback.
    """
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    import matplotlib as mpl
    ET.register_namespace('',SVG);ET.register_namespace('xlink',XLINK)
    figure=Figure(figsize=(12,5),dpi=100,facecolor='#f0f0f0')
    FigureCanvasAgg(figure);axes=figure.add_axes((.01,.01,.98,.98))
    frame_machine_axes(axes,prepared['geometry'],width_inches=12.,max_height_inches=5.)
    root=None
    duration=FRAME_COUNT*FRAME_DURATION_MS/1000
    key_times=';'.join(format(i/FRAME_COUNT,'.12g') for i in range(FRAME_COUNT+1))
    try:
        settings={'svg.hashsalt':'research-machine-animation'}
        if vector_conduits: settings['image.composite_image']=False
        with mpl.rc_context(settings):
            for index in range(FRAME_COUNT):
                draw_animation_frame(axes,prepared,index,vector_conduits=False)
                if vector_conduits:
                    from .svg_conduits import prepare_vector_conduits, replace_svg_conduits
                    ports=prepare_vector_conduits(axes)
                stream=io.BytesIO();figure.savefig(stream,format='svg',metadata={'Date':None})
                frame=ET.fromstring(stream.getvalue())
                if vector_conduits: replace_svg_conduits(frame,axes,ports)
                if root is None:
                    root=ET.Element('{'+SVG+'}svg',frame.attrib)
                    ET.SubElement(root,'{'+SVG+'}desc').text='Production machine animation: 66 frames, 85 ms/frame. Indicated motion and replay temperatures; no force/inertia simulation.'
                prefix=f'frame_{index}_'
                identifiers={node.attrib['id']:prefix+node.attrib['id'] for node in frame.iter() if 'id' in node.attrib}
                for node in frame.iter():
                    for key,value in list(node.attrib.items()):
                        if key=='id': node.set(key,identifiers[value])
                        elif value.startswith('#') and value[1:] in identifiers:
                            node.set(key,'#'+identifiers[value[1:]])
                        elif 'url(#' in value:
                            node.set(key,re.sub(r'url\(#([^)]+)\)',lambda m:'url(#'+identifiers.get(m[1],m[1])+')',value))
                group=ET.SubElement(root,'{'+SVG+'}g',id=f'frame_{index}',visibility='visible' if index==0 else 'hidden')
                for child in frame: group.append(child)
                values=['visible' if tick%FRAME_COUNT==index else 'hidden' for tick in range(FRAME_COUNT+1)]
                ET.SubElement(group,'{'+SVG+'}animate',attributeName='visibility',values=';'.join(values),
                    keyTimes=key_times,dur=f'{duration:.2f}s',repeatCount='indefinite',calcMode='discrete')
        return ET.tostring(root,encoding='utf-8',xml_declaration=True)
    finally: figure.clear()
