"""Experimental SVG conduits keep exact production ports and baseline rendering."""
import io
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import numpy as np
import pytest
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from dada_solver.research import transfer_block as block
from dada_solver.research.svg_conduits import prepare_vector_conduits, replace_svg_conduits, SVG


@pytest.mark.parametrize('layout',['UU','DD'])
def test_svg_gradients_and_exact_ports_without_preview_offset(layout):
    figure=Figure(figsize=(8,4),dpi=100);FigureCanvasAgg(figure);axes=figure.subplots()
    spec=block.LAYOUTS[layout]
    for branch,y,direction in [('Ho',.8,'left'),('Hi',-.8,'right')]:
        block.draw_row(axes,y=y,name=branch,**spec[branch],valve_direction=direction,
            valve_open=branch=='Hi',left_color=(.2,0,.8),right_color=(.8,0,.2),hx_color=(.5,0,.5))
    axes.set_xlim(1,5.5);axes.set_ylim(-2,2);axes.axis('off')
    figure.canvas.draw()
    hx=[patch for patch in axes.patches if patch.get_path() is block.HX_PATH]
    before=[patch.get_transform().transform([block.HX_VISIBLE_X0_SRC,block.HX_PORT_CENTER_Y_SRC]) for patch in hx]
    old_offset=block.HX_Y_OFFSET
    ports=prepare_vector_conduits(axes)
    assert block.HX_Y_OFFSET==old_offset==.004
    assert len(ports)==4
    after=[patch.get_transform().transform([block.HX_VISIBLE_X0_SRC,block.HX_PORT_CENTER_Y_SRC]) for patch in hx]
    for a,b in zip(before,after):
        assert a[0]==pytest.approx(b[0])
        assert a[1]>b[1]
    with mpl.rc_context({'image.composite_image':False}):
        stream=io.BytesIO();figure.savefig(stream,format='svg')
    frame=ET.fromstring(stream.getvalue())
    assert frame.findall('.//{'+SVG+'}image')
    replace_svg_conduits(frame,axes,ports)
    assert not frame.findall('.//{'+SVG+'}image')
    gradients=frame.findall('.//{'+SVG+'}linearGradient')
    assert len(gradients)==len(axes.images)
    assert all(len(g.findall('{'+SVG+'}stop'))==2 for g in gradients)
    # Exact projected exchanger port edges occur as conduit polygon vertices.
    factor=72/figure.dpi
    expected=np.array(ports).copy();expected[:,:,1]=figure.bbox.height-expected[:,:,1];expected*=factor
    paths=[p for p in frame.findall('.//{'+SVG+'}path') if p.get('id','').startswith('vector_conduit_')]
    import re
    vertices=np.concatenate([np.array([float(x) for x in re.findall(r'[-+]?\d*\.?\d+',p.get('d'))]).reshape(-1,2) for p in paths])
    for point in expected.reshape(-1,2):
        distance=vertices-point
        # Tiny axial overlap hides antialias seams without moving port heights.
        assert np.min(np.abs(distance[:,1])+np.maximum(np.abs(distance[:,0])-.2,0))<2e-6
    # A valve on the opposite row can share a port's X: never bridge the rows.
    for path in paths:
        polygon=np.array([float(x) for x in re.findall(r'[-+]?\d*\.?\d+',path.get('d'))]).reshape(-1,2)
        height=abs(polygon[0,1]-polygon[3,1])
        assert abs(polygon[0,1]-polygon[1,1])<height*.02
        assert abs(polygon[3,1]-polygon[2,1])<height*.02
    # A fresh normal draw still applies the approved default offset.
    figure.clear();axes=figure.subplots()
    block.add_exchanger(axes,1.,0.,(.5,0,.5))
    patch=next(p for p in axes.patches if p.get_path() is block.HX_PATH)
    assert patch.get_transform().transform([block.HX_VISIBLE_X0_SRC,block.HX_PORT_CENTER_Y_SRC])[1]==pytest.approx(axes.transData.transform([1.,.004])[1])


def test_cli_vector_conduits_flag_is_passed_to_viewer(monkeypatch):
    from dada_solver.research.cli import main
    from tests.test_research_inspection_tools import graph_data
    seen={}
    monkeypatch.setattr('dada_solver.research.report.compare',lambda *a,**k:graph_data())
    from dada_solver.research.graph_view import browse_graphs
    def browse(data,**kwargs):
        seen.update(kwargs);return browse_graphs(data,**kwargs)
    monkeypatch.setattr('dada_solver.research.graph_view.browse_graphs',browse)
    assert main(['graphs','campaign','--plots','positions','--svg-vector-conduits','--no-show'])==0
    assert seen['svg_vector_conduits'] is True


@pytest.mark.parametrize('visible',['both','none','small','large'])
def test_vector_machine_visibility_keeps_thermal_geometry(tmp_path,visible):
    from dada_solver.research.machine_render import candidate_mechanisms, prepare_machine_frames, draw_machine_frame, visible_machine_geometry
    from dada_solver.research.schema import load_study
    from dada_solver.research.presets import initialize_kinematics
    from tests.test_research_machine_render import cycle
    study=load_study(initialize_kinematics(tmp_path/'study.toml','slider_crank','slider_crank'))
    artifacts,_=candidate_mechanisms(study,{})
    geometry,data,low,high=prepare_machine_frames(artifacts,cycle())
    framed=visible_machine_geometry(geometry,visible)
    assert framed['sides'] is geometry['sides']
    figure=Figure(figsize=(12,5),dpi=100);FigureCanvasAgg(figure);axes=figure.subplots()
    draw_machine_frame(axes,framed,10,data,'UU',low,high,show_mechanisms=visible)
    assert not axes.images
    single=visible in ('small','large')
    assert bool(axes.collections)==(not single)
    assert len([patch for patch in axes.patches if patch.get_path() is block.HX_PATH])==(0 if single else 2)
    walls=[line for line in axes.lines if line.get_zorder()==6 and line.get_color()=='#111111']
    assert len(walls)==(3 if single else 6)
    if single:
        chosen=geometry['sides'][visible][1]
        other=geometry['sides']['large' if visible=='small' else 'small'][1]
        assert not any(np.isclose(line.get_xdata(),other['x_outer']).any() for line in walls)
        assert any(np.isclose(line.get_xdata(),chosen['x_outer']).any() for line in walls)
    # Only selected linkages contribute articulation circles; cylinders stay.
    from matplotlib.patches import Circle
    circles=[patch for patch in axes.patches if isinstance(patch,Circle)]
    assert bool(circles)==(visible!='none')
    figure.savefig(tmp_path/(visible+'.png'))
    figure.clear()


@pytest.mark.parametrize('command',['report','graphs'])
def test_cli_mechanism_visibility_forwarded(monkeypatch,tmp_path,command):
    from dada_solver.research.cli import main
    from tests.test_research_inspection_tools import graph_data
    seen={}
    def compare(*args,**kwargs):
        seen.update(kwargs)
        return graph_data()
    monkeypatch.setattr('dada_solver.research.report.compare',compare)
    monkeypatch.setattr('dada_solver.research.report.render_html',lambda *a,**k:None)
    monkeypatch.setattr('dada_solver.research.report.text_report',lambda *a,**k:'')
    args=[command,'campaign','--plots','positions','--show-mechanisms','small']
    args+=['--no-show'] if command=='graphs' else ['--html',str(tmp_path/'report.html')]
    assert main(args)==0
    assert seen['show_mechanisms']=='small'
