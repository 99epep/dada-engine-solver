"""Whole-machine glyph placement, replay state and report integration."""
import base64
import io
import json
import math
from pathlib import Path
import numpy as np
import pytest
from PIL import Image
from matplotlib.figure import Figure
from matplotlib.patches import PathPatch
from dada_solver.research import transfer_block as block
from dada_solver.research.machine_render import (
    candidate_mechanisms, machine_geometry, render_machine_webp,
    frame_series, layout_for_configuration, MACHINE_SCALE, EXCHANGER_GAP,
)
from dada_solver.research.presets import initialize_kinematics
from dada_solver.research.schema import load_study
from dada_solver.research.plot_data import candidate_plots


def assert_clear_cylinder_walls(states,links,cylinder,side):
    from dada_solver.research.machine_render import _linkage_hits_cylinder,_stroke_world,CYLINDER_LINE_WIDTHS
    outward=-1 if side=='small' else 1
    head=cylinder['x_inner']-outward*_stroke_world(CYLINDER_LINE_WIDTHS[side])
    assert not _linkage_hits_cylinder(states,links,0.,head,cylinder['x_outer'],cylinder['width'],cylinder['clearance'])


def cycle():
    angle = np.linspace(0.,2*math.pi,181)
    series = {name:(320+40*np.sin(angle+i)).tolist() for i,name in enumerate(
        ('S_temperature_K','L_temperature_K','Hi_temperature_K','Ho_temperature_K'))}
    series['Ho_to_S_valve_open'] = (angle<math.pi).astype(int).tolist()
    series['Hi_to_L_valve_open'] = (angle>=math.pi).astype(int).tolist()
    return dict(angle_rad=angle.tolist(),series=series)


@pytest.mark.parametrize('layout,ho,hi',[
    ('UU',('hx','valve'),('valve','hx')),
    ('DD',('valve','hx'),('hx','valve')),
    ('UD',('hx','valve'),('hx','valve')),
    ('DU',('valve','hx'),('valve','hx')),
])
def test_component_order_direction_and_layering(layout,ho,hi):
    assert block.LAYOUTS[layout]['Ho']['order']==ho
    assert block.LAYOUTS[layout]['Hi']['order']==hi
    fig=Figure();ax=fig.subplots()
    placed=block.PlacedAxes(ax,-2.16,2.16,.4)
    for branch,direction in [('Ho','left'),('Hi','right')]:
        end=block.draw_row(placed,y=0.,name=branch,**block.LAYOUTS[layout][branch],
            valve_direction=direction,valve_open=False,
            left_color=(0.,0.,1.),right_color=(1.,0.,0.),hx_color=(.5,0.,.5))
        assert placed.placement.transform((1.20,0.))==pytest.approx((-2.16,.4))
        assert placed.placement.transform((end,0.))==pytest.approx((2.16,.4))
        vertices,base,tip,_,_=block.valve_geometry(0.,0.,direction)
        assert (tip<base)==(direction=='left')
    hx=[p for p in ax.patches if isinstance(p,PathPatch)]
    assert len(hx)==2 and all(p.get_zorder()==1 for p in hx)
    assert all(im.get_zorder()>1 for im in ax.images)
    assert max(p.get_zorder() for p in ax.patches)==7
    # Transform applies equally to glyph y dimensions and conduit thickness.
    assert placed.placement.transform((1.20,block.PIPE_H))[1]-.4==pytest.approx(block.PIPE_H*placed.factor)
    assert block.HX_Y_OFFSET==.004
    assert block.HX_VISIBLE_X1_SRC==81.115979
    fig.clear()


def test_open_valve_has_gradient_triangle_without_bar():
    fig=Figure();ax=fig.subplots()
    block.draw_valve(ax,0.,0.,'left',True,(0.,0.,1.),(1.,0.,0.))
    assert len(ax.images)==1 and ax.images[0].get_zorder()==6
    assert len(ax.patches)==1 and len(ax.patches[0].get_xy())==4
    fig.clear()


@pytest.mark.parametrize('family',['slider_crank','four_bar','six_bar'])
def test_production_geometry_keeps_original_head_spacing(tmp_path,family):
    study=load_study(initialize_kinematics(tmp_path/'study.toml',family,family))
    artifacts,_=candidate_mechanisms(study,{})
    angles=np.linspace(0.,2*math.pi,48,endpoint=False)
    data=machine_geometry(artifacts,angles)
    assert data['heads']==dict(small=-2.16,large=2.16)
    assert data['heads']['large']-data['heads']['small']==pytest.approx(MACHINE_SCALE*EXCHANGER_GAP/2)
    for side,(states,cylinder) in data['sides'].items():
        assert cylinder['x_inner']==data['heads'][side]
        assert all(np.isfinite(list(state.values())).all() for state in states)
        assert abs(cylinder['x_outer']-cylinder['x_inner'])==pytest.approx(1.10*np.ptp([s['P'][0] for s in states]))
    assert np.ptp([s['P'][0] for s in data['sides']['small'][0]])==pytest.approx(np.ptp([s['P'][0] for s in data['sides']['large'][0]]))


def test_binary_states_are_not_interpolated_as_lift():
    query,values=frame_series(cycle(),47)
    assert set(values['Ho_to_S_valve_open'])=={0.,1.}
    assert values['Ho_to_S_valve_open'][query<math.pi].all()
    assert not values['Ho_to_S_valve_open'][query>math.pi].any()


@pytest.mark.parametrize('layout',['UU','DD'])
def test_webp_is_animated_fixed_size_and_temperature_colored(tmp_path,layout):
    study=load_study(initialize_kinematics(tmp_path/'study.toml','six_bar','six_bar'))
    artifacts,_=candidate_mechanisms(study,{})
    payload,metadata=render_machine_webp(artifacts,cycle(),layout=layout,frames=6)
    with Image.open(io.BytesIO(payload)) as image:
        assert image.format=='WEBP' and image.is_animated and image.n_frames==6
        assert image.size==(1200,metadata['height']) and image.info['loop']==0
        assert metadata['vertical_margin_px']==20 and 40<image.height<=500
        pixels=np.asarray(image.convert('RGB'))
        assert np.any((pixels[:,:,0]>180)&(pixels[:,:,2]<80))
    assert metadata['layout']==layout
    from dada_solver.research.machine_render import _stroke_world, CYLINDER_LINE_WIDTHS
    for side,sign in (('small',-1),('large',1)):
        assert metadata['transfer_block_endpoints'][side]==pytest.approx(
            metadata['cylinder_inner_faces'][side]+sign*.5*_stroke_world(CYLINDER_LINE_WIDTHS[side]))
    assert metadata['cylinder_inner_faces']['small']>-2.16 and metadata['cylinder_inner_faces']['large']<2.16
    assert len(payload)<500_000


def test_report_embeds_whole_machine_and_requests_one_replay(tmp_path,monkeypatch):
    study=load_study(initialize_kinematics(tmp_path/'study.toml','slider_crank','four_bar'))
    calls=[]
    def replay(*args,**kwargs):
        calls.append(1);return dict(machine_cycle=cycle())
    monkeypatch.setattr('dada_solver.research.plot_data.replay_thermal',replay)
    monkeypatch.setattr('dada_solver.campaign.evaluator.MachineEvaluator.evaluate',lambda *a,**k:pytest.fail('campaign reevaluation'))
    result=candidate_plots(study,dict(candidate_id='a'*64,physical={}),('mechanisms',),cache_directory=tmp_path/'cache')
    assert len(calls)==1
    motion=result['mechanisms']
    assert motion['media_type']=='image/webp'
    with Image.open(io.BytesIO(base64.b64decode(motion['data_uri'].split(',')[1]))) as image:
        assert image.n_frames==66
    from tests.test_research_cockpit import campaign
    from dada_solver.research import report
    c,_=campaign(tmp_path/'campaign',1)
    data=report.inspect(c.directory)
    data['selected']=data['records']
    data['selected'][0]['plots']=result
    page=report.render_html(data,tmp_path/'report.html').read_text()
    assert motion['data_uri'] in page
    assert 'image.src=motion.data_uri' in page and 'image.alt=' in page
    assert 'Linkage joint positions' not in page


def test_abstract_campaign_does_not_replay_for_mechanisms(tmp_path,monkeypatch):
    study=load_study(initialize_kinematics(tmp_path/'study.toml','harmonic','harmonic'))
    monkeypatch.setattr('dada_solver.research.plot_data.replay_thermal',lambda *a,**k:pytest.fail('replay'))
    result=candidate_plots(study,dict(candidate_id='a'*64,physical={}),('mechanisms',))
    assert 'complete physical mechanism' in result['mechanisms']['unavailable']


def test_continuous_opening_fades_bar_without_changing_shapes():
    fig=Figure();ax=fig.subplots()
    block.draw_valve(ax,0.,0.,'right',.25,(0.,0.,1.),(1.,0.,0.))
    assert len(ax.images)==1
    assert ax.patches[-1].get_width()==block.VALVE_BAR_WIDTH
    assert ax.patches[-1].get_alpha()==.75
    fig.clear()


@pytest.mark.parametrize('family',['six_bar'])
def test_drawn_gas_displacement_matches_production_volume_orientation(tmp_path,family):
    from dada_solver.research.machine_render import _oriented
    from dada_solver.research.mechanism_view import MechanismModel
    study=load_study(initialize_kinematics(tmp_path/'study.toml',family,family))
    artifacts,_=candidate_mechanisms(study,{})
    angles=np.linspace(0.,2*math.pi,81)
    for side in ('small','large'):
        model=MechanismModel(artifacts[side])
        q=np.array([model.state(a)['position'] for a in angles])
        p=np.array([_oriented(model,a,side)[0]['P'][0] for a in angles])
        gas_displacement=-p if side=='small' else p
        normalized=(gas_displacement-gas_displacement.min())/np.ptp(gas_displacement)
        normalized_q=(q-q.min())/np.ptp(q)
        assert normalized==pytest.approx(normalized_q,abs=1e-12)


def test_webp_missing_encoder_is_explained(tmp_path,monkeypatch):
    from PIL import features
    study=load_study(initialize_kinematics(tmp_path/'study.toml','slider_crank','slider_crank'))
    monkeypatch.setattr(features,'check',lambda feature:False)
    monkeypatch.setattr('dada_solver.research.plot_data.replay_thermal',lambda *a,**k:dict(machine_cycle=cycle()))
    result=candidate_plots(study,dict(candidate_id='a'*64,physical={}),('mechanisms',))
    assert 'WebP support is unavailable' in result['mechanisms']['unavailable']


def test_shared_crank_is_resolved_without_old_per_side_phase(tmp_path):
    study=load_study(initialize_kinematics(tmp_path/'study.toml','four_bar','four_bar',coupling='shared_crank'))
    artifacts,_=candidate_mechanisms(study,{})
    assert artifacts['small'].scientific['geometry']['phase_rad']==artifacts['large'].scientific['geometry']['phase_rad']
    machine_geometry(artifacts,np.linspace(0.,2*math.pi,48,endpoint=False))


def test_encoding_failure_preserves_other_derived_curves(tmp_path,monkeypatch):
    study=load_study(initialize_kinematics(tmp_path/'study.toml','slider_crank','slider_crank'))
    pressure=dict(series=[],unit='Pa')
    monkeypatch.setattr('dada_solver.research.plot_data.replay_thermal',lambda *a,**k:dict(machine_cycle=cycle(),pressures=pressure))
    def fail(*args,**kwargs): raise ImportError('Pillow is unavailable')
    monkeypatch.setattr('dada_solver.research.machine_render.render_machine_webp',fail)
    result=candidate_plots(study,dict(candidate_id='a'*64,physical={}),('mechanisms','pressures'))
    assert result['pressures']==pressure
    assert 'Pillow' in result['mechanisms']['unavailable']


@pytest.mark.parametrize('family',['slider_crank','four_bar'])
def test_non_six_bar_pistons_stay_inside_the_drawn_cylinders(tmp_path,family):
    study=load_study(initialize_kinematics(tmp_path/'study.toml',family,family))
    artifacts,_=candidate_mechanisms(study,{})
    data=machine_geometry(artifacts,np.linspace(0.,2*math.pi,48,endpoint=False))
    for states,cylinder in data['sides'].values():
        low,high=sorted((cylinder['x_outer'],cylinder['x_inner']))
        assert all(low<=piston<=high for piston in cylinder['piston_positions'])


def test_six_bar_presentation_preserves_every_production_bar(tmp_path):
    from dada_solver.research.mechanism_view import MechanismModel
    study=load_study(initialize_kinematics(tmp_path/'study.toml','six_bar','six_bar'))
    artifacts,_=candidate_mechanisms(study,{})
    angles=np.linspace(0.,2*math.pi,48,endpoint=False)
    data=machine_geometry(artifacts,angles)
    for side in ('small','large'):
        model=MechanismModel(artifacts[side])
        for i,angle in enumerate(angles):
            source=model.state(angle)['joints'];drawing=data['sides'][side][0][i]
            crank=np.linalg.norm(drawing['B']-drawing['A'])
            for a,b in data['links'][side]:
                assert np.linalg.norm(drawing[a]-drawing[b])/crank==pytest.approx(np.linalg.norm(np.array(source[a])-source[b]),abs=1e-12)


def test_artifact_is_kept_when_candidate_geometry_is_unchanged(tmp_path):
    study=load_study(initialize_kinematics(tmp_path/'study.toml','six_bar','six_bar'))
    artifacts,_=candidate_mechanisms(study,{})
    for side in ('small','large'):
        assert artifacts[side] is study.artifacts[side]


@pytest.mark.parametrize('ho,hi,expected',[
    ('upstream','upstream','UU'),('downstream','downstream','DD'),
    ('upstream','downstream','UD'),('downstream','upstream','DU'),
])
def test_layout_is_read_from_placements_not_operation(ho,hi,expected):
    from types import SimpleNamespace
    for motor in (False,True):
        configuration=SimpleNamespace(heat_out_valve_placement=ho,
            heat_in_valve_placement=hi,motor_operation=motor)
        assert layout_for_configuration(configuration)==expected


@pytest.mark.parametrize('family',['slider_crank','four_bar'])
@pytest.mark.parametrize('small_direction,large_direction',[(True,True),(False,False),(True,False),(False,True)])
def test_both_volume_conventions_rods_and_outboard_linkages(tmp_path,family,small_direction,large_direction):
    from dada_solver.research.artifacts import MechanismArtifact
    from dada_solver.research.mechanism_view import MechanismModel
    study=load_study(initialize_kinematics(tmp_path/'study.toml',family,family))
    artifacts,_=candidate_mechanisms(study,{})
    for side,direction in [('small',small_direction),('large',large_direction)]:
        old=artifacts[side].scientific
        artifacts[side]=MechanismArtifact.create(family,dict(old['geometry'],volume_increases_with_coordinate=direction),settings=old['settings'])
    angles=np.linspace(0,2*math.pi,361,endpoint=False)
    data=machine_geometry(artifacts,angles)
    for side in ('small','large'):
        model=MechanismModel(artifacts[side]);outward=-1 if side=='small' else 1
        states,cyl=data['sides'][side]
        volume=np.asarray(model.law.value(angles))-1.
        gas=outward*(cyl['piston_positions']-cyl['x_inner'])
        assert (gas-gas.min())/np.ptp(gas)==pytest.approx((volume-volume.min())/np.ptp(volume),abs=1e-11)
        assert gas[np.argmin(volume)]<gas[np.argmax(volume)]
        low,high=sorted((cyl['x_inner'],cyl['x_outer']))
        assert np.all((cyl['piston_positions']>=low)&(cyl['piston_positions']<=high))
        separation=np.array([outward*(state['P'][0]-piston) for state,piston in zip(states,cyl['piston_positions'])])
        assert np.all(separation>=0)
        assert separation==pytest.approx(cyl['rod_length'],abs=1e-12)
        assert_clear_cylinder_walls(states,data['links'][side],cyl,side)
        for i,state in enumerate(states):
            assert state['P'][1]==pytest.approx(0.,abs=1e-12)
            if cyl['rod_length']>0.:
                assert all(outward*(p[0]-cyl['x_outer'])>=cyl['clearance']-1e-12 for p in state.values())
            source=model.state(angles[i])['joints']
            crank=np.linalg.norm(state['B']-state['A'])
            for a,b in data['links'][side]:
                assert np.linalg.norm(state[a]-state[b])/crank==pytest.approx(np.linalg.norm(np.array(source[a])-source[b]),abs=1e-11)


def test_cadence_does_not_change_cylinders_or_transfer_block_dimensions(tmp_path):
    study=load_study(initialize_kinematics(tmp_path/'study.toml','slider_crank','four_bar'))
    artifacts,_=candidate_mechanisms(study,{})
    before=machine_geometry(artifacts,np.linspace(0,2*math.pi,48,endpoint=False))
    after=machine_geometry(artifacts,np.linspace(0,2*math.pi,66,endpoint=False))
    assert before['heads']==after['heads']==dict(small=-2.16,large=2.16)
    for side in ('small','large'):
        for key in ('x_inner','x_outer','width'):
            assert before['sides'][side][1][key]==after['sides'][side][1][key]


def test_default_webp_has_exactly_66_frames_and_slow_cycle(tmp_path):
    from dada_solver.research.machine_render import FRAME_DURATION_MS,FRAME_COUNT
    study=load_study(initialize_kinematics(tmp_path/'study.toml','slider_crank','slider_crank'))
    artifacts,_=candidate_mechanisms(study,{})
    payload,metadata=render_machine_webp(artifacts,cycle(),layout='DD')
    assert FRAME_COUNT==66 and FRAME_DURATION_MS==85
    query,_=frame_series(cycle(),FRAME_COUNT)
    assert len(query)==66 and query[0]==0 and query[-1]<2*math.pi
    assert np.diff(query)==pytest.approx(2*math.pi/66)
    with Image.open(io.BytesIO(payload)) as image:
        assert image.n_frames==66 and image.info['loop']==0
        duration=0;top=image.height;bottom=0
        for i in range(image.n_frames):
            image.seek(i);image.load()
            assert image.info['duration']==85
            duration+=image.info['duration']
            pixels=np.asarray(image.convert('RGB')).astype(int)
            rows=np.flatnonzero(np.any(abs(pixels-240)>32,axis=(1,2)))
            top=min(top,int(rows[0]));bottom=max(bottom,int(rows[-1])+1)
        assert duration==5610
        assert 18<=top<=27 and 18<=image.height-bottom<=27
    assert metadata['cycle_duration_ms']==5610
    assert len(payload)<2_000_000


@pytest.mark.parametrize('family',['slider_crank','four_bar','six_bar'])
def test_exact_extrema_mid_volume_and_six_bar_rod_clearance(tmp_path,family):
    from dada_solver.research.mechanism_view import MechanismModel
    from dada_solver.research.artifacts import MechanismArtifact
    from scipy.optimize import brentq
    study=load_study(initialize_kinematics(tmp_path/'study.toml',family,family))
    artifacts,_=candidate_mechanisms(study,{})
    if family!='six_bar':
        # Opposite conventions on the two sides, each checked independently.
        old=artifacts['large'].scientific
        artifacts['large']=MechanismArtifact.create(family,dict(old['geometry'],volume_increases_with_coordinate=False),settings=old['settings'])
    checkpoints={};angles=[]
    for side in ('small','large'):
        model=MechanismModel(artifacts[side])
        extrema=model.law.model.stationary_points(model.law.side) if family=='four_bar' else model.geometry.stationary_points()
        minimum=next(p['angle_rad'] for p in extrema if p['kind']=='minimum')
        maximum=next(p['angle_rad'] for p in extrema if p['kind']=='maximum')
        middle=brentq(lambda t:model.law.value(t)-1.5,min(minimum,maximum),max(minimum,maximum))
        checkpoints[side]=[len(angles)+i for i in range(3)]
        angles.extend((minimum,middle,maximum))
    data=machine_geometry(artifacts,angles)
    for side,indices in checkpoints.items():
        states,cyl=data['sides'][side];outward=-1 if side=='small' else 1
        distances=outward*(cyl['piston_positions'][indices]-cyl['x_inner'])
        assert distances[0]==pytest.approx(.035*cyl['width'],abs=1e-9)
        assert distances[0]<distances[1]<distances[2]
        assert distances[1]==pytest.approx(.5*(distances[0]+distances[2]),abs=1e-9)
        assert np.array([outward*(s['P'][0]-p) for s,p in zip(states,cyl['piston_positions'])])==pytest.approx(cyl['rod_length'],abs=1e-12)
        if cyl['rod_length']>0.:
            assert all(outward*(point[0]-cyl['x_outer'])>=cyl['clearance']-1e-12 for s in states for point in s.values())
        assert_clear_cylinder_walls(states,data['links'][side],cyl,side)


def test_only_cylinder_outlines_move_by_their_line_thickness(tmp_path):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from dada_solver.research.machine_render import draw_machine_frame,CYLINDER_LINE_WIDTHS,_stroke_points
    study=load_study(initialize_kinematics(tmp_path/'study.toml','slider_crank','four_bar'))
    artifacts,_=candidate_mechanisms(study,{})
    query,data=frame_series(cycle(),66);geometry=machine_geometry(artifacts,query)
    fig=Figure(figsize=(12,5),dpi=100);FigureCanvasAgg(fig);ax=fig.add_axes((.01,.01,.98,.98))
    draw_machine_frame(ax,geometry,17,data,'DD',280,380);fig.canvas.draw()
    walls=[line for line in ax.lines if line.get_zorder()==6 and line.get_color()=='#111111']
    assert len(walls)==6
    for i,side in enumerate(('small','large')):
        frames,cyl=geometry['sides'][side];inward=1 if side=='small' else -1
        head=walls[3*i];upper=walls[3*i+1];lower=walls[3*i+2]
        original_x=ax.transData.transform((cyl['x_inner'],0))[0]
        new_x=ax.transData.transform((head.get_xdata()[0],0))[0]
        assert inward*(new_x-original_x)==pytest.approx(_stroke_points(ax,CYLINDER_LINE_WIDTHS[side])*fig.dpi/72.,abs=1e-10)
        assert head.get_xdata()[0]==upper.get_xdata()[0]==lower.get_xdata()[0]
        assert upper.get_xdata()[1]==lower.get_xdata()[1]==cyl['x_outer']
        assert upper.get_ydata()==pytest.approx([cyl['width']/2]*2)
        assert lower.get_ydata()==pytest.approx([-cyl['width']/2]*2)
        piston=next(line for line in ax.lines if line.get_linewidth()==_stroke_points(ax,9.5 if side=='small' else 10.))
        assert piston.get_xdata()==pytest.approx([cyl['piston_positions'][17]]*2)
    # Geometry baselines remain fixed; transfer endpoints now overlap the walls.
    coils=[patch for patch in ax.patches if isinstance(patch,PathPatch)]
    assert len(coils)==2
    assert geometry['heads']==dict(small=-2.16,large=2.16)
    fig.clear()


def test_external_webp_is_portable_and_does_not_mutate_report_data(tmp_path,monkeypatch):
    from tests.test_research_cockpit import campaign
    from dada_solver.research import report
    c,_=campaign(tmp_path/'campaign',1);data=report.inspect(c.directory)
    payload=b'RIFF-test-WebP-bytes'
    uri='data:image/webp;base64,'+base64.b64encode(payload).decode()
    data['selected']=data['records']
    data['selected'][0]['plots']=dict(mechanisms=dict(data_uri=uri,media_type='image/webp',layout='DD'))
    destination=tmp_path/'export dir'/'report.html'
    report.render_html(data,destination,external_webp=True)
    assets=list(destination.with_suffix('.assets').glob('*.webp'))
    assert len(assets)==1 and assets[0].read_bytes()==payload
    page=destination.read_text()
    assert uri not in page and 'report.assets/mechanisms-' in page
    assert 'motion.data_uri || motion.url' in page
    assert data['selected'][0]['plots']['mechanisms']['data_uri']==uri
    embedded=report.render_html(data,tmp_path/'embedded.html').read_text()
    assert uri in embedded
    # CLI forwards the opt-in independently from curve generation.
    from dada_solver.research.cli import main
    monkeypatch.setattr(report,'compare',lambda *a,**k:data)
    output=tmp_path/'cli.html'
    assert main(['report',str(c.directory),'--plots','mechanisms','--external-webp','--html',str(output)])==0
    assert list(output.with_suffix('.assets').glob('*.webp'))


@pytest.mark.parametrize('layout',['UU','DD'])
def test_transfer_rows_overlap_walls_and_have_requested_separation(tmp_path,layout):
    from dada_solver.research.machine_render import draw_machine_frame, _transfer_endpoints, _stroke_world, CYLINDER_LINE_WIDTHS
    study=load_study(initialize_kinematics(tmp_path/'study.toml','slider_crank','four_bar'))
    artifacts,_=candidate_mechanisms(study,{})
    query,data=frame_series(cycle(),66);geometry=machine_geometry(artifacts,query)
    fig=Figure(figsize=(12,5),dpi=100);ax=fig.add_axes((.01,.01,.98,.98))
    draw_machine_frame(ax,geometry,0,data,layout,280,380)
    ends=_transfer_endpoints(ax,geometry)
    width=geometry['sides']['large'][1]['width']
    # Every conduit lies above walls, below diodes, and above the serpentine.
    conduits=[image for image in ax.images if image.get_zorder() in (8,9)]
    assert conduits
    centers=[]
    for image in conduits:
        x0,x1,y0,y1=image.get_extent()
        a,b=image.get_transform().transform(((x0,(y0+y1)/2),(x1,(y0+y1)/2)))
        a,b=ax.transData.inverted().transform((a,b))
        assert ends['small']-1e-12<=a[0]<=b[0]<=ends['large']+1e-12
        centers.append(a[1])
    assert min(centers)==pytest.approx(-.155*1.25*width)
    assert max(centers)==pytest.approx(.155*1.25*width)
    coils=[patch for patch in ax.patches if isinstance(patch,PathPatch)]
    assert all(p.get_zorder()==7 for p in coils)
    walls=[line for line in ax.lines if line.get_zorder()==6 and line.get_color()=='#111111']
    for i,side in enumerate(('small','large')):
        sign=-1 if side=='small' else 1
        assert ends[side]-walls[3*i].get_xdata()[0]==pytest.approx(sign*.5*_stroke_world(CYLINDER_LINE_WIDTHS[side]))
    assert block.SEG_OUTER>.36 and block.SEG_MID<.44
    fig.clear()


def test_assembly_strokes_scale_with_viewport_without_moving_cylinders(tmp_path):
    from dada_solver.research.machine_render import draw_machine_frame,_cylinder_head,_stroke_points,REFERENCE_PIXELS_PER_UNIT
    study=load_study(initialize_kinematics(tmp_path/'study.toml','slider_crank','slider_crank'))
    artifacts,_=candidate_mechanisms(study,{})
    query,data=frame_series(cycle(),66);geometry=machine_geometry(artifacts,query)
    heads=[];strokes=[]
    for zoom in (1.,2.):
        expanded=dict(geometry)
        xmin,xmax,ymin,ymax=geometry['limits']
        expanded['limits']=(xmin*zoom,xmax*zoom,ymin*zoom,ymax*zoom)
        fig=Figure(figsize=(12,5),dpi=100);ax=fig.add_axes((.01,.01,.98,.98))
        draw_machine_frame(ax,expanded,0,data,'DD',280,380)
        heads.append(_cylinder_head(ax,geometry['sides']['small'][1],'small'))
        strokes.append(_stroke_points(ax,6.5))
        fig.clear()
    assert heads[0]==heads[1]
    assert strokes[1]==pytest.approx(strokes[0]/2.)
    # The calibrated reference retains the original point linewidths.
    fig=Figure(figsize=(12,5),dpi=100);ax=fig.add_axes((.01,.01,.98,.98))
    ax.set_xlim(0.,1176/REFERENCE_PIXELS_PER_UNIT)
    assert _stroke_points(ax,6.5)==pytest.approx(6.5)
    assert _stroke_points(ax,10.)==pytest.approx(10.)
    fig.clear()


@pytest.mark.parametrize('direction',['left','right'])
@pytest.mark.parametrize('opening',[0.,.5])
def test_diode_bar_has_original_thickness_and_reaches_conduit_triangle_intersection(direction,opening):
    fig=Figure();ax=fig.subplots()
    vertices,base,tip,_,_=block.valve_geometry(0.,0.,direction)
    block.draw_valve(ax,0.,0.,direction,opening,(0.,0.,1.),(1.,0.,0.))
    bar=ax.patches[-1]
    assert bar.get_width()==pytest.approx(block.BAR_W)
    overlap=(bar.get_x()+bar.get_width()-tip if direction=='left' else tip-bar.get_x())
    assert overlap==pytest.approx(block.VALVE_BAR_OVERLAP)
    # At the triangle-facing bar edge, the sloped triangle boundary meets
    # the conduit top/bottom exactly, for either mirrored diode direction.
    assert .5*block.TRI_H*overlap/block.TRI_W==pytest.approx(.5*block.PIPE_H)
    fig.clear()


@pytest.mark.parametrize('side',['small','large'])
@pytest.mark.parametrize('collision',[False,True])
def test_piston_rod_can_be_zero_only_without_wall_collision(side,collision):
    from dada_solver.research.machine_render import _place, _linkage_hits_cylinder, _stroke_world, CYLINDER_LINE_WIDTHS
    sign=-1 if side=='small' else 1
    # A narrow horizontal mechanism fits inside the open-ended cylinder;
    # moving B above a wall instead makes its connecting rod cross that wall.
    samples=[dict(A=np.array((sign*3.,0.)),B=np.array((sign*2.,2. if collision else .1)),P=np.array((sign*x,0.))) for x in (0.,1.)]
    links=[('A','B'),('B','P')]
    placed,cyl=_place(samples,side,0.,2.,.1,envelope=samples,visible_stroke=1.,links=links)
    assert (cyl['rod_length']>0)==collision
    for before,after,piston in zip(samples,placed,cyl['piston_positions']):
        assert sign*(after['P'][0]-piston)==pytest.approx(cyl['rod_length'])
        assert np.linalg.norm(after['B']-after['P'])==pytest.approx(np.linalg.norm(before['B']-before['P']))
    if not collision:
        assert all(state['P'][0]==p for state,p in zip(placed,cyl['piston_positions']))
    wall_head=cyl['x_inner']-sign*_stroke_world(CYLINDER_LINE_WIDTHS[side])
    assert not _linkage_hits_cylinder(placed,links,0.,wall_head,cyl['x_outer'],cyl['width'],.1)


def test_large_cylinder_minimum_size_for_large_linkage(tmp_path,monkeypatch):
    from dada_solver.research import machine_render
    study=load_study(initialize_kinematics(tmp_path/'study.toml','six_bar','six_bar'))
    artifacts,_=candidate_mechanisms(study,{})
    angles=np.linspace(0.,2*math.pi,66,endpoint=False)
    normal=machine_geometry(artifacts,angles)
    original=machine_render._oriented
    def tall_linkage(*args):
        joints,links=original(*args)
        # Stress the display envelope independently of the piston stroke.
        return {name:point*np.array((1.,50. if name!='P' else 1.)) for name,point in joints.items()},links
    monkeypatch.setattr(machine_render,'_oriented',tall_linkage)
    geometry=machine_geometry(artifacts,angles)
    assert geometry['heads']==normal['heads']
    assert geometry['sides']['large'][1]['width']>normal['sides']['large'][1]['width']
    assert geometry['sides']['small'][1]['width']>normal['sides']['small'][1]['width']
    assert geometry['sides']['small'][1]['width']/geometry['sides']['large'][1]['width']==pytest.approx(
        normal['sides']['small'][1]['width']/normal['sides']['large'][1]['width'])
    assert np.ptp(geometry['sides']['large'][1]['piston_positions'])==pytest.approx(
        np.ptp(normal['sides']['large'][1]['piston_positions']))
    payload,metadata=render_machine_webp(artifacts,cycle(),layout='DD',frames=3)
    assert metadata['large_cylinder_width_px']>=.20*metadata['height']
    with Image.open(io.BytesIO(payload)) as image:
        assert image.n_frames==3 and image.height==metadata['height']


@pytest.mark.parametrize('output,base,edges,color',[
    ('rocker',('C','D'),[('C','H'),('D','H')],'#2ca02c'),
    ('coupler',('B','C'),[('B','H'),('C','H')],'#ff7f0e'),
])
def test_four_bar_output_body_segments_share_color(tmp_path,output,base,edges,color):
    from dada_solver.research.machine_render import draw_machine_frame
    study=load_study(initialize_kinematics(tmp_path/'study.toml','four_bar','four_bar'))
    artifacts,_=candidate_mechanisms(study,{})
    query,data=frame_series(cycle(),66)
    geometry=machine_geometry(artifacts,query)
    # Exercise both rigid output-body connectivities on production joint states.
    for side in ('small','large'):
        geometry['links'][side]=[('A','B'),('B','C'),('C','D'),('H','P'),*edges]
    figure=Figure(figsize=(12,5),dpi=100);ax=figure.subplots()
    draw_machine_frame(ax,geometry,7,data,'DD',280,380)
    for side in ('small','large'):
        state=geometry['sides'][side][0][7]
        for a,b in [base,*edges]:
            matches=[line for line in ax.lines if line.get_linewidth()==2.55 and
                np.array_equal(line.get_xdata(),[state[a][0],state[b][0]]) and
                np.array_equal(line.get_ydata(),[state[a][1],state[b][1]])]
            assert matches and all(line.get_color()==color for line in matches)
