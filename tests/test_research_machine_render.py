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
        assert image.size==(1200,500) and image.info['loop']==0
        pixels=np.asarray(image.convert('RGB'))
        assert np.any((pixels[:,:,0]>180)&(pixels[:,:,2]<80))
    assert metadata['layout']==layout and metadata['cylinder_inner_faces']==dict(small=-2.16,large=2.16)
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
        assert image.n_frames==48
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
    assert ax.patches[-1].get_width()==block.BAR_W
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
        assert all(low<=state['P'][0]<=high for state in states)


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
