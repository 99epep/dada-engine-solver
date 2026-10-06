"""Circular-cell envelopes, gas frusta and geometry-owned lossless diodes."""
from dataclasses import asdict,replace
import math
import json
import tomllib
import numpy as np
import pytest
from dada_solver.exchangers.microtube_geometry import MicrotubeBank
from dada_solver.exchangers.microtube import MicrotubeExchanger
from dada_solver.research.presets import initialize_kinematics,initialize_external_stream
from dada_solver.research.schema import load_study,compile_study
from dada_solver.research.study_io import dumps
from tests.test_exchanger_hardware import inputs


def bank(**changes):
    return MicrotubeBank(**dict(dict(tube_count=100,tube_length_m=.05,inner_diameter_m=.0005,
        wall_thickness_m=.00002,pitch_ratio=1.3,collector_half_angle_deg=30.,
        conduit_area_ratio=1.,additional_internal_volume_m3=1e-7),**changes))


@pytest.mark.parametrize('count',[1,7,100,1000])
def test_circular_geometry_and_gas_inventory(count):
    b=bank(tube_count=count);d=b.dimensions()
    pitch=1.3*(.0005+2*.00002)
    area=count*math.sqrt(3)/2*pitch**2
    db=math.sqrt(4*area/math.pi);tube_area=count*math.pi*.0005**2/4
    dc=math.sqrt(4*tube_area/math.pi);height=(db-dc)/(2*math.tan(math.pi/6))
    one=math.pi*height/12*(db*db+db*dc+dc*dc)
    assert b.effective_pitch_m==pitch
    assert d['bundle_face_area_m2']==pytest.approx(area)
    assert d['bundle_diameter_m']==pytest.approx(db)
    assert d['conduit_area_m2']==d['tube_flow_area_m2']==tube_area
    assert d['conduit_diameter_m']==pytest.approx(dc)
    assert d['collector_height_m']==pytest.approx(height)
    assert d['header_gas_volume_m3']==pytest.approx(2*one)
    assert d['working_gas_volume_m3']==pytest.approx(tube_area*.05+2*one+1e-7)
    assert d['fluid_envelope_length_m']==pytest.approx(.05+2*height)
    assert MicrotubeBank(**asdict(b)).dimensions()==d


def test_count_pitch_and_conduit_scaling():
    b=bank();d=b.dimensions()
    assert replace(b,tube_count=400).dimensions()['bundle_diameter_m']==pytest.approx(2*d['bundle_diameter_m'])
    assert replace(b,pitch_ratio=2.6).dimensions()['bundle_diameter_m']==pytest.approx(2*d['bundle_diameter_m'])
    q=replace(b,conduit_area_ratio=1.1).dimensions()
    assert q['conduit_area_m2']==pytest.approx(1.1*d['tube_flow_area_m2'])
    assert q['conduit_diameter_m']**2*math.pi/4==pytest.approx(q['conduit_area_m2'])


@pytest.mark.parametrize('change',[
    {'pitch_ratio':1},{'pitch_ratio':.9},{'pitch_ratio':float('nan')},
    {'conduit_area_ratio':.99},{'conduit_area_ratio':float('inf')},
    {'collector_half_angle_deg':0},{'collector_half_angle_deg':90},
    {'collector_half_angle_deg':float('nan')},{'conduit_area_ratio':10},
    {'header_depth_m':.001},{'pitch_m':.001},{'additional_internal_volume_m3':-1},
])
def test_invalid_geometry(change):
    with pytest.raises(ValueError): bank(**change)


@pytest.mark.parametrize('placement',['upstream','downstream'])
def test_diode_does_not_add_loss_or_restrict_area(ideal_gas,placement):
    _,props=inputs();b=bank()
    source=MicrotubeExchanger(b,props,1e-30,placement)
    assert source.outlet_valve_cda_m2==b.dimensions()['conduit_area_m2']
    components=source.build()
    for cda in (None,1e-20,1,1e10):
        changed=replace(source,outlet_valve_cda_m2=cda).build()
        for name in ('inlet','outlet'):
            link=getattr(changed,name)
            assert link.valve_cda_m2 is None
            assert link._flow_cap.effective_flow_area==b.tube_flow_area_m2
            assert link.directed_flow(200010,200000,350,ideal_gas)==getattr(components,name).directed_flow(200010,200000,350,ideal_gas)
            assert link.directed_flow(200000,200010,350,ideal_gas).mass_flow_rate==0
    a=components.inlet.directed_flow(200010,200000,350,ideal_gas)
    assert a==components.outlet.directed_flow(200010,200000,350,ideal_gas)
    assert replace(components.inlet,header_loss_coefficient=100).directed_flow(200010,200000,350,ideal_gas).mass_flow_rate<a.mass_flow_rate
    for diameter in (.00008,.0005,.0018):
        hx=replace(source,bank=replace(b,inner_diameter_m=diameter),outlet_valve_cda_m2=1e-30)
        assert hx.outlet_valve_cda_m2==hx.bank.tube_flow_area_m2
        assert hx.build().inlet._flow_cap.effective_flow_area==hx.bank.tube_flow_area_m2


def circular_study(path,version=3):
    path=(initialize_external_stream if version==3 else initialize_kinematics)(path)
    raw=tomllib.loads(path.read_text())
    raw['parameters']=[r for r in raw['parameters'] if not r['name'].endswith(('.pitch_m','.header_depth_m'))]
    raw['policies']['outlet_valve_cda']='geometry_conduit_area_v1'
    for side in ('heat_in','heat_out'):
        for field,value,unit in [('pitch_ratio',1.3,'1'),('collector_half_angle_deg',30.,'deg'),('conduit_area_ratio',1.,'1'),('additional_internal_volume_m3',1e-7,'m^3')]:
            raw['parameters'].append(dict(name=f'microtube.{side}.{field}',value=value,unit=unit))
    path.write_text(dumps(raw));return path


@pytest.mark.parametrize('version',[2,3])
def test_research_build_and_compiled_parity(tmp_path,version):
    path=circular_study(tmp_path/'study.toml',version)
    study=load_study(path);definition=compile_study(study)
    design=definition.adapter.build(definition.fixed_parameters)
    w=design.build()
    for side in ('heat_in','heat_out'):
        hx=getattr(design,side);metadata=dict(hx.build().metadata)
        assert metadata['valve_cda_m2']==metadata['conduit_area_m2']
        assert metadata['valve_model']=='ideal_diode_no_hydraulic_loss'
        assert metadata['header_gas_volume_m3']==2*metadata['single_collector_gas_volume_m3']
        assert hx.bank.header_depth_m is None and hx.bank.pitch_m is None
    from dada_solver.wall_backend import WallRHS,WallBackendSettings
    from tests.test_wall_backend import values
    pytest.importorskip('numba')
    rhs=WallRHS(w,WallBackendSettings('numba'));state=values(w);state[:2]*=1.00001
    np.testing.assert_allclose(rhs(0.,state),w.derivative(0.,state),rtol=2e-11,atol=1e-11)
    assert rhs.snapshot()['fallback_calls']==0


@pytest.mark.parametrize('field,value',[('pitch_m',.002),('header_depth_m',.001),('outlet_valve_cda_m2',.001)])
def test_derived_coordinates_not_research_variables(tmp_path,field,value):
    path=circular_study(tmp_path/'study.toml');raw=tomllib.loads(path.read_text())
    raw['parameters'].append(dict(name='microtube.heat_in.'+field,value=value,unit='m'))
    path.write_text(dumps(raw))
    with pytest.raises(ValueError,match='Unknown'): load_study(path)


def test_parameter_bounds_and_identity(tmp_path):
    path=circular_study(tmp_path/'study.toml');old=load_study(path).study_id
    raw=tomllib.loads(path.read_text());row=next(r for r in raw['parameters'] if r['name']=='microtube.heat_in.collector_half_angle_deg')
    row.pop('value');row.update(kind='continuous',initial=30.,lower=15.,upper=60.,transform='linear')
    path.write_text(dumps(raw));study=load_study(path)
    assert study.study_id!=old and len(study.space.parameters)==1
    row['upper']=90.;path.write_text(dumps(raw))
    with pytest.raises(ValueError,match='between 0 and 90'):load_study(path)


def test_external_stream_supplied_cda_is_replaced_and_reported(tmp_path):
    path=circular_study(tmp_path/'study.toml')
    d=compile_study(load_study(path));design=d.adapter.build(d.fixed_parameters)
    gas=design.configuration.gas
    hx=design.heat_in;reference=hx.build().outlet.directed_flow(200010,200000,350,gas)
    for value in (1e-30,1e10):
        changed=replace(hx,outlet_valve_cda_m2=value)
        assert changed.build().outlet.directed_flow(200010,200000,350,gas)==reference
        assert changed.outlet_valve_cda_m2==changed.bank.dimensions()['conduit_area_m2']
    metadata=dict(hx.build().metadata)
    assert 'independent_of_frustum_geometry' in metadata['header_loss_model']
    from dada_solver.research.report import HTML
    for key in ('bundle_diameter_m','conduit_area_m2','conduit_area_ratio','collector_height_m','valve_model'):
        assert key in HTML


def test_circular_rescale_keeps_additional_volume_extensive(tmp_path):
    from tests.test_research_rescale import snapshot
    from dada_solver.research.rescale import rescale
    path=circular_study(tmp_path/'study.toml')
    study,definition,record=snapshot(path,tmp_path/'source')
    source=definition.adapter.build(definition.fixed_parameters)
    scaled=rescale(tmp_path/'source',record['candidate_id'],4,tmp_path/'scaled.toml')
    d=compile_study(load_study(scaled));design=d.adapter.build(d.fixed_parameters)
    for side in ('heat_in','heat_out'):
        a=getattr(source,side);b=getattr(design,side)
        assert b.bank.tube_count==2*a.bank.tube_count
        assert b.bank.additional_internal_volume_m3==4*a.bank.additional_internal_volume_m3
        assert b.outlet_valve_cda_m2==pytest.approx(2*a.outlet_valve_cda_m2)
        assert b.bank.dimensions()['header_gas_volume_m3']==pytest.approx(2**1.5*a.bank.dimensions()['header_gas_volume_m3'])
    # Count grows by sqrt(s); circular header volume grows by s^0.75.
