"""Parity frozen from the historical compact law before its extraction."""
import json
from pathlib import Path
import tomllib
import numpy as np
import pytest
from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.hybrid_compact_kinematics import HybridCompactKinematics
from dada_solver.research.families import FAMILIES,build_side,available_metrics
from dada_solver.research.presets import initialize_kinematics,initialize_external_stream
from dada_solver.research.schema import load_study,compile_study,candidate_for_values
from dada_solver.research.study_io import dumps

FIXTURES=Path(__file__).parent/'fixtures/hybrid_compact'
REFERENCE=json.loads((FIXTURES/'reference.json').read_text())


def model(parameters):
    return HybridCompactKinematics(*(CylinderVolumeLimits(**REFERENCE['limits'][s+'_cylinder'])
        for s in ('small','large')),**parameters)


@pytest.mark.parametrize('case',[0,1])
def test_dense_historical_scalar_and_vector_parity(case):
    motion=model(REFERENCE['cases'][case])
    with np.load(FIXTURES/'trajectories.npz') as data:
        actual=np.array([motion.cylinder_volumes_and_derivatives(float(t)) for t in data['theta']])
        # Pure extraction retains the actual arithmetic, not just a curve fit.
        np.testing.assert_array_equal(actual,data[f'case_{case}'])
        vector=np.asarray(motion.cylinder_volumes_and_derivatives(data['theta'])).T
        np.testing.assert_array_equal(vector,data[f'vector_{case}'])
        np.testing.assert_allclose(vector,actual,rtol=2e-11,atol=1e-13)
        for side,index in (('small',0),('large',1)):
            parameters={k:v for k,v in REFERENCE['cases'][case].items() if k.startswith(side+'_')}
            law,_=build_side(dict(family='hybrid_compact'),parameters,side,getattr(motion,side+'_limits'))
            selected=np.array([[law.value(float(t)),law.value(float(t),1)] for t in data['theta']])
            np.testing.assert_array_equal(selected,data[f'case_{case}'][:,[index,index+2]])


def test_phase_volume_direction_and_derivative_convention():
    p=REFERENCE['cases'][0];m=model(p)
    for side,maximum,duration in (('small',p['small_max_deg'],p['small_down_duration_deg']),('large',0.,p['large_down_duration_deg'])):
        value=getattr(m,side+'_cylinder_volume');derivative=getattr(m,side+'_cylinder_volume_derivative')
        limits=getattr(m,side+'_limits')
        assert value(-np.deg2rad(maximum))==pytest.approx(limits.maximum,abs=1e-14)
        assert value(-np.deg2rad(maximum+duration))==pytest.approx(limits.minimum,abs=1e-14)
        for theta in (.1,1.3,2.7,5.4):
            assert derivative(theta)==pytest.approx((value(theta+1e-6)-value(theta-1e-6))/2e-6,rel=2e-7,abs=1e-11)
            assert value(theta)==pytest.approx(value(theta+2*np.pi),abs=1e-14)
    assert m.small_physical_stroke is None and m.large_physical_stroke is None
    assert 'maximum_absolute_second_derivative' not in available_metrics('hybrid_compact')


@pytest.mark.parametrize('key,value',[
    ('small_down_duration_deg',34.),('large_down_duration_deg',326.),
    ('large_down_rounding',0.),('small_up_rounding',.49),
    ('small_down_kink_u',.01),('small_down_kink_q',.99),
    ('large_up_kink_u',.99),('large_up_kink_q',.01)])
def test_historical_domain_guards(key,value):
    with pytest.raises(ValueError): model(dict(REFERENCE['cases'][0],**{key:value}))


def test_all_nine_coordinates_owned_independently_and_active(tmp_path):
    assert len(FAMILIES)==11 and FAMILIES[-1]=='hybrid_compact'
    path=initialize_kinematics(tmp_path/'study.toml','hybrid_compact','hybrid_compact')
    raw=tomllib.loads(path.read_text())
    for row in raw['parameters']:
        if not row['name'].startswith('kinematics.'): continue
        value=row.pop('value');row.update(initial=value,lower=.99*value,upper=1.01*value,kind='continuous',transform='linear')
    path.write_text(dumps(raw));d=compile_study(load_study(path))
    assert len(d.space.parameters)==9
    assert sum(p.name.startswith('kinematics.small.') for p in d.space.parameters)==5
    assert sum(p.name.startswith('kinematics.large.') for p in d.space.parameters)==4
    values={p.name:p.initial for p in d.space.parameters}
    first=candidate_for_values(d,values)
    key='kinematics.small.small_max_deg';values[key]+=.1
    assert candidate_for_values(d,values).candidate_id!=first.candidate_id
    d.adapter.build(dict(d.fixed_parameters,**values)).build()


@pytest.mark.parametrize('small,large',[('hybrid_compact','harmonic'),('slider_crank','hybrid_compact')])
def test_mixed_family_and_v3_schema(tmp_path,small,large):
    path=initialize_external_stream(tmp_path/'study.toml',small,large)
    study=load_study(path)
    assert study.settings['small']['family']==small
    assert study.settings['large']['family']==large
