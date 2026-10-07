"""Paired thermodynamic study generation preserves the exact source machine."""
from dataclasses import asdict
import gzip
import json
import math
import tomllib

import numpy as np
import pytest

from dada_solver.research.artifacts import MechanismArtifact, MechanismLibrary
from dada_solver.research.mechanism_adaptation import paired_thermodynamic
from dada_solver.research.presets import initialize_kinematics
from dada_solver.research.schema import load_study, compile_study
from dada_solver.research.study_io import dumps
from dada_solver.research.cli import main
from tests.test_research_synthesis_architecture import artifact
from tests.test_research_rescale import snapshot


def pair(family, *, constraints=()):
    mechanisms = {}
    for side in ('small','large'):
        original = artifact(family,side)
        mechanisms[side] = MechanismArtifact.create(family,original.scientific['geometry'],
            settings=original.scientific['settings'], constraints=constraints,
            provenance=dict(stage='synthetic_pair', side=side)).data
    return MechanismLibrary((dict(family_id='selected-pair',mechanisms=mechanisms,
                                 metadata=dict(provenance=dict(stage='opposite_local_adaptation'))),))


@pytest.fixture
def source(tmp_path,monkeypatch):
    from dada_solver.campaign.evaluator import MachineEvaluator
    monkeypatch.setattr(MachineEvaluator,'evaluate_with_control',lambda *a,**kw:pytest.fail('Generation must not integrate'))
    path = initialize_kinematics(tmp_path/'source.toml','hybrid_compact','hybrid_compact')
    raw = tomllib.loads(path.read_text())
    row = next(r for r in raw['parameters'] if r['name']=='operation.frequency_hz')
    value = row.pop('value')
    row.update(kind='continuous',initial=value,lower=value/2,upper=value*2,transform='linear')
    path.write_text(dumps(raw))
    study,definition,record=snapshot(path,tmp_path/'campaign',{'operation.frequency_hz':value*1.25})
    return path,study,definition,record


@pytest.mark.parametrize('family,count',[('slider_crank',6),('four_bar',22),('six_bar',30)])
def test_exact_pair_portable_study_and_local_search(source,tmp_path,family,count):
    path,original,definition,record=source
    library=pair(family)
    lib_path=tmp_path/'pair.json'; library.save(lib_path)
    output=tmp_path/'adapted/study.toml'
    assert main(['mechanism','adapt',str(tmp_path/'campaign'),'--candidate',record['candidate_id'][:12],
                 '--library',str(lib_path),'--family-id','selected-pair','--output',str(output)])==0
    study=load_study(output)
    values={p.name:p.initial for p in study.space.parameters}
    expected=definition.adapter.build(dict(original.fixed_parameters,**record['physical']))
    actual=compile_study(study).adapter.build(dict(study.fixed_parameters,**values))
    assert len(values)==count
    assert all(n.startswith('kinematics.') and 'branch' not in n for n in values)
    assert asdict(actual.configuration)==asdict(expected.configuration)
    assert asdict(actual.heat_in)==asdict(expected.heat_in)
    assert asdict(actual.heat_out)==asdict(expected.heat_out)
    for name,value in dict(original.fixed_parameters,**record['physical']).items():
        if not name.startswith('kinematics.'):
            assert study.fixed_parameters[name]==value
    for key in ('objective','constraints','policies','numerical','execution'):
        assert study.data[key]==original.data[key]
    for side,raw in library.member('selected-pair')['mechanisms'].items():
        assert study.artifacts[side].data==raw
        assert study.data['kinematics'][side]['sha256']==raw['content_hash']
        for name,value in raw['scientific']['geometry'].items():
            assert dict(study.fixed_parameters,**values)[f'kinematics.{side}.{name}']==value
        from dada_solver.research.families import build_side
        selected_law,_=build_side(raw['scientific']['settings'],raw['scientific']['geometry'],side,
                                  getattr(actual.configuration.machine_volumes,side+'_cylinder'))
        angles=np.linspace(0,2*math.pi,721)
        assert np.array_equal(getattr(actual.kinematics,side).value(angles),selected_law.value(angles))
    assert main(['validate',str(output)])==0
    assert study.data['search']['domain']=='local_regions_v1'
    assert study.data['search']['regions'][0]['center']==values
    from dada_solver.campaign.scheduled_search import ScheduledSobol
    # Standard center-first search uses exact values, not decoded approximations.
    search=ScheduledSobol(study.space,study.data['search'])
    search.next_point()
    assert search.last_physical==values
    provenance=study.basis.data['provenance']['paired_thermodynamic']
    assert provenance['source']['candidate_id']==record['candidate_id']
    assert provenance['artifacts']=={s:a.content_hash for s,a in study.artifacts.items()}
    assert provenance['source_assessment']['status']=='feasible'
    assert all(r['satisfied'] for r in actual.kinematics.mechanical_diagnostics)
    with pytest.raises(ValueError,match='must not exist'):
        paired_thermodynamic(path,library,'selected-pair',output)


def test_constraint_conjunction_preserves_artifacts(source,tmp_path):
    path,original,_,_=source
    raw=original.data
    raw['mechanical_constraints'].append(dict(side='small',metric='maximum_absolute_first_derivative',
                                             relation='maximum',limit=100.,unit='m^3/rad'))
    path.write_text(dumps(raw))
    library=pair('slider_crank',constraints=(dict(metric='maximum_absolute_first_derivative',
                                               relation='maximum',limit=200.,unit='m^3/rad'),))
    output=paired_thermodynamic(path,library,'selected-pair',tmp_path/'adapt.toml')
    study=load_study(output)
    assert next(r for r in study.mechanical_constraints if r['side']=='small' and r['metric']=='maximum_absolute_first_derivative')['limit']==100.
    assert study.artifacts['small'].data==library.member('selected-pair')['mechanisms']['small']
    changed=dict(study.fixed_parameters,**{p.name:p.initial for p in study.space.parameters})
    changed['kinematics.small.rod_over_crank']=.5
    from dada_solver.campaign.adapters import PreflightRejection
    with pytest.raises(PreflightRejection,match=''):
        compile_study(study).adapter.build(changed)


def test_incomplete_pair_and_constraints_fail_without_writing(source,tmp_path):
    path,*_=source
    library=pair('slider_crank')
    member=library.member('selected-pair'); member['mechanisms'].pop('small')
    with pytest.raises(ValueError,match='complete'):
        paired_thermodynamic(path,MechanismLibrary((member,)),'selected-pair',tmp_path/'absent.toml')
    raw=tomllib.loads(path.read_text())
    raw['mechanical_constraints'].append(dict(side='small',metric='minimum_primary_transmission_sine',
                                             relation='minimum',limit=.1,unit='1'))
    # The source family itself must support a declared requirement.
    raw['kinematics']['small']=dict(family='four_bar',**{k:v for k,v in artifact('four_bar').scientific['settings'].items() if k!='family'})
    raw['parameters']=[r for r in raw['parameters'] if not r['name'].startswith('kinematics.small.')]
    from dada_solver.research.families import parameter_specs
    a=artifact('four_bar'); specs=parameter_specs(a.scientific['settings'],'small')
    raw['parameters'] += [dict(name='kinematics.small.'+n,value=v,unit=specs[n].unit) for n,v in a.scientific['geometry'].items()]
    path.write_text(dumps(raw))
    with pytest.raises(ValueError,match='unavailable'):
        paired_thermodynamic(path,library,'selected-pair',tmp_path/'absent.toml')
    assert not (tmp_path/'absent.toml').exists()
    assert not (tmp_path/'absent.basis.json').exists()


def test_reference_pressure_freezes_inventory(source,tmp_path):
    path,*_=source
    from tests.test_research_reference_charge import derived
    raw=tomllib.loads(path.read_text()); derived(raw); path.write_text(dumps(raw))
    original=load_study(path)
    design=compile_study(original).adapter.build(dict(original.fixed_parameters,**{p.name:p.initial for p in original.space.parameters}))
    output=paired_thermodynamic(path,pair('six_bar'),'selected-pair',tmp_path/'adapt.toml')
    study=load_study(output)
    assert study.fixed_parameters['charge.total_mass_kg']==design.configuration.charge.total_mass
    assert study.data['policies']['charge']=='explicit_inventory'
    assert study.basis.data['provenance']['paired_thermodynamic']['inventory_policy_change'] is not None


def test_warm_state_reuses_standard_source_exact_policy(source,tmp_path):
    path,study,definition,record=source
    warm=study.basis.data['warm_start']
    assert warm is not None
    from dada_solver.campaign.evaluator import _state_record
    record['final_periodic_state']=_state_record(warm['values'],
        'four_gas_volumes_m_U_plus_H_i_H_o_wall_energy','microtube_wall_10_state','motor',
        periodic=True,wall_capacities=warm['wall_capacities_j_k'])
    record['converged']=True
    (tmp_path/'campaign/history.jsonl.gz').write_bytes(gzip.compress((json.dumps(record)+'\n').encode()))
    output=paired_thermodynamic(tmp_path/'campaign',pair('slider_crank'),'selected-pair',tmp_path/'adapt.toml',candidate=record['candidate_id'])
    new=load_study(output)
    assert new.data['warm_start']['initial_source']=='source_exact'
    assert new.basis.data['warm_start']['values']==warm['values']
    assert compile_study(new).initial_wall_state['source_candidate_id']==record['candidate_id']
    assert new.basis.data['provenance']['paired_thermodynamic']['warm_start']['selected_candidate_state']


def test_standalone_evaluation_and_cooling_objective_are_supported(tmp_path,monkeypatch):
    from dada_solver.research.presets import initialize_external_stream
    from dada_solver.campaign.evaluator import MachineEvaluator
    monkeypatch.setattr(MachineEvaluator,'evaluate_with_control',lambda *a,**kw:pytest.fail('Must not integrate'))
    path=initialize_external_stream(tmp_path/'source.toml',mode='refrigeration')
    original,definition,record=snapshot(path,tmp_path/'campaign')
    evaluation=tmp_path/'evaluation.json'
    evaluation.write_text(json.dumps(dict(artifact_type='research_evaluation_v1',name='Synthetic cooling evaluation',
        definition=dict(definition.identity,definition_id=definition.definition_id),record=record)))
    output=paired_thermodynamic(evaluation,pair('six_bar'),'selected-pair',tmp_path/'adapt.toml')
    study=load_study(output)
    assert study.data['objective']==original.data['objective']
    assert study.data['constraints']==original.data['constraints']
    assert compile_study(study).configuration.motor_operation is False
    values=dict(study.fixed_parameters,**{p.name:p.initial for p in study.space.parameters})
    expected=definition.adapter.build(dict(original.fixed_parameters,**record['physical']))
    actual=compile_study(study).adapter.build(values)
    assert asdict(actual.configuration)==asdict(expected.configuration)
    assert asdict(actual.heat_in)==asdict(expected.heat_in)
    assert asdict(actual.heat_out)==asdict(expected.heat_out)


def test_production_topology_roots_reject_candidate_before_integration(source,tmp_path,monkeypatch):
    path,*_=source
    output=paired_thermodynamic(path,pair('four_bar'),'selected-pair',tmp_path/'adapt.toml')
    study=load_study(output)
    from dada_solver.four_bar import SharedCrankFourBarVolumeKinematics
    monkeypatch.setattr(SharedCrankFourBarVolumeKinematics,'stationary_points',
        lambda *a,**kw:[dict(kind='maximum'),dict(kind='minimum')]*2)
    from dada_solver.campaign.adapters import PreflightRejection
    with pytest.raises(PreflightRejection,match='before integration') as error:
        compile_study(study).adapter.build(dict(study.fixed_parameters,**{p.name:p.initial for p in study.space.parameters}))
    assert any(r['value']==4 and not r['satisfied'] for r in error.value.diagnostics)


def test_periodic_angle_bounds_remain_centered_across_seam(source,tmp_path):
    path,*_=source
    library=pair('slider_crank')
    member=library.member('selected-pair')
    for side,data in member['mechanisms'].items():
        g=data['scientific']['geometry']; g['phase_rad']=2*math.pi+.01
        member['mechanisms'][side]=MechanismArtifact.create('slider_crank',g).data
    output=paired_thermodynamic(path,MechanismLibrary((member,)),'selected-pair',tmp_path/'adapt.toml',radius=.02)
    study=load_study(output)
    for p in study.space.parameters:
        if p.name.endswith('.phase_rad'):
            assert p.initial==2*math.pi+.01
            assert (p.lower+p.upper)/2==p.initial
            assert p.upper-p.lower==pytest.approx(2*math.pi)


def test_protocol_executes_paired_stage_with_explicit_source(source,tmp_path):
    path,study,*_=source
    from tests.test_research_synthesis_architecture import target_for
    from dada_solver.research.synthesis import SynthesisPlan,SynthesisRequest
    target=target_for('slider_crank')
    plan=SynthesisPlan(target,SynthesisRequest(target.content_hash,'slider_crank','design_exploitation',
                                              ('paired_thermodynamic',),retained_family_ids=('selected-pair',),
        mechanical_constraints=(dict(metric='minimum_rod_axis_cosine',relation='minimum',limit=.1,unit='1'),)))
    output=plan.execute(source=path,library=pair('slider_crank'),output=tmp_path/'adapt.toml')
    loaded=load_study(output)
    assert len(loaded.space.parameters)==6
    assert sum(r['metric']=='minimum_rod_axis_cosine' and r['limit']==.1 for r in loaded.mechanical_constraints)==2


@pytest.mark.parametrize('radius',[0.,math.nan,.51,True])
def test_invalid_local_radius(radius,tmp_path):
    with pytest.raises(ValueError,match='radius'):
        paired_thermodynamic('unread',pair('slider_crank'),'selected-pair',tmp_path/'adapt.toml',radius=radius)
