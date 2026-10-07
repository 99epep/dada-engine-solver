"""Frozen adapted mechanisms and portable, typed original hardware domains."""
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace
import json
import math
import tomllib

import numpy as np
import pytest

from dada_solver.research.hardware_retuning import hardware_retuning
from dada_solver.research.mechanism_adaptation import paired_thermodynamic
from dada_solver.research.presets import initialize_kinematics
from dada_solver.research.schema import load_study,compile_study,candidate_for_values
from dada_solver.research.study_io import dumps
from dada_solver.research.cli import main
from dada_solver.campaign.scheduled_search import ScheduledSobol
from tests.test_research_mechanism_adaptation import pair
from tests.test_research_rescale import snapshot


def original_source(tmp_path):
    path=initialize_kinematics(tmp_path/'original.toml','hybrid_compact','hybrid_compact')
    raw=tomllib.loads(path.read_text())
    for row in raw['parameters']:
        name=row['name']
        if name in ('volume.total_swept_m3','operation.frequency_hz','charge.total_mass_kg','microtube.heat_in.tube_length_m'):
            v=row.pop('value'); row.update(kind='continuous',initial=v,lower=v/2,upper=v*2,transform='log')
        elif name=='microtube.heat_in.tube_count':
            v=row.pop('value'); row.update(kind='integer',initial=v,lower=v//2,upper=v*2,encoding='nearest_even_v1')
    raw['parameters'].append(dict(name='valve.heat_in.placement',unit='1',kind='choice',
                                  initial='downstream',choices=['downstream','upstream']))
    # An active motion coordinate must never be restored as hardware.
    row=next(r for r in raw['parameters'] if r['name']=='kinematics.small.small_max_deg')
    v=row.pop('value'); row.update(kind='continuous',initial=v,lower=v-1,upper=v+1,transform='linear')
    path.write_text(dumps(raw))
    return path


@pytest.fixture
def forbidden_integration(monkeypatch):
    import dada_solver.campaign.evaluator as evaluator
    monkeypatch.setattr(evaluator.MachineEvaluator,'evaluate_with_control',lambda *a,**kw:pytest.fail('Generation must not integrate'))
    monkeypatch.setattr(evaluator,'solve_periodic_wall_machine',lambda *a,**kw:pytest.fail('Generation must not integrate'))


def paired_source(tmp_path,family):
    original=original_source(tmp_path)
    adapted=paired_thermodynamic(original,pair(family),'selected-pair',tmp_path/'paired.toml')
    study=load_study(adapted)
    values={p.name:p.initial for p in study.space.parameters}
    phase=next(n for n in values if n.startswith('kinematics.small.') and (n.endswith('phase_rad') or n.endswith('primary_phase')))
    values[phase]+=.01
    actual,definition,record=snapshot(adapted,tmp_path/'campaign',values)
    return original,actual,definition,record


@pytest.mark.parametrize('family',['slider_crank','four_bar','six_bar'])
def test_source_active_freezes_actual_adapted_mechanism(tmp_path,forbidden_integration,family):
    original,paired,definition,record=paired_source(tmp_path,family)
    output=tmp_path/'retuned/study.toml'
    assert main(['mechanism','retune',str(tmp_path/'campaign'),'--candidate',record['candidate_id'][:12],
                 '--scope','source-active','--output',str(output)])==0
    study=load_study(output)
    active={p.name for p in study.space.parameters}
    expected={p.name for p in load_study(original).space.parameters if not p.name.startswith('kinematics.')}
    assert active==expected
    assert [p.name for p in study.space.parameters]==[p.name for p in load_study(original).space.parameters if not p.name.startswith('kinematics.')]
    assert not any(n.startswith('kinematics.') for n in active)
    values={p.name:p.initial for p in study.space.parameters}
    rebuilt=compile_study(study).adapter.build(dict(study.fixed_parameters,**values))
    selected=definition.adapter.build(dict(paired.fixed_parameters,**record['physical']))
    assert asdict(rebuilt.configuration)==asdict(selected.configuration)
    assert asdict(rebuilt.heat_in)==asdict(selected.heat_in)
    assert asdict(rebuilt.heat_out)==asdict(selected.heat_out)
    angle=np.linspace(0,2*math.pi,721)
    for side in ('small','large'):
        assert np.array_equal(getattr(rebuilt.kinematics,side).value(angle),getattr(selected.kinematics,side).value(angle))
        a=study.artifacts[side]
        assert all(study.fixed_parameters[f'kinematics.{side}.{n}']==v for n,v in a.scientific['geometry'].items())
        assert study.data['kinematics'][side]['sha256']==a.content_hash
        if side=='small':
            assert a.content_hash!=paired.artifacts[side].content_hash
            assert a.data['provenance']['parent_artifact_hash']==paired.artifacts[side].content_hash
        else: assert a.data==paired.artifacts[side].data
    for key in ('objective','constraints','policies','numerical','execution'):
        assert study.data[key]==paired.data[key]
    assert study.mechanical_constraints==paired.mechanical_constraints
    assert main(['validate',str(output)])==0
    scheduler=ScheduledSobol(study.space,study.data['search']); scheduler.next_point()
    assert scheduler.last_physical==values
    for _ in range(32):
        scheduler.next_point()
        assert type(scheduler.last_physical['microtube.heat_in.tube_count']) is int
        assert scheduler.last_physical['charge.total_mass_kg']>0
    entry=study.basis.data['provenance']['hardware_retuning']
    assert entry['source']['candidate_id']==record['candidate_id']
    assert entry['paired_source']['study_id']==paired.study_id
    assert entry['abstract_source']['study_id']==load_study(original).study_id
    assert entry['selected_parameters']==list(values)
    assert entry['mechanisms']=={s:dict(hash=a.content_hash,parent_artifact_hash=paired.artifacts[s].content_hash) for s,a in study.artifacts.items()}
    assert load_study(output).study_id==study.study_id
    original.unlink()  # Generation and later passes must not depend on this path.
    second=hardware_retuning(output,tmp_path/'second.toml',groups=('volumes',))
    assert {p.name for p in load_study(second).space.parameters}=={'volume.total_swept_m3'}


@pytest.mark.parametrize('groups,expected',[
    (('exchangers',),{'microtube.heat_in.tube_count','microtube.heat_in.tube_length_m'}),
    (('volumes',),{'volume.total_swept_m3'}),
    (('frequency',),{'operation.frequency_hz'}),
    (('charge',),{'charge.total_mass_kg'}),
    (('valves',),{'valve.heat_in.placement'}),
    (('exchangers','volumes'),{'microtube.heat_in.tube_count','microtube.heat_in.tube_length_m','volume.total_swept_m3'})])
def test_targeted_groups_only_filter_originally_active(tmp_path,forbidden_integration,groups,expected):
    _,paired,_,_=paired_source(tmp_path,'slider_crank')
    output=hardware_retuning(paired.path,tmp_path/'retune.toml',groups=groups)
    study=load_study(output)
    assert {p.name for p in study.space.parameters}==expected
    assert all(not n.startswith('kinematics.') for n in expected)


def test_effective_domains_keep_types_and_clip_edges(tmp_path,forbidden_integration):
    original,paired,definition,_=paired_source(tmp_path,'slider_crank')
    domain=load_study(original)
    raw=paired.data
    for row in raw['parameters']:
        if row['name']=='operation.frequency_hz':
            row['value']=next(p.upper for p in domain.space.parameters if p.name==row['name'])
    # A retuning descendant's current values, not the abstract source, are authoritative.
    paired.path.write_text(dumps(raw))
    output=hardware_retuning(paired.path,tmp_path/'retune.toml',radius=.1)
    study=load_study(output); entry=study.basis.data['provenance']['hardware_retuning']
    for p in study.space.parameters:
        original_parameter=next(x for x in domain.space.parameters if x.name==p.name)
        assert type(p) is type(original_parameter)
        if p.name.startswith('valve.'):
            assert p.choices==original_parameter.choices
            continue
        assert p.lower==original_parameter.lower and p.upper==original_parameter.upper
        b=entry['effective_bounds'][p.name]
        assert p.lower<=b['lower']<=p.initial<=b['upper']<=p.upper
        if p.name.endswith('tube_count'): assert type(b['lower']) is type(b['upper']) is int
        if p.name=='operation.frequency_hz':
            assert b['normalized_center']==1. and b['normalized_upper']==1.
    sched=ScheduledSobol(study.space,study.data['search']); seen=set()
    for _ in range(64):
        sched.next_point(); seen.add(sched.last_physical['valve.heat_in.placement'])
    assert seen=={'downstream','upstream'}


def test_progressive_retuning_and_evaluation_source_preserve_lineage(tmp_path,forbidden_integration):
    _,paired,definition,record=paired_source(tmp_path,'six_bar')
    evaluation=tmp_path/'evaluation.json'
    evaluation.write_text(json.dumps(dict(artifact_type='research_evaluation_v1',name='Paired source',
        definition=dict(definition.identity,definition_id=definition.definition_id),record=record)))
    first=hardware_retuning(evaluation,tmp_path/'hx.toml',groups=('exchangers',))
    study=load_study(first)
    _,d,r=snapshot(first,tmp_path/'hx-campaign',{'microtube.heat_in.tube_count':4000})
    second=hardware_retuning(tmp_path/'hx-campaign',tmp_path/'hx-volumes.toml',candidate=r['candidate_id'],groups=('exchangers','volumes'))
    new=load_study(second); entry=new.basis.data['provenance']['hardware_retuning']
    assert entry['source']['candidate_id']==r['candidate_id']
    assert entry['paired_source']['candidate_id']==record['candidate_id']
    assert len(entry['ancestors'])==1
    assert next(p.initial for p in new.space.parameters if p.name.endswith('tube_count'))==4000
    assert {s:a.content_hash for s,a in new.artifacts.items()}=={s:a.content_hash for s,a in study.artifacts.items()}


def test_domains_without_snapshot_require_verified_explicit_source(tmp_path,forbidden_integration):
    original,paired,_,_=paired_source(tmp_path,'slider_crank')
    data=paired.basis.data; data['provenance']['paired_thermodynamic'].pop('hardware_source_domain')
    raw=paired.data; import hashlib
    text=json.dumps(data,indent=2)+'\n'; paired.path.with_suffix('.basis.json').write_text(text)
    raw['sources']['machine']['sha256']=hashlib.sha256(text.encode()).hexdigest(); paired.path.write_text(dumps(raw))
    with pytest.raises(ValueError,match='--source-study'):
        hardware_retuning(paired.path,tmp_path/'missing.toml')
    output=hardware_retuning(paired.path,tmp_path/'restored.toml',source_study=original)
    assert load_study(output).space.parameters
    with pytest.raises(ValueError,match='original paired source'):
        hardware_retuning(paired.path,tmp_path/'wrong.toml',source_study=output)


def test_invalid_selections_and_non_paired_source_are_rejected(tmp_path,forbidden_integration):
    original,paired,_,_=paired_source(tmp_path,'slider_crank')
    for kwargs in (dict(groups=('external-stream',)),dict(parameters=('kinematics.small.phase_rad',)),
                   dict(parameters=('volume.swept_ratio',)),dict(scope='source-active',groups=('volumes',)),
                   dict(radius=0),dict(radius=math.nan)):
        with pytest.raises(ValueError): hardware_retuning(paired.path,tmp_path/'invalid.toml',**kwargs)
    assert not (tmp_path/'invalid.toml').exists()
    with pytest.raises(ValueError,match='paired_thermodynamic'):
        hardware_retuning(original,tmp_path/'invalid.toml')


def test_current_external_warm_start_rescales_only_initial_guess(tmp_path,monkeypatch):
    _,paired,_,_=paired_source(tmp_path,'slider_crank')
    output=hardware_retuning(paired.path,tmp_path/'charge.toml',groups=('charge',))
    study=load_study(output); definition=compile_study(study)
    assert study.data['warm_start']['initial_source']=='source_exact'
    old=np.asarray(definition.initial_wall_state['values'])
    captures=[]
    def solve(wrapper,state,**kw):
        captures.append(state.copy())
        return SimpleNamespace(status='interrupted',message='test initial guess only',history=(),
            last_complete_state=None,converged=False,backend_statistics={})
    monkeypatch.setattr('dada_solver.campaign.evaluator.solve_periodic_wall_machine',solve)
    from dada_solver.campaign.evaluator import MachineEvaluator
    p=study.space.parameters[0]
    MachineEvaluator(definition).evaluate(candidate_for_values(definition,{p.name:p.initial}))
    assert np.array_equal(captures[-1][:8],old[:8])
    MachineEvaluator(definition).evaluate(candidate_for_values(definition,{p.name:p.initial*1.05}))
    assert captures[-1][:8:2].sum()==pytest.approx(p.initial*1.05)
    assert np.allclose(captures[-1][:8],old[:8]*(p.initial*1.05/old[:8:2].sum()),rtol=1e-14)


def test_external_stream_group_and_explicit_parameter_union(tmp_path,forbidden_integration):
    from dada_solver.research.presets import initialize_external_stream
    path=initialize_external_stream(tmp_path/'original.toml',mode='refrigeration')
    raw=tomllib.loads(path.read_text())
    names={'external_stream.heat_in.mass_flow_kg_s','operation.frequency_hz'}
    for row in raw['parameters']:
        if row['name'] in names:
            value=row.pop('value'); row.update(kind='continuous',initial=value,lower=value/2,upper=value*2,transform='log')
    path.write_text(dumps(raw))
    adapted=paired_thermodynamic(path,pair('slider_crank'),'selected-pair',tmp_path/'paired.toml')
    output=hardware_retuning(adapted,tmp_path/'stream.toml',groups=('external-stream',),parameters=('operation.frequency_hz',))
    study=load_study(output)
    assert {p.name for p in study.space.parameters}==names
    assert compile_study(study).configuration.motor_operation is False
    assert study.data['objective']==load_study(adapted).data['objective']


def test_changed_inventory_without_candidate_state_uses_standard_uniform_fallback(tmp_path,forbidden_integration):
    _,paired,_,_=paired_source(tmp_path,'slider_crank')
    first=hardware_retuning(paired.path,tmp_path/'charge.toml',groups=('charge',))
    initial=next(p.initial for p in load_study(first).space.parameters if p.name=='charge.total_mass_kg')
    _,_,record=snapshot(first,tmp_path/'charge-campaign',{'charge.total_mass_kg':initial*1.05})
    output=hardware_retuning(tmp_path/'charge-campaign',tmp_path/'next.toml',candidate=record['candidate_id'],groups=('charge',))
    study=load_study(output)
    assert study.data['warm_start']['initial_source']=='uniform'
    assert next(p.initial for p in study.space.parameters)==initial*1.05
    assert 'different inventory' in study.basis.data['provenance']['hardware_retuning']['warm_start']['reason']


def test_synthesis_protocol_executes_retuning_with_no_mechanism_release(tmp_path,forbidden_integration):
    _,paired,_,_=paired_source(tmp_path,'six_bar')
    from tests.test_research_synthesis_architecture import target_for
    from dada_solver.research.synthesis import SynthesisPlan,SynthesisRequest,release_coordinates
    target=target_for('six_bar')
    plan=SynthesisPlan(target,SynthesisRequest(target.content_hash,'six_bar','design_exploitation',('hardware_retuning',)))
    output=plan.execute(source=paired.path,output=tmp_path/'retune.toml',groups=('frequency',))
    assert release_coordinates('hardware_retuning',('small','large'),family='six_bar')==()
    assert {p.name for p in load_study(output).space.parameters}=={'operation.frequency_hz'}
