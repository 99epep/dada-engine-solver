"""Exact declaration edits, portable inputs and pure terminal controller actions."""
from dataclasses import asdict
import builtins
import copy
import json
from pathlib import Path
import tomllib
import numpy as np
import pytest
from dada_solver.research.study_editor import StudyEditor, initial_center
from dada_solver.research.schema import load_study, compile_study, candidate_for_values
from dada_solver.research.artifacts import MechanismLibrary
from dada_solver.research.mechanism_adaptation import paired_thermodynamic
from dada_solver.research.hardware_retuning import hardware_retuning
from dada_solver.research.presets import initialize_kinematics
from dada_solver.research.study_io import dumps
from dada_solver.research.cli import main
from tests.test_research_hardware_retuning import original_source
from tests.test_research_synthesis_architecture import artifact


@pytest.fixture(autouse=True)
def no_integration(monkeypatch):
    from dada_solver.campaign.evaluator import MachineEvaluator
    monkeypatch.setattr(MachineEvaluator,'evaluate',lambda *a,**k:pytest.fail('Editor must not evaluate'))
    from dada_solver.integration import CycleIntegrator
    monkeypatch.setattr(CycleIntegrator,'integrate_cycle',lambda *a,**k:pytest.fail('Editor must not integrate'))


def retuned(tmp_path,small='slider_crank',large='slider_crank'):
    source=original_source(tmp_path)
    library=MechanismLibrary((dict(family_id='pair',mechanisms={
        'small':artifact(small,'small').data,'large':artifact(large,'large').data},metadata={}),))
    paired=paired_thermodynamic(source,library,'pair',tmp_path/'paired.toml')
    return hardware_retuning(paired,tmp_path/'retuned/study.toml',scope='source-active')


def design(study):
    return compile_study(study).adapter.build(dict(study.fixed_parameters,**initial_center(study.data['parameters'])))


@pytest.mark.parametrize('small,large,count',[
    ('slider_crank','slider_crank',6),('four_bar','four_bar',22),('six_bar','six_bar',30),
    ('slider_crank','four_bar',14),('four_bar','six_bar',26),('six_bar','slider_crank',18),
])
def test_retune_release_mechanisms_keeps_initial_machine_and_hardware(tmp_path,small,large,count):
    source=retuned(tmp_path,small,large)
    editor=StudyEditor(source)
    before=editor.study
    snapshots={p:p.read_bytes() for p in source.parent.iterdir() if p.is_file()}
    inventory={p.name:p for p in editor.inventory()}
    n=f'kinematics.small.{next(iter(artifact(small,"small").scientific["geometry"]))}'
    assert inventory[n].source_of_value=='mechanism artifact'
    editor.release(groups=['mechanisms'])
    editor.configure_search('local',radius=.05)
    report=editor.review()
    assert len(report['activated_parameters'])==count
    output=editor.save(tmp_path/'combined/study.toml')
    after=load_study(output)
    assert len(after.space.parameters)==len(before.space.parameters)+count
    assert [p.name for p in after.space.parameters][:len(before.space.parameters)]==[p.name for p in before.space.parameters]
    assert initial_center(after.data['parameters'])==dict(initial_center(before.data['parameters']),
        **{n:before.fixed_parameters[n] for n in report['activated_parameters']})
    a,b=design(before),design(after)
    for field in ('configuration','heat_in','heat_out'): assert asdict(getattr(a,field))==asdict(getattr(b,field))
    angles=np.linspace(0.,2*np.pi,91)
    for side in ('small','large'):
        assert np.array_equal(getattr(a.kinematics,side).value(angles),getattr(b.kinematics,side).value(angles))
        assert before.artifacts[side].data==after.artifacts[side].data
    for key in ('objective','constraints','policies','numerical','execution'):
        assert before.data[key]==after.data[key]
    assert before.mechanical_constraints==after.mechanical_constraints
    assert after.data['search']['regions'][0]['center']==initial_center(after.data['parameters'])
    assert after.data['search']['radius_fraction']==.05
    for p,data in snapshots.items(): assert p.read_bytes()==data
    provenance=after.basis.data['provenance']
    assert provenance['paired_thermodynamic']==before.basis.data['provenance']['paired_thermodynamic']
    assert provenance['hardware_retuning']==before.basis.data['provenance']['hardware_retuning']
    assert provenance['study_edit']['source_study_id']==before.study_id
    assert after.study_id!=before.study_id
    assert main(['validate',str(output)])==0


def test_freeze_exact_initial_then_restore_saved_domain(tmp_path):
    source=original_source(tmp_path)
    editor=StudyEditor(source)
    name='operation.frequency_hz'; old=editor._rows()[name]
    editor.freeze(parameters=[name]);out=editor.save(tmp_path/'fixed/study.toml')
    fixed=load_study(out)
    assert name not in initial_center(fixed.data['parameters'])
    assert fixed.fixed_parameters[name]==old['initial']
    second=StudyEditor(out);row,origin=second.resolve_domain(name)
    assert row==old and 'recorded' in origin
    second.release(parameters=[name]);second.configure_search('global')
    active=load_study(second.save(tmp_path/'restored/study.toml'))
    assert next(p for p in active.space.parameters if p.name==name)==next(p for p in editor.study.space.parameters if p.name==name)


def test_freezing_changed_mechanical_initial_regenerates_artifact(tmp_path):
    source=retuned(tmp_path)
    editor=StudyEditor(source);name='kinematics.small.phase_rad'
    editor.release(parameters=[name])
    original=editor.study.artifacts['small']; value=editor.values()[name]+.01
    editor.edit_parameter(name,initial=value);editor.freeze(parameters=[name])
    assert editor.review()['artifact_changes']['small']['before']==original.content_hash
    generated=load_study(editor.save(tmp_path/'changed/study.toml'))
    new=generated.artifacts['small']
    assert new.scientific['geometry']['phase_rad']==value
    assert new.scientific['constraints']==original.scientific['constraints']
    assert new.content_hash!=original.content_hash
    assert generated.artifacts['large'].data==editor.study.artifacts['large'].data
    assert generated.data['kinematics']['small']['sha256']==new.content_hash
    assert new.data['provenance']['study_edit']['parent_artifact_hash']==original.content_hash


def test_active_override_does_not_relabel_reference_artifact(tmp_path):
    editor=StudyEditor(retuned(tmp_path));name='kinematics.small.phase_rad'
    editor.release(parameters=[name]);editor.edit_parameter(name,initial=editor.values()[name]+.01)
    new=load_study(editor.save(tmp_path/'active/study.toml'))
    assert new.artifacts['small'].data==editor.study.artifacts['small'].data
    assert initial_center(new.data['parameters'])[name]!=new.artifacts['small'].scientific['geometry']['phase_rad']


def test_reference_box_and_centered_synthesis_fallback_do_not_clamp(tmp_path):
    editor=StudyEditor(retuned(tmp_path));name='kinematics.small.phase_rad'
    row,origin=editor.resolve_domain(name)
    box=editor.basis['provenance']['paired_thermodynamic']['reference_boxes']['small']['phase_rad']
    assert [row['lower'],row['upper']]==box and 'paired' in origin
    editor.edit_parameter(name,value=box[1]+.2)
    row,origin=editor.resolve_domain(name)
    assert row['initial']==box[1]+.2 and 'centered synthesis' in origin
    assert row['lower']<row['initial']<row['upper']


def test_missing_domain_explicit_bounds_and_protected_categories(tmp_path):
    editor=StudyEditor(initialize_kinematics(tmp_path/'study.toml','harmonic','harmonic'))
    with pytest.raises(ValueError,match='operation.frequency_hz: no reliable'):
        editor.release(groups=['frequency'])
    value=editor.values()['operation.frequency_hz']
    editor.release(groups=['frequency'],domains={'operation.frequency_hz':dict(kind='continuous',lower=value/2,upper=value*2,transform='log')})
    editor.configure_search('local',radius=.1)
    p=next(p for p in load_study(editor.save(tmp_path/'edited/study.toml')).space.parameters if p.name=='operation.frequency_hz')
    assert p.transform=='log'
    bounds=editor.review()['effective_bounds'][p.name][0]
    assert bounds['lower']==pytest.approx(p.decode(p.encode(value)-.1))
    physical=StudyEditor(initialize_kinematics(tmp_path/'physical.toml','six_bar','six_bar'))
    with pytest.raises(ValueError,match='remain fixed'):
        physical.release(parameters=['kinematics.small.primary_branch'])
    with pytest.raises(ValueError,match='not editable'):
        physical.edit_parameter('kinematics.small.primary_branch',value=-1)


def test_integer_choice_and_source_domain_recovery(tmp_path):
    editor=StudyEditor(retuned(tmp_path))
    names=['microtube.heat_in.tube_count','valve.heat_in.placement']
    editor.freeze(parameters=names)
    editor.release(parameters=names)
    editor.configure_search('local',radius=.05)
    new=load_study(editor.save(tmp_path/'integer/study.toml'))
    from dada_solver.campaign.scheduled_search import ScheduledSobol
    scheduler=ScheduledSobol(new.space,new.data['search'])
    for _ in range(12):
        scheduler.next_point()
        assert type(scheduler.last_physical[names[0]]) is int
        assert scheduler.last_physical[names[1]] in ('upstream','downstream')
    assert new.data['search']['choice_scope']=='declared_choices'
    second=StudyEditor(hardware_retuning(tmp_path/'paired.toml',tmp_path/'limited/study.toml',groups=['frequency']))
    assert not next(p for p in second.inventory() if p.name==names[0]).active
    row,origin=second.resolve_domain(names[0]);assert row['kind']=='integer' and 'original hardware' in origin
    second.release(parameters=names);second.save(tmp_path/'expanded/study.toml')


def test_recenter_multiple_regions_requires_explicit_confirmation(tmp_path):
    editor=StudyEditor(retuned(tmp_path));raw=editor.study.data
    raw['search']['regions'].append(dict(raw['search']['regions'][0],id='second'))
    editor.source.write_text(dumps(raw));editor=StudyEditor(editor.source)
    preserved=editor.review()['search']['regions'];assert len(preserved)==2
    editor.configure_search('local',radius=.2)
    assert editor.review()['search']['regions']==preserved
    editor.release(groups=['mechanisms'])
    with pytest.raises(ValueError,match='multiple local regions'): editor.save(tmp_path/'bad/study.toml')
    assert not (tmp_path/'bad').exists()
    editor.configure_search('local',radius=.05,recenter=True)
    new=load_study(editor.save(tmp_path/'good/study.toml'))
    region,=new.data['search']['regions']
    expected=candidate_for_values(editor.definition,initial_center(editor.study.data['parameters']))
    assert region['source_candidate_id']==expected.candidate_id
    assert new.basis.data['provenance']['study_edit']['center_evidence'].startswith('unevaluated')


def test_all_fixed_and_execution_settings(tmp_path):
    editor=StudyEditor(original_source(tmp_path))
    editor.freeze(parameters=[p.name for p in editor.study.space.parameters])
    editor.configure_search('local',radius=.05)
    editor.edit_execution(budget='30m',max_candidates=512,seed=1234,scramble=False)
    report=editor.review();assert report['active_after']==0 and report['warnings']
    new=load_study(editor.save(tmp_path/'fixed/study.toml'))
    assert not new.space.parameters and new.data['search']['domain']=='fixed_global_bounds'
    assert new.data['execution']['default_budget']=='30m'
    assert new.data['search']['seed']==1234 and new.data['search']['scramble'] is False


def test_invalid_center_and_existing_destination_leave_no_partial_files(tmp_path):
    editor=StudyEditor(retuned(tmp_path))
    name='kinematics.small.rod_over_crank'
    editor.edit_parameter(name,value=.1)
    with pytest.raises((ValueError,RuntimeError)): editor.save(tmp_path/'invalid/study.toml')
    assert not (tmp_path/'invalid').exists()
    editor=StudyEditor(retuned_path:=tmp_path/'retuned/study.toml')
    with pytest.raises(ValueError,match='must be new'):editor.save(retuned_path)
    destination=tmp_path/'collision/study.toml';destination.parent.mkdir();destination.with_suffix('.basis.json').write_text('reserved')
    with pytest.raises(ValueError,match='must be new'):editor.save(destination)
    assert not destination.exists() and destination.with_suffix('.basis.json').read_text()=='reserved'


def test_write_failure_rolls_back_only_new_files(tmp_path,monkeypatch):
    editor=StudyEditor(retuned(tmp_path));output=tmp_path/'rollback/study.toml'
    original=Path.open
    def fail(path,*args,**kwargs):
        if path==output: raise OSError('simulated publication failure')
        return original(path,*args,**kwargs)
    monkeypatch.setattr(Path,'open',fail)
    with pytest.raises(OSError,match='publication failure'):editor.save(output)
    assert not list(output.parent.iterdir())


def test_batch_commands_match_engine(tmp_path):
    source=retuned(tmp_path,'slider_crank','four_bar')
    expected=StudyEditor(source);expected.release(groups=['mechanisms']);expected.configure_search('local',radius=.05)
    a=load_study(expected.save(tmp_path/'api/study.toml'))
    output=tmp_path/'cli/study.toml'
    assert main(['study','release',str(source),'--group','mechanisms','--radius','0.05','--output',str(output)])==0
    b=load_study(output)
    assert a.study_id==b.study_id
    assert main(['study','freeze',str(output),'--group','frequency','--output',str(tmp_path/'cli/fixed.toml')])==0
    assert 'operation.frequency_hz' not in initial_center(load_study(tmp_path/'cli/fixed.toml').data['parameters'])


def test_cli_explicit_domain_and_bad_input(tmp_path):
    source=initialize_kinematics(tmp_path/'study.toml','harmonic','harmonic')
    value=StudyEditor(source).values()['operation.frequency_hz']
    output=tmp_path/'out/study.toml'
    assert main(['study','release',str(source),'--group','frequency','--bounds',f'operation.frequency_hz={value/2}:{value*2}:log','--search','global','--output',str(output)])==0
    assert load_study(output).data['search']['domain']=='fixed_global_bounds'
    with pytest.raises(SystemExit) as error:main(['study','release',str(source),'--group','frequency','--output',str(tmp_path/'missing/study.toml')])
    assert error.value.code==2 and not (tmp_path/'missing').exists()


class SimulatedMenus:
    def __init__(self,actions,inputs=()): self.actions=iter(actions);self.inputs=iter(inputs);self.messages=[]
    def menu(self,*a,**k):return next(self.actions)
    def input(self,*a,**k):return next(self.inputs)
    def confirm(self,*a,**k):return True
    def message(self,text):self.messages.append(text)


def test_tui_controller_matches_batch_and_cancel_writes_nothing(tmp_path):
    from dada_solver.research.study_editor_tui import run_editor
    source=retuned(tmp_path);output=tmp_path/'interactive/study.toml'
    ui=SimulatedMenus([('toggle','kinematics-small'),('toggle','kinematics-large'),
        ('open','search'),('open','local'),('open','review'),('open','save')],['0.05',str(output)])
    assert run_editor(source,ui=ui)==output
    expected=StudyEditor(source);expected.release(groups=['mechanisms']);expected.configure_search('local',radius=.05)
    batch=load_study(expected.save(tmp_path/'batch/study.toml'))
    assert batch.study_id==load_study(output).study_id
    assert any('STUDY EDIT — REVIEW' in message for message in ui.messages)
    cancellation=SimulatedMenus([('toggle','kinematics-small'),('open','cancel')])
    before=set(source.parent.rglob('*'))
    assert run_editor(source,ui=cancellation) is None
    assert set(source.parent.rglob('*'))==before


def test_tui_non_terminal_and_optional_dependency_errors(tmp_path,monkeypatch):
    from dada_solver.research import study_editor_tui as tui
    monkeypatch.setattr(tui.sys.stdin,'isatty',lambda:False)
    with pytest.raises(ValueError,match='requires a terminal'):tui.TerminalMenus()
    monkeypatch.setattr(tui.sys.stdin,'isatty',lambda:True)
    monkeypatch.setattr(tui.sys.stdout,'isatty',lambda:True)
    original=builtins.__import__
    def missing(name,*a,**k):
        if name.startswith('prompt_toolkit'):raise ImportError('optional package absent')
        return original(name,*a,**k)
    monkeypatch.setattr(builtins,'__import__',missing)
    with pytest.raises(ValueError,match=r'dada-engine-solver\[tui\]'):tui.TerminalMenus()
    source=initialize_kinematics(tmp_path/'study.toml','slider_crank','slider_crank')
    StudyEditor(source).release(groups=['mechanisms'])  # pure engine needs no TUI
    with pytest.raises(SystemExit) as error:main(['study','--help'])
    assert error.value.code==0


def test_shared_crank_ownership_and_raw_spline_coordinates(tmp_path):
    source=initialize_kinematics(tmp_path/'shared.toml','four_bar','four_bar',coupling='shared_crank')
    editor=StudyEditor(source)
    names=editor.release(groups=['mechanisms'])
    assert len(names)==21
    assert names.count('kinematics.shared.phase_rad')==1
    assert 'kinematics.small.phase_rad' not in names
    editor.configure_search('local',radius=.05)
    new=load_study(editor.save(tmp_path/'combined/study.toml'))
    assert new.data['kinematics']['coupling']=='shared_crank'
    assert np.array_equal(design(new).kinematics.small.value(np.linspace(0,6,71)),
                          design(editor.study).kinematics.small.value(np.linspace(0,6,71)))
    spline=StudyEditor(initialize_kinematics(tmp_path/'spline.toml','free_spline','free_spline'))
    raw=[p for p in spline.inventory() if p.kind=='fixed_continuous']
    if raw:
        assert all(not p.releasable for p in raw)
        with pytest.raises(ValueError,match='remain fixed'): spline.release(parameters=[raw[0].name])


def test_user_domains_and_edits_reject_invalid_types_without_mutation(tmp_path):
    editor=StudyEditor(original_source(tmp_path))
    name='operation.frequency_hz'
    original=copy.deepcopy(editor.raw)
    with pytest.raises(ValueError): editor.edit_parameter(name,lower=-1,transform='log')
    assert editor.raw==original
    with pytest.raises(ValueError): editor.edit_parameter('microtube.heat_in.tube_count',initial=1.5)
    assert editor.raw==original
    editor.configure_search('global')
    assert editor.review()['search']==original['search']


def test_local_to_global_review_reports_removed_regions(tmp_path):
    editor=StudyEditor(retuned(tmp_path))
    editor.configure_search('global')
    report=editor.review()
    assert report['replaced_regions']
    assert any('removed' in w for w in report['warnings'])
    assert not any('one unevaluated' in w for w in report['warnings'])


def test_tui_declined_region_replacement_preserves_search(tmp_path):
    from dada_solver.research.study_editor_tui import _search
    editor=StudyEditor(retuned(tmp_path))
    editor.raw['search']['regions'].append(dict(editor.raw['search']['regions'][0],id='second'))
    before=copy.deepcopy(editor.review()['search'])
    ui=SimulatedMenus([('open','global')])
    ui.confirm=lambda *a,**k:False
    _search(editor,ui)
    assert editor.search_request is None
    assert editor.review()['search']==before


def test_tui_active_choice_menu_and_space_only_toggles_groups(tmp_path):
    from dada_solver.research.study_editor_tui import _parameter,run_editor
    source=retuned(tmp_path);editor=StudyEditor(source)
    name='valve.heat_in.placement'
    ui=SimulatedMenus([('open','value'),('open','upstream'),None])
    _parameter(editor,ui,name)
    assert editor.values()[name]=='upstream'
    cancel=SimulatedMenus([('toggle','save'),('open','cancel')])
    assert run_editor(source,ui=cancel) is None


def test_readable_review_includes_execution_changes_and_radius_scope(tmp_path):
    from dada_solver.research.study_editor_tui import review_text
    editor=StudyEditor(retuned(tmp_path))
    editor.release(groups=['mechanisms']);editor.configure_search('local',radius=.05)
    editor.edit_execution(budget='30m')
    review=editor.review();text=review_text(review)
    assert 'ALL active parameters' in text
    assert 'Newly released:' in text and 'Frozen:' in text
    assert 'Execution default_budget:' in text
    assert 'effective_bounds' not in text
    new=load_study(editor.save(tmp_path/'reviewed/study.toml'))
    assert new.basis.data['provenance']['study_edit']['execution_configuration_change']==review['execution_configuration_change']
