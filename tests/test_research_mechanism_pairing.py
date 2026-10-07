"""Human-selected homogeneous and mixed pairs preserve independent artifacts."""
import json
import pytest
from dada_solver.research.artifacts import MechanismArtifact, MechanismLibrary
from dada_solver.research.mechanism_pairing import pair_mechanisms
from dada_solver.research.cli import main
from dada_solver.research.schema import load_study
from dada_solver.research.mechanism_adaptation import paired_thermodynamic
from dada_solver.research.hardware_retuning import hardware_retuning
from tests.test_research_synthesis_architecture import artifact,target_for
from tests.test_research_hardware_retuning import original_source


def library(path,family,side):
    raw=artifact(family,side).data
    member=dict(family_id=f'{family}-{side}',mechanisms={side:raw},metadata={})
    result=MechanismLibrary((member,));result.save(path)
    return result


@pytest.mark.parametrize('small,large,count',[
    ('slider_crank','slider_crank',6),('four_bar','four_bar',22),('six_bar','six_bar',30),
    ('slider_crank','four_bar',14),('four_bar','six_bar',26),('six_bar','slider_crank',18)])
def test_cli_pair_adapt_and_retune(tmp_path,monkeypatch,small,large,count):
    from dada_solver.campaign.evaluator import MachineEvaluator
    monkeypatch.setattr(MachineEvaluator,'evaluate_with_control',lambda *a,**k:pytest.fail('No integration during pairing/adaptation/retuning'))
    sl=tmp_path/'small.json';ll=tmp_path/'large.json';output=tmp_path/'pair.json'
    sm=library(sl,small,'small');lm=library(ll,large,'large')
    if small==large:
        common=tmp_path/'common.json';MechanismLibrary(sm.members+lm.members).save(common)
        args=[str(common)];sources={'small':str(common),'large':str(common)}
    else:
        args=['--small-library',str(sl),'--large-library',str(ll)];sources={'small':str(sl),'large':str(ll)}
    assert main(['mechanism','pair',*args,'--small',sm.members[0]['family_id'],'--large',lm.members[0]['family_id'],'--output',str(output)])==0
    paired=MechanismLibrary.load(output);member=paired.members[0]
    assert member['family_id'].startswith('pair/')
    for side,original in [('small',sm),('large',lm)]:
        assert member['mechanisms'][side]==original.members[0]['mechanisms'][side]
        p=member['metadata']['provenance'][side]
        assert p['source_library']==sources[side]
        assert p['artifact_hash']==member['mechanisms'][side]['content_hash']
        assert p['physical_family']==member['metadata']['families'][side]
    assert member['metadata']['provenance']['operation']=='manual_pair_selection'
    repeat=pair_mechanisms(small_library=sources['small'],small=sm.members[0]['family_id'],large_library=sources['large'],large=lm.members[0]['family_id'],output=tmp_path/'repeat.json')
    assert repeat.to_data()==paired.to_data()
    source=original_source(tmp_path)
    adapted=paired_thermodynamic(source,paired,member['family_id'],tmp_path/'thermo.toml')
    study=load_study(adapted)
    assert len(study.space.parameters)==count
    for side,family,n in [('small',small,{'slider_crank':3,'four_bar':11,'six_bar':15}[small]),('large',large,{'slider_crank':3,'four_bar':11,'six_bar':15}[large])]:
        assert study.settings[side]['family']==family
        from dada_solver.research.synthesis import synthesis_protocol
        owned={f'kinematics.{side}.{coordinate}' for coordinate in synthesis_protocol(family).coordinates('paired_thermodynamic')}
        assert {p.name for p in study.space.parameters if p.name.startswith(f'kinematics.{side}.')}==owned
        assert len(owned)==n
        for p in study.space.parameters:
            if p.name in owned:assert p.initial==member['mechanisms'][side]['scientific']['geometry'][p.name.split('.')[-1]]
        assert study.artifacts[side].data==member['mechanisms'][side]
    retuned=load_study(hardware_retuning(adapted,tmp_path/'retuned.toml',groups=('frequency',)))
    assert not any(p.name.startswith('kinematics.') for p in retuned.space.parameters)
    for side in ('small','large'):assert retuned.artifacts[side].data==member['mechanisms'][side]
    assert main(['validate',str(retuned.path)])==0


@pytest.mark.parametrize('args',[
    ['common.json','--small-library','small.json'],
    ['common.json','--large-library','large.json'],
    ['--small-library','small.json'],['--large-library','large.json'],[],
    ['common.json','--small','one','--small','different'],
    ['--small-library','one.json','--small-library','other.json','--large-library','large.json'],
])
def test_ambiguous_cli_forms_fail_before_reading_files(tmp_path,args,capsys):
    with pytest.raises(SystemExit):main(['mechanism','pair',*args,'--small','s','--large','l','--output',str(tmp_path/'pair.json')])
    assert not (tmp_path/'pair.json').exists()
    assert capsys.readouterr().err


def test_missing_side_id_intermediate_and_invalid_inputs(tmp_path):
    sl=tmp_path/'small.json';ll=tmp_path/'large.json'
    library(sl,'slider_crank','small');library(ll,'four_bar','large')
    def pair(s='slider_crank-small',l='four_bar-large'):
        return pair_mechanisms(small_library=sl,small=s,large_library=ll,large=l,output=tmp_path/'pair.json')
    with pytest.raises(ValueError,match='Unknown'):pair(s='absent')
    # Extracting SMALL from a LARGE-only member is invalid.
    with pytest.raises(ValueError,match='no SMALL'):pair_mechanisms(small_library=ll,small='four_bar-large',large_library=ll,large='four_bar-large',output=tmp_path/'pair.json')
    from dada_solver.research.families import PRIMARY_COORDINATES
    full=artifact('six_bar','small').scientific['geometry']
    primary=MechanismArtifact.create('six_bar',{n:full[n] for n in (*PRIMARY_COORDINATES,'primary_branch')},settings=dict(component='primary'))
    MechanismLibrary((dict(family_id='primary',mechanisms={'small':primary.data},metadata={}),)).save(tmp_path/'primary.json')
    with pytest.raises(ValueError,match='complete'):pair_mechanisms(small_library=tmp_path/'primary.json',small='primary',large_library=ll,large='four_bar-large',output=tmp_path/'pair.json')
    for invalid in ('{}','[]','not JSON'):
        sl.write_text(invalid)
        with pytest.raises(ValueError,match='Invalid SMALL mechanism library'):pair()
    with pytest.raises(ValueError,match='physical'):MechanismArtifact.create('free_spline',{})
    sl.unlink();lib=library(sl,'slider_crank','small');bad=lib.to_data()
    bad['members'][0]['mechanisms']['small']['content_hash']='invalid'
    from dada_solver.campaign.candidate import content_hash
    bad['content_hash']=content_hash({k:v for k,v in bad.items() if k!='content_hash'})
    sl.write_text(json.dumps(bad))
    with pytest.raises(ValueError,match='hash'):pair()
    assert not (tmp_path/'pair.json').exists()


def test_extract_both_sides_from_complete_members(tmp_path):
    from tests.test_research_mechanism_adaptation import pair
    a=pair('slider_crank').members[0];b=pair('four_bar').members[0];a['family_id']='a';b['family_id']='b'
    path=tmp_path/'complete.json';MechanismLibrary((a,b)).save(path)
    result=pair_mechanisms(small_library=path,small='a',large_library=path,large='b',output=tmp_path/'pair.json')
    assert result.members[0]['mechanisms']=={'small':a['mechanisms']['small'],'large':b['mechanisms']['large']}


def test_catalogue_selectors_complete_and_primary(tmp_path,monkeypatch):
    from dada_solver.research.synthesis import synthesis_member
    from dada_solver.research.synthesis_catalogue import render_synthesis_catalogue
    from tests.test_research_mechanism_adaptation import pair
    target=target_for('slider_crank');artifacts={s:MechanismArtifact.from_data(a) for s,a in pair('slider_crank').members[0]['mechanisms'].items()}
    lib=MechanismLibrary((synthesis_member('quoted family',artifacts,target),))
    monkeypatch.setattr('dada_solver.research.synthesis_catalogue._preview',lambda *a:'<svg/>')
    path=render_synthesis_catalogue(lib,target,tmp_path/'custom.html',library_path=tmp_path/'source library.json')
    text=path.read_text()
    for token in ('Select as SMALL','Select as LARGE','Copy SMALL selector','Copy LARGE selector','--small-library','--large-library','Complete pair','mechanism adapt','pair-command','source library.json'):
        assert token in text
    assert 'dada-research mechanism pair' in text


def test_catalogue_selection_javascript_builds_safe_short_command(tmp_path,monkeypatch):
    import re,shlex,shutil,subprocess
    from html.parser import HTMLParser
    from dada_solver.research.synthesis import synthesis_member
    from dada_solver.research.synthesis_catalogue import render_synthesis_catalogue
    node=shutil.which('node')
    if node is None:pytest.skip('Node is optional for executing standalone catalogue controls')
    target=target_for('slider_crank');family_id="family with ' quote"
    artifacts={s:artifact('slider_crank',s) for s in ('small','large')}
    library_path=tmp_path/"source with ' quote.json"
    lib=MechanismLibrary((synthesis_member(family_id,artifacts,target),))
    monkeypatch.setattr('dada_solver.research.synthesis_catalogue._preview',lambda *a:'<svg/>')
    text=render_synthesis_catalogue(lib,target,tmp_path/'catalogue.html',library_path=library_path).read_text()
    class Buttons(HTMLParser):
        def __init__(self):super().__init__();self.buttons={};self.base=None
        def handle_starttag(self,tag,attrs):
            attrs=dict(attrs)
            if tag=='pre' and attrs.get('id')=='pair-command':self.base=attrs['data-base']
            if tag=='button' and 'data-side' in attrs:
                self.buttons[attrs['data-side']]={k.removeprefix('data-'):v for k,v in attrs.items() if k.startswith('data-')}
    parsed=Buttons();parsed.feed(text)
    script=re.search(r'<script>(.*?)</script>',text,re.S).group(1)
    runner="""const fs=require('fs'),vm=require('vm');const input=JSON.parse(fs.readFileSync(0,'utf8'));
const box={dataset:{base:input.base},textContent:''};
const context={document:{getElementById:()=>box},navigator:{}};vm.createContext(context);
vm.runInContext(input.script,context);
context.selectMechanism({dataset:input.buttons.small});context.selectMechanism({dataset:input.buttons.large});
process.stdout.write(box.textContent);"""
    result=subprocess.run([node,'-e',runner],input=json.dumps(dict(base=parsed.base,buttons=parsed.buttons,script=script)),text=True,capture_output=True)
    assert result.returncode==0,result.stderr
    assert shlex.split(result.stdout)==['dada-research','mechanism','pair',str(library_path),'--small',family_id,'--large',family_id,'--output','pair.json']
