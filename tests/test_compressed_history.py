"""Lossless member-wise journal appends, crash recovery and legacy reads."""
import gzip
import json
import pytest
from dada_solver.campaign.history import CampaignHistory,atomic_json,read_journal,journal_path
from dada_solver.campaign.candidate import canonical_json
from dada_solver.campaign.runner import OptimizationCampaign
from tests.test_research_cockpit import campaign
from tests.test_campaign import Clock,Evaluator
from tests.test_research_cli import definition


def test_gzip_roundtrip_standard_reader_and_resume(tmp_path):
    c,clock=campaign(tmp_path,3);before=c.history.path.read_bytes()
    assert c.history.path.name=='history.jsonl.gz'
    records=c.history.load()
    assert [json.loads(line) for line in gzip.decompress(before).splitlines()]==records
    ev=Evaluator(clock)
    OptimizationCampaign.resume(c.directory,evaluator=ev,clock=clock).run(100,maximum_candidates=2)
    assert c.history.path.read_bytes().startswith(before)
    after=c.history.load()
    assert after[:3]==records and len(after)==5
    assert set(ev.calls).isdisjoint(r['candidate_id'] for r in records)


@pytest.mark.parametrize('cut',[1,2,9,25,-1,-8])
def test_incomplete_final_member_recovery_is_lossless(tmp_path,cut):
    c,_=campaign(tmp_path,2);records=c.history.load()
    members=[gzip.compress((canonical_json(r)+'\n').encode(),compresslevel=6,mtime=0) for r in records]
    atomic_json(c.directory/'recovery.json',records[1])
    c.history.path.write_bytes(members[0]+members[1][:cut])
    before=c.history.path.read_bytes()
    read,_,offset=read_journal(c.history.path)
    assert read==records[:1] and offset==len(members[0])
    assert c.history.path.read_bytes()==before
    assert c.history.load()==records
    assert c.history.path.read_bytes()==b''.join(members)
    assert len(list(c.directory.glob('history_torn_tail_*.bin')))==1


@pytest.mark.parametrize('member_index',[0,1])
def test_crc_corruption_is_not_treated_as_torn_append(tmp_path,member_index):
    c,_=campaign(tmp_path,2);records=c.history.load()
    members=[bytearray(gzip.compress((canonical_json(r)+'\n').encode(),mtime=0)) for r in records]
    members[member_index][-8]^=1
    damaged=b''.join(members);c.history.path.write_bytes(damaged)
    with pytest.raises(ValueError,match='Corrupt compressed'): c.history.load()
    assert c.history.path.read_bytes()==damaged


def test_existing_plain_journal_remains_append_only(tmp_path):
    d=definition(tmp_path);directory=tmp_path/'legacy';directory.mkdir()
    clock=Clock();OptimizationCampaign(d,directory,evaluator=Evaluator(clock),clock=clock)
    path=directory/'history.jsonl';path.write_text('')
    c=OptimizationCampaign(d,directory,evaluator=Evaluator(clock),clock=clock)
    c.run(100,maximum_candidates=1);before=path.read_bytes()
    OptimizationCampaign.resume(directory,evaluator=Evaluator(clock),clock=clock).run(100,maximum_candidates=1)
    assert path.read_bytes().startswith(before)
    assert len(c.history.load())==2 and not (directory/'history.jsonl.gz').exists()
    (directory/'history.jsonl.gz').write_bytes(b'')
    with pytest.raises(ValueError,match='Both plain and compressed'): journal_path(directory)
