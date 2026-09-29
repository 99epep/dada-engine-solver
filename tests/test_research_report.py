"""Reports remain offline, escaped and read-only even with interrupted journals."""
import json

import pytest

from dada_solver.campaign.runner import OptimizationCampaign
from dada_solver.research.report import inspect, compare, render_html, select_records
from tests.test_research_cli import definition
from tests.test_campaign import Clock, Evaluator


def test_report_reads_torn_tail_and_orphans_without_mutating(tmp_path,monkeypatch):
    d=definition(tmp_path); clock=Clock()
    run=tmp_path/'run'; c=OptimizationCampaign(d,run,evaluator=Evaluator(clock),clock=clock)
    c.run(100,maximum_candidates=2)
    history=run/'history.jsonl'; lines=history.read_bytes().splitlines(keepends=True)
    (run/'recovery.json').write_bytes(lines[1])
    history.write_bytes(lines[0]+b'{"unfinished":')
    before={str(p.relative_to(run)):p.read_bytes() for p in run.rglob('*') if p.is_file()}
    monkeypatch.setattr('dada_solver.campaign.evaluator.MachineEvaluator.evaluate',lambda *a:pytest.fail('evaluation called'))
    monkeypatch.setattr('dada_solver.campaign.runner.OptimizationCampaign.run',lambda *a,**k:pytest.fail('search called'))
    data=compare([run]); assert len(data['records'])==2
    assert any('Incomplete' in w for w in data['warnings'])
    assert any('awaiting journal' in w for w in data['warnings'])
    page=render_html(data,tmp_path/'report.html').read_text()
    assert 'useful_mechanical_power_w' in page
    assert 'https://' not in page and '<script src=' not in page
    assert before=={str(p.relative_to(run)):p.read_bytes() for p in run.rglob('*') if p.is_file()}


def test_report_escapes_labels_and_preserves_missing_values(tmp_path):
    d=definition(tmp_path); clock=Clock()
    c=OptimizationCampaign(d,tmp_path/'run',evaluator=Evaluator(clock),clock=clock)
    c.run(100,maximum_candidates=2)
    data=inspect(tmp_path/'run'); data['name']='</script><script>alert(1)</script>'
    data['scientific']['study']['purpose']=data['name']
    page=render_html(data,tmp_path/'report.html').read_text()
    assert '<script>alert(1)</script>' not in page
    assert '\\u003c/script' in page
    selected=select_records(data,[data['records'][0]['candidate_id'][:12],data['records'][1]['candidate_id']])
    assert len(selected)==2
    with pytest.raises(ValueError,match='absent or ambiguous'): select_records(data,['missing'])


def test_changed_study_comparison_has_no_combined_ranking(tmp_path):
    d=definition(tmp_path); clock=Clock()
    a=OptimizationCampaign(d,tmp_path/'a',evaluator=Evaluator(clock),clock=clock); a.run(100,maximum_candidates=1)
    source=tmp_path/'study.toml'; source.write_text(source.read_text().replace('required_power = 25.0','required_power = 30.0'))
    from dada_solver.research.schema import compile_study,load_study
    b=OptimizationCampaign(compile_study(load_study(source)),tmp_path/'b',evaluator=Evaluator(clock),clock=clock); b.run(100,maximum_candidates=1)
    data=compare([tmp_path/'a',tmp_path/'b'])
    assert not data['comparison_compatible'] and not data['best']
    assert len(data['compared_sources'])==2


def test_inspection_rejects_corrupted_scientific_identity(tmp_path):
    d=definition(tmp_path); clock=Clock()
    c=OptimizationCampaign(d,tmp_path/'run',evaluator=Evaluator(clock),clock=clock)
    c.run(100,maximum_candidates=1)
    path=tmp_path/'run'/'definition.json'
    raw=json.loads(path.read_text()); raw['scientific']['fixed']['hot_air_inlet_K']=598.15
    path.write_text(json.dumps(raw))
    with pytest.raises(ValueError,match='definition does not match'): inspect(tmp_path/'run')
