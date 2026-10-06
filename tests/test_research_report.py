"""Reports remain offline, escaped and read-only even with interrupted journals."""
import json
import gzip

import pytest

from dada_solver.campaign.runner import OptimizationCampaign
from dada_solver.research.report import inspect, compare, render_html, select_records
from tests.test_research_cli import definition
from tests.test_campaign import Clock, Evaluator


def test_report_reads_torn_tail_and_orphans_without_mutating(tmp_path,monkeypatch):
    d=definition(tmp_path); clock=Clock()
    run=tmp_path/'run'; c=OptimizationCampaign(d,run,evaluator=Evaluator(clock),clock=clock)
    c.run(100,maximum_candidates=2)
    import gzip
    history=c.history.path; records=c.history.load()
    (run/'recovery.json').write_text(json.dumps(records[1]))
    history.write_bytes(gzip.compress((json.dumps(records[0])+'\n').encode(),mtime=0)+b'\x1f\x8b')
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
    source=tmp_path/'study.toml'
    from dada_solver.research.study_io import dumps
    raw=d.study.data
    raw['objective']=dict(type='maximize_motor_power', unit='W')
    source.write_text(dumps(raw))
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
    raw=json.loads(path.read_text()); raw['scientific']['study']['purpose']='changed scientific purpose'
    path.write_text(json.dumps(raw))
    with pytest.raises(ValueError,match='definition does not match'): inspect(tmp_path/'run')


def test_exchanger_total_volume_and_comparison_fields(tmp_path):
    from tests.test_research_rescale import snapshot
    from dada_solver.research.presets import initialize_external_stream
    path=initialize_external_stream(tmp_path/'study.toml')
    _,definition,record=snapshot(path,tmp_path/'source')
    record['derived']={'heat_in_gas_volume_m3':.001,'heat_out_gas_volume_m3':.002}
    journal=tmp_path/'source'/'history.jsonl.gz'
    journal.write_bytes(gzip.compress((json.dumps(record)+'\n').encode()));before=journal.read_bytes()
    data=inspect(tmp_path/'source')
    assert data['records'][0]['derived']['total_exchanger_gas_volume_m3']==pytest.approx(.003)
    assert journal.read_bytes()==before
    from dada_solver.research.report import text_report
    assert 'exchanger gas volume=0.003 m³' in text_report(data,list_candidates=True)
    page=render_html(data,tmp_path/'report.html').read_text()
    assert 'Total exchanger gas volume [m³]' in page
    assert "['composite','COP × cooling power per microtube [W/microtube]'" in page
    assert "value=composite?'composite':productivity?'microtube'" in page
    comparison=page[page.index('function compare(){'):]
    assert comparison.index('for(const name of active)parameter(name)') < comparison.index("kinematic family")
    assert "cooling_cop_times_power_per_total_microtube_w:'W/microtube'" in comparison
    record['derived'].pop('heat_out_gas_volume_m3')
    journal.write_bytes(gzip.compress((json.dumps(record)+'\n').encode()))
    assert inspect(tmp_path/'source')['records'][0]['derived']['total_exchanger_gas_volume_m3'] is None
