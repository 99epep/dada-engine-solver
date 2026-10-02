"""Research CLI and the unchanged campaign durability contracts."""
import json
from pathlib import Path

import pytest

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.campaign.runner import OptimizationCampaign
from dada_solver.research.cli import initialize, main
from dada_solver.research.schema import load_study, compile_study
from dada_solver.research.report import inspect
from tests.test_campaign import Clock, Evaluator


def definition(tmp_path):
    return compile_study(load_study(initialize(tmp_path/'study.toml')))


def test_split_resume_matches_uninterrupted_and_uses_snapshots(tmp_path):
    d = definition(tmp_path); clock = Clock()
    whole = OptimizationCampaign(d,tmp_path/'whole',evaluator=Evaluator(clock),clock=clock)
    whole.run(100,maximum_candidates=4)
    split = OptimizationCampaign(d,tmp_path/'split',evaluator=Evaluator(clock),clock=clock)
    split.run(100,maximum_candidates=2)
    (tmp_path/'study.toml').unlink(); (tmp_path/'study.basis.json').unlink()
    resumed = OptimizationCampaign.resume(tmp_path/'split',evaluator=Evaluator(clock),clock=clock)
    resumed.run(100,maximum_candidates=2)
    a,b = whole.history.load(),resumed.history.load()
    assert [r['candidate_id'] for r in a] == [r['candidate_id'] for r in b]
    assert [r['sequence_index'] for r in b] == [0,1,2,3]
    assert (tmp_path/'split'/'basis.json').exists()


def test_crash_retries_the_pending_candidate(tmp_path):
    d = definition(tmp_path); clock = Clock(); observed=[]
    class Crash:
        def evaluate(self,candidate):
            observed.append(candidate.candidate_id)
            raise RuntimeError('synthetic crash')
    campaign=OptimizationCampaign(d,tmp_path/'run',evaluator=Crash(),clock=clock)
    with pytest.raises(RuntimeError,match='synthetic'): campaign.run(100,maximum_candidates=1)
    evaluator=Evaluator(clock)
    resumed=OptimizationCampaign.resume(tmp_path/'run',evaluator=evaluator,clock=clock)
    resumed.run(100,maximum_candidates=1)
    assert evaluator.calls==observed
    assert resumed.history.state()['pending'] is None


def test_exact_duplicate_cache_is_retained(tmp_path,monkeypatch):
    d=definition(tmp_path); clock=Clock(); evaluator=Evaluator(clock)
    from dada_solver.campaign.strategy import SobolStrategy
    original=SobolStrategy.next_point
    def repeated(self):
        original(self)
        return tuple([.5]*5)
    monkeypatch.setattr(SobolStrategy,'next_point',repeated)
    campaign=OptimizationCampaign(d,tmp_path/'run',evaluator=evaluator,clock=clock)
    campaign.run(100,maximum_candidates=2)
    assert len(evaluator.calls)==1
    assert campaign.history.load()[1]['cache_hit']


def test_incompatible_source_refuses_resume_but_remains_inspectable(tmp_path,monkeypatch):
    d=definition(tmp_path); clock=Clock()
    c=OptimizationCampaign(d,tmp_path/'run',evaluator=Evaluator(clock),clock=clock)
    c.run(100,maximum_candidates=1)
    monkeypatch.setattr('dada_solver.research.schema.runtime_identity',lambda:dict(changed=True))
    with pytest.raises(ValueError,match='runtime changed'): CampaignDefinition.resume(tmp_path/'run')
    assert len(inspect(tmp_path/'run')['records'])==1
    source=tmp_path/'run'/'basis.json'; source.write_text(source.read_text()+' ')
    with pytest.raises(ValueError,match='SHA-256'): CampaignDefinition.resume(tmp_path/'run')


def test_cli_validation_and_no_budget_run(tmp_path,capsys,monkeypatch):
    path=tmp_path/'study.toml'
    assert main(['init','sixbar-thermo5d','--output',str(path)])==0
    monkeypatch.setattr('dada_solver.campaign.evaluator.solve_periodic_wall_motor',lambda *a,**k:pytest.fail('solver called'))
    assert main(['validate',str(path),'--json'])==0
    assert '"integration_started": false' in capsys.readouterr().out
    assert main(['run',str(path),'--directory',str(tmp_path/'run'),'--budget','0'])==0
    assert inspect(tmp_path/'run')['records']==[]
    assert main(['resume',str(tmp_path/'run'),'--budget','0'])==0
    assert main(['status',str(tmp_path/'run')])==0


def test_run_default_directory_and_resume_from_toml(tmp_path,monkeypatch,capsys):
    path=initialize(tmp_path/'study.toml')
    monkeypatch.setattr('dada_solver.campaign.evaluator.solve_periodic_wall_motor',lambda *a,**k:pytest.fail('solver called'))
    assert main(['run',str(path),'--budget','0'])==0
    campaign=tmp_path/'campaign'
    assert inspect(campaign)['records']==[]
    # Resume uses stored scientific inputs, never reloads the editable TOML.
    path.write_text('invalid TOML')
    assert main(['resume',str(path),'--budget','0'])==0
    with pytest.raises(SystemExit) as error:
        main(['run',str(path),'--budget','0'])
    assert error.value.code==2
    assert 'use resume' in capsys.readouterr().err


def test_resume_toml_missing_campaign_is_explicit(tmp_path,capsys):
    with pytest.raises(SystemExit) as error:
        main(['resume',str(tmp_path/'study.toml'),'--budget','0'])
    assert error.value.code==2
    message=capsys.readouterr().err
    assert 'Campaign directory not found' in message
    assert str(tmp_path/'campaign') in message
    assert not (tmp_path/'campaign').exists()
