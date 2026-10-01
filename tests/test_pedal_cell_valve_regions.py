"""Offline verification of the explicit 16D-to-18D study preparation."""
import hashlib
import shutil
import pytest
from examples.prepare_pedal_cell_valve_regions import (
    prepare, SOURCE, GLOBAL, OUTPUT, CANDIDATES, PLACEMENTS, TOPOLOGIES,
)
from dada_solver.research.schema import load_study
from dada_solver.research.report import inspect
from dada_solver.campaign.scheduled_search import ScheduledSobol


def fingerprint():
    paths=list(SOURCE.rglob('*'))+[GLOBAL,GLOBAL.parent/'study.basis.json']
    return {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()}


def test_twelve_regions_portable_exact_and_sources_untouched(tmp_path):
    before=fingerprint()
    path=prepare(tmp_path/'study.toml')
    assert fingerprint()==before
    study=load_study(path); global_study=load_study(GLOBAL)
    assert len(study.space.parameters)==18
    raw=study.data; search=raw['search']; regions=search['regions']
    assert len(regions)==12 and len({r['id'] for r in regions})==12
    assert search['radius_fraction']==.20 and search['allocation']=='round_robin'
    assert search['evaluate_centers'] is True
    assert raw['execution']['default_max_candidates']==4096
    source={r['candidate_id']:r for r in inspect(SOURCE)['records']}
    for candidate in CANDIDATES:
        group=[r for r in regions if r['source_candidate_id']==candidate]
        assert len(group)==4
        assert [tuple(r['center'][p] for p in PLACEMENTS) for r in group]==list(TOPOLOGIES.values())
        for region in group:
            assert {k:v for k,v in region['center'].items() if k not in PLACEMENTS}==source[candidate]['physical']
    for actual,expected in zip(raw['parameters'],global_study.data['parameters']):
        assert {k:v for k,v in actual.items() if k!='initial'}=={k:v for k,v in expected.items() if k!='initial'}
    for key in ('policies','charge_reference','objective','constraints','numerical','screening','warm_start'):
        assert raw[key]==global_study.data[key]
    engine=ScheduledSobol(study.space,search)
    for i in range(60):
        engine.next_point()
        region=regions[i%12]
        if i<12: assert engine.last_physical==region['center']
        assert [engine.last_physical[p] for p in PLACEMENTS]==[region['center'][p] for p in PLACEMENTS]
    relocated=tmp_path/'relocated';relocated.mkdir()
    for p in tmp_path.glob('*'):
        if p.is_file(): shutil.copyfile(p,relocated/p.name)
    assert load_study(relocated/'study.toml').study_id==study.study_id
    assert load_study(OUTPUT).study_id==study.study_id
    with pytest.raises(ValueError,match='already exist'): prepare(path)
    assert fingerprint()==before
