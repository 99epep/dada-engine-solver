"""Capacity scaling preserves ownership and explicit physical boundaries."""
from dataclasses import asdict
import json
import math
from pathlib import Path
import tomllib

import numpy as np
import pytest

from dada_solver.research.presets import initialize_v3, initialize_v2
from dada_solver.research.schema import load_study, compile_study, candidate_for_values
from dada_solver.research.study_io import dumps
from dada_solver.research.rescale import rescale
from dada_solver.research import report


def snapshot(path, directory):
    study = load_study(path); definition = compile_study(study)
    directory.mkdir()
    definition.write_snapshots(directory)
    (directory/'definition.json').write_text(json.dumps(dict(definition.identity, definition_id=definition.definition_id)))
    candidate = candidate_for_values(definition, {p.name:p.initial for p in study.space.parameters})
    record = dict(candidate.payload, candidate_id=candidate.candidate_id, status='feasible',
        metrics={}, constraints=[], duration_seconds=0, evaluation_number=0)
    (directory/'history.jsonl').write_text(json.dumps(record)+'\n')
    return study, definition, record


@pytest.fixture
def source(tmp_path):
    path = initialize_v3(tmp_path/'source.toml')
    raw = tomllib.loads(path.read_text())
    raw['policies']['outlet_valve_cda'] = 'fixed_source_cda'
    for row in raw['parameters']:
        if row['name'] in ('microtube.heat_in.tube_count', 'charge.total_mass_kg'):
            value = row.pop('value')
            row.update(initial=value, lower=value//2 if type(value) is int else value/2,
                upper=value*2, kind='integer' if type(value) is int else 'continuous')
            row.update(encoding='nearest_even_v1') if type(value) is int else row.update(transform='log')
    path.write_text(dumps(raw))
    study, definition, record = snapshot(path, tmp_path/'campaign')
    return tmp_path/'campaign', study, definition, record


def initial_design(path):
    study = load_study(path); d = compile_study(study)
    return study, d.adapter.build(dict(study.fixed_parameters, **{p.name:p.initial for p in study.space.parameters}))


def test_factor_one_exact_physics_and_no_source_mutation(source, tmp_path):
    directory, original, definition, record = source
    before = {p:p.read_bytes() for p in directory.iterdir()}
    output = rescale(directory, record['candidate_id'][:8], 1, tmp_path/'one.toml')
    scaled, actual = initial_design(output)
    expected = definition.adapter.build(dict(original.fixed_parameters, **record['physical']))
    assert asdict(actual.configuration) == asdict(expected.configuration)
    assert asdict(actual.heat_in) == asdict(expected.heat_in)
    assert asdict(actual.heat_out) == asdict(expected.heat_out)
    for angle in np.linspace(0, 2*math.pi, 73):
        assert actual.kinematics.small_cylinder_volume(angle) == expected.kinematics.small_cylinder_volume(angle)
        assert actual.kinematics.large_cylinder_volume(angle) == expected.kinematics.large_cylinder_volume(angle)
    assert scaled.study_id != original.study_id  # Explicit new provenance, not a cache hit.
    assert before == {p:p.read_bytes() for p in directory.iterdir()}


def test_capacity_five_basis_parameters_bounds_and_reversibility(source, tmp_path):
    directory, original, definition, record = source
    path = rescale(directory, record['candidate_id'], 5, tmp_path/'five.toml')
    scaled, design = initial_design(path)
    old = definition.adapter.build(dict(original.fixed_parameters, **record['physical']))
    assert design.configuration.charge.total_mass == old.configuration.charge.total_mass*5
    assert design.configuration.hydraulics.hot_to_small_valve_cda == old.configuration.hydraulics.hot_to_small_valve_cda*5
    assert design.configuration.angular_speed == old.configuration.angular_speed
    for side in ('heat_in','heat_out'):
        a,b = getattr(old,side),getattr(design,side)
        assert b.outlet_valve_cda_m2 == a.outlet_valve_cda_m2*5
        assert b.bank.tube_count == a.bank.tube_count*5
        assert b.bank.additional_internal_volume_m3 == a.bank.additional_internal_volume_m3*5
        assert b.inputs.extra_wall_capacity_j_k == a.inputs.extra_wall_capacity_j_k*5
        assert b.inputs.external_stream.mass_flow_kg_s == a.inputs.external_stream.mass_flow_kg_s*5
        assert b.inputs.external_stream.wall_conductance_w_k == a.inputs.external_stream.wall_conductance_w_k*5
        assert b.bank.tube_length_m == a.bank.tube_length_m
        assert b.bank.inner_diameter_m == a.bank.inner_diameter_m
    for a,b in zip(original.space.parameters,scaled.space.parameters):
        assert b.lower == a.lower*5 and b.upper == a.upper*5 and b.initial == a.initial*5
    assert scaled.data['constraints'] == original.data['constraints']
    assert scaled.data['screening'] == original.data['screening']
    provenance = scaled.basis.data['provenance']['capacity_scaling']
    assert provenance['factor'] == 5 and provenance['source_candidate_id'] == record['candidate_id']
    _,_,scaled_record = snapshot(path,tmp_path/'scaled_campaign')
    reversed_path = rescale(tmp_path/'scaled_campaign',scaled_record['candidate_id'],.2,tmp_path/'reversed.toml')
    _, reversed_design = initial_design(reversed_path)
    assert reversed_design.heat_in.outlet_valve_cda_m2 == pytest.approx(old.heat_in.outlet_valve_cda_m2,rel=1e-15)
    assert reversed_design.heat_in.bank == old.heat_in.bank


@pytest.mark.parametrize('factor',[0,-1,float('nan'),float('inf'),.00001,1.234567])
def test_unsafe_factors_leave_no_study(source,tmp_path,factor):
    directory,_,_,record=source
    with pytest.raises(ValueError): rescale(directory,record['candidate_id'],factor,tmp_path/'bad.toml')
    assert not (tmp_path/'bad.toml').exists()
    assert not (tmp_path/'bad.basis.json').exists()


def test_refuses_overwrite_and_invalid_selector(source,tmp_path):
    directory,_,_,record=source
    path=tmp_path/'existing.toml';path.write_text('keep')
    with pytest.raises(ValueError,match='already exists'):rescale(directory,record['candidate_id'],1,path)
    assert path.read_text()=='keep'
    with pytest.raises(ValueError,match='absent or ambiguous'):rescale(directory,'not-found',5,tmp_path/'new.toml')


def test_count_ratio_policy_scales_once_and_artifacts_are_portable(tmp_path):
    path=initialize_v2(tmp_path/'physical.toml','six_bar','six_bar')
    study,d,r=snapshot(path,tmp_path/'source')
    old=d.adapter.build(dict(study.fixed_parameters,**r['physical']))
    output=rescale(tmp_path/'source',r['candidate_id'],5,tmp_path/'scaled'/'study.toml')
    scaled,new=initial_design(output)
    assert len(scaled.artifacts)==2
    assert new.heat_in.outlet_valve_cda_m2==pytest.approx(old.heat_in.outlet_valve_cda_m2*5,rel=1e-15)
    assert new.heat_out.outlet_valve_cda_m2==pytest.approx(old.heat_out.outlet_valve_cda_m2*5,rel=1e-15)
    assert scaled.basis.data['warm_start']['values']==[v*5 for v in study.basis.data['warm_start']['values']]


def test_cli_and_offline_volume_plots(source,tmp_path,monkeypatch):
    directory,study,d,r=source
    from dada_solver.research.cli import main
    monkeypatch.setattr('dada_solver.campaign.evaluator.MachineEvaluator.evaluate',lambda *a:pytest.fail('integration'))
    output=tmp_path/'cli.toml'
    assert main(['rescale',str(directory),'--candidate',r['candidate_id'],'--factor','1','--output',str(output)])==0
    data=report.compare([directory],[r['candidate_id']],plots='volumes')
    plot=data['selected'][0]['plots']['volumes']
    assert not plot['thermodynamic_replay'] and len(plot['angle'])==721
    assert plot['angle'][0]==0 and plot['angle'][-1]==360
    model=d.adapter.build(dict(study.fixed_parameters,**r['physical'])).build().model
    for i in (0,123,720):
        volumes=model.volumes(math.radians(plot['angle'][i]))
        assert plot['small'][i]==pytest.approx(volumes.small_cylinder,rel=1e-14)
    html=report.render_html(data,tmp_path/'report.html').read_text()
    assert 'sortBy' in html and 'sortedRows' in html and 'topology_display' in html
    assert 'relative_margin' in html and 'data-candidate' in html
    assert '<script src=' not in html
    assert report.render_html(data,tmp_path/'report.html').read_text()==html
    with pytest.raises(ValueError,match='Only volume'):report.compare([directory],plots='temperatures')


def test_legacy_empty_event_display_preserves_original_record(source,tmp_path):
    directory,_,_,record=source
    record['diagnostics']=dict(valve_events=[],topology=dict(classification='non_nominal',
        reasons=['observed_event_sequence_differs_from_nominal_sequence']))
    journal=directory/'history.jsonl';journal.write_text(json.dumps(record)+'\n');before=journal.read_bytes()
    row=report.inspect(directory)['records'][0]
    assert row['topology_display']['classification']=='unavailable'
    assert row['diagnostics']['topology']['classification']=='non_nominal'
    assert journal.read_bytes()==before


def test_standalone_evaluation_reconstruction(source,tmp_path):
    directory,_,definition,record=source
    path=tmp_path/'evaluation.json'
    path.write_text(json.dumps(dict(artifact_type='research_evaluation_v1', name='Stored evaluation',
        definition=dict(definition.identity,definition_id=definition.definition_id),record=record)))
    data=report.compare([path],plots='volumes')
    assert data['selected'][0]['plots']['volumes']['small']
    scaled=rescale(path,record['candidate_id'],1,tmp_path/'copied.toml')
    _,design=initial_design(scaled)
    assert design.configuration.charge.total_mass==definition.configuration.charge.total_mass


def test_volume_sampling_motor_uses_production_direction_once(tmp_path):
    from dada_solver.research.visualization import sample_report_volumes
    path=initialize_v2(tmp_path/'motor.toml','hybrid_compact','hybrid_compact')
    study=load_study(path);d=compile_study(study)
    design=d.adapter.build(study.fixed_parameters)
    model=design.build().model
    plot=sample_report_volumes(study,{},samples=25)
    for i,angle in enumerate(plot['angle']):
        expected=model.volumes(math.radians(angle))
        assert plot['small'][i]==pytest.approx(expected.small_cylinder,rel=1e-14)
        assert plot['large'][i]==pytest.approx(expected.large_cylinder,rel=1e-14)
