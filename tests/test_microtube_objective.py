"""Manufacturing productivity uses both physical banks, regardless of ownership."""
from types import SimpleNamespace
import tomllib

import pytest

from dada_solver.campaign.objectives import (
    MaximizeCoolingPowerPerTotalMicrotube,
    MaximizeCoolingCopTimesPowerPerTotalMicrotube,
    cooling_power_per_total_microtube, cooling_cop_times_power_per_total_microtube,
)
from dada_solver.campaign.evaluator import MachineEvaluator
from dada_solver.performance import OperatingMode, ConservationReport
from dada_solver.research.presets import initialize_external_stream
from dada_solver.research.schema import load_study, compile_study, candidate_for_values
from dada_solver.research.study_schema import OBJECTIVE_UNITS
from dada_solver.research.study_io import dumps
from dada_solver.research.progress import metric_parts
from dada_solver.research import report

NAME = 'maximize_cooling_power_per_total_microtube'
COMPOSITE_NAME = 'maximize_cooling_cop_times_power_per_total_microtube'


def test_score_and_objective():
    assert cooling_power_per_total_microtube(450, 100, 200) == 1.5
    assert cooling_power_per_total_microtube(None, 100, 200) is None
    assert cooling_power_per_total_microtube(450, 0, 0) is None
    assert cooling_power_per_total_microtube(450, -100, 200) is None
    objective = MaximizeCoolingPowerPerTotalMicrotube()
    evaluation = SimpleNamespace(performance=SimpleNamespace(
        operating_mode=OperatingMode.REFRIGERATION, cooling_power=450))
    first = objective.evaluate(evaluation, heat_in_microtube_count=100, heat_out_microtube_count=200)
    assert first.available and first.value == -1.5
    evaluation.performance.cooling_power = 300
    second = objective.evaluate(evaluation, heat_in_microtube_count=50, heat_out_microtube_count=50)
    assert second.value == -3 and second.value < first.value
    from dataclasses import asdict
    from dada_solver.campaign.report import elite_records
    ranked = elite_records([
        dict(candidate_id='higher_power', status='feasible', objective=asdict(first)),
        dict(candidate_id='higher_productivity', status='feasible', objective=asdict(second)),
    ])
    assert ranked[0]['candidate_id'] == 'higher_productivity'
    assert not objective.evaluate(evaluation).available
    evaluation.performance.operating_mode = OperatingMode.MOTOR
    assert not objective.evaluate(evaluation, heat_in_microtube_count=50, heat_out_microtube_count=50).available
    evaluation.performance = None
    assert not objective.evaluate(evaluation).available


def test_composite_score_objective_unavailability_and_ranking():
    assert cooling_cop_times_power_per_total_microtube(2.0, 1.5) == 3.0
    assert cooling_cop_times_power_per_total_microtube(None, 1.5) is None
    assert cooling_cop_times_power_per_total_microtube(2.0, None) is None
    objective = MaximizeCoolingCopTimesPowerPerTotalMicrotube()
    evaluation = SimpleNamespace(performance=SimpleNamespace(
        operating_mode=OperatingMode.REFRIGERATION, cooling_power=450., cooling_cop=2.))
    result = objective.evaluate(
        evaluation, heat_in_microtube_count=100, heat_out_microtube_count=200)
    assert result.available and result.value == -3.0
    evaluation.performance.cooling_cop = None
    assert not objective.evaluate(
        evaluation, heat_in_microtube_count=100, heat_out_microtube_count=200).available
    evaluation.performance.cooling_cop = 2.
    assert not objective.evaluate(evaluation).available

    from dataclasses import asdict
    from dada_solver.campaign.report import elite_records
    def score(candidate_id, cop, productivity):
        candidate = SimpleNamespace(performance=SimpleNamespace(
            operating_mode=OperatingMode.REFRIGERATION,
            cooling_power=productivity*100, cooling_cop=cop))
        value = objective.evaluate(
            candidate, heat_in_microtube_count=50, heat_out_microtube_count=50)
        return dict(candidate_id=candidate_id, status='feasible', objective=asdict(value))

    ranked = elite_records([
        score('A', 1.5, 0.40),
        score('B', 2.2, 0.30),
        score('H_i', 3.0, 0.15),
    ])
    assert ranked[0]['candidate_id'] == 'B'


@pytest.mark.parametrize('objective,unit', OBJECTIVE_UNITS.items())
def test_schema_objective_units_and_direction(tmp_path, objective, unit):
    motor = objective in ('maximize_motor_power', 'maximize_thermal_efficiency')
    path = initialize_external_stream(tmp_path/'study.toml', mode='motor' if motor else 'refrigeration')
    raw = tomllib.loads(path.read_text())
    raw['objective'] = dict(type=objective, unit=unit)
    path.write_text(dumps(raw))
    assert compile_study(load_study(path)).objective.name == objective
    raw['objective']['unit'] = 'W' if objective == NAME else 'incorrect'
    path.write_text(dumps(raw))
    with pytest.raises(ValueError, match='objective or unit'): load_study(path)


@pytest.mark.parametrize('objective', [NAME, COMPOSITE_NAME])
def test_motor_rejects_productivity(tmp_path, objective):
    path = initialize_external_stream(tmp_path/'study.toml', mode='motor')
    raw = tomllib.loads(path.read_text())
    raw['objective'] = dict(type=objective, unit='W/microtube')
    path.write_text(dumps(raw))
    with pytest.raises(ValueError, match='operating direction'): load_study(path)


@pytest.mark.parametrize('active_sides', [(), ('heat_in',), ('heat_in', 'heat_out')])
@pytest.mark.parametrize('objective', [NAME, COMPOSITE_NAME, 'maximize_cooling_cop'])
def test_built_counts_include_fixed_coordinates(tmp_path, monkeypatch, active_sides, objective):
    path = initialize_external_stream(tmp_path/'study.toml')
    raw = tomllib.loads(path.read_text())
    raw['objective'] = dict(type=objective, unit=OBJECTIVE_UNITS[objective])
    for row in raw['parameters']:
        for side, count in [('heat_in', 100), ('heat_out', 200)]:
            if row['name'] == f'microtube.{side}.tube_count':
                row['value'] = count
                if side in active_sides:
                    row.pop('value')
                    row.update(kind='integer', encoding='nearest_even_v1', initial=count,
                               lower=count-1, upper=count+1)
    path.write_text(dumps(raw))
    definition = compile_study(load_study(path))
    from dada_solver.sizing.evaluator import EvaluationStatus
    # Production construction and assessment, replacing only the expensive integration.
    def wall(self, candidate, design, built, derived, direction, control):
        p = SimpleNamespace(operating_mode=OperatingMode.REFRIGERATION,
            cooling_power=450., cooling_cop=2., gas_power=-225., heat_in_power=450.,
            heat_out_power=-675., thermal_efficiency=None, mechanical_input_power=225.,
            conservation=ConservationReport(0., 0., 0., 0.))
        evaluation = SimpleNamespace(performance=p, usable=True, status=EvaluationStatus.CONVERGED,
                                     validity=ConservationReport(0., 0., 0., 0.))
        return self._assessment(evaluation, derived, None, None, {'cycles_completed': 1}, None)
    monkeypatch.setattr(MachineEvaluator, '_wall', wall)
    definition.constraints = ()
    result = MachineEvaluator(definition).evaluate(candidate_for_values(definition,
        {p.name:p.initial for p in definition.space.parameters}))
    assert result['derived']['heat_in_microtube_count'] == 100
    assert result['derived']['heat_out_microtube_count'] == 200
    assert result['derived']['total_microtube_count'] == 300
    assert result['metrics']['cooling_power_per_total_microtube_w'] == 1.5
    assert result['metrics']['cooling_cop_times_power_per_total_microtube_w'] == (
        result['metrics']['cooling_cop']
        * result['metrics']['cooling_power_per_total_microtube_w'])
    expected = {NAME: -1.5, COMPOSITE_NAME: -3.0, 'maximize_cooling_cop': -2.0}
    assert result['objective']['value'] == expected[objective]
    assert result['status'] == 'feasible'


def test_report_and_progress_cooling_classification(tmp_path):
    page = report.render_html(dict(name='Productivity', records=[], best=[],
        scientific=dict(objective=dict(type=NAME, unit='W/microtube'))),
        tmp_path/'report.html').read_text()
    assert '"cooling_objective": true' in page
    assert "cooling=d.cooling_objective" in page
    assert "productivity?'cooling_power_per_total_microtube_w'" in page
    assert 'Cooling power per microtube [W/microtube]' in page
    composite_page = report.render_html(dict(name='Composite', records=[], best=[],
        scientific=dict(objective=dict(type=COMPOSITE_NAME, unit='W/microtube'))),
        tmp_path/'composite.html').read_text()
    assert '"cooling_objective": true' in composite_page
    assert "composite=objective==='maximize_cooling_cop_times_power_per_total_microtube'" in composite_page
    assert 'COP×Qcold/microtube [W/microtube]' in composite_page
    parts = metric_parts(dict(
                             cooling_cop_times_power_per_total_microtube_w=3.0,
                             cooling_power_per_total_microtube_w=1.5,
                             cooling_cop=2., cooling_power_w=450.))
    assert 'COP·Q/N 3 W/µt' in parts
    assert 'Q/N 1.5 W/µt' in parts
    assert any(p.startswith('COP ') for p in parts)
    assert any(p.startswith('Qcold ') for p in parts)
    data = dict(name='Productivity', study_id='study', records=[], status_counts={},
                best=[dict(candidate_id='candidate', status='feasible', metrics=dict(
                    cooling_cop=2., cooling_power_w=450.,
                    cooling_power_per_total_microtube_w=1.5,
                    cooling_cop_times_power_per_total_microtube_w=3.0))],
                warnings=[])
    assert 'Qcold/microtube=1.5 W/microtube' in report.text_report(data)
    assert 'COP×Qcold/microtube=3.0 W/microtube' in report.text_report(data)


def test_offline_metric_uses_fixed_and_active_without_replay(tmp_path, monkeypatch):
    import json
    path = initialize_external_stream(tmp_path/'study.toml')
    definition = compile_study(load_study(path))
    # No search is needed to inspect a standalone evaluation artifact.
    candidate = candidate_for_values(definition, {})
    record = dict(candidate.payload, candidate_id=candidate.candidate_id,
                  evaluation_number=0, status='feasible', objective=dict(name='maximize_cooling_cop', value=-2., available=True),
                  constraints=[], metrics=dict(cooling_power_w=450., cooling_cop=2.), derived={})
    artifact = dict(artifact_type='research_evaluation_v1', name='Standalone evaluation',
                    definition=dict(definition.identity, definition_id=definition.definition_id), record=record)
    source = tmp_path/'evaluation.json'
    source.write_text(json.dumps(artifact))
    before = source.read_bytes()
    monkeypatch.setattr(MachineEvaluator, 'evaluate', lambda *a: pytest.fail('No replay'))
    inspected = report.inspect(source)['records'][0]
    total = sum(definition.fixed_parameters[f'microtube.{side}.tube_count']
                for side in ('heat_in', 'heat_out'))
    assert inspected['derived']['total_microtube_count'] == total
    assert inspected['metrics']['cooling_power_per_total_microtube_w'] == 450/total
    assert inspected['metrics']['cooling_cop_times_power_per_total_microtube_w'] == 2*450/total
    assert source.read_bytes() == before
