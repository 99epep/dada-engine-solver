"""Failure-only instrumentation preserves rejection and records existing values."""
import json

import pytest

from dada_solver.exchangers.failure_diagnostics import microtube_failure_snapshot
from dada_solver.exchangers.gas_correlations import MicrotubeDomainError, MicrotubeGasModel
from dada_solver.exchangers.gas_film import MicrotubeGasFilm
from tests.test_microtube_gas_model import BANK
from tests.test_solver_acceleration import variable_wrapper
from tests.test_wall_backend import values
from dada_solver.wall_backend import WallRHS, WallBackendSettings


@pytest.mark.parametrize('flow,p1,p2,criterion', [
    (.00001, 4e5, 2e5, 'large_relative_pressure_drop'),
    ((5e6+1) * BANK.tube_flow_area_m2 * MicrotubeGasModel().transport.viscosity(350.) / BANK.inner_diameter_m,
     2e5, 2e5, 'reynolds_outside_correlation_domain'),
])
def test_film_snapshot(flow, p1, p2, criterion):
    film = MicrotubeGasFilm(BANK, MicrotubeGasModel(), 0.)
    expected = film.model.diagnose(BANK, flow, p1, p2, 350.)
    with pytest.raises(MicrotubeDomainError) as caught:
        film.evaluate(350., 350., dict(frequency_hz=1., passages=((flow, p1, p2),)))
    error = caught.value
    snapshot = microtube_failure_snapshot(error)
    assert criterion in str(error)
    assert snapshot['criterion'] == str(error)
    assert snapshot['mass_flow_kg_s'] == flow
    assert snapshot['p1_pa'] == p1 and snapshot['p2_pa'] == p2
    assert snapshot['temperature_k'] == 350.
    assert snapshot['tube_count'] == BANK.tube_count
    for name in ('reynolds', 'prandtl', 'mach', 'knudsen'):
        assert snapshot[name] == getattr(expected, name)
    json.dumps(snapshot, allow_nan=False)


@pytest.mark.parametrize('backend', ['python', 'numba'])
def test_rhs_first_trial_context(backend):
    if backend == 'numba': pytest.importorskip('numba')
    wrapper = variable_wrapper()
    x = values(wrapper)
    x[:2] *= 2.
    rhs = WallRHS(wrapper, WallBackendSettings(backend))
    with pytest.raises(MicrotubeDomainError) as caught:
        rhs(.15, x)
    snapshot = microtube_failure_snapshot(caught.value)
    assert snapshot['angle_rad'] == .15
    assert snapshot['time_s'] == .15 / wrapper.model.angular_speed
    assert snapshot['exchanger'] in ('heat_in', 'heat_out')
    assert snapshot['passage'] is not None
    assert snapshot['tube_count'] > 0
    assert snapshot['p1_pa'] > 0 and snapshot['temperature_k'] > 0
    json.dumps(snapshot, allow_nan=False)


def test_research_artifact_keeps_first_failure_across_retry(monkeypatch):
    from dada_solver.campaign.evaluator import MachineEvaluator
    from dada_solver.research.schema import compile_study, load_study, candidate_for_values
    from tests.test_research_sixbar import TEMPLATE, values as candidate_values
    definition = compile_study(load_study(TEMPLATE))
    calls = []
    def fail(wrapper, state, **kwargs):
        calls.append(1)
        film = MicrotubeGasFilm(BANK, MicrotubeGasModel(), 0.)
        # Distinct retry pressures prove that the first failure remains available.
        film.evaluate(350., 350., dict(frequency_hz=1.,
            passages=((.00001, (4 + len(calls))*1e5, 2e5),)))
    monkeypatch.setattr('dada_solver.campaign.evaluator.solve_periodic_wall_motor', fail)
    result = MachineEvaluator(definition).evaluate(
        candidate_for_values(definition, candidate_values(definition)))
    assert result['status'] == 'invalid_exchanger'
    assert len(calls) == 2
    assert result['diagnostics']['first_microtube_failure']['p1_pa'] == 5e5
    assert 'large_relative_pressure_drop' in result['diagnostics']['first_microtube_failure']['criterion']
    json.dumps(result, allow_nan=False)


def test_thermal_wall_snapshot_identifies_passage_and_angle(monkeypatch):
    from dataclasses import replace
    from dada_solver.dynamics import ValveTopology
    from dada_solver.state import ThermodynamicState
    from dada_solver.valves import ValveState
    wrapper = variable_wrapper()
    x = values(wrapper)
    point = wrapper.model.instantaneous_point(0., ThermodynamicState.from_array(x[:8]),
        ValveTopology(ValveState.CLOSED, ValveState.CLOSED))
    film = wrapper.heat_in.gas_film
    flow = 3000 * film.bank.tube_flow_area_m2 * film.model.transport.viscosity(350.) / film.bank.inner_diameter_m
    pressures = point.pressures.copy(); pressures[0] *= 1.5
    point = replace(point, pressures=pressures, flows=replace(point.flows, small_to_cold=flow))
    monkeypatch.setattr(type(wrapper.model), 'instantaneous_point', lambda *args: point)
    with pytest.raises(MicrotubeDomainError) as caught:
        wrapper.derivative(.25, x)
    snapshot = microtube_failure_snapshot(caught.value)
    assert snapshot['exchanger'] == 'heat_in'
    assert snapshot['passage'] == 'small_to_cold'
    assert snapshot['angle_rad'] == .25
    assert snapshot['mass_flow_kg_s'] == flow
    assert snapshot['reynolds'] == pytest.approx(3000)


@pytest.mark.parametrize('backend', ['python', 'numba'])
@pytest.mark.parametrize('diameter,criterion', [
    (.00075, 'large_relative_pressure_drop'),
    (.00080, 'large_relative_pressure_drop'),
    (.00077, 'large_relative_pressure_drop'),
])
def test_real_first_cycle_callback_rejections(backend, diameter, criterion):
    """Exercise LSODA's nested derivative callback, not only the direct RHS."""
    from pathlib import Path
    from dada_solver.campaign.evaluator import MachineEvaluator
    from dada_solver.research.schema import compile_study, load_study, candidate_for_values
    root = Path(__file__).resolve().parents[1]
    study = load_study(root/'outputs/human_cell_stage2_microtube/human_cell_stage2_microtube.toml')
    definition = compile_study(study)
    definition.wall_backend = WallBackendSettings(backend)
    physical = {p.name: p.initial for p in study.space.parameters}
    for side in ('heat_in', 'heat_out'):
        for key, value in (('tube_count', 1000), ('tube_length_m', .8), ('inner_diameter_m', diameter)):
            physical[f'microtube.{side}.{key}'] = value
    result = MachineEvaluator(definition).evaluate(candidate_for_values(definition, physical))
    assert result['status'] == 'invalid_exchanger'
    snapshot = result['diagnostics']['first_microtube_failure']
    assert snapshot['criterion'] == criterion
    assert 'snapshot_unavailable' not in snapshot
    assert snapshot['exchanger'] in ('heat_in', 'heat_out')
    assert snapshot['passage'] is not None
    assert snapshot['angle_rad'] > 0  # Failure after preflight, inside solve_ivp.
    assert snapshot['tube_count'] == 1000
    assert snapshot['inner_diameter_m'] == diameter
    assert snapshot['tube_length_m'] == .8
    if diameter == .00077:
        assert 2300 < snapshot['reynolds'] < 4000
        assert snapshot['relative_pressure_drop'] > .2
    for key in ('mass_flow_kg_s', 'reynolds', 'prandtl', 'mach', 'knudsen',
                'pressure_ratio', 'relative_pressure_drop', 'p1_pa', 'p2_pa', 'temperature_k'):
        assert snapshot[key] is not None
    assert snapshot['hydraulic_trial_reynolds'] is None
    json.dumps(snapshot, allow_nan=False)


def test_missing_optional_frame_locals_preserve_partial_snapshot():
    from types import SimpleNamespace
    # Missing locals in any recognized frame must not erase the criterion.
    trace = None
    for module, name in (
        ('external_stream', 'derivative'), ('hardware', '_gas_flow'),
        ('gas_film', 'evaluate'), ('gas_correlations', 'require')):
        frame = SimpleNamespace(f_globals={'__name__': 'dada_solver.exchangers.'+module},
            f_code=SimpleNamespace(co_name=name), f_locals={})
        trace = SimpleNamespace(tb_frame=frame, tb_next=trace)
    error = SimpleNamespace(__traceback__=trace)
    snapshot = microtube_failure_snapshot(error)
    assert snapshot['criterion'] == str(error)
    assert 'snapshot_unavailable' not in snapshot
    for key in ('exchanger', 'passage', 'angle_rad', 'mass_flow_kg_s', 'reynolds',
                'prandtl', 'mach', 'knudsen', 'pressure_ratio', 'tube_count'):
        assert snapshot[key] is None


@pytest.mark.parametrize('outcome', ['retry_interrupted', 'rollback_then_failure'])
def test_first_failure_survives_recovery_outcomes(monkeypatch, outcome):
    from types import SimpleNamespace
    from dada_solver.campaign.evaluator import MachineEvaluator
    from dada_solver.research.schema import compile_study, load_study, candidate_for_values
    from tests.test_research_sixbar import TEMPLATE, values as candidate_values
    definition = compile_study(load_study(TEMPLATE))
    calls = []
    def solve(wrapper, state, **kwargs):
        calls.append(1)
        if outcome == 'rollback_then_failure':
            kwargs['statistics_callback'](dict(phase='microtube_failure',
                snapshot=dict(criterion='first rollback rejection', mass_flow_kg_s=None)))
            raise MicrotubeDomainError('later rejection')
        if len(calls) == 1:
            raise MicrotubeDomainError('first retry rejection')
        return SimpleNamespace(status='interrupted', message='deadline', history=(),
            last_complete_state=None, converged=False, backend_statistics={})
    monkeypatch.setattr('dada_solver.campaign.evaluator.solve_periodic_wall_motor', solve)
    result = MachineEvaluator(definition).evaluate(candidate_for_values(definition, candidate_values(definition)))
    snapshot = result['diagnostics']['first_microtube_failure']
    assert snapshot['criterion'] == ('first retry rejection' if outcome == 'retry_interrupted' else 'first rollback rejection')
    assert snapshot['mass_flow_kg_s'] is None
