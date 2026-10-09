"""Library inspection preserves recorded evidence, production geometry and angle."""
import copy
import math
from types import SimpleNamespace

import numpy as np
import pytest

from dada_solver.research.artifacts import MechanismArtifact, MechanismLibrary
from dada_solver.research.mechanism_view import (
    browser_metric, order_browser_members, browse_mechanisms, MechanismViewCache)
from tests.test_research_synthesis_architecture import artifact, target_for
from tests.test_research_six_bar_synthesis import primary_library


@pytest.fixture(autouse=True)
def headless_no_integration(monkeypatch):
    mpl = pytest.importorskip('matplotlib')
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    from dada_solver.campaign.evaluator import MachineEvaluator
    monkeypatch.setattr(MachineEvaluator, 'evaluate_with_control',
                        lambda *a, **k: pytest.fail('Unexpected thermodynamic integration'))
    yield
    plt.close('all')


def member(identifier, family='slider_crank', side='large', rms=.1, score=1., phase=None):
    a = artifact(family, side)
    if phase is not None:
        raw = a.scientific
        a = MechanismArtifact.create(family, dict(raw['geometry'], phase_rad=phase), settings=raw['settings'])
    return dict(family_id=identifier, mechanisms={side: a.data}, metadata=dict(
        evidence={side: dict(fit=dict(position_rms=rms, position_maximum_error=rms,
                                    velocity_rms_per_rad=rms), mechanical=dict(metrics={}))},
        search=dict(score=score), provenance=dict(primary_family_id='parent-primary')))


def ids(rows): return [m['family_id'] for m in rows]


def test_stable_numeric_order_and_side_filter():
    rows = [member('a', rms=.3), member('b', rms=.1), member('c', rms=.1),
            member('none', rms=None), member('nan', rms=float('nan')),
            member('inf', rms=float('inf')), member('bool', rms=True), member('missing'),
            member('other-side', side='small')]
    del rows[-2]['metadata']['evidence']['large']['fit']['position_rms']
    lib = MechanismLibrary(tuple(rows))
    assert ids(order_browser_members(lib, 'large')) == ids(rows[:-1])
    assert ids(order_browser_members(lib, 'large', 'position-rms')) == ['b','c','a','none','nan','inf','bool','missing']
    for sort in ('position-max', 'velocity-rms', 'score'):
        assert len(order_browser_members(lib, 'large', sort)) == 8
    with pytest.raises(ValueError, match='Unknown browse sort'): order_browser_members(lib, 'large', 'unknown')
    with pytest.raises(ValueError, match='No mechanism'): order_browser_members(MechanismLibrary((rows[-1],)), 'large')


def test_primary_metrics_are_never_piston_fits():
    target = target_for('six_bar')
    row = copy.deepcopy(primary_library(target).members[0])
    row['metadata']['search'] = {'score': 2.}
    other = copy.deepcopy(row); other['family_id'] = 'second'; other['metadata']['search']['score'] = 1.
    lib = MechanismLibrary((row, other))
    assert ids(order_browser_members(lib, 'large', 'score')) == ['second', row['family_id']]
    for sort in ('position-rms','position-max','velocity-rms'):
        assert browser_metric(row, 'large', sort) is None
        with pytest.raises(ValueError, match='primary E projections are not piston fits'):
            order_browser_members(lib, 'large', sort)
    mixed = MechanismLibrary((row, member('piston', rms=.5)))
    assert ids(order_browser_members(mixed, 'large', 'position-rms')) == ['piston',row['family_id']]


def test_navigation_keeps_angle_pause_and_single_animation():
    lib = MechanismLibrary(tuple(member(name, phase=i*.2, rms=value) for i,(name,value) in enumerate([('a',.3),('b',.1),('c',.2)])))
    figure, b = browse_mechanisms(lib)
    animation = b.animation
    assert animation.event_source.interval == 50 and list(animation.new_frame_seq()) == list(range(66))
    assert animation._blit is False
    animation._func(17)
    angle = b.prepared.frame_angles[b.frame_index]
    b.toggle_pause(); assert b.paused
    b.next(); assert b.current_family_id == 'b' and b.frame_index == 17 and b.paused
    assert b.prepared.frame_angles[b.frame_index] == angle
    assert b.animation is animation and b.figure is figure
    b.set_sort('position-rms'); assert b.current_family_id == 'b' and b.current_index == 0
    assert b.frame_index == 17 and b.paused
    b.previous(); assert b.current_family_id == 'a'
    b.next(); assert b.current_family_id == 'b'
    b.last(); assert b.current_family_id == 'a'
    b.first(); assert b.current_family_id == 'b'
    b.on_key(SimpleNamespace(key='right')); assert b.current_family_id == 'c'
    b.on_key(SimpleNamespace(key='left')); assert b.current_family_id == 'b'
    b.on_key(SimpleNamespace(key='end')); assert b.current_family_id == 'a'
    b.on_key(SimpleNamespace(key='home')); assert b.current_family_id == 'b'
    b.last()
    b.on_key(SimpleNamespace(key='R')); assert b.current_sort == 'library' and b.current_family_id == 'a'
    b.on_key(SimpleNamespace(key='S')); assert b.current_sort == 'position-rms'
    assert b.counter.get_text() == '3 / 3'
    assert 'parent-primary' in figure._suptitle.get_text()
    b.toggle_pause(); assert not b.paused
    figure.canvas.draw()
    assert b.frame_index == 17
    for cursor in b._update(17)[-4:-1]: np.testing.assert_allclose(cursor.get_xdata(), [angle,angle])
    animation._draw_was_started = True


def test_toolbar_keys_are_preserved(monkeypatch):
    import matplotlib as mpl
    _, b = browse_mechanisms(MechanismLibrary((member('a'), member('b'))), static=True)
    for key, config in [('p','keymap.pan'),('s','keymap.save'),('r','keymap.home')]:
        monkeypatch.setitem(mpl.rcParams, config, [key])
        b.on_key(SimpleNamespace(key=key))
        assert b.current_family_id == 'a' and b.current_sort == 'library'
    b.on_key(SimpleNamespace(key='P')); assert b.current_family_id == 'b'
    b.on_key(SimpleNamespace(key=' ')); assert b.paused and b.animation is None


def test_lru_capacity_and_revisits_do_not_reconstruct(monkeypatch):
    import dada_solver.research.mechanism_view as view
    lib = MechanismLibrary(tuple(member(str(i), phase=i*.2) for i in range(7)))
    original = view.prepare_mechanism_view; calls = []
    def tracked(*args, **kwargs):
        calls.append(args[0].content_hash)
        return original(*args, **kwargs)
    monkeypatch.setattr(view, 'prepare_mechanism_view', tracked)
    _, b = browse_mechanisms(lib, samples=11, static=True)
    b.next(); b.previous(); assert len(calls) == 2
    for _ in range(6): b.next()
    assert len(b.cache.data) == 5 and len(calls) == 7
    b.previous(); assert len(calls) == 7
    b.first(); assert len(calls) == 8
    assert len(b.cache.data) == 5
    # Artifact validation also must not rebuild production closure on a hit.
    monkeypatch.setattr(MechanismArtifact, 'from_data', lambda *a, **k: pytest.fail('Cache hit reconstructed an artifact'))
    b.first()


def test_cache_shares_identical_artifacts_and_propagates_preparation_failure():
    cache = MechanismViewCache(2)
    first = object()
    assert cache.get('a', lambda: first) is first
    assert cache.get('a', lambda: pytest.fail('prepared twice')) is first
    with pytest.raises(ValueError): cache.get('b', lambda: (_ for _ in ()).throw(ValueError('bad')))
    assert list(cache.data) == ['a']
    lib = MechanismLibrary((member('a'), member('b')))
    _, b = browse_mechanisms(lib, static=True)
    b.next(); assert len(b.cache.data) == 1


def test_variable_topologies_replace_artists_and_full_target_grid():
    target = target_for('six_bar')
    primary = copy.deepcopy(primary_library(target).members[0])
    primary['metadata']['search'] = {'score': .2}
    second_primary = copy.deepcopy(primary); second_primary['family_id'] = 'second-primary'
    raw = second_primary['mechanisms']['large']['scientific']
    second_primary['mechanisms']['large'] = MechanismArtifact.create('six_bar',
        dict(raw['geometry'], primary_e_along=2.7), settings=raw['settings']).data
    second_six = copy.deepcopy(member('second-six', 'six_bar'))
    raw = second_six['mechanisms']['large']['scientific']
    second_six['mechanisms']['large'] = MechanismArtifact.create('six_bar',
        dict(raw['geometry'], primary_phase=.6), settings=raw['settings']).data
    rows = (primary, second_primary, member('six', 'six_bar'), second_six, member('four', 'four_bar'), member('slider'))
    figure, b = browse_mechanisms(MechanismLibrary(rows), target=target, samples=11, static=True)
    for count in (5,5,9,9,6,3):
        assert len(b.prepared.frame_states[0]['joints']) == count
        assert len(b.plot_axes) == 4 and len(figure.axes) == 6
        assert len(b.plot_axes[1].lines[0].get_xdata()) == len(target.scientific['angles_rad'])
        assert len(b._update(0)) == len(b.prepared.frame_states[0]['links']) + count + 4
        figure.canvas.draw()
        b.next()
    assert 'E projection — not piston' in figure._suptitle.get_text()


def test_target_identity_and_initial_selection():
    target = target_for('slider_crank')
    row = member('a'); row['metadata']['target_hash'] = 'another-target'
    with pytest.raises(ValueError, match='supplied target differs'):
        browse_mechanisms(MechanismLibrary((row,)), target=target)
    lib = MechanismLibrary((member('a',rms=.2),member('b',rms=.1)))
    _, b = browse_mechanisms(lib, sort='position-rms', family_id='a', static=True)
    assert b.current_index == 1
    with pytest.raises(ValueError, match='No mechanism'):
        browse_mechanisms(lib, family_id='absent')


def test_dense_curves_and_pause_before_first_draw(monkeypatch):
    from dada_solver.research.motion_target import MotionTarget
    angles = np.linspace(0., 2*math.pi, 1001)
    row = dict(position=(.5+.5*np.cos(angles)).tolist(),
               first_derivative=(-.5*np.sin(angles)).tolist(), second_derivative=None, events=[])
    target = MotionTarget.create(angles, dict(small=row, large=row), source={'identity':'dense-browser-test'})
    figure, b = browse_mechanisms(MechanismLibrary((member('a'),)), target=target)
    assert len(b.prepared.curve_states) == 1001
    assert len(b.prepared.frame_states) == 66
    assert b.prepared.frame_angles[-1] < 2*math.pi
    stopped = []
    monkeypatch.setattr(b.animation.event_source, 'stop', lambda: stopped.append(True))
    b.toggle_pause()
    figure.canvas.draw()
    assert b.paused and len(stopped) >= 2
    assert 'unverified' in figure._suptitle.get_text()


def test_cli_browser_and_errors(tmp_path, capsys):
    from dada_solver.research.cli import main
    path = tmp_path/'library.json'
    MechanismLibrary((member('a',rms=.2),member('b',rms=.1))).save(path)
    original_bytes = path.read_bytes()
    args = ['mechanism','visualize',str(path),'--side','large','--no-show']
    assert main([*args,'--browse']) == 0
    assert 'current family a' in capsys.readouterr().out
    assert main([*args,'--browse','--sort','position-rms']) == 0
    assert 'current family b' in capsys.readouterr().out
    assert main([*args,'--browse','--sort','position-rms','--family-id','a','--static']) == 0
    assert 'current family a' in capsys.readouterr().out
    for extra in ([], ['--sort','position-rms']):
        with pytest.raises(SystemExit) as error: main([*args,*extra])
        assert error.value.code == 2
    artifact_path = tmp_path/'artifact.json'; artifact('slider_crank').save(artifact_path)
    with pytest.raises(SystemExit) as error:
        main(['mechanism','visualize',str(artifact_path),'--browse','--no-show'])
    assert error.value.code == 2
    assert path.read_bytes() == original_bytes
