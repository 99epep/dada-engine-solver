"""Interactive mechanism inspection from production closures, without a solver replay.

Lengths remain in the stored crank-radius frame. Presentation never mirrors,
rephases or infers a physical shaft layout. Matplotlib is imported on demand.
"""
from dataclasses import dataclass
from collections import OrderedDict
import math
import numpy as np
from dada_solver.geometry import CylinderVolumeLimits
from .artifacts import MechanismArtifact
from .families import build_side


@dataclass
class MechanismModel:
    """Reconstruct once, then share geometry between samples and drawing frames."""
    artifact: MechanismArtifact

    def __post_init__(self):
        raw = self.artifact.scientific
        self.family = raw['settings']['family']
        self.primary = raw['settings'].get('component') == 'primary'
        if self.primary:
            from .synthesis_six_bar import PrimaryProjection
            self.geometry = self.artifact.reconstruct()
            self.law = PrimaryProjection(self.geometry)
        else:
            self.law, self.geometry = build_side(raw['settings'], raw['geometry'], 'small', CylinderVolumeLimits(1., 2.))

    def state(self, angle_rad):
        if not math.isfinite(angle_rad):
            raise ValueError('Mechanism angle must be finite.')
        if self.primary:
            joints = self.geometry.joint_state(angle_rad)['joints']
            links = [('A','B'), ('B','C'), ('C','D'), ('B','E'), ('C','E')]
        elif self.family == 'four_bar':
            joints = self.law.model.joint_state(angle_rad, self.law.side)['joints']
            output = self.artifact.scientific['settings']['output']
            links = [('A','B'), ('B','C'), ('C','D'), ('H','P')]
            links += [('C','H'), ('D','H')] if output == 'rocker' else [('B','H'), ('C','H')]
        elif self.family == 'six_bar':
            joints = self.geometry.joint_state(angle_rad)['joints']
            links = [('A','B'), ('B','C'), ('C','D'), ('B','E'), ('C','E'),
                     ('E','F'), ('F','G'), ('E','H'), ('F','H'), ('H','P')]
        else:
            joints = self.geometry.joint_state(angle_rad)
            links = [('A','B'), ('B','P')]
        try:
            second = self.law.value(angle_rad, 2)
        except NotImplementedError:
            second = None
        return dict(angle_rad=angle_rad, joints=joints, links=links,
                    comparison_role='E projection; not piston P' if self.primary else 'piston',
                    position=self.law.value(angle_rad)-1.,
                    first_derivative=self.law.value(angle_rad, 1), second_derivative=second,
                    length_unit='crank_radius', angle_domain='study_angle_before_operation_transform')


def mechanism_state(mechanism, angle_rad):
    model = mechanism if isinstance(mechanism, MechanismModel) else MechanismModel(mechanism)
    return model.state(angle_rad)


def _samples(artifact, target, side, samples, *, angles=None):
    if side not in ('small', 'large'):
        raise ValueError('Choose SMALL or LARGE.')
    if type(samples) is not int or samples < 3:
        raise ValueError('At least three angle samples are required.')
    if angles is None:
        angles = np.asarray(target.scientific['angles_rad']) if target else np.linspace(0., 2*math.pi, samples)
    model = MechanismModel(artifact)
    if model.primary and target:
        from .synthesis_six_bar import primary_evidence, PrimaryProjection
        model.law=PrimaryProjection(model.geometry,axis=primary_evidence(artifact,target,side)['primary_projection']['axis'])
    states = [model.state(float(t)) for t in angles]
    return angles, states


def plot_motion_comparison(artifact, target, *, side='small', samples=361, axes=None, _prepared=None):
    """Compare position and available velocity; acceleration is not a fit score."""
    import matplotlib.pyplot as plt
    angles, states = _samples(artifact, target, side, samples) if _prepared is None else (_prepared.curve_angles, _prepared.curve_states)
    if axes is None:
        figure, axes = plt.subplots(3, 1, sharex=True)
    else:
        figure = axes[0].figure
    position = np.array([s['position'] for s in states])
    velocity = np.array([s['first_derivative'] for s in states])
    label='E projection (not piston)' if artifact.scientific['settings'].get('component')=='primary' else 'Mechanism'
    axes[0].plot(angles, position, label=label)
    axes[1].plot(angles, velocity, label='Mechanism')
    if target:
        reference = target.scientific['sides'][side]
        axes[0].plot(angles, reference['position'], '--', label='Target')
        if reference['first_derivative'] is not None:
            axes[1].plot(angles, reference['first_derivative'], '--', label='Target')
        else:
            axes[1].text(.02, .9, 'Target velocity unavailable', transform=axes[1].transAxes)
        axes[2].plot(angles, position-reference['position'], label='Position error')
    else:
        axes[2].text(.02, .5, 'No target supplied', transform=axes[2].transAxes)
    for ax, label in zip(axes, ('Position / stroke', 'Velocity / rad', 'Position error')):
        ax.set_ylabel(label)
        ax.set_xlim(0., 2*math.pi)
        ax.grid(alpha=.2)
    if _prepared is not None and _prepared.model.primary:
        axes[0].set_ylabel('E projection / span')
        axes[2].set_ylabel('E projection difference')
    axes[0].legend()
    axes[1].legend()
    axes[-1].set_xlabel('Study angle (rad)')
    return figure, axes


def _view(artifact, target, side, samples, thermodynamic, *, frame_angles=None, prepared=None, figure=None):
    import matplotlib.pyplot as plt
    from .synthesis import assess_mechanism
    from .families import side_metrics
    angles, states = _samples(artifact, target, side, samples) if prepared is None else (prepared.curve_angles, prepared.curve_states)
    figure = plt.figure(figsize=(12, 8)) if figure is None else figure
    grid = figure.add_gridspec(3, 2)
    drawing = figure.add_subplot(grid[:, 0])
    axes = [figure.add_subplot(grid[i, 1]) for i in range(3)]
    plot_motion_comparison(artifact, target, side=side, samples=samples, axes=axes, _prepared=prepared)
    names = tuple(states[0]['joints'])
    cloud = np.array([[s['joints'][name] for name in names] for s in states], dtype=float)
    for index, name in enumerate(names):
        drawing.plot(cloud[:, index, 0], cloud[:, index, 1], ':', alpha=.35)
    bars = [drawing.plot([], [], 'o-', lw=2)[0] for _ in states[0]['links']]
    labels = {name: drawing.text(0, 0, name) for name in names}
    span = max(float(np.ptp(cloud[:, :, 0])), float(np.ptp(cloud[:, :, 1])), .1)
    drawing.set_xlim(float(cloud[:,:,0].min())-.1*span, float(cloud[:,:,0].max())+.1*span)
    drawing.set_ylim(float(cloud[:,:,1].min())-.1*span, float(cloud[:,:,1].max())+.1*span)
    drawing.set_aspect('equal')
    drawing.set_xlabel('Local x / crank radius')
    drawing.set_ylabel('Local y / crank radius')
    cursors = [ax.axvline(0., color='black', alpha=.5) for ax in axes]
    if prepared is not None:
        metrics = prepared.evidence.get("mechanical", {}).get("metrics", {})
        fit = prepared.evidence.get("fit")
        heading = "Stored piston fit unavailable" if not isinstance(fit, dict) else f"Recorded position RMS {fit.get('position_rms')}; max {fit.get('position_maximum_error')}"
    elif target:
        if artifact.scientific['settings'].get('component')=='primary':
            from .synthesis_six_bar import primary_evidence
            evidence=primary_evidence(artifact,target,side)
            evidence=dict(evidence,fit=evidence['primary_projection'])
        else:
            evidence = assess_mechanism(artifact, target, side)
        metrics = evidence['mechanical']['metrics']
        heading = ('E projection (not piston): ' if artifact.scientific['settings'].get('component')=='primary' else '')+f"Position RMS {evidence['fit']['position_rms']:.4g}; max {evidence['fit']['position_maximum_error']:.4g}"
    else:
        model = MechanismModel(artifact)
        metrics = model.law.metrics if model.primary else side_metrics(artifact.scientific['settings'], model.law, model.geometry, 1440)
        heading = 'Target fit unavailable'
    display = ('stroke_over_crank', 'minimum_primary_transmission_sine', 'minimum_secondary_transmission_sine',
               'minimum_transmission_sine', 'minimum_rod_axis_cosine')
    text = '; '.join(f'{key}: {metrics[key]:.4g}' for key in display if isinstance(metrics.get(key), (int, float)))
    if thermodynamic is not None:
        if not {'candidate_id', 'metrics', 'status'} <= set(thermodynamic):
            raise ValueError('Thermodynamic display requires an identified solver record.')
        text += '\nSolver '+thermodynamic['candidate_id']+': '+str(thermodynamic['metrics'])
    figure.suptitle(artifact.scientific['settings']['family']+' / '+side+'\n'+heading+'\n'+text, fontsize=9)
    angle_label = drawing.text(.02, .98, '', transform=drawing.transAxes, va='top')
    if prepared is not None:
        angles, states = prepared.frame_angles, prepared.frame_states
    elif frame_angles is not None:
        angles, states = _samples(artifact, target, side, samples, angles=frame_angles)
    def update(index):
        state = states[index]
        for bar, (left, right) in zip(bars, state['links']):
            xy = np.array([state['joints'][left], state['joints'][right]])
            bar.set_data(xy[:,0], xy[:,1])
        for name, label in labels.items():
            label.set_position(state['joints'][name])
        for cursor in cursors:
            cursor.set_xdata([angles[index], angles[index]])
        angle_label.set_text(f'Study angle {angles[index]:.3f} rad')
        return [*bars, *labels.values(), *cursors, angle_label]
    update(0)
    if prepared is None:
        figure.tight_layout(rect=(0, 0, 1, .88))
    else:
        figure.subplots_adjust(left=.08, right=.97, bottom=.14, top=.72, hspace=.4, wspace=.4)
    return figure, update, len(states)


def plot_mechanism(artifact, *, target=None, side='small', samples=361, angle_rad=0., thermodynamic=None):
    """Static synchronized view; no file export or thermodynamic execution."""
    if not math.isfinite(angle_rad):
        raise ValueError('Mechanism angle must be finite.')
    figure, update, count = _view(artifact, target, side, samples, thermodynamic)
    angles = np.asarray(target.scientific['angles_rad']) if target else np.linspace(0., 2*math.pi, samples)
    update(int(np.argmin(abs(angles-angle_rad % (2*math.pi)))))
    return figure


def animate_mechanism(artifact, *, target=None, side='small', samples=181, interval_ms=50, frames_per_cycle=66, thermodynamic=None):
    """Return (figure, FuncAnimation); call pyplot.show() for interactive use.

    Space pauses/resumes. Keep the returned animation alive while its window is
    open. No writer, GIF or external encoder is involved. Playback uses an
    independent uniform angular grid, without duplicating the periodic seam.
    Defaults are 66 frames at 20 fps (3.3 seconds per cycle); curves and
    diagnostics keep their full resolution. Only moving artists are redrawn
    when the graphical backend supports blitting.
    """
    from matplotlib.animation import FuncAnimation
    if not math.isfinite(interval_ms) or interval_ms <= 0:
        raise ValueError('Animation interval must be positive and finite.')
    if type(frames_per_cycle) is not int or frames_per_cycle < 3:
        raise ValueError('Animation frames per cycle must be an integer >= 3.')
    frame_angles = np.linspace(0., 2*math.pi, frames_per_cycle, endpoint=False)
    figure, update, count = _view(artifact, target, side, samples, thermodynamic, frame_angles=frame_angles)
    animation = FuncAnimation(figure, update, frames=range(count), interval=interval_ms,
                              blit=figure.canvas.supports_blit)
    paused = False
    def key(event):
        nonlocal paused
        if event.key == ' ':
            animation.resume() if paused else animation.pause()
            paused = not paused
    figure.canvas.mpl_connect('key_press_event', key)
    figure._mechanism_animation = animation
    return figure, animation


def plot_library_catalogue(library):
    """Unranked family catalogue for human selection; missing evidence stays absent."""
    import matplotlib.pyplot as plt
    from .synthesis import mechanism_catalogue
    rows = mechanism_catalogue(library)
    text = []
    for row in rows:
        fit = row['fit']
        text.append([row['family_id'], row['side'], row['mechanism_family'],
                     'Unavailable' if fit is None else f"{fit['position_rms']:.5g}",
                     'Unavailable' if fit is None or fit['velocity_rms_per_rad'] is None else f"{fit['velocity_rms_per_rad']:.5g}",
                     row['artifact_hash'][:12]])
    figure, axis = plt.subplots(figsize=(12, max(3, len(rows)*.4)))
    axis.axis('off')
    axis.table(cellText=text, colLabels=['Family ID', 'Piston', 'Mechanism', 'Position RMS', 'Velocity RMS / rad', 'Artifact hash'], loc='center')
    axis.set_title('Retained families — no automatic ranking')
    return figure


BROWSE_SORTS = ('library', 'position-rms', 'position-max', 'velocity-rms', 'score')
FIT_SORT_KEYS = {'position-rms': 'position_rms', 'position-max': 'position_maximum_error',
                 'velocity-rms': 'velocity_rms_per_rad'}


def browser_metric(member, side, sort):
    """Read recorded evidence only; an E projection is never a piston fit."""
    metadata = member.get('metadata', {})
    if sort == 'score':
        return metadata.get('search', {}).get('score')
    if member['mechanisms'][side]['scientific']['settings'].get('component') == 'primary':
        return None
    fit = metadata.get('evidence', {}).get(side, {}).get('fit')
    return fit.get(FIT_SORT_KEYS.get(sort)) if isinstance(fit, dict) else None


def _finite_metric(value):
    return not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(value)


def order_browser_members(library, side, sort='library'):
    """Stable recorded-metric ordering, restricted to the requested side.

    In mixed libraries primary projections stay present, at the end of fit sorts.
    Missing, nonnumeric and nonfinite metrics retain their original relative order.
    """
    if side not in ('small', 'large'):
        raise ValueError('Choose SMALL or LARGE.')
    if sort not in BROWSE_SORTS:
        raise ValueError(f'Unknown browse sort: {sort}')
    indexed = [(i, m) for i, m in enumerate(library.members) if side in m['mechanisms']]
    if not indexed:
        raise ValueError(f'No mechanism available for side {side!r}.')
    if sort in FIT_SORT_KEYS and all(m['mechanisms'][side]['scientific']['settings'].get('component') == 'primary' for _, m in indexed):
        raise ValueError(f'{sort} is a piston-fit metric; primary E projections are not piston fits.')
    def key(item):
        i, member = item
        value = browser_metric(member, side, sort)
        return (0, float(value), i) if _finite_metric(value) else (1, 0., i)
    return [m for _, m in indexed] if sort == 'library' else [m for _, m in sorted(indexed, key=key)]


class MechanismViewCache:
    """Bounded LRU of numerical data only, never figures or animations."""
    def __init__(self, capacity=5):
        if type(capacity) is not int or capacity < 1:
            raise ValueError('View cache capacity must be a positive integer.')
        self.capacity = capacity
        self.data = OrderedDict()

    def get(self, key, prepare):
        if key in self.data:
            self.data.move_to_end(key)
            return self.data[key]
        value = prepare()
        self.data[key] = value
        while len(self.data) > self.capacity:
            self.data.popitem(last=False)
        return value


@dataclass
class PreparedMechanismView:
    artifact: MechanismArtifact
    model: MechanismModel
    curve_angles: np.ndarray
    curve_states: list
    frame_angles: np.ndarray
    frame_states: list
    evidence: dict


def prepare_mechanism_view(artifact, target, side, *, samples=181, frames_per_cycle=66, evidence=None):
    """Share one production model between full curves and animation frames.

    Recorded evidence is never reevaluated. Missing mechanical diagnostics remain
    explicitly unavailable; geometry and curves still come from production states.
    """
    if type(samples) is not int or samples < 3:
        raise ValueError('At least three angle samples are required.')
    if type(frames_per_cycle) is not int or frames_per_cycle < 3:
        raise ValueError('Animation frames per cycle must be an integer >= 3.')
    model = MechanismModel(artifact)
    evidence = evidence or {}
    if model.primary:
        from .synthesis_six_bar import PrimaryProjection
        axis = evidence.get('primary_projection', {}).get('axis')
        if axis is not None:
            model.law = PrimaryProjection(model.geometry, axis=axis)
        elif target:
            from .synthesis_six_bar import primary_target_features, primary_cadence_score, PrimaryCadencePolicy
            # Missing recorded projection orientation: choose its sign explicitly
            # from current target chronology, without recomputing mechanical fit.
            from .motion_target import PeriodicTargetSide
            policy = PrimaryCadencePolicy(**artifact.data['provenance'].get('primary_cadence_policy', {}))
            features = primary_target_features(PeriodicTargetSide(target, side), policy)
            choices = [model.law, PrimaryProjection(model.geometry, sign=-1.)]
            model.law = min(choices, key=lambda p: primary_cadence_score(p, features, policy)['score'])
    curve_angles = np.asarray(target.scientific['angles_rad']) if target else np.linspace(0., 2*math.pi, samples)
    frame_angles = np.linspace(0., 2*math.pi, frames_per_cycle, endpoint=False)
    return PreparedMechanismView(artifact, model, curve_angles,
        [model.state(float(t)) for t in curve_angles], frame_angles,
        [model.state(float(t)) for t in frame_angles], evidence)


class MechanismBrowser:
    """One figure, one animation and an unchanged study angle across selections."""
    def __init__(self, library, *, target=None, side='large', family_id=None,
                 sort='library', samples=181, frames_per_cycle=66, interval_ms=50, static=False):
        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation
        from matplotlib.widgets import Button
        if not math.isfinite(interval_ms) or interval_ms <= 0:
            raise ValueError('Animation interval must be positive and finite.')
        if type(frames_per_cycle) is not int or frames_per_cycle < 3:
            raise ValueError('Animation frames per cycle must be an integer >= 3.')
        self.library, self.target, self.side = library, target, side
        self.samples, self.frames_per_cycle = samples, frames_per_cycle
        self.current_sort = sort
        self.members = order_browser_members(library, side, sort)
        if target is not None:
            for member in self.members:
                metadata = member.get('metadata', {})
                evidence = metadata.get('evidence', {}).get(side, {})
                fit = evidence.get('fit')
                hashes = [metadata.get('target_hash'), member['mechanisms'][side]['provenance'].get('target_hash')]
                if isinstance(fit, dict): hashes.append(fit.get('target_hash'))
                if any(h is not None and h != target.content_hash for h in hashes):
                    raise ValueError(f"{member['family_id']}: supplied target differs from recorded evidence target.")
        self.current_index = 0
        if family_id is not None:
            self.current_index = self._index(family_id)
        self.frame_index, self.paused = 0, bool(static)
        self.cache = MechanismViewCache(5)
        self.figure = plt.figure(figsize=(12, 8))
        self.plot_axes = []
        self._show_current()
        self.buttons = []
        for bounds, label, callback in (([.25,.025,.16,.045], 'Previous', self.previous),
                                        ([.59,.025,.16,.045], 'Next', self.next)):
            button = Button(self.figure.add_axes(bounds), label)
            button.on_clicked(lambda event, action=callback: action())
            self.buttons.append(button)
        self.counter = self.figure.text(.5,.047,'',ha='center')
        self._update_title()
        self.animation = None if static else FuncAnimation(self.figure, self.update_frame,
            frames=range(frames_per_cycle), init_func=lambda: self.update_frame(self.frame_index),
            interval=interval_ms, blit=False, cache_frame_data=False)
        self.figure._mechanism_browser = self
        self.figure._mechanism_animation = self.animation
        self.figure.canvas.mpl_connect('key_press_event', self.on_key)
        self.figure.canvas.mpl_connect('draw_event', self._honor_pause)

    def _honor_pause(self, event):
        # FuncAnimation starts its timer on the first draw. Honor an API pause
        # requested before that draw as well as pauses during interactive use.
        if self.paused and self.animation is not None and self.animation.event_source is not None:
            self.animation.event_source.stop()

    @property
    def current_family_id(self):
        return self.members[self.current_index]['family_id']

    def _index(self, family_id):
        for i, member in enumerate(self.members):
            if member['family_id'] == family_id: return i
        raise ValueError(f'No mechanism {family_id!r} available for side {self.side!r}.')

    def _show_current(self):
        member = self.members[self.current_index]
        raw = member['mechanisms'][self.side]
        evidence = member.get('metadata', {}).get('evidence', {}).get(self.side, {})
        axis = evidence.get('primary_projection', {}).get('axis')
        key = (raw['content_hash'], self.side, self.target.content_hash if self.target else None,
               self.samples, self.frames_per_cycle, None if axis is None else tuple(axis))
        self.prepared = self.cache.get(key, lambda: prepare_mechanism_view(MechanismArtifact.from_data(raw), self.target, self.side,
            samples=self.samples, frames_per_cycle=self.frames_per_cycle,
            evidence=evidence))
        for axis in self.plot_axes: axis.remove()
        remaining = set(self.figure.axes)
        _, self._update, _ = _view(self.prepared.artifact, self.target, self.side, self.samples, None,
                                  prepared=self.prepared, figure=self.figure)
        self.plot_axes = [ax for ax in self.figure.axes if ax not in remaining]
        self._update(self.frame_index)
        self._update_title()

    def _update_title(self):
        member = self.members[self.current_index]
        raw = member['mechanisms'][self.side]['scientific']
        metadata = member.get('metadata', {})
        evidence = metadata.get('evidence', {}).get(self.side, {})
        metrics = evidence.get('mechanical', {}).get('metrics', {})
        origin = metadata.get('provenance', {})
        def display(value):
            return f'{value:.5g}' if _finite_metric(value) else 'unavailable'
        parent = origin.get('primary_family_id') or metadata.get('search', {}).get('primary_family_id') or member['mechanisms'][self.side]['provenance'].get('primary_family_id')
        title = [f"{raw['settings']['family'].upper()} / {self.side.upper()} — {self.current_index+1} / {len(self.members)}",
                 self.current_family_id, f'Sort: {self.current_sort}; parent primary: {parent or "unavailable"}']
        primary = raw['settings'].get('component') == 'primary'
        if primary:
            cadence = evidence.get('primary_cadence', {})
            title += [f"E projection — not piston; primary branch: {raw['geometry'].get('primary_branch')}; "
                      f"topology/cadence score: {display(metadata.get('search', {}).get('score', cadence.get('score')))}",
                      f"Long mirror asymmetry RMS: {display(cadence.get('long_mirror_asymmetry_rms'))}"]
            keys = ('minimum_primary_transmission_sine', 'E_span_over_crank')
        else:
            fit = evidence.get('fit') or {}
            title += [f"Recorded piston fit: RMS {display(fit.get('position_rms'))}; max {display(fit.get('position_maximum_error'))}; "
                      f"velocity RMS {display(fit.get('velocity_rms_per_rad'))}",
                      f"Secondary branch: {raw['geometry'].get('second_branch', 'not applicable')}; "
                      f"kinematic search score: {display(metadata.get('search', {}).get('score'))}"]
            keys = ('minimum_primary_transmission_sine', 'minimum_secondary_transmission_sine', 'minimum_rod_axis_cosine')
        known_target = metadata.get('target_hash') or member['mechanisms'][self.side]['provenance'].get('target_hash') or (evidence.get('fit') or {}).get('target_hash')
        target_note = 'No target supplied.' if self.target is None else ('Target identity verified.' if known_target else 'Recorded target identity unavailable; fit is unverified against supplied curves.')
        title += ['; '.join(f'{key}: {display(metrics.get(key))}' for key in keys),
                  'Recorded evidence; no thermodynamic performance. ' + target_note]
        self.figure.suptitle('\n'.join(title), fontsize=9)
        if hasattr(self, 'counter'):
            self.counter.set_text(f'{self.current_index+1} / {len(self.members)}')

    def select_family(self, family_id):
        previous_index = self.current_index
        self.current_index = self._index(family_id)
        try:
            self._show_current()
        except (ValueError, ArithmeticError):
            self.current_index = previous_index
            raise
        self.figure.canvas.draw_idle()

    def next(self): self.select_family(self.members[(self.current_index+1) % len(self.members)]['family_id'])
    def previous(self): self.select_family(self.members[(self.current_index-1) % len(self.members)]['family_id'])
    def first(self): self.select_family(self.members[0]['family_id'])
    def last(self): self.select_family(self.members[-1]['family_id'])

    def set_sort(self, criterion):
        identifier = self.current_family_id
        members = order_browser_members(self.library, self.side, criterion)
        self.members, self.current_sort = members, criterion
        self.current_index = self._index(identifier)
        self._update_title()
        self.figure.canvas.draw_idle()

    def cycle_sort(self):
        available = ['library'] + [criterion for criterion in BROWSE_SORTS[1:]
            if any(_finite_metric(browser_metric(m, self.side, criterion)) for m in self.members)]
        current = available.index(self.current_sort) if self.current_sort in available else -1
        self.set_sort(available[(current+1) % len(available)])

    def toggle_pause(self):
        if self.animation is None: return
        self.animation.resume() if self.paused else self.animation.pause()
        self.paused = not self.paused

    def update_frame(self, frame_index):
        self.frame_index = int(frame_index) % self.frames_per_cycle
        return self._update(self.frame_index)

    def on_key(self, event):
        import matplotlib as mpl
        # Preserve configured pan/zoom/save shortcuts. Shifted alternatives are
        # always available when p/s/r are assigned to the Matplotlib toolbar.
        reserved = set().union(*(set(mpl.rcParams[k]) for k in ('keymap.pan', 'keymap.zoom', 'keymap.save')))
        if 'r' in mpl.rcParams['keymap.home']: reserved.add('r')
        actions = {'left': self.previous, 'right': self.next, 'p': self.previous, 'n': self.next,
                   'home': self.first, 'end': self.last, ' ': self.toggle_pause,
                   's': self.cycle_sort, 'r': lambda: self.set_sort('library'),
                   'P': self.previous, 'S': self.cycle_sort, 'R': lambda: self.set_sort('library'),
                   'shift+p': self.previous, 'shift+s': self.cycle_sort, 'shift+r': lambda: self.set_sort('library')}
        if event.key in actions and event.key not in reserved:
            actions[event.key]()


def browse_mechanisms(library, *, target=None, side='large', family_id=None, sort='library',
                      samples=181, frames_per_cycle=66, interval_ms=50, static=False):
    """Browse recorded library members in one Matplotlib window, without exports."""
    controller = MechanismBrowser(library, target=target, side=side, family_id=family_id,
        sort=sort, samples=samples, frames_per_cycle=frames_per_cycle, interval_ms=interval_ms, static=static)
    return controller.figure, controller
