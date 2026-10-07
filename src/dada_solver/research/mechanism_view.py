"""Interactive mechanism inspection from production closures, without a solver replay.

Lengths remain in the stored crank-radius frame. Presentation never mirrors,
rephases or infers a physical shaft layout. Matplotlib is imported on demand.
"""
from dataclasses import dataclass
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
        self.law, self.geometry = build_side(raw['settings'], raw['geometry'], 'small', CylinderVolumeLimits(1., 2.))

    def state(self, angle_rad):
        if not math.isfinite(angle_rad):
            raise ValueError('Mechanism angle must be finite.')
        if self.family == 'four_bar':
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
                    position=self.law.value(angle_rad)-1.,
                    first_derivative=self.law.value(angle_rad, 1), second_derivative=second,
                    length_unit='crank_radius', angle_domain='study_angle_before_operation_transform')


def mechanism_state(mechanism, angle_rad):
    model = mechanism if isinstance(mechanism, MechanismModel) else MechanismModel(mechanism)
    return model.state(angle_rad)


def _samples(artifact, target, side, samples):
    if side not in ('small', 'large'):
        raise ValueError('Choose SMALL or LARGE.')
    if type(samples) is not int or samples < 3:
        raise ValueError('At least three angle samples are required.')
    angles = np.asarray(target.scientific['angles_rad']) if target else np.linspace(0., 2*math.pi, samples)
    model = MechanismModel(artifact)
    states = [model.state(float(t)) for t in angles]
    return angles, states


def plot_motion_comparison(artifact, target, *, side='small', samples=361, axes=None):
    """Compare position and available velocity; acceleration is not a fit score."""
    import matplotlib.pyplot as plt
    angles, states = _samples(artifact, target, side, samples)
    if axes is None:
        figure, axes = plt.subplots(3, 1, sharex=True)
    else:
        figure = axes[0].figure
    position = np.array([s['position'] for s in states])
    velocity = np.array([s['first_derivative'] for s in states])
    axes[0].plot(angles, position, label='Mechanism')
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
    axes[0].legend()
    axes[1].legend()
    axes[-1].set_xlabel('Study angle (rad)')
    return figure, axes


def _view(artifact, target, side, samples, thermodynamic):
    import matplotlib.pyplot as plt
    from .synthesis import assess_mechanism
    from .families import side_metrics
    angles, states = _samples(artifact, target, side, samples)
    figure = plt.figure(figsize=(12, 8))
    grid = figure.add_gridspec(3, 2)
    drawing = figure.add_subplot(grid[:, 0])
    axes = [figure.add_subplot(grid[i, 1]) for i in range(3)]
    plot_motion_comparison(artifact, target, side=side, samples=samples, axes=axes)
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
    if target:
        evidence = assess_mechanism(artifact, target, side)
        metrics = evidence['mechanical']['metrics']
        heading = f"Position RMS {evidence['fit']['position_rms']:.4g}; max {evidence['fit']['position_maximum_error']:.4g}"
    else:
        model = MechanismModel(artifact)
        metrics = side_metrics(artifact.scientific['settings'], model.law, model.geometry, 1440)
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
    figure.tight_layout(rect=(0, 0, 1, .88))
    return figure, update, len(states)-1


def plot_mechanism(artifact, *, target=None, side='small', samples=361, angle_rad=0., thermodynamic=None):
    """Static synchronized view; no file export or thermodynamic execution."""
    if not math.isfinite(angle_rad):
        raise ValueError('Mechanism angle must be finite.')
    figure, update, count = _view(artifact, target, side, samples, thermodynamic)
    angles = np.asarray(target.scientific['angles_rad']) if target else np.linspace(0., 2*math.pi, samples)
    update(int(np.argmin(abs(angles-angle_rad % (2*math.pi)))))
    return figure


def animate_mechanism(artifact, *, target=None, side='small', samples=181, interval_ms=40, thermodynamic=None):
    """Return (figure, FuncAnimation); call pyplot.show() for interactive use.

    Space pauses/resumes. Keep the returned animation alive while its window is
    open. No writer, GIF or external encoder is involved.
    """
    from matplotlib.animation import FuncAnimation
    if not math.isfinite(interval_ms) or interval_ms <= 0:
        raise ValueError('Animation interval must be positive and finite.')
    figure, update, count = _view(artifact, target, side, samples, thermodynamic)
    animation = FuncAnimation(figure, update, frames=count, interval=interval_ms, blit=False)
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
