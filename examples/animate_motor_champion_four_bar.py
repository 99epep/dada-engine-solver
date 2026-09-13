"""Animate the synthesized mechanism and its K2 volume laws (projection only)."""
from pathlib import Path
import math
import os
import tempfile

import numpy as np

from dada_solver.coupler_projection import load_shared_crank_coupler_projection_geometry
from dada_solver.geometry import CylinderVolumeLimits

ROOT = Path(__file__).resolve().parents[1]


def main():
    os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'dada_solver_matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    kin = load_shared_crank_coupler_projection_geometry(
        ROOT / 'outputs/motor_champion_four_bar_projection_geometry.toml',
        CylinderVolumeLimits(8.316831683168316e-6, .00084),
        CylinderVolumeLimits(9.900990099009901e-6, .001),
    )
    fig = plt.figure(figsize=(10, 7), facecolor='white')
    grid = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=.28)
    ax = fig.add_subplot(grid[0]); curve = fig.add_subplot(grid[1])
    ax.set_aspect('equal'); ax.set_xlim(-7, 9); ax.set_ylim(-5, 7)
    ax.axis('off')
    fig.suptitle('DADA motor — synthesized shared-crank four-bars', fontsize=15)
    ax.set_title('Projection-only output • geometry in crank-radius units', fontsize=10)
    colors = ['#2474b5', '#e87516']
    crank, = ax.plot([], [], 'o-', color='#333333', lw=3)
    ax.plot(0, 0, 'ks', ms=7)
    artists = []
    for side, color, name in zip([kin._small, kin._large], colors, ['Small cylinder S', 'Large cylinder L']):
        d = side.design
        pivot = np.array([d.pivot_x_ratio, d.pivot_y_ratio])
        axis = np.array([math.cos(math.radians(d.axis_angle_degrees)), math.sin(math.radians(d.axis_angle_degrees))])
        # The guide is a mathematical projection, not a finite connecting rod.
        offset = np.array([-axis[1], axis[0]]) * 3
        start = offset + axis * side._minimum
        end = offset + axis * side._maximum
        ax.plot(*np.vstack([start, end]).T, color=color, lw=7, alpha=.12)
        ax.plot(*pivot, 'ks', ms=6)
        links, = ax.plot([], [], 'o-', color=color, lw=2.5, ms=4)
        arm, = ax.plot([], [], '-', color=color, lw=2)
        projection, = ax.plot([], [], '--', color=color, lw=1, alpha=.6)
        piston, = ax.plot([], [], '|', color=color, ms=18, markeredgewidth=4)
        label = ax.text(*(end + np.array([0, .5])), name, color=color, fontsize=10, ha='center')
        artists.append((side, pivot, axis, offset, links, arm, projection, piston))
    phases = np.linspace(0, 2*math.pi, 721)
    for side, color, name in zip(['small', 'large'], colors, ['S', 'L']):
        curve.plot(np.degrees(phases), [1000*getattr(kin, side+'_cylinder_volume')(-t) for t in phases], color=color, label=name)
    cursor = curve.axvline(0, color='#333333', lw=1)
    curve.set(xlim=(0, 360), ylim=(0, 1.05), xlabel='Motor cycle angle (deg)', ylabel='Volume (L)')
    curve.legend(loc='upper right'); curve.grid(alpha=.2)
    fig.text(.5, .015, 'K2: 25 / 325 °C • 2 Hz • playback slowed 8× • indicated power 51.12 W • efficiency 20.42%', ha='center', fontsize=10)
    fig.subplots_adjust(bottom=.13, top=.89)

    def update(frame):
        phase = 2*math.pi*frame/120
        theta = -phase  # Exactly one motor-direction reversal.
        pin = np.array([math.cos(theta), math.sin(theta)])
        crank.set_data([0, pin[0]], [0, pin[1]])
        for side, pivot, axis, offset, links, arm, projection, piston in artists:
            d = side.design
            delta = pivot-pin; distance = np.linalg.norm(delta); unit = delta/distance
            along = (d.coupler_ratio**2-d.rocker_ratio**2+distance**2)/(2*distance)
            height = math.sqrt(d.coupler_ratio**2-along**2)
            joint = pin + along*unit + d.assembly_branch*height*np.array([-unit[1], unit[0]])
            coupler = (joint-pin)/d.coupler_ratio
            output = pin+d.output_along_ratio*coupler+d.output_normal_ratio*np.array([-coupler[1], coupler[0]])
            projected = offset+axis*np.dot(output, axis)
            links.set_data(*np.vstack([pin,joint,pivot]).T)
            arm.set_data(*np.vstack([pin,output,joint]).T)
            projection.set_data(*np.vstack([output,projected]).T)
            piston.set_data([projected[0]], [projected[1]])
        cursor.set_xdata([math.degrees(phase)]*2)

    animation = FuncAnimation(fig, update, frames=120, interval=1000/30)
    output = ROOT/'outputs/motor_champion_four_bar_k2.gif'
    animation.save(output, writer=PillowWriter(fps=30), dpi=110)
    update(25)
    fig.savefig(ROOT/'outputs/motor_champion_four_bar_k2_preview.png', dpi=110)
    plt.close(fig)
    print(output)


if __name__ == '__main__':
    main()
