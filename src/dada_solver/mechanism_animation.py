"""Generate publication GIFs from the configured shared-crank mechanism."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import tempfile

import numpy as np

from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_model
from dada_solver.four_bar import FourBarSliderAssembly, SharedCrankFourBarVolumeKinematics


def piston_bore(swept_volume: float, stroke: float) -> float:
    """Return the bore implied by swept volume and physical slider stroke."""

    if swept_volume <= 0.0 or stroke <= 0.0:
        raise ValueError("Swept volume and stroke must be positive.")
    return math.sqrt(4.0 * swept_volume / (math.pi * stroke))


def generate_mechanism_gif(
    configuration_path: str | Path,
    output_path: str | Path,
    *,
    frame_count: int = 120,
    frames_per_second: int = 24,
) -> None:
    """Render one kinematic revolution using the actual configured geometry."""

    if frame_count < 12:
        raise ValueError("Frame count must be at least 12.")
    if frames_per_second <= 0:
        raise ValueError("Frames per second must be positive.")
    configuration = load_simulation_configuration(configuration_path)
    kinematics = build_model(configuration).kinematics
    if not isinstance(kinematics, SharedCrankFourBarVolumeKinematics):
        raise ValueError("Animation requires shared-crank four-bar kinematics.")

    import os

    os.environ.setdefault(
        "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "dada_solver_matplotlib")
    )
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    from matplotlib.patches import Circle, Rectangle

    assemblies = (kinematics.small_assembly, kinematics.large_assembly)
    colors = ("#2474b5", "#e87516")
    names = ("Small cylinder S", "Large cylinder L")
    strokes = (kinematics.small_physical_stroke, kinematics.large_physical_stroke)
    swept = (
        configuration.machine_volumes.small_cylinder.swept,
        configuration.machine_volumes.large_cylinder.swept,
    )
    bores = tuple(piston_bore(volume, stroke) for volume, stroke in zip(swept, strokes))
    # Publication rendering deliberately enlarges cylinder diameters and places
    # both cylinders on the dry side of their piston rods. The four-bar motion
    # and normalized volume histories remain exact, but this layout is not a
    # manufacturing drawing of the configured slider constraints.
    display_bores = tuple(1.45 * bore for bore in bores)
    ground_scale = max(
        assembly.loop.rocker_pivot_x for assembly in assemblies
    )
    display_slider_ranges = tuple(
        (-0.24 * ground_scale - stroke, -0.24 * ground_scale)
        for stroke in strokes
    )

    angles = np.linspace(0.0, 2.0 * math.pi, frame_count, endpoint=False)
    all_points: list[tuple[float, float]] = [(0.0, 0.0)]
    for assembly in assemblies:
        all_points.append((assembly.loop.rocker_pivot_x, assembly.loop.rocker_pivot_y))
        for angle in angles:
            state = assembly.evaluate(*kinematics._crank(float(angle)))
            all_points.extend((state.crank_pin, state.coupler_joint, state.output_point))
    x_values = [point[0] for point in all_points]
    y_values = [point[1] for point in all_points]
    for (minimum, maximum), assembly, bore in zip(
        display_slider_ranges, assemblies, display_bores, strict=True
    ):
        clearance = 0.04 * (maximum - minimum)
        x_values.extend((minimum - clearance, maximum))
        y_values.extend(
            (
                assembly.slider.axis_origin_y - 0.55 * bore,
                assembly.slider.axis_origin_y + 0.55 * bore,
            )
        )
    span_x = max(x_values) - min(x_values)
    span_y = max(y_values) - min(y_values)
    margin = 0.08 * max(span_x, span_y)

    figure, axis = plt.subplots(figsize=(10, 6))
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlim(min(x_values) - margin, max(x_values) + margin)
    axis.set_ylim(min(y_values) - margin, max(y_values) + margin)
    axis.axis("off")

    crank_line, = axis.plot([], [], color="#333333", linewidth=4)
    link_lines = []
    rod_lines = []
    piston_lines = []
    gas_patches = []
    for index, (assembly, color, name, bore, limits) in enumerate(
        zip(
            assemblies,
            colors,
            names,
            display_bores,
            (kinematics.small_volume_limits, kinematics.large_volume_limits),
            strict=True,
        )
    ):
        pivot = (assembly.loop.rocker_pivot_x, assembly.loop.rocker_pivot_y)
        axis.add_patch(Circle(pivot, 0.012 * max(span_x, span_y), color="#222222"))
        coupler, = axis.plot([], [], color=color, linewidth=3)
        rocker, = axis.plot([], [], color=color, linewidth=3)
        output_arm, = axis.plot([], [], color=color, linewidth=1.8, linestyle="--")
        link_lines.append((coupler, rocker, output_arm))
        rod, = axis.plot([], [], color="#555555", linewidth=2)
        rod_lines.append(rod)
        piston, = axis.plot([], [], color=color, linewidth=7)
        piston_lines.append(piston)
        minimum, maximum = display_slider_ranges[index]
        clearance_length = max(
            limits.minimum / (math.pi * bore**2 / 4.0),
            0.02 * (maximum - minimum),
        )
        head_x = minimum - clearance_length
        axis.plot(
            [head_x, maximum, maximum],
            [assembly.slider.axis_origin_y - bore / 2.0] * 2
            + [assembly.slider.axis_origin_y + bore / 2.0],
            color="#777777",
            linewidth=1,
        )
        axis.plot(
            [head_x, maximum],
            [assembly.slider.axis_origin_y + bore / 2.0] * 2,
            color="#777777",
            linewidth=1,
        )
        axis.plot(
            [head_x, head_x],
            [assembly.slider.axis_origin_y - bore / 2.0, assembly.slider.axis_origin_y + bore / 2.0],
            color="#222222",
            linewidth=3,
        )
        gas = Rectangle(
            (head_x, assembly.slider.axis_origin_y - 0.45 * bore),
            0.0,
            0.9 * bore,
            facecolor=color,
            edgecolor="none",
            alpha=0.22,
        )
        axis.add_patch(gas)
        gas_patches.append((gas, head_x))

    crank_pin_marker, = axis.plot([], [], "o", color="#c62828", markersize=8)
    axis.add_patch(Circle((0.0, 0.0), 0.014 * max(span_x, span_y), color="#111111"))

    def update(frame_index: int):
        theta = float(angles[frame_index])
        states = kinematics.slider_states(theta)
        crank_pin = states[0].crank_pin
        crank_line.set_data((0.0, crank_pin[0]), (0.0, crank_pin[1]))
        crank_pin_marker.set_data((crank_pin[0],), (crank_pin[1],))
        volumes = (
            kinematics.small_cylinder_volume(theta),
            kinematics.large_cylinder_volume(theta),
        )
        artists = [crank_line, crank_pin_marker]
        for index, (assembly, state, volume) in enumerate(
            zip(assemblies, states, volumes, strict=True)
        ):
            pivot = (assembly.loop.rocker_pivot_x, assembly.loop.rocker_pivot_y)
            coupler, rocker, output_arm = link_lines[index]
            coupler.set_data(
                (state.crank_pin[0], state.coupler_joint[0]),
                (state.crank_pin[1], state.coupler_joint[1]),
            )
            rocker.set_data(
                (pivot[0], state.coupler_joint[0]),
                (pivot[1], state.coupler_joint[1]),
            )
            output_arm.set_data(
                (pivot[0], state.output_point[0]),
                (pivot[1], state.output_point[1]),
            )
            slider_y = assembly.slider.axis_origin_y
            limits = (
                kinematics.small_volume_limits,
                kinematics.large_volume_limits,
            )[index]
            minimum, maximum = display_slider_ranges[index]
            volume_fraction = (volume - limits.minimum) / limits.swept
            piston_x = minimum + volume_fraction * (maximum - minimum)
            rod_lines[index].set_data(
                (state.output_point[0], piston_x),
                (state.output_point[1], slider_y),
            )
            bore = display_bores[index]
            piston_lines[index].set_data(
                (piston_x, piston_x),
                (slider_y - 0.45 * bore, slider_y + 0.45 * bore),
            )
            gas, head_x = gas_patches[index]
            gas.set_width(max(0.0, piston_x - head_x))
            artists.extend(
                (coupler, rocker, output_arm, rod_lines[index], piston_lines[index], gas)
            )
        return artists

    animation = FuncAnimation(
        figure,
        update,
        frames=frame_count,
        interval=1000.0 / frames_per_second,
        blit=False,
    )
    animation.save(
        Path(output_path),
        writer=PillowWriter(fps=frames_per_second),
        dpi=100,
    )
    plt.close(figure)
    _rotate_and_trim_gif(Path(output_path), frames_per_second, padding=10)


def _rotate_and_trim_gif(
    path: Path,
    frames_per_second: int,
    *,
    padding: int,
) -> None:
    """Rotate counterclockwise and crop around the union of all moving frames."""

    from PIL import Image, ImageChops

    source = Image.open(path)
    frames = []
    union: tuple[int, int, int, int] | None = None
    white = Image.new("RGB", source.size, "white")
    for index in range(source.n_frames):
        source.seek(index)
        frame = source.convert("RGB")
        difference = ImageChops.difference(frame, white).convert("L")
        visible = difference.point(lambda value: 255 if value > 8 else 0)
        bounds = visible.getbbox()
        if bounds is not None:
            if union is None:
                union = bounds
            else:
                union = (
                    min(union[0], bounds[0]),
                    min(union[1], bounds[1]),
                    max(union[2], bounds[2]),
                    max(union[3], bounds[3]),
                )
        frames.append(frame.copy())
    source.close()
    if union is None:
        raise ValueError("Animation contains no visible drawing.")
    left = max(0, union[0] - padding)
    top = max(0, union[1] - padding)
    right = min(frames[0].width, union[2] + padding)
    bottom = min(frames[0].height, union[3] + padding)
    rotated = [
        frame.crop((left, top, right, bottom)).rotate(90, expand=True)
        for frame in frames
    ]
    rotated[0].save(
        path,
        save_all=True,
        append_images=rotated[1:],
        duration=round(1000.0 / frames_per_second),
        loop=0,
        disposal=2,
    )


def main(arguments: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Animate a configured shared-crank DADA mechanism."
    )
    parser.add_argument("configuration", type=Path, help="TOML configuration file")
    parser.add_argument("output", type=Path, help="Output GIF file")
    parser.add_argument("--frames", type=int, default=120)
    parser.add_argument("--fps", type=int, default=24)
    options = parser.parse_args(arguments)
    try:
        generate_mechanism_gif(
            options.configuration,
            options.output,
            frame_count=options.frames,
            frames_per_second=options.fps,
        )
        print(f"mechanism_animation = {options.output}")
        return 0
    except (OSError, ValueError) as error:
        print(f"animation_error = {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
