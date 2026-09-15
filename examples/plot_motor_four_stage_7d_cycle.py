"""Plot the retained linear 7D four-stage motor champion as stacked cycle diagnostics.

This intentionally mirrors examples/plot_best_piecewise_stageP3.py, but uses the
retained four-stage 7D motion on the K2 / 100 kPa basis.

Outputs:
    outputs/motor_four_stage_7d_best_cycle.png
    outputs/motor_four_stage_7d_best_cycle.svg

The three four-stage transition times T1/T2/T3 are drawn on every panel.
"""

from __future__ import annotations

from dataclasses import replace
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.configuration import ChargeConfiguration
from dada_solver.exchangers.air_wall import AirWallMotor
from dada_solver.exchangers.wall_cycle import (
    solve_periodic_wall_motor,
    wall_cycle_performance,
)
from dada_solver.factory import initial_valve_topology
from dada_solver.four_stage_kinematics import FourStageVolumeKinematics
from dada_solver.state import ThermodynamicState

from compare_motor_motion_laws_stage7A5 import (
    ROOT,
    A5_CAMPAIGN,
    _candidate_design_and_mass,
)
from optimize_motor_exchanger_asymmetry_stageK1 import _scaled_exchanger


THERMO3D_REPORT = ROOT / "outputs" / "motor_four_stage_thermo3d" / "report.json"

OUTPUT_CYCLE = ROOT / "outputs" / "motor_four_stage_7d_best_cycle.png"
OUTPUT_CYCLE_SVG = ROOT / "outputs" / "motor_four_stage_7d_best_cycle.svg"

CHARGE_PRESSURE_PA = 100_000.0


def _load_control() -> dict:
    if not THERMO3D_REPORT.exists():
        raise FileNotFoundError(THERMO3D_REPORT)

    report = json.loads(THERMO3D_REPORT.read_text())
    control = report.get("control_1bar")
    if control is None:
        raise RuntimeError("Thermo-3D report contains no 1-bar K2 control.")
    if not control.get("feasible"):
        raise RuntimeError("The retained 1-bar K2 control is not feasible.")
    return control


def _build_retained_design():
    control = _load_control()

    definition = CampaignDefinition(A5_CAMPAIGN)
    base_design, _legacy_mass, _saved, _a5_best, _candidate_id = (
        _candidate_design_and_mass(definition)
    )

    params = control["motion_parameters"]
    limits = base_design.configuration.machine_volumes

    kinematics = FourStageVolumeKinematics(
        limits.small_cylinder,
        limits.large_cylinder,
        **{name: float(params[name]) for name in (
            "t1", "t2", "t3", "a_l", "b_l", "a_s", "b_s"
        )},
    )

    charge = ChargeConfiguration(
        temperature=base_design.configuration.charge.temperature,
        pressure=CHARGE_PRESSURE_PA,
    )

    design = replace(
        base_design,
        configuration=replace(
            base_design.configuration,
            charge=charge,
        ),
        heat_in=_scaled_exchanger(
            base_design.heat_in,
            float(control["parameters"]["k_i"]),
        ),
        heat_out=_scaled_exchanger(
            base_design.heat_out,
            float(control["parameters"]["k_o"]),
        ),
        kinematics=kinematics,
    )

    wrapper = design.build()
    if not isinstance(wrapper, AirWallMotor):
        raise TypeError("Expected dynamic-wall AirWallMotor.")

    initial = np.asarray(control["last_complete_state"], dtype=float)

    periodic = solve_periodic_wall_motor(
        wrapper,
        initial,
        maximum_cycles=design.configuration.numerical.maximum_cycles,
        settings=definition.wall_numerical_settings,
    )
    if not periodic.converged:
        raise RuntimeError(
            f"Periodic solution did not converge: {periodic.message}"
        )

    assert periodic.angles is not None
    assert periodic.trajectory is not None

    performance = wall_cycle_performance(
        wrapper,
        periodic.trajectory,
    )

    return (
        control,
        wrapper,
        periodic.angles,
        periodic.trajectory,
        performance,
    )


def _histories(wrapper, angles, trajectory):
    n = len(angles)

    normalized_volumes = np.empty((2, n))
    pressures = np.empty((4, n))
    temperatures = np.empty((4, n))
    masses = np.empty((4, n))
    flows = np.empty((4, n))
    valve_open = np.empty((2, n), dtype=bool)

    s_limits = wrapper.model.machine_volumes.small_cylinder
    l_limits = wrapper.model.machine_volumes.large_cylinder

    for i, angle in enumerate(angles):
        state = ThermodynamicState.from_array(trajectory[:8, i])
        volumes = wrapper.model.volumes(float(angle))

        normalized_volumes[:, i] = (
            (volumes.small_cylinder - s_limits.minimum) / s_limits.swept,
            (volumes.large_cylinder - l_limits.minimum) / l_limits.swept,
        )

        pressure = state.pressures(wrapper.model.gas, volumes)
        temperature = state.temperatures(wrapper.model.gas)

        pressures[:, i] = pressure
        temperatures[:, i] = temperature
        masses[:, i] = trajectory[[0, 2, 4, 6], i]

        rates = wrapper.model.evaluate(
            float(angle),
            state,
            initial_valve_topology(),
        )
        flows[:, i] = (
            rates.flows.large_to_hot,
            rates.flows.small_to_cold,
            rates.flows.hot_to_small,
            rates.flows.cold_to_large,
        )

        # Continuous ideal-diode states.
        valve_open[0, i] = pressure[3] > pressure[0]  # H_o -> S
        valve_open[1, i] = pressure[2] > pressure[1]  # H_i -> L

    return {
        "normalized_volumes": normalized_volumes,
        "pressures": pressures,
        "temperatures": temperatures,
        "masses": masses,
        "flows": flows,
        "valve_open": valve_open,
    }


def _setup_matplotlib():
    cache = Path(tempfile.gettempdir()) / "dada_solver_matplotlib"
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _legend_below(axis, ncol=2):
    axis.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=ncol,
        frameon=False,
    )


def _decorate_angle_axis(axis, transitions_deg):
    axis.set_xlim(0.0, 360.0)
    axis.set_xticks(np.arange(0.0, 361.0, 45.0))
    axis.grid(True, alpha=0.3)
    for angle in transitions_deg:
        axis.axvline(
            angle,
            linestyle=":",
            linewidth=0.9,
            alpha=0.65,
        )


def _phase_labels(axis, transitions_deg):
    edges = [0.0, *transitions_deg, 360.0]
    labels = (
        "Low-pressure exchange",
        "Compression",
        "High-pressure exchange",
        "Expansion",
    )
    for left, right, label in zip(edges[:-1], edges[1:], labels, strict=True):
        axis.text(
            0.5 * (left + right),
            1.035,
            label,
            ha="center",
            va="bottom",
            transform=axis.get_xaxis_transform(),
            fontsize=9,
        )


def plot_cycle(control, wrapper, angles, trajectory, performance):
    plt = _setup_matplotlib()
    histories = _histories(wrapper, angles, trajectory)
    deg = np.degrees(angles)

    pressure_kpa = histories["pressures"] / 1000.0
    mass_g = histories["masses"] * 1000.0
    flow_g_s = histories["flows"] * 1000.0

    pressure_difference = np.vstack(
        (
            pressure_kpa[1] - pressure_kpa[3],  # L -> H_o
            pressure_kpa[0] - pressure_kpa[2],  # S -> H_i
            pressure_kpa[3] - pressure_kpa[0],  # H_o -> S valve
            pressure_kpa[2] - pressure_kpa[1],  # H_i -> L valve
        )
    )

    heat_in = trajectory[10] - trajectory[10, 0]
    heat_out_magnitude = -(trajectory[11] - trajectory[11, 0])
    work_out = -(trajectory[14] - trajectory[14, 0])

    motion = control["motion_parameters"]
    transitions_deg = [
        360.0 * float(motion["t1"]),
        360.0 * float(motion["t2"]),
        360.0 * float(motion["t3"]),
    ]

    fig, axes = plt.subplots(
        8,
        1,
        figsize=(12, 31),
        sharex=True,
    )

    axes[0].plot(
        deg,
        histories["normalized_volumes"][0],
        label="S",
    )
    axes[0].plot(
        deg,
        histories["normalized_volumes"][1],
        label="L",
    )
    axes[0].set_ylabel("(V - Vmin) / Vswept")
    _legend_below(axes[0])
    _phase_labels(axes[0], transitions_deg)

    for values, label in zip(
        pressure_kpa,
        ("S", "L", "H_i", "H_o"),
        strict=True,
    ):
        axes[1].plot(deg, values, label=label)
    axes[1].set_ylabel("Pressure (kPa)")
    _legend_below(axes[1], ncol=4)

    for values, label in zip(
        pressure_difference,
        ("L - H_o", "S - H_i", "H_o - S", "H_i - L"),
        strict=True,
    ):
        axes[2].plot(deg, values, label=label)
    axes[2].axhline(0.0, linewidth=1.0)
    axes[2].set_ylabel("Pressure difference (kPa)")
    _legend_below(axes[2], ncol=4)

    for values, label in zip(
        histories["temperatures"],
        ("S", "L", "H_i", "H_o"),
        strict=True,
    ):
        axes[3].plot(deg, values, label=label)
    axes[3].set_ylabel("Temperature (K)")
    _legend_below(axes[3], ncol=4)

    for values, label in zip(
        mass_g,
        ("S", "L", "H_i", "H_o"),
        strict=True,
    ):
        axes[4].plot(deg, values, label=label)
    axes[4].set_ylabel("Gas mass (g)")
    _legend_below(axes[4], ncol=4)

    for values, label in zip(
        flow_g_s,
        ("L -> H_o", "S -> H_i", "H_o -> S", "H_i -> L"),
        strict=True,
    ):
        axes[5].plot(deg, values, label=label)
    axes[5].axhline(0.0, linewidth=1.0)
    axes[5].set_ylabel("Mass flow (g/s)")
    _legend_below(axes[5], ncol=4)

    axes[6].step(
        deg,
        histories["valve_open"][0].astype(float),
        where="post",
        label="H_o -> S",
    )
    axes[6].step(
        deg,
        histories["valve_open"][1].astype(float),
        where="post",
        label="H_i -> L",
    )
    axes[6].set_ylabel("Check valve")
    axes[6].set_yticks((0, 1), ("Closed", "Open"))
    axes[6].set_ylim(-0.1, 1.1)
    _legend_below(axes[6])

    axes[7].plot(deg, heat_in, label="Heat in Q_i")
    axes[7].plot(
        deg,
        heat_out_magnitude,
        label="Heat out magnitude -Q_o",
    )
    axes[7].plot(deg, work_out, label="Indicated work out")
    axes[7].set_ylabel("Cumulative energy (J)")
    axes[7].set_xlabel("Forward motor-cycle angle (deg)")
    _legend_below(axes[7], ncol=3)

    for axis in axes:
        _decorate_angle_axis(axis, transitions_deg)

    title = (
        "Retained four-stage linear 7D motor champion — K2, 100 kPa filling\n"
        f"eta={100.0 * performance.thermal_efficiency:.5f}%   "
        f"P={performance.gas_power:.3f} W   "
        f"Qin={performance.heat_in_power:.3f} W   "
        f"T1={transitions_deg[0]:.2f}°   "
        f"T2={transitions_deg[1]:.2f}°   "
        f"T3={transitions_deg[2]:.2f}°"
    )
    fig.suptitle(title, y=0.998)
    fig.tight_layout(rect=(0, 0, 1, 0.988), h_pad=3.1)

    OUTPUT_CYCLE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_CYCLE, dpi=160)
    fig.savefig(OUTPUT_CYCLE_SVG)
    plt.close(fig)


def main():
    control, wrapper, angles, trajectory, performance = _build_retained_design()

    print("Retained linear 7D four-stage champion")
    print(json.dumps(control["motion_parameters"], indent=2))
    print(
        f"eta = {100.0 * performance.thermal_efficiency:.6f}%\n"
        f"P = {performance.gas_power:.6f} W\n"
        f"Qin = {performance.heat_in_power:.6f} W"
    )

    plot_cycle(control, wrapper, angles, trajectory, performance)

    print(f"\nWrote {OUTPUT_CYCLE.relative_to(ROOT)}")
    print(f"Wrote {OUTPUT_CYCLE_SVG.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
