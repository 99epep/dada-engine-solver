"""Plot the best Stage P3 piecewise-linear motor candidate.

Outputs:
    outputs/motor_piecewise_stageP3_best_cycle.png
    outputs/motor_piecewise_stageP3_best_pv.png

The script re-runs the best feasible P3 point with the same gas inventory and
microtube hardware used for the Stage 7A5 / harmonic comparison.

Unlike the generic WallDiagnosticCycle valve panel, the valve states plotted
here are reconstructed from the actual continuous ideal-diode pressure
conditions:
    H_o -> S open iff P_Ho > P_S
    H_i -> L open iff P_Hi > P_L
"""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path
import tempfile

import numpy as np

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.exchangers.air_wall import AirWallMotor
from dada_solver.exchangers.wall_cycle import (
    solve_periodic_wall_motor,
    wall_cycle_performance,
)
from dada_solver.factory import initial_valve_topology
from dada_solver.kinematics import IdealPiecewiseLinearVolumeKinematics
from dada_solver.state import ThermodynamicState

from compare_motor_motion_laws_stage7A5 import (
    ROOT,
    A5_CAMPAIGN,
    _candidate_design_and_mass,
    _same_inventory_design,
    _uniform_wall_initial,
)


P3_CSV = ROOT / "outputs" / "motor_piecewise_stageP3.csv"
OUTPUT_CYCLE = ROOT / "outputs" / "motor_piecewise_stageP3_best_cycle.png"
OUTPUT_PV = ROOT / "outputs" / "motor_piecewise_stageP3_best_pv.png"


def _load_best_p3() -> dict[str, float]:
    if not P3_CSV.exists():
        raise FileNotFoundError(P3_CSV)

    feasible = []
    with P3_CSV.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["feasible"].strip().lower() != "true":
                continue
            if not row["indicated_thermal_efficiency"]:
                continue
            feasible.append(row)

    if not feasible:
        raise RuntimeError("No feasible P3 candidate found.")

    best = max(
        feasible,
        key=lambda row: float(row["indicated_thermal_efficiency"]),
    )
    return {
        "small_lambda_target": float(best["small_lambda_target"]),
        "large_lambda_target": float(best["large_lambda_target"]),
        "adiabatic_sector_fraction": float(
            best["adiabatic_sector_fraction"]
        ),
        "recorded_eta": float(best["indicated_thermal_efficiency"]),
        "recorded_power_w": float(best["indicated_power_w"]),
    }


def _integrate_best():
    definition = CampaignDefinition(A5_CAMPAIGN)
    base_design, total_mass, _saved, _a5_best, _candidate_id = (
        _candidate_design_and_mass(definition)
    )
    best = _load_best_p3()

    small_limits = base_design.configuration.machine_volumes.small_cylinder
    large_limits = base_design.configuration.machine_volumes.large_cylinder

    kinematics = IdealPiecewiseLinearVolumeKinematics(
        small_limits,
        large_limits,
        best["small_lambda_target"],
        best["large_lambda_target"],
        best["adiabatic_sector_fraction"],
    )

    design = _same_inventory_design(
        base_design,
        total_mass,
        kinematics,
    )
    wrapper = design.build()
    if not isinstance(wrapper, AirWallMotor):
        raise TypeError("Expected dynamic-wall AirWallMotor.")

    initial = _uniform_wall_initial(design.configuration, wrapper)
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
    performance = wall_cycle_performance(wrapper, periodic.trajectory)

    return best, wrapper, periodic.angles, periodic.trajectory, performance


def _histories(wrapper, angles, trajectory):
    n = len(angles)

    volumes = np.empty((2, n))
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
        v = wrapper.model.volumes(float(angle))

        volumes[:, i] = (v.small_cylinder, v.large_cylinder)
        normalized_volumes[:, i] = (
            (v.small_cylinder - s_limits.minimum) / s_limits.swept,
            (v.large_cylinder - l_limits.minimum) / l_limits.swept,
        )

        p = state.pressures(wrapper.model.gas, v)
        t = state.temperatures(wrapper.model.gas)

        pressures[:, i] = p
        temperatures[:, i] = t
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

        # True continuous ideal-diode states.
        valve_open[0, i] = p[3] > p[0]  # H_o -> S
        valve_open[1, i] = p[2] > p[1]  # H_i -> L

    return {
        "volumes": volumes,
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
    import matplotlib.pyplot as plt
    return plt


def _legend_below(axis, ncol=2):
    axis.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=ncol,
        frameon=False,
    )


def _decorate_angle_axis(axis):
    axis.set_xlim(0.0, 360.0)
    axis.set_xticks(np.arange(0.0, 361.0, 45.0))
    axis.grid(True, alpha=0.3)


def plot_cycle(best, wrapper, angles, trajectory, performance):
    plt = _setup_matplotlib()
    h = _histories(wrapper, angles, trajectory)
    deg = np.degrees(angles)

    p_kpa = h["pressures"] / 1000.0
    m_g = h["masses"] * 1000.0
    flow_g_s = h["flows"] * 1000.0

    # Pressure differences in the actual hydraulic sign convention.
    dp = np.vstack(
        (
            p_kpa[1] - p_kpa[3],  # L -> H_o
            p_kpa[0] - p_kpa[2],  # S -> H_i
            p_kpa[3] - p_kpa[0],  # H_o -> S valve
            p_kpa[2] - p_kpa[1],  # H_i -> L valve
        )
    )

    qi = trajectory[10] - trajectory[10, 0]
    qo_out = -(trajectory[11] - trajectory[11, 0])
    work_out = -(trajectory[14] - trajectory[14, 0])

    fig, axes = plt.subplots(
        8,
        1,
        figsize=(12, 31),
        sharex=True,
    )

    axes[0].plot(deg, h["normalized_volumes"][0], label="S")
    axes[0].plot(deg, h["normalized_volumes"][1], label="L")
    axes[0].set_ylabel("(V - Vmin) / Vswept")
    _legend_below(axes[0])

    for values, label in zip(
        p_kpa,
        ("S", "L", "H_i", "H_o"),
        strict=True,
    ):
        axes[1].plot(deg, values, label=label)
    axes[1].set_ylabel("Pressure (kPa)")
    _legend_below(axes[1], ncol=4)

    for values, label in zip(
        dp,
        ("L - H_o", "S - H_i", "H_o - S", "H_i - L"),
        strict=True,
    ):
        axes[2].plot(deg, values, label=label)
    axes[2].axhline(0.0, linewidth=1.0)
    axes[2].set_ylabel("Pressure difference (kPa)")
    _legend_below(axes[2], ncol=4)

    for values, label in zip(
        h["temperatures"],
        ("S", "L", "H_i", "H_o"),
        strict=True,
    ):
        axes[3].plot(deg, values, label=label)
    axes[3].set_ylabel("Temperature (K)")
    _legend_below(axes[3], ncol=4)

    for values, label in zip(
        m_g,
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
        h["valve_open"][0].astype(float),
        where="post",
        label="H_o -> S",
    )
    axes[6].step(
        deg,
        h["valve_open"][1].astype(float),
        where="post",
        label="H_i -> L",
    )
    axes[6].set_ylabel("Check valve")
    axes[6].set_yticks((0, 1), ("Closed", "Open"))
    axes[6].set_ylim(-0.1, 1.1)
    _legend_below(axes[6])

    axes[7].plot(deg, qi, label="Heat in Q_i")
    axes[7].plot(deg, qo_out, label="Heat out magnitude -Q_o")
    axes[7].plot(deg, work_out, label="Indicated work out")
    axes[7].set_ylabel("Cumulative energy (J)")
    axes[7].set_xlabel("Cycle angle (deg)")
    _legend_below(axes[7], ncol=3)

    for axis in axes:
        _decorate_angle_axis(axis)

    title = (
        "Best Stage P3 piecewise-linear candidate\n"
        f"eta={100.0 * performance.thermal_efficiency:.5f}%   "
        f"P={performance.gas_power:.3f} W   "
        f"Qin={performance.heat_in_power:.3f} W   "
        f"Lambda_S={best['small_lambda_target']:.6f}   "
        f"Lambda_L={best['large_lambda_target']:.6f}   "
        f"a={best['adiabatic_sector_fraction']:.6f}"
    )
    fig.suptitle(title, y=0.998)
    fig.tight_layout(rect=(0, 0, 1, 0.988), h_pad=3.1)
    OUTPUT_CYCLE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_CYCLE, dpi=160)
    plt.close(fig)


def plot_pv(best, wrapper, angles, trajectory, performance):
    plt = _setup_matplotlib()
    h = _histories(wrapper, angles, trajectory)

    p_kpa = h["pressures"] / 1000.0
    v_l = h["volumes"][1] * 1000.0  # m3 -> L
    v_s = h["volumes"][0] * 1000.0

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(v_s, p_kpa[0], label="S cylinder")
    ax.plot(v_l, p_kpa[1], label="L cylinder")

    ax.set_xlabel("Cylinder volume (L)")
    ax.set_ylabel("Pressure (kPa)")
    ax.grid(True, alpha=0.3)
    _legend_below(ax, ncol=2)

    ax.set_title(
        "Best Stage P3 piecewise-linear P-V loops\n"
        f"eta={100.0 * performance.thermal_efficiency:.5f}%   "
        f"P={performance.gas_power:.3f} W"
    )

    fig.tight_layout(rect=(0, 0.06, 1, 1))
    OUTPUT_PV.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PV, dpi=160)
    plt.close(fig)


def main():
    best, wrapper, angles, trajectory, performance = _integrate_best()

    print("Best P3 candidate")
    print(
        f"  Lambda_S = {best['small_lambda_target']:.9f}\n"
        f"  Lambda_L = {best['large_lambda_target']:.9f}\n"
        f"  adiabatic_fraction = "
        f"{best['adiabatic_sector_fraction']:.9f}"
    )
    print(
        f"  eta = {100.0 * performance.thermal_efficiency:.6f}%\n"
        f"  P = {performance.gas_power:.6f} W\n"
        f"  Qin = {performance.heat_in_power:.6f} W"
    )

    plot_cycle(best, wrapper, angles, trajectory, performance)
    plot_pv(best, wrapper, angles, trajectory, performance)

    print(f"\nWrote {OUTPUT_CYCLE.relative_to(ROOT)}")
    print(f"Wrote {OUTPUT_PV.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
