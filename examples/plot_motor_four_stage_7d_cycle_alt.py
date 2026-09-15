"""Alternative stacked diagnostics for the retained linear 7D four-stage champion.

Compared with plot_motor_four_stage_7d_cycle.py, this version:
- removes pressure-difference traces,
- removes valve-state traces,
- removes cumulative-energy traces,
- adds instantaneous heat-transfer rates,
- adds a P-V diagram.

Outputs:
    outputs/motor_four_stage_7d_best_cycle_alt.png
    outputs/motor_four_stage_7d_best_cycle_alt.svg
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

OUTPUT_CYCLE = ROOT / "outputs" / "motor_four_stage_7d_best_cycle_alt.png"
OUTPUT_CYCLE_SVG = ROOT / "outputs" / "motor_four_stage_7d_best_cycle_alt.svg"

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
        **{
            name: float(params[name])
            for name in ("t1", "t2", "t3", "a_l", "b_l", "a_s", "b_s")
        },
    )

    charge = ChargeConfiguration(
        temperature=base_design.configuration.charge.temperature,
        pressure=CHARGE_PRESSURE_PA,
    )

    design = replace(
        base_design,
        configuration=replace(base_design.configuration, charge=charge),
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

    performance = wall_cycle_performance(wrapper, periodic.trajectory)

    return control, wrapper, periodic.angles, periodic.trajectory, performance


def _histories(wrapper, angles, trajectory):
    n = len(angles)

    volumes = np.empty((4, n))
    normalized_volumes = np.empty((2, n))
    pressures = np.empty((4, n))
    temperatures = np.empty((4, n))
    masses = np.empty((4, n))
    flows = np.empty((4, n))

    s_limits = wrapper.model.machine_volumes.small_cylinder
    l_limits = wrapper.model.machine_volumes.large_cylinder

    for i, angle in enumerate(angles):
        state = ThermodynamicState.from_array(trajectory[:8, i])
        v = wrapper.model.volumes(float(angle))

        volumes[:, i] = (
            v.small_cylinder,
            v.large_cylinder,
            v.hot_heat_exchanger,
            v.cold_heat_exchanger,
        )
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

    return {
        "volumes": volumes,
        "normalized_volumes": normalized_volumes,
        "pressures": pressures,
        "temperatures": temperatures,
        "masses": masses,
        "flows": flows,
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
        bbox_to_anchor=(0.5, -0.16),
        ncol=ncol,
        frameon=False,
    )


def _decorate_angle_axis(axis, transitions_deg):
    axis.set_xlim(0.0, 360.0)
    axis.set_xticks(np.arange(0.0, 361.0, 45.0))
    axis.grid(True, alpha=0.3)
    for angle in transitions_deg:
        axis.axvline(angle, linestyle=":", linewidth=0.9, alpha=0.65)


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
            1.03,
            label,
            ha="center",
            va="bottom",
            transform=axis.get_xaxis_transform(),
            fontsize=9,
        )


def _instantaneous_rates(wrapper, angles, trajectory):
    """Return exact exchanger-wall <-> working-gas heat-transfer rates.

    AirWallMotor.derivative() returns derivatives with respect to cycle angle.
    Multiplying by the model's positive integration angular speed recovers W.

    Quadratures:
        12: H_i wall -> working gas, positive into the gas
        13: H_o wall -> working gas, positive into the gas

    For plotting, H_o is sign-reversed so both useful transfer directions are
    positive:
        H_i -> gas
        gas -> H_o
    """
    heat_in_to_gas = np.empty(len(angles))
    gas_to_heat_out = np.empty(len(angles))

    omega = wrapper.model.angular_speed

    for i, angle in enumerate(angles):
        rates = wrapper.derivative(
            float(angle),
            np.asarray(trajectory[:, i], dtype=float),
        ) * omega

        heat_in_to_gas[i] = rates[12]
        gas_to_heat_out[i] = -rates[13]

    return heat_in_to_gas, gas_to_heat_out


def _nearest_indices(deg, transitions_deg):
    marks = [0.0, *transitions_deg]
    return [int(np.argmin(np.abs(deg - mark))) for mark in marks]

def _volumetric_flows(histories):
    """Return interface volumetric flow rates in L/s, using upstream density."""
    volumes = histories["volumes"]
    masses = histories["masses"]
    flows = histories["flows"]

    v_s = volumes[0]
    v_l = volumes[1]
    v_hi = volumes[2]
    v_ho = volumes[3]

    m_s = masses[0]
    m_l = masses[1]
    m_hi = masses[2]
    m_ho = masses[3]

    rho_s = m_s / v_s
    rho_l = m_l / v_l
    rho_hi = m_hi / v_hi
    rho_ho = m_ho / v_ho

    q_l_ho = flows[0] / rho_l
    q_s_hi = flows[1] / rho_s
    q_ho_s = flows[2] / rho_ho
    q_hi_l = flows[3] / rho_hi

    return 1000.0 * np.vstack((q_l_ho, q_s_hi, q_ho_s, q_hi_l))


def plot_cycle(control, wrapper, angles, trajectory, performance):
    plt = _setup_matplotlib()
    h = _histories(wrapper, angles, trajectory)
    deg = np.degrees(angles)

    p_kpa = h["pressures"] / 1000.0
    m_g = h["masses"] * 1000.0
    flow_g_s = h["flows"] * 1000.0
    volumetric_flow_l_s = _volumetric_flows(h)
    qidot, qodot = _instantaneous_rates(wrapper, angles, trajectory)

    motion = control["motion_parameters"]
    transitions_deg = [
        360.0 * float(motion["t1"]),
        360.0 * float(motion["t2"]),
        360.0 * float(motion["t3"]),
    ]

    fig = plt.figure(figsize=(12, 36))
    gs = fig.add_gridspec(8, 1, height_ratios=[1, 1, 1, 1, 1, 1, 1, 1.2], hspace=1.0)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(8)]
    for axis in axes[1:7]:
        axis.sharex(axes[0])

    axes[0].plot(deg, h["normalized_volumes"][0], label="S")
    axes[0].plot(deg, h["normalized_volumes"][1], label="L")
    axes[0].set_ylabel("(V - Vmin) / Vswept")
    _legend_below(axes[0])
    _phase_labels(axes[0], transitions_deg)

    for values, label in zip(
        p_kpa,
        ("S", "L", "H_i", "H_o"),
        strict=True,
    ):
        axes[1].plot(deg, values, label=label)
    axes[1].set_ylabel("Pressure (kPa)")
    _legend_below(axes[1], ncol=4)

    for values, label in zip(
        h["temperatures"],
        ("S", "L", "H_i", "H_o"),
        strict=True,
    ):
        axes[2].plot(deg, values, label=label)
    axes[2].set_ylabel("Temperature (K)")
    _legend_below(axes[2], ncol=4)

    for values, label in zip(
        m_g,
        ("S", "L", "H_i", "H_o"),
        strict=True,
    ):
        axes[3].plot(deg, values, label=label)
    axes[3].set_ylabel("Gas mass (g)")
    _legend_below(axes[3], ncol=4)

    for values, label in zip(
        flow_g_s,
        ("L -> H_o", "S -> H_i", "H_o -> S", "H_i -> L"),
        strict=True,
    ):
        axes[4].plot(deg, values, label=label)
    axes[4].axhline(0.0, linewidth=1.0)
    axes[4].set_ylabel("Mass flow (g/s)")
    _legend_below(axes[4], ncol=4)

    for values, label in zip(
        volumetric_flow_l_s,
        ("L -> H_o", "S -> H_i", "H_o -> S", "H_i -> L"),
        strict=True,
    ):
        axes[5].plot(deg, values, label=label)
    axes[5].axhline(0.0, linewidth=1.0)
    axes[5].set_ylabel("Volumetric flow (L/s)")
    _legend_below(axes[5], ncol=4)

    axes[6].plot(deg, qidot, label="H_i wall -> working gas")
    axes[6].plot(deg, qodot, label="Working gas -> H_o wall")
    axes[6].axhline(0.0, linewidth=1.0)
    axes[6].set_ylabel("Gas heat-transfer rate (W)")
    axes[6].set_xlabel("Forward motor-cycle angle (deg)")
    _legend_below(axes[6], ncol=2)

    for axis in axes[:7]:
        _decorate_angle_axis(axis, transitions_deg)

    v_s_l = h["volumes"][0] * 1000.0
    v_l_l = h["volumes"][1] * 1000.0
    pv = axes[7]
    pv.plot(v_s_l, p_kpa[0], label="S cylinder")
    pv.plot(v_l_l, p_kpa[1], label="L cylinder")

    idx_marks = _nearest_indices(deg, transitions_deg)
    labels = ("T0", "T1", "T2", "T3")
    for idx, label in zip(idx_marks, labels, strict=True):
        pv.plot(v_s_l[idx], p_kpa[0, idx], marker="o")
        pv.plot(v_l_l[idx], p_kpa[1, idx], marker="s")
        pv.annotate(label, (v_s_l[idx], p_kpa[0, idx]), xytext=(4, 4), textcoords="offset points")
        pv.annotate(label, (v_l_l[idx], p_kpa[1, idx]), xytext=(4, -10), textcoords="offset points")

    pv.set_xlabel("Cylinder volume (L)")
    pv.set_ylabel("Pressure (kPa)")
    pv.grid(True, alpha=0.3)
    _legend_below(pv, ncol=2)
    pv.set_title("P-V loops")

    title = (
        "Retained four-stage linear 7D motor champion — alternative diagnostic view\n"
        f"eta={100.0 * performance.thermal_efficiency:.5f}%   "
        f"P={performance.gas_power:.3f} W   "
        f"Qin={performance.heat_in_power:.3f} W   "
        f"T1={transitions_deg[0]:.2f}°   "
        f"T2={transitions_deg[1]:.2f}°   "
        f"T3={transitions_deg[2]:.2f}°"
    )
    fig.suptitle(title, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.988))

    OUTPUT_CYCLE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_CYCLE, dpi=160)
    fig.savefig(OUTPUT_CYCLE_SVG)
    plt.close(fig)


def main():
    control, wrapper, angles, trajectory, performance = _build_retained_design()

    print("Retained linear 7D four-stage champion — alternative view")
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
