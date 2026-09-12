"""Thermodynamic/hydraulic diagnostic for the Stage F6 motor champion.

This is a diagnostic-only script. It re-integrates the best feasible Stage F6
candidate with the same gas inventory, hardware and numerical settings used by
the optimization. It does not change any physical equation.

Outputs
-------
    outputs/motor_two_extrema_stageF6_thermo_volumes.png
    outputs/motor_two_extrema_stageF6_thermo_volume_derivatives.png
    outputs/motor_two_extrema_stageF6_thermo_pressures.png
    outputs/motor_two_extrema_stageF6_thermo_pressure_differences.png
    outputs/motor_two_extrema_stageF6_thermo_temperatures.png
    outputs/motor_two_extrema_stageF6_thermo_temperature_differences.png
    outputs/motor_two_extrema_stageF6_thermo_mass_flows.png
    outputs/motor_two_extrema_stageF6_thermo_masses.png
    outputs/motor_two_extrema_stageF6_thermo_heat_rates.png
    outputs/motor_two_extrema_stageF6_thermo_valves.png
    outputs/motor_two_extrema_stageF6_thermo_pv_S.png
    outputs/motor_two_extrema_stageF6_thermo_pv_L.png
    outputs/motor_two_extrema_stageF6_thermo_sector_mass_transfer.csv
    outputs/motor_two_extrema_stageF6_thermo_summary.json
    outputs/motor_two_extrema_stageF6_thermo_trajectory.npz

Every angle-based plot carries thin vertical lines at the sampled pressure-zero
crossings of the two passive ideal-diode valves. P-V plots use arrows for cycle
direction and theta-labelled reference points.

Important interpretation note
-----------------------------
V_S + V_L and d(V_S + V_L)/dtheta are kinematic indicators only. The actual
movement through the exchangers is diagnosed from the four simulated port
mass-flow histories and from the exchanger control-volume mass changes.

Solver state naming is S, L, C, H. In motor/exchanger plots below:
    C == heat-in exchanger gas control volume H_i
    H == heat-out exchanger gas control volume H_o
"""

from __future__ import annotations

import csv
import json
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
from dada_solver.state import ThermodynamicState

from compare_motor_motion_laws_stage7A5 import (
    ROOT,
    A5_CAMPAIGN,
    _candidate_design_and_mass,
    _same_inventory_design,
    _uniform_wall_initial,
)
from optimize_motor_two_extrema_stageF6 import (
    LARGE_MIN_DEG,
    PARAMETERS,
    _make_kinematics,
)


F6_JSON = ROOT / "outputs" / "motor_two_extrema_stageF6.json"

PREFIX = ROOT / "outputs" / "motor_two_extrema_stageF6_thermo"
OUTPUTS = {
    "volumes": Path(str(PREFIX) + "_volumes.png"),
    "volume_derivatives": Path(str(PREFIX) + "_volume_derivatives.png"),
    "pressures": Path(str(PREFIX) + "_pressures.png"),
    "pressure_differences": Path(str(PREFIX) + "_pressure_differences.png"),
    "temperatures": Path(str(PREFIX) + "_temperatures.png"),
    "temperature_differences": Path(str(PREFIX) + "_temperature_differences.png"),
    "mass_flows": Path(str(PREFIX) + "_mass_flows.png"),
    "masses": Path(str(PREFIX) + "_masses.png"),
    "heat_rates": Path(str(PREFIX) + "_heat_rates.png"),
    "valves": Path(str(PREFIX) + "_valves.png"),
    "pv_s": Path(str(PREFIX) + "_pv_S.png"),
    "pv_l": Path(str(PREFIX) + "_pv_L.png"),
    "sector_csv": Path(str(PREFIX) + "_sector_mass_transfer.csv"),
    "summary_json": Path(str(PREFIX) + "_summary.json"),
    "trajectory": Path(str(PREFIX) + "_trajectory.npz"),
}


def _load_f6_best() -> dict[str, float]:
    data = json.loads(F6_JSON.read_text())
    best = data.get("best_feasible")
    if best is None:
        raise RuntimeError("F6 contains no feasible champion.")
    return best


def _integrate_best():
    definition = CampaignDefinition(A5_CAMPAIGN)
    base_design, total_mass, _saved, _a5_best, _candidate_id = (
        _candidate_design_and_mass(definition)
    )
    best = _load_f6_best()

    values = {name: float(best[name]) for name in PARAMETERS}
    small_limits = base_design.configuration.machine_volumes.small_cylinder
    large_limits = base_design.configuration.machine_volumes.large_cylinder
    kinematics = _make_kinematics(small_limits, large_limits, values)
    design = _same_inventory_design(base_design, total_mass, kinematics)

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

    return (
        best,
        values,
        wrapper,
        periodic.angles,
        periodic.trajectory,
        performance,
    )


def _histories(wrapper, angles, trajectory):
    n = len(angles)

    volumes = np.empty((2, n))
    dvol_dtheta = np.empty((2, n))
    pressures = np.empty((4, n))
    temperatures = np.empty((4, n))
    masses = np.empty((4, n))
    flows = np.empty((4, n))
    valve_open = np.empty((2, n), dtype=bool)

    wall_temperature = np.empty((2, n))
    gas_heat = np.empty((2, n))
    air_heat = np.empty((2, n))
    air_outlet_temperature = np.empty((2, n))
    work_rate = np.empty(n)

    topology = initial_valve_topology()

    for i, angle in enumerate(angles):
        a = float(angle)
        state = ThermodynamicState.from_array(trajectory[:8, i])
        v = wrapper.model.volumes(a)
        t, p = state.temperatures_and_pressures(wrapper.model.gas, v)

        volumes[:, i] = (v.small_cylinder, v.large_cylinder)
        dvol_dtheta[:, i] = (
            wrapper.model.kinematics.small_cylinder_volume_derivative(a),
            wrapper.model.kinematics.large_cylinder_volume_derivative(a),
        )
        pressures[:, i] = p
        temperatures[:, i] = t
        masses[:, i] = trajectory[[0, 2, 4, 6], i]

        rates = wrapper.model.evaluate(a, state, topology)
        flows[:, i] = (
            rates.flows.small_to_cold,
            rates.flows.cold_to_large,
            rates.flows.large_to_hot,
            rates.flows.hot_to_small,
        )
        work_rate[i] = rates.gas_work_rate

        valve_open[0, i] = p[2] > p[1]  # H_i/C -> L
        valve_open[1, i] = p[3] > p[0]  # H_o/H -> S

        hi = wrapper.heat_in.rates(t[2], trajectory[8, i])
        ho = wrapper.heat_out.rates(t[3], trajectory[9, i])

        wall_temperature[:, i] = (
            trajectory[8, i] / wrapper.heat_in.wall_capacity_j_k,
            trajectory[9, i] / wrapper.heat_out.wall_capacity_j_k,
        )
        gas_heat[:, i] = (hi["gas_heat_w"], ho["gas_heat_w"])
        air_heat[:, i] = (hi["air_heat_w"], ho["air_heat_w"])
        air_outlet_temperature[:, i] = (
            hi["air_outlet_temperature_k"],
            ho["air_outlet_temperature_k"],
        )

    return {
        "volumes": volumes,
        "dvol_dtheta": dvol_dtheta,
        "pressures": pressures,
        "temperatures": temperatures,
        "masses": masses,
        "flows": flows,
        "valve_open": valve_open,
        "wall_temperature": wall_temperature,
        "gas_heat": gas_heat,
        "air_heat": air_heat,
        "air_outlet_temperature": air_outlet_temperature,
        "work_rate": work_rate,
    }


def _sampled_valve_events(angles, pressures):
    """Linear interpolation of ideal-diode pressure zero crossings."""
    events = []
    definitions = (
        ("H_i -> L", pressures[2] - pressures[1]),
        ("H_o -> S", pressures[3] - pressures[0]),
    )

    for valve, difference in definitions:
        for i in range(len(angles) - 1):
            if angles[i + 1] <= angles[i]:
                continue
            left = float(difference[i])
            right = float(difference[i + 1])
            if (left > 0.0) == (right > 0.0):
                continue
            denominator = right - left
            if denominator == 0.0:
                continue
            fraction = -left / denominator
            angle = float(
                angles[i] + fraction * (angles[i + 1] - angles[i])
            )
            events.append(
                {
                    "valve": valve,
                    "transition": "opening" if right > 0.0 else "closing",
                    "angle_radians": angle,
                    "angle_degrees": math.degrees(angle) % 360.0,
                    "method": (
                        "linear_interpolation_of_sampled_pressure_zero_crossing"
                    ),
                }
            )

    events.sort(key=lambda item: item["angle_degrees"])
    return events


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


def _decorate_angle_plot(axis, events):
    axis.set_xlim(0.0, 360.0)
    axis.set_xticks(np.arange(0.0, 361.0, 45.0))
    axis.set_xlabel("Study angle theta (deg)")
    axis.grid(True, alpha=0.3)

    for event in events:
        axis.axvline(
            event["angle_degrees"],
            linewidth=0.7,
            linestyle="--",
            alpha=0.28,
            label="_nolegend_",
        )


def _annotate_events_at_top(axis, events):
    ymax = axis.get_ylim()[1]
    for i, event in enumerate(events):
        label = (
            f"{event['valve']} "
            f"{event['transition']}\n"
            f"{event['angle_degrees']:.1f}°"
        )
        axis.annotate(
            label,
            xy=(event["angle_degrees"], ymax),
            xytext=(2, -4 - 20 * (i % 2)),
            textcoords="offset points",
            rotation=90,
            va="top",
            ha="left",
            fontsize=7,
            alpha=0.72,
        )


def _save_single_plot(
    path,
    title,
    ylabel,
    degrees,
    series,
    events,
    *,
    zero_line=False,
    ncol=2,
    event_labels=False,
):
    plt = _setup_matplotlib()
    fig, ax = plt.subplots(figsize=(12, 6.5))

    for values, label, linestyle in series:
        kwargs = {"label": label}
        if linestyle is not None:
            kwargs["linestyle"] = linestyle
        ax.plot(degrees, values, **kwargs)

    if zero_line:
        ax.axhline(0.0, linewidth=0.8, alpha=0.55)

    ax.set_ylabel(ylabel)
    ax.set_title(title)
    _decorate_angle_plot(ax, events)
    _legend_below(ax, ncol=ncol)
    if event_labels:
        _annotate_events_at_top(ax, events)

    fig.tight_layout(rect=(0, 0.08, 1, 1))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _interp_unique(degrees, values, angle):
    unique_deg, unique_indices = np.unique(degrees, return_index=True)
    unique_values = np.asarray(values)[unique_indices]
    return float(np.interp(angle % 360.0, unique_deg, unique_values))


def _plot_pv(
    path,
    label,
    degrees,
    volume_l,
    pressure_bar,
    extrema_angles,
    events,
    performance,
):
    plt = _setup_matplotlib()
    fig, ax = plt.subplots(figsize=(9.5, 7.2))

    ax.plot(volume_l, pressure_bar, label=f"{label} cylinder")

    for phase in (20, 65, 110, 155, 200, 245, 290, 335):
        x0 = _interp_unique(degrees, volume_l, phase)
        y0 = _interp_unique(degrees, pressure_bar, phase)
        x1 = _interp_unique(degrees, volume_l, phase + 7.0)
        y1 = _interp_unique(degrees, pressure_bar, phase + 7.0)
        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            arrowprops={"arrowstyle": "->", "lw": 1.0},
        )

    marks = []
    for name, angle in extrema_angles:
        marks.append((name, angle % 360.0))
    for event in events:
        short = (
            ("Hi-L " if event["valve"] == "H_i -> L" else "Ho-S ")
            + ("open" if event["transition"] == "opening" else "close")
        )
        marks.append((short, event["angle_degrees"]))

    dedup = []
    for name, angle in sorted(marks, key=lambda item: item[1]):
        if dedup and abs(angle - dedup[-1][1]) < 0.15:
            dedup[-1] = (dedup[-1][0] + "/" + name, dedup[-1][1])
        else:
            dedup.append((name, angle))

    for i, (name, angle) in enumerate(dedup):
        x = _interp_unique(degrees, volume_l, angle)
        y = _interp_unique(degrees, pressure_bar, angle)
        ax.plot(x, y, marker="o", markersize=4, linestyle="None")
        ax.annotate(
            f"{name}\ntheta={angle:.1f}°",
            xy=(x, y),
            xytext=(6, 7 if i % 2 == 0 else -24),
            textcoords="offset points",
            fontsize=7.5,
            ha="left",
        )

    ax.set_xlabel("Cylinder volume (L)")
    ax.set_ylabel("Pressure (bar absolute)")
    ax.grid(True, alpha=0.3)
    _legend_below(ax, ncol=1)
    ax.set_title(
        f"Stage F6 champion — {label} P-V loop\n"
        f"eta={100.0 * performance.thermal_efficiency:.5f}%   "
        f"P={performance.gas_power:.3f} W   "
        "arrows show increasing study-angle direction"
    )

    fig.tight_layout(rect=(0, 0.08, 1, 1))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _periodic_unique(degrees, values):
    """Return one unique 0..360 cycle and a 0..720 periodic extension."""
    degrees = np.asarray(degrees, dtype=float)
    values = np.asarray(values)

    unique_deg, unique_indices = np.unique(degrees, return_index=True)
    unique_values = values[..., unique_indices]

    if unique_deg[0] > 1e-12:
        raise RuntimeError("Cycle history does not start at theta=0.")
    if unique_deg[-1] < 360.0 - 1e-8:
        unique_deg = np.r_[unique_deg, 360.0]
        unique_values = np.concatenate(
            [unique_values, unique_values[..., :1]],
            axis=-1,
        )

    ext_deg = np.r_[unique_deg, unique_deg[1:] + 360.0]
    ext_values = np.concatenate(
        [unique_values, unique_values[..., 1:]],
        axis=-1,
    )
    return unique_deg, unique_values, ext_deg, ext_values


def _slice_with_interpolated_endpoints(x, y, start, end):
    mask = (x > start) & (x < end)
    xs = np.r_[start, x[mask], end]

    if y.ndim == 1:
        ys = np.r_[
            np.interp(start, x, y),
            y[mask],
            np.interp(end, x, y),
        ]
    else:
        channels = []
        for row in y:
            channels.append(
                np.r_[
                    np.interp(start, x, row),
                    row[mask],
                    np.interp(end, x, row),
                ]
            )
        ys = np.asarray(channels)

    return xs, ys


def _sector_mass_transfer(
    degrees,
    flows,
    masses,
    values,
    angular_speed,
):
    """Integrate signed/forward/reverse port traffic over the four motion sectors."""
    _, _, x, flow_ext = _periodic_unique(degrees, flows)
    _, _, x_mass, mass_ext = _periodic_unique(degrees, masses)

    smin = float(values["small_min_deg"])
    smax = float(values["small_max_deg"])
    lmax = float(values["large_max_deg"])
    while lmax <= smax:
        lmax += 360.0

    boundaries = (
        ("Smin_to_Lmin", smin, LARGE_MIN_DEG),
        ("Lmin_to_Smax", LARGE_MIN_DEG, smax),
        ("Smax_to_Lmax", smax, lmax),
        ("Lmax_to_next_Smin", lmax, smin + 360.0),
    )

    port_names = (
        "S_to_H_i",
        "H_i_to_L",
        "L_to_H_o",
        "H_o_to_S",
    )

    rows = []
    sectors = []

    for sector_name, start, end in boundaries:
        xs, fs = _slice_with_interpolated_endpoints(
            x,
            flow_ext,
            start,
            end,
        )
        dt_dtheta = 1.0 / angular_speed

        sector_record = {
            "sector": sector_name,
            "start_deg_unwrapped": start,
            "end_deg_unwrapped": end,
            "duration_deg": end - start,
            "ports": {},
        }

        for j, port_name in enumerate(port_names):
            flow = fs[j]
            forward = np.maximum(flow, 0.0)
            reverse = np.maximum(-flow, 0.0)

            scale = math.pi / 180.0 * dt_dtheta
            signed_kg = float(np.trapz(flow, xs) * scale)
            forward_kg = float(np.trapz(forward, xs) * scale)
            reverse_kg = float(np.trapz(reverse, xs) * scale)
            absolute_kg = forward_kg + reverse_kg

            record = {
                "forward_g": 1000.0 * forward_kg,
                "reverse_g": 1000.0 * reverse_kg,
                "signed_g": 1000.0 * signed_kg,
                "absolute_g": 1000.0 * absolute_kg,
            }
            sector_record["ports"][port_name] = record

            rows.append(
                {
                    "sector": sector_name,
                    "start_deg_unwrapped": start,
                    "end_deg_unwrapped": end,
                    "duration_deg": end - start,
                    "port": port_name,
                    **record,
                }
            )

        _, ms = _slice_with_interpolated_endpoints(
            x_mass,
            mass_ext,
            start,
            end,
        )
        sector_record["delta_H_i_storage_g"] = (
            1000.0 * float(ms[2, -1] - ms[2, 0])
        )
        sector_record["delta_H_o_storage_g"] = (
            1000.0 * float(ms[3, -1] - ms[3, 0])
        )
        sectors.append(sector_record)

    return rows, sectors


def main():
    best, values, wrapper, angles, trajectory, performance = _integrate_best()
    h = _histories(wrapper, angles, trajectory)

    deg = np.degrees(angles)
    events = _sampled_valve_events(angles, h["pressures"])

    volumes_l = h["volumes"] * 1000.0
    dvol_l_rad = h["dvol_dtheta"] * 1000.0
    pressure_bar = h["pressures"] / 1e5
    mass_g = h["masses"] * 1000.0
    flow_g_s = h["flows"] * 1000.0

    total_volume_l = volumes_l[0] + volumes_l[1]
    total_dvol_l_rad = dvol_l_rad[0] + dvol_l_rad[1]

    pressure_differences_bar = np.vstack(
        (
            pressure_bar[0] - pressure_bar[2],
            pressure_bar[2] - pressure_bar[1],
            pressure_bar[1] - pressure_bar[3],
            pressure_bar[3] - pressure_bar[0],
        )
    )

    temp_diff = np.vstack(
        (
            h["wall_temperature"][0] - h["temperatures"][2],
            h["wall_temperature"][1] - h["temperatures"][3],
            wrapper.heat_in.air_inlet_temperature_k - h["temperatures"][2],
            h["temperatures"][3] - wrapper.heat_out.air_inlet_temperature_k,
        )
    )

    title_prefix = (
        "Stage F6 champion "
        f"(index {int(best['index'])}, "
        f"eta={100.0 * performance.thermal_efficiency:.5f}%, "
        f"P={performance.gas_power:.3f} W)"
    )

    _save_single_plot(
        OUTPUTS["volumes"],
        title_prefix + " — cylinder volumes",
        "Volume (L)",
        deg,
        (
            (volumes_l[0], "S", None),
            (volumes_l[1], "L", None),
            (total_volume_l, "S + L", "--"),
        ),
        events,
        ncol=3,
        event_labels=True,
    )

    _save_single_plot(
        OUTPUTS["volume_derivatives"],
        title_prefix + " — angular volume derivatives",
        "dV/dtheta (L/rad)",
        deg,
        (
            (dvol_l_rad[0], "dV_S/dtheta", None),
            (dvol_l_rad[1], "dV_L/dtheta", None),
            (total_dvol_l_rad, "d(V_S+V_L)/dtheta", "--"),
        ),
        events,
        zero_line=True,
        ncol=3,
        event_labels=True,
    )

    _save_single_plot(
        OUTPUTS["pressures"],
        title_prefix + " — control-volume pressures",
        "Pressure (bar absolute)",
        deg,
        (
            (pressure_bar[0], "S", None),
            (pressure_bar[1], "L", None),
            (pressure_bar[2], "H_i / C", None),
            (pressure_bar[3], "H_o / H", None),
        ),
        events,
        ncol=4,
        event_labels=True,
    )

    _save_single_plot(
        OUTPUTS["pressure_differences"],
        title_prefix + " — hydraulic pressure differences",
        "Pressure difference (bar)",
        deg,
        (
            (pressure_differences_bar[0], "S - H_i", None),
            (pressure_differences_bar[1], "H_i - L (valve)", None),
            (pressure_differences_bar[2], "L - H_o", None),
            (pressure_differences_bar[3], "H_o - S (valve)", None),
        ),
        events,
        zero_line=True,
        ncol=4,
        event_labels=True,
    )

    _save_single_plot(
        OUTPUTS["temperatures"],
        title_prefix + " — gas and wall temperatures",
        "Temperature (K)",
        deg,
        (
            (h["temperatures"][0], "S gas", None),
            (h["temperatures"][1], "L gas", None),
            (h["temperatures"][2], "H_i gas", None),
            (h["temperatures"][3], "H_o gas", None),
            (h["wall_temperature"][0], "H_i wall", "--"),
            (h["wall_temperature"][1], "H_o wall", "--"),
            (
                np.full_like(deg, wrapper.heat_in.air_inlet_temperature_k),
                "Hot external-air inlet",
                ":",
            ),
            (
                np.full_like(deg, wrapper.heat_out.air_inlet_temperature_k),
                "Cold external-air inlet",
                ":",
            ),
        ),
        events,
        ncol=4,
        event_labels=True,
    )

    _save_single_plot(
        OUTPUTS["temperature_differences"],
        title_prefix + " — exchanger temperature differences",
        "Temperature difference (K)",
        deg,
        (
            (temp_diff[0], "H_i wall - H_i gas", None),
            (temp_diff[1], "H_o wall - H_o gas", None),
            (temp_diff[2], "Hot inlet air - H_i gas", "--"),
            (temp_diff[3], "H_o gas - cold inlet air", "--"),
        ),
        events,
        zero_line=True,
        ncol=2,
        event_labels=True,
    )

    _save_single_plot(
        OUTPUTS["mass_flows"],
        title_prefix + " — actual simulated port mass flows",
        "Signed mass flow (g/s)",
        deg,
        (
            (flow_g_s[0], "S -> H_i", None),
            (flow_g_s[1], "H_i -> L valve", None),
            (flow_g_s[2], "L -> H_o", None),
            (flow_g_s[3], "H_o -> S valve", None),
        ),
        events,
        zero_line=True,
        ncol=4,
        event_labels=True,
    )

    _save_single_plot(
        OUTPUTS["masses"],
        title_prefix + " — gas mass in each control volume",
        "Gas mass (g)",
        deg,
        (
            (mass_g[0], "S", None),
            (mass_g[1], "L", None),
            (mass_g[2], "H_i / C", None),
            (mass_g[3], "H_o / H", None),
        ),
        events,
        ncol=4,
        event_labels=True,
    )

    _save_single_plot(
        OUTPUTS["heat_rates"],
        title_prefix + " — instantaneous heat and work rates",
        "Rate (W)",
        deg,
        (
            (h["gas_heat"][0], "Heat to H_i gas", None),
            (h["gas_heat"][1], "Heat to H_o gas", None),
            (h["air_heat"][0], "External air -> H_i wall", "--"),
            (h["air_heat"][1], "External air -> H_o wall", "--"),
            (h["work_rate"], "Gas boundary-work rate", ":"),
        ),
        events,
        zero_line=True,
        ncol=3,
        event_labels=True,
    )

    _save_single_plot(
        OUTPUTS["valves"],
        title_prefix + " — passive ideal-diode states",
        "Valve state",
        deg,
        (
            (h["valve_open"][0].astype(float), "H_i -> L", None),
            (h["valve_open"][1].astype(float), "H_o -> S", None),
        ),
        events,
        ncol=2,
        event_labels=True,
    )

    extrema_angles = (
        ("Smin", float(values["small_min_deg"]) % 360.0),
        ("Lmin", LARGE_MIN_DEG % 360.0),
        ("Smax", float(values["small_max_deg"]) % 360.0),
        ("Lmax", float(values["large_max_deg"]) % 360.0),
    )

    _plot_pv(
        OUTPUTS["pv_s"],
        "S",
        deg,
        volumes_l[0],
        pressure_bar[0],
        extrema_angles,
        events,
        performance,
    )
    _plot_pv(
        OUTPUTS["pv_l"],
        "L",
        deg,
        volumes_l[1],
        pressure_bar[1],
        extrema_angles,
        events,
        performance,
    )

    sector_rows, sector_summary = _sector_mass_transfer(
        deg,
        h["flows"],
        h["masses"],
        values,
        wrapper.model.angular_speed,
    )

    OUTPUTS["sector_csv"].parent.mkdir(parents=True, exist_ok=True)
    with OUTPUTS["sector_csv"].open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=(
                "sector",
                "start_deg_unwrapped",
                "end_deg_unwrapped",
                "duration_deg",
                "port",
                "forward_g",
                "reverse_g",
                "signed_g",
                "absolute_g",
            ),
        )
        writer.writeheader()
        writer.writerows(sector_rows)

    np.savez_compressed(
        OUTPUTS["trajectory"],
        angles=angles,
        trajectory=trajectory,
        volumes=h["volumes"],
        dvol_dtheta=h["dvol_dtheta"],
        pressures=h["pressures"],
        temperatures=h["temperatures"],
        masses=h["masses"],
        flows=h["flows"],
        wall_temperature=h["wall_temperature"],
        gas_heat=h["gas_heat"],
        air_heat=h["air_heat"],
        valve_open=h["valve_open"],
    )

    summary = {
        "candidate": {
            "index": int(best["index"]),
            "parameters": values,
            "recorded_efficiency": float(
                best["indicated_thermal_efficiency"]
            ),
            "recorded_indicated_power_w": float(best["indicated_power_w"]),
        },
        "reevaluated": {
            "thermal_efficiency": performance.thermal_efficiency,
            "indicated_power_w": performance.gas_power,
            "heat_input_w": performance.heat_in_power,
            "heat_out_w": performance.heat_out_power,
        },
        "valve_events": events,
        "motion_sectors": sector_summary,
        "interpretation_note": (
            "Sector port integrals report actual simulated port traffic. "
            "Forward/reverse values on the bidirectional S-H_i and L-H_o "
            "links quantify reflux explicitly. Exchanger-control-volume mass "
            "changes show storage effects. These values are more informative "
            "than V_S+V_L alone for deciding whether gas moves through the "
            "exchangers during the short common-compression/common-expansion "
            "sectors."
        ),
    }
    OUTPUTS["summary_json"].write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    print("Stage F6 champion thermodynamic diagnostic")
    print(
        f"  index = {int(best['index'])}\n"
        f"  eta   = {100.0 * performance.thermal_efficiency:.6f}%\n"
        f"  P     = {performance.gas_power:.6f} W\n"
        f"  Qin   = {performance.heat_in_power:.6f} W"
    )

    print("\nSampled passive-valve events")
    for event in events:
        print(
            f"  {event['angle_degrees']:8.3f} deg  "
            f"{event['valve']:9s}  {event['transition']}"
        )

    print("\nMass transfer by motion sector (g)")
    for sector in sector_summary:
        print(
            f"\n  {sector['sector']}  "
            f"{sector['start_deg_unwrapped']:.3f} -> "
            f"{sector['end_deg_unwrapped']:.3f} deg"
        )
        for port, record in sector["ports"].items():
            print(
                f"    {port:10s} "
                f"forward={record['forward_g']:9.5f}  "
                f"reverse={record['reverse_g']:9.5f}  "
                f"signed={record['signed_g']:9.5f}"
            )
        print(
            f"    delta H_i storage = "
            f"{sector['delta_H_i_storage_g']:+.6f} g"
        )
        print(
            f"    delta H_o storage = "
            f"{sector['delta_H_o_storage_g']:+.6f} g"
        )

    print("\nWrote:")
    for path in OUTPUTS.values():
        print(f"  {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
