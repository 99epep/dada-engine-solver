#!/usr/bin/env python3
"""Detailed LP/HP exchange diagnostic for a Fourier motor champion.

The report is a snapshot: this script is safe to run on an in-progress 8H
campaign.  It compares the current Fourier champion with the retained linear-UU
motion rebuilt on EXACTLY THE SAME current hardware/geometry, so the plotted
differences are predominantly kinematic.

Examples
--------
6 harmonics:
    PYTHONPATH=src python3 examples/diagnose_motor_fourier_exchanges.py \
        --report outputs/motor_fourier_c2_6h_refine/report.json

8 harmonics, while the optimization is still running:
    PYTHONPATH=src python3 examples/diagnose_motor_fourier_exchanges.py \
        --report outputs/motor_fourier_c2_8h_refine/report.json

Outputs are written beside the supplied report:
    exchange_diagnostic_<Nh>_full.png
    exchange_diagnostic_<Nh>_BP.png
    exchange_diagnostic_<Nh>_HP.png
    exchange_diagnostic_<Nh>_history.csv
    exchange_diagnostic_<Nh>_summary.json

Angle convention is the same forward motor angle used by the Fourier comparison
plots: 0..360 degrees, with the retained linear phases
    LP exchange -> compression -> HP exchange -> expansion.
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import argparse
import csv
import json
import math
import os
import tempfile
import time

import numpy as np

from dada_solver.configuration import ChargeConfiguration
from dada_solver.dynamics import ValveTopology
from dada_solver.exchangers.air_wall import AirWallMotor
from dada_solver.exchangers.wall_cycle import (
    solve_periodic_wall_motor,
    wall_cycle_performance,
)
from dada_solver.state import ThermodynamicState
from dada_solver.valves import ValveState

from compare_motor_motion_laws_stage7A5 import ROOT, _uniform_wall_initial
from optimize_motor_temperature_point import (
    base_geometry,
    build_design as linear_build_design,
)
from refine_motor_four_stage_hx9d_variable_gas import _load_basis
from refine_motor_fourier_c2_260k import HOT_K, build_design as fourier_build_design


DEFAULT_REPORT = (
    ROOT / "outputs" / "motor_fourier_c2_6h_refine" / "report.json"
)


def _read_json_snapshot(path: Path, attempts: int = 30) -> dict:
    """Tolerate catching report.json during its atomic-ish rewrite."""
    last_error = None
    for _ in range(attempts):
        try:
            return json.loads(path.read_text())
        except (json.JSONDecodeError, OSError) as exc:
            last_error = exc
            time.sleep(0.1)
    raise RuntimeError(f"Could not read a stable JSON snapshot from {path}: {last_error}")


def _setup_matplotlib():
    cache = Path(tempfile.gettempdir()) / "dada_solver_matplotlib"
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _linear_parameters(current_thermo: dict, reference: dict) -> dict:
    p = dict(current_thermo)
    p.update(reference["parameters"])
    return p


def _build_pair(report: dict):
    best = report.get("best_feasible")
    if best is None:
        raise RuntimeError("Report contains no best_feasible candidate.")
    if best.get("last_complete_state") is None:
        raise RuntimeError("Current champion has no saved periodic state.")
    reference = report.get("reference_linear_UU")
    if reference is None or reference.get("parameters") is None:
        raise RuntimeError("Report contains no reference_linear_UU parameters.")

    h = int(best.get("harmonics", len(best["small_coefficients"]) // 2))
    if len(best["small_coefficients"]) != 2 * h:
        raise RuntimeError("Small-piston coefficient count does not match harmonics.")
    if len(best["large_coefficients"]) != 2 * h:
        raise RuntimeError("Large-piston coefficient count does not match harmonics.")

    campaign_definition, base = _load_basis()
    geometry = base_geometry(base)

    fourier_design = fourier_build_design(
        base,
        geometry,
        best["thermo"],
        np.asarray(best["small_coefficients"], dtype=float),
        np.asarray(best["large_coefficients"], dtype=float),
        h,
    )

    lp = _linear_parameters(best["thermo"], reference)
    linear_design = linear_build_design(base, geometry, lp, HOT_K)
    linear_design = replace(
        linear_design,
        configuration=replace(
            linear_design.configuration,
            heat_in_valve_placement="upstream",
            heat_out_valve_placement="upstream",
        ),
    )

    return best, reference, h, campaign_definition, fourier_design, linear_design



def _fixed_inventory_uniform_initial(design, total_mass_kg, wall_state_source=None):
    """Create a safe state for a new motion law at fixed total gas inventory.

    Gas is initialized uniformly at the configured charge temperature using the
    *target* motion law's theta=0 volumes, eliminating the artificial pressure
    jumps that occur when a periodic state from another kinematics is copied
    directly.  When supplied, wall energies are retained from the source
    periodic state because hardware is identical in this diagnostic.
    """
    config = replace(
        design.configuration,
        charge=ChargeConfiguration(
            temperature=design.configuration.charge.temperature,
            total_mass=float(total_mass_kg),
        ),
    )
    design = replace(design, configuration=config)
    wrapper = design.build()
    if not isinstance(wrapper, AirWallMotor):
        raise TypeError("Expected AirWallMotor while building fixed-inventory initial state.")

    initial = np.asarray(
        _uniform_wall_initial(design.configuration, wrapper),
        dtype=float,
    )

    if wall_state_source is not None:
        source = np.asarray(wall_state_source, dtype=float)
        if source.shape != (10,):
            raise ValueError("Expected ten-state wall source.")
        # Same exchanger hardware => same wall capacities, so wall energies can
        # be copied directly.  Gas states are intentionally NOT copied.
        initial[8:10] = source[8:10]

    mass = float(np.sum(initial[:8:2]))
    if not math.isclose(mass, total_mass_kg, rel_tol=1e-10, abs_tol=1e-14):
        raise RuntimeError(
            f"Fixed-inventory initialization mismatch: {mass} vs {total_mass_kg} kg"
        )

    return design, initial


def _solve(design, settings, initial, label: str):
    wrapper = design.build()
    if not isinstance(wrapper, AirWallMotor):
        raise TypeError(f"{label}: expected AirWallMotor.")

    periodic = solve_periodic_wall_motor(
        wrapper,
        np.asarray(initial, dtype=float),
        maximum_cycles=design.configuration.numerical.maximum_cycles,
        settings=settings,
    )
    if not periodic.converged:
        raise RuntimeError(f"{label}: periodic solve failed: {periodic.message}")
    if periodic.angles is None or periodic.trajectory is None:
        raise RuntimeError(f"{label}: periodic solve returned no trajectory.")

    performance = wall_cycle_performance(wrapper, periodic.trajectory)
    return wrapper, periodic.angles, periodic.trajectory, performance, periodic.last_complete_state


def _extract_raw(wrapper, angles, trajectory):
    n = len(angles)

    q = np.empty((2, n))
    dq_cycle = np.empty((2, n))
    pressure = np.empty((4, n))
    temperature = np.empty((4, n))
    mass = np.empty((4, n))
    flow = np.empty((4, n))
    gas_heat = np.empty((2, n))
    air_heat = np.empty((2, n))
    wall_temperature = np.empty((2, n))
    valve_open = np.empty((2, n), dtype=float)

    correlation = np.empty((4, n), dtype=object)
    reynolds = np.full((4, n), np.nan)
    residence = np.full((4, n), np.nan)

    limits_s = wrapper.model.kinematics.small_volume_limits
    limits_l = wrapper.model.kinematics.large_volume_limits

    closed = ValveTopology(ValveState.CLOSED, ValveState.CLOSED)

    for i, angle in enumerate(angles):
        a = float(angle)
        values = trajectory[:, i]
        state = ThermodynamicState.from_array(values[:8])

        point = wrapper.model.instantaneous_point(a, state, closed)
        volumes = point.volumes
        p = point.pressures
        t = point.temperatures
        f = point.flows

        pressure[:, i] = p
        temperature[:, i] = t
        mass[:, i] = values[[0, 2, 4, 6]]

        q[0, i] = (volumes.small_cylinder - limits_s.minimum) / limits_s.swept
        q[1, i] = (volumes.large_cylinder - limits_l.minimum) / limits_l.swept

        dvs = wrapper.model.kinematics.small_cylinder_volume_derivative(a)
        dvl = wrapper.model.kinematics.large_cylinder_volume_derivative(a)
        # motor_fraction = -theta/(2*pi), hence dq/d(motor_fraction)
        dq_cycle[0, i] = -2.0 * math.pi * dvs / limits_s.swept
        dq_cycle[1, i] = -2.0 * math.pi * dvl / limits_l.swept

        flow[:, i] = (
            f.small_to_cold,   # S -> H_i
            f.cold_to_large,   # H_i -> L
            f.large_to_hot,    # L -> H_o
            f.hot_to_small,    # H_o -> S
        )

        valve_open[0, i] = float(p[2] > p[1])  # H_i -> L
        valve_open[1, i] = float(p[3] > p[0])  # H_o -> S

        contexts = wrapper.flow_contexts(a, values)
        hi = wrapper.heat_in.thermal_point(
            t[2], values[8], context=contexts[0]
        )
        ho = wrapper.heat_out.thermal_point(
            t[3], values[9], context=contexts[1]
        )

        gas_heat[:, i] = (hi.gas_heat_w, ho.gas_heat_w)
        air_heat[:, i] = (hi.air_heat_w, ho.air_heat_w)
        wall_temperature[:, i] = (hi.wall_temperature_k, ho.wall_temperature_k)

        diag_rows = tuple(hi.film_diagnostics) + tuple(ho.film_diagnostics)
        if len(diag_rows) != 4:
            # Static-UA fallback, though current optimized cases use variable film.
            correlation[:, i] = "static_UA"
        else:
            for j, diag in enumerate(diag_rows):
                correlation[j, i] = str(diag.correlation_id)
                reynolds[j, i] = float(diag.reynolds)
                if diag.residence_time_s is not None:
                    residence[j, i] = float(diag.residence_time_s)

    return {
        "q": q,
        "dq_cycle": dq_cycle,
        "pressure": pressure,
        "temperature": temperature,
        "mass": mass,
        "flow": flow,
        "gas_heat": gas_heat,
        "air_heat": air_heat,
        "wall_temperature": wall_temperature,
        "valve_open": valve_open,
        "correlation": correlation,
        "reynolds": reynolds,
        "residence": residence,
    }


def _forward_cycle(angles, raw):
    """Convert solver theta direction to the forward motor-angle convention."""
    motor_deg = (-np.degrees(np.asarray(angles))) % 360.0

    # The integrated endpoint duplicates 0 degrees. Round only for deduplication.
    rounded = np.round(motor_deg, 10)
    unique_deg, indices = np.unique(rounded, return_index=True)

    order = np.argsort(unique_deg)
    deg = unique_deg[order]
    indices = indices[order]

    out = {}
    for key, values in raw.items():
        a = np.asarray(values)
        if a.ndim == 1:
            selected = a[indices]
            selected = np.r_[selected, selected[0]]
        else:
            selected = a[..., indices]
            selected = np.concatenate((selected, selected[..., :1]), axis=-1)
        out[key] = selected

    deg = np.r_[deg, 360.0]
    return deg, out


def _linear_motion(deg, p):
    t = np.asarray(deg) / 360.0
    knots = np.asarray([0.0, p["t1"], p["t2"], p["t3"], 1.0])
    sv = np.asarray([p["b_s"], 1.0, p["a_s"], 0.0, p["b_s"]])
    lv = np.asarray([1.0, p["b_l"], 0.0, p["a_l"], 1.0])

    qs = np.interp(t, knots, sv)
    ql = np.interp(t, knots, lv)

    dqs = np.empty_like(t)
    dql = np.empty_like(t)
    for i in range(4):
        mask = (t >= knots[i]) & (t <= knots[i + 1])
        dqs[mask] = (sv[i + 1] - sv[i]) / (knots[i + 1] - knots[i])
        dql[mask] = (lv[i + 1] - lv[i]) / (knots[i + 1] - knots[i])

    return np.vstack((qs, ql)), np.vstack((dqs, dql))


def _interp_history(deg_source, values, deg_target):
    values = np.asarray(values)
    if values.dtype == object:
        # nearest-neighbour for categorical data
        idx = np.searchsorted(deg_source, deg_target, side="left")
        idx = np.clip(idx, 0, len(deg_source) - 1)
        prev = np.clip(idx - 1, 0, len(deg_source) - 1)
        choose_prev = (
            np.abs(deg_target - deg_source[prev])
            <= np.abs(deg_target - deg_source[idx])
        )
        idx = np.where(choose_prev, prev, idx)
        return values[..., idx]

    if values.ndim == 1:
        return np.interp(deg_target, deg_source, values)

    return np.vstack(
        [np.interp(deg_target, deg_source, row) for row in values]
    )


def _common_grid(fdeg, fh, ldeg, lh):
    deg = np.linspace(0.0, 360.0, 2881)
    f = {k: _interp_history(fdeg, v, deg) for k, v in fh.items()}
    l = {k: _interp_history(ldeg, v, deg) for k, v in lh.items()}
    return deg, f, l


def _phase_markers(reference):
    p = reference["parameters"]
    return {
        "LP_end": 360.0 * p["t1"],
        "compression_end": 360.0 * p["t2"],
        "HP_end": 360.0 * p["t3"],
    }


def _candidate_extrema(deg, q):
    return {
        "S_max": float(deg[np.argmax(q[0])]),
        "S_min": float(deg[np.argmin(q[0])]),
        "L_max": float(deg[np.argmax(q[1])]),
        "L_min": float(deg[np.argmin(q[1])]),
    }


def _periodic_window(deg, values, lo, hi):
    """Slice numeric/categorical histories on an optionally wrapped window."""
    if lo >= 0.0 and hi <= 360.0:
        mask = (deg >= lo) & (deg <= hi)
        return deg[mask], values[..., mask]

    x = np.r_[deg[:-1] - 360.0, deg, deg[1:] + 360.0]
    if values.ndim == 1:
        y = np.r_[values[:-1], values, values[1:]]
    else:
        y = np.concatenate(
            (values[..., :-1], values, values[..., 1:]), axis=-1
        )
    mask = (x >= lo) & (x <= hi)
    return x[mask], y[..., mask]


def _regime_codes(*arrays):
    names = sorted(
        {
            str(v)
            for array in arrays
            for v in np.asarray(array, dtype=object).ravel()
        }
    )
    # Put the two known regimes in an intuitive order when present.
    preferred = [
        "stagnant_radial_screening",
        "hausen_constant_wall",
        "static_UA",
    ]
    ordered = [x for x in preferred if x in names]
    ordered += [x for x in names if x not in ordered]
    code = {name: i for i, name in enumerate(ordered)}
    return code


def _decorate(ax, xlim, markers, extrema, *, top=False):
    ax.set_xlim(*xlim)
    ticks = np.arange(
        math.ceil(xlim[0] / 15.0) * 15.0,
        math.floor(xlim[1] / 15.0) * 15.0 + 0.1,
        15.0,
    )
    ax.set_xticks(ticks)
    ax.grid(True, alpha=0.25)

    for name, angle in markers.items():
        candidates = [angle - 360.0, angle, angle + 360.0]
        for x in candidates:
            if xlim[0] <= x <= xlim[1]:
                ax.axvline(x, linestyle=":", linewidth=0.9, alpha=0.55)
                if top:
                    ax.annotate(
                        name.replace("_", " "),
                        xy=(x, 1.0),
                        xycoords=("data", "axes fraction"),
                        xytext=(2, -2),
                        textcoords="offset points",
                        rotation=90,
                        va="top",
                        fontsize=7,
                        alpha=0.75,
                    )

    for name, angle in extrema.items():
        candidates = [angle - 360.0, angle, angle + 360.0]
        for x in candidates:
            if xlim[0] <= x <= xlim[1]:
                ax.axvline(x, linestyle="--", linewidth=0.7, alpha=0.35)
                if top:
                    ax.annotate(
                        name,
                        xy=(x, 0.0),
                        xycoords=("data", "axes fraction"),
                        xytext=(2, 2),
                        textcoords="offset points",
                        rotation=90,
                        va="bottom",
                        fontsize=7,
                        alpha=0.65,
                    )


def _phase_labels(ax, reference, xlim):
    p = reference["parameters"]
    spans = [
        (0.0, 360.0 * p["t1"], "LP exchange"),
        (360.0 * p["t1"], 360.0 * p["t2"], "compression"),
        (360.0 * p["t2"], 360.0 * p["t3"], "HP exchange"),
        (360.0 * p["t3"], 360.0, "expansion"),
    ]
    for a, b, label in spans:
        for shift in (-360.0, 0.0, 360.0):
            aa, bb = a + shift, b + shift
            left, right = max(aa, xlim[0]), min(bb, xlim[1])
            if right > left:
                ax.text(
                    0.5 * (left + right),
                    1.03,
                    label,
                    transform=ax.get_xaxis_transform(),
                    ha="center",
                    va="bottom",
                    fontsize=8,
                )


def _legend(ax, ncol=4):
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.16),
        ncol=ncol,
        frameon=False,
        fontsize=8,
    )


def _plot_comparison(
    path,
    title,
    deg,
    fourier,
    linear,
    reference,
    extrema,
    xlim,
):
    plt = _setup_matplotlib()
    fig, axes = plt.subplots(
        8, 1, figsize=(14, 23), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.0, 1.15, 1.0, 1.1, 1.0, 0.9, 0.72]},
    )

    # Wrapped window extraction for each numeric panel.
    def window(v):
        return _periodic_window(deg, v, xlim[0], xlim[1])

    # 1 — normalized motion
    x, fq = window(fourier["q"])
    _, lq = window(linear["q"])
    axes[0].plot(x, fq[0], label="S Fourier")
    axes[0].plot(x, fq[1], label="L Fourier")
    axes[0].plot(x, lq[0], "--", label="S linear, same hardware")
    axes[0].plot(x, lq[1], "--", label="L linear, same hardware")
    axes[0].set_ylabel("Normalized\nvolume")
    _legend(axes[0], 4)

    # 2 — piston speed in normalized stroke per cycle
    x, fdq = window(fourier["dq_cycle"])
    _, ldq = window(linear["dq_cycle"])
    axes[1].plot(x, fdq[0], label="dS/dcycle Fourier")
    axes[1].plot(x, fdq[1], label="dL/dcycle Fourier")
    axes[1].plot(x, ldq[0], "--", label="dS/dcycle linear")
    axes[1].plot(x, ldq[1], "--", label="dL/dcycle linear")
    axes[1].axhline(0.0, linewidth=0.8, alpha=0.5)
    axes[1].set_ylabel("Normalized\npiston speed")
    _legend(axes[1], 4)

    # 3 — actual gas flows
    x, ff = window(1000.0 * fourier["flow"])
    _, lf = window(1000.0 * linear["flow"])
    flow_names = ("S→H_i", "H_i→L", "L→H_o", "H_o→S")
    for i, name in enumerate(flow_names):
        axes[2].plot(x, ff[i], label=f"{name} Fourier")
        axes[2].plot(x, lf[i], "--", label=f"{name} linear")
    axes[2].axhline(0.0, linewidth=0.8, alpha=0.5)
    axes[2].set_ylabel("Actual flow\n(g/s)")
    _legend(axes[2], 4)

    # 4 — cylinder pressures
    x, fp = window(fourier["pressure"] / 1e5)
    _, lp = window(linear["pressure"] / 1e5)
    axes[3].plot(x, fp[0], label="P_S Fourier")
    axes[3].plot(x, fp[1], label="P_L Fourier")
    axes[3].plot(x, lp[0], "--", label="P_S linear")
    axes[3].plot(x, lp[1], "--", label="P_L linear")
    axes[3].set_ylabel("Cylinder P\n(bar abs)")
    _legend(axes[3], 4)

    # 5 — exchanger gas/wall temperatures
    x, ft = window(fourier["temperature"])
    _, lt = window(linear["temperature"])
    _, fw = window(fourier["wall_temperature"])
    _, lw = window(linear["wall_temperature"])
    axes[4].plot(x, ft[2], label="H_i gas Fourier")
    axes[4].plot(x, ft[3], label="H_o gas Fourier")
    axes[4].plot(x, fw[0], label="H_i wall Fourier")
    axes[4].plot(x, fw[1], label="H_o wall Fourier")
    axes[4].plot(x, lt[2], "--", label="H_i gas linear")
    axes[4].plot(x, lt[3], "--", label="H_o gas linear")
    axes[4].plot(x, lw[0], "--", label="H_i wall linear")
    axes[4].plot(x, lw[1], "--", label="H_o wall linear")
    axes[4].set_ylabel("Exchanger T\n(K)")
    _legend(axes[4], 4)

    # 6 — instantaneous wall-to-gas heat
    x, fqh = window(fourier["gas_heat"])
    _, lqh = window(linear["gas_heat"])
    axes[5].plot(x, fqh[0], label="H_i gas heat Fourier")
    axes[5].plot(x, fqh[1], label="H_o gas heat Fourier")
    axes[5].plot(x, lqh[0], "--", label="H_i gas heat linear")
    axes[5].plot(x, lqh[1], "--", label="H_o gas heat linear")
    axes[5].axhline(0.0, linewidth=0.8, alpha=0.5)
    axes[5].set_ylabel("Wall→gas heat\n(W)")
    _legend(axes[5], 4)

    # 7 — exchanger gas storage
    x, fm = window(1000.0 * fourier["mass"])
    _, lm = window(1000.0 * linear["mass"])
    axes[6].plot(x, fm[2], label="m_Hi Fourier")
    axes[6].plot(x, fm[3], label="m_Ho Fourier")
    axes[6].plot(x, lm[2], "--", label="m_Hi linear")
    axes[6].plot(x, lm[3], "--", label="m_Ho linear")
    axes[6].set_ylabel("Exchanger\nmass (g)")
    _legend(axes[6], 4)

    # 8 — correlation regime map
    code = _regime_codes(fourier["correlation"], linear["correlation"])
    x, fc = window(fourier["correlation"])
    _, lc = window(linear["correlation"])
    rows = np.vstack(
        [
            np.vectorize(code.get)(fc[i]) for i in range(4)
        ] + [
            np.vectorize(code.get)(lc[i]) for i in range(4)
        ]
    ).astype(float)

    # pcolormesh expects cell edges; the minor x difference is irrelevant here.
    if len(x) > 1:
        dx = float(np.median(np.diff(x)))
    else:
        dx = 1.0
    xedges = np.r_[x - 0.5 * dx, x[-1] + 0.5 * dx]
    yedges = np.arange(9) - 0.5
    mesh = axes[7].pcolormesh(xedges, yedges, rows, shading="flat")
    axes[7].set_yticks(
        np.arange(8),
        (
            "F H_i in", "F H_i out", "F H_o in", "F H_o out",
            "L H_i in", "L H_i out", "L H_o in", "L H_o out",
        ),
    )
    axes[7].set_ylabel("Film regime")
    cb = fig.colorbar(mesh, ax=axes[7], pad=0.01)
    cb.set_ticks(list(code.values()))
    cb.set_ticklabels(list(code.keys()))

    markers = _phase_markers(reference)
    for i, ax in enumerate(axes):
        _decorate(ax, xlim, markers, extrema, top=(i == 0))
    _phase_labels(axes[0], reference, xlim)

    axes[-1].set_xlabel("Forward motor angle (deg)")
    fig.suptitle(title)
    fig.tight_layout(rect=(0, 0.01, 1, 0.985))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _phase_stats(deg, h, start, end, omega):
    mask = (deg >= start) & (deg <= end)
    x = deg[mask]
    if len(x) < 2:
        return {}

    # dt = d(angle_deg) * pi/(180*omega)
    time_scale = math.pi / (180.0 * omega)

    def integral(y):
        return float(np.trapz(np.asarray(y)[mask], x) * time_scale)

    duration_s = float(end - start) * time_scale
    result = {
        "start_deg": float(start),
        "end_deg": float(end),
        "duration_deg": float(end - start),
        "duration_s": duration_s,
        "ports": {},
    }
    names = ("S_to_Hi", "Hi_to_L", "L_to_Ho", "Ho_to_S")
    for i, name in enumerate(names):
        y = h["flow"][i]
        absolute_mass = integral(np.abs(y))
        result["ports"][name] = {
            "signed_mass_g": 1000.0 * integral(y),
            "absolute_mass_g": 1000.0 * absolute_mass,
            "mean_abs_flow_g_s": 1000.0 * absolute_mass / duration_s,
            "peak_abs_flow_g_s": 1000.0 * float(np.max(np.abs(y[mask]))),
        }

    result["heat"] = {
        "Hi_gas_heat_J": integral(h["gas_heat"][0]),
        "Ho_gas_heat_J": integral(h["gas_heat"][1]),
    }
    result["storage_change_g"] = {
        "Hi": 1000.0 * float(h["mass"][2, mask][-1] - h["mass"][2, mask][0]),
        "Ho": 1000.0 * float(h["mass"][3, mask][-1] - h["mass"][3, mask][0]),
    }

    result["correlation_time_fraction"] = {}
    labels = ("Hi.inlet", "Hi.outlet", "Ho.inlet", "Ho.outlet")
    for i, name in enumerate(labels):
        all_values = h["correlation"][i]
        values = all_values[mask]
        fractions = {}
        for regime in np.unique(values):
            indicator = (all_values == regime).astype(float)
            fractions[str(regime)] = integral(indicator) / duration_s
        result["correlation_time_fraction"][name] = fractions

    return result


def _write_csv(path, deg, fourier, linear):
    headers = ["motor_angle_deg"]
    series = []

    def add(prefix, data):
        mapping = [
            ("qS", data["q"][0]), ("qL", data["q"][1]),
            ("dqS_dcycle", data["dq_cycle"][0]), ("dqL_dcycle", data["dq_cycle"][1]),
            ("flow_S_Hi_kg_s", data["flow"][0]), ("flow_Hi_L_kg_s", data["flow"][1]),
            ("flow_L_Ho_kg_s", data["flow"][2]), ("flow_Ho_S_kg_s", data["flow"][3]),
            ("P_S_Pa", data["pressure"][0]), ("P_L_Pa", data["pressure"][1]),
            ("P_Hi_Pa", data["pressure"][2]), ("P_Ho_Pa", data["pressure"][3]),
            ("T_S_K", data["temperature"][0]), ("T_L_K", data["temperature"][1]),
            ("T_Hi_K", data["temperature"][2]), ("T_Ho_K", data["temperature"][3]),
            ("Twall_Hi_K", data["wall_temperature"][0]),
            ("Twall_Ho_K", data["wall_temperature"][1]),
            ("Qgas_Hi_W", data["gas_heat"][0]), ("Qgas_Ho_W", data["gas_heat"][1]),
            ("mass_Hi_kg", data["mass"][2]), ("mass_Ho_kg", data["mass"][3]),
            ("Re_Hi_in", data["reynolds"][0]), ("Re_Hi_out", data["reynolds"][1]),
            ("Re_Ho_in", data["reynolds"][2]), ("Re_Ho_out", data["reynolds"][3]),
        ]
        for name, values in mapping:
            headers.append(f"{prefix}_{name}")
            series.append(np.asarray(values))

        for i, name in enumerate(("corr_Hi_in", "corr_Hi_out", "corr_Ho_in", "corr_Ho_out")):
            headers.append(f"{prefix}_{name}")
            series.append(np.asarray(data["correlation"][i], dtype=object))

    add("fourier", fourier)
    add("linear_same_hw", linear)

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(headers)
        for i, angle in enumerate(deg):
            writer.writerow(
                [f"{angle:.9f}"]
                + [
                    str(values[i]) if values.dtype == object else f"{float(values[i]):.12g}"
                    for values in series
                ]
            )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    ap.add_argument(
        "--skip-linear",
        action="store_true",
        help="Only re-integrate the Fourier champion. Faster, but loses the same-hardware comparison.",
    )
    args = ap.parse_args()

    report_path = args.report if args.report.is_absolute() else ROOT / args.report
    if not report_path.exists():
        raise FileNotFoundError(report_path)

    report = _read_json_snapshot(report_path)
    best, reference, h, campaign_definition, fourier_design, linear_design = _build_pair(report)

    print(
        f"Snapshot champion: {h}H index={best.get('index')} "
        f"eta={100.0*best['result']['indicated_thermal_efficiency']:.6f}%"
    )

    saved = np.asarray(best["last_complete_state"], dtype=float)

    print("Re-integrating Fourier champion from its saved periodic state...")
    fw, fa, ft, fp, fstate = _solve(
        fourier_design,
        campaign_definition.wall_numerical_settings,
        saved,
        "Fourier",
    )
    fdeg, fh = _forward_cycle(fa, _extract_raw(fw, fa, ft))

    if args.skip_linear:
        # Keep the plotting/data schema useful: the dashed baseline becomes the
        # analytical retained linear motion for q/dq and NaN elsewhere.
        deg = np.linspace(0.0, 360.0, 2881)
        fh = {k: _interp_history(fdeg, v, deg) for k, v in fh.items()}
        qlin, dqlin = _linear_motion(deg, reference["parameters"])
        lh = {}
        for k, v in fh.items():
            if np.asarray(v).dtype == object:
                lh[k] = np.full_like(v, "not_evaluated", dtype=object)
            else:
                lh[k] = np.full_like(v, np.nan, dtype=float)
        lh["q"], lh["dq_cycle"] = qlin, dqlin
        lp = None
    else:
        print("Re-integrating retained linear UU motion on the SAME hardware...")
        total_mass_kg = float(np.sum(fstate[:8:2]))
        linear_design, linear_initial = _fixed_inventory_uniform_initial(
            linear_design,
            total_mass_kg,
            wall_state_source=fstate,
        )
        print(
            "  fixed gas inventory = "
            f"{1000.0 * total_mass_kg:.6f} g; "
            "gas state rebuilt uniformly for the linear theta=0 volumes"
        )
        lw, la, lt, lp, _ = _solve(
            linear_design,
            campaign_definition.wall_numerical_settings,
            linear_initial,
            "Linear same hardware and fixed inventory",
        )
        ldeg, lh = _forward_cycle(la, _extract_raw(lw, la, lt))
        deg, fh, lh = _common_grid(fdeg, fh, ldeg, lh)

    extrema = _candidate_extrema(deg, fh["q"])
    markers = _phase_markers(reference)

    outdir = report_path.parent
    stem = f"exchange_diagnostic_{h}h"
    full_png = outdir / f"{stem}_full.png"
    bp_png = outdir / f"{stem}_BP.png"
    hp_png = outdir / f"{stem}_HP.png"
    csv_path = outdir / f"{stem}_history.csv"
    summary_path = outdir / f"{stem}_summary.json"

    title_base = (
        f"{h}H Fourier champion index {best.get('index')} — "
        "same-hardware linear comparison"
    )

    _plot_comparison(
        full_png,
        title_base + " — full cycle",
        deg, fh, lh, reference, extrema, (0.0, 360.0),
    )
    _plot_comparison(
        bp_png,
        title_base + " — LP/BP exchange and transitions",
        deg, fh, lh, reference, extrema, (-30.0, markers["LP_end"] + 25.0),
    )
    _plot_comparison(
        hp_png,
        title_base + " — HP exchange and transitions",
        deg, fh, lh, reference, extrema,
        (markers["compression_end"] - 25.0, markers["HP_end"] + 25.0),
    )

    _write_csv(csv_path, deg, fh, lh)

    omega = float(fw.model.angular_speed)
    summary = {
        "report_snapshot": str(report_path),
        "candidate": {
            "index": best.get("index"),
            "proposal_index": best.get("proposal_index"),
            "harmonics": h,
            "recorded_efficiency": best["result"]["indicated_thermal_efficiency"],
            "recorded_power_W": best["result"]["indicated_power_w"],
            "reevaluated_efficiency": fp.thermal_efficiency,
            "reevaluated_power_W": fp.gas_power,
            "thermo": best["thermo"],
        },
        "linear_reference": {
            "comparison_mode": (
                "retained linear UU motion on Fourier champion hardware with identical total gas inventory"
                if not args.skip_linear
                else "kinematics only; thermodynamic baseline not evaluated"
            ),
            "original_linear_UU_recorded_efficiency": reference["efficiency"],
            "original_linear_UU_recorded_power_W": reference["power_W"],
            "same_hardware_efficiency": None if lp is None else lp.thermal_efficiency,
            "same_hardware_power_W": None if lp is None else lp.gas_power,
        },
        "linear_phase_boundaries_deg": markers,
        "fourier_extrema_deg": extrema,
        "LP_exchange": {
            "fourier": _phase_stats(
                deg, fh, 0.0, markers["LP_end"], omega
            ),
            "linear_same_hardware": (
                None if lp is None else
                _phase_stats(deg, lh, 0.0, markers["LP_end"], omega)
            ),
        },
        "HP_exchange": {
            "fourier": _phase_stats(
                deg, fh, markers["compression_end"], markers["HP_end"], omega
            ),
            "linear_same_hardware": (
                None if lp is None else
                _phase_stats(
                    deg, lh, markers["compression_end"], markers["HP_end"], omega
                )
            ),
        },
        "notes": [
            "Forward motor angle matches the Fourier comparison plots.",
            "Solid Fourier vs dashed retained-linear curves use identical current hardware and total gas inventory.",
            "Actual gas flow is simulated port flow; piston speed is only a kinematic driver.",
            "Correlation-regime fractions in phase summaries are trapezoidal physical-time fractions; the report-level gas diagnostics additionally provides heat-weighted domain fractions.",
        ],
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    print("\nRe-evaluated performance")
    print(f"  Fourier {h}H: eta={100*fp.thermal_efficiency:.6f}%  P={fp.gas_power:.4f} W")
    if lp is not None:
        print(f"  Linear same hardware: eta={100*lp.thermal_efficiency:.6f}%  P={lp.gas_power:.4f} W")
        print(
            "  Pure same-hardware kinematic gain: "
            f"{100*(fp.thermal_efficiency-lp.thermal_efficiency):+.6f} percentage points"
        )

    print("\nFourier extrema (forward motor angle)")
    for name, value in extrema.items():
        print(f"  {name:5s}: {value:8.3f} deg")

    print("\nWrote:")
    for p in (full_png, bp_png, hp_png, csv_path, summary_path):
        print(f"  {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
