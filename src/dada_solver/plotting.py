"""Diagnostic plots over one 0-360 degree cycle."""

from __future__ import annotations

from pathlib import Path
import os
import tempfile

import numpy as np

from dada_solver.dynamics import ThermodynamicModel
from dada_solver.integration import CycleIntegrationResult
from dada_solver.state import ThermodynamicState
from dada_solver.valves import ValveState
from dada_solver.humidity import (
    HumidityScreeningConfiguration,
    moisture_saturation_ratio_history,
)


def plot_cycle_diagnostics(
    cycle: CycleIntegrationResult,
    model: ThermodynamicModel,
    output_path: str | Path,
    humidity: HumidityScreeningConfiguration | None = None,
    charge_temperature: float | None = None,
    charge_pressure: float | None = None,
) -> None:
    """Write a multi-panel scientific diagnostic plot."""

    if not cycle.completed:
        raise ValueError("Plotting requires a completed cycle.")
    matplotlib_cache = Path(tempfile.gettempdir()) / "dada_solver_matplotlib"
    os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache))
    import matplotlib.pyplot as plt

    angles = np.degrees(cycle.angles)
    count = angles.size
    volumes = np.empty((2, count))
    pressures = np.empty((4, count))
    temperatures = np.empty((4, count))
    flows = np.empty((4, count))
    heat_rates = np.empty((2, count))
    valve_states = np.empty((2, count))

    for index, (angle, topology) in enumerate(
        zip(cycle.angles, cycle.topologies, strict=True)
    ):
        state = ThermodynamicState.from_array(cycle.states[:, index])
        instantaneous_volumes = model.volumes(float(angle))
        volumes[:, index] = (
            instantaneous_volumes.small_cylinder,
            instantaneous_volumes.large_cylinder,
        )
        pressures[:, index] = state.pressures(model.gas, instantaneous_volumes)
        temperatures[:, index] = state.temperatures(model.gas)
        rates = model.evaluate(float(angle), state, topology)
        flows[:, index] = (
            rates.flows.large_to_hot,
            rates.flows.small_to_cold,
            rates.flows.hot_to_small,
            rates.flows.cold_to_large,
        )
        heat_rates[:, index] = rates.cold_heat_rate, rates.hot_heat_rate
        valve_states[:, index] = (
            topology.hot_to_small is ValveState.OPEN,
            topology.cold_to_large is ValveState.OPEN,
        )

    panel_count = 8 if humidity is not None else 7
    figure, axes = plt.subplots(panel_count, 1, figsize=(11, 22), sharex=True)
    axes[0].plot(angles, volumes[0], label="Small cylinder")
    axes[0].plot(angles, volumes[1], label="Large cylinder")
    axes[0].set_ylabel("Volume (m³)")

    for values, label in zip(pressures, ("S", "L", "H_i", "H_o"), strict=True):
        axes[1].plot(angles, values, label=label)
    axes[1].set_ylabel("Pressure (Pa)")

    for values, label in zip(temperatures, ("S", "L", "H_i", "H_o"), strict=True):
        axes[2].plot(angles, values, label=label)
    axes[2].set_ylabel("Temperature (K)")

    for values, label in zip(
        flows,
        ("L to H_o", "S to H_i", "H_o to S", "H_i to L"),
        strict=True,
    ):
        axes[3].plot(angles, values, label=label)
    axes[3].set_ylabel("Mass flow (kg/s)")

    axes[4].plot(angles, heat_rates[0], label="Heat-in exchanger H_i")
    axes[4].plot(angles, heat_rates[1], label="Heat-out exchanger H_o")
    axes[4].set_ylabel("Heat rate (W)")

    axes[5].step(angles, valve_states[0], where="post", label="H_o to S")
    axes[5].step(angles, valve_states[1], where="post", label="H_i to L")
    axes[5].set_ylabel("Valve state")
    axes[5].set_yticks((0, 1), ("Closed", "Open"))

    axes[6].plot(angles, cycle.cold_heat - cycle.cold_heat[0], label="Q_i")
    axes[6].plot(angles, cycle.hot_heat - cycle.hot_heat[0], label="Q_o")
    axes[6].plot(angles, cycle.gas_work - cycle.gas_work[0], label="Gas work")
    axes[6].set_ylabel("Cumulative energy (J)")
    axes[6].set_xlabel("Cycle progress phi = |omega| t (deg)")

    if humidity is not None:
        if charge_temperature is None or charge_pressure is None:
            raise ValueError("Humidity plotting requires charge temperature and pressure.")
        saturation_ratios = moisture_saturation_ratio_history(
            cycle, model, humidity, charge_temperature, charge_pressure
        )
        for values, label in zip(
            saturation_ratios, ("S", "L", "H_i", "H_o"), strict=True
        ):
            axes[7].plot(angles, values, label=label)
        axes[7].axhline(1.0, color="black", linestyle="--", label="Saturation onset")
        axes[7].set_ylabel("Saturation ratio")
        axes[7].set_xlabel("Cycle progress phi = |omega| t (deg)")

    for axis in axes:
        axis.grid(True, alpha=0.3)
        axis.legend(loc="best", ncols=2)
        axis.set_xlim(0.0, 360.0)
    figure.tight_layout()
    figure.savefig(Path(output_path), dpi=150)
    plt.close(figure)
