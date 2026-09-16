"""Phase-by-phase thermal/work/hydraulic diagnostics for the retained 7D motor.

Reports, for each of the four phases:
- wall <-> working-gas heat exchanged by H_i and H_o;
- gas boundary work integral;
- transported mass and peak mass flow on each hydraulic link;
- screening hydraulic dissipation integral sum(|dp|*|m_dot|/rho_upstream dt).

The hydraulic term is a comparison proxy, not a full compressible exergy model.

Outputs:
  outputs/motor_four_stage_phase_diagnostics.json
  outputs/motor_four_stage_phase_diagnostics.csv
"""

from __future__ import annotations

import csv
import json
import math

import numpy as np

from dada_solver.factory import initial_valve_topology
from dada_solver.state import ThermodynamicState
from plot_motor_four_stage_7d_cycle_alt import ROOT, _build_retained_design


OUTPUT_JSON = ROOT / "outputs" / "motor_four_stage_phase_diagnostics.json"
OUTPUT_CSV = ROOT / "outputs" / "motor_four_stage_phase_diagnostics.csv"

PHASE_NAMES = (
    "Low-pressure exchange",
    "Compression",
    "High-pressure exchange",
    "Expansion",
)

# label, ModelRates.flows attribute, first control-volume index, second index
# State/volume order is S, L, H_i, H_o in the retained plotting convention.
LINKS = (
    ("L <-> H_o", "large_to_hot", 1, 3),
    ("S <-> H_i", "small_to_cold", 0, 2),
    ("H_o -> S", "hot_to_small", 3, 0),
    ("H_i -> L", "cold_to_large", 2, 1),
)

TRAPZ = getattr(np, "trapezoid", np.trapz)


def _integral(values, times, mask):
    x = np.asarray(times[mask], dtype=float)
    y = np.asarray(values[mask], dtype=float)
    if len(x) < 2:
        return 0.0
    return float(TRAPZ(y, x))


def _key(label):
    return (
        label.replace(" ", "_")
        .replace("<->", "to")
        .replace("->", "to")
    )


def main():
    control, wrapper, angles, trajectory, performance = _build_retained_design()

    omega = wrapper.model.angular_speed
    times = np.asarray(angles, dtype=float) / omega
    cycle_time = 2.0 * math.pi / omega

    motion = control["motion_parameters"]
    fractions = np.asarray(
        [0.0, motion["t1"], motion["t2"], motion["t3"], 1.0],
        dtype=float,
    )
    phase_times = fractions * cycle_time
    phase_degrees = fractions * 360.0

    n = len(angles)
    q_hi_to_gas = np.empty(n)
    q_ho_to_gas = np.empty(n)
    gas_work_rate = np.empty(n)

    mass_flow = {label: np.empty(n) for label, *_ in LINKS}
    hydraulic_power = {label: np.empty(n) for label, *_ in LINKS}

    for i, angle in enumerate(angles):
        values = np.asarray(trajectory[:, i], dtype=float)

        # AirWallMotor quadratures: 12/13 are wall -> working-gas heat rates;
        # 14 is the gas work term. Multiplication by omega recovers W.
        full_rates = wrapper.derivative(float(angle), values) * omega
        q_hi_to_gas[i] = full_rates[12]
        q_ho_to_gas[i] = full_rates[13]
        gas_work_rate[i] = full_rates[14]

        gas = ThermodynamicState.from_array(values[:8])
        volumes = wrapper.model.volumes(float(angle))
        pressures = gas.pressures(wrapper.model.gas, volumes)

        model_rates = wrapper.model.evaluate(
            float(angle), gas, initial_valve_topology()
        )

        masses = np.asarray(values[:8:2], dtype=float)
        volume_values = np.asarray(
            [
                volumes.small_cylinder,
                volumes.large_cylinder,
                volumes.cold_heat_exchanger,
                volumes.hot_heat_exchanger,
            ],
            dtype=float,
        )
        densities = masses / volume_values

        for label, attribute, first, second in LINKS:
            mdot = float(getattr(model_rates.flows, attribute))
            mass_flow[label][i] = mdot

            upstream = first if mdot >= 0.0 else second
            rho_up = float(densities[upstream])
            dp = abs(float(pressures[first] - pressures[second]))
            hydraulic_power[label][i] = dp * abs(mdot) / rho_up

    phases = []
    for j, name in enumerate(PHASE_NAMES):
        left = phase_times[j]
        right = phase_times[j + 1]
        if j == len(PHASE_NAMES) - 1:
            mask = (times >= left) & (times <= right)
        else:
            mask = (times >= left) & (times < right)

        row = {
            "phase": name,
            "start_deg": float(phase_degrees[j]),
            "end_deg": float(phase_degrees[j + 1]),
            "duration_deg": float(phase_degrees[j + 1] - phase_degrees[j]),
            "duration_ms": 1000.0 * (right - left),
            "Q_Hi_to_gas_J": _integral(q_hi_to_gas, times, mask),
            "Q_Ho_to_gas_J": _integral(q_ho_to_gas, times, mask),
            "Q_total_to_gas_J": _integral(q_hi_to_gas + q_ho_to_gas, times, mask),
            "gas_work_J": _integral(gas_work_rate, times, mask),
        }

        hyd_total = 0.0
        for label, *_ in LINKS:
            key = _key(label)
            transported = _integral(np.abs(mass_flow[label]), times, mask)
            dissipation = _integral(hydraulic_power[label], times, mask)
            peak = (
                float(np.max(np.abs(mass_flow[label][mask])))
                if np.any(mask)
                else 0.0
            )
            row[f"{key}_transported_mass_g"] = 1000.0 * transported
            row[f"{key}_peak_mass_flow_g_s"] = 1000.0 * peak
            row[f"{key}_hydraulic_dissipation_J"] = dissipation
            hyd_total += dissipation

        row["hydraulic_dissipation_total_J"] = hyd_total
        row["hydraulic_dissipation_average_W"] = hyd_total / (right - left)
        phases.append(row)

    totals = {
        "cycle_time_s": cycle_time,
        "frequency_hz": 1.0 / cycle_time,
        "thermal_efficiency": performance.thermal_efficiency,
        "gas_power_W": performance.gas_power,
        "heat_input_power_W": performance.heat_in_power,
        "Q_Hi_to_gas_J": float(TRAPZ(q_hi_to_gas, times)),
        "Q_Ho_to_gas_J": float(TRAPZ(q_ho_to_gas, times)),
        "gas_work_J": float(TRAPZ(gas_work_rate, times)),
        "hydraulic_dissipation_total_J": float(
            sum(TRAPZ(values, times) for values in hydraulic_power.values())
        ),
    }

    report = {
        "basis": {
            "motion_parameters": motion,
            "charge_pressure_pa": 100000.0,
            "hydraulic_note": (
                "Screening proxy: integral sum(|dp|*|m_dot|/rho_upstream dt); "
                "not a full compressible exergy-destruction calculation."
            ),
        },
        "totals": totals,
        "phases": phases,
    }

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(report, indent=2) + "\n")

    columns = list(phases[0].keys())
    with OUTPUT_CSV.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(phases)

    print()
    print(
        f"{'Phase':24s} {'deg':>8s} {'ms':>8s} "
        f"{'Q_Hi J':>10s} {'Q_Ho J':>10s} {'Wgas J':>10s} "
        f"{'Hyd J':>10s} {'Hyd Wavg':>10s}"
    )
    print("-" * 106)
    for row in phases:
        print(
            f"{row['phase']:24s} "
            f"{row['duration_deg']:8.2f} "
            f"{row['duration_ms']:8.2f} "
            f"{row['Q_Hi_to_gas_J']:10.4f} "
            f"{row['Q_Ho_to_gas_J']:10.4f} "
            f"{row['gas_work_J']:10.4f} "
            f"{row['hydraulic_dissipation_total_J']:10.5f} "
            f"{row['hydraulic_dissipation_average_W']:10.3f}"
        )

    print("\nPer-link details:")
    for row in phases:
        print(f"\n{row['phase']}")
        for label, *_ in LINKS:
            key = _key(label)
            print(
                f"  {label:10s}: "
                f"peak={row[f'{key}_peak_mass_flow_g_s']:8.3f} g/s, "
                f"transported={row[f'{key}_transported_mass_g']:8.4f} g, "
                f"hyd={row[f'{key}_hydraulic_dissipation_J']:10.6f} J"
            )

    print(f"\nWrote {OUTPUT_JSON.relative_to(ROOT)}")
    print(f"Wrote {OUTPUT_CSV.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
