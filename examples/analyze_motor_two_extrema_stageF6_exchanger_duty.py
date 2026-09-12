"""Quantify H_i / H_o duty asymmetry for the Stage F6 champion.

Requires examples/plot_motor_two_extrema_stageF6_thermo.py.

Outputs:
  outputs/motor_two_extrema_stageF6_exchanger_duty.json
  outputs/motor_two_extrema_stageF6_exchanger_duty.csv
"""

from __future__ import annotations
import csv, json, math
import numpy as np

from compare_motor_motion_laws_stage7A5 import ROOT
from plot_motor_two_extrema_stageF6_thermo import _histories, _integrate_best

OUT_JSON = ROOT / "outputs" / "motor_two_extrema_stageF6_exchanger_duty.json"
OUT_CSV = ROOT / "outputs" / "motor_two_extrema_stageF6_exchanger_duty.csv"
THRESHOLDS = (0.01, 0.05, 0.10, 0.25, 0.50)


def integ(t, y):
    return float(np.trapz(np.asarray(y, float), t))


def time_fraction(t, condition):
    return float(np.trapz(np.asarray(condition, float), t) / (t[-1] - t[0]))


def metrics(t, useful_power):
    """Power is positive in the intended heat-transfer direction."""
    q = np.asarray(useful_power, float)
    qp = np.maximum(q, 0.0)
    qr = np.maximum(-q, 0.0)
    period = float(t[-1] - t[0])

    e_use = integ(t, qp)
    peak = float(np.max(qp))

    out = {
        "net_energy_j_per_cycle": integ(t, q),
        "useful_energy_j_per_cycle": e_use,
        "reverse_energy_j_per_cycle": integ(t, qr),
        "peak_useful_power_w": peak,
        "peak_reverse_power_w": float(np.max(qr)),
        "rms_power_w": math.sqrt(integ(t, q*q) / period),
        "mean_net_power_w": integ(t, q) / period,
        "fraction_cycle_useful_direction": time_fraction(t, q > 0.0),
        "fraction_cycle_reverse_direction": time_fraction(t, q < 0.0),
        "equivalent_full_peak_duty_factor": (
            e_use / (peak * period) if peak > 0 else None
        ),
    }
    for x in THRESHOLDS:
        out[f"fraction_cycle_above_{round(100*x)}pct_of_peak"] = (
            time_fraction(t, qp > x*peak) if peak > 0 else 0.0
        )
    return out


def exchanger_record(name, role, gas_power, air_power, wall_t, air_out_t, t):
    return {
        "role": role,
        "gas_side": metrics(t, gas_power),
        "external_air_side": metrics(t, air_power),
        "wall_temperature_min_k": float(np.min(wall_t)),
        "wall_temperature_max_k": float(np.max(wall_t)),
        "wall_temperature_swing_k": float(np.ptp(wall_t)),
        "external_air_outlet_temperature_min_k": float(np.min(air_out_t)),
        "external_air_outlet_temperature_max_k": float(np.max(air_out_t)),
        "external_air_outlet_temperature_swing_k": float(np.ptp(air_out_t)),
    }


def main():
    best, values, wrapper, angles, trajectory, performance = _integrate_best()
    h = _histories(wrapper, angles, trajectory)
    t = np.asarray(angles, float) / wrapper.model.angular_speed
    period = float(t[-1] - t[0])

    # Convert to positive = intended useful heat-transfer direction.
    # AirWallExchanger.rates(): gas_heat > 0 means wall -> gas;
    # air_heat > 0 means external air -> wall.
    hi = exchanger_record(
        "H_i", "heat input to working gas",
        h["gas_heat"][0],
        h["air_heat"][0],
        h["wall_temperature"][0],
        h["air_outlet_temperature"][0],
        t,
    )
    ho = exchanger_record(
        "H_o", "heat rejection from working gas",
        -h["gas_heat"][1],
        -h["air_heat"][1],
        h["wall_temperature"][1],
        h["air_outlet_temperature"][1],
        t,
    )

    hi_g, ho_g = hi["gas_side"], ho["gas_side"]
    comparison = {
        "H_o_to_H_i_peak_useful_gas_power_ratio":
            ho_g["peak_useful_power_w"] / hi_g["peak_useful_power_w"],
        "H_o_to_H_i_rms_gas_power_ratio":
            ho_g["rms_power_w"] / hi_g["rms_power_w"],
        "H_o_to_H_i_useful_gas_energy_ratio":
            ho_g["useful_energy_j_per_cycle"] / hi_g["useful_energy_j_per_cycle"],
        "H_o_to_H_i_useful_direction_time_fraction_ratio":
            ho_g["fraction_cycle_useful_direction"] /
            hi_g["fraction_cycle_useful_direction"],
        "H_o_minus_H_i_useful_direction_cycle_fraction":
            ho_g["fraction_cycle_useful_direction"] -
            hi_g["fraction_cycle_useful_direction"],
    }

    report = {
        "candidate": {
            "index": int(best["index"]),
            "parameters": values,
            "recorded_efficiency": float(best["indicated_thermal_efficiency"]),
            "recorded_indicated_power_w": float(best["indicated_power_w"]),
        },
        "reevaluated": {
            "frequency_hz": 1.0 / period,
            "cycle_period_s": period,
            "thermal_efficiency": performance.thermal_efficiency,
            "indicated_power_w": performance.gas_power,
            "heat_input_w": performance.heat_in_power,
            "heat_out_w": performance.heat_out_power,
        },
        "sign_convention": {
            "H_i_gas": "positive = wall -> working gas",
            "H_i_external_air": "positive = hot external air -> wall",
            "H_o_gas": "positive = working gas -> wall",
            "H_o_external_air": "positive = wall -> cold external air",
        },
        "relative_peak_thresholds": THRESHOLDS,
        "exchangers": {"H_i": hi, "H_o": ho},
        "comparison": comparison,
        "note": (
            "Wall storage decouples gas-side and external-air-side peak duties; "
            "peak gas heat rate alone is not an exchanger sizing requirement."
        ),
    }

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(report, indent=2) + "\n")

    rows = []
    for ex_name, ex in (("H_i", hi), ("H_o", ho)):
        for side in ("gas_side", "external_air_side"):
            row = {"exchanger": ex_name, "side": side, **ex[side]}
            row.update({
                "wall_temperature_swing_k": ex["wall_temperature_swing_k"],
                "external_air_outlet_temperature_swing_k":
                    ex["external_air_outlet_temperature_swing_k"],
            })
            rows.append(row)

    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"Stage F6 index {int(best['index'])}: "
        f"eta={100*performance.thermal_efficiency:.6f}% "
        f"P={performance.gas_power:.3f} W"
    )
    for name, ex in (("H_i", hi), ("H_o", ho)):
        print(f"\n{name} — gas side")
        m = ex["gas_side"]
        print(f"  useful duration : {100*m['fraction_cycle_useful_direction']:.2f}% cycle")
        print(f"  useful Q        : {m['useful_energy_j_per_cycle']:.4f} J/cycle")
        print(f"  reverse Q       : {m['reverse_energy_j_per_cycle']:.4f} J/cycle")
        print(f"  peak useful     : {m['peak_useful_power_w']:.2f} W")
        print(f"  RMS             : {m['rms_power_w']:.2f} W")
        print(f"  equiv peak duty : {100*m['equivalent_full_peak_duty_factor']:.2f}%")
        print(f"{name} — external-air side")
        m = ex["external_air_side"]
        print(f"  useful duration : {100*m['fraction_cycle_useful_direction']:.2f}% cycle")
        print(f"  peak useful     : {m['peak_useful_power_w']:.2f} W")
        print(f"  RMS             : {m['rms_power_w']:.2f} W")

    print("\nH_o / H_i gas-side comparison")
    for k, v in comparison.items():
        print(f"  {k}: {v:.6f}")
    print(f"\nWrote {OUT_JSON.relative_to(ROOT)}")
    print(f"Wrote {OUT_CSV.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
