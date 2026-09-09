"""Compare fixed hardware at 2, 5 and 10 Hz; exchanger closure is uncalibrated.

Run: PYTHONPATH=src python3 examples/motor_kinematics_comparison.py
"""
from dataclasses import replace
import csv
import json
import math
from pathlib import Path

from dada_solver.configuration import load_simulation_configuration, ChargeConfiguration
from dada_solver.factory import build_model, build_initial_state, build_periodic_solver, initial_valve_topology
from dada_solver.performance import calculate_cycle_performance
from dada_solver.periodic import PeriodicStatus
from dada_solver.state import ThermodynamicState
from dada_solver.exchangers.duty import extract_exchanger_duty
from dada_solver.reporting import format_simulation_report
from dada_solver.results import extract_cycle_diagnostics
from dada_solver.validity import assess_cycle_validity

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / 'outputs/motor_kinematics_comparison'
OUTPUT.mkdir(parents=True, exist_ok=True)
piecewise = load_simulation_configuration(ROOT / 'examples/motor_demonstrator_piecewise.toml')
# Match inventory across waveforms; filling volume can depend on crank origin.
inventory = build_initial_state(piecewise, build_model(piecewise)).total_mass
rows = []
for kind in ('piecewise', 'four_bar'):
    seed = load_simulation_configuration(ROOT / f'examples/motor_demonstrator_{kind}.toml')
    for frequency in (2, 5, 10):
        config = replace(seed, angular_speed=-2*math.pi*frequency,
                         charge=ChargeConfiguration(temperature=298.15, total_mass=inventory))
        model = build_model(config)
        initial = build_initial_state(config, model)
        result = build_periodic_solver(config, model).solve(initial, initial_valve_topology())
        row = dict(kinematics=kind, frequency_hz=frequency, periodic_status=result.status.value,
                   exchanger_status='uncalibrated_constant_UA_and_orifice', total_mass_kg=inventory)
        if result.status is PeriodicStatus.CONVERGED:
            cycle = result.final_cycle
            performance = calculate_cycle_performance(cycle, ThermodynamicState.from_array(cycle.states[:,0]), config.angular_speed)
            validity = assess_cycle_validity(cycle, model, config.validity)
            report = format_simulation_report(result, performance, extract_cycle_diagnostics(cycle, model), validity, initial.total_mass)
            (OUTPUT / f'{kind}_{frequency}hz.txt').write_text(report)
            summary, columns = extract_exchanger_duty(cycle, model)
            with (OUTPUT / f'{kind}_{frequency}hz_duty.csv').open('w', newline='') as stream:
                writer = csv.writer(stream)
                writer.writerow(columns)
                writer.writerows(zip(*columns.values(), strict=True))
            (OUTPUT / f'{kind}_{frequency}hz_duty.json').write_text(json.dumps(summary, indent=2))
            row.update(indicated_gas_power_w=performance.gas_power,
                       thermal_efficiency=performance.thermal_efficiency,
                       energy_residual_j=performance.conservation.absolute_energy_residual)
        rows.append(row)
        (OUTPUT / 'summary.json').write_text(json.dumps(rows, indent=2))
        print(row, flush=True)
