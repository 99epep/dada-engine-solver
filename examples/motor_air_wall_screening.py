"""Periodic air/wall coupling smoke study, with explicit uncalibrated inputs."""
from pathlib import Path
import json
import numpy as np

from dada_solver.configuration import load_simulation_configuration
from dada_solver.factory import build_model, build_initial_state
from dada_solver.exchangers.air_wall import AirWallExchanger, AirWallMotor

ROOT = Path(__file__).resolve().parents[1]
config = load_simulation_configuration(ROOT/'examples/motor_demonstrator_piecewise.toml')
model = build_model(config)
# Numerical scenarios, not a resistance split inferred from Doty.
hot = AirWallExchanger(100, 100, 10, .05, 1005, 448.15)
cold = AirWallExchanger(100, 100, 10, .05, 1005, 298.15)
wrapper = AirWallMotor(model, hot, cold)
state = np.r_[build_initial_state(config, model).as_array(), 10*448.15, 10*298.15]
converged = False
for cycle_number in range(1, 51):
    angle, history = wrapper.integrate_cycle(state, rtol=1e-8,
        atol=np.array([1e-13, 1e-8]*4 + [1e-7]*2 + [1e-8]*5))
    final = history[:10, -1]
    error = float(np.max(np.abs(final-state)/(1e-12 + 1e-6*np.maximum(np.abs(final), np.abs(state)))))
    initial_energy = state[1:8:2].sum()+state[8:10].sum()
    final_energy = final[1:8:2].sum()+final[8:10].sum()
    residual = final_energy-initial_energy-history[10:12, -1].sum()+history[14, -1]
    print(f'cycle={cycle_number}, scaled_periodic_error={error:.6g}', flush=True)
    state = final
    if error <= 1:
        converged = True
        break
period = 2*np.pi/model.angular_speed
qi, qo = history[10:12,-1]/period
power = history[14,-1]/period
report = dict(status='converged' if converged else 'maximum_cycles', cycles=cycle_number,
              calibration='uncalibrated_wall_and_air_scenario', frequency_hz=1/period,
              scaled_periodic_error=error, total_energy_residual_j=float(residual),
              last_cycle_indicated_power_w=float(power),
              last_cycle_external_heat_input_w=float(qi), last_cycle_external_heat_out_w=float(qo),
              thermal_efficiency=float(power/qi) if converged and power>0 and qi>0 and qo<0 else None,
              useful_power_w=None, fan_power_w=None,
              final_state=state.tolist(),
              integration_relative_tolerance=1e-8,
              integration_absolute_tolerances=[1e-13, 1e-8]*4 + [1e-7]*2 + [1e-8]*5,
              wall_temperature_ranges_k=[list(map(float,(history[i].min()/10,history[i].max()/10))) for i in (8,9)])
output = ROOT/'outputs/motor_air_wall_screening.json'
output.write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
