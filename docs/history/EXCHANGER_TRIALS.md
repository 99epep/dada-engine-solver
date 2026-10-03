# Early air-wall and microtube trials

Condensed from the September 2026 geometry/air-wall trial notes. These results
used historical constant-property/Nusselt closures, rectangular collectors and
finite valve CdA. They are not predictions for the current circular-collector,
variable-transport model and must not become generic Research limits.

Current equations: [microtubes](../MICROTUBE_GAS_MODEL.md) and
[external streams](../EXTERNAL_STREAM_THERMAL_MODEL.md). Numerical acceleration
and periodic checks: [execution reference](../SOLVER_PERFORMANCE.md).

## What was learned

- The gas/wall/external-stream energy balance conserves storage and retains heat
  exchange during pauses. A uniform wall cannot generally reproduce axial
  counterflow profiles; Doty's steady UA does not identify resistance shares
  and wall capacity independently.
- More tube surface at fixed external capacity rate has diminishing benefit.
  Large external flow suppresses outlet-temperature change, not film resistance.
- Increasing count and shortening tubes at fixed surface reduces tube friction
  but can increase collector hold-up. A UA multiplier at unchanged volume,
  capacity and hydraulics is a sensitivity probe, not a realizable redesign.
- H_i and H_o need not have equal geometry. Actual event chronology and reflux
  must be measured; intended timing does not prescribe passive-valve transitions.
- Internal hydraulic losses already affect cycle work. Do not subtract them again
  as fictitious pump power. Excluded external fan consumption is not zero loss.

## Controlled sequence and artifacts

All powers below are indicated gas power. Negative values mean mechanical input,
not a negative-efficiency motor. The early comparison scale was 2 Hz and a 1 L
large-cylinder maximum volume, not the historical 66 L ceiling.

| Trial | Change / observation | Stored evidence |
|---|---|---|
| Abstract air/wall screening | Assumed 100 W/K films, 10 J/K walls and 0.05 kg/s streams; corrected tolerances gave ~23.637 W and 16.189% indicated efficiency | `outputs/motor_air_wall_screening*`, `motor_air_wall_refinement.json`; `examples/motor_air_wall_screening.py` |
| Geometry-connected four-bar | 1236 tubes, 31.75 mm; ~13.23 mL hold-up/HX and 38.25 J/K wall capacity; −39.031 W, not a functioning motor | `outputs/motor_hardware_screening*`; `examples/motor_hardware.toml` |
| Doubled count and length | 2472 tubes, 63.5 mm; surface/tube volume ×4, finite-air conductance not ×4; reference −23.70 W, lower-Lambda motion −26.94 W | `outputs/motor_doubled*`, `motor_lower_lambda*`; `examples/motor_hardware_doubled.toml` |
| Excess air, 25/175 °C | Peak-flow-based capacity sizing; opposed-Lambda result −19.211 W; doubling air changed power by ~0.078 W | `outputs/motor_excess_air_{reference,opposed,check}*`; corresponding `examples/motor_hardware_excess_air*.toml` |
| Higher source, 25/325 °C | Same inventory/motion/hardware comparison: +15.498 W, 4.440% indicated efficiency; not useful shaft power | `outputs/motor_325c*`; `examples/motor_demonstrator_325c_trial.toml`, `motor_hardware_325c.toml` |
| Constant-area parallel tubes | 3708 tubes, 42.3333 mm; tube resistance at fixed state/flow ×(1/2.25), header hold-up increased; original four-bar motion restored | `outputs/motor_parallel_325c*`; `examples/motor_hardware_parallel_325c.toml` |

`examples/motor_hardware_screening.py --help` describes the historical runner.
Configurations, reports, trajectories and refinement artifacts retain the exact
inputs. They are not overwritten or rerun by this documentation consolidation.
Frozen result values belong to their original property/source version.

## Calibration and hydraulic qualifications

The Doty nitrogen rows gave hot/cold sensible-heat imbalance around 5.1–5.3%
under equal constant Cp. This is not attribution to a particular heat leak.
Apparent-UA reduction reproduced the tabulated values within their uncertainties;
that is data-reduction consistency, not independent model prediction. Keep the
[Doty source and absolute-loss audit](../DOTY_SCREENING.md) and
[experimental comparison](../EXCHANGER_VALIDATION.md).

The legacy tube model allocates half friction and half header K to each port.
At equal port flow these sum to one exchanger; storage can make port flows
unequal. K uses total tube-passage velocity. Constant-Nusselt external longitudinal
flow, equivalent laminar passage pressure drop and fan estimates were screening
assumptions. Excess-air trials exceeded that exterior laminar regime, so their
fan figures are not validated ratings. For the new neutral external-stream
family, conductance is declared and liquid/fan/pump hydraulics are not invented.

Wall capacity follows material volume and declared additional capacity; it is
not reduced to accelerate convergence. Checkpoints/warm states remain guesses.
Refinement and small conservation residuals establish numerical consistency,
not physical calibration. Later [four-stage studies](../FOUR_STAGE_OPTIMIZATION.md)
retain the artificial-UA versus actual-hardware comparison.

## Historical reproduction entry points

Use the original files and matching implementation when reproducing exact numbers:

- `examples/motor_air_wall_screening.py`;
- `examples/motor_hardware_screening.py`, `motor_trials_refinement.py`;
- `examples/motor_parallel_comparison.py`;
- `outputs/doty_energy_balance_audit.csv`.

The [consolidation map](DOCUMENTATION_MAP.md) identifies the retired prose files.
Current model documentation, rather than old trial instructions, governs new studies.
