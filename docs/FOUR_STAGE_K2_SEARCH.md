# Seven-variable four-stage motion search on K2

The new `FourStageVolumeKinematics` implements four shared linear stages with
seven independent parameters. Time is expressed as a fraction of the forward
motor cycle; T0=0 and T4=1 closes the cycle.

| Knot | L swept-volume fraction | S swept-volume fraction |
| --- | --- | --- |
| T0 | 1 | b_s |
| T1 | b_l | 1 |
| T2 | 0 | a_s |
| T3 | a_l | 0 |
| T4 | 1 | b_s |

Zero means the original K2 minimum enclosed volume, including clearance. One
means its maximum enclosed volume. The four intermediate levels are independent;
no ordering between a and b is imposed. Linear interpolation supplies analytic,
piecewise-constant derivatives; derivative jumps are passed to the integrator
as breakpoints. The module adapts forward motor time to study angle so that
the existing factory motor reversal produces the requested chronological order.
This is an explicitly motor-oriented input representation, not a second factory
reversal. Tests check the resulting motor trajectories at every specified knot.

The thermodynamic comparison fixes K2 gas inventory, cylinder volume limits,
2 Hz speed, 25/325 deg C sources, asymmetric microtube geometry, hydraulic losses
at their reference values, wall thermal storage and numerical tolerances.
Efficiency uses external-source heat; mechanical losses remain unknown.

## First search

```sh
PYTHONPATH=src python3 examples/optimize_motor_four_stage_k2.py --budget-seconds 480 --evaluations 24
PYTHONPATH=src python3 examples/plot_motor_four_stage_k2.py
```

The first candidate uses the free target's extrema in forward motor time, with
L maximum shifted to time zero. Subsequent candidates alternate broad Sobol
exploration (one in four) and Sobol perturbations of the best feasible candidate.
All seven coordinates can change simultaneously. Initial neighborhood half-width
is 0.12, falling to 0.06 after sequence index 16. Four stages must each occupy
at least 2% of the cycle; intermediate levels lie in [0,1]. These are explicit
search bounds, not claims of mechanical feasibility. The initial target-derived
seed is not an imposed waveform shape for later candidates.

The existing screening limits are retained: at least 40 W indicated power,
at most 12 bar absolute, 850 K, 0.08 kg/s, Mach 0.2 and tube Reynolds below 2300,
plus the existing thermodynamic validity verdict. Non-converged evaluations do
not generate invented physical constraint violations. The greatest efficiency
among feasible candidates is archived independently of the latest evaluation.

Every evaluated candidate is appended to
`outputs/motor_four_stage_k2/history.jsonl`, including its seven parameters,
convergence history, diagnostics, constraint failures, warm-start source and
last complete state. `definition.json` records the search basis;
`report.json` stores the current best and phase summary. Restart continues
from the recorded sequence index, with exact completed-candidate deduplication.
Warm starts use the nearest prior converged state in the seven dimensionless
coordinates; inventory and wall capacities are fixed in this experiment.
The overall deadline is checked cooperatively inside integration. Interrupted
work is recorded as incomplete, never as an efficiency result.

This deliberately small experimental driver reuses the K2 evaluator and does
not replace the persistent campaign framework or implement SLSQP. A finite
initial exploration establishes only the best observed candidate, not a global
optimum. Piecewise-linear velocity jumps have no mechanical acceleration or
inertia-loss model attached to them.

## Initial recorded results

The first phase requested 480 seconds and finished its in-flight candidate in
506.74 seconds. Twelve candidates converged; nine met the screening constraints.
The highest feasible efficiency was candidate 6: **18.829771%**, with
**53.858662 W** indicated gas power. The initial target-derived candidate gave
16.209185%. These are initial search results, not a certified optimum.
The free K2 reference remains 20.561604% and the six-bar reference 20.381817%.

| Parameter | Best observed value |
| --- | ---: |
| T1 / period | 0.402135611 |
| T2 / period | 0.626363759 |
| T3 / period | 0.971676933 |
| a_l | 0.392564913 |
| b_l | 0.059046611 |
| a_s | 0.474203422 |
| b_s | 0.126399369 |

The last segment occupies 2.8323% of the cycle (14.16 ms at 2 Hz), near the
2% duration lower bound. This proximity should remain visible when interpreting
future improvements. No acceleration penalty is included in this kinematic
experiment. Full suite: 326 tests passed, including seven new chronology,
breakpoint, derivative, periodicity and parameter-validity checks.

Artifacts in `outputs/motor_four_stage_k2/` include `best_motion.png`, its
vector `best_motion.svg`, and `best_motion.csv` in forward motor time.
Future runs also accept `--candidate-seconds 90` to bound difficult candidates
cooperatively; callback checks may allow a short in-flight overrun.
