# Doubled exchanger and lower-Lambda trials

## Current user decisions

For this stage, exclude external-air aerodynamic losses and fan consumption
from the performance balance. Keep finite external-air heat-capacity rate and
heat exchange during working-gas flow pauses. Historical fan estimates remain
available in hardware diagnostics, but are not deducted or presented as zero
physical losses. Mechanical losses remain unknown and useful shaft power is
not claimed. Reservoir inlet temperatures remain 298.15 and 448.15 K, frequency
2–10 Hz, useful-output target approximately 100 W, large-cylinder ceiling 66 L.

The requested first trial doubles tube count from 1236 to 2472 and length from
31.75 to 63.5 mm. Diameter, wall thickness, pitch, header depth, air mass flow
(0.005 kg/s per side), properties and Nusselt scenarios remain unchanged.
Tube surface and tube volume quadruple. Total hold-up grows from 13.23 to
31.05 ml per exchanger because collectors and connections scale differently.
Static conductance increases from 3.046 to 12.088 W/K and wall capacity from
38.25 to 147.02 J/K. This does not quadruple the finite-air effective conductance:
the external capacity rate remains only 5.025 W/K per side.

## Independent mechanical trial

The reference four-bar output-point angles are -56 degrees on the small side
and +56 degrees on the large side. The second trial uses -80 and +75 degrees,
respectively. Other linkage ratios and the 119.38 mm ground scale are retained.
The independent angles were chosen through a bounded geometric scan to reduce
closure fractions at the reference valve-opening phases; this is a candidate
selection heuristic, not a motion-law or efficiency optimization.

Actual passive-valve pressure crossings are then measured in the integrated
cycle. No valve timing is imposed. The reported Lambdas below are
`(V_max-V_at_event)/V_swept` at explicitly named openings. They must not be
confused with arbitrary phase samples, prescribed piecewise-law targets, or a
universal definition of a single Lambda for non-nominal cycles. All four sampled
valve transitions, both cylinder Lambdas and event angles are retained in JSON.

Cylinder maximum volumes stay at 0.84 and 1 litre. Because linkage stroke
changes at fixed gas volume, piston bore also changes. Small/large strokes are
about 108.0/108.0 mm for the reference and 81.5/88.4 mm for the modified trial;
implied bores change from 99.0/108.0 to 114.0/119.4 mm. This effect must not be
hidden when pursuing mechanical efficiency or compactness.

## Results at 2 Hz

| Case | Indicated power (W) | Heat input (W) | Lambda_S at H_i-to-L opening | Lambda_L at H_o-to-S opening |
| --- | --- | --- | --- | --- |
| Original exchanger, reference mechanism | -39.03 | 30.78 | not extracted here | not extracted here |
| Doubled count and length, reference mechanism | -23.70 | 74.65 | 0.07661 | 0.08261 |
| Doubled count and length, independent output angles | -26.94 | 84.59 | 0.04638 | 0.05794 |

The two new cases converge, but consume mechanical work. No motor efficiency or
100 W useful-output result is claimed. The source-to-gas temperature differences
increase as suggested: hot-side time mean 16.94 to 19.19 K, cold-side time mean
22.32 to 25.31 K. That improves heat transfer but does not improve the complete
work balance for these candidates. It does not refute other lower-Lambda designs.
Tube Mach estimates fall below 0.12 in both trials; the prior 0.219 screening
failure is removed. Constant Nusselt, uniform-wall, laminar equivalent air
passage and uncalibrated Doty scaling limitations remain.

The heat-capacity rate of the imposed external air stream is now an explicit
bottleneck to investigate; increasing tube area alone has diminishing returns
at fixed air mass flow. Excluding fan penalties does not make external air a
constant-temperature reservoir throughout the exchanger.

## Runtime and periodic initial guesses

100 buffered checkpoint writes of about 1.2 kB took 0.043 s in the measured
workspace; median write time was 0.000129 s. This is a local buffered-I/O
measurement, not a forced-disk-sync benchmark. Checkpoints now default to every
10 cycles and at convergence. Measured cumulative checkpoint time was about
0.0016 s and 0.0232 s in the two runs, compared with roughly 44 s and 36 s of
integration/post-processing before plotting. I/O is not the dominant cost here.

Both trials reuse previous gas state and wall temperatures as initial guesses,
rescaling wall energies by the new heat capacity. Gas inventory is held fixed.
Optional bounded Aitken extrapolation updates wall initial guesses between
cycles every ten iterations; it never changes a state within an integrated
cycle. No extra heat or work is credited from an initial-guess jump. Acceptance
still requires the original periodic criterion on all ten states over an
unmodified cycle. Reported iteration counts are not physical warm-up times.
The trials required 32 and 23 iterations respectively.

One further unaccelerated cycle at tenfold tighter tolerances and half the
maximum angular step preserves the periodic threshold in both cases. Power
changes by about 0.00010 and 0.00011 W. These are numerical checks, not physical
validation. Unit tests cover contraction extrapolation, rejection of divergent
sequences, bounded corrections, and valid independent four-bar volume laws.

## Reproduction and plots

From the checkout:

```console
MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 examples/motor_hardware_screening.py --hardware examples/motor_hardware_doubled.toml --output-prefix outputs/motor_doubled --warm-start outputs/motor_hardware_screening.json --exclude-external-losses --accelerate-walls --plots
MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 examples/motor_hardware_screening.py --configuration examples/motor_demonstrator_lower_lambda_trial.toml --hardware examples/motor_hardware_doubled.toml --output-prefix outputs/motor_lower_lambda --warm-start outputs/motor_doubled_screening.json --exclude-external-losses --accelerate-walls --plots
PYTHONPATH=src python3 examples/motor_trials_refinement.py
```

Each output prefix has `_graphs.png`, `_screening.json`, `_ports.csv`,
`_trajectory.npz`, `_checkpoint.json` and `_refinement.json`. Graphs show all gas
and wall temperatures, inlet-air/gas temperature differences, cylinder p–V
loops, volume laws, signed port flows and separate external/gas heat rates.
Sampled valve-event angles use linear interpolation of pressure crossings;
they are diagnostics, not root-resolved lift histories.
