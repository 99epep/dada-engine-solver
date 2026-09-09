# Original motion and constant-area parallel exchangers

The user requests returning to the original shared-crank four-bar motion and
trying more parallel tubes without adding contact area. Source inlets remain
25/325 deg C, frequency 2 Hz and large-cylinder maximum volume 1 L. The original
motion is the seed with small/large output angles approximately -56/+56 degrees,
not the piecewise-linear reference and not a future symmetry requirement.

## Geometry and comparison boundary

Per exchanger, tube count increases from 2472 to 3708 (+50%) while length falls
from 63.5 to 42.3333 mm. Inner diameter, wall thickness and pitch are unchanged.
Internal area remains 0.162736887 m^2; external tube area, tube gas volume
(13.4258 mL) and tube metal heat capacity also remain constant. The estimated
two-header volume grows from 15.625 to 23.2563 mL. Total internal hold-up,
including the unchanged 2 mL connections, grows from 31.0508 to 38.6820 mL per
exchanger. The gas-side envelope changes from 62.5 x 62.5 x 67.5 mm to
76.25 x 76.25 x 46.3333 mm, excluding casing and external ducting.

Under the implemented laminar law, at fixed density and total flow, tube
resistance scales with length/count, giving 1/2.25 of its former value. The
header loss coefficient is retained with tube-passage velocity as its reference;
real header distribution remains unvalidated. Valve CdA is unchanged. The
external hydraulic diameter changes slightly with packing, so static overall
conductance changes from 12.08769 to 12.09395 W/K despite identical contact area.
Constant-property/Nusselt closures and the lack of empirical Doty calibration
remain limitations. External aerodynamic losses and fan consumption are excluded;
wall storage, pause heat transfer and finite external capacity remain active.

## Results

Gas inventory is preserved by warming both new cases from the same previous
325 deg C result. The original-motion-only control isolates the geometric effect.

| Quantity | Previous opposed motion | Original motion | Original + parallel tubes |
| --- | ---: | ---: | ---: |
| Indicated power (W) | 15.498 | 19.978 | 37.626 |
| External heat input (W) | 349.050 | 344.365 | 353.593 |
| Indicated efficiency (%) | 4.440 | 5.801 | 10.641 |
| Maximum forward H_i path pressure drop (kPa) | 9.082 | 8.305 | 3.720 |
| Maximum forward H_o path pressure drop (kPa) | 19.105 | 19.199 | 9.270 |

Path pressure differences include tube halves, headers and the outlet valve,
measured only while both adjacent ports flow forward above 1e-8 kg/s. They are
not tube-only pressure drops or closed-valve pressure differences. Hi and Ho
maxima occur in different phases. Raw signed flow histories retain local reflux.

Parallelization reduces peak path pressure differences by 55.2% for H_i and
51.7% for H_o relative to the original-motion control. At almost unchanged
thermal conductance and similar heat input, the indicated power gain is 17.65 W.
This supports pursuing hydraulic improvements before adding area in this local
comparison; it does not certify exchanger capacity for 100 W useful output.
Mechanical losses remain unknown. The minimum external/internal peak capacity
ratio is 23.26, with maximum external temperature change 1.625 K.

Both new cases satisfy the existing periodic criterion. Energy conservation
including wall storage has residual below 3e-10 J and mass drift below 3e-18 kg.
The final parallel case takes 32 iterations, about 53 seconds, including
0.0023 seconds of checkpoint writes. No core solver equations were changed.

## Reproduction and artifacts

Use `MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src` for these commands:

```sh
python3 examples/motor_hardware_screening.py --configuration examples/motor_demonstrator_original_325c.toml --hardware examples/motor_hardware_325c.toml --output-prefix outputs/motor_original_325c --warm-start outputs/motor_325c_screening.json --exclude-external-losses --accelerate-walls --plots
python3 examples/motor_hardware_screening.py --configuration examples/motor_demonstrator_original_325c.toml --hardware examples/motor_hardware_parallel_325c.toml --output-prefix outputs/motor_parallel_325c --warm-start outputs/motor_325c_screening.json --exclude-external-losses --accelerate-walls --plots
python3 examples/motor_parallel_comparison.py
```

Both prefixes contain reports, trajectories, signed flows, checkpoints, angular
graphs with valve-event markers and supplementary P–V plots. The comparison
script reconstructs pressures from saved states and configurations, writing
`outputs/motor_parallel_comparison.json`.

An unaccelerated verification cycle with rtol 1e-9, tenfold tighter absolute
tolerances and maximum angular step pi/720 yields 37.626110 W: a change of
-0.000017 W. Its scaled periodic error is 0.120, below the threshold of 1.
See `outputs/motor_parallel_325c_refinement.json`.
