# Temporary higher-temperature motor trial

The current user-selected source inlets are 25/325 deg C (298.15/598.15 K).
Optimize dimensions at this wider temperature difference, then lower the hot
source progressively. Efficiency remains the objective, with approximately
100 W useful power, 2–10 Hz and at most 66 L maximum enclosed large-cylinder
volume. This first controlled comparison remains at 2 Hz and 1 L.

## Controlled comparison

`examples/motor_demonstrator_325c_trial.toml` preserves the previous independent
four-bar configuration (small output -85 degrees, large output +40 degrees),
volumes, charge and valves. `examples/motor_hardware_325c.toml` preserves the
2472 tubes per exchanger, 63.5 mm tube length and fixed excess-air sizing rule.
Both the configuration reservoir and the hardware hot-air inlet are updated.
The external hot-air density is scaled from 0.788 kg/m^3 by 448.15/598.15 at
unchanged assumed inlet pressure. Other constant transport, heat-capacity and
material properties and Nusselt numbers are deliberately retained for this
controlled model comparison; their temperature dependence is not validated.
The warm start preserves gas inventory from the previous 175 deg C trial.

External aerodynamic losses and fan consumption remain excluded by user
choice; thermal-film resistance, finite external capacity and wall storage
remain included. Large external flow is an imposed research boundary, not a
validated ventilation design. The exchanger is not empirically calibrated to
Doty. Positive indicated work, if obtained, is not a useful-shaft-power result.

## Reproduction

```sh
MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 examples/motor_hardware_screening.py --configuration examples/motor_demonstrator_325c_trial.toml --hardware examples/motor_hardware_325c.toml --output-prefix outputs/motor_325c --warm-start outputs/motor_excess_air_opposed_screening.json --exclude-external-losses --accelerate-walls --plots
```

The report, trajectory, signed port flows, checkpoint, six angle-based panels
and separate P–V plots share the `outputs/motor_325c` prefix. Retain P–V plots
as supplementary work diagnostics; the user currently finds the angle-based
panels with common valve-event lines more useful for interpreting phases.

## Results

| Quantity | Previous 25/175 deg C | New 25/325 deg C |
| --- | ---: | ---: |
| Indicated power (W) | -19.211 | +15.498 |
| External heat input (W) | 121.407 | 349.050 |
| External heat rejection magnitude (W) | 140.691 | 333.516 |
| Indicated thermal efficiency | No positive motor output | 4.440% |
| Mean hot-inlet minus H_i gas temperature (K) | 10.323 | 29.677 |
| Mean H_o gas minus cold-inlet temperature (K) | 11.962 | 28.357 |

The higher source temperature produces positive net gas work with unchanged
geometry. It does not establish useful output or an optimized design. The
mean H_i gas temperature is 295.323 deg C and mean H_o is 53.357 deg C.
The external/internal peak capacity-rate ratios are 25.45 and 24.45; maximum
external air temperature change is 1.605 K. Internal tube peak Re is 701 and
peak Mach is 0.099 under the constant-property screening assumptions.

The periodic solver converges in 42 iterations (about 62 seconds); checkpoint
writes take 0.0066 seconds. Mass drift is below 4e-18 kg and the energy residual
including wall storage is below 3e-10 J. The residual change in cycle storage
corresponds to 0.036 W; the displayed input/output difference therefore need
not exactly equal the displayed work at this stopping tolerance.

A subsequent unaccelerated cycle with rtol 1e-9, tenfold tighter component
absolute tolerances and maximum angular step pi/720 gives 15.497574 W,
a difference of -0.000076 W, and still meets the periodic criterion (scaled
error 0.504). See `outputs/motor_325c_refinement.json`.

Doubling the capacity-rate sizing factor to 40, using
`examples/motor_hardware_325c_air_check.toml` and warm-starting from
`outputs/motor_325c_screening.json`, converges in 13 iterations. The report
`outputs/motor_325c_air_check_screening.json` gives 15.696682 W and 4.480%
indicated efficiency. The increase is 0.199 W (1.28% of the base indicated
output), with mean approach changes below 0.30 K. External capacity has a
small but nonzero residual effect; film resistance is still included.
