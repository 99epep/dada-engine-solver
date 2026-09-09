# Four-bar motor demonstrator comparison

The user requested four-bar kinematics for the motor demonstrator. The existing
`shared_crank_rocker` implementation is reused, including finite piston rods,
a common crank pin and motor-direction reversal. No valve timing is prescribed:
check valves remain pressure driven. This is imposed constant-speed kinematics,
not a flywheel, startup, friction or structural dynamics simulation.

## Reproducible seed

`examples/motor_demonstrator_four_bar.toml` uses the earlier cooling-cell
candidate 001 proportions as a starting point, not a motor optimum. Crank,
coupler and rocker lengths relative to the ground distance are 0.40, 1.00 and
1.17. The ground distance is scaled to approximately 119 mm, retaining the
previous large-piston bore/stroke ratio near one. The complete mechanism can
be substantially larger than its ground distance because of the piston rods.

The large cylinder maximum volume is 1 litre, including clearance; small/large
swept-volume ratio is 0.84. Reservoirs are 298.15 and 448.15 K. The standalone
TOML starts at 2 Hz and 1 bar uniform filling pressure. The matched comparison
instead fixes total inventory from the piecewise-linear reference so a change of waveform
or crank origin cannot silently recharge the machine.

The piecewise-linear reference is `examples/motor_demonstrator_piecewise.toml`. Cylinder and
exchanger volumes, UA, flow coefficients and working-gas properties are shared.
The comparison keeps that hardware fixed when increasing frequency; it does
not increase UA or flow area with speed to preserve an artificial similarity.

## Exchanger boundary

The approximately 49.51 W/K conductances and orifice CdA values are numerical
seeds obtained by capacity/speed scaling the historical hardware to 1 litre and
2 Hz before changing reservoir temperatures. They are
not Doty-calibrated components. Scaled exchanger hold-up is also provisional.
Doty's overall gas/gas UA and tube-side pressure loss are not sufficient to
identify the required gas/reservoir closure, shell losses and gas hold-up.
Therefore this experiment isolates waveform and speed effects; it must not be
reported as a Doty-equipped, validated 100 W demonstrator.

## Running and interpreting the comparison

Run `PYTHONPATH=src python3 examples/motor_kinematics_comparison.py`.
Reports and simultaneous-state duty summaries are saved under
`outputs/motor_kinematics_comparison/` for piecewise-linear and four-bar cases at 2, 5 and
10 Hz. Nonconverged cases are recorded without inventing performance values.
Negative gas power is retained. A converged solution is not physical validation.

Use indicated gas power and thermal efficiency to compare the two waveforms;
use the duty files for pulse duration, peak flow and local reflux. Useful shaft
power still requires losses and auxiliary consumption. The 100 W useful target
and 66 litre cylinder ceiling remain design requirements, not achieved results.

## First matched comparison results

All six cases reached a periodic state. Integrated cycle energy residuals are
below 3e-12 J. These conservation checks do not establish exchanger applicability
or numerical refinement convergence. No 100 W useful-output result is claimed.

| Frequency (Hz) | Piecewise-linear indicated power (W) | Piecewise-linear efficiency | Four-bar indicated power (W) | Four-bar efficiency |
| --- | --- | --- | --- | --- |
| 2 | 23.49 | 16.07% | 21.33 | 13.62% |
| 5 | 6.72 | 2.51% | -18.98 | unavailable |
| 10 | -288.93 | unavailable | -456.32 | unavailable |

Negative power means mechanical work is consumed. All valve chronologies are
non-nominal. At 2 Hz both applicability verdicts are indeterminate (Mach is
unavailable); the four-bar 5 Hz and both 10 Hz cases also exceed the configured
pressure-equalization limit. Keep these diagnoses with the performance values.

At 2 Hz, peak inlet flow on H_o rises from approximately 6.85 g/s for the piecewise-linear
law to 9.33 g/s for the four-bar law. Thus similar power does not imply similar
exchanger pulse requirements. The next hardware screening must use the actual
four-bar flow histories and improve the coupled thermal/hydraulic design; simply
raising speed with fixed hardware does not reach the power target here.

`outputs/motor_kinematics_comparison/volume_laws.png` compares the imposed
volume laws. Full per-case reports and duty summaries accompany `summary.json`.

## Interpretation and future search

The piecewise-linear reference is not a proven optimum or an upper performance
bound. The legacy configuration type `ideal_piecewise_linear` is retained for
compatibility only. The mirrored four-bar seed is likewise one candidate, not
a motor constraint. Future searches may vary the two loops independently while
preserving actual shared-crank compatibility. See
[MOTOR_RESEARCH_OBJECTIVES.md](MOTOR_RESEARCH_OBJECTIVES.md).
