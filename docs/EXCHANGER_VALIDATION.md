# Microtube gas-model verification and experimental comparison

## Reproduce the checks

```sh
PYTHONPATH=src python3 -m pytest -q
PYTHONPATH=src python3 examples/validate_microtube_gas_model.py
PYTHONPATH=src python3 examples/compare_microtube_gas_model.py
```

The first example writes `outputs/microtube_gas_validation.json`. The second
writes `outputs/microtube_gas_model_comparison.json`, retaining configuration,
hardware, transport settings, wall numerical settings, periodic history,
performance, domain fractions and final states. No CFD, DSMC or optimization of
new physical coefficients is performed.

## Verification is distinct from experimental validation

Analytic tests cover the incompressible limit, signed flow symmetry, equality
with the old arithmetic-mean-density Poiseuille closure at constant properties,
half-tube length, the laminar asymptote, thermal-entry behavior, turbulent Darcy
inversion and the reuse of Gnielinski. Domain tests cover Reynolds gaps,
compressibility, rarefaction, temperature bounds, missing accommodation and
incompatible mean-free-path conventions. Variable-film tests poison the legacy
UA to prove it is not being used, retain heat at rest, and verify wall energy
balance. Time/heat fractions are tested on an intentionally nonuniform grid.

The Graur/Ewart regression reproduces the published P=5 polynomial fits from
Tables 1–2. It is **not** a second independent validation against unpublished
raw measurements. Accordingly the report has `measured=null`,
`reference_kind=published_regression_not_raw_measurement`, the published fit,
computed S and coefficient uncertainty envelope. The VHS convention is retained.
Silica TMAC values and uncertainties are reference metadata only. An additional
explicit HS sensitivity test demonstrates that alternative declared A1/A2
change hydraulic flow without assigning a metal TMAC or a thermal-jump model.

## Doty absolute hydraulic comparison

Source: [Doty et al. 1991](https://dotynmr.com/download/pubs/1991_HTE_Doty_HeatExchanger.pdf),
Eq.7 and Tables 1–2; SI measurements and reported uncertainties remain in
`examples/data/doty_1991_nitrogen_reference.csv` and
`examples/data/doty_1991_helium_reference.csv`.

Assumptions are explicit: 309 equal-flow tubes, 0.33 mm internal diameter,
127 mm length; representative temperature `(T3+T4)/2`; reported pressure
interpreted as mean tube pressure; no header loss in the tube-only comparison.
Variable mu uses Sutherland N2 or the NIST helium table; ideal-gas density is
used. At fixed mean pressure, inverting the compressible formula gives exactly
the same pressure drop as Poiseuille with that mean density. This identity is
checked separately rather than advertised as agreement with experiment.

| Gas | Flow, g/s | Measured tube dP, Pa | Variable-transport prediction, Pa | Signed relative error |
| --- | ---: | ---: | ---: | ---: |
| N2 | 0.470 | 5100 | 3903 | -23.46% |
| N2 | 0.840 | 9500 | 6819 | -28.22% |
| N2 | 0.930 | 11100 | 7574 | -31.77% |
| He | 0.117 | 2300 | 3319 | +44.30% |
| He | 0.079 | 1500 | 2206 | +47.05% |
| He | 0.213 | 4400 | 5477 | +24.48% |

None of these six predictions falls inside the reported pressure-drop output
uncertainty. This is retained as a model discrepancy, not removed with a fitted
multiplier. The JSON includes `measured`, `predicted`, `absolute_error`,
`relative_error`, `inside_experimental_uncertainty`, `model_validity` and `source`
for both the fixed-300-K transport comparator and the variable-property model.
That fixed-300-K comparator is an explicitly defined screening comparison, not
an assertion about the exact properties used by Doty or the original motor.

The independent He reconstruction above is not identical to Doty's own Eq.7
numbers (3000, 2000, 5500 Pa). Doty's published predictions exceed measurements
by approximately 30%, 33% and 25%; the independent property convention produces
the larger first two discrepancies shown here. See [DOTY_SCREENING.md](DOTY_SCREENING.md).
Unreported diameter, distribution, axial temperature and property-reduction
uncertainties are not invented to make an uncertainty-overlap flag true.

Measured UA/effectiveness and four temperatures remain available. A complete
predicted gas/gas UA is deliberately **unavailable** in this tube-hydraulic
comparison: the new DADA external-air equivalent passage is not an independently
validated reconstruction of Doty's shell side. Claiming thermal validation by
feeding the measured UA back into the model would be circular. An independently
specified shell-side/axial thermal reconstruction is a genuine remaining step
for integral Doty thermal validation.

## Same champion, unchanged thermodynamics

The comparison uses the selected independent S/L six-bars, K2 volume limits,
fixed gas inventory, 2 Hz, 25/325 deg C reservoirs, asymmetric microtube geometry,
normal hydraulic loss coefficients and the existing dynamic wall states.
Only internal transport and thermal-entry treatment vary. Efficiency is
indicated gas work / **external-source** heat. Useful mechanical power remains
unavailable; mechanical losses are unknown; external-air aerodynamic losses
remain excluded from the trial balance.

| Internal model | Indicated efficiency | Indicated power | Periodic cycles |
| --- | ---: | ---: | ---: |
| Legacy replay | 20.381710% | 51.406689 W | 1 |
| Variable transport, Nu=3.66 comparison | 20.195227% | 52.032713 W | 13 |
| Variable transport + thermal entry | 20.246636% | 52.376421 W | 13 |

The production change is -0.135074 percentage point relative to the legacy
replay. More power does not imply more efficiency: external heat input changes.
The wall acceleration remains an initial-guess accelerator only; all reported
solutions satisfy an ordinary cycle's unchanged periodic convergence test.
The three evaluation times on this run were 3.49, 94.03 and 103.90 seconds;
they include different cycle counts and should not be compared as identical
workloads. Transport/correlation calculations remain scalar algebra; no field
solver is called inside the integration.

For the full-entry case, thermal bulk-state diagnostics give peak Re about
381 at Hi and 308 at Ho, peak Ma about 0.0242 and 0.0390, and maximum Kn about
3.81e-4 and 1.92e-4. Pr is about 0.690 at Hi and 0.703–0.708 at Ho. Peak Gz is
about 1.87 and 2.84. Fundamental Wo is below about 0.189 and 0.308 respectively.
Thus the selected tubes remain in the no-slip continuum domain, and the entry
correction is modest. The JSON additionally reports hydraulic extrema evaluated
at actual upstream port temperatures; these are not interchangeable with the
thermal bulk-state extrema. A valid correlation domain is not experimental
validation of the pulse response or the whole machine.

## Historical preservation

The original microtube hardware regression is reproduced with:

```sh
PYTHONPATH=src python3 examples/motor_hardware_screening.py \
  --hardware examples/motor_hardware_parallel_325c.toml \
  --configuration examples/motor_demonstrator_original_325c.toml \
  --warm-start outputs/motor_parallel_325c_screening.json \
  --maximum-cycles 1 --exclude-external-losses \
  --output-prefix outputs/microtube_gas_legacy_regression
```

Result: **37.626110566556775 W**, converged, unchanged from the established
37.6261105666 W regression. Original source JSONs and Doty measurement CSVs were
not overwritten.

## Final integration checks

The complete test suite passes: **364 tests**, including 28 focused gas-model
checks (2026-09-16). The historical screening command above still gives exactly
37.626110566556775 W after routing its hardware parsing through the shared loader.
The same command with `motor_hardware_parallel_325c_variable_gas.toml`,
`--plots`, and `--output-prefix outputs/microtube_gas_screening_smoke` also runs
through the variable-property model and saves its domain diagnostics and plots.
That deliberately one-cycle CLI smoke is not converged and is not used as a
performance result; the converged champion comparison is the table above.
