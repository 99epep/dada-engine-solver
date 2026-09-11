# Doty experimental bank screening

## Verified source and scope

The user selected Doty as the first calculation reference on 2026-09-09.
[Doty et al. (1991)](https://dotynmr.com/download/pubs/1991_HTE_Doty_HeatExchanger.pdf),
printed pages 37–38 (PDF pages 7–8), describes three modules of 103 tubes each.
Tube internal diameter is 0.33 mm, length 127 mm and wall thickness 0.1524 mm.
The complete bank measures 190 x 30 x 35 mm and weighs 0.3 kg. Its bounding-box
volume is 0.1995 litres; this is not gas hold-up. Additional manifolds, spacing,
insulation and external circuits are excluded from replicated envelope totals.

In the experiment, gas first passes through the tubes (inlet T3, outlet T4),
then through a pressure-reduction valve and heater, and returns through the
shell (inlet T1, outlet T2). The heat-capacity rates of the two streams are equal.
Reported UA is overall gas/gas conductance. Table 2 gives **tube-side** pressure
loss, not combined shell/tube loss. The tube and shell pressures must not be
assumed equal. These distinctions prevent direct reuse as a gas/liquid rating.

The three transcribed nitrogen anchors at 322 kPa remain separate: their
boundary temperatures differ. The three helium measurements span 749–825 kPa
and 79–213 mg/s flow. The reported uncertainties are retained in the source CSV;
the screening factors below are not statistical confidence bounds.

## Implemented calculation

`dada_solver.exchangers.doty.scale_bank` replicates complete banks in parallel,
assuming equal flow distribution on both sides and identical inlet conditions:

- `r = total_mass_flow / (bank_count * reference_mass_flow)`;
- `UA = bank_count * reference_UA * r**exponent * condition_factor`;
- `tube_dp = reference_dp * r * viscosity_ratio / density_ratio`.

The pressure scaling follows the laminar tube relation in the paper's Eq. 7,
anchored to measured pressure loss. It is not validated after transition or for
large pressure changes. Density is representative tube-side density. The UA
flow exponent and overall condition factor are explicit hypotheses. The
measured gas/gas UA cannot determine a gas/liquid resistance split uniquely.
The model exposes additional pressure loss separately; zero means omitted,
not proven absent. Shell-side performance remains unresolved.

No geometry resizing, real-gas correction, wall dynamics or motor-cycle coupling
has been implemented here. In particular, this module does not replace the
cycle's constant UA or orifice CdA. Zero or reverse flow is rejected rather than
inventing a steady-map thermal closure during idle intervals or local reflux.

## Helium measurements and hydraulic validation

### Absolute Poiseuille comparison (Doty Eq. 7)

Doty's Table 2 gives both measured helium tube-side pressure drops and values
calculated from Eq. 7:

| Flow (mg/s) | Measured dP (Pa) | Doty Eq. 7 (Pa) | Doty overprediction |
|---|---:|---:|---:|
| 117 | 2300 | 3000 | +30.4 % |
| 79 | 1500 | 2000 | +33.3 % |
| 213 | 4400 | 5500 | +25.0 % |

Doty explicitly states that the helium pressure drops are somewhat lower than
expected and notes minor tube-side flow maldistribution as a possible
explanation.

The project also performs an independent reproduction of Eq. 7 rather than
using the Table 2 calculated column as an input. The assumptions are:

- 309 tubes, 0.33 mm internal diameter and 127 mm length;
- representative tube temperature `(T3 + T4) / 2`;
- ideal-gas helium density at the reported pressure;
- dilute-gas helium viscosity linearly interpolated from the NIST values
  21.0 uPa.s at 325 K and 22.1 uPa.s at 350 K.

NIST source:
https://www.nist.gov/pml/sensor-science/fluid-metrology/database-thermophysical-properties-gases-used-semiconductor-9

The NIST table attributes these viscosity values to Hurly and Moldover (2000)
and gives an estimated viscosity uncertainty of 0.1 %. Doty's paper does not
state enough detail about the density/viscosity evaluation used for Table 2 to
require exact numerical reproduction of its calculated column.

| Flow (mg/s) | Measured dP (Pa) | Doty Eq. 7 (Pa) | Independent Eq. 7 (Pa) | Independent vs measured | Independent vs Doty |
|---|---:|---:|---:|---:|---:|
| 117 | 2300 | 3000 | 3319 | +44.3 % | +10.6 % |
| 79 | 1500 | 2000 | 2206 | +47.1 % | +10.3 % |
| 213 | 4400 | 5500 | 5477 | +24.5 % | -0.4 % |

Thus `MicrotubeBank.laminar_tube_loss()` is consistent with the algebraic form
of Doty's Eq. 7 and reproduces the published calculated values to about 11 %
under an explicit independent property convention. It does **not** validate
the absolute helium pressure drop experimentally. No empirical multiplier is
introduced to force agreement.

### Relative scaling between helium measurements

Using the 117 mg/s measurement as a single anchor and applying
`dp ~ mass_flow * viscosity / density`, with the same NIST viscosity
interpolation and ideal-gas representative density:

| Flow (mg/s) | Measured dP (Pa) | Predicted from 117 mg/s anchor (Pa) | Relative error |
|---|---:|---:|---:|
| 117 | 2300 | 2300 (anchor) | - |
| 79 | 1500 | 1529 | +1.9 % |
| 213 | 4400 | 3796 | -13.7 % |

The relative scaling is therefore substantially more successful than the
absolute prediction: it tracks the 79 mg/s point closely and remains within
about 14 % at 213 mg/s. This supports use of the laminar scaling as screening
physics while retaining an explicit uncertainty on absolute loss.

## Reproducible illustrative comparison

Run `PYTHONPATH=src python3 examples/doty_bank_screening.py` to regenerate
`outputs/doty_bank_screening.csv`. It uses the 0.84 g/s, 5.3 W/K, 9.5 kPa anchor,
replicates 1, 4 and 16 banks, and compares rectangular duty fractions 1, 1/2
and 1/4 at the same average flow per bank. Exponents 0 and 0.5 are sensitivity
scenarios, not fitted correlations. Properties and condition factors stay at
reference values. These are not the flows or temperatures of a solved motor.

For four banks under original steady conditions, the estimate is 21.2 W/K,
9.5 kPa tube loss, 0.798 litres summed bank envelope and 1.2 kg bank mass.
At half duty, keeping the same mean flow requires twice the active flow:
active tube loss becomes 19 kPa, while active UA is 21.2 or 30.0 W/K under the
two exponent assumptions. This is extrapolation beyond the selected anchor.
Pressure loss must not be computed from mean flow alone.

The steady surrogate has no frequency dependence at a fixed pulse duty and
amplitude. It therefore cannot distinguish 2 Hz from 10 Hz. The generated
values are active-flow estimates, not cycle-averaged heat transfer. A future
cycle comparison must use simultaneous port states, preserve reflux reporting,
and assess whether wall response or gas residence time defeats this assumption.

## Next calculation boundary

Use this reference to screen plausible thermal/hydraulic requirements for the
25/325-degree, 100 W demonstrator, then rerun the cycle with explicit candidate
hardware closures. Do not certify useful power from indicated gas power, or
compactness from cylinder volume alone. A robust ranking must survive property,
flow distribution, pulse-response and unmeasured shell-side-loss scenarios.

The selected external-air boundary and geometry sizing foundation are described
in [AIR_SOURCE_EXCHANGERS.md](AIR_SOURCE_EXCHANGERS.md).
