# Doty experimental bank screening

## Source experiment and measurement boundary

[Doty et al. (1991)](https://dotynmr.com/download/pubs/1991_HTE_Doty_HeatExchanger.pdf),
printed pages 37–38 (PDF pages 7–8), describes three modules of 103 tubes each (309 tubes total).
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
the screening factors below are not statistical confidence bounds. The SI
transcriptions are [nitrogen](../examples/data/doty_1991_nitrogen_reference.csv)
and [helium](../examples/data/doty_1991_helium_reference.csv). Their columns retain
the source sensor names, measurements and reported uncertainties.

## Screening model

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

The API requires a positive integer bank count, positive finite flow and property
ratios, `0 <= ua_flow_exponent <= 1`, and nonnegative additional pressure loss.
These input bounds are not experimental validity guarantees.

No geometry resizing, real-gas correction, wall dynamics or cycle coupling
has been implemented here. In particular, this module does not replace the
cycle's constant UA or orifice CdA. Zero or reverse flow is rejected rather than
inventing a steady-map thermal closure during idle intervals or local reflux.

## Helium scaling check

Doty's Table 2 reports tube-side measured losses of 2300, 1500 and 4400 Pa at
117, 79 and 213 mg/s, respectively. Its Eq. 7 calculated values are 3000, 2000
and 5500 Pa. Doty notes flow maldistribution as a possible explanation for the
lower measured losses; it is not used here as a fitted correction.

The project also performs an independent reproduction of Eq. 7 rather than
using the Table 2 calculated column as an input. This deliberately retained reconstruction uses a historical property convention,
not the current production helium transport closure. Its assumptions are:

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

`MicrotubeBank.laminar_tube_loss()` reproduces the algebraic form of Doty's
Eq. 7, but the independent absolute predictions remain above the measured
helium losses. No empirical multiplier forces agreement. The detailed absolute
error table and its interpretation belong to
[exchanger validation](EXCHANGER_VALIDATION.md#doty-hydraulic-comparison).

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

## Extrapolation rules

Extrapolations must retain the source conditions and disclose gas-property and
geometry scaling, Reynolds-number changes, pulse waveform/duty and sensitivity scenarios; they are
estimates, not measurements or manufacturer ratings. Frequency alone does not
justify transferring a pulsating-flow correlation: waveform, Reynolds number,
thermal boundary conditions and reflux must also be compared.

Pressure loss must not be inferred from mean flow alone. The steady surrogate
has no frequency dependence at fixed pulse duty and amplitude; its estimates
are not cycle-averaged heat transfer. Assess simultaneous port states, wall
response, gas residence time and possible reflux before using a steady map in
a pulsed application.

## Limits and related references

The [example script](../examples/doty_bank_screening.py) illustrates API sensitivity
calculations; its outputs are not measurements or model validation. Unit tests
in `tests/test_doty.py` check source transcription, scaling and the explicit
helium reconstruction without running a cycle campaign.

- [Exchanger validation](EXCHANGER_VALIDATION.md): absolute discrepancies and
  the absence of independent shell-side thermal validation.
- [Microtube model](MICROTUBE_GAS_MODEL.md): current production equations,
  transport and domains, separate from the historical interpolation above.
- [External streams](EXTERNAL_STREAM_THERMAL_MODEL.md): finite-stream and wall
  boundaries; they are not supplied by `scale_bank`.
