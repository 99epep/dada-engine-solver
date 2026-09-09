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
boundary temperatures differ. The reported uncertainties are retained in the
source CSV; the screening factors below are not statistical confidence bounds.

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
25/175-degree, 100 W demonstrator, then rerun the cycle with explicit candidate
hardware closures. Do not certify useful power from indicated gas power, or
compactness from cylinder volume alone. A robust ranking must survive property,
flow distribution, pulse-response and unmeasured shell-side-loss scenarios.

The selected external-air boundary and geometry sizing foundation are described
in [AIR_SOURCE_EXCHANGERS.md](AIR_SOURCE_EXCHANGERS.md).
