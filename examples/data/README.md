# Experimental reference data

## Nitrogen reference

`doty_1991_nitrogen_reference.csv`: three nitrogen rows at 322 kPa, transcribed
from Tables 1 and 2, printed page 38 (PDF page 8), of Doty et al. (1991),
*The Microtube Strip Heat Exchanger*, Heat Transfer Engineering 12(3).
Source: https://dotynmr.com/download/pubs/1991_HTE_Doty_HeatExchanger.pdf

Units converted to SI: mg/s to kg/s, kPa to Pa, degrees Celsius to kelvin,
percent to fraction. Uncertainty columns reproduce the reported ± values;
no confidence level or statistical distribution is inferred. Temperature
uncertainties are not reproduced. T1–T4 retain original sensor names: shell inlet/outlet T1/T2 and tube
inlet/outlet T3/T4. Pressure loss is tube-side only; UA is overall gas/gas.

These are measured prototype data under their original conditions, not a
manufacturer rating for DADA. Both fluid and inlet temperatures matter; the
temperatures vary between rows. The user authorizes extrapolation when adequate
commercial data are unavailable. Such estimates must state temperature and
property corrections, geometric scaling, flow regime, and uncertainty scenarios.
Do not silently treat other conditions as fixed or substitute these measurements
for the motor UA. See `docs/EXCHANGER_LITERATURE.md` and
`docs/MOTOR_DEMONSTRATOR.md`.

## Helium reference

`doty_1991_helium_reference.csv`: three helium rows, transcribed from Tables 1 
and 2, printed page 38 (PDF page 8), of Doty et al. (1991), same source as above.

The three measurements span 79–213 mg/s helium flow at elevated operating 
pressures (749–825 kPa). Same units and uncertainties convention as nitrogen.
Helium provides a higher-speed test of the microtube model at different fluid 
properties and Reynolds numbers.

**Important:** Doty's own Eq. 7 values overpredict the measured helium pressure
drops by about 25–33 %. An independent project reconstruction using NIST
dilute-gas helium viscosity, ideal-gas density and mean tube temperature gives
approximately 3.319, 2.206 and 5.477 kPa, versus Doty's published 3.0, 2.0 and
5.5 kPa and measured 2.3, 1.5 and 4.4 kPa. The independent reproduction is
within about 11 % of Doty's Eq. 7 values but remains clearly above experiment.
Relative scaling from a measured helium anchor is much closer. No hidden
calibration factor is introduced. See `docs/DOTY_SCREENING.md` for assumptions
and detailed comparisons.
