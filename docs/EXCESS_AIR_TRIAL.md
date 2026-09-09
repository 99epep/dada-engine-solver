# Excess-air and opposed-Lambda trial

## Scope and air boundary

External air is deliberately oversized to study the internal cycle. External
fan consumption and aerodynamic losses remain excluded from the balance by
user decision. Finite air heat capacity, gas/metal/air resistances and wall
storage remain active, including during flow pauses.

The runner reads optional `[air_sizing]` inputs and sets a fixed external mass
flow for both sources using:

`external_flow = reference_peak_internal_flow * expected_peak_margin * capacity_ratio * cp_gas / cp_air`.

The reference peak is the maximum absolute flow across all four internal ports,
not the cycle-average flow. The previous lower-Lambda result gives 0.00878144
kg/s; a 1.25 allowance and capacity ratio 20 give 0.21942688 kg/s per source.
The achieved minimum capacity ratio is checked against each exchanger's actual
sampled port peak; it is not guaranteed for untested future designs. Update the
reference peak if a new design exceeds the allowance. The check case doubles
external flow to 0.43885376 kg/s.

This suppresses air temperature change, not the air-film resistance. The
constant-Nusselt external film remains a screening assumption. These large
flows through the current narrow external passages exceed the equivalent
laminar-air model's regime (Re about 8472 at the base excess flow). Its reported
pressure loss and fan estimates are not validated design predictions. Treat
this as a thermal-boundary sensitivity experiment, not a fan or build selection.

## Results

All cases use 2472 tubes per exchanger, 63.5 mm tube length, 2 Hz, 25/175 deg C
source inlets and 1 L maximum large-cylinder volume. Gas inventory is preserved
through warm starts. This is not a 100 W functioning motor design.

| Quantity | Excess-air reference | Opposed-Lambda trial | Double-air check |
|---|---:|---:|---:|
| Small/large output angles (degrees) | -80 / +75 | -85 / +40 | -85 / +40 |
| Indicated power (W) | -20.209 | -19.211 | -19.133 |
| External heat input (W) | 131.186 | 121.407 | 121.958 |
| External heat rejection magnitude (W) | 151.375 | 140.691 | 141.108 |
| Mean inlet-to-H_i gas difference (K) | 11.154 | 10.323 | 10.229 |
| Mean H_o gas-to-inlet difference (K) | 12.871 | 11.962 | 11.835 |
| Maximum external-air temperature change (K) | 0.696 | 0.646 | 0.324 |
| Minimum external/internal peak capacity ratio | 25.29 | 23.70 | 47.41 |

At the corresponding valve openings, Lambda_S decreases from 0.046179 to
0.033763, while Lambda_L increases from 0.057903 to 0.097157. These are measured
cycle diagnostics, not imposed valve timing. The proposed mechanical direction
is achieved, but mean H_i gas temperature rises by 0.831 K rather than falling.
Mechanical input decreases by about 1 W. Doubling external flow changes
indicated power by only 0.078 W and mean temperature differences by less than
0.13 K within this model. No positive indicated efficiency is reported.

At periodic operation, rejected heat equals absorbed heat plus mechanical
input for a negative-work cycle. Larger rejection alone does not demonstrate
a superior cooler. The small remaining cycle-storage imbalance is below
0.074 W for these converged cycles; conservation including storage has a
residual below 3e-10 J. These are numerical checks, not physical validation.

## Graphs and interpretation

Each `outputs/motor_excess_air_{reference,opposed,check}_graphs.png` contains
six panels against cycle angle, including all four absolute pressures. Thin,
unlabelled dashed vertical lines show sampled passive-valve pressure-crossing
events on every panel. Close opening/closing events can appear nearly merged.

The separate `_pv.png` figures plot pressure against volume, not their ratio.
Signed area `integral(p dV)` is gas work per cycle: clockwise traversal is
positive; counterclockwise is negative. Arrows indicate progression, dots mark
valve events. Vertical phase lines cannot carry the same meaning on a volume
axis. Displayed loop-work areas use sampled trapezoidal integration; reported
indicated power uses the solver's work quadrature.

## Reproduction and verification

Run from the repository with `PYTHONPATH=src` and
`MPLCONFIGDIR=/tmp/dada-matplotlib`:

```sh
python3 examples/motor_hardware_screening.py --configuration examples/motor_demonstrator_lower_lambda_trial.toml --hardware examples/motor_hardware_excess_air.toml --output-prefix outputs/motor_excess_air_reference --warm-start outputs/motor_lower_lambda_screening.json --exclude-external-losses --accelerate-walls --plots
python3 examples/motor_hardware_screening.py --configuration examples/motor_demonstrator_opposed_lambda_trial.toml --hardware examples/motor_hardware_excess_air.toml --output-prefix outputs/motor_excess_air_opposed --warm-start outputs/motor_excess_air_reference_screening.json --exclude-external-losses --accelerate-walls --plots
python3 examples/motor_hardware_screening.py --configuration examples/motor_demonstrator_opposed_lambda_trial.toml --hardware examples/motor_hardware_excess_air_check.toml --output-prefix outputs/motor_excess_air_check --warm-start outputs/motor_excess_air_opposed_screening.json --exclude-external-losses --accelerate-walls --plots
```

All three cases satisfy the existing gas-plus-wall periodic criterion. The
relevant hardware, air/wall, motor and wall-iteration tests pass (39 tests),
including cylinder-volume bounds for both independent mechanical trials.
