# Motor operation — version 0.2.0

Reference: [English study, revision 1288](https://dada-engine.org/index.php?title=Thermodynamic_and_Mechanical_Study&oldid=1288),
read on 2026-09-08; especially sections 1, 2.1, 3, 6.5 and 8.

Version 0.2.1 uses `H_i`/`H_o` throughout current cycle plots, pressure and
temperature reports, valve and flow labels, moisture displays and coupled
exchanger reports. The current motor plot/report have been regenerated.
Historical artifacts retain their original labels; Python and TOML retain
legacy identifiers so existing configurations remain usable.

## Scientific mapping

The study now names the exchangers by physical function. Their graph and
hardware are independent of which external reservoir is connected:

| Study | Existing code and TOML branch name | Hydraulic path |
|---|---|---|
| Heat-in exchanger `H_i` | `C`, `cold_*` | `S <-> C -> L` |
| Heat-out exchanger `H_o` | `H`, `hot_*` | `L <-> H -> S` |

Legacy branch identifiers remain usable throughout the state vector, geometry,
UA, CdA, valves, validity and exchanger tools. They identify fixed hardware;
in motor operation `C` is connected to the **hot** reservoir and `H` to the
**cold** reservoir. The reservoir configuration keys `cold_temperature` and
`hot_temperature` continue to mean the actual external temperatures.

The first-level conservative state remains `(m_S,U_S,m_L,U_L,m_C,U_C,m_H,U_H)`.
Review against the new study found no need to change its mass balances,
upstream enthalpy transport, finite-UA heat transfer or cylinder boundary work.
The independent-pressure model remains the previously chosen finite-resistance
extension of the study's quasi-pressure-equalized analytical limits.

## Configuration and integration

Set the study's signed angular velocity in rad/s:

```toml
[reservoirs]
cold_temperature = 300.0
hot_temperature = 600.0

[operation]
angular_speed = -0.9520107557799403
```

Positive speed retains existing refrigeration behavior. Negative speed:

1. Reverses both prescribed cylinder laws at the same configured origin.
2. Connects `H_i` to `hot_temperature` and `H_o` to `cold_temperature`.
3. Recalculates the complete periodic state and pressure-driven valve events.

Do not also swap the temperatures or reverse `[kinematics].crank_direction`
when converting an existing refrigeration configuration: that would apply
the corresponding reversal twice. Four-bar `crank_direction` defines the
underlying installed refrigeration geometry; its effective sign is multiplied
by the sign of `[operation].angular_speed`. The mechanical adapter remains
available for animations and geometry sizing.

Integration uses increasing progress `phi = |omega| t` over `[0, 2*pi]`.
For motor operation the study angle is `theta = -phi`, so the integrator uses
`V(-phi)` and `dV/dphi = -V'(-phi)`. Piecewise discontinuities are also reflected
and sorted. Time, heat transfer and passive hydraulics always evolve forward.
`model.angular_speed` is the positive integration speed;
`model.signed_angular_speed` is the signed study speed and must be passed to
`calculate_cycle_performance`.

Stored cycle/event angles and plot axes are progress angles. The generalized
torques returned by the existing mechanical boundary are conjugate to this
progress coordinate: `sum(P*dV/dphi)`. Their product with `|omega|` is gas power.
To express torque conjugate to the study angle, multiply by the study direction.
These remain thermodynamic pressure loads, without friction or inertia.

The origin is preserved, not silently rephased. Harmonic and ideal-piecewise
examples start at `V_L,max` as in the study. Historical four-bar files expose
their own configured angular offsets, some slightly away from the exact large
cylinder maximum. Changing those offsets with pressure-defined filling would
also change inventory; this update deliberately retains that existing behavior.

## Performance and sizing

Heat is positive when received by the gas at each fixed exchanger. Successful
motor operation requires the motor reservoir assignment and
`Q_i > 0`, `Q_o < 0`, `W_cycle > 0`. Then:

```text
thermal_efficiency = W_cycle / Q_i
motor_power = W_cycle * |omega| / (2*pi)
```

`1 + Q_o/Q_i` is an independent periodic first-law check; away from periodic
closure it differs by the change in stored internal energy. The reported
efficiency uses work directly. Specifying negative speed does not guarantee
positive work, convergence, nominal topology or physical validity.

The report exposes signed `heat_in_*`, `heat_out_*` and `gas_power_W`, plus
`thermal_efficiency` and `motor_power_W` when available. COP is unavailable
under the motor reservoir assignment, even if that attempted motor still
consumes work. Signed mechanical input remains `-W_cycle` for compatibility.

The Python `CyclePerformance` retains its historical `cold_heat_per_cycle`,
`hot_heat_per_cycle`, `cooling_power` and `heating_power` fields as signed branch
quantities. In motor calculations use the neutral `heat_in_per_cycle`,
`heat_out_per_cycle`, `heat_in_power` and `heat_out_power` properties instead;
`heat_out_power` has the received-heat sign, opposite the legacy heating power.
Motor reports omit legacy cooling/heating-power labels.

Sizing now supports objectives `maximize_thermal_efficiency` and
`maximize_motor_power`, and the `minimum_motor_power` constraint with a positive
`required_power`. Negative angular-speed design bounds are allowed, but a search
interval cannot contain zero or cross between operation modes. Existing cooling
power/task constraints require actual refrigeration operation, so a motor's
heat input cannot satisfy a cooling demand. These additions provide sizing
interfaces, not an optimized motor design.

## Example and verification

From an uninstalled source checkout:

```console
PYTHONPATH=src python3 -m dada_solver examples/motor_controlled_example.toml
PYTHONPATH=src python3 -m pytest
```

After package installation, the first command can be written as
`dada-solver examples/motor_controlled_example.toml`.

The example reverses the existing ideal-piecewise cooling candidate with
explicit illustrative reservoirs of 600 K and 300 K, dry constant-property gas,
continuous ideal diodes, and 0.03 Pa orifice regularization. It is not a validated
hardware or working-fluid specification.

On 2026-09-08, the reference run converged in five cycles:

| Quantity | Result |
|---|---:|
| Heat input per cycle `Q_i` | 7127.422 J |
| Heat received at heat-out exchanger `Q_o` | -5857.849 J |
| Gas work per cycle | 1269.573 J |
| Mean motor power | 192.362 W |
| Thermal efficiency | 0.178125 (17.8125%) |
| Mass residual | -5.55e-17 kg |
| First-law residual | 5.09e-11 J |

The efficiency is below the 50% two-reservoir Carnot limit. The temperature excursion at `H_o` is 8.38%. Following the user decision,
isothermality is diagnostic-only and no longer fails core validity. Mach remains
unavailable, giving an **indeterminate** verdict, and passive valve chronology
remains non-nominal. Ideal chronology is still a design target; efficiency is
the primary objective. No experimental motor performance is claimed.

Refinement and regularization sensitivity, with the same periodic tolerances:

| Integration rtol | Maximum progress step | Regularization | Motor power | Efficiency |
|---|---:|---:|---:|---:|
| 1e-8 | 0.5 degrees | 0.03 Pa | 192.362137 W | 0.178125115 |
| 1e-9 | 0.25 degrees | 0.03 Pa | 192.362138 W | 0.178125129 |
| 1e-9 | 0.25 degrees | 0.01 Pa | 192.383485 W | 0.178086221 |

Both refined runs converge in four cycles. Reducing regularization by three
changes power by about 0.011% and efficiency by about 0.022%; this is a finite
sensitivity check, not proof of convergence to zero regularization or of exact
valve chronology.

The regression suite covers reversal and derivative signs for all five built-in
kinematic configurations, reflected piecewise breakpoints, unchanged origin and
inventory, unchanged unequal UA hardware and valve directions, signed power
and efficiency, motor sizing, and periodic mass/energy closure. The original
141 tests passed before this update.

An additional diagnostic correction recognizes cyclic rotations of the nominal
four-event valve sequence. Previously a different choice of angular origin
could produce a false `non_nominal` verdict. Simultaneous, repeated, missing,
or overlapping transitions still remain non-nominal. No valve event is imposed
or suppressed by this classification.
