# Motor operation and sign conventions

Reference: [English study, revision 1288](https://dada-engine.org/index.php?title=Thermodynamic_and_Mechanical_Study&oldid=1288),
read on 2026-09-08; especially sections 1, 2.1, 3, 6.5 and 8.

Plots and reports use `H_i`/`H_o` throughout current cycle plots, pressure and
temperature reports, valve and flow labels, moisture displays and coupled
exchanger reports. Historical artifacts retain their original labels; Python and TOML retain
legacy identifiers so existing configurations remain usable.

## Scientific mapping

The study now names the exchangers by physical function. Their branch circulation and
hardware are independent of which external reservoir is connected:

| Study | Existing code and TOML branch name | Hydraulic path |
|---|---|---|
| Heat-in exchanger `H_i` | `C`, `cold_*` | `S <-> C -> L` |
| Heat-out exchanger `H_o` | `H`, `hot_*` | `L <-> H -> S` |

The table shows downstream placement; upstream placement moves the valve to
the other half-link without changing branch circulation. Both placements are
configurable independently in Research.

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

For the external-stream wall family, stream inlet conditions are declared on
each exchanger in the basis; they are not automatically swapped by editing a
reservoir TOML example. See [external streams](EXTERNAL_STREAM_THERMAL_MODEL.md).
The reference-pressure filling policy is independent of the crank origin and
uses maximum simultaneous total gas volume.

## Performance and sizing

Heat is positive into the machine at each thermal boundary. Reservoir models
use heat into gas; wall models use external-stream heat for cycle performance. Successful
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

## Command line and verification

From an uninstalled source checkout:

```console
PYTHONPATH=src python3 -m dada_solver path/to/motor_configuration.toml
PYTHONPATH=src python3 -m pytest
```

After package installation, the first command can be written as
`dada-solver path/to/motor_configuration.toml`.

The regression suite covers reversal and derivative signs for the supported
kinematic configurations, reflected piecewise breakpoints, unchanged origin and
inventory, unchanged unequal UA hardware and valve directions, signed power
and efficiency, motor sizing, and periodic mass/energy closure.

Valve-sequence classification recognizes cyclic rotations of the nominal
four-event sequence. Simultaneous, repeated, missing or overlapping transitions
remain non-nominal when events were actually observed. An unavailable event
sequence is not itself evidence of non-nominal operation. No valve event is
imposed or suppressed by this classification.
