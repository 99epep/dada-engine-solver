# Motor operation and sign conventions

Reference: [English study, revision 1288](https://dada-engine.org/index.php?title=Thermodynamic_and_Mechanical_Study&oldid=1288),
especially sections 1, 2.1, 3, 6.5 and 8.

Plots and reports use `H_i`/`H_o` throughout current cycle plots, pressure and
temperature reports, valve and flow labels, moisture displays and coupled
exchanger reports. Python and TOML use the hydraulic port field names below.

## Scientific mapping

The study names the exchangers by physical function. Their branch circulation and
hardware are independent of which external reservoir is connected:

| Study | Existing code and TOML branch name | Hydraulic path |
|---|---|---|
| Heat-in exchanger `H_i` | `cold_*` | `S <-> H_i -> L` |
| Heat-out exchanger `H_o` | `hot_*` | `L <-> H_o -> S` |

The table shows downstream placement; upstream placement moves the valve to
the other half-link without changing branch circulation. Both placements are
configurable independently in Research.

H_i and H_o identify fixed hardware. In motor operation H_i connects to the
hot reservoir and H_o to the cold reservoir. The reservoir configuration keys
`cold_temperature` and `hot_temperature` identify actual external temperatures.
The conservative state is `(m_S,U_S,m_L,U_L,m_Hi,U_Hi,m_Ho,U_Ho)`.

## Configuration and integration

Set the study's signed angular velocity in rad/s:

```toml
[reservoirs]
cold_temperature = 300.0
hot_temperature = 600.0

[operation]
angular_speed = -1.0
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

The configured crank origin is preserved. Harmonic and ideal-piecewise laws
start at `V_L,max`; four-bar laws use their declared angular offsets.
With pressure-defined filling at the initial angle, changing the origin can
change inventory. The reference-pressure policy instead uses maximum
simultaneous total gas volume.

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
consumes work. Signed mechanical input is `-W_cycle`.

The Python `CyclePerformance` stores `cold_heat_per_cycle`,
`hot_heat_per_cycle`, `cooling_power` and `heating_power` fields as signed branch
quantities. In motor calculations use the neutral `heat_in_per_cycle`,
`heat_out_per_cycle`, `heat_in_power` and `heat_out_power` properties instead;
`heat_out_power` has the received-heat sign, opposite the heating-power sign convention.
Motor reports use indicated motor power and signed branch heat.

Sizing supports objectives `maximize_thermal_efficiency` and
`maximize_motor_power`, and the `minimum_motor_power` constraint with a positive
`required_power`. Negative angular-speed design bounds are allowed, but a search
interval cannot contain zero or cross between operation modes. Existing cooling
power/task constraints require actual refrigeration operation, so a motor's
heat input cannot satisfy a cooling demand. These are sizing interfaces, not an optimized motor design.

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
