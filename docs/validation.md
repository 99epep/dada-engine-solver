# Validation Record

This document records numerical and analytical checks of the solver. It does
not validate any physical DADA machine design or any example parameter set.

## Motor update, 2026-09-08

Version 0.2.0 adds signed motor operation against English study revision 1288.
The motor example, first-law checks, refinement and regularization sensitivity
are recorded in [MOTOR_OPERATION.md](MOTOR_OPERATION.md). It produces about
192.36 W at 17.81% thermodynamic efficiency, with an indeterminate
physical-validity verdict (Mach unavailable) and non-nominal valve topology.
Isothermality is now a diagnostic only, following the user decision.
The update also tests origin-independent classification of the nominal event
sequence and prevents motor heat input from satisfying refrigeration constraints.

## Analytical regression checks

The automated suite currently checks:

- ideal-gas reconstruction from conservative mass and internal energy;
- the closed adiabatic invariant `P * V**gamma`;
- the Appendix A.5 outflow-only adiabatic donor relation;
- exact instantaneous cancellation of internal mass and enthalpy flows;
- exact conservative-state continuity at nonzero valve events;
- zero heat transfer when gas and reservoir temperatures are equal;
- choked and unchoked compressible-orifice limits;
- upstream selection for bidirectional links;
- check-valve directionality and physical hysteresis;
- periodic convergence of a static controlled case;
- global cycle mass and energy balances.

## Mobile adiabatic integration

A controlled model with isolated cylinders, zero `UA`, closed valves and
harmonic imposed volumes was integrated over a complete revolution. Both
cylinders preserve `P * V**gamma`, and the final conservative state recovers
the initial state. A refined integration has a smaller invariant error than a
coarse integration.

## Periodic refinement audit

The controlled harmonic example was solved twice with identical physical
inputs and periodic convergence settings.

| Setting | Reference | Refined |
|---|---:|---:|
| Integration relative tolerance | `1e-9` | `2e-10` |
| Integration absolute tolerance | `1e-11` | `2e-12` |
| Maximum angular step | `0.5 deg` | `0.2 deg` |
| Cycles to convergence | 21 | 21 |

Observed refinement differences:

| Quantity | Difference |
|---|---:|
| Maximum relative final-state difference | `5.28e-9` |
| Relative cold-heat difference | `1.06e-8` |
| Relative hot-heat difference | `1.70e-8` |
| Relative gas-work difference | `5.59e-9` |
| Maximum valve-event angle difference | `2.90e-8 rad` |

Both solutions recover the same four-event sequence. The example remains a
non-refrigerating, non-nominal and physically invalid result under its selected
validity thresholds. Refinement agreement establishes numerical consistency;
it does not override those physical diagnostics.

## Controlled sizing preflight

The initial point in `examples/sizing_controlled_example.toml` reaches periodic
convergence but is infeasible. Its external 20-minute water-freezing constraint
is calculated from the configured thermal load rather than a copied power. The
resulting cooling-power margin is approximately `-410.90 W`, because the point
has negative cooling power while the ideal task requires approximately
`347.58 W`. Its pressure-equalization margin is approximately `-0.0331`.

The point requires approximately `92.8 W` mechanical input and therefore
satisfies the configured `150 W` nominal human-power limit. Pressure,
temperature, mass-flow and periodic-convergence constraints are also satisfied.
No full optimization result is claimed from this example.
# Runtime profiling

The verified high-COP candidate was profiled with CPython `cProfile` and the
instrumented Radau integration. Five cycles required approximately 106 seconds
under profiling, 166,136 solver-reported right-hand-side evaluations, 6,062
numerical Jacobian evaluations, and about 45 million Python calls. Profiling
overhead means this elapsed time is not a production benchmark.

The cumulative profile attributed about 102 seconds to Radau stepping, 81
seconds to right-hand-side callbacks, 62 seconds to `ThermodynamicModel`
evaluation, and 31 seconds to numerical Jacobian construction. These times
overlap: Jacobian construction obtains its columns by calling the same right-
hand side repeatedly. The reference compressible-orifice calculations
accounted for about 9 cumulative seconds and were not the dominant cost.

Only eight integration segments occurred in the final cycle and the observed
valve chronology contained the expected four events. The cost is therefore
not caused by event chatter or by a large number of physical phases. It is
caused primarily by the number of implicit solver steps, nonlinear iterations,
and finite-difference Jacobian calls within those segments. The square-root
orifice behavior near pressure equality is a likely mathematical source of
poor conditioning, but this causal interpretation still requires a controlled
regularity study and is not treated as proven.

One initial-cycle comparison at the same tolerances gave approximately 15.9 s
for dense numerical Jacobians and 14.9 s when Radau was supplied the exact
topology-independent Jacobian sparsity pattern. Reported right-hand-side calls
fell from 41,136 to 30,613, Jacobian evaluations from 1,351 to 1,057, and LU
decompositions from 9,232 to 6,484. Sparsity is therefore retained, but the
speed gain is modest because the system has only eleven augmented states.

BDF completed the same initial cycle in about 17.3 s. LSODA completed it in
about 9.6 s, but the subsequent periodic fixed-point run did not finish within
the bounded observation time. RK45 did not finish the initial cycle in one
minute. No integration method is changed solely from these single-cycle
measurements; Radau remains the verified reference.

The next performance work should prioritize, in order:

1. explaining the very large number of steps near pressure equalization;
2. supplying or deriving a reliable Jacobian for the smooth portions;
3. reducing repeated checked dataclass and NumPy allocations in the callback;
4. evaluating volumes and their derivatives together;
5. reconsidering how the three one-way quadratures participate in error
   control;
6. only then considering compilation or another language.

The simulation report now exposes elapsed integration time, right-hand-side
and event-function evaluations, final-cycle segment count, Jacobian
evaluations, and linear decompositions. These counters allow every future
optimization to be checked against numerical trajectories and conservation,
not wall-clock time alone.

### Segment localization and allocation reductions

Per-segment statistics were subsequently added. On the initial cycle of the
high-COP candidate, the two long transfer segments dominated:

- 54-180 degrees: 10,813 RHS evaluations and about 5.22 s;
- 234.022-360 degrees: 11,322 RHS evaluations and about 5.53 s.

Together they represented roughly 72% of measured integration time. Individual
event-ending segments cost at most about 0.33 s. This confirms that root finding
and event count are not the dominant cost; sustained pressure-driven transfer
within each open topology is.

Three implementation-only optimizations were then applied without changing an
equation: a combined volume/derivative kinematic evaluation, a single joint
temperature/pressure reconstruction, and scalar validation of the eight-state
dataclass instead of repeated tiny NumPy reductions. With Jacobian sparsity,
the same initial cycle fell from the measured dense baseline of about 15.9 s
to about 12.6 s, a roughly 21% reduction. Solver counters remained at 30,613
RHS evaluations, 1,057 Jacobians, and 6,484 LU decompositions after the
allocation changes, demonstrating that these gains affect callback cost rather
than numerical trajectory.

An optimized one-cycle profile still recorded about 9.98 million Python calls.
Radau stepping dominated, followed by thermodynamic evaluation. Sparse
Jacobian construction itself remained costly because SciPy constructs and
factorizes sparse objects for a system with only eleven augmented states. A
fully or partly analytical dense Jacobian is therefore a more promising next
experiment than additional sparsity tuning.

Use `--integration-profile` to include every final-cycle segment in the normal
simulation report:

```console
dada-solver configuration.toml --integration-profile
```

### Scalar four-bar closure

Profiling the E0 continuous-diode case showed that small NumPy arrays and
`numpy.linalg.solve` on a two-by-two velocity-closure matrix dominated the
Python-side cost. The four-bar position, velocity closure, rocker/coupler output
and slider projection were rewritten as the same explicit scalar equations.
No physical or numerical approximation was introduced.

The profiled initial cycle fell from approximately 19.1 s and 8.68 million
function calls to 6.9 s and 3.62 million calls, a 2.76-fold profile speed-up.
The complete four-cycle E0 periodic run fell from approximately 49.6 s to
19.3 s wall time on the same machine. Cooling power changed from 177.85357 W
to 177.85409 W and COP from 2.0189601 to 2.0189797 under the configured loose
exploratory integration tolerance; conservation remained near machine
precision. This difference is numerical-tolerance scale, not a physical model
change. The full 121-test suite passed after the rewrite.

### Extensive and speed similarity

Candidate 001 was scaled to the 20-minute ideal water-load power using the exact
first-level similarity factor 3.471579. The independent rerun converged in four
cycles and returned 347.582 W cooling, 103.465 W thermodynamic input and COP
3.359416. The reference values multiplied by the factor agree within integration
tolerance. This validates the model-level similarity transformation, but not the
manufacturability of its required 433.95 W/K exchanger conductances or enlarged
effective flow areas.
