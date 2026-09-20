# Periodic-map and exact-kinematics reuse diagnostic

Date: 2026-09-20. Follows [Stage 4](SOLVER_ACCELERATION_STAGE4.md).
This is a separate numerical experiment. It changes neither production periodic
iteration nor the kinematics interface, and implements no Newton/Broyden/Anderson
solver. Python remains the default backend. The experiment explicitly selects
the production Numba thermodynamic backend with ordinary Python kinematics.

## Why fixed motion still incurs evaluations

For prescribed motion, geometry, assembly branch, normalization extrema and phase
parameters are fixed for a candidate. The four-bar constructor already prepares
its normalization. A compatible thermal-only candidate can reuse that prepared
geometry. Frequencies, cylinder limits and phases must still be included wherever
they affect the cached quantity: angular derivatives and time derivatives are
not interchangeable.

The functions V(theta) and dV/dtheta nevertheless need values at the adaptive
angles requested by LSODA. Different gas states or thermal parameters generally
change those angles. One precomputed uniform table cannot answer every request
exactly without interpolation or another approximation. This does not imply
that compilation is necessary: exact repeated-angle reuse and simpler scalar
arithmetic should be measured first.

An experimental provider stores only the last angle and its four outputs
(V_S, V_L, dV_S/dtheta, dV_L/dtheta). It returns them only on exactly equal angle
requests. It is attached to one already constructed immutable candidate, after
its motor-direction transform; it never stores gas states or RHS values. There
is no rounded key, interpolated motion, cross-candidate mutable cache, or second
motor-direction transform. The experiment also counts distinct angles, which
adds some instrumentation overhead in both compared paths.

## Exact reuse: measured one-cycle results

Three serial paired repeats alternate enabled/disabled order, with compilation
already warm. Each repeat starts from the same saved tighter-periodic endpoint.
The **complete angle array and all 15 trajectory rows are bitwise identical**
with and without reuse, including the heat/work quadratures.

| One physical cycle | Calls | Consecutive identical angles | Median without reuse (s) | Median with reuse (s) |
|---|---:|---:|---:|---:|
| Production four-stage | 6,461 | 4,207 (65.1%) | 0.261 | 0.232 |
| Smooth four-bar | 4,786 | 3,365 (70.3%) | 0.440 | 0.251 |

The observed cycle reductions are approximately 11% and 43%. These are
instrumented single-cycle measurements, not full-campaign or reporting gains.
They justify testing a small generic exact cache before compiling mechanisms;
they do not yet constitute the full frozen adoption gate for such a cache.
Cross-cycle or cross-candidate cache hit rates were not measured. Reusing an
entire sampled trajectory across different adaptive grids is not established.

The chosen mechanism priorities are slider-crank, four-bar and six-bar. The
existing generic interface should remain usable for each; no compiled provider
is required merely to admit a family. This experiment benchmarks only the current
four-stage and four-bar workloads, not new slider-crank/six-bar implementations.
Solenoids remain deferred, with their eventual numerical cost to be measured.

## Definition of the periodic map

Phi advances one complete physical cycle at fixed forcing phase, with the current
LSODA method, exact kinematic breakpoints and physical checks. Its state is

    (m_S, U_S, m_L, U_L, m_Hi, U_Hi, m_Ho, U_Ho, E_wall_Hi, E_wall_Ho).

Five heat/work quadratures reset to zero each cycle and are not shooting unknowns.
Total gas mass removes one degree of freedom. At each saved anchor x*, let
S = diag(abs(x*)). An orthonormal 10-by-9 matrix B spans the null space of the
mass row after scaling by S. The nine independent coordinates are defined by

    x = x* + S B z,
    G(z) = B^T S^-1 [Phi(x* + S B z) - x*].

Every perturbation therefore conserves the same total mass to floating-point
precision. D G is the representation of D Phi on that tangent space. This avoids
the artificial neutral eigenvalue associated with changing the enclosed charge.
The output projection is only used to calculate the diagnostic matrix; integrated
states are not corrected or projected, and actual mass drift is recorded.

Anchors are the saved Stage-4 tighter-periodic production and smooth endpoints.
At ordinary ODE tolerances their next physical cycle has normalized periodic
error 0.00560 and 0.00434 respectively, well below the original threshold 1.
These are approximate periodic anchors, not exact mathematical fixed points.

Central differences use h = 1e-3, 3e-4 and 1e-4 in scaled tangent coordinates.
Each complete matrix costs 18 physical cycles, plus a shared baseline cycle.
A second process repeats the study with all ODE relative/absolute tolerances
multiplied by 0.1. This tighter ODE setting is experimental only; the production
settings and the original periodic acceptance criterion remain unchanged.

The map excludes all between-cycle extrapolation. In particular, the frozen
benchmark's existing `accelerate_walls=True` schedule is not part of Phi.
Consequently, local eigenvalues must not be used to predict the benchmark's
22/32-cycle counts as though those runs were unaccelerated fixed-point iteration.

## Measured spectrum and mode interpretation

| Anchor | Largest eigenvalue | Second eigenvalue | Next two, approximately | cond(D G - I) |
|---|---:|---:|---|---:|
| Production | 0.950029 | 0.932684 | 0.02559, 0.01082 | 90.3 |
| Smooth four-bar | 0.950092 | 0.931388 | 0.01134, 0.00271 | 75.4 |

The two dominant eigenvalues are real and remain within about 1e-6 across all
successful difference steps and both ODE settings. The much smaller eigenvalues
are more sensitive to differentiation noise; no physical interpretation is
assigned to their tiny real/imaginary parts. Matrix norm changes with step are
approximately 1.2e-4 for production and up to 1.4e-4 for smooth at ordinary ODE
tolerances; tighter smooth integrations reduce these changes to 7.3e-6–1.9e-5.
Condition numbers depend on the stated relative-state coordinate scaling.

The dominant mode couples gas and both wall temperatures in the same direction.
For production, normalizing its cold-wall relative-temperature component to 1
gives approximately:

| Relative-temperature component | S | L | Hi gas | Ho gas | Hi wall | Ho wall |
|---|---:|---:|---:|---:|---:|---:|
| Dominant production mode | 0.945 | 0.768 | 0.759 | 1.010 | 0.758 | 1.000 |

These are eigenvector ratios, not temperature changes in kelvin. Gas temperatures
are reconstructed linearly as delta U/U - delta m/m. The second mode has the two
wall temperatures moving in opposite directions, with coupled gas redistribution.
This supports a thermal-storage interpretation, rather than a purely hydraulic
slow mode or two independent wall-only unknowns.

For an isolated local mode, a factor 0.95003 means about 19.5 cycles per e-fold
and 45 cycles per tenfold reduction. This explains why repeated physical cycles
can relax slowly even when each integration is cheap. It does not certify global
convergence. Indeed, the measured operator norms are about 2.25 and 2.00: some
combinations of perturbations can initially grow in this relative-state metric
even though all computed eigenvalues lie inside the unit circle.

## Failures and numerical limits

Not every difference stencil is admissible. At ordinary ODE tolerances, the
production h=1e-4 stencil raises

    MicrotubeDomainError: Transition flow has no validated hydraulic closure.

With tenfold tighter ODE tolerances, h=1e-4 succeeds but h=3e-4 raises the same
exception. The h=1e-3 production stencil succeeds at both settings. All three
smooth stencils succeed at both settings. Each production failure is recorded
with the initial perturbation and one authoritative Python fallback; it is not
filled with a fabricated matrix column, regularized away, or counted as success.

These errors occur during RHS evaluations; this experiment does not establish
whether an accepted physical trajectory must enter the unsupported domain or
whether an internal LSODA trial excursion causes it. Tighter tolerances alone
are not a general remedy. No physical validity boundary is weakened.

Maximum relative mass drift along successful map cycles is below 6.8e-15 at
ordinary ODE tolerances and 1.7e-14 in the tighter audit. Every stored state is
positive. Existing anchor validity was established by Stage 4; a full final
report/validity replay is not performed for every perturbed diagnostic cycle.
These experiments are not an acceptance gate for optimized candidates.

Including partial failed stencils, the ordinary study attempts 42 production
and 55 smooth map cycles, taking 10.83 and 24.34 s after compilation. The tighter
audit attempts 45 and 55, taking 25.68 and 32.34 s. The exact-angle experiments
are separate. Thus ten forward-difference cycles would be a preliminary estimate,
not the cost of this step-size and tolerance audit. Fresh-process compilation
is also excluded from these map-cycle totals.

## Next experimental gate

The spectrum justifies a guarded periodic-shooting experiment. A damped Newton
reference is useful for accuracy and basin-of-attraction evidence, but rebuilding
nine columns can cost as much as many existing convergence cycles. Anderson or
Broyden should be compared on **complete physical-cycle count and elapsed time**,
including rejected steps and Jacobian preparation, against the current fixed
wall-extrapolation policy as well as plain fixed-point iteration.

Use the same nine mass-conserving scaled coordinates. Require positive admissible
states, bounded corrections, residual decrease or rollback, and a fallback to
normal cycles when a trial integration or difference stencil fails. Both slow
modes and the observed transient amplification argue against unrestricted
componentwise extrapolation. A proposed solution must be followed by a normal
physical cycle satisfying the original periodic criterion, with the unchanged
physical validity and conservation checks. No accelerator has passed that gate
in this diagnostic, and none is enabled by it.

## Reproduction and artifacts

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 examples/diagnose_periodic_map.py \
  --output outputs/periodic_map_diagnostic/new_production_tolerances
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 examples/diagnose_periodic_map.py \
  --output outputs/periodic_map_diagnostic/new_tighter_ode --ode-factor 0.1
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m pytest -q \
  tests/test_periodic_map_diagnostic.py
```

Results are in `outputs/periodic_map_diagnostic/production_tolerances/` and
`tighter_ode/`, including full matrices, complex modes, scales/bases, source and
anchor hashes, cycle records, fallbacks and exact-cache comparisons. The initial
`baseline/` directory is an incomplete first attempt, retained for provenance;
its failure is reproduced and recorded in `production_tolerances/`.
Two focused tests pass: a known nine-mode linear map verifies differentiation
and mass reduction, and adjacent floating-point angles verify exact cache keys.
Stage 4 separately retains its 506-test production validation.

The latest [study export](<Thermodynamic and Mechanical Study — dada-engine.org.html>)
was reviewed as revision 1659. Its generalized pressure-equalized reductions and
expanded phase rationale do not replace the full four-volume equations in this
experiment; see the [study update record](SOLVER_ACCELERATION_IMPLEMENTATION.md).
