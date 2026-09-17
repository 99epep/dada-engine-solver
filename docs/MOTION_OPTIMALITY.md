# Motion-law optimality and certification strategy

## Purpose

The motion-law search should not stop at the statement that no better motion was found. The long-term objective is to distinguish clearly between:

1. a **best-found motion** in a finite numerical search;
2. a motion that is **locally optimal** within a chosen parameterization;
3. a motion that is **globally optimal within a bounded finite-dimensional family**;
4. a motion that is **near-optimal among a much broader class of admissible periodic laws**.

The final two levels require explicit upper bounds on achievable performance, not only repeated searches.

This document defines a practical route toward such a certificate for the DADA motor while preserving the current thermodynamic model and passive-valve rules.

## What an optimality claim means

Let the admissible motion class be `K`, and let

```text
eta(k)
```

be the converged indicated thermal efficiency obtained from the complete periodic thermo-hydraulic simulation for motion `k` in `K`, subject to all declared validity and engineering constraints.

A candidate `k*` is **globally optimal in K** if:

```text
eta(k*) = sup eta(k),  k in K
```

For engineering purposes, an **epsilon-optimal certificate** is sufficient:

```text
sup eta(k) - eta(k*) <= epsilon
```

This is much stronger than reporting that many searches converged to the same region.

The value of `epsilon`, the admissible motion class, the constraints, the thermodynamic model, the exchanger model, the gas inventory and all numerical tolerances must be stated with the claim.

## Level 1 — best-found numerical evidence

Sobol campaigns, local refinement and repeated searches from different seeds provide useful evidence but do not prove global optimality.

Evidence should include:

- the best feasible objective value;
- the best competing regions found;
- distance from active parameter bounds;
- sensitivity around the champion;
- repeatability across independent seeds;
- numerical convergence and periodic-state tolerances;
- all active physical constraints.

The correct wording at this level is:

> best found under the stated search budget and parameter bounds

and not:

> globally optimal motion.

## Level 2 — parameterization-independence evidence

A strong intermediate test is to solve the motion problem with progressively richer and structurally different representations.

Examples include:

- the reduced two-extrema family;
- the fully independent one-sided-extrema family;
- periodic `FreeKinematics` splines with increasing control-point count;
- alternative smooth bases such as truncated Fourier series.

If richer representations converge toward the same efficiency and the same qualitative motion, this is evidence that the result is not an artifact of one chosen parameterization.

Useful convergence checks are:

```text
eta_N -> eta_infinity
```

and

```text
||V_S,N - V_S,N+1|| -> 0
||V_L,N - V_L,N+1|| -> 0
```

with analogous checks on angular derivatives.

This still does not prove global optimality, but it can establish that the remaining uncertainty from motion representation is small.

## Level 3 — certified global optimum in a finite-dimensional family

For a bounded family such as an 11-variable motion family, a true global certificate is possible in principle with deterministic global optimization.

The preferred concept is branch-and-bound:

1. start from the complete bounded parameter box;
2. compute a rigorous upper bound on efficiency for every box;
3. retain the best feasible lower bound from evaluated candidates;
4. discard any box whose upper bound cannot beat the current best;
5. subdivide the remaining boxes;
6. stop when the global upper-bound gap is below the chosen `epsilon`.

At termination:

```text
eta_best <= eta_global <= eta_upper
```

with

```text
eta_upper - eta_best <= epsilon
```

That is a numerical proof of epsilon-global optimality **within that bounded motion family**.

### Main difficulty

The DADA objective is not a simple algebraic function of the motion parameters. Each evaluation requires a periodic solution of a nonlinear, topology-varying thermo-hydraulic system with dynamic exchanger walls and passive valves.

Therefore a useful branch-and-bound implementation needs certified bounds not only on the kinematic law but also on the resulting periodic thermodynamic solution.

A naive interval evaluation of the complete solver would probably be too loose and too expensive.

A practical certification effort should therefore begin only after the best-motion region has been narrowed substantially.

## Level 4 — upper bounds beyond a chosen parameterization

The most valuable scientific result would be an upper bound applying to a broad class of admissible periodic motions rather than only to one finite parameterization.

The motion laws can be viewed as controls:

```text
V_S(theta), V_L(theta)
```

or equivalently their angular derivatives, subject to periodicity, stroke, velocity, acceleration and mechanical-feasibility constraints.

The conservative gas and wall states form the dynamical system.

This converts the search conceptually into an optimal-control problem.

Two complementary approaches are relevant.

### Necessary optimality conditions

Pontryagin-type or calculus-of-variations analysis may identify necessary conditions on an optimal motion.

Such analysis could explain features observed numerically, for example:

- location of extrema;
- rapid versus slow sectors;
- asymmetric approach and departure from extrema;
- whether some velocity or acceleration bounds should be active;
- whether certain thermo-hydraulic switching events should coincide with kinematic features.

These conditions are scientifically useful even when they do not provide a global proof.

### Relaxed upper bound

A stronger route is to define a relaxed problem whose feasible set contains all physically admissible motions:

```text
K_physical subset K_relaxed
```

Then:

```text
sup eta(k), k in K_physical
<=
sup eta(k), k in K_relaxed
```

If a computable relaxed upper bound is close to the efficiency of the best physical motion, the gap itself becomes a near-optimality certificate.

For example, if the best admissible motion gives:

```text
eta_best = 21.10%
```

while a rigorously justified relaxation proves:

```text
eta <= 21.20%
```

for every admissible motion in the declared class, then the remaining optimality gap is at most:

```text
0.10 percentage point
```

That statement is substantially stronger than an arbitrarily large numerical search.

## Separation of proof scopes

Every optimality statement must specify its scope.

A proof for an 11-parameter family does not prove optimality among arbitrary periodic motions.

A proof for unconstrained smooth motions does not prove mechanical realizability.

A thermodynamic upper bound that ignores acceleration, force or mechanism constraints may still be useful, but it is an upper bound on an enlarged problem rather than on the final machine.

The following scopes should remain distinct:

- thermodynamic motion optimum;
- smooth bounded-motion optimum;
- mechanically realizable motion optimum;
- four-bar-realisable motion optimum;
- complete machine optimum including mechanical losses.

## Agreed next project sequence

1. Finish the current piecewise-linear machine and temperature exploration,
   using the production gas model and explicit design/feasibility conventions.
2. Take the selected best linear candidate and construct a smooth C2 motion
   approximating it.
3. Allow limited thermodynamic reoptimization of that C2 law to measure whether
   smoothing loses or gains efficiency.
4. Synthesize a six-bar mechanism approximating the selected C2 target.
5. Evaluate the **actual six-bar kinematics** in the periodic thermodynamic
   solver; geometric RMS fit alone is not a performance result.

The [four-stage chronology](FOUR_STAGE_OPTIMIZATION.md) records the implemented
experiments; [MOTOR_RESEARCH_OBJECTIVES.md](MOTOR_RESEARCH_OBJECTIVES.md) defines
the planned whole-machine temperature map. Existing spline and six-bar support
are reusable capabilities. They do not establish that C2 smoothing, its
reoptimization or mechanism synthesis for the newly selected four-stage target
has already been completed.

Velocity/acceleration and mechanical bounds must be made explicit before
claiming realizability. Certification, optimal-control bounds and branch-and-
bound described above remain possible later research, not the immediate
implementation plan or an achieved result.

## Reporting language

Use the following terminology consistently.

| Evidence available | Recommended wording |
|---|---|
| finite Sobol/local search only | best found |
| repeated convergence across parameterizations | robust best-found motion |
| local first/second-order conditions only | locally optimal in the stated family |
| exhaustive/certified bounded global search | epsilon-globally optimal in the stated finite-dimensional family |
| rigorous relaxed upper bound close to best admissible motion | certified near-optimal motion in the stated admissible class |
| no explicit upper bound | do not claim global optimality |

## Current status (2026-09-17)

The project remains at the **best-found / parameterization-comparison** stage.
Shared seven-parameter four-stage timing has been explored and finely refined.
An 11D independent S/L timing search retained the shared seed; a focused 4D
relative-event test then found a small gain from sub-degree RMS separation.
Shared timing is therefore a useful reduced approximation, not exact physics.

An independent local climb from 150/30/150/30 deg found another, more balanced
timing basin. Those comparisons used the legacy gas model. The current physical
9D refinement varies shared motion and independent H_i/H_o tube counts with the
production variable-property model; its committed checkpoint is incomplete.
See [the report-based chronology](FOUR_STAGE_OPTIMIZATION.md) for exact artifacts,
parameters, filling rules, budgets and model generations. Efficiencies across
those generations must not be ranked as optimization improvements.

The four-stage law is still only a thermodynamic design probe. Its piecewise
constant velocities jump at stage boundaries with no acceleration, inertia,
stress or mechanical-loss penalty. Those jumps are mechanically unrealizable
as written, and the 2% minimum duration is only a search bound. Passive-valve
behavior and thermal coupling remain part of the physical evaluation.

No present result proves global optimality, a useful shaft efficiency or
experimental validity of transient gas heat transfer.
