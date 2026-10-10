# Motion-law optimality and certification strategy

## Purpose

The motion-law search should not stop at the statement that no better motion was found. The long-term objective is to distinguish clearly between:

1. a **best-found motion** in a finite numerical search;
2. a motion that is **locally optimal** within a chosen parameterization;
3. a motion that is **globally optimal within a bounded finite-dimensional family**;
4. a motion that is **near-optimal among a much broader class of admissible periodic laws**.

The final two levels require explicit upper bounds on achievable performance, not only repeated searches.

This document defines a practical route toward such a certificate for the experimental motor while preserving the current thermodynamic model and passive-valve rules.

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
- sensitivity around the best-found motion;
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

The optimization objective is not a simple algebraic function of the motion parameters. Each evaluation requires a periodic solution of a nonlinear, topology-varying thermo-hydraulic system with dynamic exchanger walls and passive valves.

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

## Separation of proof scopes

Every optimality statement must specify its scope.

A proof for a finite-dimensional family does not prove optimality among arbitrary periodic motions.

A proof for unconstrained smooth motions does not prove mechanical realizability.

A thermodynamic upper bound that ignores acceleration, force or mechanism constraints may still be useful, but it is an upper bound on an enlarged problem rather than on the final machine.

The following scopes should remain distinct:

- thermodynamic motion optimum;
- smooth bounded-motion optimum;
- mechanically realizable motion optimum;
- four-bar-realisable motion optimum;
- complete machine optimum including mechanical losses.

## From abstract motion to a realizable mechanism

Optimize the selected abstract motion under explicit design and feasibility
constraints. If it has corners or other unrealizable features, construct a
smooth target and evaluate its thermodynamic behavior before fitting a physical
mechanism. Evaluate the **actual six-bar kinematics** in the periodic solver;
geometric RMS fit alone is not a performance result.

[Motor research objectives](MOTOR_RESEARCH_OBJECTIVES.md) defines the
whole-machine questions. Spline and six-bar support are reusable capabilities,
not evidence that synthesis for a selected target is complete.

Velocity, acceleration and mechanical bounds must be explicit before claiming
realizability. Certification and optimal-control bounds described above are
possible research methods, not achieved results.

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

## Interpretation of the four-stage probe

The four-stage law is a thermodynamic design probe, not a proven near-optimal
or mechanically realizable motion. A prototype requires a realizable motion law.
Approximating the linear target does not require retaining its corners or freezing its shape. Its piecewise
constant velocities jump at stage boundaries with no acceleration, inertia,
stress or mechanical-loss penalty. Those jumps are mechanically unrealizable
as written. Any minimum stage duration is a study search bound. Passive-valve
behavior and thermal coupling remain part of the physical evaluation.

No present result proves global optimality, a useful shaft efficiency or
experimental validity of transient gas heat transfer.
