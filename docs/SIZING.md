# Constrained Sizing API

This lower-level API remains useful to existing scripts. New persistent studies
should normally start with [Research](DADA_ENGINE_RESEARCH.md). Its geometry-owned
parameters must not be replaced by independently optimized derived UA/hold-up.

The sizing layer evaluates immutable thermodynamic configurations. It does not
define a preferred design objective or silently repair infeasible simulations.

The API also supports motor objectives `maximize_thermal_efficiency` and
`maximize_motor_power`, the `minimum_motor_power` constraint (`required_power`
in W), and negative angular-speed bounds that stay strictly below zero.
See [MOTOR_OPERATION.md](MOTOR_OPERATION.md) for reservoir and branch conventions.
Cooling constraints require refrigeration operation; motor heat input cannot
satisfy them.

## Evaluation flow

```text
bounded DesignVariable values
    -> DesignPoint
    -> transformed SimulationConfiguration
    -> periodic thermodynamic simulation
    -> performance, diagnostics and validity
    -> objective and constraint margins
```

An evaluation is usable only when its periodic simulation converges. Expensive
periodic evaluations are cached by the exact design-point values.

## Supported design parameters

- small and large swept volumes;
- small and large clearance volumes;
- cold and hot heat-exchanger volumes;
- cold and hot `UA` values;
- all four first-level `CdA` values;
- charge pressure or total gas mass, but not both;
- angular speed.

All values and bounds use SI units.

## Objectives

The current interchangeable scalar objectives are:

- maximize cooling COP;
- minimize total swept volume;
- minimize charging pressure;
- minimize total `UA`.

Maximizing COP does not itself impose a cooling duty. A separate minimum
cooling-power constraint must be included when that operating requirement is
intended.

## Constraints

Constraints return a margin that is non-negative when satisfied. Available
constraints include:

- minimum cooling power;
- maximum pressure;
- maximum temperature;
- maximum absolute mass flow;
- maximum Mach number;
- periodic convergence;
- maximum piston gas-side force for explicitly supplied piston areas.

An unavailable quantity makes its constraint unavailable and therefore
infeasible. A Mach constraint needs an available geometric flow-area diagnostic;
microtubes supply it, whereas a bare orifice CdA alone does not. Pressure inequality
is diagnostic-only.

## SLSQP adapter

`SlsqpSizingOptimizer` maps every bounded physical variable to `[0, 1]` before
calling SciPy SLSQP. The caller must explicitly provide:

- an objective scale;
- one scale for every constraint;
- the penalty used when an objective is unavailable;
- the negative margin used when a constraint is unavailable;
- iteration and termination limits.

These numerical values are not physical assumptions, but making them explicit
prevents hidden weighting from changing the sizing result.

The SLSQP adapter is one local constrained-search strategy. The same
`SizingProblem` can later support other optimizers or Pareto exploration.

## TOML and command line

The TOML file declares a base simulation configuration, variables, one selected
objective, constraints, scales and optimizer termination settings.

After package installation, a configured run is started with:

```console
dada-sizing path/to/sizing_problem.toml
```

Before optimization, evaluate the initial point alone:

```console
dada-sizing path/to/sizing_problem.toml --initial-only
```

This reports periodic status, objective availability and every signed
constraint margin. A nonzero exit status indicates that the point is not a
feasible starting design; it does not mean the thermodynamic integration
necessarily failed.

Each function evaluation may require many complete thermodynamic cycles. The
command therefore does not silently loosen periodic tolerances or replace
failed points with plausible thermodynamic results.

## Pareto comparison

`assess_objective_set` evaluates a caller-supplied collection of design points
against several minimization objectives. `pareto_front` then removes
infeasible, unavailable and dominated points. It introduces no weights and
does not select one point from the resulting tradeoff set.

## Feasibility exploration

Before local optimization, a reproducible Latin-hypercube study can test
whether the configured bounds contain any admissible point:

```console
dada-sizing path/to/sizing_problem.toml \
    --feasibility-samples 20 \
    --seed 42 \
    --csv feasibility.csv
```

The configured initial point is included in addition to the requested sample
count. Progress is printed after each periodic simulation. The summary counts
integration statuses and every unavailable or violated constraint. The CSV
retains all design coordinates, objective values, availability flags,
satisfaction flags and signed margins.

Latin-hypercube sampling improves coverage compared with unstructured random
draws, but it does not prove that no feasible region exists when none of the
finite samples is admissible. Bounds, sample count and seed are user inputs;
the solver does not infer physically meaningful ranges.

## Mechanical boundary

`extract_thermodynamic_loads` exposes:

- small- and large-cylinder pressures;
- gas-side piston forces `P * area`;
- the generalized gas torque
  `P_S * dV_S/dtheta + P_L * dV_L/dtheta`.

This generalized torque is work-conjugate to the imposed angle. It is not a
complete shaft-load prediction. Net piston forces, linkage forces, inertia,
friction, bearing reactions and structural stresses require a future
`MechanicalResponseModel` and additional mechanical data.

## Interchangeable model preparation

See [PLUGGABLE_MODELS.md](PLUGGABLE_MODELS.md) for independent periodic free
kinematics, exchanger construction and the family-specific campaign adapters. The sizing enum is preserved. Free-motion derivative violations return
`invalid_kinematics` with separate diagnostic margins before integration.
Swept-volume/clearance changes update free-motion ranges without changing shape.
