# Interchangeable kinematics and exchanger models

## Scope and construction

Four-bar motion is not a fundamental assumption of the DADA cycle. Microtubes
are the currently best-supported geometry family in the connected motor
experiments, not a fundamental thermodynamic assumption. The conservative four
gas volumes, hydraulic graph, passive valves, heat/work signs and equations
are unchanged by this refactor. No mechanical efficiency is introduced.

`MachineDesign` in `dada_solver.machine` composes a `SimulationConfiguration`
with an optional directly injected `KinematicsModel` and independent
`ExchangerModel` designs for both sides. Injected motion uses study angle;
the factory skips family construction and applies motor reversal exactly once.
The existing TOML selector remains the default when no implementation is supplied. Its `build()`
method calls the existing factory and, if provided, the exchanger connector.
There is no loader, registry, entry-point discovery or dependency-injection
framework. A family is selected explicitly at construction time, outside the
integration loop. Family choice is fixed for an optimization campaign; compare
different families in separate campaigns, not with a continuous family variable.

The existing modules stay in place to preserve public imports. The new files
are `free_kinematics.py`, `machine.py`, `exchangers/base.py` and
`exchangers/microtube.py`. No large source tree was relocated for appearance.

## Common kinematics interface

`kinematics.KinematicsModel` is a structural `typing.Protocol` exposing:

- small/large cylinder volume and first angular derivative;
- small/large `CylinderVolumeLimits`, including swept volume;
- optional physical strokes, represented by `None` for mathematical volume laws;
- derivative-discontinuity breakpoints for integration.

The existing name `VolumeKinematics` remains an alias. A combined
`cylinder_volumes_and_derivatives` method remains an optional fast path.
`SmoothKinematicsModel` describes optional analytical second derivatives;
thermodynamics requires only the first derivative. Diagnostics specific to a
family remain on that implementation, not in the thermodynamic balance.

The canonical cycle domain is theta in radians, one cycle `[0, 2*pi)`. Free
laws wrap all finite angles periodically, including negative angles. Definitions
retain the existing study-angle convention: the factory applies `V(-theta)`
and the first-derivative sign reversal once for negative-speed motor operation.
Second derivatives reverse their argument but not their sign. The physical
time coordinate always advances. A free law need not put L at maximum volume
at zero; initial filling uses its actual volumes at the declared origin.

`dynamics.py`, `integration.py`, `periodic.py` and the wall integrator do not
import a concrete four-bar or free-motion class. Family selection is confined
to configuration/factory boundaries and implementation-specific tools.

## Four-bar compatibility

`four_bar.FourBarKinematics` aliases the existing
`SharedCrankFourBarVolumeKinematics` class. Existing circle closure, assembly
branches, finite slider rods, normalization, stroke evaluation and singularity
checks are reused without rewriting their mathematics. Existing configuration
names, animation and sizing code keep the same concrete objects. The factory's
motor-direction handling for the shared-crank type is unchanged.

A regression test compares full volume/first-derivative trajectories and slider
states with the original constructor over three revolutions. There is no
numerical change. Analytical second volume derivatives were not previously
implemented by this class; they remain unavailable rather than introducing
finite differences or unrelated new mechanics. Toggle/branch diagnostics remain
accessible through the original assembly and slider-state APIs.

Harmonic examples and the historical `ideal_piecewise_linear` family are
retained. The latter is a reference parameterization, not a proven optimum.

## Free periodic motion

`FreeKinematics` exists to reveal thermodynamically desirable motion before
trying to reproduce it with a realizable mechanism. It has no prescribed
plateaus, transfer windows, Lambdas, valve timing, mirror symmetry or preferred
number of extrema. Small and large definitions are independent and may even
use different control counts, sharing only the same angle coordinate.

`FreeMotionDefinition` is an immutable, serializable dataclass with control
values, minimum/maximum volumes and optional bounds on absolute first and
second angular derivatives. `FreeKinematicsConfiguration` contains two such
definitions. TOML type `free` reads separate `[kinematics.small]` and
`[kinematics.large]` tables. Cylinder ranges come from the existing geometry
table, a single input authority; Python configurations reject mismatched ranges.
Existing sizing and similarity transformations update those ranges consistently.

There are N equally spaced distinct controls, with N >= 4; the endpoint value
is repeated internally, not an independent variable. The example uses six
controls per cylinder. A periodic cubic spline is C2, including at the seam;
its third derivative may jump and jerk constraints are deferred. See the
[SciPy CubicSpline contract](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline.html).

Controls are canonicalized to zero mean and unit Euclidean norm, removing
positive affine offset/amplitude equivalence. Already canonical serialized
controls are preserved on reconstruction. For future optimization,
`FreeMotionDefinition.from_shape_coordinates` maps N-2 unconstrained coordinates
through a stereographic sphere chart and an orthonormal zero-mean basis to N
controls. Thus six controls need **four shape coordinates per cylinder**, eight
for independent S/L. The chart excludes one pole; alternate charts or direct
control definitions cover that representation edge. Campaign bounds should be
placed on these chart coordinates, not on N independently mutable canonical
controls. No redundant global phase parameter is added.

True spline extrema are found from polynomial derivative roots plus knots,
not from a sampling grid. The entire law is then transformed affinely to the
requested physical minimum and maximum. This is the definition of the volume
parameterization, not clipping an invalid path. Between-knot overshoot is
included in normalization. Switching extrema can make shape-parameter
sensitivities nonsmooth; the path itself stays C2. No differentiability of a
future objective with respect to design coordinates is promised.

Analytic first and second derivatives come from the spline polynomial. Only
immutable coefficient and knot tuples are stored, evaluated by Horner's rule;
there is no mutable interpolation cache. The global first-derivative bound is
computed using second-derivative roots and knots. The second-derivative bound
uses knots, since it is piecewise linear. Optional signed margins are
`configured limit - measured maximum`. Volumes remain positive by their
validated bounds. Constant/unresolved or nonfinite controls and invalid ranges
are explicit construction errors. Derivative violations retain the unchanged
trajectory and diagnostics; `require_feasible()` raises
`KinematicConstraintViolation`. The factory refuses these candidates before
integration. The sizing evaluator reports `invalid_kinematics`, preserving
individual diagnostics and an explicit 'integration not started' message.
No objective penalty hides this state.

## Exchanger boundary

`ExchangerModel.build()` produces immutable `ExchangerComponents`:

- gas hold-up volume;
- inlet and outlet hydraulic closures (outlet includes the selected valve loss);
- either a static `HeatTransferModel` or optional `LumpedWallThermalModel`;
- geometric/other metadata and explicit validity-domain descriptions.

`connect_exchangers` maps these components onto the existing graph and keeps
passive valve thresholds. `MicrotubeExchanger` derives them using the unchanged
`MicrotubeBank`, `HardwareInputs`, `build_exchanger` and `TubeHalfLink`. It
retains the same tube/header loss allocation, wall capacity, thermal resistance
and external-air assumptions. `connect_hardware` remains a compatibility
constructor with its former return shape and report keys.

`HeatTransferModel` needs only `heat_rate(Tgas)`. Optional reservoir-normalized
isothermality diagnostics are unavailable if a static closure lacks a reference
temperature. Existing reservoir closures retain their former numerical results.
The one-wall integrator consumes already evaluated heat/storage rates instead
of rebuilding a hard-coded linear film. The equations and numerical values for
current air-wall models are unchanged; a test-only nonlinear model proves that
conductance and air-side input fields are not required by this integration.
The legacy `air_heat_w` key denotes external-source heat in this capability.

A future family may supply other geometry, correlations or measured behaviour
without sharing microtube input parameters. No new production exchanger family
is implemented. Current integrators support two static exchangers or two
one-wall-energy exchangers; mixed storage or distributed multi-state exchanger
models explicitly need another state-layout integrator. This limited capability
is stated rather than forcing every future family into a microtube/one-wall
parameter set. Air-specific trial plotting and fan reporting remain explicitly
specialized utilities, not generic evaluator requirements.

## Derived quantities and sizing migration

In a microtube candidate, geometry/material/correlation choices determine gas
volume, conductances, hydraulic closures and wall capacity. They must not also
be independent search variables. With `MachineDesign` exchanger designs,
legacy UA, hold-up and passage closures in the thermodynamic configuration are
superseded seeds, not additional campaign coordinates. Valve CdA remains an
explicit component input of the current microtube outlet closure.

Likewise, free physical volume ranges are separate from canonical shape; the
four-bar's stroke and trajectories are derived from its mechanism and configured
cylinder ranges. A free law has no inferred physical stroke, bore, realizable
linkage or mechanical loss model.

`SizingProblem`, `DesignPoint`, objectives and explicit constraint margins are
preserved. The legacy `DesignParameter` enum is not extended with spline or
microtube internals. For the next campaign, use a family-specific parameter
adapter from bounded coordinates to immutable `MachineDesign` inputs. Reuse
legacy `DesignPoint` for its supported operating/volume parameters, then apply
shape and exchanger-design coordinates in their own adapters. Reject attempts
to vary superseded derived exchanger quantities in that campaign adapter.
Cache/persist the complete serialized candidate, fixed family choices, model
assumptions and numerical settings, not only the legacy enum vector.

This refactor does not implement a global optimizer, persistence controller,
human review workflow, force/stress/friction/inertia model or speculative
mechanisms. Those are not prerequisites to interchangeable evaluation.

## Verification and examples

Before refactoring, all 206 existing tests passed. Focused tests now cover
periodicity, seam derivative continuity, true extrema, derivative bounds,
independent laws, scalar/vector agreement, deterministic immutable definitions,
serialization, sizing rejection, original four-bar trajectory equivalence and
both thermodynamic backends. Exchanger tests cover the compatibility connector,
original RHS equivalence, non-microtube static and nonlinear wall fixtures,
conservation, and explicit unsupported mixed state layouts. An AST check guards
against concrete kinematics imports in the generic integration modules.

Run the suite from the checkout:

```sh
PYTHONPATH=src python3 -m pytest -q
PYTHONPATH=src python3 examples/pluggable_kinematics_smoke.py
```

The smoke example saves `outputs/pluggable_kinematics_smoke.json`; it uses the
legacy constant-UA/orifice closures to validate periodic evaluation, not the
microtube hardware trial. Its powers must not be compared as hardware gains.
The current microtube case is separately checked with an unaccelerated cycle
from its saved converged state in
`outputs/motor_parallel_architecture_regression_screening.json`.

### Recorded verification results

The full periodic smoke cases converge in 14 cycles for four-bar and
37 cycles for free motion. The arbitrary free example consumes mechanical
work; successful integration is not an efficiency or optimality claim. The
standalone constant-UA smoke cases use different exchanger closures from the
current hardware study.

The saved hardware regression produces 37.626110567 W, compared with
37.626110390 W in the previous tighter-tolerance verification cycle
(difference 1.765e-07 W). Its scaled periodic error is
0.119902, below one. This check supports backward
compatibility, not independent exchanger calibration.

Final full-suite verification: **232 tests passed** in 39.09 seconds, compared
with 206 tests before the refactor (26 new focused tests). Both original
four-bar public tests and historical sizing/exchanger tests remain included.

## Persistent campaign layer

The first orchestration layer is now implemented. See
[OPTIMIZATION_CAMPAIGN.md](OPTIMIZATION_CAMPAIGN.md) for bounded family adapters,
Sobol continuation, exact candidate caching, durable history, approximate time
budgets, human reports and the physical free-motion smoke example. The next
review will choose local refinement and compatible periodic-state warm starts.
Dynamic-wall campaigns still require their own evaluator adapter.
