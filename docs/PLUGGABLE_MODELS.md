# Interchangeable kinematics and exchanger models

## Scope and construction

Four-bar motion is not a fundamental assumption of the cycle modeled by Dada Engine. Microtubes
are a supported geometry family, not a fundamental thermodynamic assumption.
Model composition preserves the conservative four gas volumes, hydraulic graph,
passive valves, heat/work signs and equations. No mechanical efficiency is introduced.

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

Core composition files are `free_kinematics.py`, `machine.py`, `exchangers/base.py` and
`exchangers/microtube.py`.

## Common kinematics interface

`kinematics.KinematicsModel` is a structural `typing.Protocol` exposing:

- small/large cylinder volume and first angular derivative;
- small/large `CylinderVolumeLimits`, including swept volume;
- optional physical strokes, represented by `None` for mathematical volume laws;
- derivative-discontinuity breakpoints for integration.

A combined
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

## Shared-crank four-bar kinematics

`SharedCrankFourBarVolumeKinematics` provides volume and first-derivative
trajectories from the assembly and slider-state APIs. Analytical second volume
derivatives are unavailable; finite differences are not substituted.
Toggle and branch diagnostics remain accessible through those APIs.

Harmonic motion and `ideal_piecewise_linear` are supported families. The latter
is a reference parameterization, not a proven optimum.

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
is repeated internally, not an independent variable. Each cylinder has its own
control vector. A periodic cubic spline is C2, including at the seam;
its third derivative may jump and jerk constraints are deferred. See the
[SciPy CubicSpline contract](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline.html).

Controls are canonicalized to zero mean and unit Euclidean norm, removing
positive affine offset/amplitude equivalence. Already canonical serialized
controls are preserved on reconstruction. For optimization,
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
- inlet and outlet hydraulic closures, with the declared valve placement/model;
- either a static `HeatTransferModel` or optional `LumpedWallThermalModel`;
- geometric/other metadata and explicit validity-domain descriptions.

`connect_exchangers` maps these components onto the existing graph and keeps
passive valve thresholds. `MicrotubeExchanger` derives them using the unchanged
`MicrotubeBank`, `HardwareInputs`, `build_exchanger` and `TubeHalfLink`. It
retains the same tube/header loss allocation, wall capacity, thermal resistance
and external-air assumptions. Assemble these components through
`connect_exchangers` or `MachineDesign`.

`HeatTransferModel` needs only `heat_rate(Tgas)`; no reservoir reference is
needed for validity. Existing reservoir closures retain their numerical results.
The one-wall integrator consumes already evaluated heat/storage rates instead
of rebuilding a hard-coded linear film. The equations and numerical values for
current air-wall models are unchanged; a test-only nonlinear model proves that
conductance and air-side input fields are not required by this integration.
The `external_heat_w` key denotes external-source heat in this capability.

A future family may supply other geometry, correlations or measured behaviour
without sharing microtube input parameters. The neutral external-stream family reuses the same wall-state capability with
declared external conductance; it does not invent liquid correlations. Integrators support two static exchangers or two
one-wall-energy exchangers; mixed storage or distributed multi-state exchanger
models explicitly need another state-layout integrator. This limited capability
is stated rather than forcing every future family into a microtube/one-wall
parameter set. Air-specific fan estimates are specialized model outputs, not
generic evaluator requirements.

## Derived quantities and sizing

In a microtube candidate, geometry/material/correlation choices determine gas
volume, internal conductance, hydraulic closures and wall capacity. Those derived
quantities must not also be independent search variables. In the external-stream
family, declared external conductance remains a separate scenario input. With `MachineDesign` exchanger designs,
configuration UA, hold-up and passage closures in the thermodynamic configuration are
superseded seeds, not additional campaign coordinates. Valve CdA remains an
explicit component input for rectangular microtubes. Circular collectors
derive a lossless-diode section from conduit geometry; see the microtube reference.

Likewise, free physical volume ranges are separate from canonical shape; the
four-bar's stroke and trajectories are derived from its mechanism and configured
cylinder ranges. A free law has no inferred physical stroke, bore, realizable
linkage or mechanical loss model.

`SizingProblem`, `DesignPoint`, objectives and explicit constraint margins are
preserved. The `DesignParameter` enum is not extended with spline or
microtube internals. Campaigns use a family-specific parameter
adapter from bounded coordinates to immutable `MachineDesign` inputs. Reuse
`DesignPoint` for its supported operating/volume parameters, then apply
shape and exchanger-design coordinates in their own adapters. Reject attempts
to vary superseded derived exchanger quantities in that campaign adapter.
Cache/persist the complete serialized candidate, fixed family choices, model
assumptions and numerical settings, not only the sizing enum vector.

## Campaign and verification boundary

Research composes these production models with fixed/active continuous, integer
and choice parameters, exact evaluation, local/global Sobol scheduling and
compatible warm starts. See [Research families](DADA_ENGINE_RESEARCH_KINEMATICS.md)
and [campaign internals](OPTIMIZATION_CAMPAIGN.md). Family-specific geometry
belongs in adapters; conservative integration remains independent of those choices.

Tests cover spline extrema/derivatives, branch/closure validation, direct
kinematics injection, exchanger ownership, conservation and original-model parity.

## Independent six-bar integration

`dada_solver.six_bar` provides `SixBarCylinderMechanism` and
`IndependentSixBarVolumeKinematics`. They compose the primary A-B-C-D linkage,
E rigid on BC, the E-F-G RR dyad, H rigid on EF, and a finite H-P rod with the
positive slider-closure branch. Circle closures and their differentiated
constraint equations provide analytic slider velocities. The two cylinders
have independent geometry and stored phases in the common study-angle domain.
There is no extra L reflection or phase adjustment. Injection through
`MachineDesign(..., kinematics=motion)` or `build_model(..., kinematics=motion)`
uses the single motor-direction transformation.

Assembly branches, full primary rotation and closure/singularity conditions
are checked. A 1440-interval scan brackets slider-velocity roots; Brent refinement
determines continuous travel extrema. Geometry is checked on the scan and at
every evaluation. This is numerical preflight, not a proof excluding all
possible singularities between scan points.

Normalization uses `q = 1 - (slider - slider_min) / stroke`, with volume obtained
from the given `CylinderVolumeLimits`. Refined extrema define the continuous
volume law without clipping. Lengths are expressed in crank-radius units;
`IndependentSixBarVolumeKinematics` reports physical strokes as `None`.
Research can derive physical stroke from explicit `crank_radius_m` metadata.

Research mechanism artifacts carry explicit geometry and branch selections; see
[DADA_ENGINE_RESEARCH_KINEMATICS.md](DADA_ENGINE_RESEARCH_KINEMATICS.md).

Research constructs this production family from declared coordinates or
versioned mechanism artifacts without changing thermodynamic equations. See
the [Research kinematics reference](DADA_ENGINE_RESEARCH_KINEMATICS.md) for
coordinate conventions, physical scale and artifact identity.
