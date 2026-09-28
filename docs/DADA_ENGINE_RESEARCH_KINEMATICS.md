# Research V2 — kinematics and mechanism reference

Research can now select a family independently for each cylinder, declare fixed
or active coordinates, and evaluate them through the existing campaign engine.
Schema 2 adds construction and ownership; it does not replace the solver,
optimizer, history, exact cache, warm starts or Sobol continuation.
The [user guide](DADA_ENGINE_RESEARCH.md) introduces the commands.

## First configuration

From the checkout:

```sh
research() { PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m dada_solver.research "$@"; }
research init kinematics --small slider_crank --large harmonic --output outputs/my_motion/study.toml
research validate outputs/my_motion/study.toml
research evaluate outputs/my_motion/study.toml --output outputs/my_motion/reference.json --budget 3m
research report outputs/my_motion/reference.json --html outputs/my_motion/reference.html
```

Each template is entirely fixed initially. It is an executable reference, not an
optimized family or a scientifically recommended machine. The `.basis.json`
contains portable machine inputs; physical families also get a mechanism JSON
per cylinder. Neither source `examples/` nor old `outputs/` are runtime inputs.
A preset refuses to overwrite files. `validate` builds and screens geometry but
does not integrate thermodynamics.

To search a phase, replace its fixed table, rather than adding a duplicate:

```toml
[[parameters]]
name = "kinematics.small.phase_rad"
unit = "rad"
kind = "continuous"
initial = 4.345
lower = 4.30
upper = 4.40
transform = "linear"
```

The fixed equivalent is:

```toml
[[parameters]]
name = "kinematics.small.phase_rad"
unit = "rad"
value = 4.345
```

Then run and resume two bounded phases:

```sh
research run outputs/my_motion/study.toml --directory outputs/my_motion/run --budget 3m --max-candidates 1
research resume outputs/my_motion/run --budget 3m --max-candidates 1
research report outputs/my_motion/run --html outputs/my_motion/candidates.html
research compare outputs/my_motion/reference.json outputs/my_motion/run --html outputs/my_motion/comparison.html
```

`initial` controls standalone evaluation, not the first Sobol point. `--set`
overrides only declared active parameters, within their bounds. To change a fixed
coordinate, edit the study and start a new result/campaign. A fully fixed study
supports `evaluate`; `run` explains that at least one active coordinate is needed.

The same declaration works for geometry, thermal inputs, exchanger dimensions,
frequency, inventory and sizing. Active continuous coordinates require `kind`,
`initial`, `lower`, `upper`, `unit`, and `transform` (`linear` or `log`). Integer
counts require `kind = "integer"` and `encoding = "nearest_even_v1"` instead of
`transform`. Integer rounding preserves the existing inclusive-bound semantics.
Assembly branches, direction signs, orientation booleans, family selection and
representation sizes are **fixed scientific categories**, not Sobol floats.
Compare their alternatives in separate studies; categorical enumeration is deferred.

All family coordinates must be supplied by declarations or an artifact. Machine
coordinates omitted from TOML retain explicit basis defaults, included in the
scientific identity and resolved report. Active declarations remove a coordinate
from the fixed set. Unknown, duplicate, wrong-family and incorrectly typed inputs
are rejected. There is one candidate vector, in declaration order.

## Families and mathematical conventions

`[kinematics] coupling = "independent"` permits any valid small/large combination.
The family lives in `[kinematics.small]` / `[kinematics.large]`. Below, motion
parameters have prefix `kinematics.<side>.`. Angles passed to the production
model are **study angles**, before the factory applies its one motor-direction
transform. Historical motor-fraction laws retain `t = -theta/(2*pi) mod 1`.
This is part of their definition, not an extra Research direction correction.

| Family | Reused mathematical representation | Owned coordinates and units |
| --- | --- | --- |
| `harmonic` | Production cosine volume law | `phase_rad` [rad], independently on each side; large phase 0 reproduces the original reference origin |
| `slider_crank` | Finite centered/offset crank-slider, positive-root closure | `rod_over_crank`, `offset_over_crank` [crank_radius], `phase_rad` [rad], `volume_increases_with_coordinate` [boolean, unit 1] |
| `four_bar` | Production four-bar closure with finite output rod; `output = "rocker"` or `"coupler"` | See geometry table below |
| `six_bar` | Existing primary four-bar, downstream dyad, finite piston rod | The established 15 continuous coordinates plus two fixed branches |
| `free_spline` | Production periodic cubic spline, canonical normalization and stereographic chart | `count`, `representation` settings; fixed `control_0…N-1` or active/fixed `shape_0…N-3` [1]; `phase_rad` [rad] |
| `fourier_c2` | Smooth periodic Fourier law with refined extrema normalization | `harmonics = H`; `coefficient_0…2H-1` [1] per side, alternating cosine and sine coefficients |
| `structured_c2_15p` | Historical piecewise quintic C2 law, scalar PPoly fast path | Three timing/extremum parameters, four curvature magnitudes, two BP bowing coordinates, two three-coordinate HP kinks |
| `ideal_piecewise` | Existing ideal chronology law, including pauses | `small_lambda_target`, `large_lambda_target`, `adiabatic_sector_fraction` [1] per selected law |
| `four_stage` | Existing four-segment linear law with fixed global origin | `t1,t2,t3` and that side's `a_s,b_s` or `a_l,b_l` [1] |
| `independent_four_stage` | Existing independent four-segment law | Small: `t0_s,t1_s,t2_s,t3_s,a_s,b_s`; large: `t1_l,t2_l,t3_l,a_l,b_l` [1] |

The four-stage presets reproduce the historical shared seven-coordinate law
when both sides have the same stage times. V2 can release each cylinder's stage
times independently. No hidden equality ties them; if synchronized active stage
times are required, an explicit shared abstract-timing declaration is future
work. The independent representation additionally releases the small cyclic
origin. Unused other-side helper parameters never enter a candidate vector.
Derivative jumps in piecewise linear laws remain jumps; no smooth acceleration
is invented at their knots.

Fourier coefficients retain the historical ordering and normalization, including
its 2049-point bracketing followed by bounded scalar extremum refinement. This
is not an analytic global-extremum certificate for arbitrary unbounded harmonic
counts. Full-cycle screens remain explicit. Coefficient scale redundancy is
historical; no undocumented coordinate canonicalization changes stored laws.

Free-spline controls are canonicalized by `FreeMotionDefinition`. Fixed stored
controls plus an active phase reproduce the historical phase wrapper. Shape
search uses its non-redundant stereographic coordinates (`N-2`, not `N` controls);
the phase can independently be fixed or active. The chart removes affine
control offset/amplitude, not the continuous phase of the spline knot grid.
First/second derivative limits use the production spline extrema diagnostics,
not a noisy finite difference. Sharp extrema of a spline are not linkage targets
that must be copied mechanically.

### Structured C2 15p

The owned names are:

- Small: `small_max_deg`, `small_down_duration_deg`, `small_max_curvature`,
  `small_min_curvature`, `small_up_bp_mid_q`, `small_down_kink_u`,
  `small_down_kink_q`, `small_down_kink_width_rel`.
- Large: `large_down_duration_deg`, `large_max_curvature`,
  `large_min_curvature`, `large_down_bp_mid_q`, `large_up_kink_u`,
  `large_up_kink_q`, `large_up_kink_width_rel`.

Names are retained verbatim from the fitting lineage, e.g.
`kinematics.small.small_max_deg`. `_deg` uses degrees; the other coordinates use
unit `1`. Curvatures are magnitudes of normalized-position second derivatives
with respect to normalized motor time, **not radians or physical acceleration**.
The global large maximum fixes the time origin. Durations are strictly between
0 and 360 degrees; curvatures positive; kink positions and bowing values strictly
between 0 and 1; relative kink widths between 0 and 2. The template also explicitly
screens for two reversals and valid volume range. No shape clipping is performed.

```sh
research init structured-c2-3952 --output outputs/c2_reference/study.toml
research evaluate outputs/c2_reference/study.toml --output outputs/c2_reference/result.json --budget 3m
```

This preset uses **candidate 3952**, its original hardware and predecessor warm
state. It is not the separate candidate 3335 named best in the final search
report. The reference returns 40.719090181447584 W indicated and efficiency
0.23508063005202154 after three cycles on the recorded Numba backend.

### Four-bar geometry

Crank radius is the length datum, 1. Geometry coordinates are `coupler`, `rocker`,
`ground_x`, `ground_y`, `output_along`, `output_normal`, `rod_length`,
`slider_origin_x`, `slider_origin_y` [crank_radius], and `axis_angle`, `phase_rad`
[rad]. `loop_branch`, `slider_branch`, `crank_direction` are signed integers
-1/+1; `volume_increases_with_coordinate` is boolean. Output-point coordinates
retain the selected rocker/coupler frame semantics of `four_bar.py`.

`coupling = "shared_crank"` currently requires two four-bar assemblies already
expressed in a common frame. The phase and direction then have names
`kinematics.shared.phase_rad` and `kinematics.shared.crank_direction`, and a
single physical radius applies to both. Research never silently rotates an
arbitrary mechanism into a guessed shaft layout. The compact historical seed
is already reconstructed in the common frame. Its `envelope_frame_angle_rad`
preserves the original local bounding-box screen: an axis-aligned bounding-box
diagonal is not rotation invariant. This scientific screening frame is distinct
from a drawing offset.

### Six-bar geometry

Topology: A–B–C–D primary loop, E attached to BC, E–F–G secondary dyad, H
attached to EF, finite rod H–P, with P constrained to its slider axis.

| Stage | Exact coordinate names | Units |
| --- | --- | --- |
| Primary (6) | `primary_ground`, `primary_coupler`, `primary_rocker`, `primary_e_along`, `primary_e_normal`, `primary_phase` | First five crank_radius; phase rad |
| Downstream (9) | `second_pivot_x`, `second_pivot_y`, `link_ef`, `link_gf`, `h_along_over_ef`, `h_normal_over_ef`, `piston_rod`, `slider_axis_offset`, `slider_axis_angle` | Lengths crank_radius; H fractions 1; angle rad |
| Fixed branches | `primary_branch`, `second_branch` | Signed integer -1/+1, unit 1 |

`E_along` / `E_normal` are lengths in crank radii. H coordinates are fractions
of EF. Do not interchange those conventions. The volume law keeps
`q = 1 - (slider - minimum) / stroke` and the production refined extrema.

For all physical families, `crank_radius_m` is optional family metadata. Without
it, physical stroke stays unavailable; dimensionless mechanism geometry does
not determine cylinder stroke in metres. Supplying it scales the geometry's
stroke, without changing the normalized volume law or assuming a bore.

## Machine and hardware ownership

| Coordinates | Units / interpretation |
| --- | --- |
| `volume.swept_ratio`, `volume.total_swept_m3` | Small/large swept ratio [1], total swept volume [m^3] |
| `volume.small_clearance_ratio`, `volume.large_clearance_ratio` | Minimum enclosed / swept volume [1] |
| `operation.frequency_hz`, `charge.total_mass_kg` | Positive Hz, kg; operation sign comes from the basis |
| `microtube.heat_in.*`, `microtube.heat_out.*` | Integer `tube_count` [1]; `tube_length_m`, `inner_diameter_m`, `wall_thickness_m`, `pitch_m`, `header_depth_m` [m] |
| `thermal.heat_in.*`, `thermal.heat_out.*` | `air_inlet_temperature_k` [K], `air_mass_flow_kg_s` [kg/s], `metal_conductivity_w_m_k` [W/(m*K)], `metal_density_kg_m3` [kg/m^3], `metal_cp_j_kg_k` [J/(kg*K)] |

Heat transfer, hold-up, film resistance, wall capacity and losses are derived
from production hardware. UA is not an independent coordinate. Valve outlet CdA
uses the explicit `source_cda_times_count_ratio_v1` or `fixed_source_cda` policy.
Working-gas properties, numerical physics, valve placement and reservoir inputs
remain portable basis configuration, not an invented universal parameter list.
V2 currently requires explicit gas inventory; a pressure-charge study must first
resolve that inventory explicitly, as the Fourier parity fixture does.

For a refrigerator, the basis selects positive angular speed, the objective
`maximize_cooling_cop`, and constraints such as `minimum_cooling_power` and
`maximum_mechanical_input_power` [W]. The latter limits **indicated** input, not
human shaft power. Shaft losses remain unknown. V2 validated refrigeration with
the reservoir branch. V3 also validates the conservative wall/external-stream
branch in refrigeration; see [external thermal streams](EXTERNAL_STREAM_THERMAL_MODEL.md).
Both use the same kinematic and campaign interfaces. No calibrated human-power
model or cooling-cell optimization is
supplied. Generic kinematic constructors contain no motor-efficiency objective. The historical motor-time laws retain their angle
orientation; selecting refrigeration does not silently redesign their chronology.

## Artifacts, scientific identity and reconstruction

`MechanismArtifact.create(family, geometry, settings=..., constraints=...,
provenance=...)` produces a versioned scientific payload and canonical content
hash. `save`, `load`, `from_data` verify complete geometry, family, units and
branches. Artifacts record scale convention, phase/orientation, slider geometry,
constraints and provenance. Presentation layout is not scientific identity;
artifact provenance is excluded from the mechanism hash.

A minimal fixed input is:

```toml
[kinematics.small]
family = "six_bar"
artifact = "small.mechanism.json"
sha256 = "<canonical content hash from the artifact>"
```

This hash identifies canonical scientific content, not the pretty-printed file
bytes. The machine basis separately uses a byte SHA-256. A declaration can
release or override a geometry coordinate from an artifact; resolved fixed
values and the source artifact both remain in the study identity. Scale/output
settings cannot be silently relabelled: create a new artifact to change them.
`reconstruct()` returns the production geometry object. `families.build_side`
combines the artifact's geometry/settings with cylinder limits and a selected
side to build its exact `CylinderLaw` and `KinematicsModel` backend.

`MechanismLibrary` retains a set of named family entries and member mechanisms
plus metadata. It imposes no winner and never replaces primary-family diversity
with the primary search score. The documented ranks 1, 4, 12 and 50 are covered
by geometry/diagnostic round-trip tests; rank is a source identifier, not a new
taxonomy of mechanism families.

Study identity includes fixed/active ownership, families, artifact scientific
content, basis, numerical settings, constraints and policies. Candidate identity
also includes normalized and decoded active values through the existing engine.
Cache keys remain exact. Source paths, study display name and wall-clock phase
budget do not define physics. Changing scientific input starts a new campaign.
The campaign snapshots mechanisms beside its existing definition and history;
there is no second history format or implicit mutation of old records.

Schema 1 `init sixbar-thermo5d` still works. Existing histories and reports remain
readable offline. As before, source edits change the strict runtime fingerprint;
execution resume of a pre-V2 campaign requires its original matching source
version. There is no automatic conversion or identity rewrite.

## Feasibility and constraint evidence

Construction uses the original production closure/branch checks, then cheap
family-specific screens before integration. `[[mechanical_constraints]]`
contains `side`, `metric`, `relation` (`minimum`, `maximum`, `equal`), `limit`,
and `unit`. Artifact constraints are inherited; duplicate declarations are
rejected. A study can define both a lower and upper stroke bound.

| Family / metric | Meaning and evidence |
| --- | --- |
| All: `zero_crossing_count` [1] | Signed velocity reversals on the stated cyclic sample grid; no smoothing or unsigned-speed substitution |
| All: `maximum_absolute_first_derivative` [m^3/rad] | Volume derivative; spline uses analytic extrema, others a sampled bound |
| Smooth abstract / slider: `maximum_absolute_second_derivative` [m^3/rad^2] | Available analytic derivative sampled over cycle; spline uses exact polynomial extrema |
| Physical: `stroke_over_crank` [crank_radius], `minimum_rod_axis_cosine` [1] | Stroke / finite-rod alignment; slider-crank has analytic full-cycle bounds |
| Slider: `closure_margin` [crank_radius] | rod/crank − 1 − absolute offset/crank; must already be strictly positive for construction |
| Four/six: `minimum_primary_transmission_sine` [1] | Production primary closure transmission margin |
| Four: `stroke_over_envelope` [1] | Historical local envelope proxy, not a packaged machine size |
| Six: `minimum_secondary_transmission_sine`, `H_axis_lateral_rms_over_stroke`, `H_axis_lateral_span_over_stroke` [1] | Downstream transmission and transverse motion relative to piston stroke |
| Six: `EH_over_crank`, `crank_axis_to_EFH_clearance_over_crank` [crank_radius] | Output reach and crank-axis clearance to rigid EFH triangle |

Available metrics, including unconstrained ones, are retained in each result's
`derived.kinematic_metrics`. Constraint records contain current value, relation,
limit, unit, signed absolute margin, relative margin where meaningful, and
`satisfied` / `violated` / `unavailable`. Near-active satisfied constraints are
amber within 5% of a nonzero limit. This is a reading aid, not a safety factor.
A preflight rejection keeps its geometric evidence and unavailable thermal
constraints. Categorical verdicts have no relative SI margin. Legacy categorical
margins retain their original convention in read-only inspection.

The default mechanical grid has 1440 angles (minimum configurable 360). Sampled
screens are **not continuous proofs**; increasing samples cannot establish one.
Historical six-bar thresholds are explicit preset inputs, not universal defaults
for all machines. The maximum large-cylinder enclosed volume is separately
screened; motor studies retain the 0.066 m^3 ceiling. No constraints are weakened
to make a demonstration feasible. Reports offer no automatic scientific advice.

## Established synthesis workflow and current boundary

The authoritative methodology remains
[Mechanism synthesis search theory](MECHANISM_SYNTHESIS_SEARCH_THEORY.md),
[Primary four-bar families](PRIMARY_FOUR_BAR_FAMILIES.md),
[Six-bar mechanism families](SIX_BAR_MECHANISM_FAMILIES.md), and
[Compact kinematic synthesis](COMPACT_KINEMATIC_SYNTHESIS.md).
This interface carries those decisions; it does not replace those documents.

`research.synthesis.release_coordinates(stage, sides)` identifies the six-bar
coordinates to release. Fixed/active declarations and independent artifacts
represent primary-only 6-D, downstream-only 9-D, full local 15-D, opposite-side
15-D and paired 30-D studies. Branches remain fixed. Hardware retuning freezes
both artifacts and activates the relevant machine coordinates. An initialization
from the opposite mechanism becomes a new independent artifact; mirror symmetry
is not a persistent constraint. Family libraries preserve multiple basins.

`SynthesisRequest` records target scientific study identity, physical family,
ordered stages, mechanical constraints and retained family IDs. Position remains
the main reference, acceleration diagnostic-only, mirroring initialization-only,
and the paired objective the study's actual thermodynamic objective.
`FreshIslandPolicy` separately records fixed bounds, independent deterministic
seeds, a fixed named search policy and both primary branches. Saturation requests
reject retained seeds; design requests reject a saturation policy.

**Executable now:** construct, screen, evaluate, and Sobol-search declared
thermo-mechanical coordinates through `OptimizationCampaign`; retain artifacts
and compare real thermodynamic results. A researcher can create successive
studies with the appropriate coordinates fixed/released.

**Deferred:** automatic primary-family discovery/clustering, target-position
cadence fitting, downstream-dyad fitting, local numerical polish, mirror/adaptation
operators, automatic stage scheduling, and fresh-island saturation execution.
The request/policy objects are handoff contracts, not hidden optimizers. The
six-bar `release_coordinates` helper must not be interpreted as a four-bar
parameter mapping. No generic least-squares/acceleration fitter was introduced.
No physical reversal is removed by segmentation, and no PCA trajectory is
silently substituted for the final finite-rod slider motion.

## Visualization and validation

`research.visualization.sample_motion(study, active_values, samples=361)` exports
normalized position, study-angle velocity and available acceleration, physical
volume, optional stroke/scale, and physical joint coordinates/connectivity for
slider, four-bar and six-bar families. It constructs the same production objects
as evaluation and never integrates thermodynamics. Time derivatives require the
chosen angular speed; acceleration is not supplied where the backend lacks it.
Full animation and browser motion replay are deferred. Display layouts must
remain separate from these scientific frames.

See [validation](DADA_ENGINE_RESEARCH_VALIDATION.md) for measured tolerances and
bounded demonstrations, and [migration status](DADA_ENGINE_RESEARCH_MIGRATION_MATRIX.md)
for the distinction between mathematical families and historical protocols.
