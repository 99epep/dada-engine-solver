# Research kinematics and mechanism reference

## 1. Purpose and scope

This reference describes the current kinematic/mechanism layer used by schema-2
and schema-3 Research studies. Each study selects cylinder laws, declares their
fixed and active coordinates, and builds a production `KinematicsModel` for the
existing solver and campaign engine.

The authoritative registry and ownership rules live in
`dada_solver.research.families`: `FAMILIES`, `PHYSICAL_FAMILIES` and
`parameter_specs(settings, side)`. Search protocols do not define new families.

Thermal and hardware coordinates share the same fixed/active declaration and
candidate vector. Their physical meaning belongs in the
[configuration reference](DADA_ENGINE_RESEARCH_REFERENCE.md),
[microtube model](MICROTUBE_GAS_MODEL.md),
[external-stream model](EXTERNAL_STREAM_THERMAL_MODEL.md) and
[current ownership rule](PHYSICS_DECISIONS.md).

## 2. Quick start

From a source checkout:

```sh
research() {
    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m dada_solver.research "$@"
}

research init kinematics \
    --small slider_crank \
    --large harmonic \
    --output outputs/my_motion/study.toml

research validate outputs/my_motion/study.toml

research evaluate \
    outputs/my_motion/study.toml \
    --output outputs/my_motion/reference.json \
    --budget 3m
```

These paths are user-created destinations, not dependencies on committed
historical results. The generated template starts with fixed coordinates and
portable machine inputs; physical families also receive mechanism artifacts.
It is a starting configuration, not an optimized or recommended machine.

Use `evaluate` while all coordinates remain fixed. Replace selected `value`
fields with bounded active declarations before `run`. `validate` constructs and
screens geometry but does not integrate thermodynamics. Command and campaign
usage belongs in the [Research guide](DADA_ENGINE_RESEARCH.md).

## 3. Fixed and active coordinate ownership

A fixed coordinate is declared as:

```toml
[[parameters]]
name = "kinematics.small.phase_rad"
unit = "rad"
value = 4.345
```

To release it, replace that table with:

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

Active continuous declarations require `kind`, `initial`, `lower`, `upper`,
`unit` and `transform` (`linear` or `log`), in addition to `name`. Integer
coordinates use `kind = "integer"` and `encoding = "nearest_even_v1"` instead
of a continuous transform.

Kinematic assembly branches, direction signs, orientation booleans,
representation sizes and family selection remain fixed scientific categories.
This restriction does not prohibit explicitly supported categorical parameters
elsewhere in a machine study.

All family coordinates must be supplied by declarations or a mechanism artifact.
An active declaration removes its coordinate from the fixed set. Unknown,
duplicate, wrong-family and wrongly typed declarations are rejected, including
incorrect units. There is one active candidate vector in declaration order.
Search bounds describe the explored region, not an intrinsic physical domain.

## 4. Family registry

`FAMILIES` currently contains exactly these eleven representations:

| Family | Representation | Physical mechanism |
|---|---|---|
| `harmonic` | Cosine volume law | No |
| `slider_crank` | Centered or offset finite-rod slider crank | Yes |
| `four_bar` | Four-bar loop and finite output rod | Yes |
| `six_bar` | Primary loop, secondary dyad and finite piston rod | Yes |
| `free_spline` | Periodic cubic spline | No |
| `fourier_c2` | Normalized Fourier motion | No |
| `structured_c2_15p` | Structured C2 motion | No |
| `ideal_piecewise` | Ideal piecewise-linear chronology | No |
| `four_stage` | Four-stage piecewise-linear motion | No |
| `independent_four_stage` | Independently timed piecewise-linear motion | No |
| `hybrid_compact` | Compact hybrid motion with rounded and kinked branches | No |

`PHYSICAL_FAMILIES` is exactly `slider_crank`, `four_bar`, `six_bar`.
An abstract volume law does not establish linkage realizability.

## 5. Family-specific conventions

Coordinates below are local names, prefixed with `kinematics.small.` or
`kinematics.large.` in parameter declarations unless shared-crank ownership is
explicitly selected. Settings belong in `[kinematics.small]` or
`[kinematics.large]` and remain fixed.

### 5.1 `harmonic`

The production cosine volume law owns `phase_rad` in radians. Research uses the
production small-side phase convention independently for either cylinder.

### 5.2 `slider_crank`

The centered/offset finite crank slider uses the production positive-root closure.

| Coordinate | Unit / type |
|---|---|
| `rod_over_crank` | `crank_radius`, positive |
| `offset_over_crank` | `crank_radius` |
| `phase_rad` | `rad` |
| `volume_increases_with_coordinate` | Boolean, fixed; declaration unit `1` |

The orientation boolean determines the volume direction; it is not an additional
operation-direction transform. Closure is checked by the production mechanism.

### 5.3 `four_bar`

The production loop drives a finite output rod. The fixed setting
`output = "rocker"` or `output = "coupler"` selects the output-point frame.

| Continuous coordinates | Unit |
|---|---|
| `coupler`, `rocker`, `ground_x`, `ground_y` | `crank_radius` |
| `output_along`, `output_normal`, `rod_length` | `crank_radius` |
| `slider_origin_x`, `slider_origin_y` | `crank_radius` |
| `axis_angle`, `phase_rad` | `rad` |

`coupler`, `rocker` and `rod_length` must be positive. `loop_branch`,
`slider_branch` and `crank_direction` are fixed integer signs (`-1` or `+1`).
`volume_increases_with_coordinate` is a fixed boolean; all four categories use
unit `1` in declarations.

`envelope_frame_angle_rad` is optional screening-frame metadata. It defines the
frame used for the envelope proxy, not a reconstructed shaft layout.

### 5.4 `six_bar`

The topology is:

```text
A-B-C-D primary loop
E on BC
E-F-G secondary dyad
H on EF
finite H-P piston rod
P on slider axis
```

The production dataclass owns exactly 15 continuous coordinates, exposed as
`SIXBAR_CONTINUOUS`. `PRIMARY_COORDINATES` contains the first six and
`DOWNSTREAM_COORDINATES` the remaining nine.

| Group | Coordinates | Unit |
|---|---|---|
| Primary | `primary_ground`, `primary_coupler`, `primary_rocker` | `crank_radius` |
| Primary | `primary_e_along`, `primary_e_normal` | `crank_radius` |
| Primary | `primary_phase` | `rad` |
| Downstream | `second_pivot_x`, `second_pivot_y`, `link_ef`, `link_gf` | `crank_radius` |
| Downstream | `h_along_over_ef`, `h_normal_over_ef` | `1` |
| Downstream | `piston_rod`, `slider_axis_offset` | `crank_radius` |
| Downstream | `slider_axis_angle` | `rad` |

`primary_e_along` and `primary_e_normal` are lengths, not fractions of the
coupler. The two H coordinates are dimensionless fractions of EF.
`primary_branch` and `second_branch` are fixed integer signs with unit `1`.

See the [primary seed catalogue](PRIMARY_FOUR_BAR_FAMILIES.md),
[paired six-bar catalogue](SIX_BAR_MECHANISM_FAMILIES.md) and
[synthesis method](MECHANISM_SYNTHESIS_SEARCH.md) for geometry and methodology.

### 5.5 `free_spline`

A periodic cubic spline requires integer `count = N >= 4` and
`representation = "controls"` or `"shape_coordinates"`.

- `controls`: `control_0` through `control_(N-1)` are fixed continuous values.
- `shape_coordinates`: `shape_0` through `shape_(N-3)` provide `N-2`
  nonredundant coordinates, each fixed or active.
- Both representations own an independent `phase_rad` in radians.

Control and shape coordinates use unit `1`. The shape chart removes affine
control offset/amplitude, not phase. Spline extrema and derivative diagnostics
use the production spline implementation.

### 5.6 `fourier_c2`

The setting `harmonics = H >= 1` is a fixed integer. The family owns
`coefficient_0` through `coefficient_(2H-1)`, all with unit `1`.

Coefficients alternate cosine then sine for successive harmonics: indices
`2(k-1)` and `2(k-1)+1` multiply `cos(2πkt)` and `sin(2πkt)` respectively,
for `k = 1 ... H`. The production convention is
`t = (-theta / (2π)) mod 1`, followed by full-stroke normalization.
No separate phase coordinate is owned; phase is represented by coefficients.
The coefficient vector must describe nondegenerate motion.

### 5.7 `structured_c2_15p`

This structured C2 law owns 15 coordinates across the pair:

| SMALL (8) | LARGE (7) |
|---|---|
| `small_max_deg` | `large_down_duration_deg` |
| `small_down_duration_deg` | `large_max_curvature` |
| `small_max_curvature` | `large_min_curvature` |
| `small_min_curvature` | `large_down_bp_mid_q` |
| `small_up_bp_mid_q` | `large_up_kink_u` |
| `small_down_kink_u` | `large_up_kink_q` |
| `small_down_kink_q` | `large_up_kink_width_rel` |
| `small_down_kink_width_rel` | |

Coordinates ending in `*_deg` use degrees; the others use unit `1`.
Curvatures are normalized-position second derivatives with respect to normalized
motion time, not physical acceleration. Each side owns only its listed inputs.

### 5.8 `hybrid_compact`

This compact hybrid law owns nine coordinates across the pair:

| SMALL (5) | LARGE (4) |
|---|---|
| `small_max_deg` | `large_down_duration_deg` |
| `small_down_duration_deg` | `large_down_rounding` |
| `small_up_rounding` | `large_up_kink_u` |
| `small_down_kink_u` | `large_up_kink_q` |
| `small_down_kink_q` | |

The `*_deg` coordinates use degrees; rounding and kink coordinates use unit `1`.
Its name does not imply that every branch has an available second derivative:
use the metric availability contract below.

### 5.9 `ideal_piecewise`

This abstract chronology law owns `small_lambda_target`, `large_lambda_target`
and `adiabatic_sector_fraction`, all with unit `1`, on each selected side.
It does not describe a physical mechanism. Piecewise-linear derivative jumps
remain real jumps.

### 5.10 `four_stage`

Each side owns shared-origin stage timings and its own level coordinates:

- SMALL: `t1`, `t2`, `t3`, `a_s`, `b_s`;
- LARGE: `t1`, `t2`, `t3`, `a_l`, `b_l`.

All use unit `1`; timings are normalized cycle fractions. In independent
coupling these are separate declarations on each side, not hidden synchronization
between cylinders. Piecewise-linear derivative jumps remain real jumps.

### 5.11 `independent_four_stage`

The independently timed piecewise-linear law owns:

- SMALL: `t0_s`, `t1_s`, `t2_s`, `t3_s`, `a_s`, `b_s`;
- LARGE: `t1_l`, `t2_l`, `t3_l`, `a_l`, `b_l`.

All use unit `1`, with timings in normalized cycle fractions. Piecewise-linear
derivative jumps remain real jumps; no acceleration smoothing is implied.

## 6. Independent and shared-crank coupling

```toml
[kinematics]
coupling = "independent"
```

This is the normal generic mode and permits any valid SMALL/LARGE family
combination. It does not impose mirror symmetry.

`coupling = "shared_crank"` is the explicit common-shaft four-bar mode. Both
assemblies must already use a common physical frame and the same physical crank
radius if supplied. Research does not rotate arbitrary mechanisms into a common
shaft layout.

The per-side phase and direction coordinates are replaced by:

- `kinematics.shared.phase_rad`: fixed or active continuous, unit `rad`;
- `kinematics.shared.crank_direction`: fixed integer sign, unit `1`.

A shared direction remains categorical even when the phase is active.

## 7. Physical scale

For `slider_crank`, `four_bar` and `six_bar`, positive `crank_radius_m` is optional
family metadata. Without it, normalized geometry remains valid, but physical
stroke is unavailable. With it:

```text
physical_stroke = stroke_over_crank * crank_radius_m
```

This does not determine bore. Abstract families do not acquire a physical crank
scale merely by specifying a swept volume.

For finite-rod slider closures, changing the assembly branch changes the
kinematic law and potentially the stroke; it is not a rigid translation.
Crank phase is not generally piston-extremum phase, especially for offset sliders
and non-harmonic mechanisms. A common crank does not require symmetric outputs.

A sampled planar envelope is only a geometric screening proxy. It does not
certify component thickness, bearings, supports, cylinder bodies, axial layering
or collision clearance. Crossing line segments in a planar drawing alone do not
establish physical collision. Normalized geometry and acceleration/shape scores
do not determine loads, stress, fatigue, bearing forces or mechanical losses;
these require physical scale and separate engineering models.

## 8. Mechanism artifacts and libraries

`dada_solver.research.artifacts` provides `MechanismArtifact` and
`MechanismLibrary`. Artifacts are restricted to the three physical families.

- `MechanismArtifact.create(...)` validates complete family geometry, settings
  and mechanical-constraint declarations, and constructs production closure.
- `MechanismArtifact.from_data(...)` validates serialized schema and content hash.
- `MechanismArtifact.load(...)` loads an artifact, optionally checking an expected hash.
- `MechanismArtifact.reconstruct()` reconstructs the production geometry object.

An artifact contains its schema, settings, complete geometry, length and angle
conventions, mechanical constraints, content hash and provenance. Its hash
covers canonical scientific content; provenance is deliberately outside the
scientific mechanism hash. Changing geometry or settings changes identity.
Presentation/drawing layout is not scientific mechanism identity.

A minimal reference is:

```toml
[kinematics.small]
family = "six_bar"
artifact = "small.mechanism.json"
sha256 = "<canonical mechanism content hash>"
```

A declaration can release or override a geometry coordinate while the source
artifact remains part of study identity. Artifact settings cannot be relabelled:
changed scale or output settings require a corresponding artifact. Embedded
constraints remain applicable; configurable study design limits belong in
`mechanical_constraints`.

`MechanismLibrary` retains several named mechanism families and imposes no
winner. The [primary catalogue](PRIMARY_FOUR_BAR_FAMILIES.md) and
[six-bar catalogue](SIX_BAR_MECHANISM_FAMILIES.md) illustrate distinct lineages,
without defining an optimizer ranking.

## 9. Mechanical diagnostics and constraints

`side_metrics()` computes diagnostics from the production kinematic object.
`available_metrics(family)` declares what can legally be constrained for that
family; its exact availability is authoritative.

| Scope | Current metrics |
|---|---|
| All families | `zero_crossing_count`, `maximum_absolute_first_derivative` |
| Families with second derivative support | `maximum_absolute_second_derivative` |
| Physical families | `stroke_over_crank`, `minimum_rod_axis_cosine` |
| Slider crank | `closure_margin` |
| Four-/six-bar | `minimum_primary_transmission_sine` |
| Four-bar | `stroke_over_envelope` |
| Six-bar | `minimum_secondary_transmission_sine`, `EH_over_crank`, `H_axis_lateral_rms_over_stroke`, `H_axis_lateral_span_over_stroke`, `crank_axis_to_EFH_clearance_over_crank` |

Currently the second-derivative constraint is available for `harmonic`,
`slider_crank`, `free_spline`, `fourier_c2` and `structured_c2_15p`. Do not assume
it exists for other families. First and second volume derivatives use
`m^3/rad` and `m^3/rad^2`, not time derivatives.

```toml
[[mechanical_constraints]]
side = "small"
metric = "minimum_primary_transmission_sine"
relation = "minimum"
limit = 0.30
unit = "1"
```

`validate_mechanical_constraint()` in `dada_solver.research.margins` checks metric
availability, unit, relation and a finite nonnegative limit, including an integer
reversal count where applicable. These declared design limits are distinct from
intrinsic requirements such as real closure and valid branches.

The normal screening grid is 1440 samples; configurable grids require at least
360. Mechanical preflight precedes thermodynamic integration. A sampled screen
is not a continuous proof. Some diagnostics use stronger production methods,
such as analytic spline extrema; consult the recorded method rather than
assuming all metrics share the same guarantee.

## 10. Constraint margins

`margin_record()` exposes value, limit, relation, unit, availability, satisfaction,
absolute and relative margin, and presentation flags.

```text
minimum: margin = value - limit
maximum: margin = limit - value
equal:   margin = -abs(value - limit)
```

A nonnegative margin means satisfied when the value is available. Unavailable
values are reported as unavailable, not as evidence of satisfaction. Relative
margin is `margin / abs(limit)` for a nonzero numeric limit and a non-equality
relation; otherwise it is unavailable.

`near_active` flags a satisfied non-equality constraint within 5% of a nonzero
limit. It is a reading aid, not a safety factor and not a manufacturing tolerance
certificate. Margins retain their individual units and meanings; no universal
mechanical robustness score is implied.

## 11. Synthesis handoff

`dada_solver.research.synthesis` exposes the ordered `STAGES`, coordinate-group
helper `release_coordinates()`, and declarative contracts `SynthesisRequest`
and `FreshIslandPolicy`.

| Stage | Six-bar continuous coordinates released |
|---|---:|
| `primary_discovery` | 6 per selected side |
| `downstream_fit` | 9 per selected side |
| `full_local_polish` | 15 per selected side |
| `mirror_initialization` | None |
| `opposite_local_adaptation` | 15 per selected side |
| `paired_thermodynamic` | 30, requiring SMALL and LARGE |
| `hardware_retuning` | No mechanism coordinates |

Branches stay fixed. `release_coordinates()` is a grouping helper, not an
optimizer. The scientific handoff preserves:

- position as the primary synthesis reference;
- acceleration as diagnostic only;
- mirroring as initialization only;
- the actual study thermodynamic objective after paired release;
- fresh-island saturation separate from seeded design exploitation.

`FreshIslandPolicy` declares fixed bounds, distinct deterministic seeds, an
explicit search policy and both primary branches. `SynthesisRequest` validates
the ordered stages and separates the two protocols; it does not execute them.

Current Research can construct, screen, evaluate and Sobol-search declared
thermo-mechanical coordinates, persist/resume campaigns and compare actual
thermodynamic results.

It does not yet generically execute automatic primary-family discovery,
clustering, target-position cadence fitting, downstream-dyad fitting, local
polish, mirror/adaptation operators, synthesis-stage scheduling or fresh-island
saturation. Declared stages do not imply implemented synthesis operators.

See the [synthesis method](MECHANISM_SYNTHESIS_SEARCH.md) for the methodological
contract and why primary discovery, complete-mechanism fitting and thermodynamic
assessment are different tasks.

## 12. Kinematic visualization

`sample_motion()` in `dada_solver.research.visualization` compiles the study and
builds the same production kinematics used by evaluation, without integrating
thermodynamics. It samples study angle and returns:

- SMALL/LARGE physical volumes;
- normalized position and velocity per radian;
- normalized acceleration per radian squared when supported, otherwise unavailable;
- physical stroke when scale exists;
- joint positions and linkage connectivity for physical families.

Its angle domain is `study_angle_before_operation_transform`. Derivatives are
with respect to study angle, not time. Joint coordinates use stored local
mechanism frames in crank-radius units. No drawing offset or common physical
shaft layout is inferred.

`sample_report_volumes(...)` instead samples current production volumes in
solver-cycle angle, applying the operation transform exactly once. Its angle
output is in degrees and volume output in cubic metres; it also requires no ODE
integration. Report rendering and thermodynamic replay options belong in the
[cockpit reference](DADA_ENGINE_RESEARCH_COCKPIT.md).

## 13. Scope boundary and related references

This document owns kinematic family selection, coordinate ownership, physical
mechanism geometry, artifacts, mechanical diagnostics and constraint evidence,
synthesis handoff and kinematic visualization.

For other concerns, use:

- [Research guide](DADA_ENGINE_RESEARCH.md): commands and workflow;
- [configuration reference](DADA_ENGINE_RESEARCH_REFERENCE.md): machine inputs and campaign settings;
- [cockpit](DADA_ENGINE_RESEARCH_COCKPIT.md): report presentation and plots;
- [microtube model](MICROTUBE_GAS_MODEL.md): exchanger geometry and correlations;
- [external-stream model](EXTERNAL_STREAM_THERMAL_MODEL.md): external thermal boundaries;
- [working-fluid models](WORKING_FLUID_MODELS.md): thermodynamic property contracts;
- [validation and evidence](validation.md): current status and evidence boundaries;

Thermodynamic equations, hardware physics, campaign persistence internals and
historical parity reports are outside this kinematics reference.
