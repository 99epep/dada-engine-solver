# Research kinematics and mechanism reference

## 1. Purpose and scope

This reference describes the current kinematic/mechanism layer used by schema-3
Research studies. Each study selects cylinder laws, declares their
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

These paths are user-created destinations. The generated template starts with
fixed coordinates and
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

## 11. Motion targets and synthesis handoff

`dada_solver.research.motion_target.MotionTarget` is an immutable, hashed
scientific target independent of the source family. It contains both SMALL and
LARGE on an increasing study-angle grid from zero through `2*pi`, including the
periodic endpoint. Study angle precedes the single operation transform applied
by the thermodynamic factory. Position is normalized using declared cylinder
limits, `(volume - minimum_volume) / swept_volume`, without sampled rescaling or
clipping. Available derivatives are per study-angle radian; absent derivatives
are `null`, never finite-differenced into existence.

Each side carries extensible named events with angles, optional directed/wrapped
intervals and source descriptors. Extraction retains a periodic seam and
bracketed derivative-root turnarounds when available. For `hybrid_compact`, exact
source parameters additionally identify extrema, rounding intervals, cadence
joins and kinks. Event metadata also retains exact normalized position and velocity at these
source-defined features. Kink metadata carries its branch boundaries; no neighborhood
width or refit heuristic is inferred. Events describe the source, not mandatory
features of a realizable mechanism.

`target_from_study()` extracts production kinematics at configured or explicitly
supplied active coordinates without integration. `load_motion_target()` accepts
a saved target, a schema-3 study TOML, or an identified stored Research result.
Campaign extraction requires a candidate selector. Scientific source identity,
coordinate conventions, events and samples determine the target hash;
presentation provenance is stored separately. Extraction reconstructs current
production kinematics, not a thermodynamic trajectory replay.

```console
dada-research motion-target path/to/study.toml --output path/to/target.json
dada-research motion-target path/to/campaign --candidate CANDIDATE_ID \
    --output path/to/target.json
dada-research mechanism synthesize path/to/target.json --family slider_crank \
    --stage global_discovery --validate-only
```

`PROTOCOLS` and `synthesis_protocol()` centrally own ordered stages and continuous
coordinate groups. Slider-crank and four-bar protocols use `global_discovery`,
`full_local_polish`, `paired_thermodynamic` and `hardware_retuning`; neither has
`downstream_fit`. Slider-crank releases `rod_over_crank`, `offset_over_crank` and
`phase_rad`, keeping orientation fixed. Four-bar releases its eleven continuous
coordinates, keeping output type, closure branches, direction and orientation
fixed. The six-bar protocol retains its hierarchy below. A paired release
requires both independent sides; hardware retuning releases no mechanism
coordinates.



`dada_solver.research.synthesis` exposes the six-bar ordered `STAGES`, coordinate-group
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
explicit search policy and both primary branches. This primary-loop saturation
contract applies to four-bar and six-bar, not slider-crank. `SynthesisRequest` validates
the ordered stages and separates the two protocols; it does not execute them.

Current Research can construct, screen, evaluate and Sobol-search declared
thermo-mechanical coordinates, persist/resume campaigns and compare actual
thermodynamic results.

`SynthesisPlan.execute()` runs direct `global_discovery` and `full_local_polish`
for slider-crank and four-bar, writing diverse `MechanismLibrary` members.
The shared `SearchPolicy` controls independent deterministic islands, bounds,
discrete categories, sampling, clustering and budgets. Each artifact retains
its target identity, search origin and separate evidence; polish additionally
records its parent family and artifact. Six-bar `primary_discovery`,
`downstream_fit`, `full_local_polish`, `mirror_initialization` and
`opposite_local_adaptation` execute independently through the same engine.
Fresh-primary saturation is not implemented by this geometric operator.
Paired thermodynamics and hardware retuning generate ordinary Research studies
through `mechanism adapt` and `mechanism retune`, using explicit complete Research
sources rather than a geometric proxy objective.

### Paired thermodynamic study generation

`paired_thermodynamic()` in `research.mechanism_adaptation` accepts a source
study at its declared initials, a stored evaluation, or a campaign with an
explicit candidate selector. A complete selected `MechanismLibrary` member
supplies SMALL and LARGE. A `MotionTarget` alone cannot supply the machine.
The source may use `hybrid_compact` directly; no spline intermediate is required.

```console
dada-research mechanism adapt path/to/source-campaign --candidate best \
  --library path/to/mechanical-pair.json --family-id FAMILY_ID \
  --output path/to/thermodynamic/study.toml --radius 0.1
dada-research validate path/to/thermodynamic/study.toml
dada-research run path/to/thermodynamic/study.toml
```

Generation writes a portable machine basis and two unchanged `MechanismArtifact`
files referenced by relative path and scientific hash. It screens geometry but
never starts integration. Only continuous mechanical coordinates are active:
**6 for slider-crank, 22 for four-bar, 30 for six-bar**. Each side retains its
fixed categories, branches, output and orientation. Mixed-family pairs are also
allowed, with the sum of their respective coordinate counts.

All other coordinates are fixed at the exact selected source values. Objectives,
thermodynamic constraints, numerical settings and execution budgets are retained.
Reference-pressure filling is explicitly converted to fixed selected inventory:
changing motion must not refill the machine. This policy change is recorded.
Source mechanical constraints and artifact constraints both apply; overlapping
bounds use their stricter conjunction, and incompatible equalities are rejected.
An unavailable family-specific source metric fails validation rather than being
silently dropped. Exactly two piston reversals are required on each side.
Production analytic-velocity root refinement screens topology before integration;
four-/six-bar root sampling is refined at increasing resolution, not claimed as
a mathematical continuous-domain proof.

The `centered_mechanical_local_regions_v1` search policy uses the existing
`refine` center-first Sobol scheduler (`local_regions_v1`). Each coordinate's
reference box is centered on the selected geometry with the width of the
versioned synthesis default bounds; positive-length half-widths are capped at
half the current value. Periodic angles remain unwrapped about the initial angle.
`--radius` defaults to 0.1, in (0, 0.5], and specifies the normalized half-width
inside that box. Thus positive lengths move by at most 10% by default when their
half-width cap applies, and phase windows have a default half-width of 36 degrees.
These are configurable **search bounds**, not physical domains. The exact pair is
evaluated first; hardware is never active.

The basis records source study/definition/candidate IDs, source assessment when
available, library and artifact hashes, selected family lineage, initial geometry,
effective reference boxes and warm-start decisions. Compatible selected-candidate
wall states use the standard `source_exact` initial-guess path, with the source's
safe-retry policy unchanged. Without one, the source initialization policy is
retained. No periodic solution is assumed for the new motion. This provenance
supports comparison with the source and adapted Research results without any
COP/fit composite objective.

### Hardware retuning with fixed mechanisms

`hardware_retuning()` in `research.hardware_retuning` generates a study after
paired thermodynamic adaptation. Its source is an exact paired candidate, an
evaluation, a study at its declared initials, or a retuning descendant. It never
starts integration and releases **zero mechanism coordinates** for all three
physical families. Geometry comes from the selected candidate's effective values,
not merely the initial artifact still linked by the paired study. Unchanged
geometry retains its artifact hash and provenance; changed geometry gets its
correct scientific hash and an explicit link to the parent artifact/candidate.
Branches, output, orientation, scale and mechanical constraints remain fixed.

```console
dada-research mechanism retune path/to/paired-campaign --candidate best \
  --scope source-active --radius 0.1 --output path/to/retuning/study.toml
dada-research mechanism retune path/to/paired-campaign --candidate CANDIDATE_ID \
  --group exchangers --group volumes --output path/to/retuning-hx/study.toml
dada-research validate path/to/retuning/study.toml
dada-research run path/to/retuning/study.toml
```

`source-active` (also the default without a targeted selection) restores only
originally active, still available **non-kinematic** parameters. `mechanism adapt`
records their complete Research declarations and a domain content hash in its
portable machine basis; no original file path is required afterward. If that
snapshot is absent, `--source-study` supplies the original source and its study ID
is verified. An empty immediate source domain can also use an explicitly supplied
source identified in `motion_refit` provenance; no change to motion refit or
invention of bounds is involved. Fixed original parameters cannot be reopened.

Targeted selection is the union of repeated `--group` and `--parameter NAME`.
It cannot be combined with explicit `--scope source-active`. Groups only filter
parameters actually active in the original source:

| Group | Existing parameter namespaces |
|---|---|
| `exchangers` | `microtube.*`, `thermal.*` |
| `volumes` | `volume.*` |
| `frequency` | `operation.frequency_hz` |
| `charge` | `charge.*` |
| `valves` | `valve.*` |
| `external-stream` | `external_stream.*` |

Unavailable source parameters are excluded and recorded; explicitly requesting
one is an error. An empty selection is rejected. All unselected coordinates stay
fixed at the current candidate, and every reopened initial value must lie in its
original source domain. No exchanger scaling or capacity transformation occurs;
this operation is independent of `rescale`.

The `source_active_hardware_local_regions_v1` policy retains original declarations,
units, continuous transforms, integer encodings and choice sets. The existing
center-first `local_regions_v1` Sobol scheduler clips each numeric interval to
`[max(0, center-radius), min(1, center+radius)]` in the **original normalized
domain**. `--radius` defaults to 0.1 and lies in (0, 1]. Effective physical bounds
are recorded, including nearest-even integer decoding; no domain extrapolation
occurs. Choice order is not a geometric metric: explicit
`choice_scope = "declared_choices"` keeps every original category eligible,
independently of numeric radius. The exact current candidate is evaluated first.

Objectives, thermodynamic and mechanical constraints, charge policy, numerical
settings and safe-retry policy come from the current scientific source. Explicit
inventory is not changed back into reference-pressure filling. Charge can be
reopened only if originally active and selected. Compatible current wall guesses
use standard `source_exact` initialization: original-inventory gas states are
preserved, other inventories use the same initial-guess rescaling as cached
Research states, and wall energies follow current capacities. Integrated mass is
never renormalized and convergence criteria remain unchanged.

Provenance retains the abstract/domain source, paired study and candidate IDs,
frozen mechanism hashes, current parent candidate, parameter selection, radius,
effective bounds, warm-start decisions and preceding retuning passes. A later
pass can reopen volumes plus exchangers after an exchanger-only pass using the
same original domains and the new current candidate. Ordinary `report`/`compare`
identities and source assessments support the before/after comparisons; no
additional score or reporting system is introduced.

The [direct-operator reference](MECHANISM_SYNTHESIS_SEARCH.md#direct-slider-crank-and-four-bar-operators)
defines default search bounds, the position-dominant score, progressive mechanical
screens, circular geometry distance and category-stratified retention. These
numerical search policies do not redefine any mechanical domain. A target may
come from any current kinematic family, including `free_spline`.

Discovery and polish write a standalone sortable HTML catalogue, defaulting to
the library path with an `.html` suffix. It exposes IDs, linkage views, point
paths, target/motion curves, fit errors, categories and individual mechanical
metrics. No thermodynamic performance is inferred. Selected members use the
same production joints in interactive `mechanism visualize` animations.

```console
dada-research mechanism synthesize path/to/target.json --family slider_crank \
    --stage global_discovery --side both --budget 2m --output path/to/slider-library.json
dada-research mechanism synthesize path/to/target.json --family four_bar \
    --stage global_discovery --side both --islands 8 --seed 1234 \
    --budget 2m --output path/to/fourbar-library.json
dada-research mechanism synthesize path/to/target.json --family four_bar \
    --stage full_local_polish --side large --library path/to/fourbar-library.json \
    --family-id large-family-003 --output path/to/polished-library.json
```

`--config` overrides the versioned search policy with explicit bounds, categories
and mechanical constraints; CLI search options override matching config fields.
`--max-evaluations` is reproducible under a fixed numerical runtime; cooperative
wall budgets stop at machine-dependent points. `--validate-only` does no search.
Dimensional derivative constraints require explicit cylinder volume limits.
Polish preserves all parent settings, categories and embedded constraints and
optimizes each selected family independently, without reducing them to one winner.

`SynthesisPlan.artifact()` and `member()` also accept submitted production geometry.
`assess_mechanism()` keeps three levels separate:

- angle-weighted RMS and maximum position error, and optional velocity RMS;
- production mechanical metrics and individual constraint margins;
- optional identified solver results supplied separately, never inferred from fit.

Acceleration is not included in the fit score. Without explicit machine volume
limits (`volume_limits` per piston in `member()`), dimensional volume-derivative
metrics are unavailable; dimensionless
geometry metrics remain available. A normalized artifact does not imply a
physical machine volume or crank scale. `mechanism_catalogue()` preserves all
members and their order, with family IDs, mechanism family, fit, mechanical
metrics, provenance, artifact hash and optional thermodynamic evidence.

### Hierarchical six-bar synthesis

The [executable six-bar reference](MECHANISM_SYNTHESIS_SEARCH.md#executable-hierarchical-six-bar-operators)
owns the versioned bounds, structured initialization, score terms and CLI chain.
The hierarchy is strictly six-dimensional primary discovery, nine-dimensional
fixed-primary downstream discovery, then fifteen-dimensional **local** polish.
No global fifteen-dimensional fit is provided.

Intermediate primaries use the existing `MechanismArtifact` format with
`family = "six_bar", component = "primary"`; they store only the six primary
coordinates and `primary_branch`. `SixBarPrimaryMechanism` supplies the same
production A–B–C–D/E closure used by complete six-bars. It has no piston law.
Primary evidence labels a principal-axis E projection as a chronology opportunity,
not a P fit; downstream constraints remain explicitly deferred until applicable.
Both primary branches and both downstream branches are explored categorically.

Each stage consumes explicit retained family IDs, produces a reusable
`MechanismLibrary` and shared HTML catalogue, and keeps parent family/artifact
hashes in provenance. Downstream cannot merge different primary parents.
Fit, mechanical constraints, intermediate preferences and thermodynamics remain
separate. Effective parent-relative bounds and seeds make geometric searches
inspectable and reproducible under the same numerical runtime.

`mirror_initialization` selects one destination side, reflects the source and
aligns the destination maximum to its target. It retains the first side artifact
unchanged. `opposite_local_adaptation` locally releases the destination's fifteen
coordinates, preserving its branches and leaving the first side untouched;
the pair has no implicit symmetry coupling. SMALL or LARGE may be synthesized
first. All intermediate and complete artifacts use the ordinary interactive
`mechanism visualize` path; no renderer or geometry is duplicated.

### Feature-aware motion refit

`MotionRefitRequest.execute()` geometrically fits a `MotionTarget` to
`free_spline` with **15 points per piston**: uniform nodes spaced by 24 degrees,
13 `shape_coordinates` and an independent `phase_rad` for each piston.
`MotionRefitResult` retains the target identity, canonical controls, coordinates,
phases, per-event errors and numerical policy (`uniform15_feature_position_phase_v1`).
It contains no thermodynamic result. Refit is an abstract-motion initialization,
not a mechanism fit.

The production phase convention is `q(theta - phase_rad)` in study angle before
operation transformation. Nodes therefore occur at `phase_rad + j*2*pi/15`.
The deterministic phase exploration spans one node spacing; whole-node shifts
are equivalent to circular control permutations. The default `RefitPolicy`
explores 32 phases, fits each phase with local position least squares, then
polishes controls and phase together around the best grid result. It uses no
stochastic global search and does not claim a globally optimal fit.

The frozen target is periodically interpolated using position and available first
derivatives, enriched by exact source event values. Position-only targets use a
periodic cubic interpolant without claiming an available target velocity.
Target acceleration is never used or reconstructed for fitting;
`acceleration_role = diagnostic_only` remains the contract.

A uniform background mesh of 1440 samples covers the complete cycle. Additional
sample blocks use the following explicit numerical weights and neighborhoods:

| Feature | Total position weight | Neighborhood |
| --- | ---: | --- |
| Background | 1 | Entire cycle |
| Maximum / minimum | 2 each | Quarter node spacing on each side; half the weight at the event |
| Derivative-root turnaround | 1 each | Same rule as extrema |
| Rounding | 0.5 each | Full directed interval plus exterior radius `min(spacing/8, interval/4)` |
| Cadence change | 0.25 each | Eighth node spacing on each side |
| Kink | 1 each | Radius `min(spacing/2, branch_length/8)` |
| Periodic seam / other descriptor | 0.1 each | Quarter node spacing on each side |

Feature blocks use 33 interior samples by default, include the exact event, and
rounding blocks include exterior samples. Weights are normalized after assembling
all blocks. These are inspectable numerical initialization policies, not physical
domains or mandatory features of a realizable mechanism. Exact hybrid event
metadata supplies the branch limits; no acceleration heuristic is involved.

Position is primary: local least squares minimizes weighted position MSE. Phase
selection first keeps only fits within `1e-12` MSE of the best position result,
then compares `position_MSE + 1e-4 * normalized_velocity_MSE` within this numerical
tie. Velocity differences are normalized by the node spacing before squaring.
The report exposes both terms separately. Velocity cannot trade a materially
worse position fit for a lower score. Missing target velocity contributes no term.

Every accepted spline has exactly one maximum and one minimum. Topology is checked
using all roots of its piecewise polynomial derivative, classified between roots;
stationary intervals and additional stationary points are rejected. A fit without
an admissible topology fails rather than returning an oscillatory seed.
`FreeMotionDefinition.to_shape_coordinates()` inverts the production chart,
round-trips the canonical controls and explicitly rejects the excluded pole.

Diagnostics compare the naive node projection at phase zero with the fitted
motion: cycle position RMS and sampled maximum error, feature-neighborhood RMS,
optional velocity RMS, event-local errors, extrema angles and count, and analytic
maximum absolute first and second spline derivatives. Position-error maxima use
a dense evaluation mesh, not a certified continuous supremum. `plot_refit()`
shows both pistons, position, error, available velocity, source events, extrema
and all 15 shifted nodes. It produces a static Matplotlib figure, not a GIF.

### Motion-only study generation

`refit_study()` accepts a current study, stored evaluation, or campaign with an
explicit candidate selector. It reconstructs that exact machine and writes a
portable schema-3 study plus machine basis. Only the 28 spline coordinates are
active (13 shapes plus phase per side); all other candidate coordinates are fixed.
Physical constraints and numerical settings are retained. An incompatible
family-specific mechanical constraint is rejected, never silently removed.
A reference-pressure filling policy is converted to explicit inventory using the
selected candidate's computed mass, so changing motion does not refill the machine.

The default motion-only search is deliberately narrow around the refitted law.
The joint canonical-control angular radius is **0.02 radians**. Its stereographic
image is a ball; the study uses a box centered exactly on the fitted coordinates,
contained in that ball. Each coordinate half-width is the initial point's distance
to the ball boundary divided by `sqrt(13)`. Even simultaneous changes at box
corners therefore stay inside the canonical cap. This prevents the 13 independent
shape coordinates from collectively expanding the requested neighborhood.

This is a search-region policy, not a uniform physical displacement bound or a
claim that the source motion is optimal. A cap touching the excluded pole is
rejected. Machine-basis provenance records the policy, cap, box half-width and
pole distance. The independent phase search half-width defaults to **0.48 degrees**,
2% of a node spacing. `--shape-radius` and `--phase-radius-fraction` can explicitly
widen or narrow these local search policies. These bounds do not change the refit
itself or its one-cell phase exploration. They apply only to newly generated
studies; existing studies and campaign identities are not rewritten.
Other active design coordinates are not inherited. The initial refit is evaluated
before ordinary Sobol sampling when the new study is subsequently run.

```console
dada-research motion-refit path/to/campaign --candidate CANDIDATE_ID \
    --output path/to/spline/study.toml --report path/to/spline/refit.json --plot
dada-research validate path/to/spline/study.toml
dada-research run path/to/spline/study.toml --directory path/to/spline/campaign
```

A `MotionTarget` JSON alone supports `--report` fitting and plotting but cannot
provide a complete machine study. `--validate-only` checks the request without
fitting or integration. Existing destination files are not overwritten.
No thermodynamic integration or thermodynamic optimization occurs during refit.

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

### Interactive mechanism inspection

`dada_solver.research.mechanism_view` provides `mechanism_state()`,
`plot_mechanism()`, `animate_mechanism()`, `plot_motion_comparison()` and
`plot_library_catalogue()`. `MechanismModel` reconstructs a portable artifact
once. Slider-crank, four-bar and six-bar frames use production joint states;
the four-bar `joint_state(theta, side)` exposes the same closure used by its
volume law. No drawing-specific closure, phase correction or reflection is
introduced. Six-bar joints A–B–C–D–E–F–G–H–P are all available.

Views combine linkage trajectories, normalized position, available target
velocity, fit error, mechanical metrics and a shared angle cursor. Local lengths
remain in crank-radius units. A solver record may supply additional COP/power
annotations; visualization never evaluates thermodynamics.

Install the optional `dada-engine-solver[plot]` dependencies for these views.
`animate_mechanism()` returns a figure and Matplotlib `FuncAnimation`; retain the
animation and call `matplotlib.pyplot.show()`. Space pauses/resumes. No GIF,
encoder or export step is required. Libraries are explicitly selected by family
ID and piston; no implicit best member is chosen.

```console
dada-research mechanism catalogue path/to/library.json --plot
dada-research mechanism visualize path/to/mechanism.json --target path/to/target.json
dada-research mechanism visualize path/to/library.json --family-id FAMILY_ID --side large
```

Use `--static` for a still view or `--no-show` to construct a view without opening
a window. Neither changes the artifact or starts an integration.

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

Thermodynamic equations, hardware physics and campaign persistence internals
are outside this kinematics reference.
