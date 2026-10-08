# Mechanism-synthesis search for the DADA engine

## Purpose

This note describes a hierarchical method for transforming an optimized
thermodynamic piston motion into mechanically realizable slider-crank, four-bar
and six-bar linkages. Targets, bounds and objective weights belong to the selected study.

The main methodological result is:

> Mechanism synthesis should be treated as a hierarchical basin-discovery and
> thermo-mechanical adaptation problem, not as a blind fit of all linkage
> dimensions at once.

The practical sequence is:

1. identify the piston carrying the most restrictive kinematic demand;
2. discover a diverse set of primary four-bar families;
3. give selected primaries a fair downstream-dyad search;
4. locally release the complete six-bar geometry;
5. mirror a successful mechanism toward the opposite piston as an
   initialization only;
6. locally adapt the mirrored mechanism;
7. release both mechanisms simultaneously and optimize them using the real
   thermodynamic objective;
8. with the mechanisms fixed, re-tune the small set of relevant thermodynamic
   and exchanger dimensions;
9. use a separate fresh-island saturation experiment when the coverage of the
   primary search space itself must be assessed.

The retained primary seed geometries are documented separately in
[`PRIMARY_FOUR_BAR_FAMILIES.md`](PRIMARY_FOUR_BAR_FAMILIES.md), and complete
six-bar family comparisons belong in
[`SIX_BAR_MECHANISM_FAMILIES.md`](SIX_BAR_MECHANISM_FAMILIES.md).

---

## 1. Thermodynamic target and mechanical mechanism are different objects

An optimized thermodynamic motion is a design probe. It is not automatically a
mechanically meaningful trajectory.

A realizable linkage should preserve the thermodynamic effects that matter
while remaining free to smooth, redistribute or reshape features introduced by
the target parameterization itself.

The working rules are:

1. **position is the primary kinematic reference**;
2. velocity is useful for motion topology and cadence;
3. sharp acceleration features should not automatically be fitted as physical
   requirements;
4. noisy derivative information is not permission to discard the underlying
   position law;
5. once a viable mechanism exists, its final quality must be assessed by the
   thermodynamic solver.

This distinction is important because an optimizer can exploit an omitted
position interval by introducing a pause or reversal while still appearing to
match derivative-based targets well.

The target therefore guides discovery, but the final mechanism is not required
to reproduce every feature of the abstract trajectory.

---

## Family-owned software protocols

A periodic `MotionTarget` is the common input for abstract-motion refit and all
three physical synthesis families. Its normalized SMALL/LARGE positions,
available derivatives, extensible source events and scientific identity are
independent of the source representation. Extraction uses production study-angle
kinematics, without thermodynamic integration; missing derivatives remain absent.
The [Research kinematics reference](DADA_ENGINE_RESEARCH_KINEMATICS.md#11-motion-targets-and-synthesis-handoff)
owns the target conventions, artifact APIs and CLI.

The central `synthesis_protocol()` registry separates family semantics:

| Family | Ordered stages | Continuous coordinates |
|---|---|---|
| `slider_crank` | `global_discovery`, `full_local_polish`, `paired_thermodynamic`, `hardware_retuning` | Rod/crank, offset/crank, phase; orientation fixed |
| `four_bar` | `global_discovery`, `full_local_polish`, `paired_thermodynamic`, `hardware_retuning` | Eleven production coordinates; branches, output and orientation fixed |
| `six_bar` | Primary discovery, downstream fit, complete polish, mirror initialization, opposite adaptation, paired thermodynamics, hardware retuning | Six primary + nine downstream per side |

Slider-crank and four-bar have no downstream-dyad stage. The hierarchical method
below applies to six-bar, not to an artificial common topology. Categories remain
fixed scientific choices rather than continuous optimization coordinates.

`SynthesisPlan` validates family stages and target identity, creates existing
`MechanismArtifact` objects from submitted physical geometries, and records
separate evidence in `MechanismLibrary` members. The library may retain mixed
physical families and multiple basins of the same family. It never prunes to one
winner automatically. Direct slider-crank and four-bar discovery and polish are
implemented, along with the five geometric six-bar stages described below.

Fit, mechanical quality and thermodynamic performance are three distinct levels.
Fit uses angle-weighted position RMS/max error and optional velocity RMS; it does
not impose acceleration fitting. Production mechanical metrics and constraint
margins remain individually inspectable. COP, cooling power or indicated power
come only from an identified solver result, not from a combined proxy score.
Dimensional volume-derivative constraints require an explicit machine volume
scale; normalized portable geometry alone does not supply one.

Interactive Matplotlib views use production joints for linkage animation,
point trajectories, target/mechanism position and available velocity, shared
angle cursors, fit error and mechanical metrics. `FuncAnimation` and `show()` are
the normal path; no GIF is required. An unranked catalogue exposes family ID,
mechanism family, evidence, provenance and artifact hash for human selection.
Selected families can subsequently be polished or adapted using the actual
Research thermodynamic objective.

The abstract-motion refit operator converts a `MotionTarget` to
`structured_c2_15p` with fifteen active pair coordinates. A deterministic small
multistart position-led fit uses hybrid features as seeds and diagnostics, not
fixed constraints. Continuous polynomial derivative checks reject nonmonotone
solutions. Target acceleration is not fitted. The generated study freezes all
other machine coordinates. This operator does not perform mechanical synthesis
or judge thermodynamic performance. See the
[refit contract](DADA_ENGINE_RESEARCH_KINEMATICS.md#feature-aware-motion-refit)
for numerical policy, diagnostics and study generation.

### Direct slider-crank and four-bar operators

`SynthesisPlan.execute()` implements `global_discovery` and `full_local_polish`
for `slider_crank` and `four_bar`, using the shared `GeometrySearch` engine.
Fresh-primary saturation is not implemented by this geometric operator.
Paired adaptation and hardware retuning generate Research studies via
`mechanism adapt` and `mechanism retune`, as described below. No geometric stage
constructs or calls a thermodynamic evaluator.

The versioned `SearchPolicy` is `direct_geometry_islands_v1`; its default bounds
are `crank_normalized_direct_bounds_v1`. Lengths are in crank-radius units.
These are **search bounds**, not physical domains or manufacturing limits:

| Family | Coordinate | Default interval |
|---|---|---|
| Slider-crank | `rod_over_crank` | 1.2–12 |
| Slider-crank | `offset_over_crank` | −4–4 |
| Both | `phase_rad` | 0–`2*pi` |
| Four-bar | `coupler`, `rocker` | 1.2–8 |
| Four-bar | `ground_x`, `ground_y` | −6–6 |
| Four-bar | `output_along` | −6–6 |
| Four-bar | `output_normal` | −4–4 |
| Four-bar | `rod_length` | 1.5–16 |
| Four-bar | `slider_origin_x`, `slider_origin_y` | −8–8 |
| Four-bar | `axis_angle` | −`pi`–`pi` |

Discovery explicitly explores both slider-crank orientations and all 32 four-bar
combinations of output (`rocker`/`coupler`), loop branch, slider branch, crank
direction and volume orientation. Discretes remain separate categories, never
continuous coordinates. The default is four independent islands per category
and piston, population 16 and 24 generations. Every island starts with one
phase-aligned generic geometry and independent random samples, not source
mechanism dimensions. Differential evolution uses `rand/1/bin`, mutation factor
uniform in 0.5–0.9 and crossover probability 0.8; populations are interleaved
across categories, islands and pistons. Slider-crank uses the same small
infrastructure in three dimensions. Discovery then locally polishes retained
basins (24 least-squares iterations by default).

Seeds derive deterministically from the root seed, piston, category and island.
Fixed settings and evaluation counts reproduce identities in the same numerical
runtime. A cooperative wall deadline is machine-dependent; its observed stopping
point is not a promise of identical wall-budget replay. `--max-evaluations` caps
actual geometry evaluations, including local finite differences. Final topology
validation and HTML rendering occur after the search budget. Evidence records
seeds, island progress, rejection stages, actual evaluations and policy.

Geometry is rejected progressively: full crank-revolution loop closure,
piston-rod closure, nondegenerate stroke/singularity, two reversals, declared
mechanical constraints, then dense fit. Production closure and metrics are the
authority. Vectorized poses share the scalar production equations. Slider-crank
turnarounds use exact collinear dead centers. Four-bar acceptance refines roots
of analytic velocity on doubled periodic grids and screens tangencies; this is
a resolution-checked screen, **not a certified continuous root count**. No
implicit transmission floor or auxiliary mechanical preference is imposed.

The fit has 721 uniform cycle samples by default. Nine samples per recognized
extremum, turnaround, cadence join, rounding or kink reinforce coverage; their
combined weight is only 0.1 relative to the background weight 1. Directed events
cover their interval; point events use a radius of `2*pi/36`. These configurable
numerical policies prevent neglected intervals without requiring reproduction
of fine abstract-motion features. Available target velocity is secondary;
missing velocity stays absent and target acceleration is never used.

The internal discovery score is weighted position MSE plus
`0.001 * normalized_velocity_MSE / (1 + normalized_velocity_MSE)`. Velocity is
normalized by target RMS velocity, bounded below by `1/(2*pi)`; its contribution
is bounded by 0.001. Both terms are reported separately. `assess_mechanism()`
independently reports cycle position RMS/max error, optional velocity RMS,
mechanical metrics and individual constraint margins. None is thermodynamic
performance. Dimensional volume-derivative constraints require explicit
`volume_limits`; normalized geometry does not invent a cylinder size.

Diversity uses RMS differences in normalized geometry, with shortest circular
phase distance divided by `2*pi`. Four-bar descriptors remove common rigid
frame rotation and axial slider-origin gauge; categories are compared separately.
A score-ordered greedy archive retains representatives at distance at least
0.05, with category round-robin retention and up to 64 members per piston by
default. Bounds, threshold and retention are configurable. This direct-mechanism
archive is distinct from the primary-loop connected-component saturation method
below; its count is not an exhaustive geometric family count. Multiple basins
survive; the catalogue never declares an absolute winner.

`full_local_polish` requires explicit library family IDs. Each parent gets its
own bounded least-squares search (80 iterations by default), within 10% of each
configured coordinate range. Full-period angular windows wrap at the seam.
Output, branches, direction, orientation, settings and embedded constraints are
preserved. One improved artifact is retained per selected parent and piston,
with parent family/artifact identity in provenance. It does not merge parents
or select a global winner.

```console
dada-research mechanism synthesize path/to/target.json --family slider_crank \
    --stage global_discovery --side both --seed 1234 --islands 4 \
    --budget 2m --output path/to/slider-library.json
dada-research mechanism synthesize path/to/target.json --family four_bar \
    --stage global_discovery --side both --seed 1234 --islands 8 \
    --budget 2m --output path/to/fourbar-library.json
dada-research mechanism synthesize path/to/target.json --family four_bar \
    --stage full_local_polish --side large --library path/to/fourbar-library.json \
    --family-id large-family-003 --output path/to/fourbar-polished.json
```

Every execution writes a `MechanismLibrary` and a standalone sortable HTML
catalogue beside it; `--html` selects another path. Install the optional `[plot]`
dependencies. The catalogue provides family IDs, production linkage views and
point paths, target/mechanism position and velocity, fit errors, mechanical
metrics, categories and full evidence. Select a member for interactive animation:

```console
dada-research mechanism visualize path/to/fourbar-library.json \
    --family-id large-family-003 --side large --target path/to/target.json
```

`--config` accepts a TOML `SearchPolicy`, with CLI flags overriding the matching
fields. For example:

```toml
seed = 1234
islands = 4
population = 16
generations = 24
cluster_distance = 0.05
retain_per_side = 64

[bounds]
rod_over_crank = [1.5, 10.0]

[categories]
volume_increases_with_coordinate = [true, false]

[[mechanical_constraints]]
metric = "minimum_rod_axis_cosine"
relation = "minimum"
limit = 0.5
unit = "1"
```

Bounds and categories must belong to the requested family. Constraints are
unscoped declarations applied to each selected piston. Optional
`[volume_limits.small]` and `[volume_limits.large]` tables require `minimum` and
`maximum` in cubic metres. Existing output files are not overwritten.
`--validate-only` checks family protocol and policy without search or integration.
The source can be a `MotionTarget`, study, evaluation or identified campaign
candidate, including a spline source; no hybrid-specific synthesis path exists.


---

### Executable hierarchical six-bar operators

The `hierarchical_six_bar_geometry_v1` adapters reuse **the same**
`GeometrySearch.discover()` differential evolution, archive, budget control and
`GeometrySearch.polish()` bounded least-squares implementation. There is no
15-dimensional global six-bar search and no thermodynamic integration.
Execute one stage at a time and retain explicit `--family-id` parents.

| Stage | Released geometry | Categories explored or preserved |
|---|---|---|
| `primary_discovery` | Six primary coordinates | Two primary branches |
| `downstream_fit` | Nine downstream coordinates, primary fixed | Two second branches per selected primary |
| `full_local_polish` | All fifteen coordinates locally | Both parent branches fixed |
| `mirror_initialization` | No search | Branches transformed by the reflection |
| `opposite_local_adaptation` | Fifteen coordinates on the destination side only | Destination branches fixed |

A primary is a `MechanismArtifact` with settings
`family = "six_bar", component = "primary"`, exactly the six primary coordinates
and `primary_branch`. It reconstructs `SixBarPrimaryMechanism`, not a piston law;
it cannot be injected into a machine until a downstream linkage exists.
The existing `MechanismLibrary` retains these intermediate artifacts and their
family IDs without inserting dummy slider or dyad dimensions.

Primary discovery uses the principal-axis projection of the production E
trajectory as an **intermediate chronology opportunity**, following method
[6.4](#64-topology-and-cadence-proxy). Both projection signs are examined and
selected by topology and cadence, never by piston-position resemblance.
The executable policy is `primary_topology_cadence_v2`: correctly timed
turnarounds, monotonic signs, zero endpoint speeds, extra-reversal penalties,
short-branch speed ratio and displacement sharing, and explicit soft
transmission/directionality preferences. Optional long-branch symmetry compares
the candidate's own branch; it does not penalize internal speed variations.

The catalogue labels the projection as such: **E is not P**. Position and
velocity-profile RMS remain diagnostics only; piston fit is `null`.
Evidence exposes projection axis/span, axis variance fraction, primary
transmission, matched turnarounds, cadence descriptors and each weighted score
component. Constraints requiring downstream geometry are recorded as deferred,
not satisfied; applicable primary constraints are enforced. Downstream fitting
continues to compare the actual piston trajectory directly with the target.

Default primary search bounds, in crank-radius units except phase, are:
`primary_ground` 1.2–10, `primary_coupler` 1.2–12,
`primary_rocker` 1.2–10, `primary_e_along` −8–12,
`primary_e_normal` −8–8 and `primary_phase` 0–`2*pi`.
These are numerical search bounds, not new mechanical domains.

Downstream bounds are derived separately from each parent's E trajectory.
Let `d` be its bounding-box diagonal (floored by `minimum_trajectory_extent`, default 0.1 crank radius),
`c` the bounding-box center, and `r = pivot_envelope_radius*d` (default 2d).
G is searched in the coordinate box `c +/- r`; EF/GF span
`link_extent_minimum*d` through `link_extent_maximum*d` (defaults 0.25d–4d).
H coordinates span −1–2 along EF and −1–1 normal to EF.
The rod spans 0.25d–`rod_extent_maximum*d` (default 8d);
axis offset spans `+/-(norm(c)+r)`, and axis angle spans −`pi`–`pi`.
All effective bounds are stored; `[bounds]` overrides take precedence.
These relative-bound settings are configurable numerical search policies.
Local stages inherit the effective parent bounds, including the primary bounds;
mirrored seed regions are reflected and phase-centered consistently. Explicit
configuration bounds override this inheritance.

Structured downstream seeds place G around E's envelope, choose equal EF/GF
lengths exceeding half the maximum E–G distance with an extent-based margin,
place H on/near EF, and initialize the slider from H's principal axis.
The rod exceeds the transverse excursion. Both axis directions are considered
as seeds using the production closure. These heuristics supply starting points;
they impose no equality, rod/stroke ceiling or transmission floor.
The evaluation order checks primary, secondary and rod closure, stroke,
piston topology, explicit mechanical constraints, then dense piston fit.
Complete piston acceptance refines analytic velocity roots on doubled grids and
screens tangencies: exactly one maximum and one minimum are required. This is
resolution-checked, not a certified continuous proof.

Clustering uses RMS coordinate differences divided by effective search widths;
phase and axis-angle distances are circular. Primary discovery uses six primary
coordinates; downstream uses nine downstream coordinates **within one parent**.
Complete geometry balances primary and downstream group RMS equally.
Categories and parent lineages never merge. The configurable default distance
0.05 remains a normalized search-space distance, not a physical similarity law.
Polish retains one result per selected parent rather than one global winner.
Multiple selected parents receive fair shares of the remaining evaluation/time
budget. A parent with no admissible descendants is recorded in
`failed_parent_searches` and does not erase successful descendants from others.
Seeds, actual bounds, evaluations and parent artifact hashes are retained.

Mirroring reflects the local y geometry, reverses both branch signs and crank
phase, giving `q_mirror(theta) = q_source(-theta)` before phase alignment.
It aligns the destination maximum with its own target and retains the original
side artifact verbatim. The resulting pair contains independently owned geometry.
Opposite adaptation releases only the destination geometry; there is no mirror
constraint. Both sides are selectable; neither is assumed restrictive.
A local window is 10% of configured coordinate width by default, with periodic
wrapping. `local_radius` and `polish_evaluations` configure this search region.

Every stage writes a library and the shared sortable HTML catalogue. Primary
views show A–B–C–D/E and the E path beside the target/projection; full views use
production A–B–C–D–E–F–G–H–P states. Metrics, dimensions, categories, hashes and
lineage remain inspectable, and each artifact works with `mechanism visualize`.
Neither the catalogue nor mirroring infers thermodynamic performance.

```console
dada-research mechanism synthesize path/to/target.json --family six_bar \
    --stage primary_discovery --side large --islands 12 --budget 5m \
    --output path/to/primary-large.json
dada-research mechanism synthesize path/to/target.json --family six_bar \
    --stage downstream_fit --side large --library path/to/primary-large.json \
    --family-id primary-L-001 --family-id primary-L-004 --budget 5m \
    --output path/to/downstream-large.json
dada-research mechanism synthesize path/to/target.json --family six_bar \
    --stage full_local_polish --side large --library path/to/downstream-large.json \
    --family-id primary-L-001/downstream_fit-001 --output path/to/polished-large.json
dada-research mechanism synthesize path/to/target.json --family six_bar \
    --stage mirror_initialization --side small --library path/to/polished-large.json \
    --family-id primary-L-001/downstream_fit-001/full_local_polish-001 \
    --output path/to/small-seed.json
dada-research mechanism synthesize path/to/target.json --family six_bar \
    --stage opposite_local_adaptation --side small --library path/to/small-seed.json \
    --family-id primary-L-001/downstream_fit-001/full_local_polish-001/mirror-small \
    --output path/to/pair.json
```

IDs above illustrate lineage; select the actual IDs from each catalogue.
Paired thermodynamics, hardware retuning and fresh-primary saturation remain
separate, non-executed operators.

## 2. Six-bar topology and variables

The retained six-bar topology consists of:

- primary four-bar `A-B-C-D`;
- point `E` rigidly attached to coupler `BC`;
- secondary dyad `E-F-G`;
- point `H` rigidly attached to link `EF`;
- a finite piston rod connecting `H` to a slider axis.

The crank length is normalized:

\[
AB = 1.
\]

Each cylinder mechanism has 15 continuous coordinates.

### Primary four-bar

\[
(AD,\ BC,\ CD,\ E_\parallel,\ E_\perp,\ \phi)
\]

where the first three lengths and the point-E coordinates are expressed in
crank-radius units.

### Downstream mechanism

\[
(G_x,\ G_y,\ EF,\ GF,\ H_\parallel/EF,\ H_\perp/EF,
L_{rod},\ s_{slider},\ \alpha_{slider})
\]

There are also two discrete assembly branches:

- primary four-bar branch;
- secondary dyad branch.

A paired two-cylinder machine therefore contains 30 continuous mechanism
coordinates once the discrete branches have been chosen.

The production representation of this topology is implemented by
`SixBarCylinderMechanism`; the search-stage coordinate groups are exposed by
`dada_solver.research.synthesis`.

---

## 3. Why a blind 15-dimensional search is inefficient

A direct global search over all 15 continuous dimensions has several structural
disadvantages:

- the feasible region is thin relative to the raw parameter box;
- loop closure, transmission and slider constraints create large invalid
  regions;
- flat invalid penalties provide little information to a global optimizer;
- useful structure exists already in the primary four-bar;
- the downstream dyad can substantially transform the primary motion;
- randomizing primary and downstream geometry simultaneously spends a large
  fraction of the search budget rediscovering basic feasibility.

The useful decomposition is therefore hierarchical.

This is not merely a software optimization. It reflects the physical roles of
the two stages: the primary produces a temporal and planar motion structure,
while the downstream dyad acts as a genuine motion transformer.

---

## 4. Hierarchical synthesis workflow

The generic stage names used by `dada_solver.research.synthesis` are retained
below.

### 4.1 `primary_discovery`

Search only the six primary coordinates:

\[
x_p=(g,c,r,E_\parallel,E_\perp,\phi).
\]

The goal is not a final piston mechanism.

The primary should provide a useful temporal cadence and a point-E trajectory
that a downstream dyad can exploit.

Because this problem has only six continuous dimensions, many independent
global search islands are practical.

Several mechanically distinct families should be retained rather than only the
lowest-scoring primary.

### 4.2 `downstream_fit`

For each selected primary, freeze its six coordinates and search the remaining
nine continuous coordinates plus the secondary assembly branch.

The downstream dyad is deliberately allowed to reshape the primary motion.

A primary-only score is therefore not a reliable predictor of downstream
transformability.

### 4.3 `full_local_polish`

Once a viable primary and downstream dyad exist, release all 15 continuous
coordinates in a moderate local region.

This stage is important because modest primary changes can allow large
improvements in the complete piston motion that are impossible with the
primary frozen.

### 4.4 `mirror_initialization`

The mechanism synthesized for the more restrictive piston may be spatially
mirrored and phase-aligned toward the opposite piston.

The mirror is an **initialization only**.

It must not become a permanent symmetry constraint unless a separate physical
argument justifies that restriction.

### 4.5 `opposite_local_adaptation`

The mirrored mechanism receives its own independent local adaptation while its
chosen discrete branches are retained.

The two final cylinder mechanisms are therefore allowed to become different.

### 4.6 `paired_thermodynamic`

Once realizable mechanisms exist on both sides, release the two complete
mechanisms simultaneously.

At this point the optimization objective changes fundamentally:

- the actual thermodynamic objective becomes primary;
- position error relative to the abstract target becomes diagnostic;
- mechanism-family identity is preserved through a local search region and
  fixed discrete branches rather than through exact geometric imitation.

The abstract target has served its purpose once it has led the search into a
useful realizable basin.

This stage is implemented for `slider_crank`, `four_bar` and `six_bar` by
`paired_thermodynamic()` and the `mechanism adapt` command. It requires an explicit
reconstructible Research source and one library member containing both pistons.
The source may be a hybrid study/candidate directly. Geometry discovery is not
repeated and no thermodynamic integration starts during generation.

```console
dada-research mechanism adapt path/to/source-campaign --candidate best \
  --library path/to/pair.json --family-id FAMILY_ID \
  --output path/to/adaptation/study.toml --radius 0.1
dada-research validate path/to/adaptation/study.toml
dada-research run path/to/adaptation/study.toml
```

The unchanged mechanism artifacts supply exact initial coordinates and fixed
categories. The generated study releases 6, 22 or 30 continuous coordinates for
the respective paired families, freezes all hardware and other machine coordinates
at the source candidate, and reuses its objectives, constraints and numerical
policies. Reference-pressure charge is explicitly frozen to the selected inventory.
Source and artifact mechanical constraints remain jointly enforced, including
production closure and refined two-reversal topology screens before integration.

The existing center-first `local_regions_v1` Sobol policy explores a normalized
radius (default 0.1) about the exact pair. Its centered reference boxes use synthesis
default coordinate widths, capped at half the current value for positive lengths;
angles remain unwrapped locally. This `centered_mechanical_local_regions_v1`
policy is a configurable search region, not a mechanical domain. The
[paired study reference](DADA_ENGINE_RESEARCH_KINEMATICS.md#paired-thermodynamic-study-generation)
specifies bounds, compatible warm starts, immutable artifact linkage and recorded
source/pair provenance. Fit remains diagnostic; ordinary Research evaluates the
actual source objective without a COP/RMS score.

### Manual selection between synthesis and adaptation

A catalogue discovers basins; the human chooses SMALL and LARGE independently.
`mechanism pair LIBRARY --small SMALL_ID --large LARGE_ID --output pair.json`
selects from a common library. For mixed sources, use
`mechanism pair --small-library SMALL_FILE --small SMALL_ID --large-library LARGE_FILE --large LARGE_ID --output pair.json`.
Do not combine the two forms. This operation preserves complete artifact payloads
and hashes and records independent `manual_pair_selection` provenance.

Homogeneous workflow: one catalogue → select SMALL and LARGE → pair → adapt.
Heterogeneous workflow: separate catalogues → select one side in each → pair →
adapt. No family is preferred and no combination search is performed. The shared
HTML renderer provides short-command selection, copyable side selectors and
direct adaptation snippets for already complete pairs. Primary-only components
remain visible but cannot be paired. See the
[manual pairing reference](DADA_ENGINE_RESEARCH_KINEMATICS.md#manual-smalllarge-pairing)
for exact syntax and validation.

`paired_thermodynamic` releases coordinates from each artifact's own family:
3, 11 or 15 per side. Mixed counts are 14 (slider/four-bar), 18 (slider/six-bar)
and 26 (four-bar/six-bar). Subsequent hardware retuning freezes either kind of pair.

### 4.7 `hardware_retuning`

`hardware_retuning()` and `mechanism retune` implement this final stage for
slider-crank, four-bar and six-bar. They freeze the exact adapted SMALL/LARGE
mechanisms and reopen selected thermodynamic or exchanger coordinates. No
mechanism coordinate is active and generation never integrates a thermal cycle.

```console
dada-research mechanism retune path/to/paired-campaign --candidate best \
  --scope source-active --radius 0.1 --output path/to/retuning/study.toml
dada-research validate path/to/retuning/study.toml
dada-research run path/to/retuning/study.toml
```

The default source-active set restores only non-kinematic variables that were
active in the original Research source. Their declarations are stored portably
by paired adaptation, with their original domains and types. Initials come from
the exact current candidate, not the abstract source. Numeric local intervals
are clipped to the original normalized domain through the existing local Sobol
scheduler; categorical choices keep their original declared set. No new hardware
domain, capacity scaling or charge-policy conversion is introduced.

Repeated `--group exchangers`, `--group volumes`, `--group frequency`,
`--group charge`, `--group valves`, `--group external-stream`, and
`--parameter NAME` filter this set without adding absent/fixed source parameters.
Later passes can reopen additional groups using the original domains. Immutable
artifacts describe the adapted geometry with correct hashes and parent lineage;
mechanical constraints continue to apply when volume or other hardware changes.
The current source supplies objectives, evaluation policies and compatible warm
guesses. See the [hardware retuning reference](DADA_ENGINE_RESEARCH_KINEMATICS.md#hardware-retuning-with-fixed-mechanisms)
for domain recovery, bounds, typed choices and complete provenance.

This separates:

1. performance lost because of mechanical realization;
2. performance recoverable by adapting the surrounding machine to the new
   mechanism.

---

## 5. Primary four-bar search space

Normalize the crank to \(AB=1\) and search the six coordinates
\((g,c,r,E_\parallel,E_\perp,\phi)\). Declare the coordinate bounds and admissible
assembly branches for the study. Full-revolution closure is mandatory.

Primary transmission quality is measured by

\[
s_{\min}
=
\min_\theta
\frac{|BC\times DC|}{|BC||DC|}
=
\min_\theta |\sin\mu|.
\]

Declare any hard transmission floor separately from a soft preference. Likewise,
a required point-E span must be an explicit design requirement, not an implicit
universal threshold.

> A mechanical preference should not silently become a topological exclusion
> unless a real physical limit justifies it.

---

## 6. Primary-search proxy limitations

A mathematically plausible primary objective can be mechanically misleading.
The following distinctions constrain how a proxy should be interpreted.

### 6.1 Acceleration-lobe timing is not enough

The tangential acceleration associated with point-E speed magnitude is

\[
a_t
=
\frac{E'\cdot E''}{\|E'\|}
=
\frac{d\|E'\|}{d\theta},
\]

The derivative of a non-negative speed magnitude is not signed piston
acceleration.

### 6.2 Signed PCA velocity improves sign information but not physical identity

Projecting point-E motion onto its principal PCA axis supplies signed velocity,
but has two limitations:

- the PCA axis is not the final piston axis because the downstream dyad can
  rotate and transform the motion;
- event-based penalties can allow a nominally required event to become
  arbitrarily weak without a correspondingly large penalty.

### 6.3 Statistically optimal segmentation can violate motion topology

Automatic piecewise velocity segmentation can produce a low approximation
error while allowing a fitted regime to cross a real piston reversal:

> A segmentation may minimize approximation error while violating the topology
> of the physical motion.

Real turnarounds must therefore structure the cadence objective.

### 6.4 Topology-and-cadence proxy

Begin from the selected target's physical turnarounds. Measure its branch
durations, fast/slow mean-speed ratio and displacement sharing rather than
assuming one operating point's cadence is universal.

The proxy rewards:

1. two correctly timed turnarounds;
2. no additional reversal;
3. correct monotonic sign on both branches;
4. target-like fast/slow cadence on the short branch;
5. target-like displacement sharing on that branch;
6. acceptable primary transmission;
7. weak geometric directionality preferences.

A study may leave the return branch less constrained instead of copying its
detailed target velocity profile. If approximate mirror symmetry is chosen as
a study preference, the candidate's own branch asymmetry can be measured by:

\[
A_{BP}
=
\frac{
\operatorname{RMS}[v(u)-v(1-u)]
}{
\sqrt{
\operatorname{mean}
\left[
\frac{v(u)^2+v(1-u)^2}{2}
\right]
}
}.
\]

A symmetric fast-slow-fast branch is therefore allowed.

The implementation extracts extrema from `MotionTarget` events when present,
otherwise from its periodic position interpolant. The shorter ascending or
descending branch is selected cyclically. Interior `kink`/`cadence_change` events
provide candidate splits. Without such events, a deterministic two-level split
of available target velocity is accepted only when its mean-speed contrast and
explained variance meet the configured detection thresholds. Event-supported
splits also require a resolved speed contrast. No target acceleration is used.
If cadence cannot be resolved, `missing_cadence = "disable"` omits both cadence
terms and records the reason; `"error"` refuses discovery. No cadence is invented.

A search-settings TOML may contain, for example:

```toml
[primary_cadence]
missing_cadence = "disable"
monotonicity_weight = 2.0
endpoint_zero_weight = 0.40
turning_weight = 0.65
turning_scale_rad = 0.2617993877991494
extra_crossing_weight = 0.50
speed_ratio_weight = 0.55
displacement_fraction_weight = 0.80
long_symmetry_weight = 0.0
transmission_preference_weight = 0.12
directionality_preference_weight = 0.03
preferred_transmission_sine = 0.35
preferred_axis_variance_fraction = 0.70
minimum_speed_ratio = 1.5
minimum_split_explained_variance = 0.5
reject_extra_turnarounds = false
```

All weights and preferences belong to `PrimaryCadencePolicy`, not universal
mechanical limits. Long-branch symmetry is disabled unless requested explicitly.
The scalar score sums the weighted RMS wrong-sign velocity, normalized endpoint
speed, turnaround RMS divided by `turning_scale_rad`, extra-reversal count,
absolute log speed-ratio error, absolute displacement-fraction error, optional
mirror asymmetry and normalized deficits below the two preferred mechanical
values. Speeds in topology terms are normalized by the candidate's cycle RMS.
Non-reversing stationary points are not extra turnarounds.

Default sampling guards exclude 3/360 of a cycle around turns for sign checks,
8/360 for cadence means and long symmetry, and 6/360 around a cadence split;
they are capped relative to branch/subphase length. Candidate cadence means
use 65 points in each guarded subphase. A numerical split must leave at least
15% of the short branch on either side. These configurable sampling/detection
policies are recorded alongside effective weights and search/root resolutions.
Analytic production velocities and refined root searches determine turns;
this is a doubled-grid/tangency check, not a certified continuous proof.



This is still only a search proxy.

The physical question is ultimately whether the complete cylinder-volume
evolution produces good thermodynamic behavior.

---

## 7. Primary score ownership

A primary-discovery proxy may combine monotonicity, velocity at turnarounds,
zero-crossing timing, additional reversals, fast/slow mean-speed ratio,
displacement sharing and transmission quality. The study must declare the
weights, reference descriptors, hard constraints and any directionality or
symmetry preferences separately. Evaluate both possible PCA-axis signs when
using signed projected velocity.

The score should not be interpreted as:

- thermodynamic efficiency;
- complete-mechanism quality;
- manufacturing robustness;
- geometric probability density.

It exists only to make primary family discovery computationally tractable.

---

## 8. Search islands and basin-capture probability

A differential-evolution island is not one small geometric patch.

A fresh island initializes a population across the complete bounded domain,
according to its declared population policy. The population then evolves
globally to locally.

The useful statistical quantity is therefore **algorithmic basin-capture
probability**:

\[
P_i^{capture}
=
P(\text{one fresh search island terminates in family }i).
\]

For fixed:

- bounds;
- constraints;
- optimizer;
- population policy;
- evaluation budget;
- branch policy;

it can be estimated by

\[
\hat P_i=\frac{n_i}{N}.
\]

This must not be described as the literal geometric volume of a mechanism
family in continuous parameter space.

---

## 9. Operational definition of a primary family

For fresh-island saturation, use the connected-component definition below.
Executable design-exploitation discovery instead uses the score-ordered archive
and the `2*pi` circular scale specified above. Their family counts and distance
thresholds must not be treated as interchangeable.

Opposite primary assembly branches are always treated as different families.

For mechanisms on the same branch, define normalized coordinate differences.

For the five non-periodic coordinates,

\[
z_i
=
\frac{|x_i-y_i|}{x_{i,\max}-x_{i,\min}}.
\]

For phase,

\[
z_\phi
=
\frac{
|\operatorname{wrap}_{[-\pi,\pi)}(\phi_x-\phi_y)|
}{\pi}.
\]

The geometric distance is

\[
d(x,y)
=
\sqrt{
\frac{1}{6}
\sum_i z_i^2
}.
\]

If the branches differ,

\[
d=\infty.
\]

Families are connected components of the graph containing an edge whenever

\[
d < d_{\mathrm{threshold}}.
\]

Using connected components makes the result independent of insertion order,
although chaining remains a reason not to over-interpret the absolute number of
families.

Declare the clustering threshold with the search policy and assess sensitivity
to that choice before interpreting family counts.

---

## 10. Fresh-island saturation evidence

Assess saturation with independent fresh islands, fixed bounds and optimizer
budget, balanced assembly branches and no injection of design candidates.
Record the family-count distribution under the declared clustering policy.
Singletons and repeated captures can inform estimates of unseen algorithmic
capture mass, but do not certify exhaustive geometric coverage.

Concentration of capture probability in a few basins does **not** prove that:

- every geometrically possible four-bar family has been discovered;
- family counts are invariant to the clustering threshold;
- the largest capture basin produces the best complete six-bar;
- capture frequency is geometric phase-space volume.

---

## 11. Seeded design search and fresh-island saturation have different purposes

Mechanism seeds are useful during design exploitation
because they preserve expensive discoveries and accelerate convergence.

They are inappropriate when the scientific question is whether a fresh search
repeatedly rediscovers the same primary basins.

The two protocols must therefore remain explicit and separate.

### Design exploitation

May use:

- retained families;
- warm starts;
- local continuation;
- selected mechanisms.

Question:

> How efficiently can a useful design be improved?

### Fresh-island saturation

Uses:

- fixed broad bounds;
- independent deterministic random seeds;
- no design-candidate injection;
- balanced discrete branches;
- fixed optimizer policy and budget.

Question:

> Under this declared search policy, how concentrated is fresh-search capture
> among the observed families?

This distinction is represented explicitly by `SynthesisRequest` and
`FreshIslandPolicy`.

---

## 12. Downstream transformability matters

A primary four-bar is an intermediate motion generator, not the final piston
mechanism.

The downstream dyad can:

- rotate the effective motion direction;
- reshape timing;
- amplify or attenuate parts of the primary trajectory;
- change the relation between primary cadence score and piston-position fit.

Consequently:

> A mediocre or rare primary under the simplified primary proxy can still
> become valuable after downstream optimization.

This is why several mechanically different primary families should survive long
enough to receive a fair downstream search.

The four retained reference seeds are documented in
[`PRIMARY_FOUR_BAR_FAMILIES.md`](PRIMARY_FOUR_BAR_FAMILIES.md).

---

## 13. Synthesize the more restrictive piston first

For a paired machine, the most useful first target is generally the piston with
the sharper or more restrictive kinematic demand.

The recommended rule is:

> Identify the harder motion, synthesize that side first, then use its
> mechanism as an initialization for the easier side.

The identity of the harder side is not universal.

It may change with:

- motor versus refrigeration operation;
- temperature levels;
- valve timing;
- pressure ratio;
- thermodynamic target family.

The rule is therefore “hard side first”, not “always synthesize one named
cylinder first”.

---

## 14. Mirroring is initialization, not symmetry

A successful mechanism may be mirrored geometrically toward the second
cylinder to provide a strong starting point.

Afterward the second mechanism must be free to adapt.

Permanent symmetry would unnecessarily remove mechanically and
thermodynamically useful degrees of freedom.

The generic synthesis interface therefore records:

```text
mirror_policy = "initialization_only"
```

rather than enforcing equality between the two cylinder mechanisms.

---

## 15. Mechanical constraints and proxy constraints are different

Direct mechanical protections and search preferences must not be confused.

Examples of direct mechanical requirements include:

- full loop closure;
- avoiding toggles below an explicitly retained transmission limit;
- finite piston-rod closure;
- required motion topology;
- any declared packaging or mechanical-clearance limit.

Examples of proxy preferences include:

- favored transmission margin above the hard minimum;
- PCA directionality;
- cadence resemblance;
- geometric compactness preferences.

A useful rule is:

> Relax proxy constraints when evidence shows that they are artificial, but do
> not casually relax direct mechanical protections.

Study-specific mechanical limits belong to the study configuration or complete
mechanism assessment, not to an immutable universal score.

---

## 16. Primary score, piston fit and thermodynamic quality are different rankings

Synthesis distinguishes three different notions of quality.

### Primary quality

A cheap score used to discover four-bar basins with useful motion topology and
cadence.

### Complete-mechanism kinematic quality

How well the full realizable linkage produces an appropriate piston motion
after the downstream dyad and slider are included.

### Thermodynamic quality

How well the actual machine performs when driven by that realizable mechanism.

These rankings need not agree.

A primary that is not the best proxy scorer may be highly transformable.

A complete mechanism that most closely fits the abstract piston curve need not
give the best thermodynamic result.

This is the reason for retaining several families through the hierarchy.

---

Complete-mechanism fitting must score the actual finite-rod piston/slider
displacement, not rocker angle or an intermediate point projection. Failure of
one fitted geometry or one bounded search does not establish a limitation of
the complete mechanism family.

## 17. Why position error becomes diagnostic

During early synthesis, normalized piston-position error is useful because it
provides a cheap, interpretable target.

It should not remain the final objective.

Once a viable mechanism exists:

- the thermodynamic solver can directly judge the machine;
- the mechanism may beneficially smooth target artifacts;
- small deviations from the abstract piston law may improve real performance;
- exact curve imitation can become an unnecessary constraint.

The final optimization therefore uses the study thermodynamic objective while
position error remains diagnostic.

The target is a map toward useful geometry, not a trajectory that the final
machine must copy exactly.

---

## 18. Recommended reusable workflow

For a new DADA machine configuration:

### Step 1 — obtain a useful thermodynamic reference motion

Optimize the thermodynamic machine before imposing a specific realizable
linkage.

### Step 2 — identify trustworthy motion structure

Determine:

- turnarounds;
- monotonic branches;
- important cadence changes;
- derivative regions that are physically meaningful;
- derivative regions that are parameterization artifacts.

### Step 3 — identify the more restrictive piston

Begin mechanism synthesis on that side.

### Step 4 — discover primary families

Search the six primary coordinates over broad bounds and both assembly
branches.

Retain diversity, not only the best scalar score.

### Step 5 — search the downstream dyad

Give several distinct primary families a fair downstream search.

### Step 6 — release the full mechanism locally

Allow all 15 continuous dimensions to adapt around successful constructions.

### Step 7 — initialize the opposite piston

Mirror and phase-align a good mechanism as a starting point only.

### Step 8 — adapt the opposite mechanism independently

Do not enforce permanent symmetry.

### Step 9 — optimize the pair thermodynamically

Release both mechanisms simultaneously and use the actual thermodynamic
objective.

### Step 10 — re-tune the surrounding machine

With mechanisms fixed, adapt relevant thermodynamic and exchanger dimensions.

### Step 11 — evaluate robustness separately

Manufacturing sensitivity and mechanical-margin robustness are not represented
by the primary cadence score.

### Step 12 — use fresh-island saturation only when needed

Run an explicitly seed-free search when the question is search-space coverage
rather than design exploitation.

---

## 19. Generalization

The long-term objective is not to restart a completely unconstrained mechanism
search for every operating point.

Known families should serve as reusable initial structures.

For a nearby target:

1. begin from several retained families;
2. locally adapt them;
3. compare their complete-mechanism behavior;
4. use the thermodynamic solver for the final ranking;
5. add fresh broad exploration when existing families cease to cover the useful
   design space.

This applies across changes in:

- source temperatures;
- working fluid;
- machine scale;
- motor versus refrigeration operation;
- operating frequency;
- exchanger configuration.

The difficult piston and the winning mechanism family may change.

---

## 20. Simpler mechanism classes remain valid competitors

Six-bar mechanisms are not always necessary.

For applications where simplicity dominates, the same thermodynamic problem
should also be tested with simpler physical mechanism families.

A final engineering comparison may include:

- thermodynamic performance;
- mechanical efficiency;
- part count;
- joints and bearings;
- tolerance sensitivity;
- packaging;
- fabrication difficulty;
- cost;
- maintenance.

Mechanism complexity must earn its place through system-level benefit.

---

## 21. Software boundary

The reusable software abstraction belongs in:

```text
src/dada_solver/research/synthesis.py
```

It records:

- ordered synthesis stages;
- released coordinate groups;
- explicit separation between design exploitation and fresh-island saturation;
- retained family identities;
- mechanical constraints;
- position as the primary synthesis reference;
- acceleration as diagnostic only;
- mirroring as initialization only;
- the thermodynamic objective after paired release.

Fresh-island saturation retains declarative stage and ownership contracts.
The five geometric six-bar stages, and direct slider-crank and four-bar discovery
and local polish, execute through the shared geometric search engine described above. Thermodynamic evaluation remains
an independent Research operation.

---

## 22. What should not be inferred

This method does not establish that:

- one primary proxy is universally optimal;
- only a fixed number of useful primary families exists;
- optimizer capture frequency equals geometric mechanism-space volume;
- the closest kinematic fit gives the best thermodynamic machine;
- a mirrored pair should remain symmetric;
- six-bar mechanisms are always preferable to simpler linkages;
- one target cadence is universal.

The durable result is the hierarchical search architecture and the distinction
between discovery proxies, complete-mechanism quality and final
thermodynamic assessment.
