# Research technical reference

## Purpose and scope

This is the current contract for study schema, scientific identity, parameter ownership,
policies, search scheduling, persistence and compatibility. Production loaders and
campaign APIs are authoritative. For everyday commands, start with the [Research
guide](DADA_ENGINE_RESEARCH.md).

This reference does not define thermodynamic equations, mechanism families or report
layout. Their canonical references are listed at the end.

## Supported schema and presets

Study and machine basis files require `schema_version = 3`. The loader accepts
only this schema and performs no automatic schema conversion.
Generic initialization entry points are `init kinematics`,
`init external-stream-refrigeration` and `init external-stream-motor`.
They share the fixed/active parameter and campaign layer, with explicit physical
families and matching study/basis versions.

## Study anatomy

| Group | Ownership |
|---|---|
| `schema_version` | Selects the accepted study/basis generation |
| `study` | Nonempty `name`, `protocol = "machine_design"`, `purpose`; optional `parent_candidate_id` |
| `sources` | Portable machine basis at `sources.machine.path`, verified by `sha256` |
| `kinematics` | Family, coupling and settings for SMALL/LARGE |
| `parameters` | Explicit fixed/active coordinate declarations |
| `policies` | Physical and derivation conventions |
| `objective` | Ranking objective and unit |
| `constraints` | Thermodynamic/design requirements |
| `mechanical_constraints` | Side-scoped mechanical requirements |
| `screening` | Mechanical sampling and optional declared volume screen |
| `search` | Global or local Sobol definition |
| `numerical` | Candidate evaluation settings |
| `warm_start` | Initial-state policy |
| `execution` | Default/per-invocation scheduling controls |
| `charge_reference` | Required only with reference-pressure inventory |

The display name and source filesystem location are not scientific physics. Verified
basis content is. Machine coordinates omitted from explicit declarations retain basis
defaults; family coordinates must be supplied by declarations or mechanism artifacts.
Unknown keys and invalid ownership are rejected.

`numerical` contains `maximum_cycles`, `backend`, `candidate_budget_seconds` and
`domain_error_retry`; other integration inputs come from the basis. `execution` contains
`default_budget`, `default_max_candidates`, `initial_evaluation_seconds` and
`deadline_grace_seconds`. Execution controls do not define scientific study identity.
Invocation budgets and candidate caps control how much work is attempted, not the
physical acceptance criteria.

## Parameter ownership and parameter types

One declaration owns a coordinate. Active coordinates form one vector in declaration
order; fixed and active ownership cannot overlap. Names, units, family ownership, types,
bounds and initial values are validated.

A fixed input:

```toml
[[parameters]]
name = "operation.frequency_hz"
unit = "Hz"
value = 2.8
```

A `continuous` input:

```toml
[[parameters]]
name = "operation.frequency_hz"
unit = "Hz"
kind = "continuous"
initial = 2.8
lower = 1.0
upper = 4.0
transform = "linear"
```

`ContinuousParameter` supports `linear` and `log`. Bounds must be finite, ordered and
contain the initial value; logarithmic bounds must be positive. Normalized coordinates
lie in `[0, 1]`.

`IntegerParameter` uses `kind = "integer"` and `encoding = "nearest_even_v1"` instead of
`transform`. Bounds are inclusive, positive, ordered true integers, as is `initial`;
floating counts and booleans are rejected. Decoding rounds to nearest even, with the
supported range bounded by exact integer representability in the encoding.

`ChoiceParameter` uses `kind = "choice"`, a nonempty ordered `choices` list of distinct
finite scalar values, and an `initial` belonging to that list. Choices occupy equal
normalized bins. Their order is an encoding, not physical distance. Candidate
construction canonicalizes a choice to its selected bin centre.

Family-defined fixed categories, including mechanism branches, must remain fixed.
Declaring a number or a `choice` does not bypass a family's ownership rules.

## Study editor

`study edit` is a terminal assistant for an existing, valid schema-3 TOML.
Install the optional UI with `pip install "dada-engine-solver[tui]"`; it uses
`prompt_toolkit` and is imported only for interactive editing. The solver,
batch operations and CLI help do not require this extra. Redirected stdin/stdout
are refused rather than entering a blocking questionnaire.

```console
dada-research study edit path/to/retuned/study.toml
dada-research study release path/to/retuned/study.toml --group mechanisms \
  --radius 0.05 --output path/to/combined/study.toml
dada-research study freeze path/to/study.toml --group frequency \
  --output path/to/fixed-frequency/study.toml
dada-research study release path/to/study.toml \
  --parameter operation.frequency_hz --bounds operation.frequency_hz=1:4:log \
  --search global --output path/to/bounded/study.toml
```

The menu uses arrows to navigate, Space to toggle activity, Enter to open a
group/parameter and Escape to go back or cancel. `[~]` indicates partial activity;
`[-]` marks a fixed scientific category. Group activation preserves already active
hardware. Individual edits support exact values, numeric bounds, linear/log
transforms, integer counts and declared choices; displayed values are never
rounded before saving. Search, execution, review and cancellation have separate
pages. Review includes effective bounds for every local region, value changes,
artifact hash changes, source identity and warnings. Nothing is written on cancel.

`research.study_editor.StudyEditor` implements these operations without terminal
imports or thermodynamic integration. Inventory includes machine-basis values,
explicit fixed/active declarations and coordinates supplied by artifacts.
`mechanisms` selects continuous physical coordinates on both sides;
`kinematics-small` and `kinematics-large` also support abstract families.
Hardware groups are `exchangers`, `volumes`, `frequency`, `charge`, `valves` and
`external-stream`; they filter existing coordinates only. Branches, orientations,
shared crank direction, family/settings and raw spline controls remain fixed by
their current contracts. Shared-crank phase is one coordinate, not two independent
phases. Charge remains governed by its current policy; reference-inventory charge
cannot be made an independent mass parameter by the editor.

Automatic domain resolution uses the current declared domain, then recorded
original domains (including `paired_thermodynamic.reference_boxes`), then the
centered, positive-length-aware synthesis reference-box policy for physical
mechanisms. Angles stay unwrapped. A box which no longer contains the exact
current mechanical value can fall back to that centered synthesis policy, never
to clamping. Other incompatible recorded domains require explicit resolution.
These are **search bounds**, not physical validity domains. Missing domains
require interactive bounds or batch `--bounds NAME=LOW:HIGH[:linear|log]`;
`--choices 'NAME=["first", "second"]'` supplies categorical domains. Explicit
user domains override automatic resolution and undergo the same Research checks.

**Declared domain and local radius are different controls.** `--search global`
uses `fixed_global_bounds` and the full declarations; `--radius` is rejected in
that mode. Release defaults to local Sobol with radius 0.05. With
`local_regions_v1`, `radius_fraction` restricts **all active parameters**, including
previously active hardware, to `[max(0,c-r), min(1,c+r)]` in their declared
normalized coordinates. Log transforms and integer decoding use the existing
parameter APIs. `choice_scope = "declared_choices"` retains all declared choices;
category order is not a distance. The assistant shows effective physical bounds.

Recalibrate a selected group's **declared bounds** around its exact current
initials, independently of the scheduler radius:

```console
dada-research study bounds recenter path/to/active/study.toml \
  --group mechanisms --half-width 0.01 --output path/to/narrow/study.toml
```

`StudyEditor.recenter_bounds(groups=..., parameters=..., half_width=...)` uses
`encode(initial) ± half_width` in the **current** declared normalized domain,
clips to `[0, 1]`, and decodes with the existing linear/log or integer encoding.
It changes only `lower/upper`, never initials, types, units or transforms, and
never expands the previous domain. Fixed parameters and categorical choices are
reported as skipped. An integer interval collapsing to one value is rejected;
all selected rows are checked before any change is applied. Unselected exchanger
bounds remain unchanged. This sets search windows, not equal thermodynamic
sensitivities or new physical validity domains.

Each group submenu offers **Recalibrate bounds for this group**, with the old
reference intervals, proposed intervals and Apply/Back. Repeating the operation
uses the newly declared intervals and therefore narrows cumulatively; the preview
says so explicitly. **Undo last bounds recalibration** restores the preceding
intervals in the session, provided subsequent parameter edits would not be lost.
Nothing is saved until Save. Applied operations are recorded in `changed_bounds`
and `study_edit.bounds_recalibrations`, including their reference declarations.

The scheduler radius is unchanged by recalibration. The TUI separately asks
whether to set a local radius to 1.0; accepting affects **all active parameters**,
including unmodified hardware. Batch `--radius 1` requests that same separate
change explicitly. With radius 1, Sobol explores the entire newly declared
numeric domain. **Show effective bounds** groups declared bounds, exact initials,
effective bounds and clipping indicators, and displays the global radius;
choices retain `choice_scope = "declared_choices"` without numeric distances.
Group menus remain open after toggles and individual edits, preserving parameter
focus and scroll position. Back or Escape returns to the main menu.

An altered active set or changed initials creates one region at the exact current
initials; unchanged valid regions are preserved. Replacing multiple regions
requires explicit interactive confirmation or `--recenter` in batch mode.
The center is unevaluated: its provenance uses the real computed identity of the
source initial configuration, not an invented evaluated candidate. An entirely
fixed study uses global settings and is intended for `evaluate`, not an empty
search. Freezing always uses the current `initial`, never an assumed campaign best.

Save creates a new portable TOML, basis and referenced artifacts. Active overrides
can keep their reference artifacts; fixed effective geometry differing from the
reference receives a new artifact, preserving settings, embedded constraints and
parent provenance. Unchanged artifacts retain their data and hashes. Existing
objectives, constraints, scientific policies, warm-start data and prior refit/
adaptation/retuning provenance are preserved. `basis.provenance.study_edit` records
edits, domain sources, search changes and unevaluated-center evidence. All staged
inputs are loaded, compiled and production-preflighted before publication; write
failures roll back only newly created files. Source files, existing destinations
and running campaign directories are never overwritten. No integration, automatic
geometry repair or automatic campaign launch occurs.

## Objectives and constraints

The objective contains `type` and `unit`:

| Type | Unit | Operating direction |
|---|---|---|
| `maximize_thermal_efficiency` | `1` | Motor |
| `maximize_cooling_cop` | `1` | Refrigeration |
| `maximize_motor_power` | `W` | Motor |
| `maximize_cooling_power` | `W` | Refrigeration |
| `maximize_cooling_power_per_total_microtube` | `W/microtube` | Refrigeration |
| `maximize_cooling_cop_times_power_per_total_microtube` | `W/microtube` | Refrigeration |

`maximize_cooling_power_per_total_microtube` maximizes
`cooling_power_w / (microtube.heat_in.tube_count + microtube.heat_out.tube_count)`.
The denominator counts all physical microtubes in both exchangers, not hydraulic
half-links. The corresponding metric is `cooling_power_per_total_microtube_w`.
This manufacturing productivity metric imposes no minimum COP or cooling power
by itself.

`maximize_cooling_cop_times_power_per_total_microtube` maximizes
`cooling_cop * cooling_power_per_total_microtube_w`. This objective jointly
rewards cooling efficiency and cooling productivity per physical microtube. It
does not impose minimum COP or cooling power by itself. Its reported metric is
`cooling_cop_times_power_per_total_microtube_w`, in `W/microtube`.

The objective must agree with the basis operating direction. Motor power and mechanical
input are indicated quantities, not useful shaft power. Mechanical losses remain unknown
unless a separate physical model supplies them.

`PHYSICAL_CONSTRAINTS` in `research.margins` defines the accepted vocabulary:

| Constraint `type` | Numeric field | Unit |
|---|---|---|
| `minimum_motor_power` | `required_power` | `W` |
| `minimum_cooling_power` | `required_power` | `W` |
| `maximum_mechanical_input_power` | `limit` | `W` |
| `maximum_pressure` | `limit` | `Pa` |
| `maximum_temperature` | `limit` | `K` |
| `maximum_absolute_mass_flow` | `limit` | `kg/s` |
| `maximum_mach_number` | `limit` | `1` |
| `valid_thermodynamic_model` | None | `1` |
| `periodic_convergence` | None | `1` |

Each declaration includes `type` and `unit`; numeric fields must be positive and finite.
Duplicates are rejected. Generic studies require an explicit `valid_thermodynamic_model`
constraint.

Model applicability/domain, study requirements, search bounds and numerical settings
have distinct owners. A stricter study Mach limit does not redefine the exchanger's
correlation domain. See the [current ownership rule](PHYSICS_DECISIONS.md).

Mechanical requirements use separate side-scoped declarations:

```toml
[[mechanical_constraints]]
side = "small"
metric = "minimum_primary_transmission_sine"
relation = "minimum"
limit = 0.30
unit = "1"
```

Metric availability, units and relations are family-specific; see the [kinematics
reference](DADA_ENGINE_RESEARCH_KINEMATICS.md). Screening normally uses 1440 samples,
with a minimum of 360. `maximum_large_enclosed_volume_m3` is an optional explicit
screening ceiling. Sampled screens are not continuous proofs or manufacturing
certification; some individual metrics have stronger analytic implementations.

## Policies

Policies are explicit declarations:

| Key | Supported convention |
|---|---|
| `volume_partition` | `total_swept_and_clearance_ratios` |
| `mechanical_losses` | `unknown` |
| `useful_power` | `unavailable` |
| `local_reflux` | `retain_signed_flows` |

`charge` accepts `explicit_inventory` or
`reference_pressure_at_maximum_total_volume_v1`. `outlet_valve_cda` accepts
`source_cda_times_count_ratio_v1`, `fixed_source_cda` or `geometry_conduit_area_v1`.
These select existing derivation conventions; the [microtube
reference](MICROTUBE_GAS_MODEL.md) owns their geometry rules and compatibility
requirements.

The external boundary requires:

| Key | Convention |
|---|---|
| `external_loop_hydraulics` | `unmodelled` |
| `external_pump_fan_consumption` | `excluded_from_balance` |

An external fluid label does not provide a pump or hydraulic model. Excluded losses are
not asserted to be physically zero. `POLICIES_V3` in `research.study_schema`
define the supported base policy sets.

## Warm starts and domain retry

Schema 3 `warm_start.initial_source` accepts `uniform` or `source_exact`.
`source_exact` requires a stored compatible source state and its original
inventory at the study initial point. Explicit inventory can be active: candidate
guesses at other inventories use the existing cached-state rescaling, preserving
specific gas energies; original-inventory gas states remain exact. Wall energies
follow current capacities. Integrated states are never renormalized.
It is incompatible with geometry-derived reference-pressure charge. A warm
state is an initial guess, not proof of periodic convergence.

`numerical.domain_error_retry` accepts `none` or
`once_safe_uniform_state_with_source_wall_temperatures`. A safe retry changes
initialization, not the physical model or validity limits. Earlier rejected trial
evidence must remain distinct from final periodic-cycle constraint evidence;
presentation belongs in the [cockpit](DADA_ENGINE_RESEARCH_COCKPIT.md).

`dada-research mechanism retune` generates a hardware-only local study from a
paired thermodynamic candidate or retuning descendant. Original non-kinematic
active domains are restored from portable paired provenance, with exact current
initials and fixed adapted mechanism artifacts. `--scope source-active` restores
the full eligible set; `--group`/`--parameter` select a subset. Numeric radii use
the original normalized domains; `search.choice_scope = "declared_choices"`
keeps original categories eligible without interpreting their order as distance.
See [hardware retuning](DADA_ENGINE_RESEARCH_KINEMATICS.md#hardware-retuning-with-fixed-mechanisms)
for groups, artifact linkage, provenance and generation without integration.

## Scientific identity

### `study_id`

This is the hash of canonical scientific study content. For schema 3 it includes basis
content, fixed/active ownership, parameter bounds/transforms/ encodings, kinematic and
mechanism scientific content, policies, objective, constraints, mechanical constraints,
screening, search definition, numerical settings and warm-start policy. Reference-charge
settings are included too.

It excludes `study.name`, the machine source path, mechanism artifact path locations and
`[execution]`. Other retained study metadata, including purpose and parent lineage when
present, remain in the canonical content. Basis hashes are verified; moving a portable
study with unchanged contents does not redefine the science. Re-serializing a basis can
change its recorded hash, so portability means preserving verified contents, not
arbitrary reformatting.

Changing the search definition changes study identity, even when the underlying
thermodynamic machine is unchanged.

### `definition_id`

The executable definition adds current runtime identity to the scientific study.
`runtime_identity()` includes the production Python-source digest and Python, NumPy and
SciPy versions. Backend identity is also recorded. Schema 3 adds fluid/kernel execution
identity where applicable.

A source/runtime change can make execution resume incompatible even when the scientific
study is unchanged. Documentation-only edits do not alter the production source digest.

### `candidate_id`

The canonical candidate payload includes `definition_id`, normalized active coordinates,
physical active values, selected families and numerical settings. Fixed inputs are
represented through the definition. Choice coordinates are canonicalized. Exact
configured/centre values are retained by physical-point construction rather than
replaced by a floating round trip.

Local-region origin is provenance outside the candidate payload. Two regions can
therefore produce the same candidate ID within a definition. This does not create a
cache across different studies/definitions. Historical IDs must not be reused as current
cache identities. Cache equality is exact payload equality, not approximate geometric or
numerical similarity.

## Validation and exact evaluation

`research validate STUDY` parses declarations, verifies hashes, constructs production
geometry/machine objects and performs preflight checks. It does not integrate
thermodynamics and does not establish cycle feasibility.

`research evaluate STUDY --output RESULT` evaluates one exact configured point, using
active `initial` values unless overridden with `--set NAME=VALUE`. Overrides must
address declared active coordinates within their domains. It does not run a search and
never overwrites an existing result.


## Global Sobol search

```toml
[search]
type = "sobol"
domain = "fixed_global_bounds"
seed = 29092026
scramble = true
```

Optional `evaluate_initial = true` schedules the exact initial physical point first. If
absent or false, global search starts at Sobol index zero without inserting that point.
An unevaluated initial value is not an incumbent. This option participates in study
identity.

## Local refinement and explicit initial evaluations

`research refine SOURCE --candidate ID --radius 0.20 --output outputs/local/study.toml`
accepts schema-3 sources with an active parameter space. Multiple sources must have
the same `study_id`. Campaign sources require selectors; standalone single-evaluation
artifacts may supply their sole record implicitly.

Refinement creates a new portable study, copying the basis and required mechanism
artifacts without changing thermodynamic physics or global bounds. It never overwrites
destination files and does not import source evaluations as already completed results.
The new search definition makes this a new study.

The local search fields are `domain = "local_regions_v1"`, `radius_fraction` in `(0,
1]`, `allocation = "round_robin"` and `evaluate_centers = true`. Each region embeds
`id`, `source_candidate_id`, `source_study_id` and `center`. Its centre contains exactly
the active parameters in physical units; source IDs are provenance, not external file
dependencies.

For a globally normalized centre coordinate `z` and radius `r`, the local interval is
`[max(0, z-r), min(1, z+r)]`. Thus `r = 0.20` means ±20% of the global normalized search
width, not ±20% of the physical value. Log transforms retain their normal
interpretation; clipping near a global bound creates asymmetry.

Exact distinct centres are scheduled before any Sobol draw. Duplicate centres are
evaluated once and credited to all associated regions. Each region then has its own
Sobol index, with deterministic round-robin allocation and the same seed/scramble
policy. Normal integer and choice encodings remain in force; scheduled physical points
are canonically re-encoded. Region provenance is recorded separately as `search_origin`.

Resume preserves centres, region ordering, indices and pending work. There is no
adaptive basin selection, local gradient optimizer or cross-study cache.

## Capacity rescale

`research rescale` accepts schema-3 stored candidates and creates a new portable
study/basis, with parent/provenance and preserved constraints. It never integrates,
changes source histories or silently relaxes constraints. Stored global or local
candidates supply their physical coordinates; local search basins are replaced by
a global Sobol domain in the standalone output study. Detailed transformation
rules belong in [capacity scaling](DADA_ENGINE_RESEARCH_CAPACITY.md).

## Filling at a reference pressure and maximum total gas volume

The production policy is `reference_pressure_at_maximum_total_volume_v1`:

```toml
[policies]
charge = "reference_pressure_at_maximum_total_volume_v1"

[charge_reference]
pressure_pa = 100000.0
temperature_k = 293.15
volume_state = "maximum_total_gas_volume"

[warm_start]
initial_source = "uniform"
```

Under this policy `charge.total_mass_kg` must not be declared. Pressure is absolute;
pressure and temperature must be positive and finite. The current implementation
requires `CaloricallyPerfectGas`; it does not silently apply an ideal filling
approximation to a real-gas table.

Both cylinders are evaluated at the same solver angle. The final connected production
volumes include actual exchanger hold-up without summing independent piston maxima or
adding exchanger volumes twice. Inventory is `m = p Vmax / (R T)`. `source_exact` is
rejected; ordinary compatible warm guesses may be inventory-rescaled. Reference settings
participate in scientific identity.

The method identifier is
`periodic_total_volume_grid8192_breakpoints_derivative_roots_local_refinement_v1`. It
uses 8192 periodic intervals, declared kinematic breakpoints, derivative-root bracketing
and local refinement of sampled peaks. Near-equal maxima use a deterministic first-angle
tie rule; the recorded volume belongs to that angle. This is a deterministic numerical
maximum search, not an analytic global proof for arbitrary pathological motion laws. It
is independent of `screening.samples`.

`derived.charge` records policy, reference pressure/temperature, total gas volume,
reference angle, derived mass and method. The normal total-mass observable reflects that
derived inventory.

## Categorical parameters

For a supported machine coordinate:

```toml
[[parameters]]
name = "valve.heat_in.placement"
unit = "1"
kind = "choice"
initial = "downstream"
choices = ["downstream", "upstream"]
```

Bins are equal-width in declared order; encoding returns the selected bin centre.
Categorical Sobol draws canonicalize to that category. Local radius operates on these
bins, not on a physical category distance: a binary centre with radius 0.20 stays in its
own bin; a sufficiently large radius can reach another bin. Reports exclude categories
from numerical bound-pressure interpretation. Mechanism assembly branches remain fixed
unless their family interface explicitly supports another ownership contract.

## Run, resume and persistence

`run` requires at least one active coordinate, creates a campaign directory and refuses
a nonempty destination. The default directory is `campaign` beside the study TOML;
`--directory` overrides it. `resume` accepts a campaign directory or a TOML path whose
sibling `campaign` is selected. A missing campaign never silently starts a new search.

Research snapshots `definition.json`, `study.toml`, `basis.json` and `study.json`, plus
required mechanism artifacts. Edits to the original external TOML do not mutate that
stored campaign definition.

Campaigns use `history.jsonl.gz`, `recovery.json` and `state.json`. The compressed
journal is lossless and append-only: one complete canonical JSON record per independent
gzip member, with deterministic gzip timestamp. `recovery.json` is a temporary durable
completion slot, not a second history.

Campaign persistence requires `history.jsonl.gz` and the single `recovery.json`
completion slot. Plain journals and per-candidate recovery layouts are rejected.

After evaluation, completion is written atomically to recovery and synchronized, then
appended/synchronized to the journal; state is published durably before recovery is
removed. See the [cockpit contract](DADA_ENGINE_RESEARCH_COCKPIT.md) for crash-order
details and presentation.

## Cache, recovery and compatibility

Exact candidate identity drives cache reuse, with cache hits persisted as attempts. It
is not approximate physical-value deduplication. Pending candidates are saved before
evaluation, allowing interruption or power-loss recovery. `--retry-incomplete`
explicitly retries unresolved deadline-limited work; scheduled searches also preserve
interrupted pending work for continuation.

An incomplete final journal tail is recoverable under the journal contract. Non-final
corruption is an error; compressed CRC/decode corruption and invalid JSON in a complete
gzip member are also errors, not silently skipped records.

Execution resume reconstructs the definition from stored snapshots and requires the
recorded `definition_id`. Source/runtime incompatibility prevents execution resume.
Offline inspection can remain available when execution is incompatible; inspection does
not authorize continuing under a silently changed definition.

## Inspection and reporting boundary

`status`, `report` and `compare` inspect stored campaigns/evaluations without changing
scientific identity. Requested derived curves may be reconstructed using current
production code; thermodynamic curves can require an explicit report replay. Runtime
mismatch remains visible. Cross-study reports do not invent a combined ranking for
incompatible studies.

Command use belongs in the [guide](DADA_ENGINE_RESEARCH.md); report controls and
evidence presentation belong in the [cockpit](DADA_ENGINE_RESEARCH_COCKPIT.md).

## Related references

| Subject | Canonical documentation |
|---|---|
| Everyday commands | [Research guide](DADA_ENGINE_RESEARCH.md) |
| Study schema, identity, policies, search and compatibility | This file |
| Kinematic families, mechanisms and mechanical metrics | [Kinematics reference](DADA_ENGINE_RESEARCH_KINEMATICS.md) |
| Reports, persistence UX and evidence presentation | [Cockpit](DADA_ENGINE_RESEARCH_COCKPIT.md) |
| Capacity transformation | [Capacity scaling](DADA_ENGINE_RESEARCH_CAPACITY.md) |
| Microtube physics | [Microtube model](MICROTUBE_GAS_MODEL.md) |
| External streams | [External-stream model](EXTERNAL_STREAM_THERMAL_MODEL.md) |
| Working fluids | [Fluid models](WORKING_FLUID_MODELS.md) |
| Current validation status | [Validation](validation.md) |

For retained geometry and synthesis methodology, use the [primary
catalogue](PRIMARY_FOUR_BAR_FAMILIES.md), [six-bar
catalogue](SIX_BAR_MECHANISM_FAMILIES.md) and [synthesis
method](MECHANISM_SYNTHESIS_SEARCH.md).
