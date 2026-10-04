# Research technical reference

## Purpose and scope

This is the current contract for study schema, scientific identity, parameter ownership,
policies, search scheduling, persistence and compatibility. Production loaders and
campaign APIs are authoritative. For everyday commands, start with the [Research
guide](DADA_ENGINE_RESEARCH.md).

This reference does not define thermodynamic equations, mechanism families or report
layout. Their canonical references are listed at the end.

## Supported schemas and presets

| Schema | Role today |
|---|---|
| Schema 1 | Historical `fixed_pair_thermo5d` regression protocol; loadable and evaluable, resumable when runtime-compatible |
| Schema 2 | Generic `machine_design` studies with configurable kinematics, machine parameters, policies and wall/reservoir models |
| Schema 3 | `machine_design` extension for external-stream wall models and compiled-fluid support |

Schemas 2 and 3 share the fixed/active parameter and campaign layer. Schema 3 extends
the machine representation, not the campaign semantics. The study and portable basis
must have matching schema versions. Schema 1 is compatibility support, not the
recommended starting point for new research.

Current generic initialization entry points are `init kinematics`, `init
external-stream-refrigeration` and `init external-stream-motor`. The named regression
presets are described under legacy compatibility below.

## Schema 2/3 study anatomy

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

## Objectives and constraints

The objective contains `type` and `unit`:

| Type | Unit | Operating direction |
|---|---|---|
| `maximize_thermal_efficiency` | `1` | Motor |
| `maximize_cooling_cop` | `1` | Refrigeration |
| `maximize_motor_power` | `W` | Motor |
| `maximize_cooling_power` | `W` | Refrigeration |

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
constraint. Historical power, pressure, temperature and absolute-flow guards are not
universal defaults.

Model applicability/domain, study requirements, search bounds and numerical settings
have distinct owners. A stricter study Mach limit does not redefine the exchanger's
correlation domain. See the [limit-ownership audit](RESEARCH_LIMIT_OWNERSHIP_AUDIT.md).

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

Schema 2 base policies are explicit declarations:

| Key | Supported convention |
|---|---|
| `volume_partition` | `total_swept_and_clearance_ratios` |
| `mechanical_losses` | `unknown` |
| `useful_power` | `unavailable` |
| `local_reflux` | `retain_signed_flows` |
| `external_air_aerodynamic_losses` | `excluded_from_balance` |
| `fan_consumption` | `excluded_from_balance` |

`charge` accepts `explicit_inventory` or
`reference_pressure_at_maximum_total_volume_v1`. `outlet_valve_cda` accepts
`source_cda_times_count_ratio_v1`, `fixed_source_cda` or `geometry_conduit_area_v1`.
These select existing derivation conventions; the [microtube
reference](MICROTUBE_GAS_MODEL.md) owns their geometry rules and compatibility
requirements.

Schema 3 replaces the two external-air/fan keys with:

| Key | Convention |
|---|---|
| `external_loop_hydraulics` | `unmodelled` |
| `external_pump_fan_consumption` | `excluded_from_balance` |

An external fluid label does not provide a pump or hydraulic model. Excluded losses are
not asserted to be physically zero. `POLICIES` and `POLICIES_V3` in `research.schema_v2`
define the supported base policy sets.

## Warm starts and domain retry

Schema 2/3 `warm_start.initial_source` accepts `uniform` or `source_exact`.
`source_exact` requires a stored compatible source state and its original fixed
inventory. It is incompatible with geometry-derived reference-pressure charge. A warm
state is an initial guess, not proof of periodic convergence.

`numerical.domain_error_retry` accepts `none` or
`once_safe_uniform_state_with_source_wall_temperatures`. A safe retry changes
initialization, not the physical model or validity limits. Earlier rejected trial
evidence must remain distinct from final periodic-cycle constraint evidence;
presentation belongs in the [cockpit](DADA_ENGINE_RESEARCH_COCKPIT.md).

## Scientific identity

### `study_id`

This is the hash of canonical scientific study content. For schema 2/3 it includes basis
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

For schema 1 only, `--reference` selects the legacy regression point. Schema 2/3 uses
configured fixed/initial values and rejects `--reference`.

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
identity; old studies without it retain their scheduling behaviour.

## Local refinement and explicit initial evaluations

`research refine SOURCE --candidate ID --radius 0.20 --output outputs/local/study.toml`
accepts schema-2/3 sources with an active parameter space. Multiple sources must have
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

`research rescale` accepts schema-2/3 stored candidates and creates a new portable
study/basis, with parent/provenance and preserved constraints. It never integrates,
changes source histories or silently relaxes constraints. Local-region sources are
rejected: rescale the global source first, then refine. Detailed extensive/intensive
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

New campaigns use `history.jsonl.gz`, `recovery.json` and `state.json`. The compressed
journal is lossless and append-only: one complete canonical JSON record per independent
gzip member, with deterministic gzip timestamp. `recovery.json` is a temporary durable
completion slot, not a second history.

Legacy `history.jsonl` and `candidates/*.json` remain readable/recoverable. An existing
plain journal continues in plain format; simultaneous plain and compressed journals are
rejected as ambiguous. No automatic history migration or deletion of legacy result files
occurs.

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

## Legacy compatibility

Schema 1 remains the loadable `fixed_pair_thermo5d` protocol. `init sixbar-thermo5d`
creates its regression input, and `evaluate --reference` is schema-1-only. Old campaigns
remain inspectable; execution resume still requires strict runtime/source matching.

`init structured-c2-3952` is a named historical schema-2 regression seed/preset. Neither
legacy preset defines generic defaults for new studies. `refine` and `rescale` require
schema 2 or 3.

Historical extraction details belong in the [validation
ledger](history/RESEARCH_VALIDATION_LEDGER.md), [migration
matrix](DADA_ENGINE_RESEARCH_MIGRATION_MATRIX.md) and [implementation
plan](DADA_ENGINE_RESEARCH_IMPLEMENTATION_PLAN.md).

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
| Historical Research extraction evidence | [Historical ledger](history/RESEARCH_VALIDATION_LEDGER.md) |
| Limit ownership history | [Ownership audit](RESEARCH_LIMIT_OWNERSHIP_AUDIT.md) |

For retained geometry and synthesis methodology, use the [primary
catalogue](PRIMARY_FOUR_BAR_FAMILIES.md), [six-bar
catalogue](SIX_BAR_MECHANISM_FAMILIES.md) and [synthesis
method](MECHANISM_SYNTHESIS_SEARCH.md).