# Dada-Engine Research

Research V3 adds declared external thermal streams, motor/refrigerator wall
cycles and an ideal-generated compiled fluid-table validation path. It reuses
the V2 campaign and mechanism layer.

Research V2 selects kinematic families independently for the small and large
cylinders, with fixed or active coordinates in the existing campaign engine.
The V1 fixed six-bar / five-parameter study remains available as a regression
preset. Exact identities, Sobol continuation, deadlines and recovery are retained;
production thermodynamic physics is unchanged.

## Limit ownership and migration

Generic V2/V3 presets impose no inherited 25 W minimum power, 1.2 MPa pressure,
850 K temperature, 0.08 kg/s internal-flow or 66 L volume ceiling. These are
optional explicit study requirements, not universal physical limits. The named
V1 historical-parity preset keeps its original requirements. Model domains,
study requirements, search bounds and numerical settings have distinct owners;
see [the audit and migration notes](RESEARCH_LIMIT_OWNERSHIP_AUDIT.md).

Pressure equalization is diagnostic-only; isothermality fields/criteria are
removed. Old validity keys load but are ignored. Microtube Mach validity comes
from each exchanger's gas model, while `maximum_mach_number` remains an optional
study constraint. New mechanism seed quality limits live in editable
`[[mechanical_constraints]]`, outside seed artifacts. Existing embedded artifact
constraints remain enforced.

Existing results are not reclassified or overwritten. After this implementation
change, regenerate the intended study and use a new campaign directory; runtime
compatibility prevents silently resuming the old scientific computation.
`refine`/`rescale` preserve explicit source constraints, including historical caps.

## Daily cockpit workflow

```sh
dada-research report outputs/my_study/campaign
dada-research resume outputs/my_study/campaign --budget 30m
```

Reports default to `CAMPAIGN/report.html` and may overwrite derived HTML.
Terminal output is compact; use `--list-candidates` or `--json` explicitly.
Generated/rescaled studies default to 512 new attempts per invocation, with
`--max-candidates` available for a deliberate override.

On an interactive terminal, each completed evaluation gets exactly one permanent
line: its campaign evaluation number (one-based), region and center/Sobol origin,
status, and either principal performance metrics or rejection evidence. `BEST`
marks an improved feasible champion; `(cached)` identifies exact cache reuse.
A separate bottom line refreshes in place about every ten seconds during solver
progress checks and immediately after each completed evaluation. It shows phase
counts, elapsed time/budget and the champion, without repeated error counters.
Its width follows the output terminal and preserves progress, feasibility and
COP/objective ahead of secondary metrics. `Pin` is indicated mechanical input;
`Pgas` is indicated gas power, never an assumed useful shaft output.

Redirected/non-TTY output contains the start line, one line per completed
evaluation and a final `Finished` summary. It contains no periodic status spam,
carriage returns or terminal escape sequences. Normal completion, Ctrl-C and
exceptions terminate/clear the interactive status before returning control.
Detailed evidence stays in the journal; presentation neither recalculates reports
nor changes scheduling, candidate identity construction or durable resume.
The existing source/runtime compatibility checks still apply after code updates.

The HTML shows the
campaign funnel, bound pressure, factual review actions and copyable commands
that follow the selected candidates. Diagnostic suggestions filter the exact
failure category or violated constraint and open matching stored records; bound
actions show bound evidence. Only report regeneration gets a regeneration command.
The full candidate data remains in a scrollable table (about 15 visible rows),
with a Basin column and sticky sortable headers. Comparison opens on the best
two candidates under the existing objective ranking; selectors remain editable.
Constraint evidence prioritizes violated/near-boundary limits with actual values,
signed/relative margins and candidate/passage context. Missing model boundaries
are not guessed. Native solver stderr is captured around
solver calls while original Python exceptions remain visible.

New evaluations use a single `recovery.json` followed by durable journal/state
commits; old per-candidate files remain readable and recoverable. See the
[cockpit and persistence reference](DADA_ENGINE_RESEARCH_COCKPIT.md) for exact
thresholds, crash ordering, compatibility and fd 2 scope.

## External streams and refrigeration (V3)

```sh
PYTHONPATH=src python3 -m dada_solver.research init external-stream-refrigeration --output outputs/my_cooling/study.toml
PYTHONPATH=src python3 -m dada_solver.research validate outputs/my_cooling/study.toml
PYTHONPATH=src python3 -m dada_solver.research evaluate outputs/my_cooling/study.toml --output outputs/my_cooling/reference.json --budget 2m
PYTHONPATH=src python3 -m dada_solver.research report outputs/my_cooling/reference.json --html outputs/my_cooling/reference.html
```

The preset is a bounded validation fixture, not an optimized cooling cell.
`external-stream-motor` selects the motor direction. Add
`--small structured_c2_15p --large structured_c2_15p` to configure the existing
structured family; changing family is a new study and does not guarantee a
valid refrigeration cycle. Individual stream, inventory, frequency, exchanger
and motion parameters use the same fixed/active declaration as V2. Nothing is
active by default. Choose `maximize_cooling_cop` [1] or
`maximize_cooling_power` [W], with an optional indicated-input bound.

See [external stream equations, units and configuration](EXTERNAL_STREAM_THERMAL_MODEL.md)
and [working-fluid interfaces, compiled tables and limitations](WORKING_FLUID_MODELS.md).
The external fluid label implies no liquid hydraulic correlation or pump power.
The fluid table proves compiled reconstruction using analytic ideal-gas data;
it is not a cryogenic helium model. Old V1/V2 inputs and offline inspection are
retained; source/runtime changes still require a fresh execution directory.

The [bounded candidate report](../outputs/research_v3/validated/candidate_comparison.html)
contains two Sobol candidates and a persisted resume. The
[cross-boundary report](../outputs/research_v3/validated/cross_boundary_report.html)
shows the refrigerator, its table replay and an air motor side by side, with
explicit differing-study warnings. It is not a fair optimization ranking.

## Select a motion family (V2)

```sh
PYTHONPATH=src python3 -m dada_solver.research init kinematics --small slider_crank --large harmonic --output outputs/my_motion/study.toml
PYTHONPATH=src python3 -m dada_solver.research validate outputs/my_motion/study.toml
PYTHONPATH=src python3 -m dada_solver.research evaluate outputs/my_motion/study.toml --output outputs/my_motion/reference.json --budget 3m
PYTHONPATH=src python3 -m dada_solver.research report outputs/my_motion/reference.json --html outputs/my_motion/reference.html
```

Templates start with fixed values. Replace a parameter's `value` with
`initial`, `lower`, `upper`, `kind = "continuous"`, and `transform = "linear"`
to activate it, retaining its name and unit. A study with no active coordinates
uses `evaluate`; a study with active coordinates can also `run` and `resume`.
Thermal, hardware and kinematic coordinates use the same declaration and vector.

Available families are harmonic, centered/offset slider-crank, four-bar, six-bar,
periodic free spline, Fourier C2, structured C2 15p, ideal piecewise, shared-origin
four-stage, independent four-stage laws, and the nine-coordinate `hybrid_compact` law. `--small` and `--large` accept mixed
families. `init structured-c2-3952` creates the original structured historical
candidate with its own thermal basis and warm state.

The [kinematics and mechanism reference](DADA_ENGINE_RESEARCH_KINEMATICS.md)
contains executable examples, exact names/units, artifact conventions, screening
constraints, and the boundary between available thermo-mechanical search and
future hierarchical synthesis. The [bounded demonstration](../outputs/research_kinematics_v2/comparison.html)
compares recorded family evaluations; it is not a fair optimization ranking.

Every reported constraint now exposes current value, limit, absolute margin,
relative margin where meaningful, and satisfied/violated/unavailable state.
Amber highlights satisfied limits within 5%; it is not a safety factor.
Research provides evidence without automatic scientific recommendations.

## Rescale and inspect an existing candidate

`rescale SOURCE --candidate ID --mode capacity --factor 5 --output study.toml`
creates a new portable V2/V3 study and basis centered on that candidate. It
scales inventory, cylinder capacity, parallel tubes, CdA and extensive thermal
inputs together, including basis-owned inputs. Active extensive bounds follow
the same factor; intensive inputs and all constraints stay unchanged. It never
integrates, overwrites an existing study, or modifies source histories.

Reports now provide clickable column sorting and a **Sort by** control for
metrics, active coordinates and constraint margins. Topology availability,
reflux detection and worst signed flow are visible in the main table.
`report SOURCE --candidate ID --plots volumes --html comparison.html` adds
both cylinder volume curves without thermodynamic integration. Repeat
`--candidate` to compare candidates; HTML remains standalone and offline.

See [capacity scaling, exact rules and examples](DADA_ENGINE_RESEARCH_CAPACITY.md)
for limitations, valve-event availability and the measured Human Cell ×5 checks.

## Local refinement and explicit initial evaluations

`run` starts a campaign from a study; `resume` continues its saved schedule.
`rescale` changes machine capacity and its extensive inputs. `refine` instead
creates a new portable study with unchanged physical inputs, objective,
constraints, global parameter bounds and transforms. It does not run an
optimization or import a source result as an already evaluated center.

```sh
dada-research refine outputs/source/campaign --candidate abc123 --radius 0.20 --output outputs/local/study.toml
dada-research validate outputs/local/study.toml
dada-research run outputs/local/study.toml --directory outputs/local/campaign --budget 3m --max-candidates 32
dada-research resume outputs/local/campaign --budget 3m --max-candidates 32
dada-research report outputs/local/campaign
```

Use a full ID, a unique prefix, or `best`. Repeat `--candidate` for several
basins in one source campaign. Several standalone `research_evaluation_v1`
artifacts can supply centers implicitly:

```sh
dada-research refine replay_A.json replay_B.json replay_C.json --radius 0.20 --output outputs/multi/study.toml
```

Sources must have the **same scientific study identity**; runtime/source-code
fingerprints may differ. This deliberately strict compatibility check refuses
mixing different physical definitions, bounds or policies. Candidate IDs and
source study IDs are embedded provenance, never external dependencies. The
new study includes a copied basis and any required mechanism artifacts.
Existing files are never overwritten. Creation defaults to 512 new attempts
per invocation; validation runs should explicitly use a small limit.

The generated V2/V3 TOML contains the following search structure (the real
center table contains every active coordinate in physical units):

```toml
[search]
type = "sobol"
domain = "local_regions_v1"
seed = 29092026
scramble = true
radius_fraction = 0.20
allocation = "round_robin"
evaluate_centers = true

[[search.regions]]
id = "basin_1"
source_candidate_id = "<full SHA-256 candidate ID>"
source_study_id = "<full SHA-256 study ID>"

[search.regions.center]
"volume.swept_ratio" = 1.5
# All other active physical values follow.
```

For global normalized coordinate `z`, the interval is
`[max(0, z-radius), min(1, z+radius)]`. Radius is finite and in `(0, 1]`;
0.20 means **±20% of the original global normalized width**, not ±20% of the
physical value. Log transforms keep their existing meaning. Near a global
bound the interval becomes asymmetric. Sobol points are mapped into this
interval and decoded using the existing parameter implementation. Integer
counts keep nearest-even decoding: rounding can move the encoded integer by
up to half a bin beyond the continuous local interval, never beyond the global
integer bounds. Exact center values are preserved without an encode/decode
round trip. Branches and unordered categories remain fixed scientific inputs;
they are not assigned a numeric distance or made active by refinement.

Every distinct center is evaluated before any local Sobol point. The first
center also becomes the standalone `evaluate` initial. Each region then has
its own Sobol index, using the same seed/scramble pattern. Allocation is fixed
round-robin in declared region order, with no adaptive basin selection. Equal
centers are evaluated once and credited to their associated regions. Overlapping
regions reuse exact candidate results, including identical decoded integer
points. Region origin lives in `search_origin` outside the candidate payload.
Identical physical points **within the same study/runtime definition** therefore
share candidate IDs, regardless of origin. A new refinement study does have a
new study/definition ID: its centers and search protocol participate in that
identity, so IDs are not promised to match the source study.

The existing journal/recovery/state protocol persists the in-flight candidate,
center progress, next region and each region's Sobol index. An interrupted local
evaluation remains pending, including a deadline-limited center, and is retried
before the schedule advances on resume. Completed centers are never inserted
again. Reports count evaluation attempts (including retries and cache hits);
shared centers can contribute to more than one region's summary. `status` and
the offline cockpit show per-region attempts, convergence, feasibility, best
objective/COP and rejection evidence, with exact center values in HTML details.

Global V2/V3 studies can explicitly opt into one initial evaluation:

```toml
[search]
type = "sobol"
domain = "fixed_global_bounds"
seed = 29092026
scramble = true
evaluate_initial = true
```

Absent this field, historical behavior is unchanged: the campaign starts at
Sobol index zero and does not insert the initial. Existing snapshots are never
amended on resume. Adding the flag changes study identity and requires a new
campaign. Current presets retain their historical setting; local studies always
evaluate their centers. V1 remains unchanged.

There is no adaptive allocation, local gradient optimizer or cross-study cache.
Capacity scaling of an already local study is explicitly refused, because its
embedded centers would also need a scientifically consistent transformation;
rescale the global source first, then refine the resulting evaluation. No
thermodynamic equation, correlation, validity threshold or constraint is changed
by this search layer.

The [bounded Human Cell check](../outputs/human_cell_stage2c/README.md)
reproduces the atmospheric center exactly (COP 1.1220429906193439) and records
one additional feasible local Sobol point. Its 6-minute budget stopped after
7 attempts; 32 was a cap, not an achieved sample count. The small observed
feasibility rate is evidence of local sampling, not a convergence claim.

## V1 fixed-pair validation study

The following walkthrough remains valid for `sixbar-thermo5d` (schema 1).
Use the linked V2 reference to configure other families or different ownership.

## Start from the checkout

No installation is required. In the commands below, `research` is a shell helper
for the module entry point:

```sh
research() { PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m dada_solver.research "$@"; }
research init sixbar-thermo5d --output outputs/my_research/study.toml
research validate outputs/my_research/study.toml
```

An installed package also exposes the equivalent `dada-research` command.
The wheel includes both the template and reference basis; it does not need the
checkout's historical `outputs/` or `examples/` directories to run this study.
NumPy and SciPy are required. The existing optional Numba backend accelerates
evaluation when installed; its actual use or Python fallback is recorded.
Rendering HTML requires no plotting package, server or internet connection.

`init` creates an editable TOML file and a frozen `study.basis.json` beside it.
It refuses to overwrite either file. The basis contains both mechanisms,
resolved hardware and configuration, initial state, wall capacities, numerical
settings and historical source digests. The source paths are provenance, not
live dependencies. The digest in TOML detects changes to the copied basis.

`validate` checks schema, units, ownership, source integrity and cheap geometry
construction without a periodic solve. Use `validate ... --json` for the full
validation/runtime metadata.

## Configure the first experiment

Edit the five `[[parameters]]` tables in `study.toml`:

| TOML parameter | Meaning | Unit | Short name for `evaluate --set` |
| --- | --- | --- | --- |
| `volume.swept_ratio` | Small / large swept-volume ratio, at fixed total swept volume | 1 | `swept_ratio` |
| `microtube.heat_in.tube_count` | H_i parallel tubes | integer count | `n_i` |
| `microtube.heat_in.tube_length_m` | H_i tube length | m | `length_i_m` |
| `microtube.heat_out.tube_count` | H_o parallel tubes | integer count | `n_o` |
| `microtube.heat_out.tube_length_m` | H_o tube length | m | `length_o_m` |

`lower` and `upper` bound the search; `initial` defines the default standalone
evaluation and the unevaluated initial reference in campaign reports. Sobol
does not insert that initial point. Counts must be integers, not `4029.0`.
The integer decoder uses nearest-even rounding, with inclusive bounds and
half-width endpoint bins; it is not a uniform categorical sampler.

The pair, gas inventory, source temperatures, frequency, total swept volume,
clearance ratios, transport model and finite external-air flow are fixed in
this protocol. Those assumptions are visible under `[fixed]` and `[policies]`.
Changing them requires a new explicit basis/study rather than editing a label
while retaining incompatible physics. The adapter changes exchanger hold-up,
wall capacity, heat transfer and losses through geometry and scales valve CdA
with tube count. Those derived quantities are not independent parameters.

The V1 parity preset retains five explicit historical study constraints. Their thresholds are configurable;
units and meanings are validated. After changing bounds, constraints or physical
inputs, start a **new campaign directory**. `resume` reads its stored definition,
not edits to the original TOML. Presentation names, input path relocation and
per-invocation time/count budgets do not change the scientific study identity.

## Evaluate a configuration without optimization

```sh
research evaluate outputs/my_research/study.toml --reference --output outputs/my_research/reference.json
research evaluate outputs/my_research/study.toml --reference --set n_i=4100 --output outputs/my_research/trial.json
research compare outputs/my_research/reference.json outputs/my_research/trial.json --html outputs/my_research/comparison.html
```

Open `comparison.html` in a browser. `--reference` uses the stored historical
champion's five coordinates. Without it, evaluation starts from TOML initials.
Repeated `--set` accepts either short names or full parameter names. Requested
physical floats are preserved exactly; normalization is used for identity and
distance, not to perturb the requested configuration through a decode roundtrip.
An existing evaluation output is never overwritten.

Standalone evaluations are immutable JSON artifacts, not Sobol journal entries.
They retain the complete scientific definition, candidate identity, source
lineage, physical values, statuses, constraints and numerical diagnostics.
They can be inspected with `status`, `report` and `compare`. Their default time
limit is 180 s; override with `--budget` if needed. A stopped standalone
evaluation must be requested again. A deadline result is never labelled feasible.

**Physical reproduction means the same configuration gives the same physical
results within the declared numerical tolerance.** No agreement with the old
optimization order, incumbent updates, shrinking radii or final best candidate
is claimed or required. Historical candidate IDs remain lineage references;
new campaign IDs also cover normalized coordinates, definition and runtime.

## Run, stop and resume a small search

```sh
research run outputs/my_research/study.toml --directory outputs/my_research/campaign --budget 2m --max-candidates 1
research resume outputs/my_research/campaign --budget 2m --max-candidates 1
research status outputs/my_research/campaign
research report outputs/my_research/campaign --html outputs/my_research/campaign.html
```

`--max-candidates` is a cap for this invocation, not a lifetime total. Both the
wall-clock budget and that cap apply. The initial evaluation estimate defaults
to 30 s, so a smaller budget can legitimately start no candidate. Subsequent
estimates use the existing robust integration-time history.

Ctrl-C leaves the campaign's pending candidate recoverable. Budget interruption
retains the last complete cycle as an initial guess where available; partial
cycles cannot establish convergence. Use a larger budget to deliberately retry
the latest deadline-interrupted evaluation:

```sh
research resume outputs/my_research/campaign --budget 5m --retry-incomplete
```

Stopping is cooperative at solver progress boundaries. A per-candidate deadline
and a campaign deadline apply together. At most one explicitly recorded safe
initial-state retry follows a microtube-domain error. Its gas state is rebuilt
at fixed inventory and its wall temperatures retained; it does not change
states during a cycle or weaken periodic convergence.

The existing exact cache remains unchanged: normalized coordinates are part of
candidate identity. Different Sobol coordinates that decode to the same tube
count are not silently merged. Unexpected exceptions leave pending work for
recovery rather than producing a feasible synthetic result.

## Compare and inspect

The offline HTML contains status and constraint filters, efficiency versus
indicated-power and progress plots, and two candidate selectors. The comparison
table gives values and B − A differences for parameters, metrics and constraint
margins. It also shows achieved external/internal peak capacity-rate ratios and
minimum signed port flows. Negative local flow remains visible. The recorded
reflux flag uses a declared 1e-8 kg/s threshold; raw extrema are retained.

Candidate detail exposes cycle count, normalized periodic error history, warm
source, retry, backend and phase timings. The physical boundary is always
external-stream heat input. Missing metrics display as unavailable. A failed
solver or unavailable constraint cannot enter the feasible ranking. Saved
configuration and source digests remain inspectable in the provenance section.

To select two candidates explicitly, copy IDs from `status` and use full IDs
or unambiguous prefixes:

```sh
research report outputs/my_research/campaign --candidate best --candidate CANDIDATE_ID --html outputs/my_research/selected.html
```

`report` and `compare` never run physics. They read completed orphan records and
ignore a torn final append for inspection without repairing or writing to the
journal; execution resume performs the existing recovery. Runtime/source changes
are labelled while inspection stays available. Execution resume requires the
recorded compatible runtime and code. If scientific study identities differ,
cross-study values are labelled and no combined ranking is assigned.

Trajectories, synchronized animation, a graphical editor, local refinement and
other study protocols are deferred. No renderer silently replays missing data.

## Scope of this reference machine

The validation preset is historical: **298.15/558.15 K**, 2 Hz, fixed gas mass
0.0012027347024506055 kg, and a 25 W indicated-power floor. It is separate from
the current 298.15/598.15 K demonstrator research envelope and approximately
100 W useful-output goal. The raw base configuration still contains a 598.15 K
reservoir seed; the actual injected wall source uses 558.15 K inlet air. Reports
use the actual exchanger inputs. No metadata normalization changes the model.

Mechanical losses remain unknown, useful output unavailable, and external-air
aerodynamic losses/fan power excluded from the balance rather than physically
zero. Finite air flow, thermal-film resistance, pause heat transfer, wall
storage and signed local reflux remain active. The large-cylinder maximum
enclosed volume, including clearance, is checked against the 0.066 m³ ceiling.
Mechanism lengths are in crank-radius units; no physical stroke or manufacturing
scale is inferred. This ceiling belongs to the V1 study; generic V2/V3 studies
may omit or choose another ceiling. Isothermality diagnostics are no longer computed.

## Acceptance and the next review

See the [validation record](DADA_ENGINE_RESEARCH_VALIDATION.md) for measured
results and the ready-to-open first-study artifacts in
[`outputs/research_first_study/`](../outputs/research_first_study/).

The reference and two changed configurations are compared with the original
evaluator under matching inputs, warm states and numerical settings. The
wrapper tolerance uses the existing backend-equivalence scale (relative 2e-11,
absolute 1e-11), not a physical-model accuracy claim. Pure persistence tests
separately verify split-run order, exact cache and pending-candidate recovery.
Real Sobol serialization is covered in addition to mocked runner tests.

Run the checks from the source checkout:

```sh
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 -m pytest -q
```

Before adding other campaigns, the user must try this first study and assess:

1. Can the five variables and units be located and changed without Python?
2. Do validation errors make an invalid bound or count easy to correct?
3. Are exact evaluation, search, stopping and resume clearly distinguished?
4. Does comparing two candidates explain the efficiency/power change and the
   constraints or numerical failures that matter?

Record that feedback before broadening the module. Automated and agent-operated
smoke checks are not a substitute for this real researcher usability review.

## Filling at a reference pressure and maximum total gas volume

V2/V3 retain `charge = "explicit_inventory"` unchanged. To derive inventory from
each candidate's final geometry instead, remove every `charge.total_mass_kg`
parameter declaration and use:

```toml
[policies]
charge = "reference_pressure_at_maximum_total_volume_v1"
# Keep the study's other policies unchanged.

[charge_reference]
pressure_pa = 100000.0
temperature_k = 293.15
volume_state = "maximum_total_gas_volume"

[warm_start]
initial_source = "uniform"
```

Pressure is absolute. This policy currently requires the calorically perfect
working gas; it does not apply an ideal filling approximation to real-gas tables.
It derives `m = p_ref * max_theta(V_total(theta)) / (R * T_ref)` after construction
of the candidate kinematics and connection of both actual exchanger designs.
Both cylinders use the **same angle**. The production model's four gas volumes
already include its clearances, exchanger gas, headers and additional hold-up;
configuration HX seeds are replaced, not added a second time.

The basis inventory remains historical input metadata and is overwritten before
integration. The uniform initial gas temperature is `T_ref`; at the solver's
starting angle zero the filling pressure can differ from `p_ref`. The reference
pressure applies at the maximum-volume position, not at every shaft position.
Reservoir/external temperatures and all periodic-convergence criteria are unchanged.
`source_exact` is rejected. Ordinary nearby candidate guesses are rescaled to the
new inventory (gas mass and energy together), retaining specific internal energies
and, for wall models, wall temperatures.

The versioned numerical method uses 8192 equal angular intervals, declared
kinematic breakpoints, bracketed roots of the total-volume derivative and local
bounded refinement of sampled peaks. It is independent of `screening.samples`.
All families share this method, including harmonic laws whose analytic maximum
is used as a test oracle. It is a deterministic numerical search, not a proof of
the global maximum for arbitrarily narrow/pathological motion features. Roots
and local refinements use a 1e-13 rad absolute target; maxima within 2e-14 relative
volume are tied by the lowest solver angle. The reported volume is evaluated at
that reported angle. Flat total-volume laws choose angle zero. Changing this
algorithm requires a deliberate policy/method version change.

`derived.charge` records `policy`, `reference_pressure_pa`,
`reference_temperature_k`, `reference_total_gas_volume_m3`, `reference_angle_rad`,
`derived_total_mass_kg` and `reference_volume_method`. The angle is the solver
angle after the existing operation-direction transform. The mass appears in
`metrics.total_mass_kg` and the HTML inventory column even when subsequent
integration fails; it is not evidence of a converged cycle. Earlier geometry
rejections cannot have a computed reference volume.

The policy and reference settings participate in scientific study/candidate
identity. Existing explicit-inventory studies retain their scientific identity;
normal source/runtime compatibility guards still apply. Capacity rescaling keeps
the reference pressure/temperature and derives mass again from the scaled final
geometry, rather than introducing an independent mass parameter.

### Explicit design limits in new studies

A stricter Mach requirement is optional and distinct from the selected
exchanger's correlation domain:

```toml
# Example study requirement, not a generic model limit.
[[constraints]]
type = "maximum_mach_number"
limit = 0.20
unit = "1"
```

The screen only needs its numerical resolution. Declare an enclosed-volume
ceiling only if the study requires one:

```toml
[screening]
samples = 1440
# Optional study requirement; values above 0.066 m^3 are allowed.
# maximum_large_enclosed_volume_m3 = 0.100
```

New six-bar presets expose engineering defaults directly in the study, for
example the following existing row can be edited or removed without regenerating
the mechanism artifact (do not add a duplicate row):

```toml
[[mechanical_constraints]]
side = "small"
metric = "minimum_primary_transmission_sine"
relation = "minimum"
limit = 0.30
unit = "1"
```

A bound on a searched coordinate belongs in that parameter's declaration;
it does not automatically create a feasibility constraint.
