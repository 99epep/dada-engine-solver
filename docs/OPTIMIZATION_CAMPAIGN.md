# Campaign internals

For day-to-day work, use [Research](DADA_ENGINE_RESEARCH.md). This is the
lower-level implementation reference, also used by the historical `dada-optimize`
CLI. Research snapshots have their own filenames described in its reference.

## Scope

The campaign is a sequential orchestration layer above the existing physical
and sizing evaluators. It reuses `evaluate_configuration`, objective objects,
constraint objects and `SizingAssessment`. `SizingProblem`, its local SLSQP
optimizer and the existing Latin-hypercube feasibility tools remain available
and unchanged. No physical equation, sign convention, valve rule, exchanger
correlation, four-bar mathematics or periodic stopping tolerance is replaced.
No mechanical efficiency is assumed; indicated power remains distinct from
unknown useful mechanical output.

New code lives in `dada_solver.campaign`:

- `parameters`: immutable continuous, integer and choice parameters and parameter spaces;
- `candidate`: canonical JSON and SHA-256 candidate identities;
- `adapters`: family-specific construction and parameter ownership;
- `definition`: TOML loading, snapshots and runtime/model fingerprints;
- `strategy`: the minimal strategy protocol and incremental Sobol exploration;
- `evaluator`: preflight, existing physical evaluation and raw sizing assessment;
- `history`: durable records and recovery;
- `runner`: budgets, cache, sequence continuation and feasible archives;
- `report`: JSON/text phase summaries and deterministic human suggestions;
- `cli`: the `dada-optimize` command.

There is no database, plugin discovery, distributed execution, Bayesian search
or new external optimization dependency. The first strategy is Sobol, not
SLSQP. Research V2/V3 additionally use `scheduled_search.ScheduledSobol` for
explicit initial evaluations and center-first local regions; the original global
Sobol path stays unchanged. See [Research local refinement](DADA_ENGINE_RESEARCH_REFERENCE.md#local-refinement-and-explicit-initial-evaluations)
for normalization, round-robin allocation and resume semantics.

## Kinematics injection

`MachineDesign(configuration=..., kinematics=..., heat_in=..., heat_out=...)`
now accepts a `KinematicsModel` directly. Existing positional exchanger
arguments remain compatible; injection is an optional additional field.
`build_model(configuration, kinematics=...)` skips the family selector when
an implementation is supplied, checks its cylinder ranges and calls its
optional feasibility validator. Without injection, the existing configuration
and `kinematics_type` path is unchanged.

**Inject study-angle motion, before the operation-direction transformation.**
For negative-speed motor operation, the construction boundary reverses it
exactly once. Existing four-bar objects retain their concrete type through
the historical crank-direction replacement; other implementations use the
generic reversal wrapper. Do not inject a kinematics object extracted from an
already motor-transformed model. Explicit `ReversedVolumeKinematics` injection
is rejected to catch that mistake. A future mechanism can be injected without
changing the thermodynamic integrator or adding a family switch there.

`evaluate_configuration(..., model=...)` accepts the already constructed model
instead of rebuilding it. This preserves the physical evaluator while allowing
injected motion to reach it. Configuration-only charging-pressure objectives
use the evaluated model when available, so they do not rebuild another motion
family merely to obtain the initial filling volume.

## Parameters and ownership

Each `ContinuousParameter` has a stable name, lower/upper bounds, an initial
value and `linear` or `log` transform. Bounds are finite, strictly ordered and
contain the initial value. Log bounds must be positive. `ParameterSpace`
requires unique names and normalized coordinates in `[0,1]`; invalid values
are rejected, not clipped. Endpoint decoding is exact. Interior decoding is:

```text
linear: x = (1-u)*lower + u*upper
log:    x = exp((1-u)*log(lower) + u*log(upper))
```

Integer parameters use the shared encode/decode convention; choice parameters
store explicit categories, not scientific numeric distances. Candidate physical
values may therefore be strings.

Old global studies without `evaluate_initial` retain the historical Sobol-only
start. Research studies can explicitly evaluate their initial point first; local
regions always evaluate their embedded centres before Sobol. Each local region
has its own Sobol index and deterministic round-robin allocation. Search-origin
metadata does not alter physical candidate identity. The study's search definition
and durable schedule still participate in reproducible execution.

Adapters avoid expanding the legacy enum:

- `common.<existing DesignParameter value>` uses `DesignPoint` for operating,
  charge, cylinder ranges and compatible legacy exchanger inputs;
- `operation.frequency_hz` maps positive frequency to the base configuration's
  existing signed angular-speed convention;
- `free.small.<index>` and `free.large.<index>` vary independent stereographic
  shape coordinates; other coordinates remain explicitly fixed in `[free]`;
- `four_bar.<dimensionless field>` updates the shared-crank design while keeping
  discrete assembly branches fixed;
- `microtube.heat_in.<geometry field>` and `microtube.heat_out.<geometry field>`
  are owned by the separate microtube geometry adapter. Tube count can be an integer
  coordinate. Circular collectors derive pitch, header volume and diode section;
  independent legacy UA/hold-up/CdA coordinates cannot override that geometry.

Wrong-family or unknown names are rejected. Frequency and angular speed cannot
both vary; neither can pressure and inventory, or clearance volume and ratio.
Microtube ownership explicitly forbids independent legacy UA, exchanger gas
volume and equivalent inlet hydraulic-resistance coordinates. Its geometry
adapter derives components with the existing `MicrotubeExchanger`; it does not
fit independently free UA or hold-up values.

### Physical evaluators

Reservoir closures use the existing eight-state periodic solver. Microtube
exchangers construct the existing `AirWallMotor` and use the reusable ten-state
wall-cycle evaluator factored from hardware screening. The extra states remain
H_i and H_o wall energies. The wall convergence rule, correlations, valve
equations and integration tolerances are unchanged.

For wall cycles, indicated thermal efficiency is indicated gas work divided by
external heat supplied through the external stream. Wall-to-gas heat is a
separate diagnostic. Useful shaft power remains unavailable and mechanical
losses remain unknown. Peak tube Reynolds and Mach numbers use the actual
candidate geometries; actual selected correlation-domain violations prevent
feasibility; modeled transition is reported with its uncertainty.

`WallCycleNumericalSettings` records the wall integration method, relative
tolerance, fifteen-component absolute-tolerance vector, maximum angle step,
periodic relative/absolute scales, progress interval and bounded Aitken option.
The method and angle step default from unambiguous generic numerical settings;
wall-specific values remain explicit. These values are stored in candidate
payloads and the immutable campaign identity. Historical defaults remain LSODA,
relative tolerance 1e-8, the established componentwise absolute tolerances and
maximum step pi/360.

In wall-cycle diagnostics, `cold_heat_rate_extrema` and
`hot_heat_rate_extrema` mean heat into the working gas from the H_i and H_o
walls. They are evaluated with the actual wall states and exchanger `rates()`
functions. Cycle efficiency continues to use heat delivered by the external
stream. These are deliberately different thermodynamic boundaries.

The optional `wall_numerical.accelerate_walls` applies only the existing bounded
Aitken update to a subsequent wall-energy initial guess. It cannot declare
convergence; the next ordinary complete physical cycle must satisfy the same
periodic test.

## Candidate identity and reproducibility

A `Candidate` stores immutable canonical JSON containing normalized coordinates,
decoded physical values, fixed family identifiers, numerical settings and a
campaign definition identifier. Canonical JSON uses sorted keys and rejects
nonfinite numbers. SHA-256 replaces Python's randomized `hash()`. Reconstructed
candidate payloads are checked against their identities.

The definition identifier covers the parsed campaign, the complete base TOML
snapshot and Python/NumPy/SciPy versions plus a digest of the package's Python
sources. Resume refuses changed definitions or runtimes rather than silently
mixing records or treating changed physics as cached evaluations. Move deliberate
bounds, frozen-parameter, fidelity or model changes to a new campaign directory.
Compatible warm-start policies are explicit; they never bypass periodic convergence.

Candidate identity includes normalized coordinates as well as decoded values.
No approximate/geometric-equivalence cache is attempted. Exact duplicates reuse
the persisted result, including rejection results, without integration.

## Sobol continuation

`SobolStrategy` uses `scipy.stats.qmc.Sobol` with explicit dimension, seed and
scrambling. Every active coordinate varies in the same sampled vector. The
state stores the number of consumed points; resume reconstructs the same engine
and calls `fast_forward(index)`. The journal also records each sequence index,
allowing recovery when the state snapshot lags a completed evaluation.

Arbitrary budget endpoints use incremental `random(1)`. Power-of-two prefixes
have the usual Sobol balance property; stopping at an arbitrary count does not
claim that property. No initial point is skipped, thinned or replaced. Runtime
versions are pinned by the campaign fingerprint because implementation changes
must not silently alter sequence continuation. See the
[SciPy Sobol documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.qmc.Sobol.html).

## Evaluation and statuses

Each candidate is decoded and constructed through its family adapter. Cheap
parameter, geometric and free-motion derivative-limit checks happen before the
periodic solve. `MachineDesign` then builds the injected model, which is passed
to the existing evaluator. Raw objective values and every raw constraint margin,
availability flag and satisfaction flag are retained separately.

Statuses distinguish:

- `invalid_parameterization`;
- `invalid_kinematics`;
- `invalid_exchanger`;
- `invalid_fluid_domain`;
- `integration_failure`;
- `periodic_non_convergence`;
- `converged_infeasible`;
- `feasible`.

A feasible record requires periodic convergence, all requested constraints and
an available finite objective. No fabricated penalty replaces an unavailable
objective or a failed solve. Unexpected programming errors propagate, leaving
the pending candidate recoverable, rather than masquerading as physical
infeasibility. Geometry-definition errors and derivative-limit margins remain
explicit in the record.

Reusable states are persisted after convergence, maximum-cycle exhaustion and
deadline interruption whenever at least one complete cycle exists. They carry
an explicit `periodic_solution` flag. Compatible nearest-state selection and
inventory/wall-capacity rescaling are described below; convergence tolerances
are never relaxed.

## Files, durability and crash recovery

```text
campaign_directory/
    campaign.toml             # original campaign definition snapshot
    base.toml                 # complete thermodynamic configuration snapshot
    definition.json           # immutable content/runtime identity
    state.json                # Sobol index, pending candidate, phase and archive IDs
    history.jsonl.gz          # new campaigns: append-only gzip members, one full JSON record each
    recovery.json            # one transient durable completed result before state commit
    candidates/<sha256>.json  # legacy only; retained/read, never newly generated
    report.json
    report.txt
    reports/phase_0001.json
    reports/phase_0001.txt
    ...
```

Before integration, `state.json` records the pending candidate and advanced
sequence index. After evaluation, its full candidate record is written using
an atomic rename and filesystem sync to `recovery.json`, then the journal is
appended and synced, then the state/archive snapshot is updated and synced.
Only then is recovery removed and the directory synced. A completed recovery
record (or legacy candidate file) absent from the journal is recovered on restart;
an already journaled recovery is acknowledged without duplication. An unfinished pending candidate may
be retried; a completed result is not silently re-integrated.

The sole journal-edit exception is recovery of a torn final append: its bytes
are saved to `history_torn_tail_*.bin`, the incomplete tail is removed, and any
completed recovery record or legacy candidate file is recovered. Non-final corruption fails explicitly.
Thus a process crash loses at most the in-flight evaluation. This assumes a
local filesystem honoring the sync/atomic-rename operations; it is not a
replicated storage system. A POSIX advisory lock prevents concurrent writers.

Cache hits append a complete attempt referring to the original evaluation,
using the same single recovery slot. Legacy candidate files are not rewritten. Archives are reconstructed from history on resume.
They contain distinct feasible candidate identities ranked by minimum objective,
not the last iterate. The phase report also preserves the best at phase start,
best at phase end, and the best candidates evaluated in that phase.

## Time budget and human reports

The budget is a per-run approximate elapsed-time allowance. Before starting a
candidate, the runner compares remaining time with 1.1 times the 75th percentile
of the last ten uncached integrations. Until such timings exist, it uses explicit
`initial_evaluation_seconds`. Cheap preflight rejections do not make the next
periodic solve appear artificially cheap. A zero or insufficient budget starts
nothing. A cooperative deadline is checked from the RHS progress boundary.
`deadline_grace_seconds` is configurable; the default is 12 seconds and the
microtube smoke uses 3 seconds with one-second progress checks. An interrupted
partial `solve_ivp` cycle is discarded, while the last complete cycle and its
error history are retained. Its status is `budget_exhausted`, never physical
infeasibility or periodic non-convergence, and it is excluded from the exact
result cache. Exploration advances normally; `--retry-incomplete` deliberately
retries the latest unfinished candidate with a later, larger budget.

Warm-start selection searches compatible prior states in normalized space and
prefers converged states before distance. Compatibility includes state layout,
evaluator family and operating direction. Gas inventory is rescaled while
preserving specific internal energies. Wall energy is scaled by the new/old
wall-capacity ratio to preserve wall temperature. Source identity, status and
normalized distance are persisted. Every candidate must still satisfy the
original periodic criterion.

JSON and readable text reports include requested/actual duration, attempted and
cached counts, preflight rejection/integration/convergence/feasibility counts,
phase-start/end bests, objective improvement, physical metrics, individual
constraint margins, top distinct feasible candidates and timing statistics.
Parameter changes include physical start/end values and normalized deltas.
Repeated proximity to either bound uses a declared 0.05 normalized threshold.
Failure statuses, detailed reasons and violated/unavailable constraints are
counted separately.

Suggestions are deterministic observations, never campaign edits. Rules cover
few feasible points and dominant failures; improving elites at a bound; tight
elite clustering; competitive distant regions; and a phase with little
improvement and positive listed margins. Positive margins are not automatically
called comfortable engineering margins. Thresholds and evidence are exposed in
the report. Bounds, assumptions, freezing and eventual local refinement remain
human decisions. No evidence of global optimality is inferred from a small run.


## Verification and legacy entry points

`tests/test_campaign.py` and the Research tests exercise identity, recovery,
cache, interruption and deterministic local/global continuation. `dada-optimize`
accepts a campaign TOML for creation or an existing campaign directory for resume.
Use Research for new declarative multi-family studies.

No test-count snapshot or old smoke performance is a guarantee for a changed
property model. See [validation](validation.md) and
[current limit ownership](PHYSICS_DECISIONS.md) before interpreting old verdicts.
