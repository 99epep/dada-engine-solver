# Persistent optimization campaigns

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

- `parameters`: immutable named continuous parameters and parameter spaces;
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
SLSQP. A later strategy can implement the same incremental protocol, but local
refinement is not implemented here.

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

The initial values are available as a reference but are not automatically
inserted into or substituted for the Sobol sequence. If there is no feasible
candidate at phase start, parameter-change reporting explicitly uses these
unevaluated configured initial values; objective improvement is unavailable.

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
  are owned by the separate microtube geometry adapter. Integer tube count is
  fixed in this first continuous parameter representation.

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
external heat supplied through the incoming air stream. Wall-to-gas heat is a
separate diagnostic. Useful shaft power remains unavailable and mechanical
losses remain unknown. Peak tube Reynolds and Mach numbers use the actual
candidate geometries; leaving the declared laminar/Mach domain prevents
feasibility.

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
source air. These are deliberately different thermodynamic boundaries.

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
Nearest-state transfer between compatible campaigns is a later feature.

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
    history.jsonl             # append-only attempt/result journal
    candidates/<sha256>.json  # durable original result for each unique candidate
    report.json
    report.txt
    reports/phase_0001.json
    reports/phase_0001.txt
    ...
```

Before integration, `state.json` records the pending candidate and advanced
sequence index. After evaluation, its full candidate record is written using
an atomic rename and filesystem sync, then the journal is appended and synced,
then the state/archive snapshot is updated. Completed candidate files absent
from the journal are recovered on restart. An unfinished pending candidate may
be retried; a completed result is not silently re-integrated.

The sole journal-edit exception is recovery of a torn final append: its bytes
are saved to `history_torn_tail_*.bin`, the incomplete tail is removed, and any
completed candidate file is recovered. Non-final corruption fails explicitly.
Thus a process crash loses at most the in-flight evaluation. This assumes a
local filesystem honoring the sync/atomic-rename operations; it is not a
replicated storage system. A POSIX advisory lock prevents concurrent writers.

Cache hits append an attempt referring to the original evaluation without
rewriting its candidate file. Archives are reconstructed from history on resume.
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

## Running the small physical example

`examples/free_kinematics_campaign.toml` varies eight parameters simultaneously:
two shape coordinates per cylinder, speed, charge pressure, large swept volume
and small clearance volume. Each cylinder has six spline control points (four
chart coordinates), with the remaining chart coordinates explicitly fixed.
The small initial region uses the smooth independent laws from the existing
free-motion integration example; this smoke seed does not impose a preferred
shape on the free-kinematics family.

The example retains 25/325 deg C, roughly 1 L, and 2–2.2 Hz. It caps each
physical evaluation at three complete cycles to keep the architecture smoke
modest. The periodic tolerances are unchanged; failure to converge within this
cap is recorded as non-convergence, not admitted as a feasible result. Increase
the cycle cap in a new campaign definition for a substantive physical search. Its objective is
the existing indicated thermal efficiency. The 1 W minimum motor-power
constraint is only a smoke-test feasibility condition; it does not replace the
project's approximately 100 W useful-output goal. UA/CdA values are legacy
screening inputs, not the latest geometry-connected exchanger design.

Installed command:

```sh
dada-optimize examples/free_kinematics_campaign.toml --directory outputs/free_campaign_final_smoke --budget 3m --max-candidates 1
dada-optimize outputs/free_campaign_final_smoke --budget 3m --max-candidates 1
```

From a source checkout, replace `dada-optimize` with
`PYTHONPATH=src python3 -m dada_solver.campaign.cli`. General budgets can be
`30m`, `1h30m`, `45s`, or seconds. The example definition itself caps each run
at two candidates unless `--max-candidates` overrides it. An uncapped production
campaign can omit `maximum_candidates` in its definition.

## Verification boundary and next steps

Fast tests use a lightweight evaluator and simulated clock for transforms,
hashing across processes, sequence continuation, cache hits, durable history,
preflight rejection, feasibility archives, report deltas, bound pressure,
budgets, crash recovery and changed-definition rejection. Injection tests cover
both motion families in both directions. The physical example is a separate
CLI smoke run, not a slow full campaign in the ordinary test suite.

Next review: local refinement, robustness studies and optional frozen external
reference warm starts. Local SLSQP is deliberately absent.
Do not infer an optimal waveform or useful shaft efficiency from these tests.

### Recorded fast-suite verification

The complete suite passes: **287 tests**. The last full run took 46.44 seconds.
The fixture budget tests use a simulated
clock and do not sleep. Exact duplicate lookup, interrupted in-flight recovery,
completed-file recovery after a torn append, feasible-best preservation and
cross-process candidate hashing are covered explicitly. Added wall-campaign
checks cover hardware snapshots, actual exchanger composition, deadline grace,
retry semantics, last-complete-cycle retention, warm-start compatibility and
wall-capacity rescaling, external-heat efficiency and report availability.

### Recorded physical smoke verification

The source CLI completed Sobol point 0 in
`outputs/free_campaign_final_smoke`: 81.52 seconds of evaluation, three
integrated cycles, and explicit `periodic_non_convergence`. The requested
phase budget was 180 seconds. No objective or feasible optimum is claimed;
the cycle cap did not relax any convergence tolerance. The completed record,
candidate payload and first phase report are durable.

A separate process resumed this directory and selected Sobol point 1 rather
than repeating point 0. At this verification checkpoint it was still in flight
after 738 seconds; its pending payload and advanced sequence index were
persisted. Completed physical evaluation after resume is therefore not yet
verified by this run. Automated continuation and interrupted-state recovery
tests pass independently. Inspect `state.json` and the append-only history for
the subsequent outcome; `reports/phase_0001.json` retains the completed first
phase even if the latest report changes.

This second point demonstrated substantial integration-time variability and
motivated cooperative in-flight deadlines. The current runner interrupts at a
progress boundary, discards the partial cycle and retains the preceding complete
cycle. Numerical-cost diagnostics still deserve review before long physical
exploration; no solver equation was changed to hide this limitation.

### Dynamic-wall microtube smoke

`examples/microtube_free_campaign.toml` composes independent free motion with
the parallel 3708-tube hardware snapshot at 25/325 deg C. It varies two free
coordinates, both exchanger lengths, speed, charge pressure and large swept
volume. The base maximum of 100 periodic cycles is retained; wall-clock
deadlines control interactive runtime.

Two distinct Sobol candidates were evaluated across resumed processes in
`outputs/microtube_free_campaign_correction_smoke`. Point 0 completed eight
cycles in 31.19 s for 30 s requested with a 3 s grace. Its normalized periodic
error fell from 464423.79 to 8159.15. Point 1 used point 0 as a compatible
non-periodic initial guess at normalized distance 1.5255 and completed thirteen
cycles in 42.50 s for 40 s requested with the same grace. Its error started at
127870.87 and reached 4207.97 at cycle 10. The configured bounded Aitken update
then changed the wall-energy guess: the next errors were 63529.05, 1939.03 and
27.36. The ordinary cycle-13 error remains above the convergence threshold 1.
Both statuses are `budget_exhausted`; neither is feasible or presented as an
optimum. Their hashes are distinct, sequence indices are 0 and 1, and the next
persisted Sobol index is 2.

The established parallel-microtube regression was rerun from its saved state:
37.6261105666 W indicated power, 353.5945084 W external heat input and
0.1064103363 indicated thermal efficiency. This agrees with the recorded
37.626110567 W result.
