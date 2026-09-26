# Dada-Engine Research — first implementation proposal

Status: proposed API and schema for human review, not installed commands.
This plan follows the [local audit](DADA_ENGINE_RESEARCH_MIGRATION_MATRIX.md)
and preserves the guarantees of [OPTIMIZATION_CAMPAIGN.md](OPTIMIZATION_CAMPAIGN.md).
The application name is **Dada-Engine Research**, independent of temperature.

## First deliverable

One fixed independent six-bar pair, five active thermal/hardware coordinates,
explicit historical assumptions, a bounded persistent Sobol campaign, and
offline inspection/comparison of two candidates. Users edit TOML, not Python.
V1 adds no local refinement, dynamic plugin registry, database, GUI server,
assumed mechanical efficiency, new exchanger physics, or mandatory symmetry.

The regression preset reproduces the historical hlat25 rank-01 machine at
298.15/558.15 K. The current demonstrator preset will separately use
298.15/598.15 K, the 2–10 Hz envelope, and the 0.066 m³ large-cylinder enclosed
volume ceiling. Neither preset may turn the approximately 100 W useful-output
goal into an indicated-power equivalence. Useful power remains unavailable
until a mechanical-loss model or measurements exist. Isothermality remains
diagnostic-only. Independent free spline campaigns continue to work unchanged.

## Minimal module and API boundary

| Proposed entry point | Responsibility and reused implementation |
| --- | --- |
| `research.schema.load_study(path) -> StudyDefinition` | Strict versioned TOML validation; units, active/fixed ownership, explicit policies and source digests; no integration or output writes |
| `research.schema.compile_study(study) -> CampaignDefinition` | Produce a canonical existing-engine definition with a six-bar family and a normalized frozen basis artifact |
| `research.sixbar.load_basis(path, expected_sha256) -> SixBarStudyBasis` | Validate typed resolved configuration, both mechanisms, branches, source selectors, units and optional warm state; no runtime imports from `examples` |
| `research.sixbar.SixBarStudyAdapter.build(physical) -> MachineDesign` | Apply the five coordinates and declared couplings, then construct production `IndependentSixBarVolumeKinematics` and `MicrotubeExchanger` objects |
| `research.cli.main(argv=None)` | `init`, `validate`, `run`, `resume`, `status`, `report`; delegate execution to `OptimizationCampaign` |
| `research.report.render(directory, destination)` | Read saved records and study lineage; produce JSON/text/offline HTML without loading an evaluator |

Start with these few modules rather than a protocol/plugin hierarchy. A thin
explicit six-bar branch at the campaign construction boundary is acceptable;
no concrete family checks belong in thermodynamic integration code.

The registered entry point should be `dada-research`, since the checkout has
separate `dada-solver` and `dada-optimize` entry points and no `dada` dispatcher.
Keep both historical commands compatible. Proposed usage:

```sh
dada-research init sixbar-thermo5d --output study.toml
dada-research validate study.toml
dada-research run study.toml --directory outputs/research_rank01 --budget 2m --max-candidates 2
dada-research status outputs/research_rank01
dada-research resume outputs/research_rank01 --budget 2m --max-candidates 2
dada-research resume outputs/research_rank01 --budget 5m --retry-incomplete
dada-research report outputs/research_rank01 --html /tmp/rank01.html
```

`init` refuses to overwrite a study and bundles immutable reference inputs next
to it. `validate` checks data and backend availability and performs cheap
geometry construction; it must not start a periodic solve. `status` and `report`
remain available after code/runtime changes, clearly labelling those changes.
`resume` must enforce the original exact runtime/definition identity.
Ctrl-C and cooperative budget exhaustion use the existing runner's recovery
rules; document that stopping is cooperative, not instantaneous.

## Versioned TOML and parameter semantics

[proposed_rank01.toml](research_audit/proposed_rank01.toml) is a complete parseable
**schema proposal**, not an input accepted by the current CLI. It pins the
captured reference bundle digest and explicitly records:

- schema version, study protocol and historical purpose;
- the immutable source bundle and historical parent candidate;
- objective boundary and historical constraints;
- fixed mass, total swept volume, independent fixed pair and clearances;
- finite external-air boundary, loss exclusions and count-scaled outlet valves;
- ordered active coordinates, units, bounds, decoding and initial values;
- Sobol seed/scrambling, wall settings source, backend, retry and warm-state policy.

The initial parameter values are reference metadata, not a Sobol point. The
stored champion is an explicit regression input, not a silently injected seed.
Path resolution is relative to the study file. Unknown keys, unknown units,
nonfinite values, noninteger counts, duplicate ownership, inconsistent fixed
values versus the resolved source, and missing diagnostics fail validation.
No arbitrary expressions or Python evaluation are supported in TOML. Derived
policies are a small named, versioned set with explicit equations.

Source temperature validation must use the resolved exchanger air inputs. The
historical base configuration retains a 598.15 K reservoir seed although the
injected heat-in wall uses 558.15 K air. Snapshot both facts, label the actual
boundary correctly, and test that any later metadata normalization does not
alter model diagnostics or physical results.

| Coordinate | Owner | Unit / decoding |
| --- | --- | --- |
| `volume.swept_ratio` | Six-bar study volume partition | `1`, continuous linear |
| `microtube.heat_in.tube_count` | Heat-in geometry | `1`, integer nearest-even |
| `microtube.heat_in.tube_length_m` | Heat-in geometry | `m`, continuous linear |
| `microtube.heat_out.tube_count` | Heat-out geometry | `1`, integer nearest-even |
| `microtube.heat_out.tube_length_m` | Heat-out geometry | `m`, continuous linear |

Put the reusable `IntegerParameter` beside `ContinuousParameter` in
`campaign/parameters.py`. For v1, decode valid normalized `u` with
`int(round((1-u)*lower + u*upper))`, Python's ties-to-even rule. Bounds and
initial values must be strict integers (reject bool), counts positive, bounds
ordered. Encode `k` by `(k-lower)/(upper-lower)`; test endpoints and every
half-integer tie. This matches historical rounding semantics, but a bounded
global Sobol search still does not reproduce incumbent-centred proposals.
Endpoint bins are half-width; do not claim a uniform categorical distribution.

Keep existing continuous TOML decoding and payload hashes unchanged. The
decoded count is a JSON integer. Preserve normalized coordinates in
`Candidate` identities: two coordinates mapping to the same integer hardware
are still different campaign candidates under the current contract. V1 accepts
that possible redundant evaluation; it does not introduce a hidden decoded-only
or approximate cache. Existing exact duplicate reuse stays exact.

For total swept volume `Vt`, ratio `r`, and clearance ratios `cS`, `cL`:

```text
Vs = Vt*r/(1+r); Vl = Vt/(1+r)
Vmin,S = cS*Vs; Vmax,S = (1+cS)*Vs
Vmin,L = cL*Vl; Vmax,L = (1+cL)*Vl
CdA_side = source_CdA_side * count_side / source_count_side
```

Reuse `_volume_limits` as the parity reference when extracting these equations.
Reject concurrent direct swept/clearance coordinates when this partition policy
owns them. Fixed inventory excludes varying charge pressure. UA, hold-up,
wall capacity and hydraulic losses remain derived by the exchanger model;
they cannot also be independent coordinates. No physical mechanism scale or
piston force can be inferred from dimensionless crank-radius geometry.

## Extend the campaign instead of copying it

1. Extend `CampaignDefinition`'s parameter factory with explicit `kind`, defaulting
   to historical continuous behavior. Accept six-bar basis settings only for
   the six-bar family. The source bundle must be a normalized typed artifact,
   not an arbitrary historical report with a guessed `best` selector.
2. Extend adapter ownership and dispatch at construction only. The six-bar
   adapter owns geometry and coupled volume partition; exchanger count/length
   belong to microtube geometry. Avoid forcing a pair into
   `SharedCrankRockerDesign`, which has different ownership and assumptions.
3. Include basis contents/digest and study ID in **new-family** definition
   identity. Preserve legacy payload schemas and fingerprint algorithms; do
   not rewrite old histories to match current source fingerprints. Package code
   changes already invalidate execution resume by design. Old journals remain
   inspectable; resume requires a compatible recorded runtime/source checkout.
4. Add basis snapshot writing/loading to creation and resume. Snapshot required
   source artifacts before publishing the definition manifest, under the same
   campaign lock. Resume must resolve the stored snapshot, not the original
   mutable path. Validate source hash changes before expensive construction.
5. Retain `OptimizationCampaign`, `CampaignHistory`, `Candidate`,
   `SobolStrategy`, `make_report`, atomic writes, fsync, locks, cache behavior,
   orphan recovery and incomplete-attempt semantics. Do not implement another
   journal or round-robin controller for the first single-family campaign.
6. Reuse `MachineEvaluator`, objective/constraint objects and wall solver.
   Add explicit opt-in policies only where needed for parity: provenance-checked
   external initial state, per-candidate deadline, and one bounded safe-state
   retry after `MicrotubeDomainError`. A retry stays within both deadlines,
   preserves inventory, and cannot claim convergence. Persist its source and
   outcome; unrelated programming errors must propagate.
7. The replay acceptance fixture supplies the **exact historical predecessor**
   state and wall capacities. Search runs may use the existing compatible
   nearest-state policy; that is a recorded numerical-policy difference from
   the script's seed-normalized thermal distance. Never synthesize a fake
   completed campaign candidate merely to inject a warm state.

Historical `interrupted` maps to `budget_exhausted`, not physical rejection;
`converged` maps to `feasible` or `converged_infeasible` after all assessments.
Retain the raw legacy status in reference comparisons. The first constraints
are efficiency availability, 25 W indicated floor, 1.2 MPa pressure ceiling,
850 K gas-temperature ceiling, 0.08 kg/s flow ceiling and model validity.
Mechanical closure is checked before integration. The hlat25 synthesis limits
are inherited evidence for this frozen pair, not universally applied thermal
constraints or continuous feasibility proofs. If geometry becomes active in a
later study, its mechanical margins must become explicit evaluated constraints.

## Provenance and reporting

Use a separate `study_id = content_hash(scientific_definition)` and explicit
historical parent IDs. The scientific definition includes schema version,
source bytes/digests, family/policy versions, parameters, constraints, thermal
boundaries and numerical settings. Exclude report colors, titles, output paths
and per-invocation time budget from this new study identity. Preserve the
existing campaign identity semantics independently, including runtime identity.

Keep `study.toml`, `study.json`, source snapshots, and an immutable manifest
beside the existing campaign files; reports are replaceable derived outputs.
Each selected artifact manifest has schema version, candidate/study/definition
IDs, source digests, units and whether it was stored or recomputed. Portable
input bundles must include both resolved physical settings and source lineage;
do not import the long chain of examples in the production CLI.

The first HTML report embeds all data/assets locally and shows identity,
assumptions, unavailable useful power, candidate statuses/counts, feasible
ranking by the declared objective, elapsed evaluation progress, efficiency vs
indicated power, constraint margins and parameter differences for two selected
candidates. A numerical failure is visible, not a plotted zero. Report
generation never invokes `run()` or a physical evaluator. Escape report strings
and JSON embedding so names cannot inject HTML. Presentation changes must not
invalidate scientific identity.

Phase 2 adds `TrajectoryDiagnostics` (a new distinct name) and `MechanismFrames`:
angle domain/direction, array shape, named series and units, boundary/sign
definitions, availability/validity, selected geometry, named joints/links,
scale and layout disclosures. Adapt production replay/geometry; frontend code
does not reconstruct thermodynamic or linkage laws. Full trajectories are
saved only for selected candidates. Missing data shows an explicit unavailable
message; only an explicitly requested compatible-runtime replay may create it.
Detailed cylinder P–V plots remain work diagnostics; a global pressure plot
cannot silently replace individual pressure-work integrals.

## Ordered work and acceptance gates

| Step | Files to add/change | Checks and gate | Expected computation / risk |
| --- | --- | --- | --- |
| 0: review this package | Audit, matrix, snapshots and this plan only | Review historical preset choice and bounded Sobol scope | Completed; no migration |
| 1: definition and integer space | `campaign/parameters.py`, `definition.py`; `research/schema.py`, `research/sixbar.py`; compact `tests/fixtures/research/` input bundle | Strict schema/units/ownership; tie and endpoint tests; input digest mismatch; legacy candidate payload regression | Seconds, no integration; risk of identity drift |
| 2: construction and evaluator parity | `campaign/adapters.py`, opt-in hooks in `evaluator.py`; six-bar adapter | Compare volumes/derivatives on full-cycle grid, both mechanisms/branches, charge mass, hardware/valves/air conditions; reject bad geometry before solver | Cheap construction; one baseline plus two varied thermal points, allow up to 180 s each |
| 3: durable CLI | `campaign/runner.py` snapshot plumbing, `definition.py` resume; `research/cli.py`, `pyproject.toml` | Tiny two-candidate real smoke; fake evaluator split-run/order equality; kill/retry; duplicate and rejected cache; orphan/torn-tail/lock tests; no files required outside bundle | Mostly seconds; real smoke budget 2–6 min; warm-state/runtime compatibility risk |
| 4: offline inspection | `research/report.py`, README/user guide | Two saved candidate views; report must pass with solver patched to fail if called; absent metrics; incompatible runtime still readable; offline/headless HTML check | Seconds; prevents accidental expensive replay |
| 5: acceptance | `tests/test_research_schema.py`, `test_research_sixbar.py`, `test_research_cli.py`, `test_research_report.py`, existing tests | Same inputs produce matching metrics and constraint statuses; source/runtime mismatch fails resume; no legacy regressions | Full suite plus bounded parity runs; do not promise a universal run time |
| 6: phase 2+ | Diagnostics/frame export, then sweeps, staged protocols, robustness | One independent historical parity reference per protocol; animation follows actual geometry; no perturbation retuning | Separate scope; full optimization may take hours |

All tests run from the checkout with `PYTHONPATH=src python3 -m pytest`.
Separate mock persistence tests from physical parity tests. Mock success is not
physical evidence. Check missing metrics, nonconvergence, domain failure and
budget exhaustion explicitly. Pin exact initial states when comparing stopping
criteria; a different warm start can legitimately change cycle counts.
Use existing six-bar/backend equivalence tests to set tolerances, then record
the chosen per-metric thresholds and actual differences.

The measured original rank-01 replay took 9.466 s (including first-use compiled
backend overhead), versus a stored 3.380 s historical evaluation. This suggests
a small smoke can be inexpensive, but invalid or slow regions can consume the
full candidate budget. Record setup, screening, integration, replay and rendering
separately. No performance claim is inferred from one converged warm start.

## Review boundary

The initiating audit's section 10 explicitly says: “Stop for human review before
writing the new implementation or cleaning the repository.” This package makes
that review concrete. Approval would authorize steps 1–5 above: a fixed-pair
thermo-5D CLI using the existing campaign engine, with historical physical parity
and offline reports. It would not authorize deleting examples, changing physical
equations, reproducing a full local search campaign, or treating the historical
25 W indicated floor as the demonstrator's useful-power requirement.
