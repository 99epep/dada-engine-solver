# Dada-Engine Research — Architecture Audit and Implementation Plan

**Status:** design proposal for Codex; no implementation authorized by this document.  
**Repository:** `99epep/dada-engine-solver`  
**Inspected branch:** `master`, tree `5e94be1e9f1f8aa70e0bad3a618378dcc237550a` (26 September 2026).  
**Product name:** **Dada-Engine Research**. “Dada” is the inventor's surname; preserve the product's capitalization and hyphenation.

## 1. Executive decision

Build **Dada-Engine Research** as a thin, declarative research layer **on top of** the existing `dada_solver.campaign` and physical-model packages, not as another optimization framework. Retain existing numerical solvers, mechanical definitions, validity logic, cache/replay infrastructure, and the existing campaign journal. Make research strategies, study protocols, evaluators, result views, and renderers separately replaceable. Use TOML and a CLI first; generate standalone HTML reports with optional static publication exports. Defer a graphical editor, distributed scheduler, and database.

**First milestone:** reproduce one existing six-bar thermomechanical study from a declarative definition, pause/resume it, and inspect/compare its candidates without editing Python. Do not migrate or delete legacy examples before parity is established.

## 2. Evidence and audit boundaries

This audit is based on the current GitHub tree, the directory inventories of `examples/`, `src/dada_solver/`, `src/dada_solver/campaign/` and `docs/`, and direct inspection of representative implementation and documentation: `campaign/{cli,definition,runner,report,strategy,history,evaluator}.py`, `cli.py`, `plotting.py`, `mechanism_animation.py`, `docs/OPTIMIZATION_CAMPAIGN.md`, `docs/COMPACT_KINEMATIC_SYNTHESIS.md`, and representative six-bar thermo and robustness scripts. The inventory is comprehensive at filename level, but **not every example was read line by line**. Therefore the detailed migration mapping below is a prioritized hypothesis, not a claim that all scripts have been semantically audited. Codex should verify each group before extraction. The Git tree includes `__pycache__` entries; treat those as cleanup candidates only after a separate tracked-file and compatibility check.

The repository already documents its campaign subsystem in `docs/OPTIMIZATION_CAMPAIGN.md`. That document is authoritative on its original scope; this proposal extends it and must not silently reinterpret its guarantees.

## 3. Existing assets and actual gaps

| Existing asset | Verified capability | Gap for Research |
|---|---|---|
| `campaign/parameters.py`, `adapters.py`, `definition.py` | Named continuous parameter space, family ownership, TOML loading, snapshots and runtime fingerprints | Discrete/count coordinates, derived/coupled parameters, six-bar geometry ownership, reusable study presets |
| `campaign/strategy.py`, `runner.py` | Incremental Sobol; time budgets; sequential persistent execution; continuation | Strategy selection, multi-stage orchestration, local refinement, nested/multifidelity workflows |
| `campaign/candidate.py`, `history.py` | Canonical identities, persistent history and recovery | Standardized cross-campaign lineage and selected-candidate artifact registry |
| `campaign/evaluator.py` | Physical evaluation, preflight and explicit failure statuses | Mechanism-only, coupled thermo-mechanical, tolerance and sweep study evaluators |
| `campaign/report.py` | Feasible elites, time statistics, bound-pressure and deterministic human suggestions | Cross-campaign comparisons, Pareto/constraint trade-off views, provenance-rich visualization |
| `plotting.py` | Thermodynamic diagnostic plots from a completed cycle | Reusable diagnostics schema, selective replay and interactive plotting |
| `mechanism_animation.py` | Existing four-bar GIF generator | General six-bar geometry/frame export, reusable animation API, synchronized plots |
| `examples/` | Large body of experimental searches, plots, refinements and animations | Duplicated research orchestration and hard-coded paths/parameters; not yet an end-user API |

**Important:** `campaign` is not merely a prototype: it already implements important correctness contracts. Keep its candidate identities, immutable campaign definitions, exact duplicate cache, interrupted-candidate behavior, explicit feasibility statuses and reproducible Sobol continuation. Never claim that a result is feasible when a diagnostic is unavailable or a periodic solve has failed.

## 4. Inventory of experimental workflows

The following groups are derived from filenames and representative scripts. Some scripts may belong to more than one group; Codex must inspect their actual imports, I/O and algorithms before consolidation.

| Workflow family | Representative examples | Proposed reusable primitive |
|---|---|---|
| Thermal parameter exploration and temperature maps | `run_motor_temperature_map.py`, `optimize_motor_temperature_point*.py`, `evaluate_motor_temperature_map_baseline.py` | Grid/sweep study and bounded local refinement; parent/child campaign linkage |
| Free motion and law fitting | `optimize_motor_free_spline_260k*.py`, `optimize_motor_fourier_c2_260k.py`, `optimize_motor_hybrid_c2_15p_260k.py`, `fit_motor_hybrid_c2_*` | Motion-law adapters, shape constraints, fitting evaluator |
| Phase and valve studies | `optimize_motor_four_stage_*.py`, `screen_motor_valve_placements_300k.py`, `optimize_motor_valve_topologies_260k.py` | Continuous + enumerated parameter studies; comparative topology reports |
| Exchanger and hardware studies | `optimize_motor_exchanger_asymmetry_stageK*.py`, `motor_hardware_screening.py`, `compare_microtube_gas_model.py` | Exchanger parameter adapters, hardware validity and coupled evaluator |
| Four-bar synthesis and diagnostics | `search_compact_motor_motion.py`, `evaluate_motor_champion_four_bar_k2.py`, `animate_motor_champion_four_bar.py` | Mechanism synthesis adapter, mechanical screening, geometry/animation export |
| Six-bar synthesis and staged search | `search_large_primary_cadence_3952*.py`, `search_large_second_dyad_v41_panel_3952.py`, `search_small_sixbar_stage2*.py`, `optimize_sixbar_pairs_thermo_3952*.py` | Six-bar adapters, hierarchical stages, multiple families and candidate lineage |
| Coupled thermo-mechanical refinement | `optimize_sixbar_pairs_thermo5d_3952_hlat25.py`, `optimize_rank01_compact_piston_rods*.py` | Frozen/active parameter groups, constrained refinement and coupled evaluation |
| Robustness | `test_sixbar_robustness_3952.py` | Tolerance perturbation protocol with cheap mechanical screening followed by selected thermo replay |
| Diagnostics and benchmarks | `benchmark_solver_stage*.py`, `benchmark_periodic_anderson.py`, `diagnose_periodic_map.py`, `audit_wall_acceleration.py` | Separate developer benchmark suite; diagnostic export, **not** user-facing optimization strategies |
| Presentation and animation | `plot_motor_*.py`, `plot_sixbar_3952_comparison*.py`, `animate_sixbar_pair_v*.py`, `visualize_large_second_dyad_v41_panel_3952.py` | Stable renderers and reusable data views; retain publication-specific layouts as optional presets |

### 4.1 Highest-value extraction targets

1. **Six-bar paired-candidate serialization and construction.** `optimize_sixbar_pairs_thermo5d_3952_hlat25.py` reconstructs `SixBarCylinderMechanism` from a JSON candidate, builds `IndependentSixBarVolumeKinematics`, holds working-gas inventory constant, modifies exchanger tube count/length and scales outlet valve CdA with tube count. These are physical *study policies*, not generic assumptions. Express them explicitly and snapshot them.
2. **Multi-stage search protocol.** Existing scripts demonstrate primary-geometry searches, secondary-dyad searches, paired thermo searches, local refinement and multi-family comparison. Preserve stage boundaries and source-candidate lineage; do not pretend that one flat continuous optimizer reproduces these procedures.
3. **Robustness protocol.** `test_sixbar_robustness_3952.py` uses a 30-dimensional Sobol perturbation of both mechanisms, separate length/angle tolerance rules, mechanical pre-screening, then thermodynamic replay of a subset with no reoptimization. Encode these as a distinct study type with explicit perturbation distributions, seeds and acceptance definitions. Its acceptance rate is a *screening statistic*, not certified industrial reliability.
4. **Diagnostics and presentation.** Extract physics-derived arrays once; make all plotting and animation consume the same versioned diagnostic/frame representation. Do not have each plot script reimplement mechanical geometry or thermodynamic state extraction.

### 4.2 Do not extract indiscriminately

The many `animate_sixbar_pair_v3_*` and `v4_*` files are valuable layout history but should not become dozens of public renderers. First identify the latest authoritative geometry renderer and separately preserve any publication-specific layout choices. Likewise, acceleration and solver benchmark examples belong in development/test tooling unless they expose a genuine reusable diagnostic.

## 5. Proposed architecture

```text
src/dada_solver/
    campaign/                  # Existing engine; extend minimally
    research/
        schema.py              # Versioned study definition and validation
        presets.py             # Named example study templates
        protocols/             # sweep, search, refine, compare, robustness
        adapters/              # six-bar and additional study-specific adapters
        orchestration.py       # Stages, dependencies, child campaigns
        artifacts.py           # Manifest and selected-candidate artifact references
        diagnostics.py         # Typed, units-aware analysis datasets
        compare.py             # Reproducible candidate/campaign comparisons
        visualization/
            plots.py
            mechanisms.py
            report.py
        cli.py
```

This is a **responsibility map**, not a mandate to create one file per name. Codex should prefer fewer modules where existing APIs already suffice. In particular, extend `campaign` rather than copying its runner, identity, history or evaluator. `research` owns study composition, data presentation and high-level user workflows.

### 5.1 Contracts

- **Study definition:** validated TOML; named base configuration; study type; active, fixed and derived parameters; bounds, transforms and discrete domains; objective(s), constraints, physical assumptions, strategy, seed, fidelity and time budget. Include explicit units and reject ambiguous ownership.
- **Strategy:** resumable state with `ask`/`tell`-like semantics or an equally minimal incremental protocol. Wrap existing Sobol first; later add local refinement using installed SciPy facilities where appropriate. Do not change the existing Sobol sequence behavior.
- **Evaluator:** normalized request + model/geometry construction + cheap preflight + physical evaluation; return raw metrics, constraint margins, failure status and optional reusable converged state. No arbitrary penalty should turn failures into feasible candidates.
- **Protocol:** orchestrates sweeps, stages, family comparisons and tolerance replay by composing campaigns. A robustness protocol is **not** a search for a repaired optimum.
- **Artifact:** immutable manifest with definition/runtime fingerprint, parent IDs, source-file digests, unit conventions and schema version; JSON/JSONL for metadata and compact array storage for large trajectories. Avoid persisting full trajectories for every rejected candidate.
- **Renderer:** read-only, deterministic consumer of stored diagnostics; where absent, selectively replay the exact selected candidate with a compatible model/runtime and explicitly mark recomputed data.

### 5.2 Data identity and reproducibility

Preserve the existing campaign definition/candidate hash contracts. Add a separate `study_id` and explicit parent/child lineage rather than retroactively changing old candidate IDs. Distinguish immutable *run identity* from presentation choices (plot colors, layout, resolution). Keep original TOML snapshots and source-data SHA-256s, particularly for historical `outputs/` dependencies. Changing physical model, tolerances, hardware, parameterization, target law or constraint semantics creates a new scientific study identity; changing only report appearance does not.

For integer tube counts, define rounding/encoding rules once and record the **decoded integer**; otherwise nominally distinct continuous candidates may evaluate identical hardware. Avoid approximate cache reuse across nonidentical physics. Reuse existing exact kinematic caching, compatible warm-state logic and shared diagnostic replay only through their validated interfaces.

### 5.3 Physical semantics that must survive migration

- Maintain the solver's operation-direction convention: inject study-angle motion and apply motor reversal **exactly once** at the established model boundary.
- Keep indicated gas work/power distinct from unknown useful shaft power. Preserve the correct external-source heat boundary for wall-coupled thermal efficiency.
- Do not silently substitute fixed inventory for fixed charging pressure, or hold exchanger hold-up/UA constant while changing microtube geometry.
- Preserve discrete assembly branches and check full-cycle closure and mechanical transmission constraints before expensive thermodynamic integration.
- Distinguish sampled geometric screens from continuous mathematical proofs and numerical convergence from physical design feasibility.
- Preserve the precise constraints and assumptions of each historical comparison rather than globally adopting the 3952-specific limits.

## 6. User experience: CLI first, HTML second

Illustrative **proposed** commands, not current interfaces:

```bash
dada research init motor-sixbar --output study.toml
dada research validate study.toml
dada research run study.toml --budget 2h
dada research status outputs/my_study
dada research resume outputs/my_study --budget 1h
dada research report outputs/my_study --html
dada research compare outputs/study_a outputs/study_b --html
dada research animate outputs/my_study --candidate best --format mp4
```

Use an appropriate registered console entry point or an existing CLI dispatcher; do not force a breaking change to `dada-optimize`. Error messages should identify invalid parameters, absent source artifacts, incompatible physical models and missing diagnostic data before costly computation.

### 6.1 Reports and visualization

A useful first HTML report includes: study identity and assumptions; candidate count/status breakdown; best feasible candidates; objective vs evaluations/wall time; power vs efficiency scatter with constraint filters; active bounds/constraint margins; per-phase progress; parameter changes; numerical convergence diagnostics; and links to detailed candidate pages. Distinguish unavailable metrics from zeros. A report must be useful offline.

Selected-candidate detail views should expose 0–360° piston position/normalized volume and derivatives; valve states; pressures and temperatures; mass flow; heat rates; cumulative work/heat; mechanism transmission and rod angles where supported; and a synchronized mechanical animation. A renderer must use the **actual selected mechanism and geometry**, with any schematic enlargement or nonmanufacturing layout prominently disclosed. Export SVG/PNG/CSV plus GIF/MP4 where dependencies permit. A lightweight standalone HTML report can use embedded plot data; avoid a server requirement for v1.

Use a canonical `CycleDiagnostics`-like representation (angle, named arrays, units, boundary definitions, validity flags) and `MechanismFrames` (angle, named joints/links/sliders, axis orientation, geometric scale). Build adapters from existing production evaluators; never recompute physical laws inside the frontend.

## 7. Implementation sequence and acceptance gates

### Phase 0 — audit completion and baseline preservation (no refactor)

Codex should enumerate **every** Python/TOML example, its imports, CLI flags, input dependencies, output formats, optimization algorithm, feasibility logic, plotting/animation code and documentation references. Build a machine-readable migration matrix with columns: `file`, `category`, `canonical_successor`, `shared_logic`, `study_specific_logic`, `inputs`, `outputs`, `tests`, `migration_status`. Explicitly locate current historical output datasets; some examples use `Path.cwd()/outputs` and import other scripts as modules. Snapshot at least one reproducible representative from each major family before moving files.

**Gate:** migration matrix reviewed; no existing scripts deleted; reference commands and checksums recorded; clear distinction between verified code facts and filename-derived hypotheses.

### Phase 1 — smallest end-to-end vertical slice

Implement a declarative study based on one existing six-bar thermomechanical workflow (prefer a small, representative version of the final 3952 thermo-5D study, not an expensive full reproduction initially). Use the existing campaign journal/runner wherever possible; add only the necessary six-bar adapter and discrete tube-count handling. Run short, deterministic tests with a tiny candidate budget, then verify selected evaluations against the original script under identical solver settings.

**Gate:** a non-Python user can configure, validate, run, stop and resume the study; same candidate parameters give matching physical metrics, status and constraints; provenance identifies all required source reports.

### Phase 2 — diagnostic schema and report

Extract reusable diagnostics from existing `plotting.py` and six-bar plotting/animation examples. Generate offline HTML for one campaign and two selected candidates. Render animations without coupling them to optimizer internals.

**Gate:** report generation does not rerun optimization; missing trajectory data triggers only an explicitly requested, provenance-checked selective replay; figure labels and physical units are validated; animation uses the same geometry as the selected candidate.

### Phase 3 — protocols and cross-study analysis

Add sweep/temperature-map, staged search/refinement and robustness protocols in that order. Preserve independent families and their history. Add campaign comparison with explicit checks for model/hardware/temperature comparability; label incompatible comparisons rather than plotting misleading direct rankings.

**Gate:** parity against one reference study per protocol; reproducible seeded sampling and deterministic resume; robustness never silently retunes a perturbed candidate.

### Phase 4 — consolidation and cleanup

Only after parity tests: migrate canonical examples to `examples/research/` or declarative presets, archive historically significant experiment scripts, retire truly redundant render variants, and remove generated/tracked caches with an explicit file-by-file review. Update documentation and entry points. Never delete historical output dependencies before snapshots/migration.

## 8. Verification plan

1. **Unit tests:** TOML schema, parameter ownership and unit validation; discrete encoding; objective/constraint availability; strategy state serialization; diagnostic schema.
2. **Determinism:** same seed/definition/runtime yields the same Sobol candidates; split run plus resume matches uninterrupted candidate order; exact duplicate cache remains exact.
3. **Physical parity:** identical geometry, gas inventory, exchanger hardware, valve placement and numerical tolerances produce matching evaluations between original example and Research wrapper. Set tolerances according to existing tests and numerical backend, not invented universal decimal thresholds.
4. **Failure and interruption:** invalid geometry rejected before thermal integration; numerical failure is not disguised as physical infeasibility; interrupted candidates can be retried and reports remain readable.
5. **Provenance:** source-file hash changes are detected; old studies remain inspectable; new incompatible studies do not silently share cached physical results.
6. **Visualization:** independently checked angles/units/signs, endpoint periodicity, mechanism branch consistency, static export, animation-frame correspondence, and a headless CI smoke test.
7. **Performance:** benchmark against the existing script, separately recording setup, cheap screening, physical integration, diagnostics replay and rendering. Preserve acceleration gains; do not make every candidate generate full trajectories.

## 9. Risks and non-goals

- **Overengineering:** no initial GUI framework, service, SQL database, task queue or arbitrary Python plugin system.
- **Premature generality:** do not make every historical experiment a first-class strategy. Consolidate common operations and keep rare studies as explicitly documented specialized protocols.
- **Historical output dependency:** a script that imports another example or reads `outputs/...` cannot be migrated safely by moving it alone.
- **Numerical equivalence:** changed solver settings or cached/replayed states may change results; record them and verify parity rather than assuming it.
- **Scientific communication:** multiobjective trade-offs, robustness screens and physically incomparable studies must not be collapsed into a single unexplained 'champion' score.

## 10. Immediate instructions for Codex

**First assignment: complete and validate this audit; do not implement the full module yet.**

1. Read this document and `docs/OPTIMIZATION_CAMPAIGN.md`; inspect the current tree, campaign package, relevant tests and all `examples/` sources. Correct any hypotheses here that do not survive code inspection.
2. Produce `docs/DADA_ENGINE_RESEARCH_MIGRATION_MATRIX.md` (or CSV plus short MD), mapping every example to a canonical workflow or justified historical/developer category. Mark shared logic and the exact implementation to reuse; flag import/path and output-file dependencies.
3. Propose a minimal concrete API and versioned TOML schema for the first six-bar thermo-5D vertical slice. Explicitly explain how to extend rather than duplicate the current campaign infrastructure, and where integer coordinates and six-bar ownership belong.
4. Identify exact reference input artifacts and available regression outputs for the first slice. If they are not in the repository, report this as a blocker and define a small checked-in fixture rather than inventing historical results.
5. Provide an ordered implementation plan with files to touch, tests to add, acceptance gates, expected computational cost and risks. **Stop for human review before writing the new implementation or cleaning the repository.**

The end-user product is **Dada-Engine Research**. Its primary purpose is to let a researcher formulate, run, understand, compare and extend rigorous machine studies without editing Python, while retaining full scientific traceability and the user's control over engineering decisions.
