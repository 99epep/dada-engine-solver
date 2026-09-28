# Dada-Engine Research — local migration audit

Status: V1 completed and user-tested; V2 kinematics/mechanism extension implemented.
The local checkout, including local outputs and untracked inputs, is authoritative.
The CSV/inventory have been regenerated from this working tree. Hashes identify
actual inspected content; the source commit alone does not identify uncommitted work.
The original phase-0 workflow decisions below remain historical design context.
No historical example or output has been deleted.
The compact-hybrid follow-up refreshes the two affected inventory rows and adds
the eleventh Research family; unrelated local work is outside that update.

## Reviewed V2 mathematical families versus search protocols

The explicit [reviewed status map](research_audit/kinematics_v2_migrations.json)
is applied by the inventory generator. It distinguishes extracted representation,
existing production reuse and retained search/plotting code. AST imports alone
never mark a script as parity-verified. The [V2 reference](DADA_ENGINE_RESEARCH_KINEMATICS.md)
documents the parameter schema and available execution boundary.

| Historical lineage | Mathematical representation | V2 status / protocol distinction |
| --- | --- | --- |
| Harmonic production law | Cosine volume vs study angle | Reused; independent phases |
| `compare_slider_crank_motion`, `evaluate_slider_crank_k2` | Centered/offset finite crank-slider, explicit inverted/conventional volume direction | K2 physics extracted to `slider_crank.py`; fitting/search remains historical |
| `compact_coupler_geometry`, finite-coupler K2 evaluation | Four-bar finite-rod slider output in a reconstructed common crank frame | Production closure reused; portable normalized seed and retained local envelope frame; projection-only legacy variant remains outside this V2 family |
| `six_bar`, documented rank 1/4/12/50 geometries | Primary four-bar + dyad + finite rod, 15 coordinates per piston | Declarative artifacts/ownership and mechanical screens; discovery, saturation and local fitting execution deferred |
| Free local F1 and free-spline V1/V2/V3, spline refinement/from-linear | Same canonical periodic cubic spline plus phase | One `free_spline` family. V2/V3 scalar phase wrapper extracted; V1 slower wrapper retained. Radii, seeds, initial projection, safe starts and staged refinement are search protocols |
| Fourier optimization and refinement 4H/6H/8H | Same smooth Fourier representation with variable harmonic count | One `fourier_c2` family. Class extracted; harmonic count is a representation setting, successive local bounds/radii are protocols |
| Piecewise P1/P3 | Same `IdealPiecewiseLinearVolumeKinematics` | One `ideal_piecewise` family; P1/P3 bounds and feasibility policies are not new motion laws |
| Four-stage K2, symmetry4D, thermo3D, HX9D, refinements | Four linear segments, shared historical timing/origin | Existing `four_stage_kinematics` exposed per cylinder; dimensional/hardware passes and symmetry constraints remain study choices |
| Independent four-stage 11D | Independent side timing plus small cyclic origin | Distinct `independent_four_stage` representation already in production, now exposed |
| Hybrid compact optimization | `HybridCompactKinematics`: rounded linear BP branches and quintic HP kinks (9 coordinates) | Extracted to `hybrid_compact_kinematics.py` as the 11th Research family; frozen dense and candidate-501 thermal parity verified; historical search retained |
| C2 11p fitting | `CompactC2FitKinematics`: extrema curvature with split quintic kink branches | Genuinely distinct 11-coordinate fitting representation; retained in examples, extraction deferred |
| Source16 C2 15p fit and thermodynamic optimization | `StructuredMotion15` + scalar `StructuredKinematics15`: BP bowing and variable-width HP kink with global extrema curvature | Extracted to `structured_kinematics.py`; historical candidate 3952 reproduced; least-squares fitting and historical local search loops retained |

Changing V1/V2/V3 search implementation, Sobol seed, bounds, local radius or seed
source does not create a new mathematical family. Conversely, the two older
compact hybrid laws above are not aliases for the extracted 15p representation.
The nine-coordinate compact hybrid law is now extracted; the distinct 11p
fitting representation remains example-only. No production Research module imports `examples`. Remaining example-only logic
includes the older 11p compact fitting law, legacy fitting/search protocols, saturation
statistics, plotting/layout and old basis-construction utilities. Their presence
is recorded explicitly rather than claiming complete campaign migration.

## Inventory and evidence

The [CSV matrix](DADA_ENGINE_RESEARCH_MIGRATION_MATRIX.csv) has one row for every
Python/TOML example: **241 files, 188 Python and 53 TOML**, including the locally
untracked bounded V2 demonstration script. It supplies all nine requested migration columns,
plus hashes, tracking state, and documentation references.

The [detailed inventory](research_audit/example_inventory.json) records full-source
AST/TOML evidence: imports, transitive imports of other examples, CLI declarations
with defaults, input/output call sites, path expressions and artifact literals,
optimization calls, feasibility functions, rendering calls, and direct test
references. Line numbers refer to the hashed files. Empty test lists mean no
direct reference was found, not that production functionality has no coverage.

All files were parsed without executing them. Categories and successors are
proposals; static analysis does not establish scientific equivalence, resolve
arbitrary dynamic paths, or prove which branch a run used. `open()` is listed
under both I/O categories when direction is not statically classified. Inherited
CLI and feasibility logic must be followed through `transitive_example_imports`.
The original semantic review covered the first slice; the explicit V2 status
map now adds the reviewed motion and mechanism lineages above. No row is
marked migrated or parity verified merely because its filename looks relevant.

Regenerate the static inventory from any directory:

```sh
python3 docs/research_audit/generate_inventory.py
```

The script uses its own checkout location. Generated review documents are excluded
from documentation-reference discovery, so the inventory does not cite itself.
This is developer audit tooling, separate from the executable Research application.

## Workflow decisions

| Family | Canonical reusable implementation | Study policy to retain | Migration gate |
| --- | --- | --- | --- |
| Thermal exploration / temperature maps | `campaign.runner.OptimizationCampaign`, `optimize_motor_temperature_point.build_design` as extraction reference | Separate five-coordinate hardware and seven-coordinate motion passes; charging-pressure policy and per-temperature power floor | Deferred sweep protocol; preserve source-temperature lineage |
| Free motion / fitting | `free_kinematics`, family-specific motion classes; SciPy fit calls inventoried per file | Independent S/L laws, normalization, extrema rules, parameter bounds | No universal symmetry or shape assumptions |
| Phase / valve studies | Existing valve construction and physical evaluators | Enumerated topology, source-dependent constraints and distinct histories | Compare compatible physics before ranking |
| Hardware trials | `exchangers.hardware`, `microtube`, `wall_cycle`, `gas_diagnostics` | Finite inlet-air flow, gas-property closure, valve CdA rule, fan exclusion, inventory policy | Preserve thermal boundary and geometry-derived quantities |
| Four-bar synthesis | `four_bar`, `mechanical_optimization`, `compact_coupler_geometry` as reference | Target curve, scale, assembly choices, sampled screening thresholds | Verify geometry/derivative parity before extraction |
| Six-bar synthesis | `six_bar.SixBarCylinderMechanism`, `IndependentSixBarVolumeKinematics` | Separate primary/dyad stages, branch choices, family-specific mechanical limits | No flat optimizer substitution |
| Coupled refinement | Existing geometry and wall evaluators | Which mechanical/thermal coordinates are active; compact-rod and lateral-motion limits | Deferred; not a first-slice strategy |
| Fixed-pair thermo-5D | `optimize_sixbar_pairs_thermo5d_3952_hlat25.build_design` as parity oracle; campaign runner as persistence engine | Fixed mass and independent geometry, total swept volume, clearances, five coordinates, count-scaled outlet CdA | First slice proposed below |
| Robustness | `test_sixbar_robustness_3952.perturb_pair`, `screen_pair`, `stable_subset` as protocol references | Shared 30-D Sobol perturbations, distinct dimensional/angular errors, no retuning | Separate acceptance statistics; not certified reliability |
| Diagnostics / animation | `results.extract_cycle_diagnostics`, `diagnostic_replay.replay_wall_trajectory`, production geometry | Coordinate conventions, valve signs, heat boundary, layout disclosure | Read-only reports first; explicit replay only |
| Solver acceleration | Existing backend, exact kinematic cache and shared replay | Numerical-equivalence cases and timing instrumentation | Remain developer benchmarks |
| Refrigerator / load / sizing TOML | Existing configuration and sizing/load loaders | Historical operating conditions and objectives | Retain; do not force into motor presets |

For every individual file, the CSV lists exact production imports and example
dependencies. These are existing reuse points; proposed Research destinations
are explicitly identified as future work.

## Corrections to the original audit

1. **The physical six-bar implementation already exists.** Reuse `six_bar.py`,
   including branch/closure checks and refined stroke normalization. The missing
   piece is campaign ownership/serialization, not new linkage mathematics.
2. **Warm states, deadlines and numerical acceleration are already integrated.**
   `EvaluationControl`, `select_warm_start`, `rescale_wall_state`,
   `budget_exhausted`, `--retry-incomplete`, backend selection and diagnostic
   replay must survive. The runner currently constructs Sobol internally and
   its resume classmethod constructs `CampaignDefinition` directly; evaluator
   injection alone is not enough for an independently serialized Research study.
3. **Diagnostics already have two complementary representations.**
   `results.CycleDiagnostics` stores extrema/events/topology, while
   `WallTrajectoryReplay` holds evaluation-local arrays and sampled physical
   points. It checks object identity and is not a portable persisted dataset.
   Add a serializable trajectory export around these APIs; do not create a
   second same-named `CycleDiagnostics` with conflicting meaning.
4. **The first reference is historical 25/285 deg C**, i.e. 298.15/558.15 K,
   at 2 Hz and a 25 W **indicated** floor. It is not the current 25/325 deg C
   research envelope or a demonstrated 100 W useful-output machine. A new
   598.15 K study must have a separate identity; parity retains 558.15 K.
5. **Thermo-5D is incumbent-centred local search**, using round-robin families,
   shrinking radii, an explicitly evaluated seed, rounded tube counts, skipped
   proposals and a finite Sobol array. The persistent campaign uses incremental
   Sobol over fixed bounds and does not insert its configured initial values.
   V1 must claim evaluation parity, not historical search-trajectory parity.
   Local refinement remains deferred under the project conventions.
6. **The original evaluator has policies absent from the generic adapter**:
   historical initial-state selection and one safe uniform-state retry after a
   `MicrotubeDomainError`. It uses a 180-second candidate deadline as well as
   the global deadline. These require explicit settings or an explicitly
   labelled difference; they cannot disappear in a purported parity test.
7. **The historical result ID is not a campaign ID.** The script hashes family,
   pair ID and decoded parameters; `Candidate` also hashes normalized coordinates,
   definition and numerical settings. Preserve both identities through lineage;
   do not reuse legacy IDs as campaign cache keys.
8. **Production microtube validity is model-dependent.** Current evaluator code
   uses `gas_domains` when available, not a universal Re < 2300 rule. The older
   prose in `OPTIMIZATION_CAMPAIGN.md` describes the laminar screening closure.
   Retain the selected correlation domain and unavailable-constraint handling.
9. **Examples are not isolated entry points.** The chosen script imports a
   transitive chain of example helpers, several of which import plotting
   dependencies. It also reads Stage 7A5 and thermo-3D reports to construct its
   basis. Moving the top-level script alone would break it.
10. **Later compact results and local presentation work exist.** The tracked
    `sixbar_rank01_compact_*` studies are distinct from the four-family hlat25
    reference. The untracked `export_final_compact_thermodynamic_cycle.py`
    patches globals in `export_rank01_thermodynamic_cycle`; the untracked
    `animate_sixbar_graphs_svg_v5.py` is a later presentation candidate. Neither
    is automatically the canonical renderer. Existing renderers reconstruct
    linkage geometry themselves; parity with the selected production mechanism
    must precede extraction. No user work was overwritten.
11. **There are 339 tracked cache files** under `__pycache__` at this audit point.
    This is evidence for a later file-by-file cleanup review, not authorization
    to delete them now.

The resolved historical configuration still carries a legacy 598.15 K reservoir
seed, while the actual heat-in exchanger air inlet is 558.15 K. Injected wall
hardware determines the physical source boundary. A report must display the
exchanger inlet rather than blindly reading the legacy configuration field;
normalizing that metadata must not silently change the parity evaluation.

The original script's docstring also names a nonexistent
`optimize_sixbar_pairs_thermo5d_3952.py`; the actual executable has the
`_hlat25.py` suffix. Preserve the original during this audit and use the actual
filename in new documentation.

## Reference inputs and preservation

The [rank-01 snapshot](research_audit/rank01_reference_snapshot.json) contains
the original pair and selected result, resolved machine/exchanger configuration,
all wall settings, initial/final state, source hashes, backend/runtime identity,
and the twelve actual JSON/JSONL/TOML inputs observed during reconstruction.
All twelve are present locally and tracked. There is **no missing-input blocker**
for this first slice. The snapshot is a small review artifact, not an invented
historical result or a substitute for its source lineage.

The construction dependencies are:

- `examples/motor_mechanics_stage7A5_edge.toml`, its base machine TOML and
  `motor_hardware_parallel_325c.toml`;
- Stage 7A5 `report.json` and the exact candidate file selected by it;
- `outputs/motor_four_stage_thermo3d/report.json` (`control_1bar`);
- `outputs/motor_spline_from_linear_260k/report.json` (fixed gas inventory);
- `outputs/motor_fourier_c2_8h_refine/report.json` (hardware/volume seed);
- hlat25 rank-01 `best_pair.json` and matching `report.json`;
- thermo-5D rank-01 `report.json` and `history.jsonl` for the selected result
  and its precise warm-start predecessor.

The [family reference manifest](research_audit/family_reference_manifest.json)
also preserves checksums, representative stored result excerpts, script hashes,
and reference commands across the major workflow families. These other studies
were inspected/snapshotted, **not rerun**. The commands describe the historical
entry points; some write default historical locations and must be redirected
before future execution. No example or output was moved in this phase.

## Verification performed

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_campaign.py tests/test_six_bar.py tests/test_machine_injection.py
MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:examples python3 docs/research_audit/capture_reference.py
```

The focused tests passed: **61 passed in 2.61 seconds**. The original rank-01
evaluator replay converged in 13 cycles, with the same saved warm-start state:

| Metric | Stored result | Fresh original-evaluator replay |
| --- | ---: | ---: |
| Indicated thermal efficiency | 0.23255335337157246 | 0.23255335337157246 |
| Indicated gas power | 41.4337050897658 W | 41.4337050897658 W |
| External heat input | 178.16859868523701 W | 178.16859868523701 W |
| Normalized periodic error | 0.5907169380737661 | 0.5907169380737661 |
| Feasible under historical constraints | true | true |

The six compared physical metrics have zero observed difference in this runtime.
Setup took 0.151 s after imports; evaluation took 9.466 s, including first-use
Numba overhead (no fallback calls). This is one measured baseline, not an
all-platform bitwise guarantee or evidence of Research-wrapper parity. The
wrapper has not been written. Tolerances for future cross-backend parity must
follow the existing numerical regression tests and measured error budgets.

Final artifact checks also passed: exact 240-file CSV/JSON coverage and hashes,
all captured input/source digests, all twelve family artifact digests, proposed
TOML syntax/bounds/reference digest, equality of the full compact replay result
and feasibility flag with the stored reference, and local documentation links.
`git diff --check` passed. Only audit documentation/data/utilities were added;
the original proposal received a follow-up link. Test-generated bytecode from
this run was removed; pre-existing tracked and untracked caches were preserved.

The proposed API, schema and ordered acceptance gates are in
[DADA_ENGINE_RESEARCH_IMPLEMENTATION_PLAN.md](DADA_ENGINE_RESEARCH_IMPLEMENTATION_PLAN.md).
