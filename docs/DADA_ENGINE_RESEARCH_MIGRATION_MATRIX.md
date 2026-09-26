# Dada-Engine Research — local migration audit

Status: phase 0 review package, 27 September 2026. No Research implementation,
legacy migration, deletion, or cleanup has been performed. The working source
is the current checkout, not the GitHub tree named in the original proposal.
The audited HEAD is `7c7c7c69b9b2ecf17dfffc10b3867cad09f1fa02`; content hashes,
rather than HEAD alone, identify the files inspected.

## Inventory and evidence

The [CSV matrix](DADA_ENGINE_RESEARCH_MIGRATION_MATRIX.csv) has one row for every
Python/TOML example: **240 files, 187 Python and 53 TOML**, including two locally
untracked Python scripts. It supplies all nine requested migration columns,
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
Complete semantic review is concentrated on the first slice and its construction,
evaluation, feasibility, persistence and presentation dependencies. No row is
marked migrated or parity verified merely because its filename looks relevant.

Regenerate the static inventory from any directory:

```sh
python3 docs/research_audit/generate_inventory.py
```

The script uses its own checkout location. Generated review documents are excluded
from documentation-reference discovery, so the inventory does not cite itself.
This is developer audit tooling; it is not the proposed Research application.

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
