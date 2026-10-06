# Numerical execution and performance

This reference consolidates the September 2026 acceleration reviews and
experiments. It describes implementation boundaries, not a new numerical
policy. Scientific inputs, tolerances and convergence criteria remain owned by
the study. See [working-fluid backends](WORKING_FLUID_MODELS.md) for tabulated
properties and [campaign internals](OPTIMIZATION_CAMPAIGN.md) for orchestration.

## What is prepared, reused and compiled

| Work | Lifetime / implementation |
|---|---|
| Geometry, tube sections, collector volumes, wall capacities | Prepared once per constructed candidate object |
| Instantaneous pressures, temperatures, flows and film context | Computed for the current trial state; no approximate state cache |
| Conservative mass/energy assembly | Shared physical balance implementation |
| Repeated kinematic angle | Candidate-local exact last-angle cache; no rounded key/interpolation |
| Supported wall RHS | Optional Numba kernel using shared numerical primitives |
| Final-trajectory diagnostics | One shared authoritative replay, then aggregation |
| Exact candidate cache | Persistent identity-based reuse; not a nearby-geometry approximation |

`ThermodynamicModel.instantaneous_point()` separates reconstruction and hydraulic
facts from conservative rate assembly. Wall heat enters that same balance.
Prepared geometry is reconstructed when its immutable inputs change; derived
values are not recomputed inside every RHS. Diagnostic replay distinguishes
upstream hydraulic temperature/half-tube length from bulk film temperature/full
thermal entry length. Those contexts are not interchangeable.

`kinematics_cache.ExactAngleKinematics` caches the four combined volume/rate
values for one exactly equal angle, after the motor-direction transform. Adjacent
floating-point angles remain different. It stores no gas state or RHS result;
failed evaluations do not poison the last successful entry. The wall solve uses
it by default; the direct solver can disable it for equivalence checks.

## Backend selection

The lower-level wall API accepts `WallBackendSettings(name="numba")`; Research
exposes backend settings in its [reference](DADA_ENGINE_RESEARCH_REFERENCE.md).
Python remains the reference implementation; optional Numba imports are not a
requirement for Python operation. Use the **actual backend and fallback counters**
in the result, not merely the requested setting, when interpreting timings.

The compiled path receives prepared primitive arrays, uses float64 without
fastmath and calls no Python property library for a supported state. Unsupported
families/states use explicit authoritative Python fallback. Valid no-slip
microtube states use the compiled laminar, transition and turbulent closures.
The existing transition endpoint interpolation is shared with Python; a bracketed
flow solve runs inside Numba. Bennett heat transfer and cumulative Shah
entrance loss also run in the shared compiled primitives; laminar flow now
needs a scalar root solve instead of only the former analytic quadratic.
No regime is relabelled or silently rejected to
keep a benchmark fast. Both valve placements are supported
by the prepared source/destination and one-way flags. The ideal and tabulated
validation backends are separate dispatch paths.

The optional ideal-kernel disk cache uses a content/runtime/compiler/CPU namespace
and a source-derived numerical module, not a separately maintained physics copy.
Report cold compilation, populated-cache startup and warmed calls separately.
The table prototype uses in-process JIT only. Settings/source changes belong to
runtime identity; do not resume an incompatible campaign by suppressing a check.

## Periodic iteration is a separate cost

Wall initial-guess extrapolation acts between cycles and must be explicitly
selected where applicable. It never modifies a physical trajectory or credits
heat/work from an initial-guess jump. Only a normal cycle satisfying the original
periodic test can certify a solution. Compatible warm starts reduce work without
relaxing that test.

Deadlines are cooperative. A JIT compilation or individual scientific operation
can delay the next callback. Interrupted records retain completed physical-cycle
endpoints; partial cycles and numerical guesses are not periodic solutions.
LSODA's native failure text can be captured while Python exceptions and rejection
status remain visible in technical diagnostics.

## Frozen evidence, not current speed guarantees

Measurements below are historical local observations on matched frozen workloads,
not universal performance claims or timings for the current transport revision.
The input/runtime/source manifests remain under `outputs/solver_acceleration_stage2/`.
Do not mix instrumented timings with ordinary latency or add nested profile times.

| Complete evaluation (seconds) | Original Python median | Clean Python median | Production Numba median | Numba + exact angle cache + shared replay |
|---|---:|---:|---:|---:|
| Warm | 17.544 | 6.387 | 4.361 | 1.687 |
| Nearby compatible warm | 29.855 | 12.122 | 4.727 | 1.964 |
| Cold | 126.889 | 68.913 | 8.925 | 6.046 |
| Smooth four-bar | 170.978 | 85.760 | 14.007 | 7.424 |

Columns come from successive recorded batches, not paired confidence intervals.
Python cleanup and exact reuse/shared replay were checked against complete
trajectories and reports. Compiled-vs-Python comparison used the original stated
numerical tolerances: Stage 4's largest reported relative power/efficiency
changes were 3.63e-7 / 3.04e-7; supported states had no fallback. These measurements
precede the later dilute-transport update and do not certify an old fixture
against changed physics.

Useful artifacts (local-only output directories may be absent elsewhere):

| Evidence | Location |
|---|---|
| Initial profile and exact warm state | `outputs/solver_performance_stage1_20260917.json` |
| Frozen inputs, Python cleanup and adaptive guesses | `outputs/solver_acceleration_stage2/manifest.json`, `p6_equivalence.json`, `final_equivalence.json` |
| Isolated Numba prototype | `outputs/solver_acceleration_stage3/` |
| Production adapter/parity/cache | `outputs/solver_acceleration_stage4/final_comparison.json`, `final_compiled_pointwise.json`, `tight_comparison.json` |
| Exact cache/shared replay | `outputs/solver_acceleration_stage5/p1_exact.json`, `p2_exact.json`, `campaign_contract.json`, `tight_exact.json` |

## Negative experiments worth retaining

**Periodic-map diagnostic.** A ten-state map with fixed total gas mass was
projected onto nine scaled coordinates: `x = x* + S B z`, where B spans the
scaled mass-row null space. Central differences found two slow modes near
0.950 and 0.932–0.933 in the tested cold/four-bar cases. Modes coupled gas and wall
temperatures; they were not two independent wall-only states. Some perturbations
failed the then-unsupported transition domain. Failed columns were not filled
or regularized. This is local evidence, not a global convergence guarantee.
Artifacts: `outputs/periodic_map_diagnostic/{production_tolerances,tighter_ode}/`.

**Approximate periodic interpolation.** Cubic volume/rate interpolation converged
geometrically and passed aggregate state/power comparisons, but failed the full
report gate. Small adaptive-grid shifts changed extrema near stagnant flow.
At 2048 intervals, 62 four-bar and 77 six-bar report fields exceeded the original
comparison scale. Production retained exact reuse only. Experimental scripts:
`experimental_periodic_interpolation.py`, `study_stage5_interpolation.py` under
`examples/`; evidence under `outputs/solver_acceleration_stage5/`.

**Bounded Anderson acceleration.** The direct experimental API uses mass-conserving
scaled coordinates, bounded least-squares coefficients, positivity checks and
rollback. It remains opt-in and has no Research campaign selection. All nine cold
screens exhausted 100 cycles versus 22 for the reference policy. The fastest
converged four-bar screen took 59 versus 32 cycles and failed the separate output
accuracy gate. No common setting was adopted. A coefficient L1 cap of 10 excluded
the cancellation coefficients (~39) needed for a scalar 0.95 contraction; that
explains stagnation but does not justify increasing the cap without a new study.
Evidence: `outputs/periodic_anderson/`, including `screen/`, `residuals.csv`,
`default_comparison.json` and `interruption_comparison.json`.

No Julia migration, new default integrator or weakened scientific gate follows
from these experiments.

## Reproduce a bounded benchmark

Run from the checkout, with new output directories and no concurrent optimization
or test workload when measuring time:

```sh
PYTHONPATH=src python3 examples/check_wall_backend.py
PYTHONPATH=src python3 examples/benchmark_solver_stage4.py --backend numba --output outputs/my_backend_benchmark --repeats 3
PYTHONPATH=src python3 examples/benchmark_solver_stage5.py --output outputs/my_shared_replay_benchmark
PYTHONPATH=src python3 examples/benchmark_research_v3_backends.py --calls 10000 --output outputs/my_fluid_benchmark.json
```

These scripts use frozen input hashes and may require their local manifests.
For equivalence, compare matched settings and source versions; do not loosen a
fixture because a later physical model changed. Record build, integration,
periodic iterations, diagnostics, persistence, JIT and fallback costs separately.
Current targeted tests include `test_wall_backend.py`, `test_diagnostic_replay.py`,
`test_kinematics_cache.py` and the periodic-map/Anderson tests. The
[validation map](validation.md) gives the complete-suite command.

## Compiled transition benchmark

The October 2026 DD10 investigation found substantial Python fallback cost,
not repeated collector geometry calculation. A saved evaluation spent 46.98 s
in integration with 47,126 fallbacks out of 120,736 RHS calls. Different geometry
can change flow regimes and periodic convergence, so cross-campaign durations
alone do not isolate a geometry implementation regression.

The no-slip transition/turbulent hydraulic solve and thermal correlations now run
inside Numba, sharing formulas with the Python reference. The laminar algebraic
path remains unchanged. Invalid states retain authoritative Python diagnostics;
no Mach, pressure-drop, entry, transport or slip guard is relaxed.

A matched 128-state DD10 benchmark measured 45.19 ms for the previous kernel
versus 11.02 ms for the new one (median of three warmed repetitions, **4.10x**).
Fallbacks fell from 168/384 to zero. Python-only took 83.01 ms. RHS parity passed
`atol=1e-11, rtol=2e-11`; a separate single-cycle comparison remained within the
existing trajectory tolerance. These local timings are not universal campaign
speed claims. [Frozen inputs, results and reproduction command](../outputs/research_compiled_transition/README.md)
use `tools/benchmark_microtube_transition_backend.py`, without starting a search.


The 2026-10-04 developing-entry review found that retaining the old transition
endpoint beside the new laminar entrance loss caused numerical switching.
The authorized revision matches the complete laminar endpoints to the turbulent
endpoints over 2300–4000, without relaxing integration tolerances. Both branches
and the interpolation execute inside Numba. Bounded DD13 evidence, including
superseded diagnostic attempts and final converged replays, is kept under
`outputs/research_microtube_developing_entry/`; final results use its
`continuous_transition_v1/` subdirectory. Historical timings above are not a
benchmark of this changed physical model.
