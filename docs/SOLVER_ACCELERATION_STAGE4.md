# Stage 4: production Numba backend, cache experiment and residual costs

Date: 2026-09-20. Follows the [Stage-3 prototype](SOLVER_ACCELERATION_STAGE3.md).
Python remains the default and the authoritative physical reference. This stage
introduces no Julia implementation, ODE-method change, kinematics compilation,
optimization campaign, or change to the convergence policy.

## Implemented production architecture

`WallBackendSettings` selects `python` or `numba`. A `WallRHS` instance belongs to
one candidate and is passed explicitly to the existing air-wall integration
loop. It persists across that candidate's physical cycles. No model method is
replaced and no process-global monkey-patching is used by the backend or the
Stage-4 benchmark runner. The Stage-3 monkey-patched runner remains a historical
research artifact, not the production entry point.

```python
from dada_solver.wall_backend import WallBackendSettings
from dada_solver.exchangers.wall_cycle import solve_periodic_wall_motor

result = solve_periodic_wall_motor(
    wrapper, initial_state, maximum_cycles=maximum_cycles, settings=settings,
    backend=WallBackendSettings(name="numba"),
)
print(result.backend_statistics)
```

A microtube campaign can opt in with:

```toml
[wall_backend]
name = "numba"
```

Omitting the table preserves the default Python policy and the previous identity
schema. An explicit table enters campaign identity together with backend source,
Python, NumPy, Numba, llvmlite/LLVM and CPU/compiler configuration. The existing
whole-project runtime/source identity remains in force. As before, changed source
or settings require a new campaign identity; no existing history is rewritten.
Reservoir campaigns reject a wall-backend table. `pip install .[numba]` is an
optional dependency choice, not a requirement for Python operation.

The optional runtime is imported only for an explicit Numba selection. Missing
Numba or disabled JIT produces a recorded Python fallback. Unsupported object
families use the bound authoritative Python derivative. Unsupported instantaneous
regimes fall back for that entire RHS evaluation, preserving public domain
validation and exception messages. Counters distinguish requested/actual backend,
model-level fallback reason and state-level fallback calls. Completed and
interrupted results expose a backend snapshot; side-channel records also retain
counters on failed integrations. Opt-in campaign results retain those records.

Progress/deadline callbacks retain the Python ordering: initial-state validation
and derivative preflight precede the first segment callback, followed by the
existing checks inside long segments and at subsequent segment boundaries.
An invalid initial state therefore raises the same physical error even when a
callback would request interruption. Compilation itself, like an individual
derivative call, is not asynchronously preemptible and can delay that first
check. This is cooperative cancellation, not a hard wall-clock guarantee.
Retained physical endpoints and interruption statuses keep the existing semantics.

## One source for numerical physics

`src/dada_solver/numerical_primitives.py` contains the ordinary Python numerical
functions and the compilable kernel. Public model classes call these functions
while keeping their validation, protocols and descriptive diagnostics. Numba
registers compilation support for the functions without replacing their Python
bindings. The shared locations cover:

- existing Sutherland/Shomate gas properties and the existing helium data interpolation;
- primitive temperature/pressure reconstruction and enthalpy;
- pressure-squared Poiseuille, header/valve minor losses and the sonic orifice cap;
- tube dimensionless numbers, correlation-domain flags and laminar heat transfer;
- local film conductance, series resistance and wall/air/gas heat rates;
- conservative transfer accumulation and piston work.

There is one kernel parameterized by prepared species codes and gas constants.
Air, nitrogen, argon and helium use their existing authoritative transport laws;
no new helium physics or approximate lookup table was introduced. The pre-existing
helium interpolation is retained as published in the Python transport model.
The current compiled subset remains continuum/laminar, with explicit fallback
for slip, transition, turbulence or other unsupported validity regimes.

Source/destination indices, one-way flags and exchanger-port associations are
prepared arrays rather than hardcoded branches inside the kernel. The adapter
still supports only the existing downstream-valve model. Alternative placement
remains a separate physical feature; the numerical architecture does not require
a parallel kernel for it.

The Python object layer remains authoritative for unsupported models and errors.
The object and compiled paths still have different orchestration and validation
boundaries; zero duplicate lines was not pursued at the expense of generic model
interfaces. All shared compiled arithmetic remains float64 without fastmath,
parallel loops, GPU operations, state interpolation or approximate caching.
LSODA, tolerances and exact kinematic breakpoints are unchanged.

## Validation and frozen benchmarks

Measurements are stored under `outputs/solver_acceleration_stage4/`. The original
Stage-2 manifest and source-input hashes are reused. Timing runs are serial and
never overlap the test suite. Each run records source hashes, environment,
backend identity, final states/trajectories, reports and solver statistics.

The initial shared-primitives refactor passed 81 focused tests and the existing
455-test suite. The first production adapter added 47 checks, including all four
species, boundary temperatures, conservation, reversal/zero flow, breakpoint
neighbors, fallback errors, absent optional dependencies, identity and retained
endpoints. The final suite passes **506 tests** (68 existing NumPy `trapz`
deprecation warnings). Missing-Numba behavior is checked in a fresh subprocess.
The final control replay also verifies invalid-state precedence over interruption.

The initial no-cache integration is measured before cache changes. Three Numba
repetitions cover warm, nearby warm, cold, smooth four-bar, historical fallback,
invalid-domain and controlled-interruption cases. A fresh Python replay checks
that extracting shared primitives preserves the historical physical results.
The Stage-3 pointwise bound remains `1e-11 + 2e-11 * abs(reference_rhs)`; final-state
and output comparisons retain the original Stage-3 scale and are not relaxed.

The final Python replay reproduces the frozen Python reports and trajectories
exactly. The three final compiled replays also reproduce the Stage-3 compiled
reports and trajectories exactly. The 1,110 saved-state RHS comparisons pass;
the largest discrepancy consumes 0.01438 of the allowed pointwise bound, with
zero fallback on supported states. Full-run maximum relative differences from
Python are 3.63e-7 for power, 3.04e-7 for efficiency and 1.91e-14 for mass.
Maximum scaled final-state distance is 0.02182 (acceptance threshold 1).

| Frozen workload | Cycles | Python, fresh replay (s) | Numba, three-run median (s) | Numba maximum (s) |
|---|---:|---:|---:|---:|
| Production warm | 1 | 6.911 | 4.361 | 10.032 |
| Nearby compatible warm | 3 | 13.296 | 4.727 | 4.855 |
| Production cold | 22 | 72.002 | 8.925 | 8.960 |
| Smooth four-bar | 32 | 87.382 | 14.007 | 14.601 |
| Historical unsupported | 1 | 1.783 | 1.805 | 1.821 |

For historical before/after context, Stage-3 prototype medians were 4.554 s
warm, 4.375 s nearby, 8.040 s cold, 13.802 s smooth and 2.036 s fallback.
These were separate sessions, not paired timing trials; production integration
retains comparable throughput without claiming an additional speedup.

The first warm Numba run includes 5.534 s of first-call compilation. Subsequent
cases share the compiled structure. Python numbers are single fresh replays,
not three-run medians. Before the separate cache refactor, production Numba
medians were 4.882, 4.766, 8.228, 14.181 and 1.835 s respectively; this variation
is not evidence of a cache-off speed improvement or regression.

Supported workloads have no RHS fallback. The historical case explicitly uses
Python; the invalid-domain case invokes authoritative validation with the same
exception type/message. Controlled interruption retains zero completed cycles;
tests additionally interrupt after a retained cycle. A final control replay
(`control_confirmation`) validates preflight/progress ordering after the last
callback-order correction; physical computations in the timed replays are unchanged.

With periodic tolerances tightened by 100, cold/smooth need 33/52 cycles, exactly
as in the saved Python comparison. Relative power differences are 8.83e-8 and
2.54e-7; scaled final-state distances are 0.01064 and 0.04481. These pass the
**original** output comparison scale. ODE tolerances were not tightened, so this
is independent supporting evidence, not certification at 1e-8 state accuracy.

Artifacts: `final_python_exact.json`, `prototype_exact.json`,
`final_compiled_pointwise.json`, `final_comparison.json`, `tight_comparison.json`
and `final_pytest.txt` in the measurement directory.

## Separate disk-cache experiment

The cache is opt-in and is measured after the initial backend validation.
Numba documents limitations in cross-file and global-value invalidation; see
[Numba caching](https://numba.readthedocs.io/en/stable/developer/caching.html).
The production policy therefore uses a content-addressed copy of the complete
shared numerical module, with runtime/compiler/CPU identity in the namespace.
That copy is generated from the authoritative source, not a separately maintained
physical implementation. Numerical helpers and their compiled caller live in the
same module. Candidate parameters remain explicit arrays, never cached globals.

Fresh empty-cache, fresh populated-cache and warmed-process results must be
reported separately. Cache-enabled selection remains optional, including its
path and effective cache status in reproducibility metadata. No stale binary is
accepted merely because source size or modification time happens to match.

| Process/cache condition | First supported RHS (s) | Complete first warm evaluation (s) |
|---|---:|---:|
| Fresh process, cache disabled | 5.534 | 10.032 |
| Fresh process, empty disk cache | 6.183 | 13.240 |
| Fresh process, populated disk cache | 0.378 | 5.250 |
| Warm process, populated cache | 0.000048 | 4.805 |

Each fresh condition is one independent process, not a distribution of process
startups. The empty-cache run also spent unusually long in diagnostics; the full
13.240 s must not be attributed entirely to cache construction. Both cache runs
pass the physical comparison gate. A different nearby candidate needs only
about 47–50 microseconds on its first RHS after compilation/loading; candidate
values do not trigger specialization. Cache hit/miss counters are cumulative
per process. Cache namespace invalidation and unwritable-cache fallback are tested.
The cache remains **off by default**. First-call time in this refactor exceeds
the historical prototype's approximately 2.9 s and is recorded rather than hidden.

## Residual-cost measurements and next-stage decisions

Full-evaluation section timers are supplemented by bounded native-kernel and
Python-dispatch microbenchmarks. Per-RHS timers are opt-in and their results are
separate from uninstrumented latency. The native loop includes allocation/native
call overhead; subtracting its time from Python dispatch is an estimate, not an
exact partition. Final-cycle samples do not reproduce every cold transient state.

The following seconds come from the second instrumented replay, after compilation.
Indented conceptual components are presented as separate rows; native kernel and
crossing estimates partition dispatch and must not be added again to dispatch.

| Component (s) | Warm | Cold | Smooth four-bar |
|---|---:|---:|---:|
| Complete evaluation | 4.476 | 8.458 | 15.444 |
| Candidate/model construction | 0.00038 | 0.00038 | 0.03212 |
| Complete periodic solve | 0.2395 | 4.6076 | 12.8045 |
| Kinematics inside RHS | 0.0394 | 0.7488 | 7.7933 |
| Compiled dispatch, total | 0.0767 | 1.4571 | 2.0072 |
| Native kernel, estimated portion | 0.0261 | 0.5407 | 0.5012 |
| Python/compiled crossing, estimated portion | 0.0506 | 0.9164 | 1.5060 |
| Python RHS adapter | 0.0371 | 0.7051 | 0.8903 |
| Outer RHS/timer bookkeeping | 0.0106 | 0.2003 | 0.2678 |
| Integrator/callback/segment-report residual | 0.0732 | 1.4538 | 1.8078 |
| Periodic control outside segments/preflight/preparation | 0.00190 | 0.04161 | 0.03757 |
| Generic diagnostics | 1.7816 | 1.5865 | 1.1262 |
| Microtube diagnostics | 2.3228 | 2.1521 | 1.3617 |
| Validity | 0.1309 | 0.1109 | 0.1188 |
| Report materialization/other | 0.00087 | 0.00076 | 0.00076 |
| JSON serialization, outside evaluation timer | 0.00051 | 0.00175 | 0.00131 |

Preparation/preflight and all unrounded components are in `cost_summary.json`.
The residual is **not isolated LSODA internal work**. Timers perturb the execution;
use the uninstrumented table for end-to-end latency. Preflight and per-RHS timer boundaries overlap slightly; the residual table
is a measured attribution aid, not an exact accounting identity.

| Integration count | Warm | Cold | Smooth four-bar |
|---|---:|---:|---:|
| RHS calls including preflight | 7,060 | 134,565 | 153,650 |
| nfev | 7,059 | 134,543 | 153,618 |
| njev / nlu | 212 | 3,773 | 5,817 |
| Accepted steps | 2,255 | 46,222 | 42,713 |

**Recommendations, distinct from implemented changes:**

1. Prioritize shared final-trajectory reconstruction for warm evaluations:
   generic/microtube diagnostics and validity consume 94.6% of warm elapsed time.
   The current paths repeatedly reconstruct hydraulics and film rates. Preserve
   all acceptance-critical checks, samples, weights and reporting semantics.
2. Investigate exact kinematics reuse before compilation. Four-stage kinematics
   cost about 5.6 microseconds/call; smooth four-bar about 50.7. The latter accounts
   for 60.9% of integration and 50.5% of full evaluation. Geometry normalization
   is already prepared once, but adaptive-angle evaluations remain. Measure exact
   repeated-angle reuse and cheaper scalar arithmetic before introducing an
   optional compiled provider. Keep slider-crank, four-bar and six-bar interfaces
   generic; solenoids are deferred. This is not authorization to duplicate all
   mechanisms in compiled form.
3. A bounded matched-accuracy ODE comparison is justified after those measurements,
   especially for cold workloads. Collect nfev/njev/nlu, accepted steps and complete
   physical outputs; no default-method change follows from the present profile.
4. Defer Julia. Neither the 1.45 s cold nor 1.81 s smooth residual isolates a
   removable language/integrator cost. The larger measured opportunities are
   diagnostics and four-bar evaluation. Any later Julia trial must move the whole
   integration loop and establish end-to-end value.
5. Study the local mass-conserving periodic map separately before selecting a
   shooting accelerator. Production periodic convergence and experimental adaptive
   wall extrapolation remain unchanged.

No acceptance-critical validity check is skipped or deferred in this stage.

## Reproduction

```sh
PYTHONPATH=src python3 -m pytest
PYTHONPATH=src python3 examples/check_wall_backend.py
PYTHONPATH=src python3 examples/benchmark_solver_stage4.py \
  --backend numba --output outputs/solver_acceleration_stage4/new_numba --repeats 3
PYTHONPATH=src python3 examples/benchmark_solver_stage4.py \
  --backend python --output outputs/solver_acceleration_stage4/new_python --repeats 3
PYTHONPATH=src python3 examples/compare_numba_wall.py \
  outputs/solver_acceleration_stage4/new_python \
  outputs/solver_acceleration_stage4/new_numba \
  outputs/solver_acceleration_stage2/manifest.json
PYTHONPATH=src python3 examples/benchmark_solver_stage4.py \
  --backend numba --profile --cases production_warm production_cold smooth_four_bar \
  --output outputs/solver_acceleration_stage4/new_profile --repeats 2
PYTHONPATH=src python3 examples/profile_wall_backend.py \
  outputs/solver_acceleration_stage4/new_microprofile.json
```

Use the existing `manifest_tight_periodic.json` for the independent tighter
periodic experiment. Keep the original manifest as the comparison scale argument
when reproducing the stated Stage-3 output/state gate. Output directories must be
new. Cache experiments and instrumented runs are not pooled into ordinary timing
medians. Three repeats provide observed medians/maxima, not population tail estimates.

A subsequent [periodic-map and exact-kinematics reuse diagnostic](PERIODIC_MAP_DIAGNOSTIC.md)
measures the next questions separately from this production integration.
