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
endpoints. Final counts and full benchmark gates are recorded below once the
separate cache experiment is validated.

The initial no-cache integration is measured before cache changes. Three Numba
repetitions cover warm, nearby warm, cold, smooth four-bar, historical fallback,
invalid-domain and controlled-interruption cases. A fresh Python replay checks
that extracting shared primitives preserves the historical physical results.
The Stage-3 pointwise bound remains `1e-11 + 2e-11 * abs(reference_rhs)`; final-state
and output comparisons retain the original Stage-3 scale and are not relaxed.

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

## Residual-cost measurements and next-stage decisions

Full-evaluation section timers are supplemented by bounded native-kernel and
Python-dispatch microbenchmarks. Per-RHS timers are opt-in and their results are
separate from uninstrumented latency. The native loop includes allocation/native
call overhead; subtracting its time from Python dispatch is an estimate, not an
exact partition. Final-cycle samples do not reproduce every cold transient state.

The report will separately identify candidate/model construction, Python
kinematics, compiled dispatch, approximate native kernel cost, Python marshalling,
LSODA/callback overhead, periodic overhead, generic and microtube diagnostics,
validity and serialization/report materialization. No acceptance-critical validity
check is skipped or deferred in this stage.

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
