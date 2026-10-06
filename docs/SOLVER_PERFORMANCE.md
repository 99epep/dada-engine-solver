# Numerical execution and performance

This reference describes numerical execution and performance measurement.
Scientific inputs, tolerances and convergence criteria remain owned by
the study. See [working-fluid backends](WORKING_FLUID_MODELS.md) for tabulated
properties and [campaign internals](OPTIMIZATION_CAMPAIGN.md) for orchestration.

## Prepared and reused work

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

Candidate-result reuse requires exact candidate identity within a compatible
scientific study and runtime definition. A nearby design or changed source is
not a cache hit; see [Research identity](DADA_ENGINE_RESEARCH_REFERENCE.md).

## Backend selection and fallback

The lower-level wall API accepts `WallBackendSettings(name="numba")`; Research
exposes backend settings in its [reference](DADA_ENGINE_RESEARCH_REFERENCE.md).
Python remains the reference implementation; optional Numba imports are not a
requirement for Python operation. Use the **actual backend and fallback counters**
in the result, not merely the requested setting, when interpreting timings.

The compiled path receives prepared primitive arrays, uses float64 without
fastmath and calls no Python property library for a supported state. Unsupported
families/states use explicit authoritative Python fallback. Supported microtube closures use the same numerical primitives
in Python and
Numba, including scalar flow solves. No regime is relabelled or physical guard
relaxed to keep execution fast. Both valve placements use prepared
source/destination and one-way flags. See the
[microtube model](MICROTUBE_GAS_MODEL.md) for physical domains and the
[working-fluid reference](WORKING_FLUID_MODELS.md#compiled-execution) for the
separate ideal and tabulated dispatch paths.

The optional ideal-kernel disk cache uses a content/runtime/compiler/CPU namespace
and a source-derived numerical module, not a separately maintained physics copy.
Report cold compilation, populated-cache startup and warmed calls separately.
The table prototype uses in-process JIT only. Settings/source changes belong to
runtime identity; do not resume an incompatible campaign by suppressing a check.

## Periodic iteration and interruption

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

## Performance measurement

Measure in a quiet environment with matched physical inputs, initialization,
numerical settings and source versions. Distinguish cold compilation,
populated-cache startup and warmed execution; warm every exercised path before
steady timing. Record the actual backend, fallbacks and runtime identity.
Cross-candidate elapsed times also reflect flow regimes and periodic convergence,
so they do not isolate implementation speed.

Separate construction, integration, periodic iteration, diagnostics, persistence
and JIT costs where counters are available. Do not add nested profile times or
compare instrumented timings directly with ordinary latency. Check instantaneous
RHS parity, complete trajectories and sampled diagnostics separately: agreement
in aggregate power alone is insufficient. Never loosen tolerances or model-domain
checks to improve a benchmark.

The bounded benchmark is available from a source checkout:

- [Saved-state RHS benchmark](../tools/benchmark_microtube_transition_backend.py):
  accepts a compatible study, raw candidate record and trajectory containing
  `angles` and `trajectory`. It warms the exercised paths, checks Python parity
  and records backend/fallback counters without integration or optimization.
  Input compatibility is the caller's responsibility.

```sh
PYTHONPATH=src python3 tools/benchmark_microtube_transition_backend.py --help
```

Use a new output destination. The benchmark generates measurements; its outputs
are not normative speed guarantees. The [validation map](validation.md) owns
current test instructions. Focused execution checks include
`tests/test_wall_backend.py`, `tests/test_wall_backend_cache.py`,
`tests/test_kinematics_cache.py`, `tests/test_diagnostic_replay.py` and
`tests/test_solver_acceleration.py`.

## Rejected numerical shortcuts

Gas and wall storage are coupled in the periodic map. A local map analysis can
help diagnose slow convergence but cannot establish a global convergence guarantee.

Approximate periodic interpolation is not a substitute for exact kinematic reuse
based only on aggregate agreement. Changes in adaptive sampling can change extrema
and near-stagnant-flow diagnostics. Production retains exact angle reuse.

Anderson acceleration remains an opt-in experimental direct wall-cycle API, not
an adopted Research policy or a replacement for the ordinary convergence test.
No common setting is established for general production use. Initial-guess
experiments do not authorize a new default integrator, a language migration or
weaker scientific acceptance criteria.
