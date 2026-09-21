# Stage 5: exact kinematics reuse, shared diagnostic replay and interpolation study

Completed 2026-09-21. This stage follows
[Stage 4](SOLVER_ACCELERATION_STAGE4.md) and the
[periodic-map diagnostic](PERIODIC_MAP_DIAGNOSTIC.md).
Implementation order was P1 exact reuse, its validation/benchmark, P2 shared
replay, its validation/benchmark, then P3 interpolation. The additional six-bar
measurements supplement the frozen four-stage/four-bar gate.

**Adopted:** candidate-local exact-angle reuse in wall periodic solves, and shared
final-trajectory diagnostic reconstruction in wall candidate evaluation.
**Not adopted:** approximate periodic interpolation. It passes the original
state/power/efficiency comparison scales in the tested runs but fails the broader
scientific-report comparison described below. It remains in experimental scripts.
No Anderson, Broyden, Newton, Julia, ODE-method change, new periodic policy,
thermodynamic approximation or physical-validity relaxation is included.

## P1: exact production cache

`kinematics_cache.ExactAngleKinematics` wraps an existing prepared combined
`cylinder_volumes_and_derivatives(angle)` provider. It stores one exact input
angle and the corresponding immutable four-value tuple. There is no rounding,
interpolation, thermodynamic state, RHS value or cross-candidate mutable cache.
Failed evaluations do not replace the last successful entry.

`solve_periodic_wall_motor` creates a fresh adapter and a private model/wrapper
copy once per solve. The cache persists over that solve's cycles, is attached
**after** the existing motor-direction transform, and never changes the supplied
model. Both Python and Numba wall paths use it. Providers without the optional
combined method are returned unchanged; the generic scalar methods delegate to
the original provider. Existing standalone integration entry points remain usable.
Prepared kinematics must be deterministic for fixed geometry and angle, as assumed
by the current prescribed-motion interface.

The cache is enabled by default for this exact wall-solver optimization.
`exact_kinematics_cache=False` remains available for comparison/debugging.
Calls, hits, misses and hit fraction appear in `backend_statistics` and Numba
backend progress records. No new campaign parameter or warm-start rule is added;
the existing project source identity accounts for the implementation change.

| Frozen solve | Combined calls | Exact hits | Hit fraction |
|---|---:|---:|---:|
| Production warm, four-stage | 7,060 | 4,656 | 65.95% |
| Nearby warm, four-stage | 20,289 | 13,405 | 66.07% |
| Production cold, four-stage | 134,565 | 85,640 | 63.64% |
| Smooth four-bar | 153,650 | 108,015 | 70.30% |
| Additional six-bar | 132,695 | 92,600 | 69.78% |

These are complete periodic-solve counts, including preflight calls, rather than
isolated-cycle hit rates. The exact-cache comparisons require identical complete
angle arrays and all 15 trajectory rows, final states and public result reports.
Both Python and Numba interruption tests retain precisely the same completed
endpoint. Adjacent representable angles are treated as distinct. No dedicated
slider-crank implementation exists in this checkout; adding one was outside this
acceleration stage. The available six-bar family was practical to evaluate.

## P2: shared final-trajectory replay

`diagnostic_replay.replay_wall_trajectory` builds evaluation-local typed facts:

- one authoritative `InstantaneousPoint` per stored sample, containing volumes,
  angular rates, temperatures, pressures, signed flows and effective topology;
- two `WallThermalPoint` records, retaining wall temperature, gas/source heat,
  wall energy rate and the film diagnostic objects already calculated;
- separate hydraulic upstream diagnostics for the four microtube passages.

The replay also exposes read-only pressure/temperature arrays. Explicit custom
heat-rate providers still override the diagnostic heat rates when supplied. It is bound to
one wrapper, cycle and trajectory, and rejects mismatched consumers. It is a
completed-trajectory reconstruction, not a cache used by the integrator.

The existing `extract_cycle_diagnostics`, `assess_cycle_validity` and
`cycle_microtube_diagnostics` accept an optional replay. Their standalone paths
remain available as the comparison reference. The wall campaign evaluator and
bounded evaluation scripts prepare one shared replay for a converged cycle and
pass it to all consumers. Legacy constant-film tube validity can also consume
its hydraulic point instead of reconstructing it again.

`AirWallExchanger.thermal_point` is the sole wall-rate implementation;
`rates()` materializes its original dictionary contract. The film calculation
and diagnostic objects are now retained together rather than recomputed for the
report. Custom wall implementations retain their public `rates` protocol and a
fallback diagnostic path. No diagnostic-specific copy of a property, flow or
film equation is introduced: these paths continue to call the authoritative
object APIs and Stage-4 numerical primitives.

Hydraulic and film diagnostics are deliberately distinct. The former use the
actual signed-flow upstream temperature and half-tube hydraulic length. Film
diagnostics use exchanger gas temperature and the original full thermal
entrance length. Their equality must not be assumed, especially at reverse or
zero flow. The original zero-flow pressure treatment is retained.

The numerical facts required for acceptance are separate from report aggregation:
consumers produce extrema, validity and microtube ranges/fractions from the shared
records. Detailed dictionaries are built at the reporting boundary. The original
trapezoidal time weighting, absolute-heat weighting and half-film heat split are
unchanged. No report-only calculation is deferred or omitted. Pressure,
temperature, flow, Mach, Reynolds, domain issues, topology and piston-force
inputs retain their contracts, including existing unavailable classifications.

The current wall diagnostic cycle has no explicitly localized valve events;
its existing empty event sequence/topology classification is preserved, rather
than inventing new events during this refactor. Effective instantaneous valve
states remain part of the shared point. The interpolation study separately
records pressure-based sampled valve transitions as an experimental diagnostic.

## Exact validation and feasibility

The final full source-checkout suite passes **520 tests**. The 580 warnings are
the existing `np.trapz` deprecation exercised more often by equivalence tests.
P1 separately passed 81 targeted tests; P2 passed 51. The replay test counts one
hydraulic reconstruction per sample and compares every report field exactly.

The following complete benchmark comparisons pass exact equality, without
numerical tolerances:

- P1 cache off versus on: every frozen case, two repetitions;
- P2 legacy versus shared replay: every frozen case, three repetitions;
- additional six-bar case at both steps, two repetitions;
- historical unsupported backend, invalid-domain exception/message, and
  controlled interruption cases remain unchanged.

Equality covers complete trajectories and heat/work quadratures, periodic state,
cycle counts, performance, pressure/temperature/flow/heat extrema, valve-event and
topology output, pressure regularization fraction, validity, all microtube
passage ranges, hydraulic upstream ranges, domain time fractions, heat fractions,
failed criteria and total absolute gas heat. Therefore no floating-point change
is accepted for these exact production steps.

A read-only configured assessment of the saved results additionally reproduces
all campaign constraints, objectives, conservation metrics and statuses exactly:
production warm, nearby warm, production cold and smooth four-bar remain feasible.
This check does not launch a campaign or write campaign histories/caches.
The saved tighter-periodic manifest also reproduces the Stage-4 complete reports,
final states and trajectories exactly: 33 production-cold and 52 smooth cycles
(`tight_exact.json`). The corresponding warmed complete evaluations took 7.423
and 10.748 s; these single runs are an accuracy check, not pooled timing medians.

Integration cancellation and deadline checkpoints remain in their existing order;
the shared replay uses the same coarse evaluation-control boundaries. Compilation
and individual scientific operations remain cooperatively, not asynchronously,
interruptible, as in Stage 4.

Artifacts under `outputs/solver_acceleration_stage5/` include `p1_exact.json`,
`p2_exact.json`, `p1_sixbar_exact.json`, `p2_sixbar_exact.json`,
`campaign_contract.json`, and `final_tests.txt`. A final replay after the
protocol/custom-heat-provider boundary checks also passes `final_exact.json`
(warm, smooth, invalid-domain and interruption cases); its single-run timings
are kept separate from the earlier medians.

## Exact benchmark setup and elapsed time

All runs use the production Numba backend. One supported RHS is compiled before
the timing loop; the excluded compilation time is recorded in each environment
file. Timed runs are serial and do not overlap tests. Source identities, frozen
input hashes, states, full trajectories and reports are saved. The original
Stage-2 manifest, LSODA settings, tolerances, exact breakpoints and fixed wall
extrapolation policy are retained.

The additional six-bar candidate uses the frozen production thermodynamic
parameters, existing independently synthesized small/large mechanisms, their
fixed branches and a uniform compatible initial charge. Its exact initial state
is frozen in `manifest_sixbar.json`. It converges in 27 cycles. The original four
workloads retain 1, 3, 22 and 32 cycles respectively.

| Complete evaluation, seconds | Cache off, old replay | Exact cache, old replay | Exact cache + shared replay | Overall observed ratio |
|---|---:|---:|---:|---:|
| Production warm | 4.632 | 4.167 | 1.687 | 2.75x |
| Nearby warm | 4.637 | 4.200 | 1.964 | 2.36x |
| Production cold | 10.129 | 8.263 | 6.046 | 1.68x |
| Smooth four-bar | 14.832 | 9.192 | 7.424 | 2.00x |
| Additional six-bar | 9.930 | 7.375 | 5.503 | 1.80x |

The first two columns use two repetitions, the last three for frozen cases and
two for six-bar. These are observed medians from successive batches, not paired
confidence intervals. Small differences can include runtime variability. In
particular, exact caching does not change the old diagnostic path, so changes
in its measured time cannot be credited to the cache.

| Periodic solve only, seconds | Cache off | Exact cache, P1 |
|---|---:|---:|
| Production warm | 0.269 | 0.198 |
| Nearby warm | 0.691 | 0.563 |
| Production cold | 5.551 | 4.022 |
| Smooth four-bar | 12.272 | 6.390 |
| Additional six-bar | 6.920 | 4.607 |

The strongest cache evidence is the large smooth-mechanism integration reduction
with identical grids and trajectories. Warm replay work decreases from about
4.36 s before these changes to 1.42 s with shared facts/aggregation: roughly 67%
of that cost is removed without dropping scientific information.

## Residual profile after accepted changes

The table uses medians of two separate instrumented replays. Per-RHS timers add
overhead, so use the preceding uninstrumented table for throughput. Kinematics
is a component of periodic solve, not an additional elapsed-time term.

| Component, seconds | Warm | Nearby | Cold | Four-bar |
|---|---:|---:|---:|---:|
| Periodic solve | 0.243 | 0.674 | 4.467 | 6.470 |
| Kinematics inside RHS, included above | 0.0234 | 0.0656 | 0.4480 | 2.3122 |
| Shared replay core | 1.4274 | 1.2480 | 1.2633 | 0.7346 |
| Generic diagnostic aggregation | 0.0172 | 0.0164 | 0.0162 | 0.0076 |
| Microtube report aggregation | 0.1621 | 0.1487 | 0.1516 | 0.0724 |
| Validity aggregation | 0.00044 | 0.00023 | 0.00034 | 0.00020 |
| Report materialization and unassigned overhead | 0.0123 | 0.0124 | 0.0124 | 0.0062 |
| JSON serialization, outside evaluation timer | 0.00095 | 0.00061 | 0.00180 | 0.00139 |
| Complete instrumented evaluation | 1.863 | 2.100 | 5.912 | 7.325 |

Microtube aggregation includes its final reporting dictionaries; its aggregation
and dictionary-construction costs are not independently isolated. The residual
materialization row includes result conversion and other small unassigned work.
Geometry inside the replay is included in its core timer, not in the RHS
kinematics timer. Exact unrounded values are in `cost_summary.json`.

Warm evaluations still spend most time reconstructing physical diagnostic facts.
The shared replay is ordinary Python calling shared authoritative primitives;
this stage adds no separately compiled diagnostic model. Compiling an adapter
that returns all required diagnostic facts could be a later bounded experiment,
but the present RHS kernel does not expose that complete record. Report-only
aggregation is now much smaller; deferring it offers limited gain and has not
been shown safe for the reporting contract. Nothing is deferred here.
For cold/smooth evaluations, periodic integration is again the dominant cost.

## P3: separate periodic interpolation experiment

`examples/experimental_periodic_interpolation.py` contains the experimental
provider; production modules and campaign configuration do not select it.
For each smooth exact mechanism, a periodic cubic spline is constructed from
both volume laws at uniform knots. Its angular derivatives come from those same
polynomials, not separately fitted derivative samples. Coefficients are generated
automatically by SciPy `CubicSpline`; scalar Horner evaluation avoids per-call
array allocation. Exact-cache reuse remains enabled around the resulting provider.

The experiment rejects declared breakpoints, nonperiodic endpoint data and
nonpositive interpolated extrema. The two selected mechanisms have fixed
assembly branches and undergo their ordinary exact closure checks. This does
not certify arbitrary mechanisms merely because they declare no breakpoints.
No interpolation is applied to the piecewise four-stage law. No new slider-crank
geometry or parallel family-specific implementation is introduced.

Grids of 128, 256, 512, 1024 and 2048 intervals are compared on an independent
half-offset grid of 32,771 angles, augmented by derivative-root extrema, endpoint
neighborhoods, high-curvature neighborhoods and four-bar public closure-margin
neighborhoods. The selected six-bar synthesis records have minimum primary/secondary
transmission sines 0.587/0.764 (small) and 0.661/0.753 (large), with minimum rod
axis cosines 0.954 and 0.950; these are not near-toggle examples. Six-bar
dense/high-curvature coverage is not a global proof of singularity clearance; no branch-near family beyond these fixed geometries is
claimed validated. Volume errors are normalized by swept volume; angular-rate
errors by the maximum exact angular rate, avoiding division by zero at extrema.
Absolute errors for all four components, worst angles and boundary errors are
saved in each `*_geometry.json`.

| Intervals | Four-bar max normalized volume error | Four-bar rate error | Six-bar volume error | Six-bar rate error |
|---|---:|---:|---:|---:|
| 128 | 1.66e-8 | 1.87e-6 | 2.16e-7 | 2.18e-5 |
| 256 | 1.04e-9 | 2.34e-7 | 1.34e-8 | 2.70e-6 |
| 512 | 6.49e-11 | 2.92e-8 | 8.32e-10 | 3.37e-7 |
| 1024 | 4.06e-12 | 3.65e-9 | 5.19e-11 | 4.21e-8 |
| 2048 | 2.55e-13 | 4.53e-10 | 3.23e-12 | 5.26e-9 |

Each entry is the maximum of the two cylinders. Periodic closure of interpolated
volume and derivative is exact at 0/2pi. At 2048 intervals, observed preparation
cost is about 0.107 s for four-bar and 0.059 s for six-bar.

### Full thermodynamic gate and rejection

All resolutions with maximum normalized geometric error below 1e-6 undergo a
complete periodic evaluation: four-bar at 256–2048 intervals and six-bar at
512–2048. Exact and interpolated runs share the same initial state and numerical
settings; compilation is warm. Interpolant construction is included in elapsed
time. Each experimental timing is a single run, not a stable speedup estimate.

| Mechanism / intervals | Complete time (s) | Final-state scaled distance | Largest relative power/efficiency/heat/mass difference | Full-report gate |
|---|---:|---:|---:|---|
| Four-bar exact | 7.670 | 0 | 0 | reference |
| Four-bar 256 | 5.573 | 0.03254 | 4.60e-7 | fail |
| Four-bar 512 | 5.255 | 0.01861 | 3.24e-7 | fail |
| Four-bar 1024 | 5.426 | 0.03175 | 5.09e-7 | fail |
| Four-bar 2048 | 5.388 | 0.02729 | 3.88e-7 | fail |
| Six-bar exact | 6.326 | 0 | 0 | reference |
| Six-bar 512 | 6.031 | 0.01862 | 1.55e-7 | fail |
| Six-bar 1024 | 5.045 | 0.01059 | 5.04e-7 | fail |
| Six-bar 2048 | 4.957 | 0.01415 | 4.36e-7 | fail |

The original Stage-3/4 state/output scales are 1e-6 relative and 1e-12 absolute;
scaled state acceptance is <=1. Every listed run passes that state gate and the
aggregate output gate. All preserve 32 four-bar / 27 six-bar cycles and physical
validity/domain classifications. Relative conservation residual comparisons also
pass the same absolute/relative scale; near-zero absolute energy residuals are
recorded rather than judged by relative error about zero.

For this stage's broader reporting contract, the same scale is additionally
applied field by field to extrema, microtube ranges, fractions and experimental
sampled valve diagnostics. Those checks fail. At 2048 intervals, 62 four-bar and
77 six-bar fields exceed that scale. For example, four-bar minimum Hi inlet
Reynolds changes from 0.000581 to 0.000260, and maximum residence time from
504.7 to 1127.0 s. Such extrema near stagnant flow are highly sensitive to the
adaptive sample locations. Heat-rate extrema and other ranges also differ.

These differences are not evidence that power or conservation is physically
wrong. The tiny geometric errors alter the adaptive integration grid; sampled
near-zero-flow quantities then have unstable extrema. This interpretation is
consistent with the saved unequal angle arrays and the lack of monotone output
improvement as interpolation density increases. Geometry converges strongly;
thermodynamic report equality does not. This stage does not redefine reporting
samples, remove sensitive fields or loosen a gate to admit the approximation.
Pressure-based valve opening fractions and estimated crossing angles are stored
separately; each valve still has two detected transitions in these runs.

**Decision:** interpolation has measurable potential for the four-bar case, but
has not met the complete adoption gate. Keep only exact reuse in production.
The result does not justify declaring periodic interpolation universally useless;
it establishes that it is not currently a transparent replacement under the
requested scientific-output contract. No approximate campaign path is enabled.

## Reproduction and next boundary

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m pytest -q
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 examples/benchmark_solver_stage5.py \
  --output outputs/solver_acceleration_stage5/new_baseline --no-exact-cache --legacy-replay
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 examples/benchmark_solver_stage5.py \
  --output outputs/solver_acceleration_stage5/new_cache --legacy-replay
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 examples/benchmark_solver_stage5.py \
  --output outputs/solver_acceleration_stage5/new_shared
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 examples/compare_stage5_exact.py \
  outputs/solver_acceleration_stage5/new_baseline outputs/solver_acceleration_stage5/new_shared
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 examples/study_stage5_interpolation.py \
  --output outputs/solver_acceleration_stage5/new_interpolation
```

Use `--profile` for separate instrumented runs and `--manifest
outputs/solver_acceleration_stage5/manifest_sixbar.json` for the additional
mechanism. Output directories must be new. `check_stage5_campaign_contract.py`
checks the saved feasibility contract without running an optimization campaign.
`summarize_solver_stage5.py` emits the exact timing breakdown.
The subsequent Anderson experiment remains a separate stage and report, using
these exact production changes as its numerical baseline.
