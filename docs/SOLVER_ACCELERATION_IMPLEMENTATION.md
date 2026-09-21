# Solver acceleration: implementation and measured gates

Follow-up: the [Stage-3 compiled prototype](SOLVER_ACCELERATION_STAGE3.md)
records the subsequent Numba experiment; the results below remain the Python
cleanup and adaptive-wall evidence.

Date: 2026-09-20. Implements the Python cleanup proposed in
[SOLVER_ACCELERATION_STAGE2.md](SOLVER_ACCELERATION_STAGE2.md), followed by a
separately authorized experimental adaptive wall initial guess.

## Outcome and adoption decision

The default Python path is faster with **exactly unchanged completed-cycle
trajectories, states and physical reports** on the frozen workloads. No ODE
method, physical equation, property correlation, motion law, or tolerance was
changed. No Numba or Julia backend was introduced.

Adaptive wall extrapolation is available **only by explicit opt-in**. It reduces
cycles in some trials, but the advantage is not universal at tighter periodic
accuracy. It is not promoted to the campaign default. Both positive and negative
results are retained below. The existing fixed cycle-10 extrapolation remains
controlled by `wall_numerical.accelerate_walls` when adaptive mode is absent.

The subsequently supplied [local study export](<Thermodynamic and Mechanical Study — dada-engine.org.html>)
(revision 1656, last modified 2026-09-20 04:27) has now been read. Sections 1,
5–7 and Appendix A allow either upstream or downstream placement independently
for the two check valves, while preserving branch circulation directions.
Placement changes which cylinder remains connected to each exchanger when the
valve closes and which local pressures govern opening. It therefore affects
physical results, not just nomenclature.

The four-volume mass/energy equations retain their general form. The current
implementation nevertheless hardcodes downstream valves and directed enthalpy
transport on those links. Supporting alternative placements will require
configurable hydraulic links, corresponding valve events and upwind enthalpy
handling on whichever links become bidirectional, plus consistent diagnostics.
That is a separate implementation task; this acceleration series retains the
frozen reference topology. The study also explicitly keeps heat/mass exchange
active through compression/expansion and separates valve events from kinematic
stage boundaries, consistent with the existing full-cycle approach. Its basic
reservoir-UA closure does not replace the current finite-air/wall model.

## Frozen baseline and exact cleanup results

All baseline and new measurements use Debian 13.7, Python 3.13.5, NumPy 2.2.4,
and SciPy 1.15.3. Julia 1.13.0 and Numba 0.61.2 are installed but unused.
The unchanged checkout passed 364 tests. The final checkout passes **430 tests**;
the 68 existing `np.trapz` deprecation warnings are unrelated to this work.
See [test output](../outputs/solver_acceleration_stage2/completed_pytest.txt).

The [frozen manifest](../outputs/solver_acceleration_stage2/manifest.json)
records exact parameters, initial states, tolerances, input hashes, source commit
and runtime identity. The production candidate comes from the Stage-1 artifact,
not a mutable current champion. Saved trajectories and full reports support
rechecking equivalence. Original campaign histories and source artifacts are
unchanged; new results live under `outputs/solver_acceleration_stage2/`.

The following are three serial repetitions, reporting median and observed
maximum wall time. The maximum is not a p99 estimate. Baseline and P6 use the
same machine/runtime; no gain is inferred from earlier Debian-12 timings.

| Frozen workload | Cycles, unchanged | Baseline median / max (s) | Clean Python median / max (s) | Median speedup |
| --- | ---: | ---: | ---: | ---: |
| Production warm | 1 | 17.544 / 17.676 | 6.387 / 6.570 | 2.75× |
| Nearby compatible warm | 3 | 29.855 / 30.120 | 12.122 / 12.123 | 2.46× |
| Production cold | 22 | 126.889 / 127.224 | 68.913 / 72.470 | 1.84× |
| Smooth four-bar | 32 | 170.978 / 196.143 | 85.760 / 87.082 | 1.99× |
| Historical regression | 1 | 2.647 / 2.698 | 1.766 / 1.804 | 1.50× |

Sources: [baseline](../outputs/solver_acceleration_stage2/p0_baseline/summary.json),
[P6](../outputs/solver_acceleration_stage2/p6_clean_python/summary.json),
[exact comparison](../outputs/solver_acceleration_stage2/p6_equivalence.json).
The invalid transport-domain case retains its exception and message; its
sub-millisecond timing is not a speedup claim. Controlled interruption checks
status and the retained physical endpoint, rather than equal elapsed time or
cycle count. P6 timings precede the boundary callback correction described below.
The [final exact replay](../outputs/solver_acceleration_stage2/final_equivalence.json)
verifies the default after adding that correction and the opt-in adaptive branch.
It overlapped the test suite for part of its run, so its timings are excluded
from performance claims. The [final reference check](../outputs/solver_acceleration_stage2/final_reference.json)
also repeats all 379 point comparisons.

The warm median decreased through 13.105 s (one microtube report), 9.145 s
(one hydraulic evaluation), 8.030 s (prepared constants), and 6.615 s (typed
diagnostic aggregation). These are whole-evaluation measurements. The separate
2048-call RHS microbenchmark improved from 1.082 to 0.949 s when preparing
constants; it is supporting evidence, not a substitute for complete evaluations.

## Implementation boundaries

- `ThermodynamicModel.instantaneous_point()` computes geometry, temperatures,
  pressures, topology and signed flows once. `assemble_rates()` is the sole
  conservative balance implementation; `evaluate()` remains the compatibility
  facade. Read-only primitive arrays are local to one trial state.
- The air-wall derivative supplies actual gas-wall heat rates to that assembler,
  removing repeated hydraulics and per-call model copying. Closed-valve tube
  contexts preserve equal passage pressures and the separate port pressure ratio.
- Immutable candidate-owned geometry, areas, external-air conductance and flow
  caps are prepared at construction. Replacement reconstructs them; serialized
  geometry fields and public independent report dictionaries remain compatible.
  Variable gas properties, upstream/bulk temperature distinctions and domain
  checks remain evaluation-local.
- Microtube diagnostics are calculated once per evaluation and aggregated from
  typed records, avoiding recursive per-sample `asdict()`. All samples, extrema,
  physical-time weights, absolute-heat weights and failure criteria are retained.
- The pre-refactor implementation exists only as a test oracle. Reference checks
  cover 379 sampled points from five saved trajectories, plus focused valve,
  reversal, zero-flow, near-equality, breakpoint and invalid-domain cases.

No approximate state cache, rounded temperature cache, property lookup table,
search algorithm change, independent UA optimization, or mechanical-efficiency
assumption is included.

## Adaptive wall guesses: algorithm and limitations

Four retained cycle endpoints supply three temperature increments and two
successive contraction ratios for each wall. A proposal requires significant,
same-sign contracting increments and sufficiently stable ratios. The estimated
remaining temperature change is damped according to ratio drift. Each correction
is bounded by both 30 K and 20 times the latest increment, so slowing or bending
trends limit the jump. Defaults require confidence at least 0.95 and contraction
ratios below 0.98; increments below 1e-7 K are ignored.

Only the two wall-energy entries change in the next **initial guess**. All eight
gas state entries, including all masses, are copied unchanged. No physical state
is modified during a cycle. The guessed wall energy is not counted as physical
heat input: heat/work and conservation refer to the following ordinary trajectory.

Every proposed guess must survive ordinary integration and the original
all-state periodic criterion. A plausible gas adjustment transient is provisional
until another ordinary cycle contracts sufficiently. Excessive residual growth
or a trial domain/numerical error causes rollback to the last retained endpoint,
followed by at least four ordinary recovery cycles. Interrupted or exhausted
runs return a retained completed trajectory, never an unverified guess.
Trials consume the cycle-attempt budget; rejected completed trials remain marked
in history, and failed integrations/rollbacks remain in solver statistics.
Only retained endpoints notify `cycle_callback`.

Single-run measurements at the original periodic tolerances:

| Workload | Fixed schedule cycles | Adaptive cycles | Adaptive time (s) |
| --- | ---: | ---: | ---: |
| Production cold | 22 | 14 | 44.371 |
| Smooth four-bar | 32 | 25 | 65.832 |

These are promising stopping-time results, **not proof of equal final output
accuracy**. The early less selective experiment actually increased cold-start
cycles from 22 to 26; its artifacts remain in `adaptive_trial/`. The retained
settings are recorded in `adaptive_stable_trial/environment.json`.

A separate common manifest tightened both periodic tolerances 100-fold, while
keeping the ODE settings and starting states unchanged:

| Workload | Fixed cycles / time (s) | Adaptive cycles / time (s) |
| --- | ---: | ---: |
| Production cold | 33 / 107.659 | 35 / 115.457 |
| Smooth four-bar | 52 / 128.449 | 35 / 100.311 |

These too are one run per method, not medians. At the original state-comparison
scale, the final state distances are 0.142 and 0.030 (threshold 1). Relative
indicated-power differences are 4.16e-7 and 2.51e-7. Relative indicated-efficiency
differences are 1.36e-6 and 6.46e-7: **the cold efficiency difference exceeds a
1e-6 output comparison threshold**, despite the tighter periodic residual.
Validity classifications agree. Relative total-mass differences between methods
are below 2.2e-14. See the reproducible
[accuracy summary](../outputs/solver_acceleration_stage2/adaptive_accuracy.json).

Thus periodic residuals are not output-error guarantees. The adaptive mode
passes conservation/safeguard tests and can save cycles, but does not meet a
universal speed/accuracy adoption gate. It remains experimental. Decision
arithmetic costs only a few milliseconds per full evaluation; extra physical
validation cycles, not proposal arithmetic, dominate its possible overhead.

Python use:

```python
from dada_solver.exchangers.wall_iteration import AdaptiveWallAccelerationSettings

result = solve_periodic_wall_motor(
    wrapper, initial_state, maximum_cycles=maximum_cycles, settings=settings,
    adaptive_acceleration=AdaptiveWallAccelerationSettings(),
)
```

For an explicitly separate microtube campaign, add the following TOML table:

```toml
[adaptive_wall_acceleration]
damping = 0.99
minimum_confidence = 0.95
```

Omitting the table preserves the previous numerical-settings schema and policy.
Including it replaces the fixed cycle-10 schedule, even if `accelerate_walls`
is true. Controls belong to exact campaign/cache identity; changing source or
controls requires a new campaign identity under the existing rules. Reservoir
campaigns reject this wall-only option. No persistent campaign was launched.

## Instrumentation, interruptions and next gate

Optional `statistics_callback` records per-segment method, `nfev`, `njev`, `nlu`,
actual RHS calls, accepted steps/sample counts, angular step extrema/median,
kinematic-breakpoint status, elapsed time and solver status. `measure_rhs_time`
adds an opt-in RHS timer. `EvaluationControl` also exposes section timing and
warm-start/status metadata as a side channel; persisted default campaign record
schemas are unchanged. Interrupted calls report unavailable SciPy counters as
`None`, never fabricated zeros.

On P6, the warm case takes 6080 RHS evaluations, 165 Jacobian evaluations /
factorizations and 2126 accepted steps. Cold takes 134000 / 3803 / 45730 over
22 cycles; smooth takes 154080 / 5872 / 42638 over 32 cycles. Integration remains
the main cold/smooth cost. Warm evaluation retains substantial diagnostic cost
(about 2.83 s integration versus 3.55 s diagnostics/validity in one replay).
An explicitly timed warm replay records 3.189 s inside the RHS out of 3.310 s
inside solver segments (96.3%, excluding preflight/reporting). Its physical
results remain exact; per-call timer overhead means these times are not part of
the speedup table. See
[instrumented replay](../outputs/solver_acceleration_stage2/final_rhs_instrumented/production_warm_0.json).
This supports a later targeted RHS/ODE prototype and comparison; it does not yet
justify migration or a claimed compiled-backend speedup.

Two interruption fixes are separately regression-tested:

1. After fixed-schedule extrapolation, return the last actual physical endpoint
   on interruption rather than the extrapolated guess.
2. Check cooperative progress/deadlines at each segment boundary. Otherwise,
   acceleration can make every segment shorter than the callback interval and
   starve cancellation for many cycles. The controlled P6 interruption exposed
   this (~19 s); the final path interrupts before entering the first segment.

Checks within long segments still use the configured interval. Diagnostic
postprocessing has no cooperative deadline, so this is not an end-to-end latency
or hard real-time guarantee.

## Reproduction

Run from the source checkout. Measurement directories must be new; frozen
manifests and input hashes may not be overwritten silently.

```sh
PYTHONPATH=src python3 -m pytest
PYTHONPATH=src python3 examples/benchmark_solver_acceleration.py \
  --output outputs/solver_acceleration_stage2/new_default --repeats 3
PYTHONPATH=src python3 examples/compare_solver_acceleration.py \
  outputs/solver_acceleration_stage2/p0_baseline \
  outputs/solver_acceleration_stage2/new_default
PYTHONPATH=src python3 examples/check_solver_acceleration_reference.py \
  outputs/solver_acceleration_stage2/p0_baseline
PYTHONPATH=src python3 examples/benchmark_solver_acceleration.py \
  --output outputs/solver_acceleration_stage2/new_adaptive --repeats 3 \
  --cases production_cold smooth_four_bar --adaptive-wall
PYTHONPATH=src python3 examples/audit_wall_acceleration.py \
  outputs/solver_acceleration_stage2/tight_fixed \
  outputs/solver_acceleration_stage2/tight_adaptive \
  outputs/solver_acceleration_stage2/manifest.json
```

`--manifest outputs/solver_acceleration_stage2/manifest_tight_periodic.json`
reproduces the tighter study for either method. `--adaptive-controls FILE.json`
records alternative settings as a separate experiment. `--measure-rhs` enables
per-call timing; keep such timings separate from uninstrumented latency.


## Study export update, 2026-09-20 evening

The replacement [local export](<Thermodynamic and Mechanical Study — dada-engine.org.html>)
is revision 1659 (19:02), compared with the previously reviewed revision 1656
(04:27). Section 4 generalizes pressure-equalized pair reductions to connected
sets, including inlet/outlet enthalpy and local exchanger mass redistribution.
Section 6 expands the physical rationale for exchange and distinguishes spatial
pressure equalization from pressure constancy over time. The swept-volume versus
temperature ratio is a first-order sizing guide, not an imposed geometric rule.
Motor chronology is recalculated with reversed kinematics and unchanged physical
valve orientations, not obtained by reversing every hydraulic flow.

These changes clarify reference reductions and chronology; they do not require
replacing the full four-volume balances used by the acceleration benchmarks.
The existing solver retains heat exchange and passive pressure-driven valves
throughout all kinematic phases. Independent upstream/downstream placement for
the H_i and H_o check valves is now supported by the Python and production Numba
paths for `continuous_ideal_diode`; the valve CdA moves with the directed half-link.
The historical DD placement remains the default. Alternative placements remain
unsupported for `discrete_hysteretic`.
No new pressure-equalization or exact-isobaric constraint is introduced here.
