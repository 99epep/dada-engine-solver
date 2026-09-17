# Solver acceleration — stage 1 evidence and handoff

Date: 2026-09-17. Status: **first review, not an implementation decision**.
The user requests a bounded initial investigation, followed by an independent
review with Chat, before a later implementation pass. No production solver,
physical equation, convergence tolerance or optimizer was changed here. No
optimization campaign was launched. Three bounded saved-state replays were used
for timing, rather than a large benchmark or language-port experiment.

## Recommendation for the second review

First eliminate duplicated work in the present Python evaluator, measure again,
then decide whether a compiled numerical kernel is necessary. Keep Julia as a
candidate for a measured prototype, not an already-approved rewrite. There are
concrete avoidable costs before reaching that decision, particularly the repeated
final diagnostics and double thermo-hydraulic evaluation per RHS call.

Separate four objectives:

1. lower the cost of one right-hand-side evaluation;
2. lower the number of RHS evaluations per physical cycle without degrading accuracy;
3. lower the number of cycles needed for the unchanged periodic criterion;
4. lower diagnostics/orchestration cost and improve campaign throughput.

A useful cost decomposition is:

```text
candidate time = build + sum(cycle integration costs)
               + final diagnostics + serialization/persistence
cycle integration cost ~= RHS count * RHS cost + integrator overhead
```

A faster warm replay is not automatically a faster cold or nearby-candidate
solve. A faster compiled RHS cannot remove unnecessary Python post-processing.
Parallel throughput is not the same as single-candidate latency.

## Evidence collected on the current model

The source commit is `58e62b71b0b054222be6ad5ade386ca7a527c7dd`. The input is the
production variable-property four-stage campaign, not the earlier constant-UA
motor profile in [EXCHANGER_NEXT_STEPS.md](EXCHANGER_NEXT_STEPS.md).
The [measurement artifact](../outputs/solver_performance_stage1_20260917.json)
contains report/history hashes, candidate identity, its parameters, the exact
saved ten-state initial guess, numerical settings, runtime versions, probe
source code, timing results and cProfile tables. No campaign history was edited.

### Campaign costs already recorded

At this review the [production campaign history](../outputs/motor_four_stage_hx9d_variable_gas/history.jsonl)
contains **312 attempts**, a later checkpoint than the 216-attempt documentation
snapshot. Statistics below are derived from those recorded evaluations:

| Status | Count | Median elapsed (s) | 90th percentile (s) | Maximum (s) |
| --- | ---: | ---: | ---: | ---: |
| Converged | 256 | 108.4273 | 168.2619 | 189.8058 |
| Interrupted | 20 | 180.6643 | 181.5057 | 181.8960 |
| Invalid exchanger | 36 | 47.2542 | 95.8850 | 149.6690 |

Converged records have median 13 cycles, 90th percentile 22, maximum 24.
These are mixed workloads, not controlled benchmarks; interrupted durations
are censored by their deadlines. They do show that both per-cycle cost and
periodic iteration count matter. Some exchanger-domain rejections are expensive;
this is motivation for diagnosis, not permission to predict rejection from an
unjustified surrogate or to remove instantaneous validity checks.

The recorded best candidate (index 299) used 12 cycles. Its normalized periodic
error fell from 511.49 to 8.00 to 2.59, then only to 1.73 by cycle 10. The existing
wall initial-guess acceleration was applied at cycle 10; errors were 15.96 and
0.319 on cycles 11 and 12. This illustrates a slow late convergence mode but
does not identify its eigenstructure or establish the best acceleration method.

### Bounded replay and profiling

All replays start from that same candidate's saved final state and use the
existing numerical settings and `_evaluate` path. Each completes one cycle and
passes the original periodic test, with normalized error 0.09299627416662144.
This is a warm replay, not the cost of searching for a new periodic solution.
The unprofiled replay took **20.9601 s**; cProfile took **46.8491 s**.
Both produced exactly the same reported efficiency 0.22070634462673688 and
indicated power 52.40224307474156 W. The small difference from the source
candidate's stored performance is a further-cycle result within the existing
periodic criterion, not an optimization gain.

The instrumented call tree reports:

| Operation | Calls | Cumulative profiled time (s) |
| --- | ---: | ---: |
| Full evaluation | 1 | 46.849 |
| Periodic cycle integration | 1 | 14.844 |
| Microtube cycle diagnostics | 2 | 26.874 |
| Conservative model `evaluate` (including diagnostics) | 24942 | 19.991 |
| Instantaneous microtube `diagnose` | 159894 | 12.854 |
| `dataclasses.asdict` | 34082 | 11.519 |

These are **nested cumulative times**: never add their rows. cProfile more than
doubles elapsed time on this workload; its percentages cannot be directly
converted into an unprofiled speedup promise. There are 6080 integrator RHS
calls, four `solve_ivp` segments and 2126 accepted steps. The preflight RHS call
adds one. No Jacobian/factorization counters are retained by the wall wrapper;
this profile alone cannot diagnose stiffness or its origin.

A third replay uses only coarse function-boundary timers to distinguish
integration and post-processing without per-call profiling. Its measured
breakdown is recorded below and in `stage_timing` in the artifact.

| Section | Elapsed (s) | Share of this replay |
| --- | ---: | ---: |
| Periodic integration | 6.6184 | 32.6% |
| Generic diagnostics | 2.5872 | 12.7% |
| Generic validity | 0.1248 | 0.6% |
| Microtube diagnostics, first pass | 5.6668 | 27.9% |
| Microtube diagnostics, repeated pass | 5.3051 | 26.1% |
| Other evaluation overhead | 0.0037 | 0.0% |
| Total | 20.3060 | 100% |

The second microtube pass alone takes about 5.3 s (26% of this warm replay).
Eliminating exactly that work, if all other costs stayed unchanged, would imply
about a 1.35x end-to-end speedup **for this case only**. This is an arithmetic
opportunity estimate, not a measured optimized implementation. The combined
diagnostics dominate this one-cycle case; periodic integration will dominate
more strongly when many cycles are necessary.


Only one unprofiled run and one coarse-timed repeat were made; machine load,
cache effects and variance have not been characterized. No hardware-independent
speedup claim is made. The probe deadlines apply to integration callbacks;
post-processing currently has no equivalent cooperative deadline. A later
implementation must include that cost when promising end-to-end phase budgets.

## Concrete code findings and proposed order

### 1. Share final diagnostic work — first implementation candidate

In [the current evaluation helper](../examples/compare_motor_motion_laws_stage7A5.py),
`_tube_validity()` calls `cycle_microtube_diagnostics()` to obtain extrema, then
`_evaluate()` calls the same complete diagnostic routine again to retain domain
fractions. The [campaign evaluator](../src/dada_solver/campaign/evaluator.py)
has the same pattern. This is verified by two calls in the profile.

Compute that report once and use it for both extrema and the final domain
verdict. Next share the reconstructed state/flows/thermal rates among generic
and microtube diagnostics. Within
[gas_diagnostics.py](../src/dada_solver/exchangers/gas_diagnostics.py), contexts,
thermal rates and gas-film diagnostics currently cause further repeated work.
`asdict` recursively materializes diagnostic dictionaries for every sampled
port; aggregate from typed/scalar values and serialize at the boundary instead.

Retain all extrema, failure reasons, physical-time-weighted domain fractions
and absolute-heat fractions. Do not decimate trajectories or remove safety
checks merely to obtain a speedup. In the first patch, reusing the identical
report is simpler and lower risk than redesigning diagnostics wholesale.
Its gain is especially relevant to one-cycle warm replays; for a many-cycle
candidate, integration occupies a larger share. Re-measure both cases.

### 2. Evaluate the instantaneous physical state and flows once per RHS

[AirWallMotor.derivative](../src/dada_solver/exchangers/air_wall.py) currently:

1. obtains thermal rates; production flow contexts call `model.evaluate`;
2. builds a replacement model with prescribed wall-to-gas heat;
3. calls `evaluate` again for the conservative state derivative and work.

Both evaluations recompute hydraulics for the same angle and conservative state.
State reconstruction, temperature/pressure calculation and some kinematics are
also repeated. `dataclasses.replace` constructs a new model per RHS call.

Proposed boundary: a small instantaneous evaluation result with volumes/rates,
thermodynamic state and the four signed flows. Use it first to calculate the
contextual gas film, then assemble conservative mass/energy rates using those
flows and the actual gas-wall heat. The generic thermodynamic layer must remain
independent of concrete kinematics and exchanger families. Keep one implementation
of the balances; do not create a divergent 'fast physics' copy.

**Do not memoize by angle alone.** Adaptive integration revisits angles with
different trial states, including rejected steps. Reuse within one evaluation
or from an exact completed sample, not an approximate state lookup.

### 3. Prepare immutable candidate constants; reuse properties locally

The hot paths repeatedly calculate tube geometry, packing dimensions, flow
areas, resistance prefactors, species constants and external-air effective
conductance. These are invariant within one candidate. Prepare them once at
construction, with clear invalidation when geometry/material/air inputs change.
Do not cache values across candidates by an incomplete geometry key.

Variable transport must stay variable. Within one instantaneous evaluation,
reuse exact mu/k/cp values at the same temperature rather than recomputing them
inside several diagnostics. Do not round temperature or introduce property
interpolation tables in this first optimization pass. Distinguish hydraulic
upstream temperature from exchanger bulk temperature; they are not interchangeable.
Separate reporting-object allocation from the algebra used for film, hydraulic
closure and mandatory instantaneous domain checks.

### 4. Add cost counters before changing numerical algorithms

The eight-state integrator already records `nfev`, `njev`, `nlu`; the wall path
currently discards the corresponding `solve_ivp` counters. Retain them per
segment/cycle, plus accepted samples, integration time, diagnostics time,
serialization time, warm-start source and periodic error history. Include costs
for failed and interrupted evaluations, not only successful ones.

The current maximum step is pi/360 radians: at least roughly 720 steps per full
cycle even before adaptivity adds steps. Four-stage slope jumps split integration
into four solves; preserve exact breakpoints. The measured 2126 steps show that
this step cap alone does not explain the workload. Assess counter data before
trying LSODA/Radau/BDF comparisons, Jacobians, sparsity or step-limit changes.
Keep these experiments separate from algebra-preserving refactoring.

### 5. Periodic convergence and compatible warm starts

Warm starts and bounded wall Aitken acceleration are already implemented, not
new proposed features. Assess their real cycle savings and failure modes before
adding complexity. Nearness in normalized design coordinates is not proof of
nearness of the periodic state. Preserve inventory rescaling and wall-temperature
preservation when capacities change.

A possible later experiment is earlier/adaptive bounded extrapolation or a
safeguarded fixed-point accelerator. It must preserve positivity, inventory and
compatibility, allow fallback, and be followed by an **ordinary unmodified cycle**
that passes the original periodic criterion. Do not accept the extrapolated
state as a converged solution or lower maximum cycles to disguise poor convergence.
Keep a last physically completed-cycle state distinct from any extrapolated
initial guess, especially when interruption occurs after extrapolation.

Direct shooting/Newton methods could eventually reduce cycle count, but require
conditioning, mass-invariant and non-smooth valve/kinematic analysis. They are
not the first low-risk acceleration step.

## Python kernel, Julia, or parallel evaluation?

| Option | Proposed role | Main qualification |
| --- | --- | --- |
| Remove duplicate Python work | First changes after review | Measured call duplication; preserve APIs and equations |
| Compiled numerical kernel within Python | Small prototype if RHS still dominates | Compile coherent scalar algebra, not many tiny crossings between runtimes |
| Julia/SciML evaluation backend | Competing prototype if justified | Move the numerical loop/RHS together; retain config, persistence and reporting contracts |
| Process-level candidate parallelism | Later throughput experiment | Does not accelerate a single periodic solve; define deterministic scheduling |
| GPU or whole-project rewrite | Defer | No evidence from this small-state workload justifies the added scope |

Numba's documented no-Python compilation is one prototype option, but object
wrappers need a suitable numerical boundary. Start without `fastmath`; it can
relax floating-point semantics. Do not infer speed from unrelated benchmark
examples. See [Numba performance tips](https://numba.readthedocs.io/en/stable/user/performance-tips.html).

A Julia prototype should measure complete evaluations, including compilation,
first-call latency and steady repeated throughput separately. In-place RHS and
allocation discipline are relevant there too; a literal port of duplicated work
is not a sound benchmark. See [SciML ODE optimization](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/faster_ode_example/).
There is no measured Julia speedup for this repository and no migration decision.

SciPy already provides compiled LSODA integration; Python remains in the RHS
and orchestration. Its documentation exposes solver counters and warns that
vectorization is not universally faster. Changing `vectorized=True` is not a
batch-candidate solver. See [solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html).

For parallel campaigns, isolate LSODA in worker processes rather than assume
thread reentrancy. Incumbent-centered proposals depend on earlier results:
parallel completion order changes centers, warm-start choices and possibly the
search trajectory. Use explicit reproducible batches/ordered assimilation if
parallelism is approved, with one persistence owner and bounded worker budgets.
Do not claim the same search merely because the Sobol seed is unchanged.

JSON/checkpoint I/O was not established as a bottleneck here: the replay measures
physical evaluation, not whole-campaign disk writes. Measure it separately before
weakening immediate history persistence or crash recovery.

## Acceptance contract for the later implementation

Before modifying code, select a small frozen benchmark set: legacy 37.6261105666 W
regression, current production four-stage candidate, nearby candidate from a
compatible warm start, cold start, one invalid-domain/interrupted case and an
existing smooth/six-bar case. The production input identity in this note is
reproducible even if the campaign later advances.

For algebra-preserving changes, compare RHS vectors over saved real trajectories,
including valve closure/pressure equality, zero flow, reversed link flow and
stage boundaries. Preserve domain classification and exception/status semantics.
Then compare complete-cycle and periodic results at **unchanged tolerances**:
external heat, work, efficiency, mass/energy balance, pressure/temperature/flow
extrema and domain fractions. Derive acceptance tolerances from existing tests
and measured baseline repeatability, not from a convenient percent error.
Avoid changing operation order unnecessarily; investigate changed feasibility
or candidate rankings rather than writing a looser assertion.

Numerical algorithm changes require a separate accuracy study and review.
Smoothing the four-stage law changes the candidate motion and cannot count as
an unchanged-problem solver speedup. Never remove pause heat, variable properties,
wall storage, model-domain checks, motor reversal or the external-heat efficiency
boundary. Preserve the ten physical states and five diagnostic quadratures.

Report median and tail latency, cycles, RHS/Jacobian counts, memory and complete
evaluation throughput across repeated runs. Include slow/invalid cases and
compilation cost. Use Amdahl's law to set expectations from **non-overlapping**
measured costs; do not multiply optimistic speedups for overlapping operations.
Run the full regression suite after implementation and retain the reference path
until equivalence is demonstrated.

## Questions and deliverable for Chat's second review

Please review this note and the measurement artifact before proposing code:

1. Is sharing one diagnostic pass and one instantaneous flow evaluation the
   correct first boundary? Sketch the smallest API that preserves generic
   kinematics/exchanger composition and avoids duplicating conservative balances.
2. Which work belongs in immutable prepared candidate data, evaluation-local
   data, and reporting? Identify risks around signed flows, blocked valves,
   validity and static versus dynamic thermal conductance.
3. Which additional measurements distinguish RHS overhead from stiffness and
   periodic convergence? Keep the benchmark set bounded and representative.
4. Should a compiled Python kernel or Julia backend be prototyped after the
   first refactor, and what measured criterion would justify it? Specify fair
   warm/cold/compilation comparisons, not an assumed speedup.
5. Is earlier bounded wall acceleration justified, or should it wait? Define
   safeguarding and proof of ordinary-cycle convergence.
6. Produce a prioritized implementation plan with small reviewable changes,
   regression requirements, go/no-go gates and explicit deferred work.

Place that independent response in `docs/SOLVER_ACCELERATION_STAGE2.md` when it
exists. The later implementation should read both reviews, reconcile disagreements
explicitly, and update decisions based on measurements. No optimization algorithm,
physical approximation or language migration is approved by this stage-1 note.

## Scope and verification of this review

Only this note, discovery links and the profiling evidence were added. The three
replays agree exactly on their reported efficiency, power and periodic error.
The probe source is embedded in the JSON for reproduction at the recorded code
commit; use a separate output filename/directory for later measurements and
verify input hashes before comparing them. The coarse-timer probe temporarily
wraps functions in its own process; it does not edit repository code.

Local Markdown links, evidence arithmetic and `git diff --check` were checked.
The full test suite was not rerun for this documentation-only reflection; it
must be run after the implementation phase. No compiled backend, Jacobian,
new warm-start algorithm or numerical setting was installed or tested here.

## Literature-based Julia assessment while temperature campaigns run

Follow-up, 2026-09-17: the user defers all new benchmarking to avoid competing
with the active temperature campaign. This addition uses published evidence and
the measurements already recorded above. No new numerical evaluation is run.

[The SciML cross-language work-precision comparison](https://docs.sciml.ai/SciMLBenchmarksOutput/stable/MultiLanguage/ode_wrapper_packages/)
compares stiff and non-stiff examples at measured solution accuracy. Its stated
versions are Julia 1.7, SciPy 1.6.1, MATLAB 2019B and deSolve 1.3.0: this is
historical evidence, not a benchmark of our current installation. The SciPy
wrapper uses a Julia-JIT ODE function; the authors explicitly distinguish that
from a normal SciPy/Python setup. Consequently neither plotted speed ratios
nor their wrapper-specific acceleration claim can be transferred to DADA.
This is a benchmark maintained by the SciML authors, not an independent
measurement of our solver.

[Rackauckas and Nie, *Confederated Modular Differential Equation APIs for
Accelerated Algorithm Development and Benchmarking*](https://arxiv.org/abs/1807.06430)
explains the modular solver/comparison approach. Together with the
[SciML optimization tutorial](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/faster_ode_example/),
it supports comparing complete specialized ODE implementations rather than
assuming that the programming language alone determines runtime. The tutorial
emphasizes allocation discipline and in-place numerical work. It does not
supply a DADA-specific speedup.

**Inference for this repository:** a compiled, allocation-conscious kernel is a
credible opportunity because the physical system is small (10 states plus 5
quadratures), yet requires thousands of scalar RHS calls per cycle. Much of
our measured workload is Python reconstruction, dispatch and diagnostic-object
creation, not a large matrix factorization. Moving the whole numerical loop
and RHS together could remove these overheads. A several-fold gain against the
current object-heavy implementation is plausible, but remains a hypothesis,
not a measured result or a factor inferred from another paper. The incremental
gain over an equally cleaned, compiled Python-accessible kernel is unknown.

Distinguish cleanup gains, compilation/runtime gains, integrator/step-selection
gains and periodic-iteration gains. Julia alone does not guarantee fewer cycles
to periodic convergence. Identical nominal tolerances across different solvers
also do not guarantee identical actual error: compare work versus measured
accuracy and preserve the physical validity decisions.

A useful bound follows from our existing coarse timing, without new computation:
if only the 6.618 s integration became arbitrarily fast, the unchanged 13.688 s
remainder would limit the 20.306 s warm evaluation to about 1.48x. An illustrative
5x integration acceleration alone would give approximately 15.012 s total,
only 1.35x. These are Amdahl scenarios, **not Julia predictions**. They explain
why the user's planned cleanup of diagnostics and redundant calculations is
essential. For a candidate requiring many cycles, the balance will differ.

The later Julia/Python comparison must cover both smooth natural motion and
four-stage motion. Four-stage kinematics temporarily simplifies features found
in freer-motion searches; it is not the permanent mechanical or numerical
architecture. Smooth motion may alter step selection, but that is a different
motion problem and must not be counted as a language-only acceleration.

Practical conclusion: literature and the measured call structure justify giving
Julia a serious, bounded prototype comparison later. They do not justify either
an immediate whole-project port or a promised 10x/100x speedup. Keep configuration,
history and human reporting in their existing layer if only the computational
backend needs replacement. No new language decision is made here.
