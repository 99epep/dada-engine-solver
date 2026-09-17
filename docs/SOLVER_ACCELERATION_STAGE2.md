# Solver acceleration — stage 2 independent review and implementation gates

Date: 2026-09-18. Status: **independent architectural review; no production code change approved by this document alone**.

This review follows [SOLVER_ACCELERATION_STAGE1.md](SOLVER_ACCELERATION_STAGE1.md) and checks its main conclusions against the current solver structure, especially:

- [dynamics.py](../src/dada_solver/dynamics.py);
- [air_wall.py](../src/dada_solver/exchangers/air_wall.py);
- [wall_cycle.py](../src/dada_solver/exchangers/wall_cycle.py);
- [wall_iteration.py](../src/dada_solver/exchangers/wall_iteration.py);
- [gas_film.py](../src/dada_solver/exchangers/gas_film.py);
- [gas_diagnostics.py](../src/dada_solver/exchangers/gas_diagnostics.py);
- [hardware.py](../src/dada_solver/exchangers/hardware.py);
- [campaign/evaluator.py](../src/dada_solver/campaign/evaluator.py);
- the Stage-1 measurement artifact
  [solver_performance_stage1_20260917.json](../outputs/solver_performance_stage1_20260917.json).

The purpose of this document is to define a small, reviewable acceleration path before implementation. It deliberately separates:

1. removal of duplicated algebra and reporting work;
2. measurement of the resulting Python implementation;
3. numerical-method experiments;
4. optional compiled-kernel experiments;
5. optional Julia/SciML experiments;
6. periodic fixed-point acceleration;
7. campaign-level parallelism.

The first implementation pass should change **none** of the physical equations, validity domains, convergence tolerances, search logic, or accepted/rejected candidate semantics.

---

## 1. Executive decision

Stage 1 identified the correct first direction. The first implementation should be a **structural cleanup inside the existing Python solver**, not a language migration and not a new ODE algorithm.

The two first targets are:

1. **compute final microtube diagnostics once per completed trajectory and reuse the result;**
2. **compute the instantaneous thermodynamic/hydraulic point once per RHS call, then assemble the conservative derivative using the actual wall heat rates without recomputing the flows.**

The second item is the more important architectural change. The current wall path effectively performs:

```text
state
  -> model.evaluate(...) to obtain flows for gas-film context
  -> gas-film / wall heat calculation
  -> dataclasses.replace(model, prescribed heat transfers)
  -> model.evaluate(...) again to obtain conservative rates
```

The replacement should instead be:

```text
state
  -> one instantaneous point:
       geometry + dV/dt
       temperatures + pressures
       effective valve topology
       signed hydraulic flows
  -> gas-film / wall heat calculation using that point
  -> one conservative balance assembly using:
       the same point
       actual cold/hot gas heat rates
```

There must remain **one implementation of the conservative balances**. This is not permission to introduce a separate "fast physics" path.

The existing public `ThermodynamicModel.evaluate()` should remain available and should be reimplemented as a compatibility façade over the same split primitives.

---

## 2. What the Stage-1 measurements establish

The Stage-1 warm replay is not a full campaign benchmark, but it proves enough to prioritize the first changes.

For the recorded production candidate:

- one unprofiled warm evaluation took about **20.96 s**;
- the coarse-timed repeat took about **20.31 s**;
- integration took about **6.62 s** in the coarse-timed run;
- generic diagnostics took about **2.59 s**;
- microtube diagnostics took about **5.67 s**, then the same complete diagnostic routine was run again for about **5.31 s**;
- the profiled run recorded **6080 integrator RHS calls**;
- `ThermodynamicModel.evaluate()` was called **24942 times**;
- `cycle_microtube_diagnostics()` was called twice;
- `MicrotubeGasModel.diagnose()` was called **159894 times**;
- recursive `dataclasses.asdict()` activity was extremely large and expensive.

The exact profile percentages are distorted by cProfile overhead and must not be interpreted as unprofiled speedup factors. The call counts, duplicated paths and coarse boundary timings are nevertheless actionable.

A particularly important conclusion is that there are **two different kinds of duplication**:

### 2.1 Duplication during integration

`AirWallMotor.thermal_rates()` requests flow contexts. For the variable gas film those contexts call `ThermodynamicModel.evaluate()` to obtain flows. `AirWallMotor.derivative()` then constructs a replacement thermodynamic model with prescribed heat rates and calls `evaluate()` again.

This repeats geometry/state reconstruction and, most importantly, the hydraulic network solution for the same `(theta, state)`.

### 2.2 Duplication after integration

The completed trajectory is replayed through multiple diagnostic consumers. In the campaign evaluator, `_tube_validity()` invokes the complete variable-film microtube diagnostic report, then `_wall()` invokes it again to retain domain information. Generic cycle diagnostics also reconstruct thermal states and heat rates separately.

The first implementation should address these categories independently so that their gains remain measurable.

---

## 3. Recommended numerical boundary

### 3.1 Smallest useful API split

The current `ThermodynamicModel.evaluate(theta, state, topology)` combines:

- kinematics;
- primitive thermodynamic reconstruction;
- valve state selection;
- four hydraulic flows;
- mass and enthalpy transport;
- heat transfer;
- boundary work;
- conservative derivative assembly;
- residual reporting.

That combination is convenient for reservoir/static-heat models but forces the air-wall model to calculate the hydraulic state before the wall heat rates are known, then calculate it again.

The smallest useful refactor is therefore a two-stage model API.

Conceptually:

```python
point = model.instantaneous_point(theta, state, topology)
rates = model.assemble_rates(
    point,
    cold_heat_rate=...,
    hot_heat_rate=...,
)
```

and the current API becomes:

```python
def evaluate(theta, state, topology):
    point = instantaneous_point(theta, state, topology)
    q_cold = cold_heat_transfer.heat_rate(point.temperatures[COLD])
    q_hot = hot_heat_transfer.heat_rate(point.temperatures[HOT])
    return assemble_rates(point, q_cold, q_hot)
```

Names are not mandated by this review. The separation is.

### 3.2 Proposed instantaneous result

A small frozen typed object is preferable during the first Python refactor because it makes the semantic boundary explicit and testable.

It should contain only values already calculated at that instant, for example:

```text
theta
volumes: S, L, C, H
volume rates: S, L
temperatures: S, L, C, H
pressures: S, L, C, H
effective valve topology
four signed flows
four choking flags
```

It may also contain small directly reusable primitive quantities if measurement later shows that they avoid material duplicate work.

It should **not** contain:

- reporting dictionaries;
- serialized diagnostic structures;
- campaign metadata;
- candidate IDs;
- accumulated extrema;
- approximate caches keyed by angle;
- wall heat rates that have not yet been calculated.

The signed flow convention must remain exactly the convention of `NetworkFlows`.

### 3.3 Conservative assembly

The existing balance logic in `ThermodynamicModel.evaluate()` should be factored, not copied.

The balance assembler should consume:

- the instantaneous point;
- the two actual gas heat rates.

It should then produce the same `ModelRates` semantics as today:

- eight gas-state derivatives;
- signed flows;
- cold/hot heat rates;
- gas work;
- mass residual;
- energy residual.

The functions `_apply_bidirectional_link()` and `_apply_directed_link()` can remain the authoritative transport logic initially. If later compilation requires flatter scalar algebra, that should be a later refactor with equivalence tests.

### 3.4 Wall path after the split

`AirWallMotor.derivative()` can then become conceptually:

```text
gas state from first 8 state components
        |
        v
model.instantaneous_point(...)
        |
        +----> build Hi/Ho flow contexts from the same point
        |            |
        |            v
        |       dynamic gas films
        |            |
        |            v
        |       two gas heat rates
        |
        v
model.assemble_rates(point, q_i, q_o)
        |
        v
gas derivatives + wall derivatives + quadratures
```

No `dataclasses.replace()` of the model is required per RHS call.

### 3.5 Flow-context construction

`AirWallMotor.flow_contexts()` currently reconstructs pressures and calls `model.evaluate()`. Add an internal path that builds the same contexts **from the instantaneous point**.

For compatibility, a public or legacy helper may still accept `(angle, values)` and internally obtain the point once. The RHS must use the already-created point.

The existing important closed-valve behavior must remain unchanged:

> a closed valve can support a pressure jump, but zero flow must not be interpreted as an axial tube pressure gradient.

Therefore the current distinction between stored port pressure ratio and the zero-flow passage pressures must survive the refactor exactly.

### 3.6 Effective topology

`ThermodynamicModel.effective_topology()` currently reconstructs pressures by calling `state.pressures(..., self.volumes(theta))`. If `instantaneous_point()` has already calculated volumes and pressures, it must not call this path and repeat them.

The clean solution is either:

- a private topology helper operating on already known pressures; or
- topology selection directly inside instantaneous-point construction.

The existing public helper can remain for callers that do not have a prepared point.

---

## 4. Prepared candidate data, evaluation-local data and reporting

The solver should distinguish three lifetimes.

## 4.1 Candidate-invariant prepared data

Anything that is mathematically invariant during one candidate evaluation can be prepared once.

Good candidates include:

### Machine and gas constants

- angular speed and frequency;
- constant control-volume sizes;
- gas constant and calorically perfect `cp`, `cv`, `gamma` where applicable;
- any immutable kinematic constants already owned by the selected motion model.

### Microtube geometry

Values currently repeatedly obtained from `MicrotubeBank.dimensions()` or reconstructed algebraically:

- inner and outer diameter;
- tube length and half length;
- tube count;
- total tube flow area;
- internal heat-transfer area;
- relevant header/core geometric constants.

Do not create a second formula for these values. Preparation should obtain them from one authoritative geometry calculation.

### Wall and external-air constants

For one candidate with constant external air inputs:

- wall capacities;
- half-wall resistances;
- air capacity rates `m_dot * cp`;
- the external-air effective conductance

```text
C_dot * (1 - exp(-UA / C_dot))
```

when `C_dot > 0`.

This quantity is currently recomputed in `AirWallExchanger.rates()` although the air-side conductance, air mass flow and air `cp` are candidate constants.

### Hydraulic fixed factors

Where exact equivalence permits:

- geometric areas;
- half lengths;
- fixed loss multipliers;
- header-loss constants;
- valve CdA;
- pressure-independent prefactors.

Do **not** precompute values that depend on instantaneous temperature, density, pressure ratio, flow direction, Reynolds regime, slip regime or gas properties.

### Preparation invalidation rule

Prepared data belongs to one fully built candidate. Do not introduce a global cache indexed by a hand-written subset of geometry fields. A missing key would silently reuse the wrong physics.

The safest first implementation is preparation during candidate/model construction with immutable ownership.

---

## 4.2 Evaluation-local values

These change with `(theta, state)` and must be calculated exactly from the current adaptive integrator trial state:

- cylinder volumes and derivatives;
- masses and internal energies;
- temperatures;
- pressures;
- effective valve states;
- four signed flows;
- upstream side selected by flow sign;
- density;
- variable viscosity/conductivity/`cp` where applicable;
- Reynolds, Mach, Knudsen and other instantaneous domain quantities when required;
- instantaneous internal gas-film conductance;
- gas-wall heat rates;
- conservative state derivatives.

Reuse is encouraged **within the same exact instantaneous evaluation**.

Do not memoize by `theta` alone. Adaptive ODE solvers can evaluate the same or nearly the same angle at different trial states, including rejected steps.

Do not round state or temperature to increase cache hits in the first acceleration pass.

---

## 4.3 Reporting-only data

Rich diagnostic objects and serialization should be kept out of the integration hot path unless a value is physically required for validity enforcement at that instant.

Reporting-only work includes:

- converting diagnostic dataclasses recursively to dictionaries;
- constructing large nested domain summaries;
- human-readable issue collections;
- time/heat weighted fractions over a completed cycle;
- extrema aggregation;
- JSON-ready structures.

A useful rule is:

> **hot-path physics returns typed/scalar numerical facts; reporting materializes dictionaries at the boundary.**

This is particularly relevant because Stage 1 recorded very large `asdict()` and `deepcopy()` costs.

This must not weaken instantaneous model-domain rejection. If an instantaneous physical closure currently calls `require()` or raises `MicrotubeDomainError`, the optimized path must retain the same check at the same physical state.

---

## 5. Gas-film and hydraulic property reuse

The variable-property microtube path deserves special care because an apparently harmless cache can alter the model.

### 5.1 Safe local reuse

Within one exact instantaneous evaluation, it is valid to reuse a gas property evaluated at exactly the same temperature and for the same role.

Examples:

- conductivity used twice for the two half-port heat-transfer estimates when the bulk temperature is the same by model definition;
- viscosity used repeatedly in diagnostics for an identical upstream temperature;
- fixed gas constants.

### 5.2 Unsafe conflations

Do not merge quantities merely because they are numerically close.

In particular:

- the exchanger lumped bulk gas temperature and the hydraulic upstream temperature are not generally interchangeable;
- reverse flow changes the upstream node;
- zero flow has special pressure-gradient treatment;
- dynamic internal conductance must not be replaced by the legacy static gas conductance;
- the external-air effective conductance may be static while the gas-side film remains dynamic.

The current documentation correctly distinguishes the legacy `overall_static_conductance_w_k` from the dynamic internal-film model. Preserve this distinction.

---

## 6. Final diagnostics: one trajectory replay, not two

The first low-risk patch should remove the duplicated complete call to `cycle_microtube_diagnostics()`.

Current campaign flow effectively does:

```text
_tube_validity()
    -> cycle_microtube_diagnostics()

then

_wall()
    -> cycle_microtube_diagnostics()
```

For the variable-film path, calculate the report once and derive from it:

- maximum hydraulic-upstream Reynolds;
- maximum hydraulic-upstream Mach;
- model-domain failure criteria;
- physical-time-weighted domain fractions;
- absolute-heat-weighted domain fractions;
- retained detailed passage report.

The generic/static fallback in `_tube_validity()` may remain separate if it is needed for models without variable flow-context diagnostics.

### 6.1 Next diagnostic refactor, only after the first patch

After measuring the single-pass change, consider a shared completed-trajectory replay object that reconstructs for each accepted sample:

- thermodynamic state;
- instantaneous point;
- wall heat rates;
- gas-film diagnostics.

Multiple aggregators could consume that shared sample without recomputing it.

Do not combine this larger diagnostic redesign with the first patch unless necessary. The purpose is to retain reviewability and make the measured gain attributable.

---

## 7. Instrumentation required before numerical-method experiments

Stage 1 correctly notes that the current wall wrapper discards information that SciPy already exposes.

Retain per segment and per cycle:

- solver method;
- `nfev`;
- `njev`;
- `nlu`;
- number of returned/accepted samples;
- angular interval;
- integration elapsed time;
- RHS elapsed time if a low-overhead aggregate timer is retained;
- minimum/maximum/representative accepted angular step statistics;
- whether the segment ends at a kinematic breakpoint;
- completion/interruption status.

Retain per candidate:

- number of physical cycles completed;
- periodic normalized error after each cycle;
- whether an extrapolated wall initial guess was used;
- integration time;
- final generic-diagnostic time;
- final microtube-diagnostic time;
- reporting/serialization time;
- total evaluation time;
- warm-start source and normalized design distance;
- build/preflight time;
- status, including failed/invalid/interrupted evaluations.

Counters should be available without enabling cProfile.

### 7.1 What distinguishes Python/RHS overhead from stiffness

No single metric proves stiffness. Use the following evidence together:

- `nfev` per accepted step;
- `njev` and `nlu` for methods that expose them;
- clustering of very small accepted steps away from mandatory stage boundaries;
- comparison of LSODA against at least one explicit and one implicit SciPy method on frozen cases, **only after** algebra-preserving cleanup;
- measured work versus final physical accuracy, not nominal tolerance alone.

A high RHS time with modest solver bookkeeping suggests a compiled scalar kernel may help.

A large implicit-factorization count or severe step collapse suggests a numerical-method/stiffness problem rather than only Python dispatch.

The four-stage breakpoints are real non-smooth locations and should remain exact segment boundaries. Small steps near those boundaries are not by themselves evidence of physical stiffness.

---

## 8. Frozen benchmark set

Before the first production modification, create a small benchmark manifest with hashes/identities rather than relying on "current best" names that will move as campaigns continue.

The set should contain at least:

1. **legacy regression** producing the established approximately
   `37.6261105666 W` result;
2. **Stage-1 production candidate** using the recorded saved ten-state state and settings;
3. **nearby compatible warm-start candidate** requiring several cycles;
4. **cold-start version** of a representative production candidate;
5. **slow or invalid-domain case**;
6. **interrupted/deadline case**;
7. **smooth/natural or six-bar motion case**, so acceleration is not tuned only to the current four-stage law.

The benchmark should record:

- source commit;
- input/report hashes;
- Python, NumPy and SciPy versions;
- OS/platform;
- candidate parameters;
- initial state;
- numerical settings.

### 8.1 Debian migration

The laptop is about to move from Debian 12 to Debian 13. Therefore:

- finish this architectural review first;
- install Debian 13;
- run the full test suite;
- rerun the frozen baseline with **unchanged code** on Debian 13;
- treat those measurements as the performance baseline for implementation.

Do not compare a post-refactor Debian-13 timing directly with the old Debian-12 Stage-1 timing and call the difference an optimization gain.

The Stage-1 artifact remains valuable for call structure and historical comparison.

---

## 9. Acceptance contract for algebra-preserving refactors

The first patches are intended to be algebra-preserving. Their acceptance should be stronger than "similar efficiency".

### 9.1 Instantaneous equivalence tests

Build test points from saved real trajectories and explicit edge cases.

For each point compare old/reference and refactored paths for:

- volumes;
- volume rates;
- temperatures;
- pressures;
- effective topology;
- all four signed flows;
- choking flags;
- cold/hot gas heat rates;
- eight conservative derivatives;
- gas work rate;
- mass residual;
- energy residual.

Include explicitly:

- positive and reverse bidirectional flow;
- zero/near-zero flow;
- open and blocked passive valves;
- pressure equality;
- valve reversal threshold neighborhood;
- four-stage breakpoints on both sides;
- low/high temperature states from real trajectories;
- valid laminar microtube state;
- supported turbulent state if present in current production model;
- a state that must raise the current domain error.

The test should use tight tolerances derived from operation-order changes. If the new structure can preserve operation order, prefer exact or near-machine equivalence rather than writing a broad percentage allowance.

### 9.2 Completed-cycle equivalence

At unchanged ODE and periodic tolerances compare:

- final ten physical states;
- five quadratures;
- external heat in/out;
- gas work;
- indicated power;
- efficiency;
- conservation residuals;
- pressure, temperature and flow extrema;
- microtube domain classification;
- domain fractions;
- maximum Reynolds and Mach;
- status and exception semantics.

### 9.3 Campaign semantics

Verify that the refactor does not change:

- candidate feasibility;
- constraint availability;
- rejection reason classes;
- warm-start compatibility;
- retained last complete physical state after interruption;
- checkpoint/history schema unless intentionally versioned.

A changed candidate ranking or feasibility boundary must be investigated, not accepted by simply loosening tests.

---

## 10. Prioritized implementation plan

Each step should be committed and benchmarked separately where practical.

### P0 — establish the new-machine baseline

After Debian 13 installation:

1. install the project dependencies;
2. run the full regression suite;
3. freeze the benchmark manifest;
4. run repeated baseline measurements;
5. record median and tail latency and numerical results.

**Go gate:** tests pass and frozen cases reproduce expected physics.

No solver change before this gate.

---

### P1 — eliminate duplicate final microtube report

Change `_wall()` / `_tube_validity()` so the variable-film diagnostic report is produced once and reused.

Do not redesign the diagnostic object yet.

Measure:

- completed diagnostic time;
- end-to-end warm replay;
- representative many-cycle evaluation.

**Go gate:** numerical/report equivalence and a measurable reduction in final diagnostic cost.

---

### P2 — introduce the instantaneous-point / conservative-assembly split

Refactor `ThermodynamicModel.evaluate()` internally while keeping its current external semantics.

Add:

- one authoritative instantaneous thermodynamic/hydraulic calculation;
- one authoritative conservative assembler.

Make `evaluate()` call both.

Initially switch no wall code to the new split until direct equivalence tests pass.

**Go gate:** instantaneous and reservoir-path equivalence tests pass.

---

### P3 — switch `AirWallMotor.derivative()` to one hydraulic evaluation

Use the instantaneous point to:

- create both exchanger flow contexts;
- calculate both dynamic gas films;
- calculate wall/air rates;
- assemble the conservative derivative.

Remove per-RHS `dataclasses.replace()` and the second hydraulic evaluation.

Measure RHS count and cost separately from accepted step count.

**Go gate:** identical physical/domain semantics and complete-cycle equivalence; clear RHS cost reduction.

---

### P4 — prepare immutable candidate constants

Move clearly invariant calculations out of hot paths, starting with the highest measured repeat counts.

Examples:

- tube geometry/dimensions;
- areas and half lengths;
- wall resistance constants;
- external-air effective conductance;
- fixed hydraulic prefactors.

Keep variable transport variable.

**Go gate:** per-RHS improvement without changed physical results.

---

### P5 — reporting allocation cleanup

Remove recursive `asdict()` from per-sample aggregation.

Aggregate typed/scalar values and serialize once at the output boundary.

If useful, introduce a single shared completed-trajectory replay to serve generic and microtube aggregators, but only after the simpler duplication fixes have been measured.

**Go gate:** diagnostic outputs equivalent, including weighting and failure strings/semantics where part of public output.

---

### P6 — retain solver counters and run the post-refactor benchmark

At this point establish where time remains.

Produce for every frozen case:

```text
build
periodic integration
  cycles
  accepted steps
  nfev / njev / nlu
  RHS time
final diagnostics
serialization
total
```

Only now decide whether the next bottleneck is:

- RHS scalar algebra;
- solver algorithm/step selection;
- periodic fixed-point convergence;
- final diagnostics;
- campaign orchestration.

This is the principal **go/no-go gate for compiled backends**.

---

## 11. Compiled Python kernel decision

A compiled Python-accessible kernel is worth prototyping only if, after P1–P6, a coherent numerical kernel still represents a large part of complete-evaluation cost.

A practical prototype gate is:

- numerical integration/RHS work is still roughly **half or more** of representative multi-cycle end-to-end time; and
- Amdahl analysis shows that accelerating this part could plausibly improve complete candidate latency by at least about **1.5x**.

The exact threshold is an engineering rule, not a physical claim. Its purpose is to avoid optimizing a minority cost.

### 11.1 Numba prototype shape

If Numba is selected, compile a **coherent scalar kernel**, not many tiny Python↔JIT crossings.

Prefer:

- no-Python mode;
- explicit numerical arrays/scalars;
- prepared immutable parameter arrays/structures;
- no reporting dictionaries;
- no `fastmath` initially;
- no approximate memoization.

Numba documentation continues to distinguish nopython execution and optional `fastmath`; relaxing floating-point semantics should be a separate experiment:
https://numba.readthedocs.io/en/latest/user/performance-tips.html

### 11.2 Retention criterion

A compiled Python kernel should not become production complexity for a marginal gain.

Suggested retention criterion:

- substantial end-to-end improvement on the frozen benchmark, not only a microbenchmark;
- no important tail-latency regression;
- unchanged physics/status semantics;
- maintainable reference Python path retained until confidence is high.

A median complete-evaluation speedup around **1.5x or better** would justify serious integration work. Below that, maintenance cost should be weighed carefully.

This number is intentionally a project-management gate, not a promised achievable speedup.

---

## 12. Julia/SciML decision

Julia remains a credible **competing prototype**, not a predetermined rewrite.

The current SciML interfaces explicitly support in-place ODE functions and specialization choices. Their documentation notes that in-place forms can avoid allocations and that full specialization trades compilation time for runtime performance, particularly relevant to long optimization loops:

- https://docs.sciml.ai/SciMLBase/stable/interfaces/Differential_Equation_Problem_Types/
- https://docs.sciml.ai/DiffEqDocs/stable/api/ordinarydiffeq/api/common_interface/
- https://docs.sciml.ai/DiffEqDocs/stable/tutorials/faster_ode_example/

The older cross-language SciML benchmark remains useful as motivation for measurement but not as a DADA speedup estimate:
https://docs.sciml.ai/SciMLBenchmarksOutput/stable/MultiLanguage/ode_wrapper_packages/

### 12.1 What to port for a fair prototype

Do not port the complete project.

A meaningful Julia prototype should include together:

- prepared candidate numerical parameters;
- the complete instantaneous RHS physics required by the selected frozen case;
- the ODE integration loop across exact phase breakpoints;
- enough periodic-cycle iteration to measure realistic repeated throughput.

Keep initially in Python:

- campaign configuration;
- candidate generation;
- checkpoint/history format;
- human reporting;
- high-level optimization orchestration.

A Julia function called thousands of times from Python for individual scalar operations would be an unfair architecture.

### 12.2 Measurements

Report separately:

1. environment/process startup;
2. first compilation/first solve latency;
3. second identical solve;
4. repeated same-structure candidates;
5. complete periodic evaluation;
6. cold and compatible warm starts.

Compare against the **cleaned Python reference**, not the current duplicated implementation.

Compare physics at measured accuracy, not only equal nominal tolerances.

### 12.3 Julia production gate

Because a second language and FFI/backend boundary impose more maintenance than an internal Python optimization, Julia should clear a stronger bar.

A reasonable production-integration gate is approximately:

- **2x or better complete-evaluation improvement** over cleaned Python on representative multi-cycle work, or another clearly material advantage such as substantially better robust numerical behavior;
- equivalent physical outputs/domain decisions;
- acceptable first-call/compilation behavior for the intended campaign workflow;
- a clean backend API that does not duplicate reporting/orchestration logic.

If Julia is only slightly faster than a cleaned/compiled Python kernel, prefer the simpler architecture.

Again, this is a project decision threshold, not a prediction that Julia will or will not reach it.

---

## 13. ODE algorithm experiments

Do not mix an ODE-method change with P1–P5.

After instrumentation exists, compare methods on frozen cases if the counters justify it.

Possible questions:

- Is LSODA doing expensive switching/factorization work?
- Does an explicit method take excessive steps?
- Does Radau/BDF materially reduce evaluations for the same measured physical accuracy?
- Is the fixed `pi/360` maximum step actually binding often?
- Are small steps concentrated at real stage discontinuities or throughout smooth portions?

Preserve exact stage breakpoints.

Do not count a changed motion law, smoothed discontinuity or looser tolerance as solver acceleration.

A solver-method change is accepted only after an accuracy/work study.

---

## 14. Periodic convergence acceleration

### 14.1 Current wall Aitken behavior

The current `extrapolate_wall_energies()` is already conservative in several useful ways:

- wall energies only;
- componentwise Aitken form;
- only when consecutive steps have the same sign and are shrinking;
- correction bounded to ±30 K equivalent per wall;
- positive wall energy required;
- applied only between complete cycles;
- the extrapolated state itself is never accepted as convergence;
- a following ordinary physical cycle must satisfy the original periodic criterion.

These safeguards should remain.

### 14.2 Recommendation

**Do not make it earlier or more aggressive in the first acceleration implementation.**

Stage 1 shows a candidate whose error decreases slowly, then whose cycle-10 wall extrapolation causes a large transient normalized error before convergence on the following cycle. That proves the mechanism can be useful but does not establish the optimal trigger.

First measure after the RHS cleanup:

- wall-only endpoint residuals;
- gas-state endpoint residuals;
- error contraction ratios per cycle;
- effect of warm-start distance;
- cycles and wall time saved by the existing cycle-10 trigger.

### 14.3 Later bounded experiment

If periodic cycles remain the dominant cost, an earlier/adaptive experiment is justified.

Minimum safeguards:

- operate only on initial guesses between complete cycles;
- retain the last physically completed state separately;
- preserve total gas inventory;
- preserve positivity;
- wall temperatures remain inside a physically defensible bound around observed completed-cycle values;
- reject extrapolation when residual direction oscillates or contraction evidence is poor;
- maximum correction remains bounded;
- if the accelerated next cycle worsens markedly, fall back to the last ordinary completed state;
- convergence is certified only by an unmodified complete physical cycle with the original criterion.

Do not lower `maximum_cycles` to turn slow convergence into apparent success.

Direct shooting/Newton remains deferred until there is evidence that fixed-point iteration is the limiting factor after cheaper measures.

---

## 15. Campaign parallelism

Parallel candidate evaluation can improve throughput but is independent of single-candidate latency.

Do not introduce it in the first acceleration patches.

The current incumbent-centered/warm-start search has stateful dependencies: completion order can alter the incumbent, warm-start source and subsequent trajectory.

If later approved, use explicit deterministic batches or ordered assimilation, one persistence owner and independent worker budgets.

Do not claim an identical optimization search merely because the Sobol seed is unchanged.

LSODA/thread safety assumptions should not be relied upon without verification; process-level isolation is the conservative initial experiment.

---

## 16. Explicitly deferred work

The following are **not** part of the first implementation:

- whole-project Julia rewrite;
- GPU work;
- campaign multiprocessing;
- direct shooting/Newton periodic solve;
- Jacobian derivation;
- sparse Jacobian engineering;
- tolerance relaxation;
- maximum-step relaxation;
- trajectory decimation that affects diagnostics;
- approximate property lookup tables;
- angle-only memoization;
- motion-law smoothing;
- validity-domain simplification;
- removal of instantaneous domain checks;
- altered optimizer/search logic;
- changed physical models.

Any of these may become sensible later, but none is justified as a prerequisite for removing known duplicate work.

---

## 17. Expected sequence for the next implementation pass

The later coding pass should read Stage 1 and Stage 2 together and follow this sequence unless new measurements contradict it:

```text
Debian 13 migration
        |
        v
frozen baseline + full tests
        |
        v
P1 one final microtube diagnostic pass
        |
        v
measure
        |
        v
P2 split instantaneous point / conservative assembly
        |
        v
direct equivalence tests
        |
        v
P3 one hydraulics evaluation per wall RHS
        |
        v
measure
        |
        v
P4 prepare invariant candidate constants
        |
        v
measure
        |
        v
P5 reporting/allocation cleanup
        |
        v
P6 complete counters + benchmark
        |
        +----------------------+----------------------+
        |                      |                      |
        v                      v                      v
RHS still dominant?      cycles dominant?      diagnostics dominant?
        |                      |                      |
        v                      v                      v
Numba/Julia gate      bounded periodic        further replay/
and ODE study          acceleration            reporting cleanup
```

This ordering is intentionally conservative. It obtains low-risk gains first and creates a much cleaner numerical boundary for any later compiled backend.

---

## 18. Reconciliation with Stage 1

### Agreed

This review agrees with Stage 1 that:

- duplicated final microtube diagnostics should be removed first;
- the same instantaneous state/flows should not be solved twice per RHS;
- immutable candidate constants should be prepared;
- property reuse should be local and exact;
- solver counters must be retained before numerical-method changes;
- warm starts and wall acceleration must remain safeguarded;
- Julia deserves a bounded prototype only after cleanup;
- parallel throughput and single-candidate latency are separate problems;
- GPU/whole-project rewrite is not justified;
- unchanged-model refactors require strong regression checks.

### Refinements added by Stage 2

Stage 2 makes the implementation boundary more explicit:

1. split `ThermodynamicModel.evaluate()` into one instantaneous thermodynamic/hydraulic point and one conservative assembler;
2. retain `evaluate()` as the compatibility façade, so there is only one balance implementation;
3. construct wall flow contexts from that already evaluated point;
4. distinguish candidate-invariant, evaluation-local and reporting lifetimes;
5. prepare the external-air effective conductance but never freeze the dynamic internal gas film;
6. define staged implementation gates rather than combining cleanup, compilation and numerical changes;
7. require the Debian-13 unchanged-code baseline before attributing future speedups;
8. use stronger production justification for Julia than for a local Python compilation optimization;
9. defer earlier wall extrapolation until post-refactor counters establish that cycle count is still a major limiter.

No material disagreement with Stage 1 was found. Its recommended direction is supported by direct inspection of the present code.

---

## 19. Definition of success

The first acceleration campaign is successful if it produces a solver that is:

- physically equivalent at the existing tolerances and domain rules;
- simpler in its numerical data flow;
- measurably faster on both warm and multi-cycle representative cases;
- instrumented well enough to explain where remaining time goes;
- structured so that a compiled Python kernel or Julia backend can be tested without rewriting campaign/reporting logic.

The immediate goal is **not** to maximize the headline speedup in one patch. It is to remove known waste, preserve confidence in the physics, and expose a clean numerical kernel whose remaining cost can be measured honestly.
