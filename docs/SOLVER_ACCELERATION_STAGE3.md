# Stage 3: compiled wall-RHS prototype

Date: 2026-09-20. Follows the completed
[Python cleanup](SOLVER_ACCELERATION_IMPLEMENTATION.md) and the compiled-kernel
gates in [Stage 2](SOLVER_ACCELERATION_STAGE2.md).

## Scope and comparison method

This is an explicit research prototype in `examples/numba_wall_prototype.py`,
not a production default or a second campaign implementation. The production
Python source is unchanged in this stage. The prototype compiles one coherent
kernel for the current air/continuum/laminar microtube wall model: variable gas
properties, hydraulic flow, local thermal films, conservative gas/wall balances
and five heat/work quadratures. It uses one Python-to-JIT call per RHS evaluation.

Kinematics remain behind the existing model interface. The prototype uses its
volume/derivative provider and preserves all exact integration breakpoints.
LSODA, integration tolerances, maximum step, periodic tolerances, fixed-cycle
wall acceleration, diagnostics and reports remain the existing Python paths.
The adaptive wall-guess experiment is not combined with compilation.

The port follows `gas_transport.py` (Sutherland/Shomate properties),
`hardware.py` and `hydraulics.py` (pressure-squared Poiseuille, minor losses and
sonic cap), `gas_correlations.py` / `gas_film.py` (local film/domain checks), and
`dynamics.py` / `air_wall.py` (conservative accumulation and angular conversion).
The reference remains authoritative for every unsupported state.

The compiled arithmetic uses float64, `njit`, no `fastmath`, no parallel loop
and no disk JIT cache. Explicit immutable candidate arrays hold parameters;
there is no state/property approximation or interpolation. The first supported
call compiles lazily; subsequent candidates of the same structure reuse the
compiled signature, not a previous candidate's numerical parameters. See the
[Numba JIT documentation](https://numba.readthedocs.io/en/stable/user/jit.html).

The numerical port intentionally supports only exact built-in model classes,
dilute air, no explicit slip model, the present downstream ideal diodes and
variable microtube gas films. Other families, including the historical constant
film case, use the authoritative Python derivative. At an unsupported state
(transport range, rarefaction, transition/turbulence, thermal development or
compressibility guards), the entire RHS falls back to Python. Python retains
all domain decisions, exception types/messages and any exploratory report-only
behavior. Fallback counts are recorded, so a fast result cannot silently omit
costly regimes. No new valve-placement capability is implied.

The standalone context manager patches the RHS only within its benchmark scope
and restores it afterward. It is not suitable as a concurrent production backend.
The experiment duplicates a limited numerical port of the reference equations;
that maintenance cost must be resolved before promoting it into production.

## Validation and measurements

**455 tests pass**, including 25 optional prototype checks. Coverage includes
conservation, valve/reflux and temperature boundaries, invalid-domain errors,
large pressure gradients, candidate ownership, unsupported families, malformed
input and context restoration. The final run is recorded in
[the test artifact](../outputs/solver_acceleration_stage3/pytest.txt).
The compiled RHS was also compared at **1110 saved physical states** covering
four workloads, state extrema and breakpoint neighbors. No sampled supported
state fell back to Python. Differences are below the specified pointwise bound
`1e-11 + 2e-11 * abs(reference_rhs)`; the largest scaled discrepancy is 0.0144
of that bound. See [point checks](../outputs/solver_acceleration_stage3/real_state_checks.json).

Three serial complete-evaluation repetitions use the same frozen Stage-2
manifest, runtime and machine. Python columns below are the existing P6 serial
medians; a fresh single-run Python confirmation is retained separately. No test
suite runs alongside these timing trials. Values include integration and final
reports. Maximum means observed maximum, not a population tail estimate.

| Workload | Clean Python median (s) | Numba median / max (s) | Median speedup | Cycles in both |
| --- | ---: | ---: | ---: | ---: |
| Production warm | 6.387 | 4.554 / 7.830 | 1.40× | 1 |
| Nearby compatible warm | 12.122 | 4.375 / 4.591 | 2.77× | 3 |
| Production cold | 68.913 | 8.040 / 8.425 | 8.57× | 22 |
| Smooth four-bar | 85.760 | 13.802 / 14.448 | 6.21× | 32 |
| Historical, Python fallback | 1.766 | 2.036 / 2.128 | 0.87× | 1 |

Sources: [timings](../outputs/solver_acceleration_stage3/numba/summary.json),
[comparison](../outputs/solver_acceleration_stage3/comparison.json),
[backend metadata](../outputs/solver_acceleration_stage3/numba/prototype.json).
The [fresh Python confirmation](../outputs/solver_acceleration_stage3/python_confirmation/summary.json)
returns 6.803 s warm, 12.842 s nearby, 70.025 s cold, 92.792 s smooth and
2.024 s historical, with [exact reference outputs](../outputs/solver_acceleration_stage3/python_confirmation_equivalence.json).
The historical Numba-scope median of 2.036 s is essentially the same as this
fresh Python run: the older 1.766 s median does not establish a wrapper-overhead
cause by itself. Single-run confirmation timings are not substituted for the
three-repeat baseline medians.
The observed latency advantage is restricted to supported workloads; the
prototype does not meet a universal first-call/unsupported-family latency gate.

The warm maximum includes compilation; only its second and third runs are fully
warmed (4.554 and 4.200 s). The first supported RHS call took **2.933 s** including
lazy compilation. Backend import took 0.251 s, separate from the runner's common
0.913 s import time. These are measured in-process components, not a complete
operating-system process-startup benchmark. No disk cache is used, so a fresh
process pays compilation again. Preflight compilation is not cooperatively
interruptible; warmed interruption semantics remain unchanged.

All four compiled workloads retain status, cycle count and physical validity
classifications. Final state differences are at most 0.022 of the original
periodic comparison scale. Maximum relative differences against Python are
3.63e-7 for indicated power, 3.04e-7 for indicated thermal efficiency, and
1.91e-14 for total mass. Invalid-domain exception/message and interruption status
are preserved. The unsupported historical case retains exact states and outputs.
No supported valid trajectory used the fallback, across all three repetitions;
the invalid-temperature case deliberately used it once per repetition.

These compiled results are **not bitwise identical** to Python. Ordinary
floating-point rounding changes LSODA's accepted grid and counters. For example,
warm `nfev` increases from 6080 to 7059, yet integration becomes much faster.
Tolerances, maximum step, physical equations and periodic acceptance are unchanged.
Periodic state tolerance is not a rigorous bound on output error; a separate
tighter-periodic comparison is retained as an accuracy check.

The common tighter-periodic experiment divides both periodic tolerances by 100,
leaving the ODE tolerances unchanged, and again starts from the frozen initial
states. One run per method gives:

| Workload | Python cycles / time (s) | Numba cycles / time (s) |
| --- | ---: | ---: |
| Production cold | 33 / 107.659 | 33 / 14.099 |
| Smooth four-bar | 52 / 128.449 | 52 / 22.263 |

The first Numba cold run includes compilation. At the original comparison scale,
state distances are 0.0107 and 0.0449 or less. Relative indicated-power differences
remain below 2.54e-7 and efficiency differences below 3.05e-7. Both meet the stated
original 1e-6 relative output-comparison threshold; validity classifications agree.
This is supporting evidence at tighter periodic convergence, not a proof of
absolute ODE/discretization accuracy or identical fixed points at 1e-8. See
[the tighter comparison](../outputs/solver_acceleration_stage3/tight_comparison.json).

The remaining cost is now workload-dependent. In a warmed run, integration takes
0.262 s of a 4.554 s warm evaluation; diagnostics/validity take about 4.29 s.
Cold integration takes 4.26 s, and smooth-motion integration 11.19 s. The generic
Python kinematics interface is intentionally retained, so smooth mechanisms have
more Python work per compiled call. Compiling reporting or replacing the ODE
algorithm is not included in this experiment.

## Adoption decision

**Go for integration work on the supported multi-cycle workloads:** measured
complete-evaluation gains comfortably exceed the Stage-2 1.5× prototype gate.
**No automatic replacement of the production default:** one-cycle warm work has
a smaller gain, fresh-process compilation can increase latency, and unsupported
families gain nothing. The independently ported equations and benchmark-only
patching interface also need a production design before campaign adoption.

The usable deliverable here is the isolated, reproducible prototype and its
validation evidence. Production defaults, campaign identities/histories and the
physical reference remain unchanged. Promoting this kernel should provide an
explicit backend selection owned by campaign identity, include the Numba/runtime
identity, avoid process-global patching, and prevent the reference and compiled
formulas from drifting apart. Julia is not needed to establish this gain and
has not been prototyped in this stage.

## Reproduction

Numba is optional and already installed in the measurement environment. It is
not added as a mandatory package dependency. Optional prototype tests skip when
Numba is absent.

```sh
PYTHONPATH=src python3 -m pytest
PYTHONPATH=src python3 examples/check_numba_wall.py
PYTHONPATH=src python3 examples/benchmark_numba_wall.py \
  --backend python --output outputs/solver_acceleration_stage3/new_python \
  --repeats 3
PYTHONPATH=src python3 examples/benchmark_numba_wall.py \
  --backend numba --output outputs/solver_acceleration_stage3/new_numba --repeats 3
PYTHONPATH=src python3 examples/compare_numba_wall.py \
  outputs/solver_acceleration_stage3/new_python \
  outputs/solver_acceleration_stage3/new_numba \
  outputs/solver_acceleration_stage2/manifest.json
```

Use `--manifest outputs/solver_acceleration_stage2/manifest_tight_periodic.json`
and `--cases production_cold smooth_four_bar --repeats 1` to reproduce the
separate tighter-periodic replay. Keep the original manifest as the third
argument to `compare_numba_wall.py` when applying the documented original
accuracy comparison scale. Both replay directories must have the same numerical
manifest hash; the comparison script exits nonzero on a failed stated gate.

The comparison requires matching cases in both directories; select matching
`--cases` when creating a restricted baseline. All measurements are serial;
run no tests or other simulations alongside timing trials. Output directories
must be new. Frozen inputs are hash-checked by the Stage-2 runner. Each directory
contains environment/source identities, a prototype source hash, per-run reports,
trajectories, solver statistics, section timings and the additional
`prototype.json` compilation/fallback metadata.
