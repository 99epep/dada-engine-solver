# Exchanger design lessons and next steps

## Current status and design lessons (2026-09-17)

The implemented motor path now couples microtube geometry, hydraulic passages,
dynamic wall storage and the opt-in production variable-property internal gas
film. See [MICROTUBE_GAS_MODEL.md](MICROTUBE_GAS_MODEL.md) for implemented
closures and [EXCHANGER_VALIDATION.md](EXCHANGER_VALIDATION.md) for their evidence
and limitations. Quasi-steady pulse response is not experimentally validated;
there is no empirical pulsating heat-transfer multiplier.

The [four-stage experiment chronology](FOUR_STAGE_OPTIMIZATION.md) separates
legacy-gas-model controlled tests from the production local 9D campaign and
links each conclusion to its exact report. Preserve these design lessons:

1. Multiplying thermal conductance at fixed dead volume and wall capacity is a
   sensitivity experiment, not a realizable redesign. The saturation study
   multiplied both gas-wall and air-wall conductances at zero geometric cost.
2. Physical tube-count changes simultaneously alter conductance, flow area,
   hydraulic loss, tube hold-up, headers and wall thermal capacity. These
   geometry-derived quantities must not become independent free coordinates.
3. A larger physical H_o can help shorten LP exchange yet reduce efficiency,
   even where an artificial multiplier improves it. The 120/60 and 150/30
   comparisons demonstrate the coupled tradeoff; they do not separately identify
   every storage/volume contribution.
4. H_i and H_o need not have equal best-found geometry. The production 9D
   checkpoint uses a smaller H_i and slightly larger H_o than its reference.
5. Static UA is only reference metadata in the production model; the actual
   internal gas-side conductance follows instantaneous flow and gas properties.
6. External-air fan power is excluded by explicit project decision, not
   physically zero. The objective is indicated thermal efficiency, not useful
   shaft or complete-system efficiency.
7. External-air-side modelling is still a screening model and a major future
   hardware-design uncertainty. Improving the internal film does not validate
   the external-air passage or its predicted fan requirements.

Compression and expansion retain gas-wall heat exchange. The phase ablations
measure local sensitivity of frozen cycles; they cannot establish global
benefit when motion and compression magnitude can also be redesigned.

The next study is the [whole-machine temperature map](MOTOR_RESEARCH_OBJECTIVES.md),
with fixed total swept volume, variable S/L ratio and candidate-specific 100 kPa
filling. Each temperature receives hardware and then motion adaptation. Preserve
full domain verdicts (including supported turbulent states) rather than importing
historical laminar-only campaign guards. No new campaign was run for this audit.

## Historical assessment (2026-09-08)

The following assessment preserves the earlier rectangular-channel/fixed-point
workflow and profiling evidence. It is not the current microtube architecture
or a fresh authorization for model changes. Dynamic-wall coupling and internal
variable-property closures have since been implemented as described above.

### Interpreting the diagnostics in the earlier assessment

The conservative solver already integrates variable temperatures and separate
cylinder/exchanger pressures. Failure to meet the configured isothermality or
pressure-equalization target therefore does not itself invalidate those balance
equations. Separate three kinds of result in a future reporting revision:

- numerical reliability: periodic convergence, conservation, refinement;
- model applicability: fluid properties, phase changes, correlation domains,
  omitted hydraulic/thermal dynamics;
- design behavior: temperature excursion, nominal phase order, useful power,
  efficiency, pressure and force requirements.

The existing isothermality metric is `(T_max-T_min)/T_reservoir`. It measures
excursion, not the mean gas/reservoir temperature difference or lost work.
The subsequent user decision makes it a diagnostic only for the current work.
It no longer fails the core validity verdict. Historical explicitly configured
isothermality sizing constraints remain opt-in; none is added to the motor work.

Non-nominal valve order is likewise a behavior diagnostic. Ideal diodes prevent
reverse valve flow by construction, which does not validate actual valve lift,
inertia or leakage. Additional valve dynamics are justified if event/flow
sensitivity materially changes useful performance or hardware loads.

The current power, volume, speed and temperature requirements are recorded in
[MOTOR_DEMONSTRATOR.md](MOTOR_DEMONSTRATOR.md).

### Earlier exchanger priority

At that assessment, the code provided channel geometry, thermal resistances,
pressure-drop closures, an external liquid-side screening model, and a fixed-
point cycle/geometry calculation. These should be reused.

The main unresolved step is predicting transfer over the actual pulsating cycle.
The historical fixed-point calculation combines an adjacent-port peak flow with
representative pressure/temperature extrema and sends a constant UA back to the
0D cycle. That is a screening approximation; the selected extrema need not occur
simultaneously. A more detailed geometry calculation alone does not remove it.

User-confirmed sequence (2026-09-08):

1. Export simultaneous cycle pressure, temperature and port-flow histories.
   Specify an **intermittent unidirectional** intended flow, cycle frequency,
   pulse duration and active duty; do not substitute zero-mean oscillatory flow.
   Any local reflux calculated by existing bidirectional links is a reported
   discrepancy, never silently removed or used to redefine the intended flow.
2. Compare published/manufacturer evidence for plate-fin and parallel microtube
   designs. The motor search is not constrained to low-tech fabrication.
3. Use qualified experimental/manufacturer maps with explicit fluid, geometry,
   boundary conditions, uncertainty and validity range. Under the subsequent
   user decision, use explicitly documented Doty extrapolation if adequate
   commercial data are unavailable; distinguish estimates from measurements.
4. Before any dynamic exchanger extension, briefly explain its benefit and
   computational cost to the user. No such extension is authorized by this list
   alone without that explanation.

Motor efficiency is the primary design objective. Ideal chronology remains the
aspirational design target. Isothermality is deferred and has no validity veto;
its existing excursion diagnostic remains available. See
[EXCHANGER_LITERATURE.md](EXCHANGER_LITERATURE.md) for the first duty envelope and
source comparison. The transient-model-first ordering above has been superseded.

### Current acceleration review

The [stage-1 solver-acceleration note](SOLVER_ACCELERATION_STAGE1.md) profiles
the production variable-property motor and separates integration, periodic
convergence and duplicated diagnostics. It supersedes the old profile below
as the starting point for performance discussion, without approving a rewrite.
An independent second review will precede implementation.

### Historical profiling and language discussion

The current motor example was profiled using:

```console
PYTHONPATH=src python3 -m cProfile -o /tmp/dada-motor-profile.pstats \
    -m dada_solver examples/motor_controlled_example.toml
```

In this instrumented run, the periodic solver took approximately 15.0 s over
five cycles and 31,766 right-hand-side evaluations. About 9.6 s were spent in
the right-hand-side path and 3.0 s in continuous-diode event reconstruction.
These cumulative timings have profiling overhead and are not a Julia comparison.

First remove repeated state/geometry reconstruction and unnecessary allocations
in those measured paths, preserve the regression cases, and reuse periodic
solutions between nearby design points. Process-level parallel evaluation is
appropriate for independent candidates; shared LSODA instances should not be
assumed thread-safe.

Then benchmark a compiled numerical kernel on representative cases at matched
error tolerances. A Julia prototype was discussed as a possible later comparison, not approved
or implemented. Any such comparison would need the right-hand side and event
handling, with compilation time separated from repeated-solve time. No speedup
factor is established. The current decision is to discuss migration separately
after the temperature experiment.

SciPy's [LSODA interface](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)
already wraps a Fortran integrator. SciML's
[ODE optimization guidance](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/faster_ode_example/)
likewise emphasizes allocations, right-hand-side implementation and solver
selection. A language change can improve throughput; it cannot validate the
exchanger heat-transfer closure.
