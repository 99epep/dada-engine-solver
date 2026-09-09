# Proposed next steps after the motor update

Assessment dated 2026-09-08. These are recommendations, not implemented physics.

## Interpreting the existing diagnostics

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

## Exchanger priority

The current code already provides channel geometry, thermal resistances,
pressure-drop closures, an external liquid-side screening model, and a fixed-
point cycle/geometry calculation. These should be reused.

The main unresolved step is predicting transfer over the actual pulsating cycle.
The fixed-point calculation currently combines an adjacent-port peak flow with
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

## Speed and language choice

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
error tolerances. A Julia prototype is reasonable, especially for a future
larger exchanger state, but needs to include the right-hand side and event
handling rather than retaining frequent Python callbacks. Include compilation
time separately from repeated-solve time. No speedup factor is established yet.

SciPy's [LSODA interface](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)
already wraps a Fortran integrator. SciML's
[ODE optimization guidance](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/faster_ode_example/)
likewise emphasizes allocations, right-hand-side implementation and solver
selection. A language change can improve throughput; it cannot validate the
exchanger heat-transfer closure.
