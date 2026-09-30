# Research terminal UX demonstration

This is a **synthetic display test**, not a physical Human Cell study.
`demo.py` runs the production campaign runner with four artificial evaluator
results and a simulated clock (52 seconds). No thermodynamic integration or
optimization was run. Source provenance hashes label synthetic fixtures.
No existing Human Cell campaign was read, stopped or modified by this demo.

The actual runner uses two center-first local regions and produces:

```text
[0001] basin_1:center  feasible  BEST  COP 1.1472  Qcold 151.8 W  Pin 132.3 W  mdot 0.061
[0002] basin_2:center  converged_infeasible  maximum_absolute_mass_flow  mdot 0.084 > 0.080
[0003] basin_1:sobol  invalid_exchanger  large_relative_pressure_drop
[0004] basin_2:sobol  feasible  BEST  COP 1.1892  Qcold 286.4 W  Pin 240.8 W  mdot 0.078
```

On the 80-column PTY used for this check, the last transient line is:

```text
4/4 · 52s/5m · feasible 2 · COP 1.1892 · best 397c54c0d1b4 · Qcold 286.4 W
```

The indicated input and convergence count are omitted from this transient line
because of terminal width. They remain available in the permanent evaluation
and final summary. The budget is 5 minutes and the test explicitly stops at
four evaluations. No sleeps are used to produce progress events.

- `run.log`: readable non-TTY output, without carriage returns or ANSI escapes.
- `terminal_capture.txt`: raw interactive writes, intentionally including terminal
  controls. Inspect as escaped text rather than printing it into an active shell.
- `events.json`: the actual runner callback stream.
- `check.json`: capture assertions: four evaluation lines, six permanent lines
  including start/finish, and transient status fitting within 79 display cells.
- `campaign/`: isolated synthetic history and state, using normal persistence.

The terminal and log sinks received the same callbacks during one runner
execution. For reproduction, copy `demo.py` to a new empty directory and run it
from the checkout with `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 .../demo.py`.
Its initializer refuses overwriting an existing study.

Scheduling, Sobol allocation, candidate payloads, cache and persistence code paths
are unchanged. The normal strict source/runtime fingerprint check for execution
resume remains in force; this display change does not bypass that check.

Validation: 67 targeted tests passed; the complete relevant Research/campaign
suite passed 290 tests in 244.62 s, with zero warnings. See `targeted_tests.log`
and `research_tests.log`. Tests exercise terminal refresh/width, permanent lines,
scientific reasons and limits, CLI interruption cleanup, and unchanged candidate
IDs/ordering with and without presentation for global and local searches.
