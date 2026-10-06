# Research cockpit and execution UX

This is the current reference for terminal progress, stored-result inspection,
derived reports and interruption recovery. For everyday setup and commands, use
the [user guide](DADA_ENGINE_RESEARCH.md). Study schema, scientific identity,
parameter ownership and search semantics belong to the
[technical reference](DADA_ENGINE_RESEARCH_REFERENCE.md).

## Commands and candidate inspection

```sh
research status outputs/my_study/campaign
research report outputs/my_study/campaign
research report outputs/my_study/campaign --candidate best
research report outputs/my_study/campaign --plots none
research compare outputs/a/campaign outputs/b/campaign
research resume outputs/my_study/campaign --budget 30m
```

`research` is the source-checkout shorthand defined in the guide; the installed
entry is `dada-research`. These paths are illustrative destinations.

`status` inspects stored evidence without thermodynamic integration. Terminal
inspection normally shows a summary and the best feasible result rather than
every candidate. `--list-candidates` expands the selected candidate listing;
`--json` returns the complete structured inspection dataset, without the HTML
detail cap. Use `--plots none` with `report` to avoid curve reconstruction.

Selectors accept `best`, `second`, a complete ID or an unambiguous prefix;
repeat `--candidate` for multiple selections. Unknown or ambiguous selectors
fail explicitly. Comparison resolves selection across the supplied sources but
does not invent a common objective ranking for incompatible studies.

## Terminal progress

Each completed evaluation produces one permanent line with status, search origin,
useful metrics or rejection reasons, and `BEST` when the feasible objective
improves. Cached attempts are labelled. Metrics distinguish motor indicated gas
power from refrigeration cooling power and indicated input.

On a TTY, a separate transient status line displays attempts, elapsed time,
feasibility and best-result information. It is cleared before permanent lines
and cleaned on completion or interruption. Non-TTY output retains permanent
lines and start/finish summaries without ANSI refresh sequences or periodic
status spam.

The runner sends progress snapshots at start and finish, for completed evaluations,
and at solver progress boundaries. Ordinary periodic notifications are throttled:
within ten seconds of the previous notification they require at least ten more
attempts. Evaluation notifications are not subject to that throttle.
An evaluation that returns to solver progress boundaries can therefore report
progress before finishing; an uninterrupted native call cannot.
This is synchronous reporting, not a monitoring thread or parallel scheduler.

## Reports, curves and replay

`report` regenerates `CAMPAIGN/report.html` by default. For a standalone evaluation,
the default is a sibling `.html`; `--html PATH` overrides it. `compare` only writes
HTML when requested with `--html`. Reports are standalone/offline artifacts with
no server or network dependency. Regeneration does not change stored scientific
verdicts or metrics.

Default report curves are positions, gas-to-wall heat exchange, pressures and
temperatures for the two best feasible candidates. Explicit selectors change the
targets; `--plots` chooses curves or supported mechanism animations, and
`--plots none` omits them. Comparison does not request curves by default.
See the [guide examples](DADA_ENGINE_RESEARCH.md#reports-curves-and-animations).

Thermodynamic curves can require a separately announced one-cycle replay from the
saved periodic endpoint. This does not restart optimization or full periodic
convergence. Missing state or replay failure makes the curve unavailable, with an
explanation, without changing the historical result.

The derived `.research-plot-cache` beside the HTML avoids identical thermal
reconstructions. Its identity includes the saved state and current runtime.
A reconstruction under another checkout/runtime is new derived evidence, not the
exact historical trajectory. Runtime compatibility and replay metadata remain
visible. Kinematic sampling and local mechanism-frame conventions are documented
in the [kinematics reference](DADA_ENGINE_RESEARCH_KINEMATICS.md).

## Deliberate HTML detail limit

The presentation policy `best_distinct_feasible_tenth_v1` retains detailed records
for the best distinct feasible candidates up to `ceil(0.10 * attempts)` for
compatible comparisons. Required plot targets and explicit selections remain
included even outside that base. In incompatible comparisons the renderer retains
feasible distinct selections up to the cap without imposing a joint ranking.

Global statistics, the funnel and progression still cover all attempts. Complete
scientific records remain in the journal and structured JSON inspection. This is
a presentation policy, not deletion of scientific data or search sampling.
HTML filters operate on the embedded selection: an empty filtered table does not
prove that the full journal contains no matching candidate.

## Cockpit evidence and review actions

The cockpit presents an attempts/integrated/converged/feasible funnel, cache-hit
and distinct-candidate counts, failure categories and violated constraints.
These are stored-outcome counts; a cached outcome does not imply fresh integration.
Original reasons and available constraint values remain the underlying evidence.

Bound pressure counts distinct coordinates within the outer 5% of each normalized
numerical interval, including rejected candidates. Choice categories are excluded
from numerical bound-pressure interpretation. Up to five distinct feasible elites
are selected using the existing objective. At least three elites with 60% near a
bound trigger a bound-review signal. Otherwise fewer than 512 distinct candidates
or five elites can expose a continuation suggestion. These are UX conventions,
not evidence of saturation, optimality or periodic convergence.

Actions have distinct meanings:

- Shell commands, such as `resume` or report regeneration, are copyable suggestions
  for the operator; the HTML never executes them.
- A `bounds` action displays an inspection button that scrolls to bound pressure.
- A `filter` action displays an inspection button selecting a failure category or
  violated constraint, scoped to its source.

The renderer does not replace `bounds` or `filter` actions with report commands.
Review actions never change bounds, constraints or model domains automatically.
Cross-study evidence remains separated by source.

## Journal-first persistence

New campaigns use `history.jsonl.gz`. Each appended result is an independent gzip
member containing one canonical JSON record terminated by a newline. Compression
is lossless and uses a deterministic gzip timestamp.

A completed evaluation, including a cache hit, follows this durable order:

1. Atomically write and fsync `recovery.json`.
2. Append the result to the journal and fsync it.
3. Atomically update and fsync `state.json`.
4. Remove `recovery.json` after durable state acknowledgement.

Directory updates are also synchronized. A verified recovery result absent from
the journal is appended exactly once; a matching already-journaled result is not
duplicated. Conflicting records or invalid identities fail explicitly. Pending
work is reconciled with durable completions before a new evaluation.

Read-only inspection can include unjournaled recovery without changing storage.
Execution resume performs repairs under the exclusive writer lock. An incomplete
last gzip member can be preserved as `history_torn_tail_*.bin` and truncated back
to the last complete member. CRC corruption of a complete member is an error,
not an interrupted append; non-final corruption is not silently repaired.

Campaign persistence requires `history.jsonl.gz` and the single `recovery.json`
completion slot. Plain journals and per-candidate recovery layouts are rejected.

Interruptions preserve pending work for controlled resume. `--retry-incomplete`
explicitly retries unresolved deadline-limited work. Execution compatibility is
checked against stored snapshots, while offline inspection may remain available
when resume is incompatible. See the [technical reference](DADA_ENGINE_RESEARCH_REFERENCE.md)
for the complete cache and compatibility contract.

## Native solver stderr

The integration wrappers capture native solver stderr around the relevant solver
calls. Captured text can remain in technical diagnostics without changing the
scientific status or reason. The descriptor is restored even after an exception
or `KeyboardInterrupt`; normal Python error reporting remains available.

This is a POSIX, process-wide facility. Wrapped calls are serialized, but unrelated
threads must not write to stderr during a solve. Capture is not a change to solver
arguments, tolerances, validity limits or retry rules.

## Execution limits are not scientific criteria

Generated and capacity-rescaled studies default to 512 attempts per invocation.
The CLI interprets the historical configured default 16 as 512 without rewriting
snapshots; explicit `--max-candidates 16` still requests 16. Other configured caps
remain effective. This is not a total campaign-size limit.

The attempt cap, invocation time budget and individual evaluation deadline have
different scheduling roles. None establishes physical validity or convergence.
Capacity transformation rules belong to the [capacity reference](DADA_ENGINE_RESEARCH_CAPACITY.md);
current evidence boundaries and known failures belong to [validation](validation.md).

Terminal productivity labels use `Q/N` and `COP·Q/N`, with `W/µt` meaning watts
per physical microtube; `N` counts both exchangers. HTML temperature curves show
only SMALL/LARGE gas temperatures. Candidate comparisons list active parameters
first and report total exchanger gas hold-up (both exchangers), not external
packaging volume.
