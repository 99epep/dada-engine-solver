# Research cockpit and execution UX

This pass changes presentation, execution defaults and durable storage only.
Physical equations, numerical tolerances, scientific study inputs, candidate
payload construction and Sobol generation are unchanged. Execution options do
not enter the scientific study identity. The existing source/runtime identity
check still applies to resume; incompatible source versions are not silently
accepted or rewritten.

## Short commands

```sh
dada-research report outputs/my_study/campaign
dada-research resume outputs/my_study/campaign --budget 30m
dada-research report outputs/my_study/campaign --candidate a16ed2c7 --plots volumes
```

`report CAMPAIGN` writes `CAMPAIGN/report.html`. For a standalone evaluation
JSON it writes the sibling `.html`. `--html PATH` overrides the destination.
HTML is a derived artifact and may be regenerated in place; studies, bases,
evaluation JSON and histories retain their existing protection rules.

Normal report/status output contains the campaign summary and one best feasible
result, not every candidate. Report output also identifies the HTML path.
`--list-candidates` explicitly lists selected candidates; `--json` prints the
complete inspection dataset as valid JSON. A candidate selector accepts `best`,
a full ID or an unambiguous prefix. Comparison resolves prefixes across the
supplied sources. Unknown or ambiguous prefixes fail explicitly.

Generated and capacity-rescaled Research studies now set
`execution.default_max_candidates = 512`. Tests and smoke studies may explicitly
set a different cap. The CLI interprets the historical default of 16 as 512
without editing existing snapshots. `--max-candidates 16` remains an explicit
16-attempt override. Other explicitly configured caps remain effective. The cap
limits new attempts in this invocation, not total campaign size. Budget and
per-candidate deadlines remain unchanged.

## Lightweight progress

The runner emits a side-channel snapshot at start, about every ten seconds or
ten completed attempts, on objective improvement, and at finish. Solver progress
boundaries also allow updates while a candidate is still integrating. There is
no timer thread and no report recomputation. An uncooperative long native call
cannot refresh until it returns to a progress boundary.

The snapshot contains phase/index, attempts/cap, elapsed/budget, converged and
feasible counts, principal failure statuses and the best feasible objective.
Refrigeration metrics include COP, cold power and indicated input when present.
The best may come from an earlier phase; counters describe the current phase.
TTY output updates with carriage return and clears the previous line; non-TTY
output uses periodic complete lines. Completion emits one compact summary.

## Native solver stderr

`solver_output.capture_native_stderr` duplicates POSIX fd 2 and redirects it to
a temporary file immediately around each `solve_ivp` call. Python and C buffers
are flushed before redirection/restoration. The saved descriptor is restored
and closed in `finally`, including callback exceptions and `KeyboardInterrupt`.
A file avoids bounded-pipe deadlocks. The solver receives identical arguments,
and the original Python exception object/type/message is re-raised.

Captured text is attached as `native_solver_stderr` to a failing exception or a
solver result. Wall-evaluation failures can retain it in `technical_diagnostics`;
the scientific failure status and reason remain unchanged. It does not suppress
Python exception reporting after the descriptor is restored. The generic gas
and external-wall integration paths share this wrapper. No tolerance, validity
rule, retry policy or solver algorithm is changed.

Descriptor replacement is process-wide, so these solver calls are serialized
with a reentrant lock. Research already evaluates sequentially. Unrelated
threads must not write to fd 2 during a solve. This is a POSIX facility, matching
the existing filesystem locking contract.

## Cockpit evidence

The HTML adds:

- A funnel of attempted, integrated, converged and feasible stored outcomes,
  cache hits and distinct candidates. Counts describe records, including cache
  outcomes; they do not imply a new integration was performed on cache hits.
- Stable rejection categories, with exact original reasons retained in records.
  Transport temperature failures use the declared transport domain to classify
  below/above-domain cases. When that evidence is unavailable, the category
  remains `transport_temperature_outside_domain`; no cutoff is guessed.
- For each active parameter and each bound, the number/percentage of distinct
  candidate coordinates within the outer 5% of the normalized interval, plus
  elite counts. This includes rejected coordinates and honors log transforms.
  The latest record per candidate ID is used. Elites are up to five distinct
  feasible candidates ranked by the existing minimized objective.
- Deterministic factual review actions and copyable commands. At least three
  elites with 60% near a bound marks pressure and takes precedence over a
  continuation prompt. Otherwise fewer than 512 distinct points or five elites
  exposes a continuation action. These are presentation thresholds, not search
  saturation or scientific convergence criteria. Dominant domains and violated
  constraints are surfaced without proposing relaxed limits.

Scientific domains and bounds are never edited automatically. Selected-candidate
volume commands update when either comparison selector changes. IDs use a short
prefix checked against the whole source campaign. Paths and arguments use POSIX
shell quoting. Copy uses the clipboard API with a selection-copy fallback. The
standalone HTML executes no shell command and needs no server or network.
Cross-study reports retain separate funnels, bounds and commands per source.
Existing sorting, filters, margins, topology/reflux and volume graphs remain.

## Journal-first persistence

New results no longer produce `candidates/<sha256>.json`. Instead:

1. Atomically write and fsync `recovery.json`, including its directory rename.
2. Append the complete result to `history.jsonl` and fsync it.
3. Atomically update and fsync `state.json` with its existing schema.
4. Remove `recovery.json` and fsync the directory.

The runner clears recovery only after the state commit, including when resuming
from a recovered result. Cache-hit attempts follow the same durable sequence.
A recovery record absent from the journal is verified and appended exactly once.
A residual recovery matching a journal record is not duplicated. Conflicting
records or invalid identities fail explicitly. The pending Sobol point is
reconciled with durable completed records before any new evaluation.

A torn final append is preserved in the existing `history_torn_tail_*.bin`
artifact, truncated and recovered from the completion slot. Non-final
corruption remains an error. Read-only inspection includes unjournaled recovery
without modifying the journal, state or recovery file. Execution resume repairs
under the existing exclusive writer lock.

Legacy `candidates/*.json` files remain readable and recoverable. No automatic
migration or deletion occurs, and new evaluations in those campaigns also use
the single recovery slot. Exact cache, original IDs, warm states and Sobol order
are retained. A large campaign now has one growing result journal, one transient
completion file and existing per-phase summaries, rather than a file per result.
