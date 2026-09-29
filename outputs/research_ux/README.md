# Research cockpit UX validation

The demonstration reads the existing local Human Cell Stage 0 campaign without
integration, optimization or journal mutation. It selects original candidates
`a16ed2c7` and `256fdb49` for the existing volume comparison while retaining
campaign-wide funnel and bound evidence.

- [Standalone cockpit](cockpit.html).
- [Browser assertions](browser_checks.json).
- [Visual check](cockpit_check.png).

Nine Chromium checks cover the funnel, all 28 parameter-bound rows, copy buttons,
unique short IDs, selection-driven command updates, both selected candidates,
existing four volume curves, no command-execution links, and the exact command
passed to the clipboard handler. The clipboard API was stubbed for the last
assertion so validation did not alter the user's clipboard.

The targeted suite passed **93 tests** (`pytest -q`): Research cockpit, campaign,
CLI, report and rescale tests. Tests use deterministic fake evaluations for
progress, cache and crash recovery. Actual fd 2 writes, a real failing LSODA
callback, exception identity and descriptor restoration are exercised without
an optimization campaign.

The single final full-suite run passed **782 tests, 0 warnings in 286.82 s**
with `PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 -m pytest -q`.

The existing source/runtime compatibility warning remains visible on the old
campaign: this pass does not silently rewrite execution identity or histories.
Physical science identity remains separate from execution defaults and derived
HTML. See [the implementation contract](../../docs/DADA_ENGINE_RESEARCH_COCKPIT.md).
