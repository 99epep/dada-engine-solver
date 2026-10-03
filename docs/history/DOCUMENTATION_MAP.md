# Documentation consolidation map

Consolidated on 2026-10-03 from the local working tree. The pre-consolidation
tracked source is local Git commit `0baf7fd6a79a12e05793070fd88b8ed65538344a`.
This is provenance, not a dependency on GitHub. To read a retired document in a
clone retaining this history:

```sh
git show 0baf7fd6a79a12e05793070fd88b8ed65538344a:docs/SOLVER_ACCELERATION_STAGE4.md
```

The standalone documentation preserves the adopted contracts, significant
negative results, primary-source references and artifact locations. Intermediate
handoff prompts, duplicate pending plans, repeated test-count snapshots, profiler
listings and withdrawn numerical sweeps were removed from current guides.
The complete physics decision ledger is retained because the reasons for past
decisions and withdrawals remain scientifically useful.

## Merged notes

| Former path under `docs/` | Current reference |
|---|---|
| `SOLVER_ACCELERATION_STAGE1.md` | [SOLVER_PERFORMANCE.md](../SOLVER_PERFORMANCE.md) |
| `SOLVER_ACCELERATION_STAGE2.md` | [SOLVER_PERFORMANCE.md](../SOLVER_PERFORMANCE.md) |
| `SOLVER_ACCELERATION_IMPLEMENTATION.md` | [SOLVER_PERFORMANCE.md](../SOLVER_PERFORMANCE.md) |
| `SOLVER_ACCELERATION_STAGE3.md` | [SOLVER_PERFORMANCE.md](../SOLVER_PERFORMANCE.md) |
| `SOLVER_ACCELERATION_STAGE4.md` | [SOLVER_PERFORMANCE.md](../SOLVER_PERFORMANCE.md) |
| `SOLVER_ACCELERATION_STAGE5.md` | [SOLVER_PERFORMANCE.md](../SOLVER_PERFORMANCE.md) |
| `PERIODIC_MAP_DIAGNOSTIC.md` | [SOLVER_PERFORMANCE.md](../SOLVER_PERFORMANCE.md) |
| `PERIODIC_ANDERSON_EXPERIMENT.md` | [SOLVER_PERFORMANCE.md](../SOLVER_PERFORMANCE.md) |
| `AIR_SOURCE_EXCHANGERS.md` | [history/EXCHANGER_TRIALS.md](../history/EXCHANGER_TRIALS.md) |
| `AIR_WALL_COUPLING.md` | [history/EXCHANGER_TRIALS.md](../history/EXCHANGER_TRIALS.md) |
| `GEOMETRIC_MOTOR_COUPLING.md` | [history/EXCHANGER_TRIALS.md](../history/EXCHANGER_TRIALS.md) |
| `EXCHANGER_NEXT_STEPS.md` | [history/EXCHANGER_TRIALS.md](../history/EXCHANGER_TRIALS.md) |
| `WORKING_FLUIDS.md` | [WORKING_FLUID_MODELS.md](../WORKING_FLUID_MODELS.md) |
| `DOUBLED_EXCHANGER_TRIAL.md` | [history/EXCHANGER_TRIALS.md](../history/EXCHANGER_TRIALS.md) |
| `EXCESS_AIR_TRIAL.md` | [history/EXCHANGER_TRIALS.md](../history/EXCHANGER_TRIALS.md) |
| `HIGHER_TEMPERATURE_TRIAL.md` | [history/EXCHANGER_TRIALS.md](../history/EXCHANGER_TRIALS.md) |
| `PARALLEL_EXCHANGER_TRIAL.md` | [history/EXCHANGER_TRIALS.md](../history/EXCHANGER_TRIALS.md) |

Four trial filenames explicitly referenced by the project instructions remain
as short forwarding notes. Other retired files have no production-code/test
references and were removed after documentation links were updated.

## Rewritten references

- Root README: current capabilities and a short route into the documentation.
- `PHYSICS_DECISIONS.md`: current contracts; original text moved to
  [the historical ledger](PHYSICS_DECISION_LEDGER.md).
- `validation.md`: verification map and dated current check record, without
  obsolete optimization/profile recommendations.
- `MOTOR_RESEARCH_OBJECTIVES.md`: ordered objectives and links to completed work,
  replacing the stale “temperature map planned” instructions.
- `COOLING_CELL.md`: external load and power boundary; withdrawn results removed
  from the guide and kept in the ledger.
- Model/campaign/sizing references: obsolete rollout plans and validity claims
  corrected against current code; compatibility and source evidence retained.

## Deliberately retained

Mechanism methodology and family geometry, correlation equations and citations,
license notices, historical source artifacts, TOML/JSON/CSV inputs, example code,
tests and all campaign histories are unchanged in substance. Dated study reports
are visibly labelled so their constraints/results cannot become universal defaults.
No physical calculation or optimization was run for this consolidation.
