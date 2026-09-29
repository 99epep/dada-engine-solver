# Human Cell local-refinement validation

This is a bounded search-tool validation, not a completed optimization.
All scientific physical inputs, constraints and original global bounds are
unchanged from the atmospheric-charge source evaluation. Frequency is 2.8 Hz.
No microtube correlation, domain limit, objective or solver tolerance was changed.

## Reproduction

From the source checkout, prefix commands with
`PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src`.
The creation command refuses existing outputs; use a new directory to repeat it.

```sh
python3 -m dada_solver.research refine outputs/human_cell_stage2_reference_charge/replay_old_top5/f0de738d73e0.json --candidate best --radius 0.20 --output outputs/human_cell_stage2c/study.toml
python3 -m dada_solver.research validate outputs/human_cell_stage2c/study.toml
python3 -m dada_solver.research run outputs/human_cell_stage2c/study.toml --directory outputs/human_cell_stage2c/campaign --budget 6m --max-candidates 32
python3 -m dada_solver.research report outputs/human_cell_stage2c/campaign
```

Only this single bounded phase was integrated. It stopped after 303.258 s and
7 attempts: the estimated next evaluation cost exceeded the remaining budget.
The cap of 32 was **not** reached. There was no automatic extension or follow-up
optimization. Snapshot reconstruction was additionally checked without integration:
the next point is basin_1, local Sobol index 6 (campaign sequence index 7).

## Results

| Quantity | Result |
| --- | ---: |
| Attempts | 7 |
| Converged | 6 |
| Feasible, including center | 2 / 7 (28.57%) |
| New local Sobol points | 6 |
| Feasible local Sobol points, excluding center | 1 / 6 (16.67%) |
| Converged but infeasible | 4 |
| Rejected exchanger domain | 1 (hydrodynamic entry domain) |

All generated coordinates lie in the declared clipped local region. The previous
global campaign's user-reported rate was 1/861 (0.116%). These small local counts
are higher, including after excluding the supplied center, but do not establish
a general basin-capture probability or a new optimum.

The center was evaluated first, not copied from the source result. It exactly
reproduced the following stored floating-point metrics:

| Quantity | Center |
| --- | ---: |
| Cooling COP | 1.1220429906193439 |
| Cooling power | 135.56729236333018 W |
| Indicated mechanical input power | 120.82183436527677 W |
| Derived inventory | 0.010357424172130394 kg |

The new feasible neighbor `353205fccae9` has COP 0.7677956131302766,
cooling power 278.43499920037243 W and indicated input 362.6420813544407 W.
It satisfies the **existing declared constraints**; no human-power ceiling or
mechanical efficiency was added. The center remains the best COP candidate.
Its new campaign ID starts `35873e5e84f0`; source provenance retains
`67950d2e6a43986af9b83c9350a659198bdc36d6b870661ef733b8e5f1bcc8ff`.
The new study/search identity explains the different candidate ID.

## Evidence

- [Portable study](study.toml) and [basis](study.basis.json).
- [Offline cockpit](campaign/report.html), including region counts, center,
  objective/COP, rejection rates and candidate origins.
- [Input parity](input_parity.json) and [exact center parity](center_parity.json).
- [Validation summary](validation_summary.json), [CLI validation](validation.txt),
  [run log](run.log) and [text report](report.txt).
- [Test result](tests.log): 267 passed, zero warnings, in 263.68 s, using
  `python3 -m pytest -q tests/test_research*.py tests/test_campaign.py`.
  This includes 22 new local-search tests and the existing rescale, global Sobol,
  persistence, kinematics, exchanger/charge integration and report coverage.
  JavaScript syntax was also checked with `node --check`.

The journal and state remain resumable with this exact runtime/source identity.
No source artifact or historical campaign was modified.
