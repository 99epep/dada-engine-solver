# Human Cell capacity-scaling validation

Source of truth: the local `examples/human_cell_stage0/campaign` snapshots and
candidate records. No optimization was launched; source files and histories
were not rewritten. The existing untracked `examples/human_cell_stage0/best.html`
was left untouched.

Source A: `a16ed2c7f3fb24dbef2f4e5497e6141af04b17acdb3df106593511e97649243b`.
Source B: `256fdb49f8ee5bdf87212eac149964f40f016ce0874cc8627d473d7bec3a4e95`.

## Measured results

| Case | Cooling [W] | Indicated mechanical input [W] | Cooling COP | Cycles | Model / study |
| --- | ---: | ---: | ---: | ---: | --- |
| A historical | 46.378565 | 26.056771 | 1.779905 | — | valid / feasible |
| A×1, uniform start | 46.381582 | 26.056739 | 1.780023 | 15 | valid / feasible |
| A×5 | 231.885186 | 130.273861 | 1.779982 | 15 | valid / feasible |
| B historical | 53.620433 | 31.355965 | 1.710055 | — | valid / feasible |
| B×5 | 268.295029 | 156.899192 | 1.709983 | 15 | valid / feasible |

The signs of indicated gas power are negative for these refrigerators. Positive
numbers above are indicated mechanical **input**, not useful shaft output or
human power including unmodelled losses.

A×5 has input ratio 4.99961652 and cooling ratio 4.99983528 relative to stored A.
B×5 has input ratio 5.00380682 and cooling ratio 5.00359682 relative to stored B.
Tube lengths/diameters and all constraints were retained. Header packing is
recomputed by the production geometry and is not exact volume similarity.
Maximum pressures are 129089.64 Pa / 170765.58 Pa and maximum temperatures
361.49424 K / 397.78978 K for A×5 / B×5. Tube Reynolds maxima are 266.496 /
294.742 and Mach maxima 0.065253 / 0.049303. Both validity verdicts are `valid`,
with no failed or unavailable model criteria. All declared constraints pass.

Normalized periodic errors are 0.673921 / 0.504180 (required <= 1). Absolute
conservative energy residuals are 2.89e-10 J / -1.29e-10 J. Finite periodic
storage tolerance remains present; no heat/work values were forced to exact
similarity or adjusted after integration.

## Factor-one parity

`../human_cell_A1/study.toml` and its basis reconstruct the selected A physical
configuration exactly. `../human_cell_A1/evaluation.json` and
`../human_cell_A1/source_uniform.json` have exactly equal complete `metrics`
dictionaries on this runtime. Both use uniform initialization. Their small
difference from the historical warm-started campaign record is therefore an
initialization/convergence effect, not scaling physics.

## Artifacts

- [Original A/B with volume curves](comparison_volumes.html).
- [A×5/B×5 with volume curves](comparison_scaled_x5.html).
- [Machine-readable validation](validation_summary.json).
- [Browser interaction checks](browser_checks.json).
- [Volume figure screenshot](volume_plot_check.png).
- [Table screenshot](candidate_table_check.png).
- [A×1 study](../human_cell_A1/study.toml), basis and two evaluation artifacts.
- [A×5 study](../human_cell_A5/study.toml), basis and evaluation artifact.
- [B×5 study](../human_cell_B5/study.toml), basis and evaluation artifact.

Each generated basis contains `provenance.capacity_scaling`, including parent
identities, factor, policies and before/after values. All three studies retain
the original active parameter ownership; `evaluate` uses their new initials.

The original A/B report displays topology as `unavailable`, preserves the
original legacy diagnostic in the detailed evidence, and displays reflux for
both. The new evaluations explicitly report
`wall_integrator_does_not_record_valve_events`. Local reflux remains detected:
the worst signed flows for A×5 / B×5 are -0.002899807 / -0.006735210 kg/s.

Chromium checks exercised both sorting directions, column-header toggling,
active-coordinate sorting, hidden pressure sorting, margin keys, null-last
ordering, four volume curves, unavailable chronology and retained reflux.
The table and curves were also visually inspected. Screenshots focus the
relevant sections; the delivered HTML retains the complete report.

## Reproduction

Use unused destination paths; `rescale`, `evaluate` and `report` refuse overwrite.

```sh
research rescale examples/human_cell_stage0/campaign --candidate a16ed2c7 --factor 1 --output outputs/my_A1/study.toml
research rescale examples/human_cell_stage0/campaign --candidate a16ed2c7 --factor 5 --output outputs/my_A5/study.toml
research rescale examples/human_cell_stage0/campaign --candidate 256fdb49 --factor 5 --output outputs/my_B5/study.toml
research evaluate outputs/my_A5/study.toml --output outputs/my_A5/evaluation.json --budget 3m
research evaluate outputs/my_B5/study.toml --output outputs/my_B5/evaluation.json --budget 3m
research report examples/human_cell_stage0/campaign --candidate a16ed2c7 --candidate 256fdb49 --plots volumes --html outputs/my_comparison.html
```

No shaft-efficiency, pump-power, new liquid correlation or large campaign was
introduced. Cross-study reports do not assign a combined optimization ranking.

Final full-suite result: **760 passed, 0 warnings in 258.43 s** using
`PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 -m pytest -q`.
The old static-cycle assertion was updated to require unavailable chronology
for an empty event sequence, as requested; pressure checks and nonempty-event
classification tests remain intact.
