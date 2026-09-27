# Dada-Engine Research validation

Validation date: 27 September 2026. Steps 1–5 of the
[implementation plan](DADA_ENGINE_RESEARCH_IMPLEMENTATION_PLAN.md) are implemented
for the fixed independent six-bar pair and five thermal/hardware parameters.
This validates the first study; it does not establish an optimum or complete
the user's required researcher usability review.

## Checks

- Full source-checkout suite: **595 passed**, 94.35 s. The 964 warnings were
  existing NumPy `trapz` deprecations. Command:
  `PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 -m pytest -q`.
- After the final report-only change to distinguish close axis tick values,
  the 30 schema, CLI and report tests passed again in 3.87 s.
- Physical checks cover the stored reference and two changed configurations
  against the original evaluator. Both mechanisms, volume derivatives,
  configuration, exchangers, warm states and numerical settings are matched.
  Six physical metrics use relative tolerance 2e-11 and absolute tolerance
  1e-11; statuses and constraint failures must agree. This is wrapper parity,
  not an estimate of physical-model accuracy. The two original-script checks
  require historical checkout artifacts; the packaged reference is portable.
- Persistence checks cover continuation, pending-candidate recovery, exact
  cache identity, integer decoding, runtime/source mismatch and real Sobol
  result serialization. Reports are checked without a callable physical solver.
- Chromium exercised the final offline HTML: candidate selectors changed the
  comparison values, status and constraint filters worked, both charts rendered,
  and small efficiency differences had distinct axis ticks. A 1400 × 1200
  screenshot was visually inspected.
- A wheel was built using the installed setuptools backend. Its console entry
  point and both bundled data files were verified. From an extracted wheel
  outside the checkout, `init` and `validate` succeeded without historical
  `examples/` or `outputs/` dependencies.

## Saved first-study evaluations

The editable [study](../outputs/research_first_study/study.toml), frozen basis,
two standalone results and [interactive comparison](../outputs/research_first_study/comparison.html)
are saved together. These use the historical **298.15/558.15 K** inlet-air
conditions, 2 Hz and fixed gas inventory. They do not replace the current
298.15/598.15 K demonstrator envelope.

| Quantity | Historical reference | H_i count changed to 4100 |
| --- | ---: | ---: |
| H_i tube count | 4029 | 4100 |
| Indicated gas power [W] | 41.4337050897658 | 41.53214627179981 |
| Indicated thermal efficiency [1] | 0.23255335337157246 | 0.23249272480082311 |
| External-air heat input [W] | 178.16859868523701 | 178.6384769991425 |
| Maximum pressure [Pa] | 419804.92714453617 | 418803.3687248253 |
| Maximum gas temperature [K] | 554.0162928282651 | 554.0799658328734 |
| Maximum absolute mass flow [kg/s] | 0.007420304051590775 | 0.007405498969694076 |
| Periodic cycles | 13 | 13 |
| Status | feasible | feasible |

The reference's six recorded metrics differ from the frozen historical result
by **zero** in this runtime. The changed-count evaluation illustrates a small
power increase accompanied by a small efficiency decrease, without optimization.
Useful mechanical power remains unavailable because mechanical losses are
unknown. External-air aerodynamic losses and fan consumption are excluded from
this trial balance, not assumed physically zero.

The separate [campaign report](../outputs/research_first_study/campaign.html)
contains a real one-candidate run followed by a one-candidate resume:

| Sequence | Phase | Outcome |
| --- | --- | --- |
| 0 | 1 | Feasible after 22 cycles; 27.20932066597371 W indicated, efficiency 0.1898589002032464 |
| 1 | 2 | `invalid_exchanger`: `MicrotubeDomainError: large_relative_pressure_drop`, retained after the permitted safe retry |

The second outcome is a recorded model-domain rejection, with no objective or
synthetic physical metrics. The journal contains sequence indices 0 and 1,
and no pending candidate remains. All saved evaluations and the campaign
passed identity verification and match the final implementation runtime.
This short Sobol search verifies operation and continuation; it makes no claim
about historical candidate order or search quality.

## Researcher review still required

Use the [guide](DADA_ENGINE_RESEARCH.md) to modify bounds, validate inputs,
evaluate an exact configuration and compare candidates. Before other campaign
families are added, record the user's experience of locating parameters,
understanding validation errors, distinguishing evaluation from optimization,
and interpreting efficiency, power, constraints and failures. Agent-operated
checks above do not substitute for that real first-study review.
