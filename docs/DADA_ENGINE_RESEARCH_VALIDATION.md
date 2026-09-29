# Dada-Engine Research validation

V2 validation is recorded below under [V2 kinematics and mechanism validation](#v2-kinematics-and-mechanism-validation). The first sections preserve the V1 review record.

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
passed identity verification and matched the implementation runtime when generated.
This short Sobol search verifies operation and continuation; it makes no claim
about historical candidate order or search quality.

## Initial researcher feedback and further trials

The user validated the first hands-on trials. The interface is not yet intuitive,
but the documentation makes it usable. The user then authorized the V2 kinematics/mechanism extension described below.
Further hands-on trials remain the usability gate before large research campaigns.

Use the [guide](DADA_ENGINE_RESEARCH.md) to modify bounds, validate inputs,
evaluate an exact configuration and compare candidates. Before other campaign
families are added, record the user's experience of locating parameters,
understanding validation errors, distinguishing evaluation from optimization,
and interpreting efficiency, power, constraints and failures. Include the
specific steps that still require consulting the guide. Agent-operated
checks above do not substitute for that real first-study review.

## Diagnostic integration maintenance

The deprecated NumPy `trapz` call in microtube diagnostics was replaced with
`scipy.integrate.trapezoid`. SciPy is already required, and this keeps support
for the declared NumPy 1.24 minimum. The physical-time trapezoidal rule and
weighting are unchanged. Direct comparisons on four nonuniform-grid inputs
(smooth, indicator, zero and signed values) produced exactly equal results.
The full suite then passed with deprecations promoted to errors: **595 passed,
no warnings**, 103.18 s, using
`PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 -m pytest -q -W error::DeprecationWarning`.

As with any source edit, this changes the campaign runtime fingerprint.
Previously saved results remain readable; execution resume requires the
original matching source version. No saved campaign or identity was rewritten.

## V2 kinematics and mechanism validation

The user authorized this extension after validating the first hands-on V1 tests.
The [V2 reference](DADA_ENGINE_RESEARCH_KINEMATICS.md) describes the implemented
schema, fixed/active ownership, families, artifacts, screening and deferred work.
No large optimization, primary-family saturation or mechanism-fitting campaign
was run. Thermodynamic physics and the existing campaign engine are retained.

### Historical parity evidence

`docs/research_audit/capture_kinematics_v2.py` captured dense historical
trajectories **before extraction**, with source file digests in
`tests/fixtures/research_v2/reference.json` and arrays in `trajectories.npz`.
There are 4321 angles from -2*pi to 4*pi for 11 motion configurations: harmonic,
slider, finite-rod four-bar, four six-bar pairs (ranks 1/4/12/50), spline, Fourier,
structured 3952 and ideal piecewise. Volume and first-derivative comparisons use
`rtol=2e-11, atol=1e-13`; available spline/structured second derivatives use
`rtol=2e-11, atol=1e-11`. Harmonic/slider/Fourier analytic acceleration also has
finite-difference checks. Phase, orientation and extrema remain those of the
original implementations. These tolerances test mathematical implementation
parity, not experimental physical accuracy.

The 15p helper and scalar path, Fourier class, slider classes and spline phase
wrapper were moved into production. Historical consumers import those classes;
search loops were not migrated. Production Research has no import of examples.
Four-bar and six-bar closure is reused directly. Six-bar joint export was
factored from the existing closure calculation, not reconstructed in a renderer.

`capture_kinematics_v2_thermo.py` separately froze machine inputs, original warm
starts and already recorded metrics into a compact portable thermal fixture,
without running any new integration. Its runtime regression tests do not import
historical examples or read large campaign histories.

| Historical thermal case | Result / comparison |
| --- | --- |
| Inverted-offset slider-crank, Stage K2 | Same 24 convergence cycles; maximum relative difference across six recorded metrics 2.45e-6. Power remains approximately 22.4627 W, below the explicit 25 W constraint: correctly `converged_infeasible` |
| Compact finite-rod coupler four-bar, Stage K2 | Same 23 cycles; maximum relative metric difference 5.01e-6; approximately 43.11695 W and efficiency 0.18497750 |
| Fourier 8H retained candidate | Original machine, resolved inventory and predecessor warm state; six recorded metrics pass 5e-8 relative / 1e-10 absolute; same 10 cycles |
| Retained source16 spline | Original fixed machine and predecessor warm state; same metric tolerance; same 3 cycles |
| V1 rank-01 six-bar through V2 declarations | Six metrics pass 2e-11 relative / 1e-11 absolute; same 13 cycles |
| Structured C2 candidate 3952 | Same stored values on this backend, including 40.719090181447584 W, efficiency 0.23508063005202154 and 3 cycles; regression tolerance 2e-11 relative / 1e-11 absolute |

The K2 integrated regressions use a measured 1e-5 relative tolerance. Roundoff in
volume repartition and common-frame normalization perturbs adaptive LSODA steps
and sampled extrema. The existing 2e-11 RHS equivalence tolerance is not a bound
on accumulated integration or extrema differences. No production numerical
setting or physical acceptance constraint was relaxed. Historical Python and
Numba backends are selected explicitly per fixture; Numba-specific stored tests
skip only when that optional backend is unavailable.

Known six-bar families round-trip through serialized artifacts, reproduce the
dense trajectories and preserve mechanical diagnostics. Historical stroke was
sampled whereas production refines its extrema; diagnostic tolerance is 2e-5
relative/absolute. The four-bar compact screen is checked on its original 2880
angle grid and original envelope frame; the numerical results are stored in
[four-bar screen parity](../outputs/research_kinematics_v2/four_bar_screen_parity.json).
These are sampled screens, not continuous mechanical guarantees.

### Bounded acceptance artifacts

Run the complete demonstration in a new directory with:

```sh
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 examples/research_kinematics_demo.py --output outputs/research_kinematics_v2_new
```

The recorded working artifacts under
[`outputs/research_kinematics_v2/`](../outputs/research_kinematics_v2/summary.json)
include editable studies, portable bases/mechanisms, exact evaluation records,
per-study HTML and a [cross-family comparison](../outputs/research_kinematics_v2/comparison.html).
They use the historical 298.15/558.15 K air inputs, not the current 598.15 K hot
inlet demonstrator target. The first seven evaluations share a thermal basis;
3952 has its own original thermal hardware. Reference motions were not optimized
under common bounds, so these results do not establish a fair family ranking.

| Family | Indicated power [W] | Indicated efficiency | Outcome |
| --- | ---: | ---: | --- |
| Harmonic | 40.262075 | 0.19750494 | Feasible |
| Slider-crank | 28.628784 | 0.17604739 | Feasible |
| Four-bar | 39.017753 | 0.20993799 | Feasible |
| Six-bar | 41.433705 | 0.23255335 | Feasible |
| Free spline | 40.945586 | 0.23536605 | Feasible |
| Fourier C2 | 38.942669 | 0.23243461 | Feasible |
| Ideal piecewise | 0.776594 | 0.00244381 | Converged, below minimum indicated power |
| Structured C2 3952 | 40.719090 | 0.23508063 | Feasible |

The [active-phase campaign](../outputs/research_kinematics_v2/active_phase/report.html)
has one active harmonic phase, two distinct Sobol points and a real disk-based
resume between them. Sequence indices are 0/1, phase attempts are 1/1, both
converge (22/13 cycles), and no pending state remains. Indicated powers are
40.03216852 W and 40.39353090 W. This verifies persistence, not search-order
parity with a historical optimizer. Warm starts remain initial guesses followed
by full periodic validation; exact caching has not been approximated.

Motion exports for the six-bar and 3952 studies contain normalized position,
angle derivatives, volume and available geometry, without thermodynamic replay.
The separate [retained mechanism library](../outputs/research_kinematics_v2/mechanism_families.json)
contains all four historical pairs, preserving diversity rather than overwriting
it with one winner. Its entries were exported from the frozen parity inputs;
no new mechanism search was performed.

Browser verification used headless Chromium against the actual offline report.
It checked candidate A/B selection, the status filter, cross-study warning,
constraint values/limits/margins/states and near-active amber rows. The rendered
constraint table was inspected visually. All stored evaluations and the resumed
campaign pass exact identity and current runtime checks. Older V1 reports remain
inspectable even when their strict execution runtime fingerprint differs.

### Complete test result

```sh
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 -m pytest -q -W error::DeprecationWarning
```

**679 passed, no warnings, 258.02 s.** This includes the unchanged V1 regression
coverage and new schema/ownership, mixed-family, fixed/active, integer/branch,
artifact, exact evaluation, Sobol/resume, offline report/margin, historical
kinematic/thermal parity and synthesis-contract tests. The hierarchical data
model covers 6/9/15/30 active coordinates, independent sides, fixed branches and
multiple retained families. No existing test was weakened or removed.

The cooling objective and indicated-input constraint compile through the same
engine with the reservoir model. At the V2 milestone, the air-wall microtube wrapper still
required motor operation. V3 removes that restriction through the separately
validated external-stream model described below. The human-powered cell was not optimized. Automatic mechanism fitting,
local polish, primary-family discovery and saturation remain deferred contracts;
see the [architecture boundary](DADA_ENGINE_RESEARCH_KINEMATICS.md#established-synthesis-workflow-and-current-boundary).
At the V2 milestone, the two older compact hybrid laws remained example-only.
The nine-coordinate `HybridCompactKinematics` has since been extracted as
`hybrid_compact`; the distinct 11p fitting law remains example-only. Further human review should focus on
configuration usability, artifact-derived defaults, candidate comparison and
study-specific mechanical thresholds before a large research run.


## Research V3 — external boundaries and compiled property tables

Baseline: 679 passing tests, no warnings. V1/V2 tests remain in place. The new
suites are `test_external_stream_v3.py`, `test_tabulated_fluid_v3.py` and
`test_research_v3.py`. They cover neutral/legacy rate and outlet parity, a frozen
pre-V3 one-cycle trajectory, signed refrigeration, fixed/active stream ownership,
table identity and immutable storage, domain/cell-mask rejection, prohibited
ideal hydraulics, state-enthalpy transport, compiled lookup without Python
properties, periodic table replay, candidate construction, resume and reports.

`tests/data/pre_v3_air_wall.json` was generated from the local pre-V3 committed
source exported to an isolated temporary directory. It captures the rank01
source warm state, derivatives, 65-point interpolated trajectory, endpoint and
heat/work quantities. The original stored-reference Research test independently
retains its 2e-11 relative / 1e-11 absolute metric and periodic-count checks.
No physical constraint or transport-domain threshold was relaxed.

The bounded [acceptance summary](../outputs/research_v3/validated/summary.json)
records all cases as feasible. The refrigerator converges in 11 cycles at about
5.36590 W cooling, 7.21388 W heating, 1.84798 W indicated input, cooling COP
2.90366 and heating COP 3.90365. Its stored-energy-aware absolute cycle energy
residual is approximately 6.1e-13 J. The neutral-field air motor converges in 13
cycles. These are numerical fixtures, not physical hardware validation.

The compiled ideal-generated table refrigerator also converges in 11 cycles.
Relative differences versus the exact ideal path are approximately 3.50e-8 for
cooling power, 1.11e-7 for indicated input and 1.46e-7 for COP. The regression
uses 2e-6 relative / 1e-7 absolute for integrated metrics, allowing adaptive
trajectory sensitivity; instantaneous RHS checks remain at 2e-10 or tighter.
The table benchmark's same-state maximum absolute RHS error is about 9.12e-13.
These are numerical parity tolerances, not real-fluid accuracy estimates.

Two tiny Sobol candidates (one before and one after true persisted resume) are
feasible. Reports preserve stream definitions and constraint margins; the
structured-C2 refrigerator input is validated for configuration and mechanical
preflight only, without claiming a successful cell design.

The initial full-suite run exposed a diagnostic test double that has no volume
provider. Ideal temperature diagnostics now keep their volume-independent fast
path; general-fluid diagnostics require the actual volumes. The existing test
was preserved. See the final V3 test result in the [V3 artifact README](../outputs/research_v3/README.md).


## Compact hybrid extraction

`hybrid_compact` is the eleventh Research family. Dense scalar and vector
position/first-derivative references were frozen from the original example
before replacing its implementation with production imports. Tests cover two
parameter sets over 4097 angles spanning three cycles, phase/extrema/direction,
original domain guards, all nine active coordinates, mixed families and V3
schema construction. No analytic second derivative is invented.

The packaged candidate 501 machine and original preceding warm state reproduce
all recorded thermal metrics at 2e-11 relative / 1e-11 absolute tolerance:
40.555304211350325 W indicated output, efficiency 0.22525540307940708,
180.04142700654225 W external input, and 7 cycles. Dense trajectories are
bitwise equal to the frozen original implementation on the validation runtime.
The bounded [Research replay](../outputs/research_hybrid_compact/champion.json)
contains the full evaluation; this is no new optimization.

Post-extraction validation: **744 passed, 0 warnings in 293.79 s** using
`PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 -m pytest -q`.

## Capacity scaling, chronology availability and offline plots

The [capacity-scaling demonstration](../outputs/human_cell_stage0/README.md)
records local Human Cell candidates A and B, factor-one reconstruction and
bounded factor-five evaluations. With identical uniform initialization, the
original A and A×1 complete metric dictionaries are exactly equal. A×5 / B×5
produce 231.885 / 268.295 W cooling with 130.274 / 156.899 W indicated input
and COP 1.77998 / 1.70998. Both converge in 15 cycles and are model-valid and
study-feasible without constraint relaxation.

Tests cover basis-owned CdA, count-ratio policy without double scaling, inventory,
extra dead volume/capacity, external streams, active bounds, integer rejection,
portable mechanism artifacts, reversal, non-mutation and exact factor-one input
parity. Reporting tests cover offline production volume reconstruction, single
motor-direction application, standalone evaluation reconstruction and legacy
empty-event presentation without rewriting stored diagnostics. V3 evaluation
checks explicitly distinguish unavailable wall events from reflux.

Twelve browser checks additionally exercise real HTML sorting and the two
historical candidates; table and four volume curves were visually inspected.
See [the exact scaling and diagnostic rules](DADA_ENGINE_RESEARCH_CAPACITY.md).

Final full-suite result: **760 passed, 0 warnings in 258.43 s** using
`PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 -m pytest -q`.
The old static-cycle assertion was updated to require unavailable chronology
for an empty event sequence, as requested; pressure checks and nonempty-event
classification tests remain intact.

## Research cockpit and journal-first UX

The targeted UX regression set passes 93 tests covering default/overwritable
HTML, compact CLI output, explicit listing/full JSON, unique/ambiguous IDs,
real fd 2 writes and real LSODA callback failure with restoration, TTY/non-TTY
progress, 512 execution defaults, unchanged scientific/candidate identity when
only execution caps change, deterministic Sobol continuation and exact cache.

Crash cases include recovery-only completion, journal plus residual recovery,
torn append, legacy candidate files, state-publication failure, invalid recovery
identity and a live writer acknowledging recovery during read-only inspection.
No new per-candidate result file is produced. Original empty-event and physical
parity tests remain in the suite.

Nine browser checks and a visual review cover the live selection/copy commands
and retained plots in the [standalone cockpit demonstration](../outputs/research_ux/cockpit.html).
No long optimization was used for validation.

The single final full-suite run passed **782 tests, 0 warnings in 286.82 s**:
`PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 -m pytest -q`.
