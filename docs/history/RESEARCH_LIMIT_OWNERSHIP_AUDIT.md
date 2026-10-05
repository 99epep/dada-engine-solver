# Research limit ownership audit — 2026-09-30

> Historical cleanup audit dated 30 September 2026, retained as migration provenance.
> The current ownership rule is defined in [Physics decisions](../PHYSICS_DECISIONS.md).

This audit uses the local checkout, including saved outputs. It does not infer
physical ratings from historical optimizer guards. No optimization was run.

## Ownership rule

| Kind | Owner | Examples |
|---|---|---|
| Model domain | Selected mathematical/empirical model | Transport temperature domain, gas-model Mach and pressure-drop domains, closure/real-valued linkage geometry |
| Study/design constraint | Explicit study | Minimum power, maximum pressure/temperature/flow, transmission quality, maximum enclosed volume |
| Search bound | Search definition | Normalized/physical parameter bounds and local-region radius; being outside a search interval is not physical invalidity |
| Numerical setting | Evaluation protocol | Integration tolerances, convergence criteria, screening samples; these are not machine ratings |

Generic V2/V3 initializers no longer add the inherited motor guards of 25 W,
1.2 MPa, 850 K or 0.08 kg/s. All four constraint types remain supported when
explicitly declared. The external-stream refrigeration validation preset still
explicitly declares a 100 W indicated-input requirement for that named fixture;
it is not a solver or general refrigeration limit.

`screening.maximum_large_enclosed_volume_m3` is optional and may exceed 0.066.
The 66 L value remains a historical demonstrator design requirement, never a
universal motor restriction. Generic presets omit it. Numerical positivity,
geometry closure and actual exchanger domains remain enforced.

## Absolute-flow inventory and decisions

The machine-readable [inventory](../research_audit/absolute_mass_flow_inventory.json)
records file-level matches from the full local audit. Searches included
`maximum_absolute_mass_flow`, `MAXIMUM_ABSOLUTE_MASS_FLOW`, `MAX_FLOW`, decimal
and exponent spellings of 0.08, and possible 80 g/s spellings. Numeric matches
alone include geometry, plot settings, search radii, external stream inputs and
recorded measurements; they are not evidence of an internal-flow limit.

| Important location | Classification and action |
|---|---|
| `src/dada_solver/research/presets.py` | Generic inherited constraint removed; V2/V3, mixed families and champion geometry presets inherit the cleanup |
| `src/dada_solver/research/data/sixbar-thermo5d.toml` | Explicitly named V1 historical physical-parity fixture; retains 25 W, 1.2 MPa, 850 K and 0.08 kg/s as historical study requirements |
| `sizing/constraints.py`, `sizing/configuration.py`, Research schema/margins | Optional explicit mass-flow constraint retained; no implicit cap |
| Campaign evaluator, progress, report, diagnostics | Absolute peak flow observable retained, including `maximum_absolute_mass_flow_kg_s` |
| Tests | Explicit-cap fixtures retained; added no-default and >0.08 acceptance tests with real domain guards |
| Saved studies, JSON/JSONL results, reports and campaign definitions | Historical evidence retained, not relabelled or reclassified |

The following historical executable searches define 0.08 directly, without a
component-based justification found in their source or documentation:

- `examples/optimize_motor_two_extrema_stageF6_11d.py`
- `examples/optimize_motor_exchanger_asymmetry_stageK1.py`
- `examples/optimize_motor_piecewise_stageP1.py`
- `examples/optimize_motor_piecewise_stageP3.py`
- `examples/refine_motor_four_stage_hx9d_variable_gas.py`
- `examples/optimize_motor_temperature_point.py`
- `examples/optimize_motor_temperature_point_fixed.py`

These remain historical protocols, not recommended defaults. Imports propagate
some guards into Fourier/C2, free-spline/mobile-spline, hybrid, sigmoid/beta,
valve-screening and four-stage K2/11D refinements. Importing only a geometry helper
does not necessarily apply its feasibility guard. Those scripts were not silently
changed; use a newly declared Research study for new work.

Historical TOML occurrences include `motor_four_bar_thermo_stage4/4b`,
`motor_mechanics_stage5_{A_reference,B_low_ratio,C_high_ratio}`,
`motor_thermo_stage6_{A,B,C}`, `motor_mechanics_stage7A{,2,3,4,5}_edge`,
`examples/human_cell_stage0/study.toml`,
`docs/research_audit/proposed_rank01.toml`, and saved Research campaign/output
TOMLs listed individually in the inventory. They retain explicit inherited
requirements, not validated component ratings. Other explicit 0.05 and 0.02
example/controlled-test limits are not replacements for the removed default.

Historical discussion is retained in `FOUR_STAGE_K2_SEARCH.md`,
`FOUR_STAGE_OPTIMIZATION.md`, `VALVE_PLACEMENT_SCREENING.md`,
`MOTOR_TEMPERATURE_AND_VALVE_CAMPAIGN.md`, the Research architecture audit,
implementation plan and validation results. Their feasibility counts describe
the constraints used then.

No exchanger, transport, hydraulic or compiled primitive was found to assume
0.08 kg/s was safe. Those paths depend on geometry, local state and their own
model domains. None of their correlations, Mach/pressure-drop/Knudsen/transport
bounds or integration tolerances was changed by this cleanup.

## Validity cleanup and compatibility

Pressure inequality is resolved by the four-volume hydraulic/thermodynamic
balances. `maximum_pressure_equalization_error` remains a measured report
quantity, but no longer contributes to thermodynamic validity. Its old
configuration threshold is ignored. Old explicit sizing pressure-equalization
constraints and corresponding optimizer scales are also ignored.

Isothermality excursion is removed: no hot/cold isothermality fields are computed
or serialized, and no isothermality sizing constraint is constructed. Old TOML
validity keys and old sizing constraint/scale entries remain loadable and are
ignored. The legacy `ValidityThresholds` constructor accepts the old positional
or keyword arguments solely for compatibility; its dataclass state now contains
only Cp-variation and compressibility-deviation tolerances. Newly generated V2/V3
basis files omit obsolete keys. Historical basis bytes and stored records are
not rewritten; their raw diagnostic JSON may still contain obsolete fields.

Cp/property variation and compressibility-factor deviation remain active model
approximation checks. Microtube Mach validity uses each selected film and
hydraulic passage's own gas model (typically `maximum_mach = 0.3`), including
asymmetric per-side declarations. The previous generic 0.2 veto is removed.
An explicit Research `maximum_mach_number` constraint can independently impose
0.2 or another design requirement on the observed maximum. Such a rejection is
a study constraint, not a correlation-domain failure. Legacy constant-property
microtubes retain their laminar Reynolds screen and sonic hydraulic closure;
no undocumented universal Mach threshold is substituted for the removed veto.

Stored scientific definitions/IDs remain inspectable without migration. A fresh
evaluation with this code can change validity classification, though its solved
trajectory is unchanged. The existing source/runtime fingerprint changes and
prevents unsafe resume/cache reuse across the implementation change. For an
exact historical verdict use its stored result or historical implementation.
Do not resume an old campaign under a newly edited definition. Regenerate the
Human Cell study and start a new campaign directory after reviewing its explicit
constraints. `refine` and `rescale` deliberately preserve source constraints;
they do not silently clean historical caps from source studies.

## Mechanism limits

New four-/six-bar seed artifacts contain geometry, branches, conventions and
provenance; preset quality limits are emitted as editable study-level
`[[mechanical_constraints]]`. The six-bar defaults retain the former stroke,
transmission, rod-axis, reach, clearance, lateral-motion and reversal values.
Changing/removing those rows changes study identity without rebuilding the seed.
Intrinsic closure and real-valued construction remain in production mechanisms.
Old artifacts with embedded constraints still load and enforce them exactly;
duplicate study/artifact constraints remain rejected.

## Report evidence and actions

The identical-command bug was in `campaign_evidence`: each diagnostic suggestion
was assigned the report-regeneration command, which the renderer displayed as
supplied. Failure suggestions now offer an offline filter for their stable
failure category; constraint suggestions filter the exact violated constraint.
Both open stored candidate details. Bound-pressure suggestions scroll to the
bound evidence. Only the regeneration suggestion carries a regeneration command;
existing exploration and selected-volume commands retain their real purposes.
No unsupported CLI inspection flag was invented. Selection/filtering does not
run a solver or local command.

The table retains all rows in a roughly 15-row scroll viewport, with sticky
headers and an explicit Basin column. ID abbreviation now targets the ID cell
instead of accidentally overwriting the origin cell. Comparison defaults use
the existing feasible objective ranking; incompatible studies get no combined
ranking. Manual selection is unchanged.

Constraint evidence is ordered by violated/unavailable/near-boundary/satisfied
state. It shows actual value, declared boundary, signed/relative margin, a 5%
near-boundary display cue and candidate/passage context. It includes stored
model-domain and first-rejected-state evidence when available; missing limits
remain unavailable rather than guessed. Sampled cycle maxima are labelled as
sampled evidence, not continuous guarantees. Historical verdicts remain stored
verdicts, including old cap or old thermodynamic-validity failures.

The V1-only `research/sixbar.py` adapter still enforces its historical 66 L
requirement; it is not used as the V2/V3 generic mechanism adapter. The remaining
0.0668 in the Hausen equation is a correlation coefficient, not a volume limit.
No other global flow rating was discovered. The provenance of the historical
25 W/1.2 MPa/850 K/0.08 guards as technical component ratings remains unsupported.

Three historical Python callers needed compatibility edits: the motion-law
comparison now uses current model-domain validity; hardware screening reports
Mach-domain status only when the gas model supplies it; the saved Stage 5
contract checker reads its former guard explicitly from its historical source
and drops obsolete isothermality fields when constructing a current report.
It compares saved benchmark assessments, not the validity of a newly designed
machine. These changes do not rewrite benchmark files.

Kinematic diagnostics are recorded even when a study declares no mechanical or
volume constraints. Removing a design requirement must not remove the associated
observable. When fewer than two feasible candidates exist, the comparison's
remaining slot uses stored order without assigning a new infeasible ranking.

The packaged bases retain explicit 0.01 Cp-variation and compressibility-factor
approximation tolerances. Their provenance as universal real-fluid acceptance
limits has not been established; future fluid studies must review these declared
values against the chosen approximation. They are not component safety ratings.
For the current calorically perfect gas both measured deviations are zero.
