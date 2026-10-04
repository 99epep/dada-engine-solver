# DADA Research

Run this helper from the repository root; no package installation is required:

```sh
research() { PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m dada_solver.research "$@"; }
```

Research configures, evaluates, searches and compares machines using the existing
physical solver. Indicated gas power is not useful shaft power: mechanical losses
remain unknown. The local study and its portable basis define the science.

## Daily commands

```sh
research validate outputs/my_study/study.toml
research evaluate outputs/my_study/study.toml --output outputs/my_study/seed.json
research run outputs/my_study/study.toml --budget 30m
research resume outputs/my_study/study.toml --budget 6h
research status outputs/my_study/campaign
research report outputs/my_study/campaign
```

- `validate` checks declarations and builds geometry without integration.
- `evaluate` evaluates fixed/initial values once, without a search.
- `run` creates `campaign/` beside the TOML; `--directory PATH` overrides it.
- `resume` accepts a campaign directory or a TOML path (its sibling `campaign/`).
  It uses stored snapshots, not subsequent edits to the original TOML.
- `status` reads existing results without integrating.
- `report` writes **`CAMPAIGN/report.html`** by default. Regeneration overwrites
  only this derived report. `--html PATH` selects another destination.

A nonempty campaign is never overwritten by `run`. Missing resume directories
are explicit errors. `--max-candidates` limits new attempts for this invocation;
generated studies default to 512. Budget and attempt limits both apply. Use an
explicit small cap for checks. Ctrl-C preserves the pending candidate for resume;
`--retry-incomplete` deliberately retries a deadline-limited attempt.

Interactive output has one permanent line per evaluation and a status line
refreshed in place. Redirected output has no ANSI/control-character updates.
Use `--list-candidates` for a detailed terminal report, or `--json` for complete
structured inspection; replay notifications go to stderr, not into JSON stdout.

## Reports, curves and animations

By default, `report` adds **positions, heat-transfer rates, pressures and temperatures for the
two best feasible candidates**, in that order, using the campaign's existing objective ranking.
For incompatible studies, comparisons never invent a combined ranking.
The HTML is standalone/offline; opening it performs no local command or integration.

```sh
# Default plots, default output: campaign/report.html
research report outputs/my_study/campaign

# Choose candidates explicitly (full ID, unique prefix, best, or second)
research report outputs/my_study/campaign --candidate best --candidate second

# Explicit plot list replaces the default list; repeated --plots also works
research report outputs/my_study/campaign --plots default,temperatures,flows
research report outputs/my_study/campaign --plots positions,velocity,acceleration
research report outputs/my_study/campaign --plots default,mechanisms
research report outputs/my_study/campaign --plots all

# Fast inspection: no curve reconstruction or thermodynamic integration
research report outputs/my_study/campaign --plots none
research report outputs/my_study/campaign --plots none --json

# Standalone evaluation artifacts and multiple studies are also supported
research report outputs/my_study/seed.json --html outputs/my_study/seed.html
research compare outputs/a/campaign outputs/b/campaign --plots default --html outputs/comparison.html
```

| Plot | Meaning | Thermodynamic replay |
|---|---|---|
| `positions` | Normalized cylinder position, `(V-Vmin)/Vswept`, small and large | No |
| `volumes` | Actual small/large gas cylinder volumes | No |
| `pressures` | Small, large and both exchanger gas pressures | One cycle |
| `heat` | Signed gas → wall heat rates on each exchanger; gas → reservoir for reservoir models | One cycle |
| `temperatures` | Four gas temperatures and both walls when present | One cycle |
| `flows` | Signed mass flow through all four network passages | One cycle |
| `velocity` | Derivative of normalized position per solver-cycle radian | No |
| `acceleration` | Second derivative per solver-cycle radian², where the family supplies it | No |
| `mechanisms` | SVG joint/link animation for slider-crank, four-bar and six-bar | No |

Cylinder curves use the production operation-direction transform exactly once.
Positions and derivatives do not imply a physical stroke when no scale exists.
Signed flows retain the production port names and expose reflux. Heat is positive
from gas to wall (negative when the wall heats the gas). Only gas-side heat is
plotted; it is not the external-boundary cooling load when wall storage varies.
The heat panel includes a zero reference line, and ordinate ticks use round values
without rounding the scientific curve samples.

Thermodynamic plots are **new derived curves**, not a replacement for stored
metrics or validity. Research announces the replay and integrates exactly one
cycle from the saved final periodic state, with the study's numerical settings
and a bounded evaluation deadline. It does not rerun periodic convergence,
optimization or the campaign's warm-start/retry search. The report displays the
endpoint drift and actual backend. A missing state or failed replay is marked
unavailable; the historical candidate verdict remains unchanged.

A source/runtime mismatch is visible: reconstructed curves use the current
checkout. They must not be mistaken for an exact recovery of an older trajectory.
The compressed `.research-plot-cache/` beside the HTML holds disposable replay
data keyed by study, candidate, saved state and runtime. Matching cache entries
avoid repeat integration. One replay supplies every thermal plot. Campaign
journals and metrics are never modified by report generation.

SVG mechanisms have play/pause and an angle slider. They use precomputed
production joint positions and linkage connectivity, with equal geometric axis
scales. Each piston retains its local coordinate frame; side-by-side views are
not a reconstructed common-shaft assembly. The angle is explicitly the **study
angle before operation reversal**, distinct from the solver-angle curve panels.
Animation speed is illustrative, not physical machine speed. Dimensionless
geometry stays labelled as such. Abstract laws show no invented mechanism.

Detailed HTML records normally retain the best distinct feasible 10% of journal
attempts (rounded up). The two default plot targets and explicitly selected
candidates are retained even outside that cap. Funnel, rejection counts, bounds
and progress still cover all attempts. Journals and JSON inspection retain the
complete scientific records. Sorting, candidate comparison and constraint margins
remain available; recovered trial failures are separate from final-cycle limits.

## Define a study

```sh
research init kinematics --small hybrid_compact --large hybrid_compact --output outputs/new/study.toml
research init external-stream-refrigeration --small harmonic --large harmonic --output outputs/cooling/study.toml
```

A parameter with `value` is fixed. A parameter with `initial`, `lower`, `upper`
and its transform/encoding is active. Categorical valve placements use
`kind = "choice"` and explicit `choices`. The basis supplies the other machine
inputs; do not independently optimize derived exchanger UA, hold-up or losses.

Four kinds of limits have different owners:

- **Model domain:** genuine correlation/EOS applicability, such as each
  microtube model's temperature, Mach and pressure-drop domain.
- **Study requirement:** chosen pressure, temperature, power, flow or mechanism
  quality constraints. Generic presets do not impose inherited motor values.
- **Search bound:** the optimizer's exploration interval, not a physical limit.
- **Numerical setting:** integrator, convergence and sampling controls.

See [kinematics, mechanisms and ownership](DADA_ENGINE_RESEARCH_KINEMATICS.md),
[external streams](EXTERNAL_STREAM_THERMAL_MODEL.md),
[microtube models, circular collectors and triangular pitch](MICROTUBE_GAS_MODEL.md), and
[working-fluid models](WORKING_FLUID_MODELS.md).
The [reference](DADA_ENGINE_RESEARCH_REFERENCE.md) details schema, scientific
identity, policies, parameter types, search scheduling, persistence and compatibility.

## Refine locally or change capacity

```sh
research refine outputs/source/campaign --candidate abc123 --radius 0.20 --output outputs/local/study.toml
research refine outputs/source/campaign --candidate abc123 --candidate def456 --radius 0.20 --output outputs/multi/study.toml
research refine replay_A.json replay_B.json --radius 0.20 --output outputs/multi_artifacts/study.toml
research rescale outputs/source/campaign --candidate abc123 --factor 5 --output outputs/larger/study.toml
```

`refine` creates a portable study; it does not run it. Radius 0.20 means ±20% of
the **global normalized interval**, clipped at its boundaries, not ±20% of each
physical value. Exact centers run first; independent Sobol regions then alternate
in deterministic round-robin order. Integer/category encodings and global bounds
remain unchanged. Small local regions can retain their center's categorical
choice; choices have no physical distance. See the reference for bin semantics.
Each `search.regions.center` must contain exactly the active coordinates. When
fixing parameters manually, remove them from every center table.

Sources must share a scientific identity. Resume preserves centers, pending work,
region order and Sobol indices. Region provenance does not enter physical
candidate identity. Global studies evaluate their initial only when explicitly
configured with `search.evaluate_initial = true`; historical studies keep their
stored schedule.

`rescale` changes extensive machine inputs and their bounds; it is not local
refinement or strict geometric similarity. Always evaluate the result. Rescale
the global source before refining it. See [capacity rules](DADA_ENGINE_RESEARCH_CAPACITY.md).

## Persistence and compatibility

New campaigns use lossless append-only `history.jsonl.gz`, a temporary
`recovery.json` and `state.json`. Old uncompressed journals remain readable and
continue in their original format; no automatic history migration occurs.
Exact caches, candidate IDs and resume checks remain scientific contracts.
Changes to physics, study definitions or runtime can prevent execution resume;
offline inspection remains available. Never combine results under a silently
changed definition.

## Further reading and future mechanism synthesis

- [Research technical reference](DADA_ENGINE_RESEARCH_REFERENCE.md)
- [Validation and measured limitations](validation.md)
- [Cockpit reference](DADA_ENGINE_RESEARCH_COCKPIT.md)
- [Limit ownership audit](RESEARCH_LIMIT_OWNERSHIP_AUDIT.md)
- [Historical migration matrix](DADA_ENGINE_RESEARCH_MIGRATION_MATRIX.md)

Fitting a physical mechanism to an abstract/hybrid law is a separate future
Research workflow; animation does not implement synthesis. The existing design
knowledge remains the starting point: primary-family discovery, downstream dyad,
local mechanism polish, opposite-side adaptation, then paired thermodynamic
optimization. Preserve several useful families and judge final machines through
thermodynamic replay; do not make acceleration matching the universal objective.
See [synthesis method](MECHANISM_SYNTHESIS_SEARCH.md),
[primary families](PRIMARY_FOUR_BAR_FAMILIES.md) and
[six-bar families](SIX_BAR_MECHANISM_FAMILIES.md).
