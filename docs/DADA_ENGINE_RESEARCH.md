# Dada-Engine Research

The first working release lets a researcher configure, evaluate, search and
compare a fixed independent six-bar pair with five thermal/hardware coordinates.
This is the first validation study, not a permanent restriction on the future
research module. The existing campaign engine owns persistence, exact identities,
Sobol continuation, deadlines and recovery. Production physics is unchanged.

## Start from the checkout

No installation is required. In the commands below, `research` is a shell helper
for the module entry point:

```sh
research() { PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m dada_solver.research "$@"; }
research init sixbar-thermo5d --output outputs/my_research/study.toml
research validate outputs/my_research/study.toml
```

An installed package also exposes the equivalent `dada-research` command.
The wheel includes both the template and reference basis; it does not need the
checkout's historical `outputs/` or `examples/` directories to run this study.
NumPy and SciPy are required. The existing optional Numba backend accelerates
evaluation when installed; its actual use or Python fallback is recorded.
Rendering HTML requires no plotting package, server or internet connection.

`init` creates an editable TOML file and a frozen `study.basis.json` beside it.
It refuses to overwrite either file. The basis contains both mechanisms,
resolved hardware and configuration, initial state, wall capacities, numerical
settings and historical source digests. The source paths are provenance, not
live dependencies. The digest in TOML detects changes to the copied basis.

`validate` checks schema, units, ownership, source integrity and cheap geometry
construction without a periodic solve. Use `validate ... --json` for the full
validation/runtime metadata.

## Configure the first experiment

Edit the five `[[parameters]]` tables in `study.toml`:

| TOML parameter | Meaning | Unit | Short name for `evaluate --set` |
| --- | --- | --- | --- |
| `volume.swept_ratio` | Small / large swept-volume ratio, at fixed total swept volume | 1 | `swept_ratio` |
| `microtube.heat_in.tube_count` | H_i parallel tubes | integer count | `n_i` |
| `microtube.heat_in.tube_length_m` | H_i tube length | m | `length_i_m` |
| `microtube.heat_out.tube_count` | H_o parallel tubes | integer count | `n_o` |
| `microtube.heat_out.tube_length_m` | H_o tube length | m | `length_o_m` |

`lower` and `upper` bound the search; `initial` defines the default standalone
evaluation and the unevaluated initial reference in campaign reports. Sobol
does not insert that initial point. Counts must be integers, not `4029.0`.
The integer decoder uses nearest-even rounding, with inclusive bounds and
half-width endpoint bins; it is not a uniform categorical sampler.

The pair, gas inventory, source temperatures, frequency, total swept volume,
clearance ratios, transport model and finite external-air flow are fixed in
this protocol. Those assumptions are visible under `[fixed]` and `[policies]`.
Changing them requires a new explicit basis/study rather than editing a label
while retaining incompatible physics. The adapter changes exchanger hold-up,
wall capacity, heat transfer and losses through geometry and scales valve CdA
with tube count. Those derived quantities are not independent parameters.

Keep the five physical constraints explicit. Their thresholds are configurable;
units and meanings are validated. After changing bounds, constraints or physical
inputs, start a **new campaign directory**. `resume` reads its stored definition,
not edits to the original TOML. Presentation names, input path relocation and
per-invocation time/count budgets do not change the scientific study identity.

## Evaluate a configuration without optimization

```sh
research evaluate outputs/my_research/study.toml --reference --output outputs/my_research/reference.json
research evaluate outputs/my_research/study.toml --reference --set n_i=4100 --output outputs/my_research/trial.json
research compare outputs/my_research/reference.json outputs/my_research/trial.json --html outputs/my_research/comparison.html
```

Open `comparison.html` in a browser. `--reference` uses the stored historical
champion's five coordinates. Without it, evaluation starts from TOML initials.
Repeated `--set` accepts either short names or full parameter names. Requested
physical floats are preserved exactly; normalization is used for identity and
distance, not to perturb the requested configuration through a decode roundtrip.
An existing evaluation output is never overwritten.

Standalone evaluations are immutable JSON artifacts, not Sobol journal entries.
They retain the complete scientific definition, candidate identity, source
lineage, physical values, statuses, constraints and numerical diagnostics.
They can be inspected with `status`, `report` and `compare`. Their default time
limit is 180 s; override with `--budget` if needed. A stopped standalone
evaluation must be requested again. A deadline result is never labelled feasible.

**Physical reproduction means the same configuration gives the same physical
results within the declared numerical tolerance.** No agreement with the old
optimization order, incumbent updates, shrinking radii or final best candidate
is claimed or required. Historical candidate IDs remain lineage references;
new campaign IDs also cover normalized coordinates, definition and runtime.

## Run, stop and resume a small search

```sh
research run outputs/my_research/study.toml --directory outputs/my_research/campaign --budget 2m --max-candidates 1
research resume outputs/my_research/campaign --budget 2m --max-candidates 1
research status outputs/my_research/campaign
research report outputs/my_research/campaign --html outputs/my_research/campaign.html
```

`--max-candidates` is a cap for this invocation, not a lifetime total. Both the
wall-clock budget and that cap apply. The initial evaluation estimate defaults
to 30 s, so a smaller budget can legitimately start no candidate. Subsequent
estimates use the existing robust integration-time history.

Ctrl-C leaves the campaign's pending candidate recoverable. Budget interruption
retains the last complete cycle as an initial guess where available; partial
cycles cannot establish convergence. Use a larger budget to deliberately retry
the latest deadline-interrupted evaluation:

```sh
research resume outputs/my_research/campaign --budget 5m --retry-incomplete
```

Stopping is cooperative at solver progress boundaries. A per-candidate deadline
and a campaign deadline apply together. At most one explicitly recorded safe
initial-state retry follows a microtube-domain error. Its gas state is rebuilt
at fixed inventory and its wall temperatures retained; it does not change
states during a cycle or weaken periodic convergence.

The existing exact cache remains unchanged: normalized coordinates are part of
candidate identity. Different Sobol coordinates that decode to the same tube
count are not silently merged. Unexpected exceptions leave pending work for
recovery rather than producing a feasible synthetic result.

## Compare and inspect

The offline HTML contains status and constraint filters, efficiency versus
indicated-power and progress plots, and two candidate selectors. The comparison
table gives values and B − A differences for parameters, metrics and constraint
margins. It also shows achieved external/internal peak capacity-rate ratios and
minimum signed port flows. Negative local flow remains visible. The recorded
reflux flag uses a declared 1e-8 kg/s threshold; raw extrema are retained.

Candidate detail exposes cycle count, normalized periodic error history, warm
source, retry, backend and phase timings. The physical boundary is always
external-air heat input. Missing metrics display as unavailable. A failed
solver or unavailable constraint cannot enter the feasible ranking. Saved
configuration and source digests remain inspectable in the provenance section.

To select two candidates explicitly, copy IDs from `status` and use full IDs
or unambiguous prefixes:

```sh
research report outputs/my_research/campaign --candidate best --candidate CANDIDATE_ID --html outputs/my_research/selected.html
```

`report` and `compare` never run physics. They read completed orphan records and
ignore a torn final append for inspection without repairing or writing to the
journal; execution resume performs the existing recovery. Runtime/source changes
are labelled while inspection stays available. Execution resume requires the
recorded compatible runtime and code. If scientific study identities differ,
cross-study values are labelled and no combined ranking is assigned.

Trajectories, synchronized animation, a graphical editor, local refinement and
other study protocols are deferred. No renderer silently replays missing data.

## Scope of this reference machine

The validation preset is historical: **298.15/558.15 K**, 2 Hz, fixed gas mass
0.0012027347024506055 kg, and a 25 W indicated-power floor. It is separate from
the current 298.15/598.15 K demonstrator research envelope and approximately
100 W useful-output goal. The raw base configuration still contains a 598.15 K
reservoir seed; the actual injected wall source uses 558.15 K inlet air. Reports
use the actual exchanger inputs. No metadata normalization changes the model.

Mechanical losses remain unknown, useful output unavailable, and external-air
aerodynamic losses/fan power excluded from the balance rather than physically
zero. Finite air flow, thermal-film resistance, pause heat transfer, wall
storage and signed local reflux remain active. The large-cylinder maximum
enclosed volume, including clearance, is checked against the 0.066 m³ ceiling.
Mechanism lengths are in crank-radius units; no physical stroke or manufacturing
scale is inferred. Isothermality remains diagnostic-only.

## Acceptance and the next review

See the [validation record](DADA_ENGINE_RESEARCH_VALIDATION.md) for measured
results and the ready-to-open first-study artifacts in
[`outputs/research_first_study/`](../outputs/research_first_study/).

The reference and two changed configurations are compared with the original
evaluator under matching inputs, warm states and numerical settings. The
wrapper tolerance uses the existing backend-equivalence scale (relative 2e-11,
absolute 1e-11), not a physical-model accuracy claim. Pure persistence tests
separately verify split-run order, exact cache and pending-candidate recovery.
Real Sobol serialization is covered in addition to mocked runner tests.

Run the checks from the source checkout:

```sh
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 -m pytest -q
```

Before adding other campaigns, the user must try this first study and assess:

1. Can the five variables and units be located and changed without Python?
2. Do validation errors make an invalid bound or count easy to correct?
3. Are exact evaluation, search, stopping and resume clearly distinguished?
4. Does comparing two candidates explain the efficiency/power change and the
   constraints or numerical failures that matter?

Record that feedback before broadening the module. Automated and agent-operated
smoke checks are not a substitute for this real researcher usability review.
