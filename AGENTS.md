# Repository instructions

## Language and scope

Write repository content in English: documentation, comments, docstrings,
identifiers, messages, reports and figures. Preserve scientific symbols,
SI units, proper names and original source URLs. Conversation may be in French.

Current user instructions take precedence over historical project decisions.
Study definitions own objectives, constraints, operating conditions, parameter
bounds and evaluation policies. Do not promote a campaign's requirements or
research roadmap into universal defaults. Do not launch optimization campaigns
or thermodynamic calculations merely to edit or validate configuration files.

## Scientific and software boundaries

Before changing the physical model, read [physical decisions](docs/PHYSICS_DECISIONS.md)
and the relevant model reference listed in the [documentation index](docs/README.md).
Disclose assumptions, evidence and model limitations; distinguish numerical
validation from experimental validation and best-found results from proven optima.

Preserve production kinematics and exchanger interfaces, family-specific parameter
ownership and explicit validity margins. Reuse production geometry and validation
rather than implementing independent approximations. Do not independently vary
quantities derived from exchanger geometry. Do not impose symmetry or shared
mechanical families unless the study explicitly declares that coupling.

Preserve scientific hashes, artifact constraints, portable input links and existing
provenance. Campaign history is append-only; preserve exact identity caching,
rejection states and resume compatibility. Do not silently weaken convergence or
feasibility checks, change scientific inputs, or reinterpret historical artifacts.

## Research command cheat sheet

Run from the repository root. The installed equivalent is `dada-research`.
Paths, budgets, radii and factors below are illustrative, not scientific defaults.
Use `research COMMAND --help` for all options and the
[Research guide](docs/DADA_ENGINE_RESEARCH.md) for detailed workflows.

```sh
research() { PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m dada_solver.research "$@"; }
```

### Prepare studies without thermodynamic integration

```sh
# Validate configuration, geometry and linked inputs.
research validate study.toml

# Edit interactively (requires the optional tui extra).
research study edit study.toml

# Reopen mechanisms while retaining existing active hardware parameters.
research study release study.toml --group mechanisms --radius 0.05 --output work/combined/study.toml
research study freeze study.toml --group frequency --output work/fixed/study.toml

# Narrow declared mechanism bounds; this does not change the scheduler radius.
research study bounds recenter study.toml --group mechanisms --half-width 0.01 --output work/narrow/study.toml

# Create a local study around a selected stored candidate.
research refine campaign --candidate best --radius 0.1 --tidy --output work/refined/study.toml

# Explicit capacity scaling, not local refinement or hardware retuning.
research rescale campaign --candidate CANDIDATE_ID --factor 2 --output work/rescaled/study.toml

# Portable export; standalone intentionally removes historical provenance.
research study export study.toml --tidy --standalone --output work/standalone/study.toml
```

For `refine`, the radius is relative to the declared normalized domains, not a
percentage of each physical value. Study-editor local radius settings apply to
all active parameters. `rescale` changes scientific inputs under the capacity
scaling policy; validate and evaluate the resulting study before comparing it.
Keep associated basis/artifact files with generated TOMLs and use new destinations.

### Evaluate and search (thermodynamic calculations)

Run these only within the requested task and budget. `evaluate` evaluates the
study's initial configuration unless explicit overrides are supplied; `run`
requires active parameters. Resume uses the stored campaign definition.

```sh
research evaluate study.toml --budget 5m --output work/initial-evaluation.json
research run study.toml --directory work/campaign --budget 30m --max-candidates 32
research resume work/campaign --budget 30m --max-candidates 32
```

### Inspect results and mechanisms

`status` reads stored results. Reports and comparisons with cycle plots can replay
one cycle from a saved state; they are not necessarily calculation-free.
`best` follows the source study's objectives and constraints. Use an exact candidate
ID or unique prefix when a specific result is required.

```sh
research status campaign --list-candidates
research report campaign --candidate best
research report campaign --candidate best --plots mechanisms --external-webp
research report campaign --plots none
research compare campaign-a campaign-b --candidate best --html work/comparison.html

# Browse saved geometry only, with no thermodynamic integration.
research mechanism visualize library.json --side large --browse
research mechanism visualize library.json --side large --browse --sort position-rms --target target.json
```

Position-RMS sorting concerns complete piston mechanisms; primary E projections
are not piston fits. For mechanical synthesis, pairing, adaptation and retuning,
consult the [kinematics workflows](docs/DADA_ENGINE_RESEARCH_KINEMATICS.md).

### Motion refit and mechanical synthesis (no thermodynamic integration)

`motion-refit` initializes the unique destination `structured_c2_15p`; it does not
run its resulting study. `motion-target` can extract motion directly from a hybrid,
structured, spline or other compatible source; a spline is not a prerequisite.
Replace uppercase IDs below with actual IDs from each produced library/catalogue;
polish and mirror steps may generate new IDs rather than retaining the parent ID.

```sh
research motion-refit path/to/hybrid/campaign --candidate best \
  --output work/structured/study.toml --report work/structured/refit.json --plot
research motion-target path/to/source/campaign --candidate best --output work/target.json

research mechanism synthesize work/target.json --family slider_crank \
  --stage global_discovery --side both --seed 1234 --islands 4 --budget 2m \
  --output work/slider-library.json
research mechanism synthesize work/target.json --family four_bar \
  --stage global_discovery --side both --seed 1234 --islands 8 --budget 2m \
  --output work/fourbar-library.json
research mechanism synthesize work/target.json --family four_bar \
  --stage full_local_polish --side large --library work/fourbar-library.json \
  --family-id LARGE_FAMILY_ID --output work/fourbar-polished.json
research mechanism visualize work/fourbar-polished.json \
  --family-id POLISHED_FAMILY_ID --side large --target work/target.json
research mechanism visualize work/fourbar-library.json \
  --side small --target work/target.json --browse --sort position-rms

research mechanism synthesize work/target.json --family six_bar \
  --stage primary_discovery --side large --islands 12 --budget 5m \
  --output work/six-primary.json
research mechanism synthesize work/target.json --family six_bar \
  --stage downstream_fit --side large --library work/six-primary.json \
  --family-id PRIMARY_ID --budget 5m --output work/six-downstream.json
research mechanism synthesize work/target.json --family six_bar \
  --stage full_local_polish --side large --library work/six-downstream.json \
  --family-id DOWNSTREAM_ID --output work/six-polished.json
research mechanism synthesize work/target.json --family six_bar \
  --stage mirror_initialization --side small --library work/six-polished.json \
  --family-id POLISHED_ID --output work/six-small-seed.json
research mechanism synthesize work/target.json --family six_bar \
  --stage opposite_local_adaptation --side small --library work/six-small-seed.json \
  --family-id MIRROR_SEED_ID --output work/six-pair.json
```

Six-bar primary discovery retains intermediate E trajectories, not complete piston
mechanisms. Its downstream, polish and opposite-side stages remain distinct.
The first side can be SMALL or LARGE. Mirroring only creates an independent seed.

### Select a pair, adapt it, then retune hardware

Pair selection is explicitly human, with no automatic ranking of combinations.
Select complete SMALL and LARGE artifacts; their physical families may differ.
Choose either the common-library or separate-library form:

```sh
research mechanism pair work/slider-library.json \
  --small SMALL_ID --large LARGE_ID --output work/selected-pair.json
research mechanism pair \
  --small-library work/slider-library.json --small SMALL_ID \
  --large-library work/fourbar-library.json --large LARGE_ID \
  --output work/mixed-pair.json

research mechanism adapt path/to/source/campaign --candidate SOURCE_ID \
  --library work/selected-pair.json --family-id PAIR_ID \
  --output work/paired/study.toml --radius 0.1
research validate work/paired/study.toml
research run work/paired/study.toml --budget 30m --max-candidates 32

research mechanism retune work/paired/campaign --candidate best \
  --scope source-active --radius 0.1 --output work/retuned/study.toml
research validate work/retuned/study.toml
research run work/retuned/study.toml --budget 30m --max-candidates 32

# On a study where both mechanisms and hardware are already active:
research study bounds recenter path/to/active/study.toml \
  --group mechanisms --half-width 0.01 --radius 1 \
  --output work/combined-narrow/study.toml
```

`pair`, `adapt` and `retune` only prepare portable inputs; the `run` commands above
perform thermodynamic searches and require task authorization. `adapt` frees the
selected mechanical coordinates with hardware fixed. `retune` fixes the adapted
mechanisms and reopens selected original hardware domains. Use the exact source
candidate intended for the comparison, independently of the motion target's family.
A complete pair already present in a library can go directly to `adapt`.

In the last command, `half-width` narrows only the selected declared domains;
`--radius 1` separately applies to ALL active parameters of the local scheduler.
It does not equalize physical sensitivities. If replacing incompatible existing
local regions is required, review the change before explicitly adding `--recenter`.

## References and checks

Use the [Research reference](docs/DADA_ENGINE_RESEARCH_REFERENCE.md) for study
schema and policies, [kinematics reference](docs/DADA_ENGINE_RESEARCH_KINEMATICS.md)
for family conventions, and [synthesis method](docs/MECHANISM_SYNTHESIS_SEARCH.md)
for mechanical search stages. For architecture and persistence, consult
[model interfaces](docs/PLUGGABLE_MODELS.md) and
[campaign internals](docs/OPTIMIZATION_CAMPAIGN.md) when relevant.
Keep scientific details in these references rather than duplicating them here.

Run relevant checks from the source checkout with
`PYTHONPATH=src python3 -m pytest`. Use the [validation map](docs/validation.md)
to select checks appropriate to the change. Check documentation links after
editing documentation and run `git diff --check` before delivery.
