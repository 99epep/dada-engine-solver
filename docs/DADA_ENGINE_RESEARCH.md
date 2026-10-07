# DADA Research user guide

Research configures, evaluates, searches and compares machines using the existing
physical solver. Indicated gas power is not automatically useful mechanical
shaft power: mechanical losses remain unknown.

From the repository root, define the shorthand used throughout this guide:

```sh
research() { PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m dada_solver.research "$@"; }
```

The installed command is `dada-research`; use it instead of `research` outside a
source checkout. All paths below are illustrative destinations for new work.
Use `research COMMAND --help` for the complete command options.

## Create a study

Choose one of the current starting templates:

```sh
research init kinematics --small slider_crank --large harmonic \
    --output outputs/my_study/study.toml
research init external-stream-refrigeration \
    --output outputs/my_refrigerator/study.toml
research init external-stream-motor \
    --output outputs/my_motor/study.toml
```

Templates start with fixed coordinates. Edit the generated study for the intended
machine and operating conditions. A fixed parameter has a `value`; to search a
coordinate, replace that value with an active declaration including its initial
value and bounds. Do this before `run`, which needs at least one active parameter.

See the [technical reference](DADA_ENGINE_RESEARCH_REFERENCE.md) for parameter
types, policies, objectives and constraints, and the
[kinematics reference](DADA_ENGINE_RESEARCH_KINEMATICS.md) for family selection.

## Validate and evaluate one point

```sh
research validate outputs/my_study/study.toml
research evaluate outputs/my_study/study.toml \
    --output outputs/my_study/reference.json --budget 3m
```

`validate` checks declarations, hashes and production geometry without integrating
thermodynamics. `evaluate` integrates one exact configured point, using active
initial values where applicable; it does not run a search. It never overwrites an
existing result, so choose a new output filename for another evaluation.

## Run, resume and inspect

After declaring active parameters:

```sh
research run outputs/my_study/study.toml --budget 30m --max-candidates 32
research resume outputs/my_study/campaign --budget 2h --max-candidates 128
research status outputs/my_study/campaign
```

`run` creates `campaign` beside the study TOML by default. Use `--directory` to
choose another destination. A nonempty destination is refused; use `resume` to
continue it. Budget and `--max-candidates` bound the current invocation.

`resume` uses the campaign's stored snapshots, not later edits to the original
TOML. It also accepts a TOML path as shorthand for its sibling `campaign` directory.
Resume checks execution compatibility; a missing campaign does not start a new
search. `status` inspects stored results without integrating thermodynamics.

For persistence and compatibility rules, see the
[technical reference](DADA_ENGINE_RESEARCH_REFERENCE.md). Progress display,
interruption recovery and result evidence are covered in the
[cockpit guide](DADA_ENGINE_RESEARCH_COCKPIT.md).

## Reports, curves and animations

```sh
research report outputs/my_study/campaign
research report outputs/my_study/campaign --candidate second
research report outputs/my_study/campaign --candidate best --plots default,flows,mechanisms
research report outputs/my_study/campaign --plots none
research compare outputs/my_study/campaign outputs/my_motor/campaign \
    --html outputs/comparison.html
```

By default, `report` writes `campaign/report.html` and requests positions, heat
exchange, pressures and temperatures for the two best feasible candidates.
`--candidate` also accepts an exact ID or unique prefix; repeat it to select
several candidates. `--plots none` omits curves, while `mechanisms` requests
animations for supported physical mechanisms. Use `--html PATH` for another HTML
destination.

Requested thermodynamic curves are reconstructed from saved states; Research
announces the thermodynamic replay when needed. The HTML is a derived artifact.
Comparison can inspect several campaigns or saved evaluations without inventing
a joint ranking for incompatible studies. Plot choices, evidence interpretation
and animation details belong in the [cockpit guide](DADA_ENGINE_RESEARCH_COCKPIT.md).

## Refine a promising result

```sh
research refine outputs/my_study/campaign --candidate best --radius 0.20 \
    --output outputs/local/study.toml
research validate outputs/local/study.toml
```

`refine` creates a portable schema-3 study around selected stored candidates;
it does not launch a campaign or copy their results as completed evaluations.
Sources must have the same scientific study identity and an active parameter
space. Repeat `--candidate` to retain several centers. The radius is a fraction
of the original normalized search width, not of each physical parameter value.

Review the generated bounds and then use `run` when ready. Exact scheduling and
radius semantics are in
[Local refinement and explicit initial evaluations](DADA_ENGINE_RESEARCH_REFERENCE.md#local-refinement-and-explicit-initial-evaluations).

## Change capacity

```sh
research rescale outputs/my_study/campaign --candidate best --factor 2 \
    --output outputs/larger/study.toml
research validate outputs/larger/study.toml
```

`rescale` creates a new portable schema-3 study without integrating or modifying
the source. It preserves declared constraints and never overwrites destination
files. Global and local candidates are accepted; source local basins are not reused
in the new study. Review the result before evaluation or search.

The [capacity reference](DADA_ENGINE_RESEARCH_CAPACITY.md) defines which
quantities scale and which remain fixed.

## Further reading

- [Research technical reference](DADA_ENGINE_RESEARCH_REFERENCE.md): schema,
  scientific identity, policies, search, persistence and runtime compatibility.
- [Cockpit and reports](DADA_ENGINE_RESEARCH_COCKPIT.md): progress, recovery,
  result inspection and evidence presentation.
- [Kinematics and mechanisms](DADA_ENGINE_RESEARCH_KINEMATICS.md): motion families,
  physical mechanisms and mechanical diagnostics.
- [Capacity scaling](DADA_ENGINE_RESEARCH_CAPACITY.md): transformation rules.
- [Current validation and limitations](validation.md): what the available
  software and physical evidence establishes.
