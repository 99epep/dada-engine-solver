# Documentation

Start with the [Research guide](DADA_ENGINE_RESEARCH.md). It begins with the
source-checkout shell helper and covers the daily commands. This index separates
current model references from dated experiments: an old candidate, limit or
benchmark is not a default for a new study.

## Use Research

| Need | Read |
|---|---|
| Configure, evaluate, run, resume, refine, rescale | [Research guide](DADA_ENGINE_RESEARCH.md) |
| Full schema, policies, persistence and compatibility | [Research reference](DADA_ENGINE_RESEARCH_REFERENCE.md) |
| Families, mechanism artifacts and mechanical constraints | [Kinematics](DADA_ENGINE_RESEARCH_KINEMATICS.md) |
| Interpret candidate evidence, margins and reports | [Cockpit](DADA_ENGINE_RESEARCH_COCKPIT.md) |
| Change machine capacity | [Capacity scaling](DADA_ENGINE_RESEARCH_CAPACITY.md) |
| Current validation | [Validation and evidence](validation.md) |

## Understand the physical model

| Topic | Reference |
|---|---|
| State, signs, charge, limits and omissions | [Physical decisions](PHYSICS_DECISIONS.md) |
| Motor direction and branch conventions | [Motor operation](MOTOR_OPERATION.md) |
| Microtube geometry, correlations and domains | [Microtube model](MICROTUBE_GAS_MODEL.md) |
| Finite external stream, wall storage and refrigeration COP | [External streams](EXTERNAL_STREAM_THERMAL_MODEL.md) |
| EOS, conservative reconstruction and compiled property tables | [Working fluids](WORKING_FLUID_MODELS.md) |
| Other exchanger screening/fixed-point interfaces | [Exchanger screening](HEAT_EXCHANGERS.md) |
| Empirical provenance and limits | [Exchanger validation](EXCHANGER_VALIDATION.md), [literature](EXCHANGER_LITERATURE.md), [Doty screening](DOTY_SCREENING.md) |
| Model domain versus design requirement | [Current physics decisions](PHYSICS_DECISIONS.md) |

## Design motion and mechanisms

The established method is **primary families → downstream dyad → local six-bar
polish → opposite-piston adaptation → paired thermodynamic evaluation → hardware
retuning**. Do not replace it with an acceleration-only fit or a blind global
15-dimensional search.

- [Mechanism synthesis search](MECHANISM_SYNTHESIS_SEARCH.md): methodology,
  failed proxies, hierarchical stages and saturation versus exploitation.
- [Primary four-bar families](PRIMARY_FOUR_BAR_FAMILIES.md) and
  [six-bar families](SIX_BAR_MECHANISM_FAMILIES.md): retained geometry, conventions,
  reproducible data and nominal mechanical margins.
- [Research kinematics](DADA_ENGINE_RESEARCH_KINEMATICS.md): production mechanisms,
  branches, physical scale and mechanical diagnostics.
- [Motion optimality](MOTION_OPTIMALITY.md): best-found evidence versus proof.
- [Motor objectives](MOTOR_RESEARCH_OBJECTIVES.md): ordered research questions.

## Develop and validate

- [Model interfaces](PLUGGABLE_MODELS.md): injection and ownership boundaries.
- [Campaign internals](OPTIMIZATION_CAMPAIGN.md): identity, recovery and scheduling.
- [Sizing API](SIZING.md): lower-level objectives, constraints and mechanical loads.
- [Numerical execution](SOLVER_PERFORMANCE.md): Python/Numba, exact reuse,
  diagnostic replay, benchmark evidence and rejected acceleration experiments.
- [Validation map](validation.md): tests to run and the meaning of their results.

## Application boundaries and historical evidence

- [Cooling cell](COOLING_CELL.md): load calculation and indicated-input boundary.
- [Domestic refrigerator](DOMESTIC_REFRIGERATOR.md): separate appliance boundary;
  dated standards/product references require review for a new comparison.
- [Motor demonstrator](MOTOR_DEMONSTRATOR.md): application-specific design brief,
  not global solver limits.

Saved `outputs/`, examples, audit JSON/CSV, fixtures and source artifacts are not
edited by documentation maintenance. Some large output directories are local-only
and may be absent from a distributed checkout. Their absence is not a reason to
invent replacement results.
