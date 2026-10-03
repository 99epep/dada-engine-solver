# Historical studies and evidence

These documents preserve results under their **original** definitions, source
versions, domains and study constraints. They are not a recommended parameter
set for a new machine. Current guidance starts at the [documentation index](../README.md).
Old test counts are dated records, not the present test status.

## Condensed records

- [Decision ledger](PHYSICS_DECISION_LEDGER.md): original decisions, reversals and
  withdrawn early results. [Current decisions](../PHYSICS_DECISIONS.md) take precedence.
- [Early exchanger trials](EXCHANGER_TRIALS.md): finite-air/wall, geometry, doubled
  banks, excess air, temperature and constant-area parallel-tube comparisons.
- [Numerical execution](../SOLVER_PERFORMANCE.md): consolidated acceleration
  evidence, adopted optimizations and negative interpolation/Anderson experiments.
- [Documentation map](DOCUMENTATION_MAP.md): old names, destinations and local
  Git provenance for the full pre-consolidation prose.

## Motor and mechanism studies retained in detail

| Record | Why retain it |
|---|---|
| [Four-stage optimization](../FOUR_STAGE_OPTIMIZATION.md), [K2 search](../FOUR_STAGE_K2_SEARCH.md) | Distinguishes artificial conductance sensitivity from actual hardware changes and different motion basins |
| [Valve screening](../VALVE_PLACEMENT_SCREENING.md), [temperature/topology campaign](../MOTOR_TEMPERATURE_AND_VALVE_CAMPAIGN.md) | Matched DD/UD/DU/UU comparisons and temperature-specific evidence |
| [Four-bar comparison](../MOTOR_FOUR_BAR.md), [champion integration](../MOTOR_CHAMPION_FOUR_BAR.md), [six-bar K2](../SIX_BAR_K2_EVALUATION.md) | Motion substitution at fixed thermodynamic inputs and sign/scale conventions |
| [Initialization synthesis](../MOTOR_KINEMATIC_INITIALIZATION_SYNTHESIS_V1.md) | Explicit confidence levels and transferable versus gas/application-specific hypotheses |
| [Compact synthesis](../COMPACT_KINEMATIC_SYNTHESIS.md) | Finite rods, output frames, compactness and branch lessons |
| [Primary families](../PRIMARY_FOUR_BAR_FAMILIES.md), [six-bar families](../SIX_BAR_MECHANISM_FAMILIES.md) | Exact geometry, discrete branches, source artifacts and robustness tradeoffs |

The [synthesis theory](../MECHANISM_SYNTHESIS_SEARCH_THEORY.md) remains a primary
methodological reference, not discarded chronology. Its rejected proxies and
hierarchical search evidence matter for future human-cell mechanisms.

## Research implementation provenance

[Initial architecture audit](../DADA_ENGINE_RESEARCH_ARCHITECTURE_AUDIT.md),
[implementation record](../DADA_ENGINE_RESEARCH_IMPLEMENTATION_PLAN.md),
[migration matrix](../DADA_ENGINE_RESEARCH_MIGRATION_MATRIX.md) and
[validation record](../DADA_ENGINE_RESEARCH_VALIDATION.md) retain source extraction,
parity cases and compatibility decisions. Use the current Research guide for
commands; early limitations in dated sections do not describe the current CLI.
Machine-readable `research_audit/` inventories remain intact.

Results under `outputs/` and `examples/` are not relocated, regenerated, cleaned
or relabelled by this documentation pass. A saved feasible verdict may include
historical guards now absent from generic presets. Exact reproduction needs the
matching source/runtime and input definition, not just a matching candidate name.
