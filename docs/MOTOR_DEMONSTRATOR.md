# Motor demonstrator design brief

> Application-specific design requirements. These limits belong to the motor
> demonstrator study, not to generic Research or intrinsic model domains.

The intended machine is a small experimental demonstrator, not a commercial
product or a validated engine. Its name is independent of source temperature.

## Design requirements

| Quantity | Requirement |
|---|---|
| Useful mechanical output | Approximately 100 W; no acceptance tolerance specified |
| Cycle frequency | 2–10 Hz |
| Heat-input air inlet | 598.15 K (325 °C) |
| Heat-rejection air inlet | 298.15 K (25 °C) |
| Large-cylinder maximum enclosed volume | At most 0.066 m³ (66 L), including clearance |
| Primary objective | Maximize efficiency subject to power and size requirements |
| External source medium | Air with finite flow |

The 66 L ceiling is not a target displacement or a total-machine volume limit.
It applies only when explicitly declared by the demonstrator study; generic
Research does not impose it. Report exchanger hold-up and packaging separately.
After improving dimensions, the project brief calls for progressively reducing
the hot-source temperature; each study must state its actual inlet conditions.

The solver reports indicated gas power. With an explicit system boundary,
`useful_power = indicated_power - mechanical_losses - auxiliary_power`, without
double counting. Mechanical losses remain unknown: **a 100 W indicated result
does not establish the 100 W useful-output requirement**. External-air aerodynamic
losses and fan consumption are excluded from the current trial balance, not
assumed physically zero.

Prefer documented commercial exchanger hardware. Where suitable data are
unavailable, Doty-based experimental extrapolation is authorized with explicit
assumptions and uncertainty. This policy does not qualify any particular product;
thermal performance, pressure loss, gas hold-up and operating limits require
traceable evidence for the selected hardware.

## Interpretation

These are application requirements, not correlation domains or universal solver
limits. Evaluate the actual realizable motion and assembled hardware under the
declared study conditions. Numerical convergence, a kinematic fit or a screening
correlation does not establish experimental machine performance or useful shaft
output. Model applicability and evidence boundaries remain separate checks.

## Current model references

- [Physical decisions](PHYSICS_DECISIONS.md): model boundaries and limit ownership.
- [Research reference](DADA_ENGINE_RESEARCH_REFERENCE.md): study configuration and constraints.
- [Kinematics](DADA_ENGINE_RESEARCH_KINEMATICS.md): realizable motion and mechanisms.
- [Microtube model](MICROTUBE_GAS_MODEL.md): working-gas exchange and geometry.
- [External streams](EXTERNAL_STREAM_THERMAL_MODEL.md): finite-air and wall model.
- [Exchanger validation](EXCHANGER_VALIDATION.md), [literature](EXCHANGER_LITERATURE.md)
  and [Doty screening](DOTY_SCREENING.md): empirical evidence and extrapolation limits.
- [Motor research objectives](MOTOR_RESEARCH_OBJECTIVES.md): ordered research questions.
- [Validation](validation.md): numerical evidence and physical limitations.
