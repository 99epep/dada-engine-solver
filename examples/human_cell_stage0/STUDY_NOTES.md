# Human Cell — Stage 0

This is the first coupled refrigeration optimization study.

## Fixed operating conditions

- machine frequency: **2.8 Hz**;
- cold external inlet: **273.15 K (0 °C)**;
- hot external inlet: **298.15 K (25 °C)**;
- objective: **maximize cooling COP**;
- indicated mechanical input power: **diagnostic only**.

## Active coordinates

The study has **14 active coordinates**:

- 9 coordinates from `hybrid_compact`;
- swept-volume ratio;
- H_i tube count and tube length;
- H_o tube count and tube length.

The hybrid initial point is historical candidate 501.

## Deliberately fixed for Stage 0

- total swept volume;
- working-gas inventory;
- clearance ratios;
- frequency;
- external liquid-like stream mass flows and cp;
- external wall conductances.

The first goal is to find a good coupled kinematic/thermodynamic basin at the
actual intended temperature lift.  Only afterwards will machine scale, gas
inventory and exchanger capacity be moved progressively toward roughly **150 W**
of indicated mechanical input, with re-optimization after each scale change.

## Why no power constraint yet

Cooling COP is:

    COP = cooling_power / indicated_mechanical_input_power

Cooling power and indicated input remain visible in every result.  If the COP
optimum collapses toward negligible power, that is a result to diagnose rather
than something hidden by an arbitrary initial power constraint.

## First run

Place this directory at:

    outputs/human_cell_stage0/

Then run:

```sh
research() { PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m dada_solver.research "$@"; }

research validate outputs/human_cell_stage0/study.toml

research evaluate outputs/human_cell_stage0/study.toml   --output outputs/human_cell_stage0/initial.json   --budget 5m

research report outputs/human_cell_stage0/initial.json   --html outputs/human_cell_stage0/initial.html
```

If the declared initial point converges and remains valid:

```sh
research run outputs/human_cell_stage0/study.toml   --directory outputs/human_cell_stage0/campaign   --budget 30m   --max-candidates 16

research status outputs/human_cell_stage0/campaign

research report outputs/human_cell_stage0/campaign   --html outputs/human_cell_stage0/campaign.html
```

Resume with:

```sh
research resume outputs/human_cell_stage0/campaign   --budget 30m   --max-candidates 16
```

## What to inspect

After the first phase, inspect COP together with:

- cooling power;
- indicated mechanical input;
- cold/hot external heat rates;
- pressure and temperature extrema;
- maximum mass flow;
- exchanger-domain validity;
- local reflux;
- constraint margins;
- coordinates accumulating on search bounds.

The next search box should be chosen from those diagnostics.
