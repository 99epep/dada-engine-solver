# Human Cell — Stage 0 notes

## Scientific purpose

This is the first coupled optimization study for the human-powered cooling cell.

The intended workflow is deliberately different from "ask the optimizer to make
150 W".  Stage 0 maximizes refrigeration COP and **observes**:

- indicated mechanical input power;
- cooling power;
- hot-side rejected heat;
- pressure / temperature / mass-flow margins;
- exchanger-domain margins.

The nominal speed is fixed at **2.8 Hz**, corresponding to the chosen
pedalling/transmission concept.

Machine scale is not optimized in this first campaign.  Once a thermodynamically
interesting region is found, later studies will rescale total swept volume,
working-gas inventory and exchanger capacity to move the indicated input toward
roughly 150 W, then re-optimize locally.

## Stage-0 active vector

Current file: 15 active coordinates:

- 10 `four_stage` kinematic coordinates:
  - SMALL: `t1`, `t2`, `t3`, `a_s`, `b_s`
  - LARGE: `t1`, `t2`, `t3`, `a_l`, `b_l`
- swept-volume ratio;
- H_i tube count and length;
- H_o tube count and length.

Frequency, total swept volume, charge, clearance ratios and external streams are
fixed.

## Important: historical 7D law versus current Research ownership

The historical `FourStageVolumeKinematics` law had seven independent parameters:

    t1, t2, t3, a_l, b_l, a_s, b_s

with `t1/t2/t3` shared by both cylinders.

Research V3 currently owns `four_stage` timing independently per cylinder.
Therefore a fully active Research study has ten kinematic coordinates, although
both sides start from the same timing.

This Stage 0 file is valid for the current architecture and intentionally accepts
that extra freedom.

If we decide that the exact historical 7D topology is scientifically important,
the clean solution is a small Research extension for shared abstract timing
coordinates.  Do not fake the equality by optimizing only one side while leaving
the other fixed.

## External temperatures

Stage 0 keeps the already validated V3 boundary:

- cold external inlet: 278.15 K;
- hot external inlet: 298.15 K.

This is a mild refrigeration benchmark, **not yet a literal ice-production
boundary condition**.  It is useful for finding and debugging the coupled optimum
before lowering the cold-side temperature.

After Stage 0 is stable, make a new scientific study for the lower-temperature
water/ice scenario.  Do not edit an existing campaign definition in place.

## Suggested first commands

From the source checkout, assuming this directory is copied to
`outputs/human_cell_stage0/`:

```sh
research() { PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m dada_solver.research "$@"; }

research validate outputs/human_cell_stage0/study.toml

# Evaluate the declared initial point before starting Sobol:
research evaluate outputs/human_cell_stage0/study.toml \
  --output outputs/human_cell_stage0/initial.json \
  --budget 5m

# If it converges and remains in-domain:
research run outputs/human_cell_stage0/study.toml \
  --directory outputs/human_cell_stage0/campaign \
  --budget 30m \
  --max-candidates 8

research status outputs/human_cell_stage0/campaign
research report outputs/human_cell_stage0/campaign \
  --html outputs/human_cell_stage0/campaign.html
```

If the 2.8 Hz initial point is invalid, inspect the reported limiting constraint
before changing the study.  In particular, do not automatically relax the
microtube-domain checks merely to obtain a result.

## Next stages

1. Find a credible COP basin at 2.8 Hz.
2. Inspect whether independent S/L timings actually separate significantly.
3. If they remain close, consider adding the exact shared 7D timing model.
4. Move the cold boundary toward the intended ice-production temperatures.
5. Rescale machine / charge / exchanger capacity toward ~150 W indicated input.
6. Re-optimize around each new scale.
7. Only then compare `four_stage` against structured C2 15p and later physical
   4-/6-bar mechanisms.
