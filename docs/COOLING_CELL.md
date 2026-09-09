# Human-Powered Cooling Cell Reference

## Ideal external load

The first external load is one kilogram of water cooled from 293.15 K to
273.15 K and fully frozen near 273.15 K. The simple load closure is:

```text
E_load = m * cp_liquid * (T_initial - T_target)
       + m * latent_energy * phase_change_fraction
```

The controlled configuration supplies the approximate properties explicitly:

```text
cp_liquid = 4180 J/(kg K)
latent_energy = 333500 J/kg
```

This gives:

| Quantity | Value |
|---|---:|
| Sensible cooling | 83.6 kJ |
| Solidification | 333.5 kJ |
| Ideal total | 417.1 kJ |
| Mean cooling power for 20 min | 347.6 W |
| Mean cooling power for 15 min | 463.4 W |
| Required COP at 150 W for 20 min | 2.32 |
| Required COP at 150 W for 15 min | 3.09 |

These are ideal load values. They exclude the container, insulation, ambient
heat leaks, contact resistances, temperature gradients, supercooling and any
inefficiency between the DADA cold exchanger and the water.

## Human-power interpretation

The exploratory human input range is 75-250 W with 150 W nominal. At 150 W,
the machine operating point must simultaneously provide sufficient cooling
power and require no more than 150 W mechanical input.

For a fixed simulated operating point:

```text
ideal_time = ideal_load_energy / machine_cooling_power
```

This estimate is reported only for a refrigerating point with positive cooling
power. If the simulated point requires more than the available human power, it
is marked inoperable at that power. Cooling power is not scaled linearly; speed,
pressure or another operating variable must be changed and the machine must be
simulated again.

`required_pedaling_power_for_target_time` divides the task cooling power by the
COP of the evaluated point. It is a local operating-point diagnostic, not a
claim that COP remains constant when the machine is rescaled.

## Exploratory machine bounds

The user-selected provisional bounds convert to SI as follows:

| Parameter | Lower | Upper |
|---|---:|---:|
| Angular speed | 3.1416 rad/s | 157.0796 rad/s |
| Charge pressure | 100000 Pa | 2000000 Pa |
| Swept volume, each cylinder | 1e-5 m³ | 5e-3 m³ |
| Heat-exchanger volume, each | 1e-6 m³ | 2e-3 m³ |
| `UA_C`, `UA_H` | 0.1 W/K | 500 W/K |
| Each effective `CdA` | 1e-8 m² | 1e-3 m² |

These bounds are exploratory, not validated physical limits.

## Selected exploratory machine boundary conditions

> **Withdrawn results:** all piecewise-linear numerical results recorded before
> correction of the angular origin used the wrong theta-zero filling volume and
> therefore the wrong total gas inventory. They are retained below only as an
> audit trail and must not be used for sizing or physical conclusions. The
> corrected law starts at `V_L = V_L_max` as required by the study.

The current machine sensitivity pass uses a user-selected constant cold-reservoir
temperature of 268.15 K. This is distinct from the water target temperature of
273.15 K. It supplies a 5 K nominal external approach at the end of freezing,
but the gas in `C` must still become colder than 268.15 K for positive heat
transfer under the finite-`UA` closure. Lower cold temperature may increase a
thermal driving difference while reducing thermodynamic COP.

Cylinder clearance is parameterized separately for each cylinder by
`V_min / V_swept`, with exploratory bounds from 0.005 to 0.10. Heat-exchanger
volumes are separate control volumes and are never included in this cylinder
clearance ratio. Pipe volumes remain unmodeled.

One-at-a-time studies can be run before constrained optimization. The selected
ideal reference uses piecewise-linear volume laws, a common Lambda target, and
initially four equal angular sectors. The following command varies the
duration assigned to each adiabatic sector while holding all other configured
variables at their initial values:

```console
python -m dada_solver.sizing.cli examples/cooling_cell_sensitivity_space.toml \
  --sensitivity adiabatic_sector_fraction --sensitivity-points 5
```

Each transfer-sector fraction is `0.5 - adiabatic_sector_fraction`. The current
0.05 to 0.25 interval is diagnostic, not a validated mechanism design range.

The first cardinal-phase sweep at the initial exploratory point converged at
all five samples from 0 to 2 pi. Nevertheless, cold-reservoir heat-transfer
power remained negative, from approximately -38 W to -107 W. Every sampled
cycle was non-nominal and invalid under the configured validity thresholds.
Consequently, this is not a refrigeration operating point and no cooling COP or
freezing-time claim is available. Scaling this point is not a valid sizing step.

With the piecewise-linear law and common Lambda target 0.7, the equal-sector
point also converged but produced approximately -32.1 W of signed cold heat
transfer and required 75.7 W of mechanical input. Shortening each adiabatic
sector from 0.25 to 0.05 of a cycle, thereby lengthening each transfer sector
from 0.25 to 0.45, improved these values monotonically over the five sampled
points to approximately -19.6 W and 51.6 W. This locally supports the proposed
direction of the duration effect, but does not reverse the cold-heat sign. All
five points remained non-nominal and invalid.

After raising the cold reservoir by 5 K to 268.15 K, a common-Lambda sweep from
0.70 to 0.95 with equal sectors improved signed cold heat-transfer power from
approximately -30.7 W to a local best of -23.8 W at Lambda 0.90, before a small
degradation at 0.95. Combining Lambda 0.90 with the duration sweep gave a local
best of approximately -20.1 W at an adiabatic-sector fraction of 0.15. Longer
transfers beyond this point degraded the result. No sampled combination was a
refrigeration point; all remained non-nominal and invalid.

Warning: a later Lambda sweep initially described as using the doubled-cylinder,
reduced-HX geometry was invalidated because sizing initial values overrode the
updated base geometry. Its reported values were withdrawn. Sensitivity files
must keep every non-swept initial design value synchronized with the intended
base point; this is now checked explicitly for the cooling-cell configuration.
