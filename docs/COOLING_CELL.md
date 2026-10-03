# Human-powered cooling-cell boundary

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

The historical exploratory human-input scenario was 75–250 W with 150 W nominal.
These are study choices, not physiological guarantees or solver defaults. The
solver reports indicated input; actual shaft input requires a separate loss model.
A 150 W indicated-input constraint must be met at the same operating point as
the cooling duty; it does not certify operation with 150 W available at the pedals.

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

## Configure the machine separately from the load

Use [Research](DADA_ENGINE_RESEARCH.md) for the current machine definition:
external cold/hot streams, kinematics, geometry, charge policy, frequency,
objective and explicit constraints. Either maximize cooling COP or cooling power
under a declared indicated-input bound. The load model does not choose between
them or infer a human mechanical efficiency.

A finite external stream is not a constant-temperature water reservoir.
Its fluid label selects no liquid correlation. Declare inlet temperature, flow,
Cp and conductance; document heat leaks, tank/container load, coupling and any
excluded pump/transmission losses separately. See
[external streams](EXTERNAL_STREAM_THERMAL_MODEL.md).

Atmospheric reference filling uses the assembled machine's maximum simultaneous
gas volume, including exchanger hold-up once. Cylinder clearance excludes that
separate exchanger volume. Circular collector geometry changes hold-up and the
derived charge; an old capacity-similarity result is not exact proof for that
new geometry. See [capacity scaling](DADA_ENGINE_RESEARCH_CAPACITY.md).

## Historical results that must not be reused

Early piecewise-linear sensitivity results were withdrawn because the angular
origin gave the wrong filling volume/inventory. Another sweep was withdrawn
because sizing initials overrode the edited base geometry. Their numerical
listings are omitted here; the reasons and original record remain in the
[decision ledger](history/PHYSICS_DECISION_LEDGER.md). Old examples and outputs
are unchanged and are not promoted to current recommended studies.

The [domestic refrigerator comparison](DOMESTIC_REFRIGERATOR.md) has a different
load, temperature and auxiliary boundary. It requires its own design study.
