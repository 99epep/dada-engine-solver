# Domestic Refrigerator Comparison Boundary

## Purpose

A domestic appliance comparison is a separate sizing problem, not a reuse of
the human-powered cooling-cell geometry. Changed reservoir temperatures alter
the useful pressure ratio, valve chronology, exchanger requirements, cylinder
ratio, kinematics and optimum operating speed.

Two comparisons are required:

1. evaluate an existing DADA machine away from its design point;
2. redesign a DADA machine for the same cabinet duty and test temperatures.

The second comparison is the fair technology comparison. The first measures
off-design robustness.

## Standards boundary

European Regulation (EU) 2019/2016 determines appliance energy consumption at
ambient temperatures of 289.15 K and 305.15 K (16 and 32 degrees Celsius), with
compartment temperatures interpolated to the targets in its Table 3. The fresh-
food target is 277.15 K (4 degrees Celsius) and the frozen-compartment target is
255.15 K (-18 degrees Celsius). Annual energy consumption uses the mean of the
daily results at the two ambient temperatures, subject to the regulation's load
and auxiliary terms.

Primary sources:

- [Commission Delegated Regulation (EU) 2019/2016](https://eur-lex.europa.eu/legal-content/EN/ALL/?uri=CELEX:32019R2016)
- [IEC 62552-2 household refrigerating appliance performance tests](https://webstore.iec.ch/en/publication/68116)

IEC 62552 contains the detailed appliance test procedure but is not freely
available in full. The project must not claim IEC conformity from the public
summary or the EU regulation alone.

## Energy boundaries

Published annual electrical consumption does not directly reveal cabinet heat
leak or refrigerating COP. A product label therefore cannot be converted into a
machine COP without a measured or modelled thermal load.

The solver keeps four quantities separate:

```text
cabinet cooling load
thermodynamic machine input
shaft or actuator input
total appliance input including auxiliaries
```

The simple cabinet model is:

```text
Q_cabinet = K_cabinet * (T_ambient - T_compartment) + Q_internal
```

`K_cabinet`, internal loads, door-opening loads, defrost energy, fan and pump
power must be supplied or measured. No default values are invented.

Annual electricity and instantaneous running power are connected through an
explicit duty fraction:

```text
P_average = duty * P_on + (1 - duty) * P_off
E_annual = P_average * hours_per_year
```

`AnnualEnergyBudget` and `annualize_cycling_appliance` implement this boundary.
An exploratory 200 kWh/year budget corresponds to 22.83 W continuous-equivalent
electrical input. It is a normalization point, not a claimed market average or
a substitute for a product-specific energy label.

The reversible COP at compartment and ambient boundary temperatures is useful
only as an upper bound. Real DADA reservoir temperatures require finite approach
differences below the compartment and above ambient, so their Carnot bound is
lower. The first comparison matrix will contain fresh-food and frozen cases at
both 16 and 32 degrees Celsius ambient temperatures.

## Required redesign procedure

For each test point:

1. define measured cabinet load and auxiliary boundary;
2. select explicit cold-side and hot-side approach temperatures;
3. optimize DADA cylinder ratio, displacement, charge, speed and kinematics;
4. size both exchanger geometries and hydraulic passages;
5. reject unavailable or violated Mach, pressure-drop and property criteria;
6. calculate total appliance input after transmission and auxiliaries;
7. compare service delivered, not thermodynamic COP alone.

Transient pull-down, door opening, automatic defrost and thermostat cycling are
later appliance-level models. They are not silently represented by the present
periodic steady-state cycle.

## Initial off-design screening

Candidate 001 was evaluated without redimensioning at a declared 5 K approach
on both sides. These are sensitivity points, not IEC appliance tests and not
optimized DADA geometries:

| Compartment | Ambient | Cold reservoir | Hot reservoir | Cooling | Input | COP |
|---|---:|---:|---:|---:|---:|---:|
| fresh food, 277.15 K | 289.15 K | 272.15 K | 294.15 K | 133.09 W | 36.54 W | 3.643 |
| fresh food, 277.15 K | 305.15 K | 272.15 K | 310.15 K | 73.50 W | 24.51 W | 2.998 |
| frozen, 255.15 K | 289.15 K | 250.15 K | 294.15 K | 33.58 W | 16.05 W | 2.093 |
| frozen, 255.15 K | 305.15 K | 250.15 K | 310.15 K | -26.19 W | 4.05 W | unavailable |

All four results have indeterminate physical validity because Mach remains
unavailable, and all have non-nominal passive-valve chronology. The warm frozen
case is not a refrigerator at all. This directly demonstrates that the
cooling-cell kinematics cannot be extrapolated to a domestic freezer and that a
temperature-specific redesign is required.

## First helium sensitivity pass

The inherited one-bar helium seed at 250.15/310.15 K gives 77.92 W cooling,
43.59 W thermodynamic input and COP 1.788 at 9.09 rpm. A one-at-a-time search,
with one-percent cylinder clearances, found:

- a local cylinder swept-volume ratio near `S/L = 0.90`;
- COP 1.809 at 9.09 rpm, with 76.91 W cooling;
- COP 2.000 at 6 rpm, with 56.83 W cooling;
- COP 2.155 at 2 rpm, with 21.95 W cooling and 75 W/K per exchanger;
- COP 2.127 at 2 rpm with only 35 W/K per exchanger and 20.24 W cooling.

The lower-speed point illustrates an energy-versus-size tradeoff, not a universal
optimum. With explicitly illustrative values of 85 percent transmission
efficiency and 2 W running auxiliaries, continuous operation would consume about
116 kWh/year. A real cabinet load, cycling losses and verified auxiliaries remain
required before this becomes an appliance prediction.

An exact first-level pressure-volume similarity then reduced every volume and
CdA by a factor of five while raising charge pressure from one to five bar. The
result is 5.94 and 6.60 litre swept volumes, 0.132 litre gas volume per exchanger,
unchanged 35 W/K conductances, and the same 20.24 W / 9.52 W / COP 2.127 result.
Maximum pressure is 9.78 bar absolute. This compact candidate remains
`indeterminate` because real hydraulic areas are not yet attached for Mach, and
its passive-valve chronology is non-nominal. Helium property and real-gas checks
are also required over its 217.9 to 361.2 K temperature range.
