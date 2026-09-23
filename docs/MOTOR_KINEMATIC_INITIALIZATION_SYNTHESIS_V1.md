# First kinematic initialization synthesis for DADA motor optimization

## Status and purpose

This note defines a **first practical initialization law** for future DADA motor optimizations.

It combines two complementary bodies of numerical evidence:

1. the multi-temperature and passive-valve campaign documented in
   [`MOTOR_TEMPERATURE_AND_VALVE_CAMPAIGN.md`](MOTOR_TEMPERATURE_AND_VALVE_CAMPAIGN.md);
2. the later 260 K search for a smooth compact C2 motion law, culminating in the
   15-parameter hybrid family and the retained near-optimal candidate **3952**.

The objective is not to declare candidate 3952 to be a universal optimum. The objective is to
use it as a **high-quality local prior** so that future searches start close to the physically
relevant region instead of rediscovering the same structure from a broad parameter space.

This document therefore proposes a **continuation strategy**:

> predict the main thermodynamic and kinematic trends from the operating conditions, deform
> the 3952 reference motion accordingly, optimize only the parameters known to be sensitive,
> and release the remaining degrees of freedom only if the reduced search stops improving.

This is version 1. It currently applies to **motor operation with air** in the temperature range
already explored. It should later be extended with helium data and then tested in refrigeration
mode.

---

## 1. Evidence used by this synthesis

### 1.1 Temperature and valve campaign

The temperature campaign established several trends that are sufficiently regular to be useful
as initialization laws.

For the preferred **UU passive-valve topology**, with a cold source fixed at 298.15 K and
temperature differences from 160 to 360 K:

| ΔT (K) | S/L swept ratio | LP exchange (deg) | Compression (deg) | HP exchange (deg) | Expansion (deg) |
|---:|---:|---:|---:|---:|---:|
| 160 | 0.9584 | 180.6 | 36.3 | 117.4 | 25.7 |
| 200 | 0.9050 | 178.0 | 46.4 | 101.2 | 34.4 |
| 240 | 0.8763 | 179.5 | 38.8 | 105.6 | 36.1 |
| 260 | 0.8614 | 177.9 | 54.8 | 92.8 | 34.4 |
| 270 | 0.8528 | 173.0 | 50.8 | 100.6 | 35.6 |
| 300 | 0.8351 | 174.7 | 47.9 | 104.7 | 32.7 |
| 330 | 0.8285 | 169.1 | 50.2 | 104.6 | 36.0 |
| 360 | 0.8242 | 169.1 | 52.6 | 100.6 | 37.8 |

The strongest trends are:

- the optimum **small/large swept-volume ratio decreases** as ΔT increases;
- **LP exchange becomes shorter**;
- **compression becomes longer overall**;
- **expansion tends to become longer**;
- **HP exchange tends to become shorter**, although this trend contains more local numerical noise.

The four-stage campaign also found the following intermediate normalized piston levels:

| ΔT (K) | `a_l` | `b_l` | `a_s` | `b_s` |
|---:|---:|---:|---:|---:|
| 160 | 0.4467 | 0.0964 | 0.4213 | 0.1348 |
| 200 | 0.4098 | 0.1052 | 0.3916 | 0.1505 |
| 240 | 0.3789 | 0.1121 | 0.3703 | 0.1469 |
| 260 | 0.3723 | 0.0996 | 0.3609 | 0.1135 |
| 270 | 0.3666 | 0.1115 | 0.3590 | 0.1456 |
| 300 | 0.3590 | 0.1141 | 0.3504 | 0.1584 |
| 330 | 0.3618 | 0.1191 | 0.3486 | 0.1584 |
| 360 | 0.3482 | 0.1241 | 0.3295 | 0.1570 |

The decrease of `a_l` and `a_s` is clear and coherent. The variations of `b_l`, and especially
`b_s`, are weaker and noisier.

The same campaign showed an important thermodynamic consistency: the optimized cylinder ratio
and exchanged stroke fractions independently reproduce, within a few percent, the approximate
quasi-isobaric exchange relation

`dV_receiver / (-dV_donor) ≈ T_receiver / T_donor`.

This gives a physical reason to preserve the observed volume-ratio and exchange-stroke trends
rather than treating them as arbitrary optimizer output.

### 1.2 Smooth compact C2 motion search at ΔT = 260 K

The later smooth-motion campaign replaced the idealized four-stage law by progressively more
flexible C2 representations and finally by a structured **15-parameter hybrid C2 family**.

The search showed that a good smooth motion does not need the large number of free controls of
the original spline representation. The compact family can reproduce the important
thermodynamic features with a much smaller parameter set.

The retained reference candidate is **3952**.

The important observation is not only its efficiency. A family of very close runners-up reached
essentially the same thermodynamic performance, and their velocity and acceleration curves near
the HP transition had the **same qualitative shape**, with only modest amplitude differences.

This suggests that the pronounced HP-side feature is not an isolated numerical accident of one
candidate. It appears to belong to the local high-performance motion family.

Candidate 3952 should therefore be regarded as a **representative of a near-optimal basin**, not
as a unique mathematical point.

---

## 2. What should be transferred from candidate 3952

The 15 numerical coordinates of the hybrid parameterization are convenient for the optimizer,
but they are not the most robust description to transfer between machine configurations.

The transferable object should instead be described in **physical kinematic space**.

The reference motion should preserve, as far as possible:

1. the angular locations of the small- and large-cylinder extrema;
2. the durations of the four thermodynamic portions of the cycle;
3. the finite C2 curvature at each piston extremum;
4. the broad shape of the LP-side branches;
5. the location, sign and relative strength of the HP-side curvature concentration;
6. the relative phasing between the two pistons;
7. the fact that each piston has one maximum and one minimum per cycle, without parasitic
   oscillations.

This distinction is important. Future motion families may use different internal parameters.
The **physical shape** is the prior; the current 15-parameter coordinates are only one encoding.

---

## 3. Features of the 260 K compact reference that should initially be preserved

The final 15-parameter search produced several useful qualitative conclusions.

### 3.1 LP-side motion

The LP-side bowing parameters settled close to the middle of their admissible branches.

For the best converged compact solution, representative values were approximately:

- `small_up_bp_mid_q ≈ 0.478`;
- `large_down_bp_mid_q ≈ 0.512`.

There is therefore no evidence that the LP portions must be forced to be exactly linear.
However, there is also no evidence that they require complicated oscillatory structure.

The first initialization law should consequently use a **simple, smooth, weakly bowed LP
branch**, close to the 3952 shape, rather than reopen a high-dimensional LP search.

### 3.2 HP-side motion

The HP-side branches contain the strongest localized change in curvature.

This feature is:

- finite;
- C2-continuous;
- relatively sharp;
- repeated across the near-optimal candidates;
- not eliminated by comparing efficiency/power runners-up.

It should therefore be preserved in the initial seed.

The mechanical mechanism search may later soften it naturally. The thermodynamic optimizer
should not artificially remove it before the realizable mechanism is studied.

### 3.3 Extremum curvature

The optimum did not drive the extremum curvatures to zero or to the allowed bounds.

The initial seed should therefore retain **finite nonzero curvature at piston extrema**.

A sinusoidal extremum is not assumed, and a flat plateau is not imposed unless a future campaign
provides evidence for it.

---

## 4. First temperature-dependent initialization law

For version 1, the most reliable independent variable is the source temperature difference

`ΔT = T_hot_source - T_cold_source`

with the important restriction that the present evidence comes from a campaign with the cold
source fixed at 298.15 K.

Within the explored 160–360 K interval, the following quantities should be initialized by
interpolation of the existing UU campaign rather than reset to generic defaults.

### 4.1 High-confidence temperature-dependent quantities

These should remain free in future optimization and should receive a ΔT-dependent initial value:

- small/large swept-volume ratio;
- LP exchange duration;
- compression duration;
- expansion duration;
- hot-side exchanger size;
- `a_l`;
- `a_s`.

The expected directions with increasing ΔT are:

| Quantity | Expected trend as ΔT increases | Confidence |
|---|---|---|
| `V_S / V_L` | decreases | high |
| LP exchange duration | decreases | high |
| compression duration | increases overall | medium-high |
| expansion duration | increases overall | medium |
| `a_l` | decreases | high |
| `a_s` | decreases | high |
| hot-side exchanger area | decreases strongly | high |

### 4.2 Medium-confidence quantities

These should remain adjustable but should be searched locally around the interpolated/reference
value:

| Quantity | Expected behavior | Confidence |
|---|---|---|
| HP exchange duration | tends to decrease, but with local noise | medium |
| cold-side exchanger area | comparatively weak variation | medium |
| `b_l` | weak upward tendency | low-medium |
| extremum timing offsets within the smooth 15p law | likely temperature-dependent | medium |

### 4.3 Low-confidence quantities

These should initially stay close to the 3952 reference and be released only after the reduced
optimization has converged:

- LP bowing parameters;
- detailed HP kink position within each branch;
- detailed HP kink width;
- detailed HP kink amplitude/shape;
- individual extremum curvature magnitudes, except where required to follow a changed branch
  duration;
- other shape parameters for which the temperature campaign provides no independent trend.

This is a deliberate reduction of search dimensionality.

---

## 5. Bridging the old four-stage law and the 15-parameter C2 law

The temperature campaign and the compact C2 campaign use different motion parameterizations.

They should **not** be connected by attempting to identify old optimization coordinates with
new coordinates one by one.

Instead, both laws should be converted to a common set of physical descriptors.

A useful descriptor vector is:

`K = {`
`extremum angles,`
`LP duration,`
`compression duration,`
`HP duration,`
`expansion duration,`
`exchange stroke fractions,`
`extremum curvatures,`
`HP curvature-feature location/width/strength`
`}`

The temperature campaign predicts the components of `K` that have demonstrated thermal trends.
Candidate 3952 supplies the remaining local shape information.

A fitting step then converts the target descriptor vector back into the 15-parameter hybrid C2
coordinates.

This allows the initialization model to survive future changes of kinematic representation.

---

## 6. Recommended optimization procedure for a new air-motor configuration

### Stage 0 — choose the nearest validated reference

Use the existing library of converged champions.

For the present database:

- use candidate 3952 as the primary smooth-motion reference around ΔT = 260 K;
- for a different ΔT, use the closest existing temperature point as the thermodynamic reference;
- when several optimized smooth champions become available, interpolate between them rather than
  extrapolating from 260 K whenever possible.

### Stage 1 — predict the thermodynamic geometry

From the temperature campaign, initialize:

- cylinder swept-volume ratio;
- hot and cold exchanger geometry;
- the main phase durations or their smooth-motion equivalents.

Simple piecewise-linear interpolation is sufficient for version 1.

There is no present justification for fitting a high-order polynomial through the sparse data,
especially for the noisier phase durations.

Outside 160–360 K, use extrapolation only as a starting guess and enlarge the search radius.

### Stage 2 — deform the 3952 motion in physical space

Start from the normalized 3952 motion and modify only the physical descriptors for which the
temperature campaign provides evidence.

In particular:

- shorten the LP portion as ΔT rises;
- lengthen compression overall;
- lengthen expansion overall;
- allow HP duration to follow its weaker decreasing trend;
- modify the exchanged stroke fractions consistently with the interpolated `a_l`, `a_s`,
  cylinder ratio, and the quasi-isobaric exchange relation;
- preserve the basic HP curvature feature unless the optimizer demonstrates a benefit from
  changing it.

The deformation should remain C2.

### Stage 3 — reduced local optimization

The first numerical optimization should vary only:

- the thermodynamic geometry known to depend on operating point;
- global/branch timing variables;
- the strongest temperature-sensitive exchange-stroke descriptors.

Keep the detailed 3952 shape parameters fixed or tightly bounded.

This stage is intended to move rapidly onto the local optimum manifold.

### Stage 4 — release medium-confidence parameters

If progress stalls, release:

- HP duration;
- extremum curvature magnitudes;
- `b_l` / `b_s` equivalents;
- moderate local motion-shape adjustments.

The optimizer should still remain centered on the continuation prediction.

### Stage 5 — release the full 15-parameter motion only when justified

A full local 15p refinement should be the final step, not the starting point.

The existing 15p search already provides a useful hierarchy of search radii. For a continuation
run close to a known operating point, start from the fine/intermediate levels rather than from
the widest global probes.

If the optimum consistently migrates toward the edge of the local trust region, enlarge the
search radius by one level and continue.

A full wide search should be reserved for cases where:

- the operating point lies far outside the validated temperature interval;
- a new gas is introduced;
- the valve topology or machine architecture changes;
- the reduced continuation solution loses significant performance;
- the optimum repeatedly reaches imposed local bounds.

---

## 7. Why continuation should be more efficient than repeated global searches

The current results suggest that the optimum is not an isolated point in a featureless
15-dimensional space.

Instead, there appears to be a relatively smooth **manifold of good machine configurations**:

- machine proportions evolve systematically with source temperature;
- several motion descriptors evolve systematically;
- the compact smooth optimum lies in a basin containing many near-equivalent candidates;
- the HP-side acceleration structure remains qualitatively stable across those candidates.

A future optimizer should therefore search primarily **along this manifold**, not repeatedly
rediscover it.

In practical terms, this should reduce:

- the number of broad Sobol/global proposals;
- the number of structurally poor motion candidates;
- the time spent rediscovering the same HP feature;
- sensitivity to arbitrary initialization;
- the number of human interventions needed before reaching the fine-refinement stage.

---

## 8. Proposed machine-readable seed model

The Markdown synthesis should eventually be backed by a small machine-readable database.

A possible structure is:

```text
working_gas
operating_mode
T_cold_source
T_hot_source
delta_T
valve_topology
machine_geometry
exchanger_geometry
physical_kinematic_descriptors
motion_family
motion_parameters
efficiency
power
validity_metrics
candidate_id
campaign_id
```

The important design choice is to store **both**:

1. the native optimizer coordinates;
2. the representation-independent physical descriptors.

This makes future interpolation and comparison possible even if the motion law changes.

---

## 9. Extension to helium

The next important test is to repeat selected operating points with helium.

The helium campaign should answer two different questions.

### 9.1 Which trends are truly thermodynamic?

If the same normalized motion trends persist with helium, they are likely linked primarily to:

- source/gas temperature ratios;
- compression and expansion thermodynamics;
- exchange-volume matching.

Those trends may eventually be expressible with dimensionless thermal variables rather than by
gas identity.

### 9.2 Which trends are gas-dependent?

Helium changes:

- heat capacity;
- heat-capacity ratio;
- thermal conductivity;
- viscosity;
- density for a given pressure and temperature;
- speed of sound;
- exchanger Reynolds/Mach behavior.

The seed law should therefore not assume in advance that `ΔT` alone remains sufficient.

Until helium data demonstrate otherwise, air and helium should use **separate local seed
families**, linked through common physical descriptors.

After enough data exist, the independent variables may be upgraded from

`K* = f(ΔT)`

to a more general relation such as

`K* = f(T_cold, T_hot, gas properties, machine scale, operating speed, ...)`.

---

## 10. Extension to refrigeration mode

Refrigeration mode should be treated as a new optimization branch rather than obtained by
blindly reversing the motor law.

The same framework remains useful:

- physical kinematic descriptors;
- continuation from nearby validated points;
- reduced optimization first;
- full motion release only if necessary.

The refrigeration campaign should determine:

- which motor trends reverse;
- which remain invariant;
- whether the same HP-side curvature feature remains beneficial;
- whether cylinder-ratio scaling follows the same exchange-volume/temperature relation;
- whether the preferred passive-valve topology remains unchanged;
- whether the best motion belongs to the same compact C2 family.

Until this is demonstrated, motor and refrigerator seeds should remain distinct.

---

## 11. Version-1 initialization rules

For a new **air / motor / UU** optimization inside or near the explored temperature range:

1. **Do not start from a generic motion.**
   Start from candidate 3952 or from the nearest later validated smooth champion.

2. **Interpolate the thermodynamic geometry from the temperature campaign.**
   In particular, initialize cylinder ratio and exchanger geometry from the existing map.

3. **Predict the main phase timing trends from ΔT.**
   LP shorter with increasing ΔT; compression and expansion generally longer; HP somewhat
   shorter but with lower confidence.

4. **Use `a_l` and `a_s` trends as strong constraints on exchanged stroke fractions.**
   Preserve consistency with the quasi-isobaric exchange-volume relation.

5. **Preserve the 3952 detailed branch shape by default.**
   Especially preserve the qualitative HP-side curvature/acceleration feature.

6. **Keep weakly supported shape parameters fixed during the first optimization stage.**
   Do not pay again for a broad search unless the local continuation fails.

7. **Release parameters progressively.**
   Geometry and timings first; medium-confidence motion parameters second; all 15p last.

8. **Use the optimizer as a corrector, not as a rediscovery engine.**
   The prediction supplies the neighborhood; the optimizer determines the precise local optimum.

9. **Store every new converged champion in the seed library.**
   The initialization model should improve after every campaign.

10. **Do not interpret version 1 as a universal law.**
    It is an empirical continuation model for the current solver, air, motor operation, and the
    present validated temperature range.

---

## 12. Main hypothesis to test in future campaigns

The working hypothesis behind this synthesis is:

> The optimum DADA motor does not require an entirely new piston motion for every operating
> condition. It belongs to a compact family whose main thermodynamic proportions and phase
> timings evolve smoothly with operating conditions, while much of the detailed normalized C2
> branch shape remains reusable.

If this hypothesis is confirmed by helium and refrigeration campaigns, future optimization can
move from repeated high-dimensional searches toward a **predictor-corrector design method**:

1. predict the optimum from a compact physical model;
2. correct it with a short local solver campaign;
3. add the new optimum to the database;
4. improve the predictor.

That would turn the present numerical optimization history into a progressively more general
design law rather than a collection of isolated champions.
