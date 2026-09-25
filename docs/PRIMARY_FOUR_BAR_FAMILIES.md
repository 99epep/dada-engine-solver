# Primary four-bar families from the candidate-3952 synthesis

## Purpose

This note is a compact catalogue of the four primary four-bar mechanisms that were deliberately carried into the downstream six-bar synthesis for structured C2 candidate 3952.

It is not a ranking of all four-bar mechanisms found by the saturation campaign. The historical identifiers **1, 4, 12 and 50** came from the V4.1 design-library ranking used when the downstream panel was assembled. They were retained because they represented useful mechanical diversity and different downstream opportunities.

The complete search method and the statistical saturation analysis are documented in `MECHANISM_SYNTHESIS_SEARCH_THEORY.md`.

---

## 1. Geometry convention

All dimensions are normalized by the crank radius:

\[
AB = 1.
\]

The primary mechanism is the four-bar loop `A-B-C-D`.

- `A = (0, 0)` is the crank axis.
- `D = (AD, 0)` is the fixed rocker pivot.
- `AB = 1` is the crank.
- `BC` is the coupler.
- `CD` is the rocker.
- the crank angle used by the mechanism is the study angle plus `primary_phase`;
- the assembly branch is `+1` or `-1`.

Point `E` is rigidly attached to the coupler. If `B` and `C` are the coupler endpoints, the solver uses

\[
E =
B
+
rac{E_{along}}{BC}(C-B)
+
rac{E_{normal}}{BC}R_{90}(C-B),
\]

where `R90` rotates the coupler vector by +90°.

Therefore `E_along` and `E_normal` are absolute lengths in crank-radius units, not fractions of `BC`.

A physical mechanism with crank radius `R` is obtained by multiplying every length in this document by `R`. Angles and assembly branches are unchanged.

---

## 2. Reference cadence

The primary search was built around the LARGE piston of candidate 3952 because it carried the more restrictive cadence change.

The important target features were:

- high-position turnaround: approximately **358.812°**;
- low-position turnaround: approximately **154.144°**;
- short high-to-low branch: approximately **155.332°**;
- long low-to-high branch: approximately **204.668°**;
- fast/slow speed ratio on the short branch: approximately **5.027**;
- fast-branch displacement fraction: approximately **57.63%**.

The long low-pressure branch was not forced to reproduce the detailed candidate-3952 velocity profile. Its detailed waviness was treated as a secondary proxy issue rather than a physical requirement.

---

## 3. Saturation context

A dedicated V4.1 saturation campaign later ran **512 fresh differential-evolution islands** over the full six-dimensional primary search box.

At family-distance threshold `0.04`:

- 16 families were observed;
- the largest family captured **47.66%** of the islands;
- the second captured **44.53%**;
- the first four together captured **96.68%**;
- the Good-Turing unseen capture-mass estimate was **1.56%**.

At threshold `0.05`, the first four families captured **98.05%** and the unseen-mass estimate fell to approximately **1.17%**.

The saturation result therefore shows two extremely dominant basins followed by a small tail of rarer families.

The four historical representatives in this document must **not** be interpreted as a one-to-one list of the four largest saturation clusters. Their purpose was different: preserve mechanically different candidates long enough to test downstream transformability.

---

## 4. Four retained primary representatives

| Historical ID | Working description | Branch | AD/AB | BC/AB | CD/AB | E along | E normal | Phase | min |sin μ| | V4.1 score | BE/BC |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Compact balanced | -1 | 1.568637 | 1.614189 | 1.413477 | 1.903279 | 0.415508 | -21.011° | 0.346720 | 0.325752 | 1.207 |
| 4 | Compact high-stroke | -1 | 1.748507 | 2.232955 | 1.947585 | 2.292977 | 1.392213 | -18.033° | 0.327220 | 0.474417 | 1.201 |
| 12 | Long-ground alternate | +1 | 3.632829 | 2.447078 | 2.260630 | 4.570318 | 1.339373 | 178.363° | 0.349942 | 1.176907 | 1.946 |
| 50 | Long-link robust | +1 | 5.387160 | 5.295486 | 1.212507 | 4.146499 | -3.932955 | -175.151° | 0.477448 | 2.188842 | 1.079 |

The scalar V4.1 score is useful only inside the primary-cadence search. It is not a prediction of final six-bar efficiency.

---

## 5. Family 1 — compact balanced primary

Family 1 was the most productive starting point of the downstream panel.

Normalized geometry:

```text
AB = 1
AD = 1.568637064
BC = 1.614188554
CD = 1.413477246

E_along  = 1.903278626
E_normal = 0.415508230

primary_phase  = -0.366712672 rad  = -21.011 deg
assembly_branch = -1
```

Key descriptors:

- V4.1 cadence score: **0.325752**;
- minimum primary transmission sine: **0.346720**;
- `BE / BC`: approximately **1.207**.

This is a relatively compact geometry with all three primary moving-link lengths of the same order as the crank. It became the source of the best final six-bar family.

The important lesson is not that this primary had the absolute best primary score. Its importance is that it combined a good cadence proxy with very strong **downstream transformability**.

---

## 6. Family 4 — compact high-stroke primary

Normalized geometry:

```text
AB = 1
AD = 1.748506860
BC = 2.232955203
CD = 1.947584850

E_along  = 2.292976503
E_normal = 1.392213452

primary_phase  = -0.314732756 rad  = -18.033 deg
assembly_branch = -1
```

Key descriptors:

- V4.1 cadence score: **0.474417**;
- minimum primary transmission sine: **0.327220**;
- `BE / BC`: approximately **1.202**.

Compared with family 1, this primary uses a larger coupler and rocker and a more strongly offset coupler point `E`.

Its descendants remained compact but eventually pushed close to the retained stroke and transmission limits. This makes it useful as an example of a family whose nominal geometry is attractive but whose downstream optimum becomes mechanically constrained.

---

## 7. Family 12 — long-ground alternate primary

Normalized geometry:

```text
AB = 1
AD = 3.632828600
BC = 2.447077658
CD = 2.260629683

E_along  = 4.570318453
E_normal = 1.339372687

primary_phase  = 3.113027812 rad  = 178.363 deg
assembly_branch = +1
```

Key descriptors:

- V4.1 cadence score: **1.176907**;
- minimum primary transmission sine: **0.349942**;
- `BE / BC`: approximately **1.946**.

This family is geometrically much farther from the compact rank-1/rank-4 region. Its crank phase is almost opposite, the fixed ground is much longer, and point `E` lies far beyond the nominal coupler length.

Its primary score is much worse than families 1 and 4, yet it was deliberately retained because the downstream dyad can transform a primary motion in ways that the simplified cadence score does not predict.

It later became the clearest demonstration that thermodynamic adaptation can move a mechanism substantially away from the original target motion while still increasing machine performance.

---

## 8. Family 50 — long-link robust historical primary

Normalized geometry:

```text
AB = 1
AD = 5.387159727
BC = 5.295485951
CD = 1.212507384

E_along  = 4.146499490
E_normal = -3.932955365

primary_phase  = -3.056968382 rad  = -175.151 deg
assembly_branch = +1
```

Key descriptors:

- V4.1 cadence score: **2.188842**;
- minimum primary transmission sine: **0.477448**;
- `BE / BC`: approximately **1.079**.

This family came from the older successful six-bar solution and was kept as a control.

Its primary cadence score is poor compared with the newer V4.1 winners, but its transmission margin is much larger and its descendants remained mechanically robust.

The family was therefore essential evidence that:

> a primary that looks mediocre under a reduced cadence proxy can still be a strong carrier for a complete six-bar mechanism.

Its final six-bar descendant became the second-best thermodynamic family after re-tuning.

---

## 9. How these primaries should be used

These four mechanisms are best treated as **seed families**, not immutable designs.

For a new target:

- normalize the crank to `AB = 1`;
- use the stored geometry as a warm start;
- allow moderate deformation of all primary dimensions;
- retain a fresh-island component to detect topology changes;
- do not assume the historically best family will remain best.

The appropriate starting side is the piston with the more restrictive kinematic demand. In the 3952 motor this was the LARGE piston. In another machine, especially a refrigeration cycle, it may be the SMALL piston.

---

## 10. What should not be inferred

This catalogue does not establish that:

- only four primary families exist;
- these four are the four largest saturation basins;
- lower V4.1 score guarantees a better six-bar;
- the primary geometry alone predicts thermodynamic efficiency;
- any one family is globally optimal.

The completed saturation campaign supports the narrower statement that the dominant primary basins are well mapped under the declared V4.1 search policy.

The downstream campaign supports a different statement: the four representatives above span useful differences in transformability and mechanical character.

---

## 11. Source artifacts

Primary panel and cadence data:

```text
outputs/large_primary_cadence_3952_v4_1.json
outputs/large_second_dyad_v41_panel_3952.json
```

Saturation:

```text
examples/search_large_primary_cadence_3952_v4_1_saturation.py
outputs/large_primary_cadence_3952_v4_1_saturation.json
```

Complete methodology:

```text
docs/MECHANISM_SYNTHESIS_SEARCH_THEORY.md
```

The final paired six-bar descendants are catalogued separately in:

```text
docs/SIX_BAR_MECHANISM_FAMILIES.md
```
