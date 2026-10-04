# Primary four-bar seed families

## Purpose

This note catalogues four primary four-bar geometries retained as mechanically
distinct seeds while converting a structured-C2 thermodynamic reference motion
at a 260 K source-temperature difference into realizable linkages.

These mechanisms are **seed families**, not a ranking and not an exhaustive
classification of the primary four-bar search space.

Their identifiers are descriptive:

- `compact_balanced`;
- `compact_offset`;
- `long_ground`;
- `long_coupler_short_rocker`.

The search method, cadence objective and basin-discovery methodology belong in
[`MECHANISM_SYNTHESIS_SEARCH.md`](MECHANISM_SYNTHESIS_SEARCH.md).
This catalogue intentionally does not duplicate saturation statistics or
historical optimizer rankings.

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
- the crank angle used by the mechanism is the study angle plus `phase_rad`;
- the assembly branch is `+1` or `-1`.

Point `E` is rigidly attached to the coupler. If `B` and `C` are the coupler
endpoints,

\[
E =
B
+
\frac{E_{\parallel}}{BC}(C-B)
+
\frac{E_{\perp}}{BC}R_{90}(C-B),
\]

where \(R_{90}\) rotates the coupler vector by +90°.

Therefore `e_along` and `e_normal` are lengths in crank-radius units, not
fractions of `BC`.

A physical mechanism with crank radius \(R\) is obtained by multiplying every
length by \(R\). Angles and assembly branches are unchanged.

The primary transmission indicator used here is

\[
s_{\min}
=
\min_\theta
\frac{|BC \times DC|}{|BC||DC|}
=
\min_\theta |\sin\mu|.
\]

It is a geometric indicator of proximity to four-bar toggle conditions. It is
not by itself a complete manufacturing-robustness metric.

---

## 2. Retained seed geometries

The exact coordinates are stored in
[`repro/primary_four_bar_families/families.toml`](repro/primary_four_bar_families/families.toml).

| Seed | Branch | AD/AB | BC/AB | CD/AB | E along | E normal | Phase | min \|sin μ\| | BE/BC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `compact_balanced` | -1 | 1.568637 | 1.614189 | 1.413477 | 1.903279 | 0.415508 | -21.011° | 0.346720 | 1.206864 |
| `compact_offset` | -1 | 1.748507 | 2.232955 | 1.947585 | 2.292977 | 1.392213 | -18.033° | 0.327220 | 1.201339 |
| `long_ground` | +1 | 3.632829 | 2.447078 | 2.260630 | 4.570318 | 1.339373 | 178.363° | 0.349907 | 1.946213 |
| `long_coupler_short_rocker` | +1 | 5.387160 | 5.295486 | 1.212507 | 4.146499 | -3.932955 | -175.151° | 0.477443 | 1.079227 |

The tabulated values are rounded for readability. Use the TOML data for
calculation.

---

## 3. `compact_balanced`

This seed keeps the ground, coupler and rocker lengths of comparable magnitude.

Point `E` extends moderately beyond the coupler endpoint and has a relatively
small normal offset.

Its value during the synthesis campaign was not that it represented a globally
best four-bar mechanism. It provided a compact primary motion that proved
highly transformable by the downstream dyad.

This is the natural first seed when a compact primary mechanism is desired and
there is no evidence that a more unusual topology is required.

---

## 4. `compact_offset`

This seed remains compact but uses a larger coupler and rocker than
`compact_balanced`.

Point `E` also has a much larger normal offset.

The geometry is useful because it represents a different way of obtaining a
similar temporal cadence from a compact mechanism. Its downstream descendants
demonstrated that an attractive primary geometry can still approach mechanical
limits once the complete mechanism is optimized.

It should therefore remain distinct from `compact_balanced` rather than being
treated as a minor perturbation of it.

---

## 5. `long_ground`

This seed occupies a substantially different region of four-bar geometry.

Its fixed-ground distance is much larger than the crank radius, its phase is
nearly opposite the compact seeds, and point `E` lies far beyond the nominal
coupler length.

It was retained because the downstream dyad can transform primary motion in
ways that a reduced primary-only cadence metric does not predict.

The broader lesson is important:

> A primary four-bar is an intermediate motion generator. Its usefulness must
> ultimately be judged after downstream mechanical and thermodynamic
> adaptation.

---

## 6. `long_coupler_short_rocker`

This seed combines a long fixed ground and coupler with a comparatively short
rocker.

Point `E` is strongly offset from the coupler axis.

Among these four seeds it has the largest minimum primary transmission sine,
but that fact alone must not be interpreted as a complete robustness ranking.

Its downstream descendants provided an important control against an
overly narrow preference for compact primary geometries.

---

## 7. How these seeds should be used

These geometries are starting families, not immutable designs.

For a new target:

- normalize the crank to `AB = 1`;
- use one or more stored geometries as warm starts;
- allow moderate deformation of the primary dimensions;
- retain fresh global exploration when the target differs substantially from
  the reference motion;
- do not assume that the seed with the most attractive primary-only metric will
  produce the best complete mechanism.

The appropriate starting piston is the one carrying the more restrictive
kinematic demand. That choice depends on the thermodynamic cycle and is not a
permanent property of these four-bar geometries.

---

## 8. Primary geometry versus complete-mechanism robustness

Manufacturing tolerance and final mechanical robustness cannot be assigned from
the primary four-bar alone.

A complete six-bar descendant adds:

- a downstream dyad;
- a second assembly branch;
- a transformed output point;
- finite piston-rod geometry;
- slider alignment;
- thermodynamic adaptation.

Those additions can amplify or attenuate sensitivity inherited from the
primary.

Consequently this document records only primary geometry and its primary
transmission indicator.

Complete-mechanism comparisons belong in
[`SIX_BAR_MECHANISM_FAMILIES.md`](SIX_BAR_MECHANISM_FAMILIES.md).

---

## 9. What should not be inferred

This catalogue does not establish that:

- only four useful primary geometries exist;
- these seeds correspond to the largest optimizer basins;
- one of them is globally optimal;
- a larger primary transmission sine guarantees a better complete mechanism;
- a compact primary necessarily produces a compact or robust final mechanism;
- primary geometry alone predicts thermodynamic efficiency.

The four seeds are retained because they span useful differences in geometry
and downstream transformability.

---

## 10. Reproducing this catalogue

The reproducible source data are:

```text
docs/repro/primary_four_bar_families/families.toml
docs/repro/primary_four_bar_families/verify.py
```

Run:

```bash
python docs/repro/primary_four_bar_families/verify.py
```

or, for a regression-style check:

```bash
python docs/repro/primary_four_bar_families/verify.py --check
```

The script:

1. loads the exact normalized geometries;
2. verifies full-revolution four-bar closure on a dense periodic scan;
3. evaluates the selected assembly branch;
4. finds the minimum primary transmission sine over the cycle;
5. computes `BE/BC`;
6. checks the reproduced values against the small reference values stored in
   the TOML.
