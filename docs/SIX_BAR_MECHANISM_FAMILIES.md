# Paired six-bar mechanism families

## Purpose

This document catalogues four mechanically distinct paired SMALL/LARGE six-bar
mechanism lineages for use as synthesis seeds.

The catalogue specifies normalized geometry and nominal mechanical margins,
independently of a thermodynamic operating point.

Each family contains two independently adapted mechanisms:

- one for the SMALL cylinder;
- one for the LARGE cylinder.

The exact normalized coordinates are stored in
[`repro/six_bar_mechanism_families/families.toml`](repro/six_bar_mechanism_families/families.toml).

The four retained lineages are:

- `compact_balanced`;
- `compact_offset`;
- `long_ground`;
- `long_coupler_short_rocker`.

These identifiers describe the primary lineage and mechanical character.
They do not imply that the final primary four-bar geometry is unchanged from
the corresponding seed in
[`PRIMARY_FOUR_BAR_FAMILIES.md`](PRIMARY_FOUR_BAR_FAMILIES.md).

The mechanisms were allowed to deform during complete six-bar adaptation.

Search methodology belongs in
[`MECHANISM_SYNTHESIS_SEARCH.md`](MECHANISM_SYNTHESIS_SEARCH.md).

---

## 1. Topology and coordinate convention

Each piston is driven by:

```text
primary four-bar:  A-B-C-D
coupler point:     E fixed on BC
secondary dyad:    E-F-G
output point:      H fixed on EF
piston rod:        H-P
slider:           P constrained to a straight axis
```

All absolute lengths are expressed in crank-radius units:

\[
AB=1.
\]

The fixed primary pivots are

\[
A=(0,0),
\qquad
D=(AD,0).
\]

The crank point is

\[
B=
(\cos(\theta+\phi),\sin(\theta+\phi)).
\]

Point `C` follows from primary four-bar closure using the stored assembly
branch.

Point `E` is rigidly attached to the coupler:

\[
E=
B+
\frac{E_{\parallel}}{BC}(C-B)
+
\frac{E_{\perp}}{BC}R_{90}(C-B).
\]

The secondary fixed pivot is

\[
G=(G_x,G_y).
\]

Point `F` follows from the `E-F-G` RR closure using the stored secondary
branch.

Point `H` is rigidly attached to the `EF` frame:

\[
H=
E+
h_{\parallel}(F-E)
+
h_{\perp}R_{90}(F-E),
\]

where `h_along_over_ef` and `h_normal_over_ef` are dimensionless fractions.

The piston rod joins `H` to the positive slider solution used by the solver.

The slider axis is defined by:

- `slider_axis_angle`;
- `slider_axis_offset`.

The production implementation is
`dada_solver.six_bar.SixBarCylinderMechanism`.

---

## 2. Scaling to a physical mechanism

The catalogue geometry is dimensionless.

For a chosen physical crank radius \(R\), multiply by \(R\):

- `AD`;
- `BC`;
- `CD`;
- `E_along`;
- `E_normal`;
- `G_x`;
- `G_y`;
- `EF`;
- `GF`;
- piston-rod length;
- slider-axis offset.

Do not scale:

- `H_along / EF`;
- `H_normal / EF`;
- phases;
- slider-axis angles;
- assembly branches.

Physical piston stroke is

\[
stroke
=
(stroke/crank)\,R.
\]

Cylinder bore is then determined from the required swept volume and physical
stroke.

The normalized linkage does not itself determine cylinder bore or machine
power.

---

## 3. Mechanical diagnostics

The catalogue uses the production diagnostic implementation:

```text
dada_solver.mechanism_diagnostics.six_bar_metrics
```

The reported full-cycle indicators are:

- `stroke_over_crank`;
- `minimum_primary_transmission_sine`;
- `minimum_secondary_transmission_sine`;
- `minimum_rod_axis_cosine`;
- `EH_over_crank`;
- `H_axis_lateral_rms_over_stroke`;
- `H_axis_lateral_span_over_stroke`;
- `crank_axis_to_EFH_clearance_over_crank`;
- `zero_crossing_count`.

The documentary verifier evaluates the standard 1440-sample screen and repeats
it on a denser 5760-sample grid.

This is a deterministic numerical screen of the nominal geometry.

It is not a certified continuous-angle proof and not a manufacturing-tolerance
probability.

Stroke uses the production mechanism's velocity-root extrema, bracketed on its
1440-interval construction scan; increasing the diagnostic grid does not change
that root search. `EH_over_crank` is a direct geometric length. The other
full-cycle indicators use the stated diagnostic sampling grid.

---

## 4. Retained design screen

The current catalogue screen is:

| Metric | Requirement |
|---|---:|
| `stroke_over_crank` | 1.0 .. 3.0 |
| `minimum_primary_transmission_sine` | ≥ 0.30 |
| `minimum_secondary_transmission_sine` | ≥ 0.30 |
| `minimum_rod_axis_cosine` | ≥ 0.95 |
| `EH_over_crank` | ≤ 7.0 |
| `crank_axis_to_EFH_clearance_over_crank` | ≥ 0.5 |
| `H_axis_lateral_rms_over_stroke` | ≤ 0.25 |
| `H_axis_lateral_span_over_stroke` | ≤ 0.65 |
| `zero_crossing_count` | = 2 |

These are retained design constraints, not literal material-failure thresholds.

In particular, stroke windows and H-lateral limits include packaging and design
judgement.

The machine-readable definition is stored in
[`repro/six_bar_mechanism_families/design_screen.toml`](repro/six_bar_mechanism_families/design_screen.toml).

All eight exact mechanisms pass all ten constraints at both resolutions. No
pass/fail result changes on the dense grid. Tables below show canonical metric
values and individual signed margins at both resolutions. Equality is reported
as pass/fail rather than as a continuous safety margin.

---

## 5. `compact_balanced`

This lineage retains relatively compact primary geometry with ground, coupler
and rocker lengths of similar order.

The SMALL and LARGE mechanisms are independent final mechanisms, not exact
mirrors.

The final LARGE geometry lies very close to the retained
`H_axis_lateral_span_over_stroke <= 0.65` limit.

That small nominal margin immediately shows that packaging or local mechanism
re-optimization should be revisited before freezing a physical design.

| Metric | SMALL | LARGE |
|---|---:|---:|
| `stroke_over_crank` | 2.540934 | 2.388543 |
| `minimum_primary_transmission_sine` | 0.366020 | 0.323713 |
| `minimum_secondary_transmission_sine` | 0.350129 | 0.361433 |
| `minimum_rod_axis_cosine` | 0.962550 | 0.962112 |
| `EH_over_crank` | 2.436754 | 2.438924 |
| `H_axis_lateral_rms_over_stroke` | 0.139669 | 0.201255 |
| `H_axis_lateral_span_over_stroke` | 0.472256 | 0.649979 |
| `crank_axis_to_EFH_clearance_over_crank` | 0.777609 | 0.804459 |
| `zero_crossing_count` | 2 | 2 |

Signed margins (canonical / dense); equality rows report pass/fail.

| Constraint | Unit | SMALL margin | LARGE margin |
|---|---|---:|---:|
| `stroke_over_crank` ≥ 1.0 | crank_radius | +1.54093433 / +1.54093433 | +1.38854271 / +1.38854271 |
| `stroke_over_crank` ≤ 3.0 | crank_radius | +0.45906567 / +0.45906567 | +0.611457286 / +0.611457286 |
| `minimum_primary_transmission_sine` ≥ 0.3 | 1 | +0.0660202985 / +0.0660178938 | +0.0237131698 / +0.0237123031 |
| `minimum_secondary_transmission_sine` ≥ 0.3 | 1 | +0.0501289688 / +0.0501282358 | +0.0614329434 / +0.061430797 |
| `minimum_rod_axis_cosine` ≥ 0.95 | 1 | +0.0125497585 / +0.0125497415 | +0.0121117859 / +0.0121117859 |
| `EH_over_crank` ≤ 7.0 | crank_radius | +4.56324552 / +4.56324552 | +4.56107582 / +4.56107582 |
| `crank_axis_to_EFH_clearance_over_crank` ≥ 0.5 | crank_radius | +0.277608863 / +0.27760861 | +0.304459114 / +0.304458452 |
| `H_axis_lateral_rms_over_stroke` ≤ 0.25 | 1 | +0.110330541 / +0.110330541 | +0.0487454353 / +0.0487454353 |
| `H_axis_lateral_span_over_stroke` ≤ 0.65 | 1 | +0.177743998 / +0.17774224 | +2.09422996e-05 / +1.34657081e-05 |
| `zero_crossing_count` = 2 | 1 | pass / pass | pass / pass |

---

## 6. `compact_offset`

This lineage remains compact but uses a larger and more strongly offset
primary geometry than `compact_balanced`.

It achieves a large stroke relative to crank radius.

Its main nominal weakness is clear from the production diagnostics: primary
and secondary transmission minima approach the retained `0.30` floor.

The family is therefore useful as a compact/high-stroke reference but should
not be interpreted as having large toggle margin.

| Metric | SMALL | LARGE |
|---|---:|---:|
| `stroke_over_crank` | 2.979288 | 2.962992 |
| `minimum_primary_transmission_sine` | 0.301339 | 0.301172 |
| `minimum_secondary_transmission_sine` | 0.300053 | 0.300339 |
| `minimum_rod_axis_cosine` | 0.959982 | 0.954711 |
| `EH_over_crank` | 2.869658 | 3.048954 |
| `H_axis_lateral_rms_over_stroke` | 0.167176 | 0.163360 |
| `H_axis_lateral_span_over_stroke` | 0.623755 | 0.637737 |
| `crank_axis_to_EFH_clearance_over_crank` | 1.589543 | 1.698068 |
| `zero_crossing_count` | 2 | 2 |

Signed margins (canonical / dense); equality rows report pass/fail.

| Constraint | Unit | SMALL margin | LARGE margin |
|---|---|---:|---:|
| `stroke_over_crank` ≥ 1.0 | crank_radius | +1.97928755 / +1.97928755 | +1.96299157 / +1.96299157 |
| `stroke_over_crank` ≤ 3.0 | crank_radius | +0.0207124537 / +0.0207124537 | +0.0370084324 / +0.0370084324 |
| `minimum_primary_transmission_sine` ≥ 0.3 | 1 | +0.00133854775 / +0.00133755901 | +0.00117206401 / +0.00117021208 |
| `minimum_secondary_transmission_sine` ≥ 0.3 | 1 | +5.32269021e-05 / +5.23019036e-05 | +0.000339283642 / +0.000337987292 |
| `minimum_rod_axis_cosine` ≥ 0.95 | 1 | +0.00998206916 / +0.00998200884 | +0.00471078097 / +0.00471078097 |
| `EH_over_crank` ≤ 7.0 | crank_radius | +4.13034229 / +4.13034229 | +3.95104632 / +3.95104632 |
| `crank_axis_to_EFH_clearance_over_crank` ≥ 0.5 | crank_radius | +1.08954259 / +1.08954259 | +1.19806785 / +1.19806776 |
| `H_axis_lateral_rms_over_stroke` ≤ 0.25 | 1 | +0.0828237161 / +0.0828237161 | +0.0866400476 / +0.0866400476 |
| `H_axis_lateral_span_over_stroke` ≤ 0.65 | 1 | +0.0262453272 / +0.0262405458 | +0.0122626274 / +0.0122626274 |
| `zero_crossing_count` = 2 | 1 | pass / pass | pass / pass |

---

## 7. `long_ground`

This lineage occupies a substantially different region of design space.

Its primary ground distance is larger and its downstream geometry contains long
links.

The final pair illustrates that mechanical constraints may become active in
different places on the two cylinder mechanisms:

- the SMALL mechanism approaches its rod-angle and H-lateral-span limits;
- the LARGE primary approaches the transmission-sine floor.

This family is therefore useful for separating the effect of primary topology
from downstream transformability.

| Metric | SMALL | LARGE |
|---|---:|---:|
| `stroke_over_crank` | 1.238708 | 1.900432 |
| `minimum_primary_transmission_sine` | 0.503546 | 0.300351 |
| `minimum_secondary_transmission_sine` | 0.530886 | 0.317379 |
| `minimum_rod_axis_cosine` | 0.950473 | 0.981588 |
| `EH_over_crank` | 4.090995 | 4.387379 |
| `H_axis_lateral_rms_over_stroke` | 0.227536 | 0.184521 |
| `H_axis_lateral_span_over_stroke` | 0.649800 | 0.549655 |
| `crank_axis_to_EFH_clearance_over_crank` | 3.065353 | 3.471589 |
| `zero_crossing_count` | 2 | 2 |

Signed margins (canonical / dense); equality rows report pass/fail.

| Constraint | Unit | SMALL margin | LARGE margin |
|---|---|---:|---:|
| `stroke_over_crank` ≥ 1.0 | crank_radius | +0.238708358 / +0.238708358 | +0.900431971 / +0.900431971 |
| `stroke_over_crank` ≤ 3.0 | crank_radius | +1.76129164 / +1.76129164 | +1.09956803 / +1.09956803 |
| `minimum_primary_transmission_sine` ≥ 0.3 | 1 | +0.203546456 / +0.203546097 | +0.000350926379 / +0.000350926379 |
| `minimum_secondary_transmission_sine` ≥ 0.3 | 1 | +0.230885521 / +0.230885521 | +0.0173793481 / +0.0173786398 |
| `minimum_rod_axis_cosine` ≥ 0.95 | 1 | +0.000473174532 / +0.000473174532 | +0.0315877154 / +0.0315877154 |
| `EH_over_crank` ≤ 7.0 | crank_radius | +2.90900458 / +2.90900458 | +2.61262088 / +2.61262088 |
| `crank_axis_to_EFH_clearance_over_crank` ≥ 0.5 | crank_radius | +2.56535316 / +2.56535316 | +2.97158922 / +2.97158869 |
| `H_axis_lateral_rms_over_stroke` ≤ 0.25 | 1 | +0.0224643057 / +0.0224643057 | +0.0654791609 / +0.0654791609 |
| `H_axis_lateral_span_over_stroke` ≤ 0.65 | 1 | +0.000199848105 / +0.000199589066 | +0.100345171 / +0.100343222 |
| `zero_crossing_count` = 2 | 1 | pass / pass | pass / pass |

---

## 8. `long_coupler_short_rocker`

This lineage combines:

- a long primary ground;
- a long coupler;
- a short rocker;
- large downstream link distances.

It retains comparatively large nominal transmission and crank-clearance
margins.

The SMALL mechanism instead lies relatively near the retained upper
`stroke_over_crank` limit.

The family is an important counterexample to the idea that compact primary
geometry is automatically preferable.

| Metric | SMALL | LARGE |
|---|---:|---:|
| `stroke_over_crank` | 2.973136 | 2.559605 |
| `minimum_primary_transmission_sine` | 0.463840 | 0.378746 |
| `minimum_secondary_transmission_sine` | 0.363427 | 0.374122 |
| `minimum_rod_axis_cosine` | 0.976569 | 0.958564 |
| `EH_over_crank` | 4.308265 | 3.309179 |
| `H_axis_lateral_rms_over_stroke` | 0.158978 | 0.131434 |
| `H_axis_lateral_span_over_stroke` | 0.559453 | 0.499738 |
| `crank_axis_to_EFH_clearance_over_crank` | 5.217632 | 4.811296 |
| `zero_crossing_count` | 2 | 2 |

Signed margins (canonical / dense); equality rows report pass/fail.

| Constraint | Unit | SMALL margin | LARGE margin |
|---|---|---:|---:|
| `stroke_over_crank` ≥ 1.0 | crank_radius | +1.97313568 / +1.97313568 | +1.55960466 / +1.55960466 |
| `stroke_over_crank` ≤ 3.0 | crank_radius | +0.026864317 / +0.026864317 | +0.440395339 / +0.440395339 |
| `minimum_primary_transmission_sine` ≥ 0.3 | 1 | +0.16383969 / +0.16383969 | +0.0787464849 / +0.0787442425 |
| `minimum_secondary_transmission_sine` ≥ 0.3 | 1 | +0.063426948 / +0.0634267648 | +0.0741223422 / +0.0741223422 |
| `minimum_rod_axis_cosine` ≥ 0.95 | 1 | +0.0265685166 / +0.0265685166 | +0.00856356692 / +0.00856316931 |
| `EH_over_crank` ≤ 7.0 | crank_radius | +2.69173534 / +2.69173534 | +3.69082126 / +3.69082126 |
| `crank_axis_to_EFH_clearance_over_crank` ≥ 0.5 | crank_radius | +4.71763157 / +4.71762996 | +4.31129598 / +4.31129598 |
| `H_axis_lateral_rms_over_stroke` ≤ 0.25 | 1 | +0.0910218379 / +0.0910218379 | +0.118566189 / +0.118566189 |
| `H_axis_lateral_span_over_stroke` ≤ 0.65 | 1 | +0.090546661 / +0.0905456162 | +0.150262294 / +0.150257715 |
| `zero_crossing_count` = 2 | 1 | pass / pass | pass / pass |

---

## 9. Mechanical comparison

The four lineages should not be ranked by one scalar mechanical score.

They span different trade-offs.

| Family | Mechanical character | Main nominal concern |
|---|---|---|
| `compact_balanced` | compact, balanced proportions | LARGE H-lateral-span margin |
| `compact_offset` | compact, high stroke | primary/secondary transmission margins |
| `long_ground` | long-ground / long-dyad topology | SMALL rod/H span and LARGE primary transmission |
| `long_coupler_short_rocker` | long-link geometry with large clearances | SMALL upper stroke margin |

The useful comparison is the vector of mechanical margins, not a single
ranking.

---

## 10. Mechanical margin versus manufacturing robustness

A nominal mechanism can be characterized deterministically by its distance from
declared mechanical limits.

For a minimum constraint,

\[
m=x-x_{\min}.
\]

For a maximum constraint,

\[
m=x_{\max}-x.
\]

Positive margin means that the nominal mechanism passes the declared screen.

Small positive margin identifies an active or nearly active design constraint.

Margins retain their individual units and meanings. They are not summed or
ranked across unlike metrics. Both stroke bounds are reported separately.

This is useful engineering information, but it must not be confused with a
manufacturing-tolerance proof.

A true bounded tolerance certification would require proving the constraints
over an entire multidimensional perturbation domain, not testing an arbitrary
finite set of random or quasi-random samples.

If certified manufacturing robustness becomes necessary, it should be added as
a generic solver/research capability with explicit uncertainty bounds.

---

## 11. How these families should be reused

These families are starting structures, not immutable mechanisms.

For a new machine:

1. choose several mechanically different lineages;
2. preserve their discrete branches initially;
3. adapt the continuous geometry to the new piston motion;
4. inspect the mechanical margins with the production diagnostics;
5. allow both cylinder mechanisms to adapt independently;
6. compare the resulting machines using the actual thermodynamic objective.

Do not assume that:

- the most compact family is best;
- the family with the largest nominal transmission margin is best;
- a lineage optimal for one operating point is optimal for another;
- the original seed coordinates should remain fixed.

---

## 12. Scaling toward a demonstrator

After selecting a family and a crank radius \(R\):

1. scale every absolute linkage length by \(R\);
2. compute physical SMALL/LARGE strokes from `stroke_over_crank`;
3. choose bore from required swept volume;
4. design bearings and pivots from real force and friction requirements;
5. verify physical clearances in the actual packaging;
6. revisit any nearly active normalized design margins;
7. evaluate manufacturing tolerance using a proper bounded analysis if needed.

The catalogue should remain normalized.

A particular demonstrator geometry is a separate engineering instantiation.

---

## 13. What should not be inferred

This catalogue does not establish that:

- only four useful six-bar families exist;
- one of these lineages is globally optimal;
- nominal mechanical margin alone predicts thermodynamic performance;
- a 0.30 transmission sine is a universal hardware limit;
- the H-lateral limits are universal packaging requirements;
- the deterministic nominal screen is a manufacturing-tolerance proof;
- a family should remain geometrically unchanged when reused.

The catalogue records useful normalized mechanisms and their present
mechanical design margins.

---

## 14. Reproducing the catalogue

The reproducible material is:

```text
docs/repro/six_bar_mechanism_families/families.toml
docs/repro/six_bar_mechanism_families/design_screen.toml
docs/repro/six_bar_mechanism_families/verify.py
```

Run:

```bash
PYTHONPATH=src python docs/repro/six_bar_mechanism_families/verify.py
```

Regression-style check:

```bash
PYTHONPATH=src python docs/repro/six_bar_mechanism_families/verify.py --check
```

Generate the diagnostic Markdown tables:

```bash
PYTHONPATH=src python docs/repro/six_bar_mechanism_families/verify.py --markdown
```

The verifier:

1. loads the eight exact normalized mechanisms;
2. constructs production `SixBarCylinderMechanism` objects;
3. computes diagnostics through production `six_bar_metrics()`;
4. applies the declared standard mechanical screen at 1440 samples;
5. repeats the screen at 5760 samples;
6. reports individual constraint margins;
7. evaluates the stored geometries directly through production mechanisms.

A failure prints both grid values and their verdicts and returns a nonzero exit
status. All mechanisms are evaluated before reporting screen failures; neither
geometry nor limits are adjusted to obtain a pass.
