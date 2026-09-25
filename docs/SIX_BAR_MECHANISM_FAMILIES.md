# Paired six-bar mechanism families from the candidate-3952 synthesis

## Purpose

This document is a mechanism catalogue, separate from the thermodynamic machine description.

It records the four paired SMALL/LARGE six-bar families carried to the end of the candidate-3952 campaign, using the final geometry obtained after simultaneous 30-dimensional thermo-mechanical adaptation.

The subsequent five-dimensional thermodynamic re-tuning changed only the cylinder swept-volume ratio and exchanger dimensions. **It did not change the six-bar geometry.**

The historical family identifiers `1`, `4`, `12` and `50` are retained for reproducibility. They should eventually be replaced or supplemented by stable descriptive family names.

---

## 1. Topology and coordinate convention

Each piston is driven by:

```text
primary four-bar:  A-B-C-D
coupler point:     E fixed on BC
secondary dyad:    E-F-G
output point:      H fixed on EF
piston rod:        H-P
slider:            P constrained to a straight axis
```

All lengths are in **crank-radius units**, with

\[
AB = 1.
\]

The fixed primary pivots are:

\[
A=(0,0),\qquad D=(AD,0).
\]

The crank point is

\[
B=(\cos(\theta+\phi),\sin(\theta+\phi)).
\]

Point `C` is obtained from the primary four-bar closure using the stored primary branch.

Point `E` is a rigid coupler point:

\[
E =
B
+
\frac{E_{along}}{BC}(C-B)
+
\frac{E_{normal}}{BC}R_{90}(C-B).
\]

The second fixed pivot is

\[
G=(G_x,G_y).
\]

Point `F` is obtained from the `E-F-G` RR closure using the stored secondary branch.

Point `H` is a rigid point in the `EF` frame:

\[
H =
E
+
h_{along}(F-E)
+
h_{normal}R_{90}(F-E),
\]

where the stored `H_along / EF` and `H_normal / EF` are dimensionless fractions of the vector `EF`.

The piston rod connects `H` to the positive slider solution used by the solver. The slider axis is defined by:

- `slider_axis_angle`;
- `slider_axis_offset`.

The exact implementation is in `src/dada_solver/six_bar.py`.

---

## 2. Scaling to a physical mechanism

The catalogue geometry is dimensionless.

If a physical crank radius `R` is chosen:

- multiply `AD`, `BC`, `CD`, `E_along`, `E_normal`, `G_x`, `G_y`, `EF`, `GF`, piston-rod length and slider-axis offset by `R`;
- leave `H_along / EF`, `H_normal / EF`, phases, axis angles and branch signs unchanged;
- physical piston stroke is

\[
stroke = (stroke/crank)\,R.
\]

Cylinder bore is then chosen from the required swept volume and physical stroke. The normalized linkage itself does not determine bore.

This makes each family reusable across different physical machine sizes.

---

## 3. Reference machine and final performance

The abstract reference is candidate 3952 at ΔT = 260 K:

- indicated efficiency: **23.5081%**;
- indicated power: **40.719 W**.

After the full mechanism search, the paired six-bars were first evaluated with the **unchanged candidate-3952 thermodynamic machine**. They were then frozen and given a 192-evaluation five-dimensional thermodynamic re-tuning over:

\[
V_S/V_L,\quad n_i,\quad L_i,\quad n_o,\quad L_o.
\]

Final comparison:

| Family | Working description | Pair RMS vs 3952 after 30D | Efficiency with 3952 thermo | Power with 3952 thermo | Efficiency after thermo 5D | Final power |
|---|---|---:|---:|---:|---:|---:|
| 1 | High-efficiency compact | 1.292% | 23.1552% | 41.383 W | **23.2553%** | **41.434 W** |
| 50 | Robust long-link | 3.133% | 22.5328% | 41.142 W | **22.6516%** | **44.204 W** |
| 4 | Compact high-stroke | 2.434% | 22.2229% | 40.793 W | **22.4499%** | **42.695 W** |
| 12 | Alternate long-dyad | 4.777% | 21.2072% | 40.871 W | **21.6567%** | **45.347 W** |

The final thermodynamic ranking remained:

\[
1 > 50 > 4 > 12.
\]

Family 1 retains approximately **98.9% of the indicated efficiency** of abstract candidate 3952 while producing slightly greater indicated power.

---

## 4. Family 1 — high-efficiency compact pair

Family 1 is the strongest result of the campaign.

It descended from the compact balanced primary family and remained relatively close to the 3952 motion while allowing enough thermo-mechanical deformation to increase power.

### SMALL mechanism

| Parameter | Value |
|---|---:|
| `AB` | 1.000000000 |
| `AD` | 1.617198584 |
| `BC` | 1.676540573 |
| `CD` | 1.494053414 |
| `E_along` | 1.762050571 |
| `E_normal` | -0.234670317 |
| primary phase | -2.999365297 rad (-171.851°) |
| primary branch | +1 |
| `G_x` | 3.645801455 |
| `G_y` | 5.166424699 |
| `EF` | 1.224201630 |
| `GF` | 4.387179697 |
| `H_along / EF` | 1.965266625 |
| `H_normal / EF` | -0.315841720 |
| piston rod | 12.222402682 |
| slider-axis offset | 0.459386223 |
| slider-axis angle | 2.638295047 rad (151.163°) |
| secondary branch | -1 |

Mechanical indicators:

| Indicator | Value |
|---|---:|
| piston position RMS vs 3952 | 1.369% |
| stroke / crank | 2.5409 |
| min primary transmission sine | 0.3660 |
| min secondary transmission sine | 0.3501 |
| min rod/axis cosine | 0.9625 |
| H lateral RMS / stroke | 0.1397 |
| H lateral span / stroke | 0.4723 |
| crank-axis clearance / crank | 0.7776 |

### LARGE mechanism

| Parameter | Value |
|---|---:|
| `AB` | 1.000000000 |
| `AD` | 1.545494496 |
| `BC` | 1.658759773 |
| `CD` | 1.473334389 |
| `E_along` | 1.775086675 |
| `E_normal` | 0.324249174 |
| primary phase | -0.432944219 rad (-24.806°) |
| primary branch | -1 |
| `G_x` | 4.252421026 |
| `G_y` | -4.318082045 |
| `EF` | 1.189940540 |
| `GF` | 4.128862386 |
| `H_along / EF` | 2.045669719 |
| `H_normal / EF` | 0.127167224 |
| piston rod | 11.838247226 |
| slider-axis offset | -0.354303468 |
| slider-axis angle | -2.559725498 rad (-146.661°) |
| secondary branch | +1 |

Mechanical indicators:

| Indicator | Value |
|---|---:|
| piston position RMS vs 3952 | 1.211% |
| stroke / crank | 2.3885 |
| min primary transmission sine | 0.3237 |
| min secondary transmission sine | 0.3614 |
| min rod/axis cosine | 0.9621 |
| H lateral RMS / stroke | 0.2013 |
| H lateral span / stroke | 0.6500 |
| crank-axis clearance / crank | 0.8045 |

### Thermodynamic re-tuning

```text
V_S / V_L = 0.877401425
n_i       = 4029
L_i       = 38.337 mm
n_o       = 4555
L_o       = 23.608 mm
```

Final result:

- efficiency: **23.2553%**;
- indicated power: **41.434 W**.

The LARGE H lateral span is essentially at the retained `0.65 stroke` ceiling, while the relaxed H-lateral RMS remains below 0.25. This family therefore deserves a tolerance / robustness study before being selected as the demonstrator.

---

## 5. Family 4 — compact high-stroke pair

This family is visually compact and produces nearly three crank radii of piston stroke on both sides.

Its main limitation is mechanical margin: both primary and secondary transmission sines sit almost exactly on the retained `0.30` floor.

### SMALL mechanism

| Parameter | Value |
|---|---:|
| `AB` | 1.000000000 |
| `AD` | 1.788104778 |
| `BC` | 2.457775120 |
| `CD` | 2.074103828 |
| `E_along` | 2.378854821 |
| `E_normal` | -1.023122787 |
| primary phase | -3.035510851 rad (-173.922°) |
| primary branch | +1 |
| `G_x` | 4.902442164 |
| `G_y` | 2.263273163 |
| `EF` | 1.382108631 |
| `GF` | 2.528781398 |
| `H_along / EF` | 2.060544569 |
| `H_normal / EF` | 0.255214366 |
| piston rod | 6.215058044 |
| slider-axis offset | -2.301688411 |
| slider-axis angle | 2.371003621 rad (135.849°) |
| secondary branch | -1 |

Mechanical indicators:

| Indicator | Value |
|---|---:|
| piston position RMS vs 3952 | 1.658% |
| stroke / crank | 2.9793 |
| min primary transmission sine | 0.3013 |
| min secondary transmission sine | 0.3001 |
| min rod/axis cosine | 0.9600 |
| H lateral RMS / stroke | 0.1672 |
| H lateral span / stroke | 0.6238 |
| crank-axis clearance / crank | 1.5895 |

### LARGE mechanism

| Parameter | Value |
|---|---:|
| `AB` | 1.000000000 |
| `AD` | 1.747675773 |
| `BC` | 2.343788684 |
| `CD` | 1.988484905 |
| `E_along` | 2.406882286 |
| `E_normal` | 1.219215768 |
| primary phase | -0.338533420 rad (-19.397°) |
| primary branch | -1 |
| `G_x` | 4.994620641 |
| `G_y` | -1.902805284 |
| `EF` | 1.503362371 |
| `GF` | 2.221524054 |
| `H_along / EF` | 1.991078657 |
| `H_normal / EF` | -0.385685686 |
| piston rod | 6.379183253 |
| slider-axis offset | 2.273241501 |
| slider-axis angle | -2.436954743 rad (-139.627°) |
| secondary branch | +1 |

Mechanical indicators:

| Indicator | Value |
|---|---:|
| piston position RMS vs 3952 | 3.017% |
| stroke / crank | 2.9630 |
| min primary transmission sine | 0.3012 |
| min secondary transmission sine | 0.3003 |
| min rod/axis cosine | 0.9547 |
| H lateral RMS / stroke | 0.1634 |
| H lateral span / stroke | 0.6377 |
| crank-axis clearance / crank | 1.6981 |

### Thermodynamic re-tuning

```text
V_S / V_L = 0.888562498
n_i       = 4669
L_i       = 38.288 mm
n_o       = 4417
L_o       = 26.227 mm
```

Final result:

- efficiency: **22.4499%**;
- indicated power: **42.695 W**.

The family gained thermodynamically when the exchanger dimensions were re-tuned, but the geometry is already close to several mechanical limits. It is therefore a useful compactness benchmark rather than an obvious prototype choice.

---

## 6. Family 12 — alternate long-dyad pair

Family 12 is the most geometrically different of the four retained mechanisms.

The SMALL side uses a very long `EF` link and a relatively short piston stroke. During thermo-mechanical adaptation it moved much farther from the candidate-3952 motion than the other families.

This family was especially useful for identifying an artificial geometric restriction: relaxing `H lateral RMS / stroke` from 0.20 to 0.25 produced a substantial efficiency gain. The continuation then encountered more direct mechanical limits instead.

### SMALL mechanism

| Parameter | Value |
|---|---:|
| `AB` | 1.000000000 |
| `AD` | 3.410671043 |
| `BC` | 2.461702291 |
| `CD` | 2.106069045 |
| `E_along` | 3.916466385 |
| `E_normal` | -1.090131795 |
| primary phase | -0.288816327 rad (-16.548°) |
| primary branch | -1 |
| `G_x` | 6.908484590 |
| `G_y` | -0.004320118 |
| `EF` | 6.241298210 |
| `GF` | 1.883153359 |
| `H_along / EF` | 0.591559383 |
| `H_normal / EF` | 0.282313429 |
| piston rod | 8.671128856 |
| slider-axis offset | 7.435944155 |
| slider-axis angle | -1.817280951 rad (-104.123°) |
| secondary branch | -1 |

Mechanical indicators:

| Indicator | Value |
|---|---:|
| piston position RMS vs 3952 | 5.579% |
| stroke / crank | 1.2387 |
| min primary transmission sine | 0.5035 |
| min secondary transmission sine | 0.5309 |
| min rod/axis cosine | 0.9505 |
| H lateral RMS / stroke | 0.2275 |
| H lateral span / stroke | 0.6498 |
| crank-axis clearance / crank | 3.0654 |

### LARGE mechanism

| Parameter | Value |
|---|---:|
| `AB` | 1.000000000 |
| `AD` | 3.655982759 |
| `BC` | 2.493001630 |
| `CD` | 2.217482623 |
| `E_along` | 4.295814511 |
| `E_normal` | 1.241403635 |
| primary phase | 3.132932234 rad (179.504°) |
| primary branch | +1 |
| `G_x` | 6.354986000 |
| `G_y` | 0.024209397 |
| `EF` | 5.617216794 |
| `GF` | 2.042105441 |
| `H_along / EF` | 0.741397336 |
| `H_normal / EF` | -0.245730621 |
| piston rod | 9.075494111 |
| slider-axis offset | -7.163890758 |
| slider-axis angle | 1.694698020 rad (97.099°) |
| secondary branch | +1 |

Mechanical indicators:

| Indicator | Value |
|---|---:|
| piston position RMS vs 3952 | 3.810% |
| stroke / crank | 1.9004 |
| min primary transmission sine | 0.3004 |
| min secondary transmission sine | 0.3174 |
| min rod/axis cosine | 0.9816 |
| H lateral RMS / stroke | 0.1845 |
| H lateral span / stroke | 0.5497 |
| crank-axis clearance / crank | 3.4716 |

### Thermodynamic re-tuning

```text
V_S / V_L = 0.947951242
n_i       = 5039
L_i       = 41.990 mm
n_o       = 4889
L_o       = 26.007 mm
```

Final result:

- efficiency: **21.6567%**;
- indicated power: **45.347 W**.

The SMALL side is almost on both the rod-angle and H-span limits, while the LARGE primary transmission is almost on the `0.30` floor. The mechanism therefore appears genuinely constrained rather than merely under-optimized.

---

## 7. Family 50 — robust long-link pair

Family 50 descends from the older successful six-bar lineage.

It is much larger geometrically than family 1 but retains substantially better transmission margins and large crank-axis clearances.

It finished second in thermodynamic efficiency after the final re-tuning.

### SMALL mechanism

| Parameter | Value |
|---|---:|
| `AB` | 1.000000000 |
| `AD` | 4.906167614 |
| `BC` | 4.863312784 |
| `CD` | 1.149903793 |
| `E_along` | 4.199153543 |
| `E_normal` | 4.585414961 |
| primary phase | -0.252757411 rad (-14.482°) |
| primary branch | -1 |
| `G_x` | 9.198575716 |
| `G_y` | -5.704989087 |
| `EF` | 2.042563560 |
| `GF` | 9.446423386 |
| `H_along / EF` | 1.628000624 |
| `H_normal / EF` | 1.341090474 |
| piston rod | 8.474157093 |
| slider-axis offset | 0.462288981 |
| slider-axis angle | 0.618786648 rad (35.454°) |
| secondary branch | +1 |

Mechanical indicators:

| Indicator | Value |
|---|---:|
| piston position RMS vs 3952 | 3.530% |
| stroke / crank | 2.9731 |
| min primary transmission sine | 0.4638 |
| min secondary transmission sine | 0.3634 |
| min rod/axis cosine | 0.9766 |
| H lateral RMS / stroke | 0.1590 |
| H lateral span / stroke | 0.5595 |
| crank-axis clearance / crank | 5.2176 |

### LARGE mechanism

| Parameter | Value |
|---|---:|
| `AB` | 1.000000000 |
| `AD` | 5.077061603 |
| `BC` | 4.993700452 |
| `CD` | 1.153574664 |
| `E_along` | 3.971069226 |
| `E_normal` | -4.242849155 |
| primary phase | -3.083369538 rad (-176.664°) |
| primary branch | +1 |
| `G_x` | 8.871255866 |
| `G_y` | 6.087378747 |
| `EF` | 1.837243263 |
| `GF` | 9.769509441 |
| `H_along / EF` | 1.335066412 |
| `H_normal / EF` | -1.209046718 |
| piston rod | 7.338631363 |
| slider-axis offset | -0.839069328 |
| slider-axis angle | -0.703987242 rad (-40.335°) |
| secondary branch | -1 |

Mechanical indicators:

| Indicator | Value |
|---|---:|
| piston position RMS vs 3952 | 2.680% |
| stroke / crank | 2.5596 |
| min primary transmission sine | 0.3787 |
| min secondary transmission sine | 0.3741 |
| min rod/axis cosine | 0.9586 |
| H lateral RMS / stroke | 0.1314 |
| H lateral span / stroke | 0.4997 |
| crank-axis clearance / crank | 4.8113 |

### Thermodynamic re-tuning

```text
V_S / V_L = 0.889357327
n_i       = 4403
L_i       = 44.708 mm
n_o       = 4622
L_o       = 27.030 mm
```

Final result:

- efficiency: **22.6516%**;
- indicated power: **44.204 W**.

This family is a particularly useful counterexample to ranking mechanisms by primary cadence score alone. Its original primary proxy score was poor, yet the complete six-bar remained mechanically robust and thermodynamically competitive.

---

## 8. Evolution from kinematic fit to thermo-mechanical optimum

Before direct thermo-mechanical optimization, the mirrored and polished pairs had combined position RMS errors approximately:

| Family | Combined RMS after mirror + local polish |
|---|---:|
| 1 | **0.984%** |
| 4 | 2.389% |
| 12 | 2.596% |
| 50 | 1.930% |

After simultaneous 30-D optimization against the fixed 3952 thermodynamics:

| Family | Combined RMS after thermo-mechanical adaptation |
|---|---:|
| 1 | 1.292% |
| 4 | 2.434% |
| 12 | 4.777% |
| 50 | 3.133% |

The thermodynamic search was explicitly allowed to move away from the target motion.

This is why family 12 could become thermodynamically better while its geometric RMS became much worse.

The result supports the intended hierarchy:

\[
\text{position fit} \rightarrow \text{synthesis guide},
\]

then

\[
\text{periodic thermodynamic efficiency} \rightarrow \text{machine objective}.
\]

---

## 9. Comparison of mechanical character

The four families occupy noticeably different regions of mechanical design space.

| Family | Main strength | Main weakness / active margin |
|---|---|---|
| 1 | highest efficiency, compact, excellent target reproduction | LARGE H lateral span near 0.65 |
| 4 | compact, very large stroke/crank | both transmission stages near 0.30 |
| 12 | radically different topology, demonstrates thermo adaptability | SMALL rod angle + H span, LARGE primary transmission |
| 50 | robust transmissions, large clearances, strong power | larger overall linkage geometry |

This table is descriptive, not a prototype selection.

Prototype choice should also include:

- dimensional-tolerance sensitivity;
- bearing loads and inertial forces;
- collision and axial-layering layout;
- manufacturability;
- friction and useful shaft efficiency.

---

## 10. Completed robustness campaign

The four final pairs were subjected to the same manufacturing-tolerance experiment.

For every family and each relative tolerance level (**±0.1%, ±0.5%, ±1%**):

- 2048 scrambled-Sobol 30-D perturbations were generated;
- normalized perturbation coordinates were identical across families;
- dimension-like variables received bounded independent uniform errors;
- phase and slider-axis angles reached **±0.5° at the ±1% level**;
- crank radius and assembly branches remained fixed;
- the perturbed mechanism was **not re-optimized**;
- mechanically accepted variants were replayed with that family's final 5-D thermodynamic machine.

The mechanical-screen survival fractions were:

| Family | ±0.1% | ±0.5% | ±1% |
|---|---:|---:|---:|
| 1 | **50.3%** | **36.8%** | **17.4%** |
| 4 | **24.7%** | **9.3%** | **3.8%** |
| 12 | **26.0%** | **19.6%** | **12.6%** |
| 50 | **69.5%** | **37.6%** | **28.9%** |

The main rejection mechanisms were consistent with the nominal active margins:

- **family 1:** LARGE `H` lateral span first, then secondary transmission;
- **family 4:** both transmission stages and stroke limits;
- **family 12:** SMALL `H` span / rod angle and LARGE primary transmission;
- **family 50:** SMALL stroke at tight tolerance, then closure / transmission at larger perturbations.

These are failures of the retained **design screen**, not necessarily physical breakages. In particular, `H`-span and stroke-window limits are partly packaging/design choices.

### 10.1 Thermodynamic robustness after the mechanical screen

Every selected thermodynamic replay converged and remained feasible.

Mean loss of indicated efficiency, in percentage points relative to the nominal family:

| Family | ±0.1% | ±0.5% | ±1% |
|---|---:|---:|---:|
| 1 | **0.017 pp** | 0.104 pp | 0.306 pp |
| 4 | 0.054 pp | 0.238 pp | 0.427 pp |
| 12 | **0.016 pp** | **0.047 pp** | **0.085 pp** |
| 50 | 0.034 pp | 0.168 pp | 0.562 pp |

Fraction remaining within **0.1 percentage point** of nominal efficiency:

| Family | ±0.1% | ±0.5% | ±1% |
|---|---:|---:|---:|
| 1 | **100%** | 55.5% | 15.6% |
| 4 | 97.7% | 14.1% | 3.9% |
| 12 | **100%** | **92.2%** | **62.5%** |
| 50 | 96.1% | 37.5% | 7.8% |

The experiment separates two notions that were previously conflated:

- **mechanical survival:** staying inside the retained closure, transmission, rod-angle, stroke and packaging margins;
- **thermodynamic survival:** retaining efficiency once the perturbed mechanism is mechanically acceptable.

Family 50 is the strongest mechanical survivor under this screen, but family 12 is the least thermodynamically sensitive among the surviving variants at ±0.5% and ±1%.

Family 1 remains the nominal efficiency champion. Its main weakness in this robustness test is not thermodynamic fragility at small error, but the fact that the optimized LARGE mechanism sits almost directly on the `H` lateral-span ceiling.

Family 4 is confirmed as the most tolerance-critical of the four because both transmission stages were already optimized close to the 0.30 floor.

### 10.2 Implication for the demonstrator

The robustness test does not produce a single automatic winner.

For the demonstrator, the next engineering step is to decide which active limits are real hardware requirements and which can be relaxed through packaging, bearing layout or modest mechanism re-optimization.

The current data suggest that a small nominal-efficiency sacrifice could buy substantial manufacturing margin, but that trade should be made only after replacing the present normalized screening limits by actual dimensional tolerances and component constraints.

Source artifacts:

```text
examples/test_sixbar_robustness_3952.py
outputs/sixbar_robustness_3952/definition.json
outputs/sixbar_robustness_3952/report.json
```

---

## 11. How to create a scaled demonstrator

Once a family and a crank radius `R` are chosen:

1. scale every absolute mechanism length by `R`;
2. compute physical S/L strokes from the stored stroke/crank values;
3. use the selected thermodynamic swept volumes to compute piston areas and bores;
4. place the cylinders along the stored slider axes;
5. design bearings, shafts, link thickness and axial layering from the actual forces;
6. verify collision clearances with physical component thickness;
7. generate the final animation from the scaled geometry.

The mechanism catalogue should remain dimensionless; the wiki demonstrator can then present one specific physical instantiation.

---

## 12. Source artifacts

Paired kinematic synthesis:

```text
outputs/sixbar_pairs_3952/mirror_family_summary.json
outputs/sixbar_pairs_3952/rank_01_small.json
outputs/sixbar_pairs_3952/rank_01_large.json
outputs/sixbar_pairs_3952/rank_04_small.json
outputs/sixbar_pairs_3952/rank_04_large.json
outputs/sixbar_pairs_3952/rank_12_small.json
outputs/sixbar_pairs_3952/rank_12_large.json
outputs/sixbar_pairs_3952/rank_50_small.json
outputs/sixbar_pairs_3952/rank_50_large.json
```

Final simultaneous thermo-mechanical geometry:

```text
outputs/sixbar_thermo_coupled_3952_fine_hlat25/rank_01/best_pair.json
outputs/sixbar_thermo_coupled_3952_fine_hlat25/rank_04/best_pair.json
outputs/sixbar_thermo_coupled_3952_fine_hlat25/rank_12/best_pair.json
outputs/sixbar_thermo_coupled_3952_fine_hlat25/rank_50/best_pair.json
```

Final thermodynamic re-tuning:

```text
outputs/sixbar_thermo5d_3952/report.json
```

Robustness:

```text
examples/test_sixbar_robustness_3952.py
outputs/sixbar_robustness_3952/definition.json
outputs/sixbar_robustness_3952/report.json
```

Mechanism implementation:

```text
src/dada_solver/six_bar.py
```

Search methodology:

```text
docs/MECHANISM_SYNTHESIS_SEARCH_THEORY.md
```

Primary four-bar catalogue:

```text
docs/PRIMARY_FOUR_BAR_FAMILIES.md
```
