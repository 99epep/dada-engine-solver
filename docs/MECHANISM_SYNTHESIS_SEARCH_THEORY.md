# Mechanism-synthesis search theory for the DADA engine

## Purpose

This note records the mechanism-synthesis method developed while transforming an optimized thermodynamic piston law into mechanically realizable four-bar / six-bar linkages.

The reference case is structured C2 candidate 3952 at a source-temperature difference of 260 K. The numerical values in this document belong to that case, but the method is intended to be reusable for other temperatures, machine sizes, working gases and operating points.

The main result is methodological:

> Mechanism synthesis should be treated as a **hierarchical basin-discovery and thermo-mechanical adaptation problem**, not as a blind fit of all linkage dimensions at once.

The completed 3952 campaign established a practical sequence:

1. identify the piston carrying the most restrictive kinematic demand;
2. discover a diverse library of primary four-bar families;
3. give selected primaries a fair downstream-dyad search;
4. locally release the complete six-bar geometry;
5. mirror the successful mechanism toward the opposite piston as an initialization only;
6. locally adapt the mirrored mechanism;
7. release both six-bar mechanisms simultaneously and optimize them directly on thermodynamic efficiency;
8. with the mechanisms fixed, re-tune the small set of thermodynamic / exchanger dimensions;
9. use a separate fresh-island saturation experiment to quantify how well the primary search space has been mapped.

This sequence is now supported by a completed 512-island primary saturation campaign and by full thermodynamic replay of four six-bar families.

---

## 1. Thermodynamic target and mechanical mechanism are different objects

The optimized thermodynamic motion is a design probe. It is not automatically a mechanically meaningful trajectory.

The realizable linkage should preserve the thermodynamic effects that matter while remaining free to smooth, redistribute or reshape features that are artifacts of the chosen target parameterization.

This distinction became especially important for candidate 3952. Its large-cylinder target contains a severe high-pressure (HP) acceleration feature near the HP transition. Automatic detection identified the unreliable region approximately as:

- start: **34.0°**;
- end: **44.25°**;
- dominant jerk fronts: approximately **34.25°, 36.5°, 38.75°, 39.5°, 41.75°, 43.75°**.

The position law remains meaningful through this interval. The detailed acceleration shape does not.

The corresponding working rules are:

1. **position is the primary kinematic truth**;
2. velocity is useful as a cadence descriptor and guarded diagnostic;
3. target acceleration should not be directly fitted through an artificial HP feature;
4. the linkage is allowed, and expected, to round such features naturally;
5. the actual linkage must finally be judged by the thermodynamic solver.

A previous downstream experiment showed why the first point matters. Removing the HP interval from the position objective allowed the optimizer to insert a pause and slight reversal while retaining an apparently good score. A noisy derivative interval is therefore not permission to ignore the underlying position law.

---

## 2. Reference thermodynamic case: candidate 3952

The mechanism campaign was anchored to candidate 3952 from the 260 K structured C2 search.

Reference performance:

- indicated thermal efficiency: **23.508063%**;
- indicated power: **40.7191 W**;
- heat input: approximately **173.21 W**.

The source machine and thermodynamic state were kept fixed during the first thermo-coupled mechanism search. In particular:

- source-temperature difference remained **260 K**;
- exchanger model and working-fluid model remained unchanged;
- total swept volume remained unchanged;
- total gas inventory remained fixed;
- the power floor remained **25 W**;
- only mechanism geometry was initially allowed to move.

This distinction is important. The first thermo-mechanical search asked:

> How much of candidate 3952 can be recovered by actual six-bar mechanisms without changing the machine thermodynamics?

Only afterward was a separate five-dimensional thermodynamic re-tuning performed.

---

## 3. Six-bar topology and variables

The selected topology consists of a primary four-bar `A-B-C-D`, with point `E` rigidly attached to the coupler `BC`, followed by a secondary dyad `E-F-G`. Point `H` is rigidly attached to link `EF`, and the piston is connected from `H` by a finite rod to a slider axis.

The crank length `AB` is normalized to 1.

Each piston mechanism has 15 continuous variables.

### Primary four-bar and coupler point E

1. `AD / AB`;
2. `BC / AB`;
3. `CD / AB`;
4. `E_along`;
5. `E_normal`;
6. primary crank phase.

### Secondary dyad, point H and piston slider

7. second fixed-pivot x;
8. second fixed-pivot y;
9. `EF`;
10. `GF`;
11. `H_along / EF`;
12. `H_normal / EF`;
13. piston-rod length;
14. slider-axis offset;
15. slider-axis angle.

There are also two discrete assembly branches:

- primary four-bar branch;
- secondary dyad branch.

For a paired SMALL/LARGE machine this gives **30 continuous mechanism dimensions** once branches are fixed.

---

## 4. Why blind 15-D synthesis was inefficient

A long blind 15-dimensional six-bar search consumed millions of evaluations without improving the previously known large-cylinder mechanism.

The failure exposed several structural facts:

- the feasible region is thin relative to the raw 15-D box;
- closure, transmission and slider constraints create large invalid regions;
- flat penalties give differential evolution little useful direction;
- the primary four-bar contains structure that can be discovered much more cheaply on its own;
- the downstream dyad can substantially transform the primary motion, so fitting everything simultaneously from random geometry wastes search effort.

This led to the hierarchical architecture described below.

---

## 5. Hierarchical synthesis architecture

### 5.1 Stage A — primary four-bar library

Search only the six primary variables

\[
x_p=(g,c,r,E_\parallel,E_\perp,\phi).
\]

This stage is cheap enough to support hundreds of independent global islands.

Its objective is not to reproduce the final piston law exactly. The primary should generate a mechanically useful temporal cadence and a trajectory for point `E` that the downstream dyad can exploit.

### 5.2 Stage B — downstream dyad

For each selected primary, freeze the six primary dimensions and optimize the remaining nine continuous variables:

\[
(G_x,G_y,EF,GF,H_\parallel/EF,H_\perp/EF,L_{rod},s_{slider},\alpha_{slider})
\]

plus the secondary branch.

The downstream dyad is intentionally allowed to reshape the primary motion.

### 5.3 Stage C — full 15-D local polish

Once a viable primary + dyad combination exists, release all 15 continuous variables within a moderate local trust region.

This stage proved essential. Small changes in the primary geometry produced large improvements in the final piston fit that the fixed-primary downstream stage could not obtain.

### 5.4 Stage D — mirror toward the opposite piston

The successful mechanism for the harder piston is physically mirrored and phase-aligned toward the other piston.

The mirror is only an **initialization**. It is not retained as a symmetry constraint.

### 5.5 Stage E — local adaptation of the mirrored mechanism

The mirrored mechanism receives its own 15-D local polish while preserving its discrete branches.

### 5.6 Stage F — paired 30-D thermo-mechanical optimization

Both SMALL and LARGE six-bar mechanisms are then released **simultaneously** in all 30 continuous dimensions.

At this stage:

- the thermodynamic machine remains fixed to candidate 3952;
- mechanism branches remain fixed within each family;
- the objective is **actual indicated thermal efficiency**;
- RMS distance from the 3952 piston law becomes diagnostic only.

This is a critical transition. Once a realizable mechanism exists, geometric fidelity is no longer the final optimization objective.

### 5.7 Stage G — thermodynamic 5-D re-tuning

Finally, freeze both six-bar mechanisms and re-optimize only:

\[
V_S/V_L,\quad n_i,\quad L_i,\quad n_o,\quad L_o.
\]

This measures how much performance was lost because a new mechanism was still using thermodynamic hardware tuned for the abstract 3952 motion.

---

## 6. Primary four-bar feasible space

The normalized primary search uses `AB = 1` and six continuous parameters:

\[
(g,c,r,E_\parallel,E_\perp,\phi).
\]

The broad bounds used by the completed saturation experiment were:

| Parameter | Bound |
|---|---:|
| `AD / AB` | 0.50 .. 12.0 |
| `BC / AB` | 0.50 .. 12.0 |
| `CD / AB` | 0.35 .. 6.0 |
| `E_along` | -12.0 .. 18.0 |
| `E_normal` | -15.0 .. 15.0 |
| phase | -π .. +π |

Both primary assembly branches are admissible.

Full-revolution closure is checked.

Primary transmission quality is measured by

\[
s_{min}
=
\min_\theta
\frac{|BC\times DC|}{|BC||DC|}
=
\min_\theta |\sin \mu|.
\]

A hard floor of **0.30** was retained, corresponding to approximately **17.5° from toggle**.

Later primary searches used a soft preference near 0.35 instead of tightening the hard floor. This distinction is deliberate:

> A mechanical preference should not silently become a topological exclusion unless a real physical limit justifies it.

---

## 7. Evolution of the primary-cadence objective

The sequence V1 → V4.1 is important because it records several ways in which an apparently reasonable mathematical proxy can be physically wrong.

### 7.1 V1 — acceleration-lobe timing

V1 used the magnitude of `E` velocity and tangential acceleration

\[
a_t
=
\frac{E'\cdot E''}{\|E'\|}
=
\frac{d\|E'\|}{d\theta},
\]

then attempted to compare those events with signed piston acceleration.

This was conceptually inconsistent: the derivative of a non-negative speed magnitude is not signed piston acceleration.

V1 nevertheless generated useful geometric diversity and exposed a strong optimizer tendency to exploit the transmission floor.

### 7.2 V2 — signed PCA velocity and local acceleration events

V2 projected `E` motion onto its PCA principal axis and used signed projected velocity.

This fixed the sign problem but revealed two further limitations:

- the PCA axis is not physically equivalent to the final piston axis because the second dyad can rotate and transform motion;
- event-based penalties could nearly eliminate a nominally required event without producing a large objective penalty.

### 7.3 V3 — automatic three-regime velocity segmentation

V3 replaced isolated acceleration peaks with three automatically discovered piecewise-linear velocity regimes.

This was structurally better, but the statistically optimal segmentation crossed the real piston reversal near 154°, mixing opposite motion directions into one fitted phase.

The lesson was fundamental:

> A segmentation can minimize approximation error while violating the topology of the physical motion.

### 7.4 V4 — preserve the actual monotonic branches

V4 starts from the two true velocity zero crossings.

For candidate 3952, after bridging the unreliable HP velocity feature:

- high-position turnaround: **358.812°**;
- low-position turnaround: **154.144°**;
- short high-to-low branch: **155.332°**;
- long low-to-high branch: **204.668°**.

The short branch naturally divides into a fast and slow cadence, with an automatic split near **38.062°**.

The guarded target metrics are approximately:

- fast/slow mean-speed ratio: **5.027**;
- fraction of short-branch displacement completed in the fast sector: **57.63%**.

V4 therefore rewards:

1. two correctly timed turnarounds;
2. no extra reversal;
3. monotonicity of both branches;
4. target-like fast/slow cadence on the short branch;
5. target-like displacement sharing on the short branch;
6. acceptable transmission quality;
7. only weak geometric directionality preferences.

### 7.5 V4.1 — do not over-constrain the BP branch

Earlier motion searches showed that a visibly non-uniform low-pressure exchange can still produce good thermodynamic efficiency.

The individual piston speed during BP exchange therefore should not be forced to match the detailed 3952 profile.

V4.1 removed the long-branch target-profile penalty. The long branch may be approximately sinusoidal, fast-slow-fast, or otherwise non-uniform, provided the important motion topology remains correct.

The retained BP proxy is approximate mirror symmetry of the candidate's own long branch:

\[
A_{BP}
=
\frac{\operatorname{RMS}[v(u)-v(1-u)]}
{\sqrt{\operatorname{mean}((v(u)^2+v(1-u)^2)/2)}}.
\]

This remains only a search proxy.

The real physical question is whether the **combined cylinder-volume evolution** keeps the BP exchange acceptably close to the desired pressure behavior.

---

## 8. Search islands and basin-capture probability

A differential-evolution island should not be interpreted as one small geometric patch.

With six variables and `popsize = 16`, an island begins with approximately

\[
6\times16=96
\]

Latin-hypercube individuals distributed across the complete bounded domain.

The population then evolves globally to locally.

Therefore the useful statistical quantity is not literal geometric basin volume but **algorithmic basin-capture probability**:

\[
P_i^{capture}
=
P(\text{one fresh DE island terminates in family }i).
\]

For a fixed algorithm, bounds, constraints and budget, this can be estimated empirically by

\[
\hat P_i=\frac{n_i}{N}.
\]

This quantity is search-policy dependent. It must not be described as the literal fraction of geometric mechanism space occupied by a family.

---

## 9. Operational definition of a primary family

A primary family is defined from normalized geometric distance in the six primary variables, with opposite assembly branches kept separate.

Because any clustering threshold is conventional, the saturation experiment reports nearby distance thresholds:

- **0.03**;
- **0.04**;
- **0.05**.

Connected-component clustering is used so that the result does not depend on candidate insertion order.

The possibility of single-link chaining remains a reason not to over-interpret the absolute family count.

The stable quantity of greatest practical interest is the concentration of capture probability into the dominant competitive families.

---

## 10. Completed 512-island primary saturation experiment

The final V4.1 saturation campaign used:

- **512 fresh islands**;
- no historical seeds;
- strict alternation of the two primary branches;
- 260 DE generations per island;
- population size 16;
- approximately 96 initial individuals per island;
- the complete primary bounds listed above.

All **512/512** islands produced dense winners.

The campaign consumed approximately **8.09 cumulative island wall-hours**.

### 10.1 Best primary found

The final best score was:

\[
\boxed{0.29863621895}
\]

with branch `-1`.

Its main descriptors are:

- minimum primary transmission sine: **0.34998**;
- `BE / BC`: **1.19639**;
- turnarounds: **154.099°** and **358.781°**;
- fast/slow speed ratio: **5.03485**;
- fast displacement fraction: **0.55685**.

The improvement relative to the earlier best V4.1 result was small. The large saturation campaign therefore increased confidence in the map of the search landscape much more than it changed the best objective value.

### 10.2 Family concentration at distance threshold 0.04

At threshold **0.04**:

- observed families: **16**;
- singletons: **8**;
- doubletons: **3**;
- Good-Turing unseen capture-mass estimate: **1.5625%**;
- top-1 capture share: **47.66%**;
- top-4 capture share: **96.68%**;
- top-10 capture share: **98.83%**.

The four largest families captured:

1. **244 / 512 = 47.66%**;
2. **228 / 512 = 44.53%**;
3. **19 / 512 = 3.71%**;
4. **4 / 512 = 0.78%**.

The first two families alone therefore account for approximately **92.19%** of fresh-island outcomes.

### 10.3 Threshold sensitivity

At threshold **0.05**:

- observed families: **12**;
- Good-Turing unseen mass: **1.17%**;
- top-4 capture share: **98.05%**.

The qualitative conclusion is unchanged:

> The useful primary landscape is dominated by a very small number of large basins, especially two major families, followed by a long tail of rare families.

The campaign does **not** prove that only 12 or 16 physical four-bar families exist. It shows that under the declared search policy, nearly all discoverable capture probability is concentrated in a few basins.

### 10.4 Meaning of the four downstream families

The four six-bar families later retained as historical ranks **1, 4, 12 and 50** should not be interpreted as four statistically equal primary basins.

They were retained because they provided useful **mechanical diversity and downstream transformability**.

This is a central methodological distinction:

> Primary-family capture frequency and final six-bar usefulness are different properties.

A rare or mediocre primary under the simplified cadence proxy can still become valuable after the second dyad is optimized.

---

## 11. Seeds and fresh islands serve different purposes

Historical seeds are useful in a production design campaign because they:

- preserve expensive discoveries;
- accelerate convergence;
- allow changed objectives to reuse earlier geometry;
- protect narrow known basins from being lost.

They are inappropriate for measuring capture probabilities because they bias the initial distribution.

Two campaign types should therefore remain separate.

### Exploitation / design campaign

Use a mixture of:

- historical seeds;
- deliberately diverse known families;
- fresh global populations.

Goal: find useful mechanisms efficiently.

### Saturation / methodology campaign

Use:

- fresh islands only;
- fixed complete bounds;
- fixed DE policy;
- deterministic independent seeds;
- balanced discrete branches.

Goal: measure search saturation and family-capture statistics.

The completed 512-island run belongs to the second category.

---

## 12. Downstream dyad experiment

Several physically distinct primaries were sent to a nine-dimensional downstream search.

The important result was not merely a better fit. It was the demonstration that downstream transformability varies strongly between primary families.

Examples from the panel included:

- a compact, visually attractive family with excellent downstream potential;
- compact candidates near transmission limits;
- a robust historical family with poorer primary proxy score but good six-bar behavior.

This confirmed that the primary objective should identify **promising carrier motions**, not attempt to predict the complete six-bar piston law by itself.

The second dyad is not a minor correction. It is a genuine motion transformer.

---

## 13. Full 15-D LARGE polish

After the fixed-primary downstream stage, all 15 LARGE-side dimensions were released locally.

The most important result was family 1:

- full-cycle LARGE piston-position RMS reduced from approximately **3.30%** to **1.095%**;
- primary transmission sine: **0.3237**;
- secondary transmission sine: **0.3809**;
- rod-axis cosine: **0.9629**;
- stroke/crank: **2.389**.

This large improvement demonstrated that the correct architecture is neither:

- blind 15-D global search from scratch, nor
- permanent freezing of the primary after Stage A.

Instead:

> Discover the primary globally, then release it locally once a viable downstream geometry exists.

Other retained families remained mechanically distinct, which was useful for the later thermo-coupled comparison.

---

## 14. Which piston should be synthesized first?

The 3952 campaign was most successful when synthesis began with the LARGE piston, because its desired law carried the sharper and more restrictive cadence change.

The successful sequence was:

1. solve the more difficult LARGE side;
2. mirror it physically toward SMALL;
3. phase-align the mirror;
4. locally adapt all SMALL dimensions.

The mirror operation changes the transverse signs and reverses the relevant assembly branches while preserving a scalar piston law equivalent to reverse crank traversal.

The physical mirror was preferred because it allows both cranks to appear rotating in the same direction in mechanical visualization.

However, the general lesson is **not** “always synthesize LARGE first.”

For another machine, especially a refrigerator, the correct rule is:

> Identify the piston carrying the sharper or more restrictive kinematic demand, synthesize that side first, then mirror toward the easier side.

---

## 15. Mirrored SMALL synthesis

The polished LARGE families were physically mirrored and adapted to the SMALL target.

The mirror worked unusually well.

Initial phase-only mirrored seeds had SMALL RMS errors of roughly 3–4%, but local 15-D adaptation reduced them substantially.

Final paired kinematic results before thermo-coupled optimization were:

| Historical family | LARGE RMS | SMALL RMS | Combined pair RMS |
|---|---:|---:|---:|
| 1 | 1.095% | **0.859%** | **0.984%** |
| 4 | 2.949% | 1.648% | 2.389% |
| 12 | 3.203% | 1.793% | 2.596% |
| 50 | 2.397% | 1.305% | 1.930% |

All retained SMALL and LARGE mechanisms had exactly two piston reversals and no extra reversal.

The rank-1 pair was the strongest pure kinematic result, but the next stage deliberately stopped using this RMS as the optimization objective.

---

## 16. Mechanical feasibility constraints

The local six-bar searches retained the following broad physical screens:

- stroke/crank between approximately 1 and 3;
- minimum primary transmission sine ≥ **0.30**;
- minimum secondary transmission sine ≥ **0.30**;
- minimum piston-rod / slider-axis cosine ≥ **0.95**;
- `EH` bounded to avoid unbounded lever-arm exploitation;
- minimum crank-axis clearance;
- limits on lateral motion of point H relative to piston stroke.

The H lateral-RMS constraint deserves special mention.

The first thermo-coupled search used:

\[
H_{lat,RMS}/stroke \le 0.20.
\]

Family 12 clearly saturated this boundary. A controlled continuation relaxed only this diagnostic-style geometric limit to:

\[
H_{lat,RMS}/stroke \le 0.25
\]

while preserving the more physically direct lateral-span, rod-angle and transmission limits.

Family 12 then improved significantly, proving that the 0.20 RMS condition had been restrictive.

The continuation eventually encountered other, more physically meaningful boundaries simultaneously:

- H lateral span;
- rod-axis cosine;
- primary transmission sine.

At that point further relaxation was not justified merely to improve thermodynamic performance.

This is a useful general principle:

> Relax proxy constraints when evidence shows they are artificial, but do not casually relax direct mechanical protections.

---

## 17. Paired 30-D thermo-mechanical optimization

After the paired mechanisms existed, the optimization problem changed.

The exact candidate-3952 thermodynamic machine was frozen, and both six-bar mechanisms were varied simultaneously in all **30 continuous geometry dimensions**.

Branches remained fixed within each family.

The objective was:

\[
\max \eta_{indicated}.
\]

The motion RMS relative to 3952 was recorded but **not penalized**.

This prevented the optimizer from being forced to imitate details of 3952 that were not actually thermodynamically necessary.

### 17.1 Fine local trust region

The final stage used deliberately small dimensional adjustments around each established family, for example:

- primary lengths: roughly ±1.5%;
- E coordinates: roughly ±2.5% with a small absolute floor;
- primary phase: approximately ±2°;
- second fixed pivot: approximately ±0.20 crank radius;
- EF/GF: approximately ±3%;
- H ratios: approximately ±0.05;
- piston rod: approximately ±4%;
- slider offset: approximately ±0.20 crank radius;
- slider-axis angle: approximately ±2°.

Infeasible 30-D Sobol directions were not repaired coordinate by coordinate. Their full direction was retained while scalar amplitude was backed off.

This preserved the intended simultaneous S/L motion of the search.

### 17.2 Result after the complete fine / Hlat25 continuation

The final paired mechanism results, still using the original 3952 thermodynamics, were:

| Family | Efficiency | Power | Combined RMS vs 3952 |
|---|---:|---:|---:|
| **1** | **23.1552%** | 41.383 W | 1.292% |
| **50** | **22.5328%** | 41.142 W | 3.133% |
| **4** | **22.2229%** | 40.793 W | 2.434% |
| **12** | **21.2072%** | 40.871 W | 4.777% |

The ranking was stable:

\[
1 > 50 > 4 > 12.
\]

The family-12 result is especially instructive. Its RMS distance from 3952 increased strongly while its efficiency improved by more than one percentage point relative to its initial fixed-thermo mechanism.

This proves that:

> Once a realizable mechanism exists, the thermodynamic optimum within that mechanism family does not necessarily correspond to the closest geometric copy of the abstract target.

---

## 18. Thermodynamic 5-D re-tuning with mechanisms frozen

Each final paired six-bar mechanism was then frozen completely.

Only five machine / exchanger coordinates were released:

\[
V_S/V_L,\quad n_i,\quad L_i,\quad n_o,\quad L_o.
\]

Fixed quantities included:

- ΔT = 260 K;
- total swept volume;
- cylinder clearance ratios;
- total gas inventory;
- exchanger model;
- valve-scaling policy;
- exact six-bar geometry and branches.

Each family received **192 non-seed evaluations** in a persistent incumbent-centered Sobol search.

### 18.1 Final results

| Family | Before thermo 5D | After thermo 5D | Gain | Final power |
|---|---:|---:|---:|---:|
| **1** | 23.1552% | **23.2553%** | +0.1001 pt | 41.434 W |
| **50** | 22.5328% | **22.6516%** | +0.1188 pt | 44.204 W |
| **4** | 22.2229% | **22.4499%** | +0.2270 pt | 42.695 W |
| **12** | 21.2072% | **21.6567%** | +0.4495 pt | 45.347 W |

The ranking again remained unchanged:

\[
1 > 50 > 4 > 12.
\]

The second block from 96 to 192 evaluations produced only small additional gains:

- family 1: about **+0.0146 percentage point**;
- family 50: about **+0.0098 point**;
- family 4: about **+0.0015 point**;
- family 12: **no improvement** beyond its earlier champion.

The thermo-5D search is therefore considered sufficiently converged for this campaign.

### 18.2 Final family-1 comparison with abstract candidate 3952

Candidate 3952:

- efficiency: **23.5081%**;
- power: **40.719 W**.

Final family-1 six-bar machine after thermo re-tuning:

- efficiency: **23.2553%**;
- power: **41.434 W**.

The remaining efficiency gap is only:

\[
23.5081 - 23.2553
\approx
\boxed{0.253\text{ percentage point}}.
\]

The realizable six-bar machine therefore retains approximately **98.9% of the indicated efficiency** of the abstract candidate-3952 law while producing slightly greater indicated power.

This is the strongest validation obtained from the synthesis campaign.

---

## 19. What the four retained families mean

The historical labels **1, 4, 12 and 50** come from earlier ranking positions and should not be treated as permanent mechanism names.

They are useful experimental identifiers, but future documentation should give the mechanism families stable descriptive identities.

The four were retained because they span useful mechanical behavior:

- family 1: highest thermodynamic performance and best kinematic match;
- family 50: mechanically robust, clearly different geometry, second-best thermodynamic result;
- family 4: compact but repeatedly close to transmission and stroke boundaries;
- family 12: strongly different geometry, lower efficiency, but particularly informative because thermo-mechanical adaptation substantially reshaped it.

The four-family panel should therefore be understood as a **diversity panel**, not as the four largest equal-probability primary basins.

The dedicated primary-family document should separately present:

- basin-capture frequency;
- normalized geometry;
- mechanical descriptors;
- downstream transformability.

---

## 20. Primary score, six-bar fit and thermodynamic quality are different rankings

The campaign established three distinct concepts.

### Primary quality

Measured by the cadence proxy, transmission quality and useful E motion.

### Downstream transformability

Measured by how effectively the second dyad can turn the primary motion into the desired piston law.

### Machine quality

Measured by the actual periodic thermodynamic solver using the realizable paired mechanisms.

These rankings are correlated but not identical.

Therefore a synthesis workflow should never prune all diversity based on one early scalar score.

A modest set of mechanically different primaries should survive long enough to receive a fair downstream search.

---

## 21. Why position RMS must become diagnostic rather than objective

During early synthesis, normalized piston-position RMS is extremely useful because it provides a cheap and physically interpretable target.

It should not remain the final objective.

The 30-D thermo-coupled search showed why:

- some families improved thermodynamically while moving farther from 3952;
- the rank ordering remained stable even when RMS distances changed considerably;
- a different realizable mechanism may require a different optimum volume ratio and exchanger sizing.

Thus the correct hierarchy is:

\[
\text{target position}
\rightarrow
\text{synthesis guide}
\]

but

\[
\text{actual thermodynamic efficiency}
\rightarrow
\text{final machine objective}.
\]

The target is a map toward useful geometry, not a constraint that the final machine must obey exactly.

---

## 22. Recommended reusable workflow

For a new DADA machine configuration, the recommended workflow is now:

### Step 1 — obtain a high-quality thermodynamic motion

Use a sufficiently expressive motion representation to identify the useful piston behavior.

Do not assume that all high-order derivative features are physical.

### Step 2 — identify the harder piston

Compare the two piston laws and determine which side carries the sharper cadence changes or more restrictive geometric demand.

Synthesize that side first.

### Step 3 — discover primary families

Run a broad six-dimensional primary search using fresh global islands.

Preserve:

- correct reversals;
- correct branch topology;
- main cadence structure;
- acceptable mechanical transmission.

Do not overfit detailed BP waviness or target acceleration artifacts.

### Step 4 — retain diversity

Select multiple primary families based on both objective quality and mechanical diversity.

### Step 5 — optimize the downstream dyad

Give the selected families comparable nine-dimensional downstream budgets.

### Step 6 — release the complete 15-D mechanism locally

Allow the primary and secondary portions to adapt together.

### Step 7 — mirror to the easier piston

Use the successful harder-side family as a physical mirrored initialization.

Do not impose mirror symmetry afterward.

### Step 8 — locally adapt the second 15-D mechanism

Preserve only justified discrete branches and physical constraints.

### Step 9 — release both mechanisms simultaneously

Run a fine 30-D local thermo-mechanical search with the thermodynamic machine fixed.

Optimize actual indicated efficiency.

### Step 10 — re-tune the thermodynamic hardware

Freeze both mechanisms and release a small, interpretable thermo / exchanger subset.

### Step 11 — perform saturation separately

If the purpose is methodological confidence rather than direct design exploitation, run fresh islands only and measure family capture statistics.

This last step need not be repeated at full 512-island depth for every machine. Once the method is mature, smaller fresh-island checks can be used unless there is evidence that the family structure has changed.

---

## 23. Generalization to other temperatures, gases and applications

The long-term goal is not to restart an unconstrained global mechanism search from zero for every new operating point.

Known mechanism families should become reusable initial structures.

A future machine search should therefore:

1. optimize or predict the thermodynamic target;
2. compare it with known family signatures;
3. warm-start from relevant primary / six-bar families;
4. include a moderate fresh-island component to detect topology changes;
5. deepen global synthesis only when the new target genuinely escapes the known family map.

This is especially relevant for:

- helium engines at higher operating frequency;
- human-powered refrigeration;
- higher-power-density machines;
- temperature sweeps where the optimum motion changes continuously.

For refrigeration, the difficult piston may not be the same side as in the 3952 motor. The “hard side first” rule should be retained, but the identity of that side must be rediscovered.

---

## 24. Comparison with simpler mechanism classes

The six-bar campaign should not be interpreted as proof that six bars are always necessary.

For future applications, especially the low-tech human-powered cooling cell, the same thermodynamic problem should also be optimized with simpler mechanism classes:

1. ideal / free thermodynamic motion;
2. paired six-bar mechanisms;
3. paired four-bar mechanisms;
4. simple slider-crank mechanisms.

Each class should be optimized on its own thermodynamic performance rather than merely fitted to the six-bar solution.

This will quantify the real value of mechanical complexity.

A four-bar or slider-crank may be preferred even with lower ideal performance if it offers sufficient gains in:

- manufacturability;
- friction;
- robustness;
- tolerance to dimensional error;
- mass;
- cost;
- maintenance.

---

## 25. Robustness should precede prototype selection

The present search mainly optimizes nominal geometry.

Before selecting a demonstrator mechanism for fabrication, a separate robustness analysis should perturb:

- link lengths;
- fixed-pivot positions;
- crank phase / indexing;
- slider-axis position and angle;
- plausible bearing and assembly tolerances.

The useful outputs are not only mean efficiency but:

- efficiency dispersion;
- power dispersion;
- probability of crossing a closure or transmission limit;
- sensitivity of the critical clearances;
- required manufacturing tolerances.

A mechanism slightly below the nominal champion may be the better demonstrator if its performance is much less sensitive to manufacturing error.

---

## 26. Documentation strategy

The method, the mechanism families and a specific machine should be documented separately.

### This document

`MECHANISM_SYNTHESIS_SEARCH_THEORY.md`

Purpose:

- explain the reusable search architecture;
- record why the proxies were chosen;
- document saturation evidence;
- separate methodological conclusions from one machine geometry.

### Primary-family catalogue

A separate document should present the reusable primary four-bar families with:

- normalized geometry;
- assembly branch;
- cadence characteristics;
- transmission metrics;
- basin-capture statistics;
- visual / qualitative signature.

### Six-bar-family catalogue

A second separate document should present the paired six-bar families independent of a specific physical scale:

- normalized S and L geometry;
- branch choices;
- stroke/crank;
- transmission metrics;
- H lateral metrics;
- rod-angle metrics;
- normalized piston laws;
- thermodynamic performance in the 3952 reference machine.

### Demonstrator / wiki machine

A selected demonstrator should then be presented as one **scaled instance** of a mechanism family, with:

- actual crank radius;
- all link lengths;
- cylinder bores and strokes;
- exchanger dimensions;
- gas charge;
- operating conditions;
- animated mechanism geometry.

This separation prevents one prototype drawing from becoming confused with the general synthesis method.

---

## 27. Repository artifacts for the completed 3952 campaign

The principal artifacts are:

### Primary search and saturation

- `examples/search_large_primary_cadence_3952_v4_1.py`
- `examples/search_large_primary_cadence_3952_v4_1_saturation.py`
- `outputs/large_primary_cadence_3952_v4_1.json`
- `outputs/large_primary_cadence_3952_v4_1_saturation.json`

### Downstream and six-bar synthesis

- `examples/search_large_second_dyad_v41_panel_3952.py`
- `outputs/large_second_dyad_v41_panel_3952.json`
- `examples/polish_large_sixbar_v41_families_3952.py`
- `outputs/large_sixbar_v41_family_polish_3952.json`

### Mirrored SMALL synthesis

- `examples/mirror_polish_small_from_large_v41_3952.py`
- `outputs/sixbar_pairs_3952/mirror_family_summary.json`
- paired mechanism JSON files under `outputs/sixbar_pairs_3952/`

### Thermo-coupled mechanism optimization

- fine simultaneous S/L mechanism search outputs under:
  `outputs/sixbar_thermo_coupled_3952_fine/`
- relaxed-H-lateral continuation under:
  `outputs/sixbar_thermo_coupled_3952_fine_hlat25/`

### Thermodynamic re-tuning

- `examples/optimize_sixbar_pairs_thermo5d_3952_hlat25.py`
- `outputs/sixbar_thermo5d_3952/`

These outputs provide the reproducible numerical record behind the methodological conclusions in this document.

---

## 28. Stopping logic validated by the campaign

The primary search can now be considered sufficiently saturated for engineering use because all of the intended stopping signals were observed:

1. the best score improved only marginally as the island count approached 512;
2. the dominant capture fractions stabilized;
3. the top four families accounted for roughly 97–98% of capture probability depending on nearby clustering threshold;
4. the Good-Turing unseen capture-mass estimate fell to approximately 1–1.6%;
5. additional rare families did not alter the competitive best score or the set of useful downstream concepts;
6. the qualitative conclusion was stable across family-distance thresholds 0.03–0.05.

This is strong saturation evidence, not a proof of mathematical global optimality.

The proper wording remains:

> The dominant primary basins are well mapped under the stated V4.1 search policy and budget.

It should not be stated that every physically possible four-bar family has been enumerated.

---

## 29. Current methodological principles

The completed campaign supports the following rules.

1. **Do not fit target artifacts more accurately than their physical meaning warrants.**
2. **Preserve motion topology before fitting detailed shape.**
3. **Use position as the primary mechanical truth during synthesis.**
4. **Use derivatives as guarded descriptors, not unquestioned targets.**
5. **Do not convert a convenient smoothness assumption into a physical requirement.**
6. **Search the primary four-bar globally before releasing the full six-bar geometry.**
7. **Do not freeze the primary permanently; release all dimensions locally after a viable downstream dyad exists.**
8. **Select the kinematically harder piston first, not automatically SMALL or LARGE.**
9. **Use physical mirroring as an initialization, not as a final symmetry constraint.**
10. **Retain mechanically diverse families even when one scalar score is slightly worse.**
11. **Treat primary capture probability, downstream transformability and final thermo performance as different quantities.**
12. **Use fresh islands to measure saturation; use seeds to exploit known families.**
13. **Interpret repeated convergence as basin-capture evidence, not literal geometric volume.**
14. **Relax proxy constraints only when evidence shows they are artificial.**
15. **Keep direct transmission, closure and rod-angle protections unless there is a physical reason to change them.**
16. **Once paired mechanisms exist, optimize both simultaneously rather than alternately freezing one side.**
17. **At the thermo-coupled stage, make target-motion RMS diagnostic only.**
18. **Let full periodic thermodynamic replay select among realizable mechanisms.**
19. **After mechanism adaptation, allow the thermodynamic hardware to re-tune around the new motion.**
20. **Document normalized mechanism families separately from specific scaled machines.**

---

## 30. Final conclusion from candidate 3952

The 3952 campaign achieved its original practical objective.

A difficult abstract C2 piston law was converted into a mechanically realizable paired six-bar architecture without requiring the final mechanism to reproduce every artifact of the target.

The dominant primary search basins were mapped with 512 fresh global islands.

Four mechanically diverse downstream families were carried through:

- downstream synthesis;
- 15-D local polish;
- physical mirroring;
- opposite-side adaptation;
- paired 30-D thermo-mechanical optimization;
- final 5-D thermodynamic re-tuning.

The best realizable family reached:

- **23.2553% indicated thermal efficiency**;
- **41.434 W indicated power**;

compared with abstract candidate 3952 at:

- **23.5081%**;
- **40.719 W**.

The remaining efficiency difference is approximately **0.253 percentage point**.

The central result is therefore not only one successful six-bar mechanism.

It is a reusable synthesis method that connects:

\[
\text{thermodynamic target}
\rightarrow
\text{primary family discovery}
\rightarrow
\text{downstream mechanism synthesis}
\rightarrow
\text{paired mechanical adaptation}
\rightarrow
\text{thermodynamic validation}.
\]

That method can now be applied to simpler mechanism classes, helium engines, refrigeration cycles and other DADA machine configurations without restarting the conceptual search from zero.
