# Mechanism-synthesis search theory for the DADA engine

## Purpose

This note records the search methodology developed while trying to transform an optimized thermodynamic piston law into a mechanically realizable four-bar / six-bar linkage.

The immediate reference case is the large-cylinder motion of structured C2 candidate 3952, but the purpose is broader: define a search architecture that can later be reused for other temperatures, machine sizes, working gases and operating points without treating every mechanism search as an unrelated global-optimization problem.

The central conclusion is that mechanism synthesis should be treated as a **hierarchical basin-discovery problem**, not as a brute-force fit of all linkage dimensions simultaneously.

---

## 1. Thermodynamic target and mechanical mechanism are not the same object

The optimized thermodynamic motion is a design probe. It is not automatically a mechanically meaningful trajectory.

The realizable linkage should preserve the thermodynamic effects that matter while remaining free to smooth, redistribute or reshape features that are artifacts of the chosen target parameterization.

This distinction became particularly important for candidate 3952. Its large-cylinder target contains a severe high-pressure (HP) acceleration feature near the HP transition. Automatic detection identified the unreliable region approximately as:

- start: **34.0°**;
- end: **44.25°**;
- dominant jerk fronts: approximately 34.25°, 36.5°, 38.75°, 39.5°, 41.75°, 43.75°.

The position law remains meaningful through this interval. The extreme acceleration shape does not.

Therefore:

1. **position remains the primary kinematic truth**;
2. velocity is useful as a cadence descriptor and weak tie-breaker;
3. target acceleration should not be directly fitted in the artificial HP feature;
4. the realizable linkage is allowed, and expected, to round the HP transition;
5. final thermodynamic replay of the actual linkage remains the decisive test.

A previous downstream search demonstrated why this matters: completely removing the HP interval from the position objective allowed the optimizer to insert an unphysical pause and slight reversal while still obtaining an apparently good score. A noisy acceleration interval is not permission to ignore position.

---

## 2. Why the blind 15-dimensional six-bar search was inefficient

The complete six-bar parameterization contains 15 continuous variables.

### Primary four-bar and coupler point E

1. primary ground length `AD / AB`;
2. primary coupler length `BC / AB`;
3. primary rocker length `CD / AB`;
4. `E_along`;
5. `E_normal`;
6. primary crank phase.

### Secondary dyad, point H and piston slider

7. second fixed-pivot x;
8. second fixed-pivot y;
9. link `EF`;
10. link `GF`;
11. `H_along / EF`;
12. `H_normal / EF`;
13. piston rod length;
14. slider-axis offset;
15. slider-axis angle.

There are also discrete primary and secondary assembly branches.

A long blind 15-D search consumed roughly twelve million evaluations and did not improve the existing large-cylinder six-bar fit. Many islands terminated immediately because invalid initial populations saw an almost constant penalty. The best mechanically attractive new solution was much worse kinematically than the existing solution.

The failure was useful evidence:

- the feasible region is thin relative to the raw 15-D box;
- closure and transmission constraints create a difficult landscape;
- flat invalid penalties destroy differential-evolution guidance;
- a monolithic optimizer spends too much effort rediscovering primary-linkage structure that can be studied much more cheaply on its own.

This motivated a hierarchical architecture.

---

## 3. Hierarchical six-bar synthesis

The preferred decomposition is:

### Stage A — primary four-bar library

Search only the six primary variables

\[
x_p=(g,c,r,E_\parallel,E_\perp,\phi).
\]

This search is cheap enough to run hundreds of independent global islands.

Its purpose is **not** to reproduce the final piston motion exactly. The primary should provide a useful temporal cadence and a mechanically attractive carrier motion for point E.

### Stage B — downstream synthesis

For each selected primary, optimize the remaining nine continuous variables

\[
(G_x,G_y,EF,GF,H_\parallel/EF,H_\perp/EF,L_{rod},s_{slider},\alpha_{slider})
\]

plus the secondary branch.

The downstream dyad is allowed to reshape the primary cadence substantially.

### Stage C — optional full local polish

Only after a successful primary + dyad combination exists should all 15 continuous variables be released together for local refinement.

### Stage D — thermodynamic replay

The actual linkage trajectory is replayed through the thermodynamic solver. Indicated efficiency, power, phase pressure behavior, exchanger behavior and all validity constraints decide whether the mechanism is actually useful.

This ordering prevents the kinematic proxy from becoming the final objective.

---

## 4. Primary four-bar feasible space

The present normalized primary search uses `AB = 1` and six continuous parameters:

\[
(g,c,r,E_\parallel,E_\perp,\phi).
\]

Current broad bounds are:

| Parameter | Bound |
|---|---:|
| `AD / AB` | 0.50 .. 12.0 |
| `BC / AB` | 0.50 .. 12.0 |
| `CD / AB` | 0.35 .. 6.0 |
| `E_along` | -12.0 .. 18.0 |
| `E_normal` | -15.0 .. 15.0 |
| phase | -π .. +π |

Both primary assembly branches are admissible.

Full-revolution closure is checked. Primary transmission quality is measured by

\[
s_{min}=\min_\theta\frac{|BC\times DC|}{|BC||DC|}=\min_\theta|\sin\mu|.
\]

A hard floor of 0.30 has been retained, corresponding to roughly 17.5° from toggle. Later searches added a **soft preference** around 0.35 rather than tightening the hard floor. This distinction matters: a mechanical preference should not silently remove potentially useful kinematic families.

---

## 5. Evolution of the primary-cadence objective

The sequence V1 → V4.1 is important because it shows several ways in which a plausible mathematical proxy can be physically wrong.

### 5.1 V1 — acceleration-lobe timing

V1 used the magnitude of E velocity and its tangential acceleration

\[
a_t=\frac{E'\cdot E''}{\|E'\|}=\frac{d\|E'\|}{d\theta},
\]

then tried to match signed piston-acceleration lobes.

This was conceptually inconsistent: the derivative of a non-negative speed magnitude is not the signed piston acceleration.

V1 nevertheless produced useful geometric diversity and exposed the tendency of the optimizer to exploit the transmission floor.

### 5.2 V2 — signed principal-axis velocity and three acceleration events

V2 projected E motion onto its PCA principal axis and searched three signed acceleration events plus linear velocity between them.

Two problems appeared.

First, the PCA sign and direction are not physically equivalent to the final piston direction, because the downstream dyad can rotate and transform motion.

Second, the objective could almost eliminate one nominally required acceleration event because event strength was only weakly penalized.

The best candidates again tended to sit at the minimum transmission sine.

### 5.3 V3 — three piecewise-linear velocity regimes

V3 replaced isolated acceleration peaks by three automatically discovered piecewise-linear velocity regimes.

This was a methodological improvement: it evaluated broad temporal structure rather than three local peaks.

However, the statistically optimal three-line partition was physically wrong. For candidate 3952 it placed boundaries near 349°, 31° and 181°. One fitted regime crossed the real piston turnaround near 154°, mixing motion in opposite directions inside one nominal phase.

This was a critical lesson:

> An automatic segmentation can minimize approximation error while violating the topology of the physical motion.

### 5.4 V4 — physical monotonic branches

V4 starts with the two real velocity zero crossings and therefore preserves the motion topology.

For candidate 3952, after bridging and smoothing the unreliable HP velocity feature:

- high-position turnaround: **358.812°**;
- low-position turnaround: **154.144°**;
- short high-to-low branch: **155.332°**;
- long low-to-high branch: **204.668°**.

The short branch naturally divides into two cadences. The automatic two-level speed split lies near **38.062°**.

With guarded core intervals, the target has approximately:

- fast/slow mean-speed ratio: **5.027**;
- fraction of short-branch displacement completed in the fast part: **57.63%**.

V4 therefore rewards:

1. two correctly timed turnarounds;
2. no additional reversal;
3. monotonicity of both branches;
4. target-like fast/slow cadence ratio on the short branch;
5. target-like displacement sharing on the short branch;
6. acceptable transmission and weak directionality preference.

Its first version also compared the long branch against the exact normalized 3952 velocity profile.

### 5.5 V4.1 — BP symmetry instead of BP profile matching

Earlier Fourier motion searches showed that a visibly wavy low-pressure (BP) exchange can still give good thermodynamic efficiency. A perfectly regular individual piston velocity is therefore not the physical objective.

During BP exchange the more relevant requirement is coordinated cylinder-volume change that keeps pressure approximately constant.

For that reason V4.1 removes the requirement that the long branch resemble the 3952 velocity profile.

The long branch may be:

- approximately sinusoidal;
- `fast → slow → fast`;
- multiply wavy;
- otherwise strongly non-uniform.

No smoothness, curvature or local-extrema-count penalty is applied.

The remaining proxy is approximate mirror symmetry of the candidate's own long branch:

\[
A_{BP}=\frac{\operatorname{RMS}[v(u)-v(1-u)]}
{\sqrt{\operatorname{mean}((v(u)^2+v(1-u)^2)/2)}}.
\]

The comparison is performed between the candidate's own two matched turnarounds, not at the exact target angles. Turnaround timing and BP symmetry are therefore separate criteria.

A perfectly symmetric `fast → slow → fast` long branch can obtain \(A_{BP}=0\).

Ultimately even this symmetry is only a search proxy. The final criterion is the pressure ripple obtained when both actual piston motions are replayed together through the thermodynamic solver.

---

## 6. Differential-evolution islands do not have a simple geometric “volume”

A common intuitive picture is that each optimization island explores one local piece of the six-dimensional box. That is not how the present differential evolution (DE) search works.

With six variables and `popsize = 16`, an island begins with approximately

\[
6\times16=96
\]

Latin-hypercube individuals distributed across the **entire** parameter box.

The island then evolves that population by mutation, crossover and selection. Its trajectory is global-to-local rather than a fixed small hypervolume.

Therefore two independent islands ending in the same mechanism family do not necessarily mean that two predefined local search volumes overlapped.

The useful quantity is instead the **basin capture probability**.

---

## 7. Geometric basin volume versus algorithmic capture probability

Two different quantities must be distinguished.

### 7.1 Geometric basin volume

Let \(V_i\) denote the literal measure of the region of feasible parameter space associated with a mechanism family or attraction basin.

This quantity depends on the chosen coordinate system and on what is considered feasible. It is difficult to estimate directly in six dimensions and is not what DE samples uniformly after the first generation.

### 7.2 Algorithmic capture probability

Define \(P_i^{capture}\) as the probability that one fresh DE island, under a fixed algorithm, population, bounds, constraints and budget, terminates in family \(i\).

This is an operational property of the **search method + objective landscape**.

For practical optimization it is often more useful than raw geometric volume. A geometrically narrow but strongly attractive valley may have a large capture probability, while a large flat region may have little practical attraction.

The capture probability is estimated by independent fresh-island frequency:

\[
\hat P_i=\frac{n_i}{N}.
\]

Historical seeds must not be used when estimating this quantity because they deliberately bias the initial condition toward known families.

---

## 8. Empirical evidence from the V4 campaign

The V4 campaign used 128 islands, of which 112 were fresh and 16 were seeded.

Considering only the 112 fresh islands and the current operational family distance, the four most frequent families captured approximately:

- **24.1%** of fresh islands;
- **20.5%**;
- **17.0%**;
- **9.8%**.

Together these four families captured about **71%** of the independent searches. The ten most frequent families captured roughly **90%**.

This is strong evidence that the major basins were not accidental discoveries. They are repeatedly recovered from independent global initial populations.

However, the same 112 fresh islands produced 21 observed winning families, of which:

- 11 were seen only once;
- 2 were seen twice.

Therefore the large basins appear well recognized while the long tail of narrow basins is clearly not saturated.

This is precisely why “the best family was found many times” and “the mechanism space is exhaustively explored” are not equivalent statements.

---

## 9. Operational definition of a mechanism family

A family is currently defined by normalized geometric distance in the six primary variables, with opposite assembly branches kept separate.

For ordinary library diversification, a representative distance threshold of approximately 0.04 has been used.

Any family definition introduces a scale. Saturation conclusions should therefore not depend on one arbitrary threshold.

The dedicated saturation experiment reports results at several nearby thresholds:

- 0.03;
- 0.04;
- 0.05.

If the qualitative saturation conclusion changes radically across this modest range, the family metric itself requires refinement.

Connected-component clustering is useful for saturation statistics because it is independent of candidate insertion order. The possibility of single-link chaining must nevertheless be monitored.

---

## 10. Saturation experiment

A proper saturation study should use only fresh islands and accumulate them in blocks, for example:

\[
64\rightarrow128\rightarrow192\rightarrow256\rightarrow\dots\rightarrow512.
\]

At every checkpoint record:

1. best objective value;
2. number of observed families;
3. family capture-frequency distribution;
4. number of singleton families \(f_1\);
5. number of doubleton families \(f_2\);
6. fraction of islands captured by the largest 1, 4 and 10 families;
7. evolution of the mechanical / kinematic Pareto front;
8. sensitivity to family-distance threshold.

A 512-island run with the present V4/V4.1 cost is suitable as an overnight experiment on the development machine.

---

## 11. Estimating undiscovered search mass

Two ecological-style estimators are useful as diagnostics. They are not proofs of global optimality.

### 11.1 Good–Turing unseen capture mass

A simple estimate of probability mass belonging to unseen families is

\[
\hat p_0\approx\frac{f_1}{N},
\]

where \(f_1\) is the number of families seen exactly once among \(N\) islands.

Interpretation:

- a high singleton fraction means appreciable capture probability is still entering newly discovered families;
- a low singleton fraction means another random island is increasingly likely to rediscover an already known family.

This estimates **unseen optimizer-outcome probability**, not unseen geometric volume.

### 11.2 Chao1 family-count estimate

A useful rough lower-bound-style estimator is

\[
\hat S_{Chao1}=S_{obs}+\frac{f_1^2}{2f_2},
\]

when \(f_2>0\).

If \(f_2=0\), a finite fallback form can be used as a diagnostic.

Again, the result depends on the operational family definition and on the DE policy. It should be used to compare saturation across checkpoints, not quoted as the literal number of physically existing four-bar families.

---

## 12. Seeds have a different purpose from fresh islands

Historical seeds are valuable in a production optimization campaign. They preserve expensive discoveries, accelerate convergence, allow a changed objective to reuse earlier geometry, and reduce the chance of losing a narrow but already known good basin.

They are inappropriate for measuring capture probabilities because they alter the starting distribution.

Therefore two campaign types should remain separate.

### Exploitation / design campaign

Use a mixture of:

- historical seeds;
- diverse known families;
- fresh Latin-hypercube islands.

Goal: find the best mechanism efficiently.

### Saturation / methodology campaign

Use:

- fresh islands only;
- fixed complete bounds;
- fixed DE settings;
- deterministic independent seeds;
- balanced treatment of discrete branches.

Goal: estimate discovery saturation and capture probabilities.

Mixing these two purposes makes both conclusions weaker.

---

## 13. Proposed stopping logic

No single statistic proves global optimality. A practical stopping rule should combine several observations.

A primary search can be considered sufficiently saturated for engineering use when successive large additions of fresh islands show all of the following:

1. the best objective no longer improves materially;
2. the useful mechanical/kinematic Pareto front is stable;
3. the dominant family capture fractions are stable;
4. the number of new **competitive** families per added block becomes very small;
5. Good–Turing unseen capture mass becomes small and continues to decrease;
6. the conclusion is stable for nearby family-distance thresholds.

“Competitive” is important. It is unnecessary to enumerate every poor local minimum if new islands no longer alter the set of mechanisms worth sending to the expensive downstream dyad search.

The final stopping thresholds should be calibrated from actual campaigns rather than chosen in advance solely for mathematical neatness.

---

## 14. Mechanical diversity must be retained explicitly

A scalar cadence score is not enough.

Two primaries with almost equal kinematic score may differ greatly in:

- minimum transmission angle;
- link-length ratios;
- distance of E from the coupler;
- trajectory directionality;
- overall packaging;
- downstream closure opportunities.

The primary library should therefore retain diverse families and expose mechanical descriptors such as:

- minimum primary transmission sine;
- `BE / BC`;
- `AD`, `BC`, `CD`;
- PCA directionality;
- E lateral span;
- E path length.

A useful downstream panel should deliberately contain different trade-offs, not only the first rows of one scalar ranking.

---

## 15. Why the downstream dyad remains essential

The primary four-bar is only an intermediate motion generator.

Evidence from the historical six-bar solution is important here: a primary that scores poorly under some simplified cadence metrics can still participate in a six-bar giving a much better final piston-position fit.

Therefore the primary objective should identify **promising carrier motions**, not pretend to know the complete downstream transform.

The decisive experiment is to take several physically distinct primaries and give each a fair nine-variable downstream-dyad optimization budget.

A useful panel should contain:

- the best compact fast/slow candidate;
- the scalar-score champion even if E is less compact;
- a candidate emphasizing BP symmetry;
- another mechanically plausible trade-off;
- the historical primary from the known ~2.08%-RMS large six-bar as a control.

If the new primary families consistently allow better downstream fits than the historical control, the hierarchical primary metric is validated.

If the control remains superior, the discrepancy identifies information missing from the primary proxy.

---

## 16. Final thermodynamic selection

Neither primary score nor final piston-position RMS is the end objective.

A realizable six-bar should be replayed through the thermodynamic model and evaluated on at least:

- indicated efficiency;
- indicated power;
- HP and BP pressure ripple;
- heat input and rejection;
- exchanger flow / Mach / Reynolds validity;
- valve behavior;
- cylinder extrema and clearance;
- mechanical transmission constraints.

For the BP exchange in particular, a wavy but symmetric individual piston law may be entirely acceptable if the **combined volume evolution of both cylinders** keeps pressure close to constant.

Thus

\[
\text{primary symmetry}\rightarrow\text{useful search proxy},
\]

whereas

\[
\text{actual BP pressure history}\rightarrow\text{physical validation}.
\]

---

## 17. Generalization to other machine configurations

The longer-term objective is not to rerun an unconstrained global mechanism search from zero for every temperature or machine configuration.

A reusable workflow should be:

1. optimize or predict the thermodynamic target;
2. identify robust physical motion features rather than exact high-order derivatives;
3. deform or warm-start from known mechanism families;
4. run a moderate fresh-island component to detect topology changes;
5. preserve family identities across operating conditions;
6. optimize downstream dyads;
7. replay thermodynamics;
8. deepen the search only when family saturation or physical replay justifies it.

This should eventually allow mechanism families themselves to be mapped as a function of source temperatures, volume ratio, gas and machine scale.

---

## 18. Current methodological principles

The campaign supports the following working rules:

1. **Do not fit target artifacts more accurately than their physical meaning warrants.**
2. **Preserve motion topology before fitting detailed shape.**
3. **Use position as the primary mechanical truth; use derivatives as guarded diagnostics.**
4. **Do not turn a convenient smoothness assumption into a physical requirement.**
5. **Separate primary cadence generation from downstream motion synthesis.**
6. **Use soft mechanical preferences when the true acceptable boundary is not known.**
7. **Use fresh islands to measure discovery saturation; use seeds to build the best machine.**
8. **Interpret repeated convergence as basin-capture evidence, not literal percentage of geometric volume.**
9. **Retain diverse mechanical families even when their scalar scores are slightly worse.**
10. **Let full thermodynamic replay decide among good realizable mechanisms.**

These principles are more valuable than the numerical optimum of any single campaign because they define a search method that can survive changes in the thermodynamic target and machine configuration.
