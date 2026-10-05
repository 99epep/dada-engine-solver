# Mechanism-synthesis search for the DADA engine

## Purpose

This note records the reusable mechanism-synthesis method developed for
transforming an optimized thermodynamic piston motion into mechanically
realizable four-bar and six-bar linkages.

The reference study used a structured-C2 motor motion at a 260 K
source-temperature difference, but the method is intended to apply to other
temperatures, machine sizes, working fluids and operating points.

The main methodological result is:

> Mechanism synthesis should be treated as a hierarchical basin-discovery and
> thermo-mechanical adaptation problem, not as a blind fit of all linkage
> dimensions at once.

The practical sequence is:

1. identify the piston carrying the most restrictive kinematic demand;
2. discover a diverse set of primary four-bar families;
3. give selected primaries a fair downstream-dyad search;
4. locally release the complete six-bar geometry;
5. mirror a successful mechanism toward the opposite piston as an
   initialization only;
6. locally adapt the mirrored mechanism;
7. release both mechanisms simultaneously and optimize them using the real
   thermodynamic objective;
8. with the mechanisms fixed, re-tune the small set of relevant thermodynamic
   and exchanger dimensions;
9. use a separate fresh-island saturation experiment when the coverage of the
   primary search space itself must be assessed.

The retained primary seed geometries are documented separately in
[`PRIMARY_FOUR_BAR_FAMILIES.md`](PRIMARY_FOUR_BAR_FAMILIES.md), and complete
six-bar family comparisons belong in
[`SIX_BAR_MECHANISM_FAMILIES.md`](SIX_BAR_MECHANISM_FAMILIES.md).

---

## 1. Thermodynamic target and mechanical mechanism are different objects

An optimized thermodynamic motion is a design probe. It is not automatically a
mechanically meaningful trajectory.

A realizable linkage should preserve the thermodynamic effects that matter
while remaining free to smooth, redistribute or reshape features introduced by
the target parameterization itself.

The working rules are:

1. **position is the primary kinematic reference**;
2. velocity is useful for motion topology and cadence;
3. sharp acceleration features should not automatically be fitted as physical
   requirements;
4. noisy derivative information is not permission to discard the underlying
   position law;
5. once a viable mechanism exists, its final quality must be assessed by the
   thermodynamic solver.

This distinction is important because an optimizer can exploit an omitted
position interval by introducing a pause or reversal while still appearing to
match derivative-based targets well.

The target therefore guides discovery, but the final mechanism is not required
to reproduce every feature of the abstract trajectory.

---

## 2. Six-bar topology and variables

The retained six-bar topology consists of:

- primary four-bar `A-B-C-D`;
- point `E` rigidly attached to coupler `BC`;
- secondary dyad `E-F-G`;
- point `H` rigidly attached to link `EF`;
- a finite piston rod connecting `H` to a slider axis.

The crank length is normalized:

\[
AB = 1.
\]

Each cylinder mechanism has 15 continuous coordinates.

### Primary four-bar

\[
(AD,\ BC,\ CD,\ E_\parallel,\ E_\perp,\ \phi)
\]

where the first three lengths and the point-E coordinates are expressed in
crank-radius units.

### Downstream mechanism

\[
(G_x,\ G_y,\ EF,\ GF,\ H_\parallel/EF,\ H_\perp/EF,
L_{rod},\ s_{slider},\ \alpha_{slider})
\]

There are also two discrete assembly branches:

- primary four-bar branch;
- secondary dyad branch.

A paired two-cylinder machine therefore contains 30 continuous mechanism
coordinates once the discrete branches have been chosen.

The production representation of this topology is implemented by
`SixBarCylinderMechanism`; the search-stage coordinate groups are exposed by
`dada_solver.research.synthesis`.

---

## 3. Why a blind 15-dimensional search is inefficient

A direct global search over all 15 continuous dimensions has several structural
disadvantages:

- the feasible region is thin relative to the raw parameter box;
- loop closure, transmission and slider constraints create large invalid
  regions;
- flat invalid penalties provide little information to a global optimizer;
- useful structure exists already in the primary four-bar;
- the downstream dyad can substantially transform the primary motion;
- randomizing primary and downstream geometry simultaneously spends a large
  fraction of the search budget rediscovering basic feasibility.

The useful decomposition is therefore hierarchical.

This is not merely a software optimization. It reflects the physical roles of
the two stages: the primary produces a temporal and planar motion structure,
while the downstream dyad acts as a genuine motion transformer.

---

## 4. Hierarchical synthesis workflow

The generic stage names used by `dada_solver.research.synthesis` are retained
below.

### 4.1 `primary_discovery`

Search only the six primary coordinates:

\[
x_p=(g,c,r,E_\parallel,E_\perp,\phi).
\]

The goal is not a final piston mechanism.

The primary should provide a useful temporal cadence and a point-E trajectory
that a downstream dyad can exploit.

Because this problem has only six continuous dimensions, many independent
global search islands are practical.

Several mechanically distinct families should be retained rather than only the
lowest-scoring primary.

### 4.2 `downstream_fit`

For each selected primary, freeze its six coordinates and search the remaining
nine continuous coordinates plus the secondary assembly branch.

The downstream dyad is deliberately allowed to reshape the primary motion.

A primary-only score is therefore not a reliable predictor of downstream
transformability.

### 4.3 `full_local_polish`

Once a viable primary and downstream dyad exist, release all 15 continuous
coordinates in a moderate local region.

This stage is important because modest primary changes can allow large
improvements in the complete piston motion that are impossible with the
primary frozen.

### 4.4 `mirror_initialization`

The mechanism synthesized for the more restrictive piston may be spatially
mirrored and phase-aligned toward the opposite piston.

The mirror is an **initialization only**.

It must not become a permanent symmetry constraint unless a separate physical
argument justifies that restriction.

### 4.5 `opposite_local_adaptation`

The mirrored mechanism receives its own independent local adaptation while its
chosen discrete branches are retained.

The two final cylinder mechanisms are therefore allowed to become different.

### 4.6 `paired_thermodynamic`

Once realizable mechanisms exist on both sides, release the two complete
mechanisms simultaneously.

At this point the optimization objective changes fundamentally:

- the actual thermodynamic objective becomes primary;
- position error relative to the abstract target becomes diagnostic;
- mechanism-family identity is preserved through a local search region and
  fixed discrete branches rather than through exact geometric imitation.

The abstract target has served its purpose once it has led the search into a
useful realizable basin.

### 4.7 `hardware_retuning`

Finally, freeze the chosen mechanisms and re-tune the small set of
thermodynamic or exchanger coordinates that had previously been optimized for
the abstract piston law.

This separates:

1. performance lost because of mechanical realization;
2. performance recoverable by adapting the surrounding machine to the new
   mechanism.

---

## 5. Primary four-bar search space

The reference primary search normalizes the crank to:

\[
AB=1.
\]

It searches:

\[
(g,c,r,E_\parallel,E_\perp,\phi).
\]

The retained broad bounds are:

| Coordinate | Bound |
|---|---:|
| `AD / AB` | 0.50 .. 12.0 |
| `BC / AB` | 0.50 .. 12.0 |
| `CD / AB` | 0.35 .. 6.0 |
| `E_along` | -12.0 .. 18.0 |
| `E_normal` | -15.0 .. 15.0 |
| phase | -π .. +π |

Both primary assembly branches are admissible.

Full-revolution closure is mandatory.

Primary transmission quality is measured by

\[
s_{\min}
=
\min_\theta
\frac{|BC\times DC|}{|BC||DC|}
=
\min_\theta |\sin\mu|.
\]

The reference search uses a hard floor:

\[
s_{\min}\ge 0.30
\]

and only a soft preference toward approximately `0.35`.

This distinction is deliberate:

> A mechanical preference should not silently become a topological exclusion
> unless a real physical limit justifies it.

The retained search proxy also requires a minimum point-E span of `0.75` crank
radii.

The exact documentary policy is stored in
[`repro/mechanism_synthesis_search/policy.toml`](repro/mechanism_synthesis_search/policy.toml).

---

## 6. Lessons from rejected primary-search proxies

Several successive proxy designs were useful because they exposed ways in which
a mathematically plausible objective can be mechanically misleading.

Their historical version names are not important. The conceptual lessons are.

### 6.1 Acceleration-lobe timing is not enough

An early proxy used the magnitude of point-E velocity and its tangential
acceleration,

\[
a_t
=
\frac{E'\cdot E''}{\|E'\|}
=
\frac{d\|E'\|}{d\theta},
\]

then tried to relate those events to signed piston acceleration.

This is not physically equivalent: the derivative of a non-negative speed
magnitude is not signed piston acceleration.

The experiment nevertheless exposed a recurring optimizer tendency to exploit
the transmission floor.

### 6.2 Signed PCA velocity improves sign information but not physical identity

A later proxy projected point-E motion onto its principal PCA axis and used
signed projected velocity.

This repaired the sign problem but revealed two other limitations:

- the PCA axis is not the final piston axis because the downstream dyad can
  rotate and transform the motion;
- event-based penalties can allow a nominally required event to become
  arbitrarily weak without a correspondingly large penalty.

### 6.3 Statistically optimal segmentation can violate motion topology

Automatic piecewise velocity segmentation was then tested.

It produced a low approximation error but allowed a fitted regime to cross a
real piston reversal.

The resulting lesson is general:

> A segmentation may minimize approximation error while violating the topology
> of the physical motion.

Real turnarounds must therefore structure the cadence objective.

### 6.4 Retained topology-and-cadence proxy

The retained reference proxy begins from the two physical turnarounds.

The reference motion has approximately:

- high-position turnaround: `358.812°`;
- low-position turnaround: `154.144°`;
- short high-to-low branch: `155.332°`;
- long low-to-high branch: `204.668°`.

The short branch contains a fast and a slow cadence.

Reference descriptors are approximately:

- fast/slow mean-speed ratio: `5.027`;
- fast-sector displacement fraction: `0.5763`.

The proxy rewards:

1. two correctly timed turnarounds;
2. no additional reversal;
3. correct monotonic sign on both branches;
4. target-like fast/slow cadence on the short branch;
5. target-like displacement sharing on that branch;
6. acceptable primary transmission;
7. weak geometric directionality preferences.

The long return branch is deliberately less constrained.

It is not required to copy the detailed target velocity profile and receives no
smoothness, curvature or local-extrema penalty.

Instead, the retained search proxy uses approximate mirror symmetry of the
candidate's own long branch:

\[
A_{BP}
=
\frac{
\operatorname{RMS}[v(u)-v(1-u)]
}{
\sqrt{
\operatorname{mean}
\left[
\frac{v(u)^2+v(1-u)^2}{2}
\right]
}
}.
\]

A symmetric fast-slow-fast branch is therefore allowed.

This is still only a search proxy.

The physical question is ultimately whether the complete cylinder-volume
evolution produces good thermodynamic behavior.

---

## 7. Reference primary score

For documentary reproducibility, the retained proxy can be summarized as

\[
\begin{aligned}
S ={}&
2.00\,M
+0.40\,Z
+0.65\,\frac{T}{15^\circ}
+0.50\,N_{\mathrm{extra}} \\
&+0.55\,\left|\ln\frac{R}{R_*}\right|
+0.80\,|F-F_*|
+0.45\,A_{BP}
+P_s
+P_a .
\end{aligned}
\]

where:

- `M` is normalized RMS wrong-sign velocity over the two monotonic branch
  cores;
- `Z` is normalized RMS projected velocity at the target turnarounds;
- `T` is RMS angular error of the two matched zero crossings;
- `N_extra` counts additional zero crossings;
- `R` is the candidate fast/slow mean-speed ratio;
- `R*` is the reference ratio;
- `F` is candidate fast-sector displacement fraction;
- `F*` is the reference fraction;
- `A_BP` is long-branch mirror asymmetry;
- `P_s` is a soft primary-transmission penalty;
- `P_a` is a weak PCA-axis preference.

For the transmission preference,

\[
P_s=0
\]

when

\[
s_{\min}\ge0.35,
\]

and otherwise

\[
P_s
=
0.12
\frac{0.35-s_{\min}}{0.35-0.30}.
\]

Candidates below the hard `0.30` transmission floor are rejected rather than
merely penalized by this soft term.

For PCA directionality,

\[
P_a=0
\]

when the principal-axis position-variance fraction is at least `0.70`;
otherwise

\[
P_a
=
0.03
\frac{0.70-a}{0.70}.
\]

The two possible signs of the PCA axis are evaluated and the lower score is
retained.

The score should not be interpreted as:

- thermodynamic efficiency;
- complete-mechanism quality;
- manufacturing robustness;
- geometric probability density.

It exists only to make primary family discovery computationally tractable.

---

## 8. Search islands and basin-capture probability

A differential-evolution island is not one small geometric patch.

With six variables and population size 16, one fresh island begins with

\[
6\times16=96
\]

Latin-hypercube individuals distributed across the complete bounded domain.

The population then evolves globally to locally.

The useful statistical quantity is therefore **algorithmic basin-capture
probability**:

\[
P_i^{capture}
=
P(\text{one fresh search island terminates in family }i).
\]

For fixed:

- bounds;
- constraints;
- optimizer;
- population policy;
- evaluation budget;
- branch policy;

it can be estimated by

\[
\hat P_i=\frac{n_i}{N}.
\]

This must not be described as the literal geometric volume of a mechanism
family in continuous parameter space.

---

## 9. Operational definition of a primary family

Opposite primary assembly branches are always treated as different families.

For mechanisms on the same branch, define normalized coordinate differences.

For the five non-periodic coordinates,

\[
z_i
=
\frac{|x_i-y_i|}{x_{i,\max}-x_{i,\min}}.
\]

For phase,

\[
z_\phi
=
\frac{
|\operatorname{wrap}_{[-\pi,\pi)}(\phi_x-\phi_y)|
}{\pi}.
\]

The geometric distance is

\[
d(x,y)
=
\sqrt{
\frac{1}{6}
\sum_i z_i^2
}.
\]

If the branches differ,

\[
d=\infty.
\]

Families are connected components of the graph containing an edge whenever

\[
d < d_{\mathrm{threshold}}.
\]

Using connected components makes the result independent of insertion order,
although chaining remains a reason not to over-interpret the absolute number of
families.

The reference experiment inspected thresholds:

- `0.03`;
- `0.04`;
- `0.05`.

The central retained documentary result uses `0.04`.

---

## 10. Fresh-island saturation evidence

The completed reference saturation experiment used:

- 512 fresh islands;
- no historical seed;
- strict alternation of the two primary assembly branches;
- 260 differential-evolution generations per island;
- population size 16;
- approximately 96 initial individuals per island;
- the complete broad primary bounds described above.

At family-distance threshold `0.04`, the 512 captured winners were distributed
among 16 observed families with counts:

```text
244, 228, 19, 4, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1
```

Therefore:

- observed families: `16`;
- singletons: `8`;
- doubletons: `3`;
- Good-Turing unseen capture-mass estimate:

\[
\frac{8}{512}=1.5625\%;
\]

- top-1 capture share:

\[
\frac{244}{512}=47.65625\%;
\]

- top-4 capture share:

\[
\frac{495}{512}=96.6796875\%;
\]

- top-10 capture share:

\[
\frac{506}{512}=98.828125\%.
\]

The useful conclusion is deliberately narrow:

> Under the declared search policy and budget, most algorithmic capture
> probability was concentrated in a small number of dominant primary basins.

This does **not** prove that:

- every geometrically possible four-bar family has been discovered;
- the family counts are invariant to clustering threshold;
- the largest capture basin produces the best complete six-bar;
- capture frequency is geometric phase-space volume.

The minimal retained counts and their derived statistics are stored in
[`repro/mechanism_synthesis_search/saturation_reference.toml`](repro/mechanism_synthesis_search/saturation_reference.toml).

---

## 11. Seeded design search and fresh-island saturation have different purposes

Historical or retained mechanism seeds are useful during design exploitation
because they preserve expensive discoveries and accelerate convergence.

They are inappropriate when the scientific question is whether a fresh search
repeatedly rediscovers the same primary basins.

The two protocols must therefore remain explicit and separate.

### Design exploitation

May use:

- retained families;
- warm starts;
- local continuation;
- previous champions.

Question:

> How efficiently can a useful design be improved?

### Fresh-island saturation

Uses:

- fixed broad bounds;
- independent deterministic random seeds;
- no historical candidate injection;
- balanced discrete branches;
- fixed optimizer policy and budget.

Question:

> Under this declared search policy, how concentrated is fresh-search capture
> among the observed families?

This distinction is represented explicitly by `SynthesisRequest` and
`FreshIslandPolicy`.

---

## 12. Downstream transformability matters

A primary four-bar is an intermediate motion generator, not the final piston
mechanism.

The downstream dyad can:

- rotate the effective motion direction;
- reshape timing;
- amplify or attenuate parts of the primary trajectory;
- change the relation between primary cadence score and piston-position fit.

Consequently:

> A mediocre or rare primary under the simplified primary proxy can still
> become valuable after downstream optimization.

This is why several mechanically different primary families should survive long
enough to receive a fair downstream search.

The four retained reference seeds are documented in
[`PRIMARY_FOUR_BAR_FAMILIES.md`](PRIMARY_FOUR_BAR_FAMILIES.md).

---

## 13. Synthesize the more restrictive piston first

For a paired machine, the most useful first target is generally the piston with
the sharper or more restrictive kinematic demand.

The recommended rule is:

> Identify the harder motion, synthesize that side first, then use its
> mechanism as an initialization for the easier side.

The identity of the harder side is not universal.

It may change with:

- motor versus refrigeration operation;
- temperature levels;
- valve timing;
- pressure ratio;
- thermodynamic target family.

The rule is therefore “hard side first”, not “always synthesize one named
cylinder first”.

---

## 14. Mirroring is initialization, not symmetry

A successful mechanism may be mirrored geometrically toward the second
cylinder to provide a strong starting point.

Afterward the second mechanism must be free to adapt.

Permanent symmetry would unnecessarily remove mechanically and
thermodynamically useful degrees of freedom.

The generic synthesis interface therefore records:

```text
mirror_policy = "initialization_only"
```

rather than enforcing equality between the two cylinder mechanisms.

---

## 15. Mechanical constraints and proxy constraints are different

Direct mechanical protections and search preferences must not be confused.

Examples of direct mechanical requirements include:

- full loop closure;
- avoiding toggles below an explicitly retained transmission limit;
- finite piston-rod closure;
- required motion topology;
- any declared packaging or mechanical-clearance limit.

Examples of proxy preferences include:

- favored transmission margin above the hard minimum;
- PCA directionality;
- cadence resemblance;
- geometric compactness preferences.

A useful rule is:

> Relax proxy constraints when evidence shows that they are artificial, but do
> not casually relax direct mechanical protections.

Study-specific mechanical limits belong to the study configuration or complete
mechanism assessment, not to an immutable universal score.

---

## 16. Primary score, piston fit and thermodynamic quality are different rankings

The synthesis campaign exposed three different notions of quality.

### Primary quality

A cheap score used to discover four-bar basins with useful motion topology and
cadence.

### Complete-mechanism kinematic quality

How well the full realizable linkage produces an appropriate piston motion
after the downstream dyad and slider are included.

### Thermodynamic quality

How well the actual machine performs when driven by that realizable mechanism.

These rankings need not agree.

A primary that is not the best proxy scorer may be highly transformable.

A complete mechanism that most closely fits the abstract piston curve need not
give the best thermodynamic result.

This is the reason for retaining several families through the hierarchy.

---

Complete-mechanism fitting must score the actual finite-rod piston/slider
displacement, not rocker angle or an intermediate point projection. Failure of
one fitted geometry or one bounded search does not establish a limitation of
the complete mechanism family.

## 17. Why position error becomes diagnostic

During early synthesis, normalized piston-position error is useful because it
provides a cheap, interpretable target.

It should not remain the final objective.

Once a viable mechanism exists:

- the thermodynamic solver can directly judge the machine;
- the mechanism may beneficially smooth target artifacts;
- small deviations from the abstract piston law may improve real performance;
- exact curve imitation can become an unnecessary constraint.

The final optimization therefore uses the study thermodynamic objective while
position error remains diagnostic.

The target is a map toward useful geometry, not a trajectory that the final
machine must copy exactly.

---

## 18. Recommended reusable workflow

For a new DADA machine configuration:

### Step 1 — obtain a useful thermodynamic reference motion

Optimize the thermodynamic machine before imposing a specific realizable
linkage.

### Step 2 — identify trustworthy motion structure

Determine:

- turnarounds;
- monotonic branches;
- important cadence changes;
- derivative regions that are physically meaningful;
- derivative regions that are parameterization artifacts.

### Step 3 — identify the more restrictive piston

Begin mechanism synthesis on that side.

### Step 4 — discover primary families

Search the six primary coordinates over broad bounds and both assembly
branches.

Retain diversity, not only the best scalar score.

### Step 5 — search the downstream dyad

Give several distinct primary families a fair downstream search.

### Step 6 — release the full mechanism locally

Allow all 15 continuous dimensions to adapt around successful constructions.

### Step 7 — initialize the opposite piston

Mirror and phase-align a good mechanism as a starting point only.

### Step 8 — adapt the opposite mechanism independently

Do not enforce permanent symmetry.

### Step 9 — optimize the pair thermodynamically

Release both mechanisms simultaneously and use the actual thermodynamic
objective.

### Step 10 — re-tune the surrounding machine

With mechanisms fixed, adapt relevant thermodynamic and exchanger dimensions.

### Step 11 — evaluate robustness separately

Manufacturing sensitivity and mechanical-margin robustness are not represented
by the primary cadence score.

### Step 12 — use fresh-island saturation only when needed

Run an explicitly seed-free search when the question is search-space coverage
rather than design exploitation.

---

## 19. Generalization

The long-term objective is not to restart a completely unconstrained mechanism
search for every operating point.

Known families should serve as reusable initial structures.

For a nearby target:

1. begin from several retained families;
2. locally adapt them;
3. compare their complete-mechanism behavior;
4. use the thermodynamic solver for the final ranking;
5. add fresh broad exploration when existing families cease to cover the useful
   design space.

This applies across changes in:

- source temperatures;
- working fluid;
- machine scale;
- motor versus refrigeration operation;
- operating frequency;
- exchanger configuration.

The difficult piston and the winning mechanism family may change.

---

## 20. Simpler mechanism classes remain valid competitors

The success of a six-bar synthesis workflow does not prove that six bars are
always necessary.

For applications where simplicity dominates, the same thermodynamic problem
should also be tested with simpler physical mechanism families.

A final engineering comparison may include:

- thermodynamic performance;
- mechanical efficiency;
- part count;
- joints and bearings;
- tolerance sensitivity;
- packaging;
- fabrication difficulty;
- cost;
- maintenance.

Mechanism complexity must earn its place through system-level benefit.

---

## 21. Software boundary

The reusable software abstraction belongs in:

```text
src/dada_solver/research/synthesis.py
```

It records:

- ordered synthesis stages;
- released coordinate groups;
- explicit separation between design exploitation and fresh-island saturation;
- retained family identities;
- mechanical constraints;
- position as the primary synthesis reference;
- acceleration as diagnostic only;
- mirroring as initialization only;
- the thermodynamic objective after paired release.

The historical primary cadence score does **not** belong in `src/`.

It is retained only as documentary reproduction material under:

```text
docs/repro/mechanism_synthesis_search/
```

This keeps the solver architecture general while preserving the reasoning that
led to the current synthesis method.

---

## 22. What should not be inferred

This method does not establish that:

- one primary proxy is universally optimal;
- only a fixed number of useful primary families exists;
- optimizer capture frequency equals geometric mechanism-space volume;
- the closest kinematic fit gives the best thermodynamic machine;
- a mirrored pair should remain symmetric;
- six-bar mechanisms are always preferable to simpler linkages;
- the reference 260 K cadence is universal.

The durable result is the hierarchical search architecture and the distinction
between discovery proxies, complete-mechanism quality and final
thermodynamic assessment.

---

## 23. Reproducing the documentary search policy

The minimal reproducible material is:

```text
docs/repro/mechanism_synthesis_search/policy.toml
docs/repro/mechanism_synthesis_search/saturation_reference.toml
docs/repro/mechanism_synthesis_search/verify.py
```

Run:

```bash
python docs/repro/mechanism_synthesis_search/verify.py
```

or:

```bash
python docs/repro/mechanism_synthesis_search/verify.py --check
```

The reproduction intentionally:

- does not import historical search scripts;
- does not read campaign outputs;
- does not retain the 512 individual winners;
- does not retain a historical best scalar score;
- stores no JSON;
- reproduces the important saturation statistics directly from the minimal
  family-count reference;
- records the final primary-search proxy as documentation rather than as a
  generic solver capability.
