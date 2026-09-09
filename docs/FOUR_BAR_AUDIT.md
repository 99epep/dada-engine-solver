# Audit of the earlier four-bar optimizer

## Status and provenance

The repository inspected at
`/home/edada/Documents/BB/Theorie/4barres/dada-engine-4bar-optimizer` is useful
prior work, but it is not a mechanical or thermodynamic specification for this
solver. Its 12 tests pass under Python 3.11 with the available NumPy/SciPy
environment. The tests mainly protect its own scoring and refinement behavior;
they do not validate a realizable DADA piston mechanism.

## Reusable elements

- planar crank, coupler, rocker and ground-link representation;
- circle-intersection position closure;
- Grashof prefilter within the particular search bounds;
- local frames attached to the rocker and coupler;
- support-point reconstruction;
- coarse enumeration followed by local refinement;
- explicit storage of rejected candidates and score components;
- high-resolution post-processing distinct from coarse search resolution.

These elements must be reimplemented behind the solver's kinematic volume-law
interface and independently tested. Copying the earlier objective would also
copy assumptions that are no longer accepted.

## Assumptions that must not enter the coupled solver

### Prescribed target waveform

The earlier objective scores agreement with a constructed plateau/rapid/slow
waveform and fixed angular windows. Those are kinematic preferences, not a
thermodynamic optimum. Passive valve events must remain pressure-driven.

### Transfer-synchronization symmetry test

The second motion is constructed as `f(-theta mod 360)` and never optimized.
This was intended as a heuristic test that a complete-cylinder transfer starts
and ends synchronously in both pistons, not as evidence that both complete
motions must be thermodynamically symmetric. The importance of exact endpoint
synchronization has not been established. The thermodynamic study has already
shown useful cases with independent `Lambda_S` and `Lambda_L`; cylinder size,
exchanger geometry and reservoir temperatures can also shift the best timing.
The synchronization metric may be retained as an observable or optional
constraint, but it must not prescribe passive valve events or global symmetry.

### Family E scoring mismatch

Family E ranks rocker angular displacement. Only after optimization does it
construct a separate rocker-mounted point and its Cartesian motion. The
physical piston input can therefore differ from the quantity that received the
score. A coupled optimizer must score the actual slider displacement and
volume law.

### Family F slider approximation

Family F equates the support-point X coordinate with piston displacement while
allowing nonzero Y motion. A pin connected to a translating piston through a
finite connecting rod does not generally have that displacement. The missing
rod length, piston-axis position and assembly branch must be explicit. If a
different follower or slot is intended, it requires its own constraint model.

### Precompression and angular filters

The minimum 20-percent precompression, plateau centers, A3 windows and score
weights are design heuristics. They must not become valve conditions,
thermodynamic validity limits or universal manufacturing constraints.

### Incomplete general Grashof classification

The earlier prefilter checks the Grashof inequality but relies on its search
bounds to make the selected crank the shortest link. A reusable mechanism API
must explicitly classify the inversion and verify continuous crank rotation.

### Assembly and singularity handling

Selecting the upper circle intersection gives one assembly branch for the
existing geometry convention. A production implementation must make the
branch explicit, preserve it continuously, report toggle/singularity margins,
and reject circuit or branch defects. Exact equality and near-tangency require
scale-aware tolerances.

### Derivatives and physical units

Finite differences in degrees are adequate for drawing and heuristic scoring,
but the thermodynamic interface requires derivatives with respect to radians:
`dV/dtheta`, and future actuator checks require `d2V/dtheta2`. Angular speed can
then convert them to time derivatives. No acceleration score may be presented
as a force or stress calculation.

## Integration sequence

1. Implement an SI four-bar geometry with explicit assembly branch.
2. Add an explicit slider/follower model, including its connecting geometry.
3. Return piston displacement and first/second angular derivatives.
4. Map displacement monotonically to configured cylinder volume limits.
5. Verify closure, periodicity, derivative consistency and singularity margin.
6. Run the thermodynamic solver with a fixed candidate geometry.
7. Optimize a free admissible volume law as an upper bound.
8. Fit and then directly optimize four-bar parameters against thermodynamic
   objectives and constraints.

Detailed stress, bearing, fatigue, inertia and friction calculations remain
outside this scope.

## Implemented clean reference

`dada_solver.four_bar` now implements the first reference output family. Two
independent coupler-rocker loops share exactly one crank radius and crank pin.
Each rocker carries an explicit local output point connected by a finite-length
rod to a declared prismatic-slider axis. Circle closure, assembly branches,
four-bar toggles and slider toggles are checked explicitly.

The slider coordinate and its analytical derivative with respect to crank
angle are mapped to configured thermodynamic cylinder-volume limits. Physical
slider stroke remains available separately, so normalization does not erase
the piston area required to realize a swept volume. The module is not yet a
configuration-file option and no geometry has yet been optimized or claimed as
a DADA candidate.

## Selected low-tech mechanism family

The primary search family uses one shaft, one rigid crank plate and one shared
crank pin `B`. Both couplers are attached at that same pin. Downstream of `B`,
the small- and large-cylinder assemblies are independent and may have different
coupler lengths, rocker lengths, fixed rocker-pivot positions, support points,
connecting rods and piston-axis locations. They therefore share the same input
circle without being required to produce symmetric output motions.

The planar kinematic model may solve both loops independently around the common
`A-B` crank. Real construction will require an explicit axial layering or forked
joint at the shared pin; collision and bearing-stack details are packaging
questions and are not silently certified by a two-dimensional solution.

A common shaft with distinct crank throws or pins remains a future, less
constrained comparison family. It is not included in the first low-tech search.
