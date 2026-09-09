# Physics Decisions and Model Handoff

This document records decisions that must not be silently reinterpreted when
the solver is continued in Codex, ChatGPT Work, ChatGPT or another development
environment. It distinguishes the published physical specification from
explicit project decisions and numerical implementation choices.

The project is distributed under `GPL-3.0-or-later`. This replaces the earlier
CC0 dedication. Third-party dependencies and referenced scientific material
retain their own terms.

## Authority and scope

### Interchangeable motion and exchanger architecture

The user authorizes an architecture-only refactor and immediate independent
free periodic spline laws. Four-bar kinematics is one family, not a cycle
assumption; free motion must not inherit plateau, Lambda, symmetry or valve-
timing targets. The factory retains existing direction conventions, while a
free law may choose arbitrary cylinder volumes at the zero-angle origin.
Uniform initial filling uses those actual volumes. No conservative equation or
physical loss assumption is changed. Microtube geometry becomes one explicit
exchanger implementation, with derived quantities owned by that implementation.
See [PLUGGABLE_MODELS.md](PLUGGABLE_MODELS.md) for interfaces, feasibility,
compatibility, tests and the next campaign's parameter-adapter boundary.

### Original four-bar motion and constant-area parallelization

The user requests restoring the original shared-crank four-bar seed (small
output approximately -56 degrees, large +56 degrees), retaining 25/325 deg C.
Test 50% more tubes with lengths divided by 1.5, retaining contact area and tube
hold-up while allowing the geometry-derived headers to grow. Compare against
an original-motion-only control so the effects can be distinguished. This does
not impose mirror symmetry on future optimization. See
[PARALLEL_EXCHANGER_TRIAL.md](PARALLEL_EXCHANGER_TRIAL.md).

### Temporary 25/325 deg C source temperatures

The user raises the hot-air inlet to 598.15 K while retaining the cold-air
inlet at 298.15 K. Optimize dimensions first, then reduce hot-source temperature
progressively. The 100 W useful-power objective, 2–10 Hz range and 66 L large-
cylinder ceiling remain unchanged. Retain P–V diagrams as optional diagnostics;
angle-based plots remain the primary phase-reading view. See
[HIGHER_TEMPERATURE_TRIAL.md](HIGHER_TEMPERATURE_TRIAL.md) for the controlled
comparison with the previous 175 deg C trial.

### Excess external air and independent opposed Lambdas

The user asks to keep external-air capacity from limiting this research step.
Size fixed external flow from expected peak internal flow with an explicit
capacity-rate margin, verify achieved ratios and check sensitivity to doubling
flow. Retain thermal-film resistance and pause heat exchange. The independent
lower-S/higher-L Lambda trial and restored pressure/event plots are recorded in
[EXCESS_AIR_TRIAL.md](EXCESS_AIR_TRIAL.md). This is a thermal-boundary experiment,
not a validated external-flow hardware selection.

### External-loss exclusion and doubled-exchanger trial

The latest user decision excludes external-air aerodynamic losses and fan
consumption from the current optimization balance, while retaining finite air
flow, temperature change and thermal resistances. This supersedes earlier
requirements to deduct fans for these trials; excluded losses are not physically
zero. The requested doubled tube count and length and the independent lower-
Lambda mechanical trial are recorded in [DOUBLED_EXCHANGER_TRIAL.md](DOUBLED_EXCHANGER_TRIAL.md).
Checkpoint timing, optional wall initial-guess acceleration and full diagnostic
plots accompany the new results. Motion-law optimality remains unproven.

### Geometry-connected exchanger model

The user authorized connecting tube/header geometry, gas/metal/air thermal
resistance, wall storage, hydraulic loss and fans. The opt-in hardware connector
now replaces exchanger volumes and all four adjacent passage closures, splitting
tube and header resistance equally around each storage node and retaining
outlet valve losses. External airflow and fan power have explicit assumptions.
Constant-property/Nusselt and laminar closures remain screening models, not
Doty calibration. See [GEOMETRIC_MOTOR_COUPLING.md](GEOMETRIC_MOTOR_COUPLING.md).

### Opt-in air/wall motor dynamics

The authorized pause-heat extension now adds two wall energies in a separate
motor wrapper using continuous ideal diodes. Existing eight-state operation
remains unchanged. External air has finite heat-capacity rate; gas heat remains
active during pauses. Gas plus wall energy is conserved against external heat
and boundary work. Conductance split and wall capacity remain explicit uncalibrated
inputs. See [AIR_WALL_COUPLING.md](AIR_WALL_COUPLING.md) for equations, verification
and the counterflow calibration limitation.

### Air-source demonstrator exchanger decision

The user selected external air for heating and cooling the demonstrator and
authorized accounting for heat exchange during flow pauses. Use finite external
heat-capacity rates and count fan power separately; source temperatures are
external-air inlet conditions. Optimize tube geometry and collectors jointly,
including working-gas hold-up. More shorter tubes at constant diameter and
surface do not inherently increase tube volume; collectors can dominate the
increase. See [AIR_SOURCE_EXCHANGERS.md](AIR_SOURCE_EXCHANGERS.md) for implemented
geometry, assumptions and the planned wall-storage extension. No dynamic cycle
equations have been changed by the geometry implementation.

### Motion-law terminology and long-term objectives

The user clarified that the piecewise-linear reference has not been shown to
be optimal. Desired chronology does not imply that this motion law is ideal.
Mirror symmetry is a historical refrigerator choice, not a motor constraint;
allow independent loop geometry subject to shared-crank compatibility.
The ordered long-term objectives are (1) search for optimal motion laws and
(2) map independent Lambdas and cylinder swept-volume ratio against reservoir
temperatures and useful-power requirements, with other design assumptions
explicit. See [MOTOR_RESEARCH_OBJECTIVES.md](MOTOR_RESEARCH_OBJECTIVES.md).
These decisions supersede any contrary interpretation of historical text.

### Current motor design brief, 2026-09-09

All project content must be English. The motor is a small experimental
demonstrator, not a commercial product; commercial refers to purchased
exchangers. Its name does not prescribe a permanent temperature difference. The motor targets approximately 100 W of
useful mechanical output at 2–10 Hz, between reservoirs at 298.15 and 448.15 K.
The large cylinder must not exceed 0.066 m^3 maximum enclosed gas volume,
interpreted conservatively as including clearance. This is an upper limit,
not a prescribed swept volume. Overall machine envelope remains unspecified;
exchanger, manifold and mechanism size must also be reported.

Maximize efficiency subject to useful power and size requirements. Current
`motor_power` is indicated gas work per unit time; friction, transmission and
auxiliary losses are absent, so 100 W indicated does not establish 100 W useful.
No assumed loss allowance is a confirmed user requirement.

Prefer a documented commercial exchanger and adapt the motor to it. If adequate
commercial data are unavailable, extrapolation from Doty's measurements is
explicitly authorized, with source conditions, scaling assumptions and
uncertainty scenarios reported. This supersedes earlier prohibitions on
extrapolation. See [MOTOR_DEMONSTRATOR.md](MOTOR_DEMONSTRATOR.md).

### User decisions following the motor update, 2026-09-08

Motor efficiency is the primary objective. Ideal chronology remains the design
target. Isothermality is deferred and its excursion is now diagnostic-only in
the core validity verdict; no isothermality constraint is added to motor sizing.
Historical explicitly configured sizing constraints remain explicit opt-ins.

The intended exchanger through-flow is discontinuous and unidirectional. The
existing bidirectional internal links are retained; any modeled local reflux is
exported as a discrepancy, not erased. Literature comparison must use pulse
shape, active duration, flow and cycle frequency rather than assume sinusoidal
zero-mean oscillation. Plate-fin and parallel microtubes are both candidates.
The motor case has no low-tech fabrication requirement inherited from the
human-powered cell. Literature/manufacturer data precede additional transient
physics; explain the purpose and cost before any such extension. See
[EXCHANGER_LITERATURE.md](EXCHANGER_LITERATURE.md).

As of 2026-09-08, the scientific source of truth is English revision 1288 of
the DADA Thermodynamic and Mechanical Study:

https://dada-engine.org/index.php?title=Thermodynamic_and_Mechanical_Study&oldid=1288

The motor update is recorded in [MOTOR_OPERATION.md](MOTOR_OPERATION.md).
It supersedes the refrigeration-only reservoir naming and positive-only
configuration speed in the historical decisions below. The original reference
was revision 1277 of the French study:

https://dada-engine.org/index.php?title=Thermodynamic_and_Mechanical_Study/fr&oldid=1277

The initial implementation date for the decisions below is 2026-09-02. The first-level
scope is a single-phase, calorically perfect ideal gas in four well-mixed 0D
control volumes. Any change to the physical equations must be discussed and
documented before implementation.

The DADA flow is described as internal, periodic, strongly pulsating and
topology-varying. It is intermittent in check-valve branches and potentially
reversible in bidirectional branches. It must not be reduced globally to either
steady flow or the continuously oscillating flow commonly associated with a
Stirling regenerator.

## Confirmed physical decisions

### Four independent conservative control volumes

The solver uses option B from the design discussion: `S`, `L`, `C` and `H`
always remain independent conservative control volumes with their own masses,
internal energies, temperatures and pressures.

The fundamental state is always:

```text
X = (m_S, U_S, m_L, U_L, m_C, U_C, m_H, U_H)
```

Temperature and pressure are reconstructed from this state. They are not
integrated as independent thermodynamic coordinates.

The reduced quasi-pressure-equalized pair equations in the study are retained
as analytical limits and validation references. They are not imposed as the
simulation state equations. Consequently, pressure differences across the
internal `S-C` and `L-H` links are resolved rather than forced to zero.

### Hydraulic graph

The first-level graph contains four configured links:

- `L <-> H`: bidirectional internal link;
- `S <-> C`: bidirectional internal link;
- `H -> S`: passive one-way check valve;
- `C -> L`: passive one-way check valve.

Every link requires its own configured effective flow area `CdA`. No hydraulic
area is inferred or supplied as a hidden default.

The first-level flow closure is the compressible-orifice equation from the
study, including choked and unchoked branches. The actual pressure gradient
selects the upstream state of a bidirectional link. Transported enthalpy is
always evaluated at the actual upstream temperature. A check valve never
permits reverse flow.

The hydraulic closure is replaceable. A future exchanger, pipe or detailed
valve model must implement the hydraulic interface without changing the mass
and energy balances.

### Passive valve behavior and topology

Valve state is determined only by pressure difference and configured physical
hysteresis:

```text
open when  delta_P >= delta_P_open
close when delta_P <= delta_P_close
```

with `delta_P_close <= delta_P_open`.

Angles and nominal `Lambda*` values never enable, disable or schedule a valve.
They are diagnostic or design quantities only.

All physically generated topologies are simulated. Simultaneous openings,
repeated openings, missing phases and unexpected event orders are not blocked.
They are classified after integration as `nominal` or `non_nominal`.

At a valve event, the complete conservative state and imposed geometry remain
continuous. No instantaneous pressure equalization or state reset occurs.

### Event-root numerical treatment

When the event solver locates an oriented threshold root, the corresponding
discrete valve transition is applied exactly. Re-evaluating the same inequality
can differ from the located root by a few floating-point ulps, so it is not used
to decide whether the already-located transition occurred.

The other valve is still evaluated passively at the same angle and state. This
preserves simultaneous physical events. This is a numerical representation of
the threshold crossing, not an added pressure tolerance or dead band.

The result retains two samples at every nonzero event: one with the preceding
topology and one with the following topology. Their values of `X`, cumulative
heat and cumulative work must be exactly identical.

### Heat transfer

The only first-level reservoir closure is:

```text
Q_dot = UA * (T_reservoir - T_HX)
```

for both `C` and `H`. Heat is positive when received by the gas. A nominally
isothermal phase does not impose a temperature and does not disable finite
temperature differences.

The heat exchangers remain thermally coupled during every topology unless a
configured `UA` is zero.

### Work and energy signs

Gas work is positive when supplied by the gas:

```text
W_dot = P_S * V_dot_S + P_L * V_dot_L
```

Because the four-volume model resolves each cylinder pressure, the pressure
used in its boundary-work term is that cylinder's own pressure.

The intended refrigeration signs are:

```text
Q_C > 0
Q_H < 0
W_cycle < 0
```

Cooling and heating COP values are reported only when their required heat and
mechanical-input signs are meaningful. An algebraic negative refrigeration COP
is reported as `unavailable`, while the signed heat and work values remain
available for diagnosis.

### Charge and initial condition

At `theta = 0`, the large cylinder is at maximum volume. The initial numerical
state is a uniform filling equilibrium at the configured charge temperature
and either:

- a configured charge pressure; or
- a configured total gas mass, from which the uniform filling pressure is
  derived.

Mass in each volume is assigned from its volume at `theta = 0`. This state is a
filling condition and is never described as the periodic operating state.

No renormalization of mass is performed between cycles. Total-mass drift is
measured and reported.

### Periodic steady state

The initial periodic solver uses successive complete cycles. Its continuous
fixed-point condition is `X(theta + 2*pi) = X(theta)` with componentwise mass
and energy scaling.

The discrete valve topology must also be recovered. Equal continuous states
with different valve states are not the same dynamical state when hysteresis is
present.

The solver distinguishes convergence, maximum cycle count, invalid physical
state and numerical integration failure. A future shooting or Newton method
must preserve the same physical state and diagnostic definitions.

### Validity verdicts

Validity thresholds are configuration inputs, not embedded universal
constants. The solver reports pressure-equalization and isothermality errors,
ideal-gas compressibility deviation and property-variation hooks.

At first level, `Z = 1` and constant `Cp` produce zero modeled deviations. This
does not prove that a chosen real working fluid is adequately ideal over the
simulated pressure-temperature domain.

`CdA` is not a geometric flow area. Without a separately supplied real section,
gas velocity and Mach number are `unavailable`. If an essential criterion is
unavailable and no other criterion fails, the overall verdict is
`indeterminate`, not `valid`. Any violated configured criterion produces an
`invalid` verdict.

### Units

All stored and calculated physical quantities use SI units. In particular,
configuration angular speed is in `rad/s`, not rpm. Plot angles may be displayed
in degrees as a diagnostic representation.

### Kinematics boundary

The thermodynamic solver depends only on volume laws and their angular
derivatives:

```text
small_cylinder_volume(theta)
small_cylinder_volume_derivative(theta)
large_cylinder_volume(theta)
large_cylinder_volume_derivative(theta)
```

The built-in harmonic law is a controlled numerical example only. Its TOML
configuration must explicitly contain `example_data = true`. It is not a DADA
mechanism model or a validated design.

## Numerical verification decisions

Near-equal floating-point pressures are not made equal using an arbitrary
hydraulic tolerance. If a validation case requires an exactly isolated volume,
the test uses an explicitly zero-flow test link. This keeps the analytical
hypotheses distinct from the finite-resistance physical model.

The appendix relations are implemented as independent verification residuals,
including the closed adiabatic `P * V**gamma` invariant and the outflow-only
adiabatic donor relation. Global mass and energy residuals are calculated from
the conservative state and integrated boundary exchanges.

Numerical convergence establishes consistency of the discretization only. It
does not override a `non_nominal`, `invalid` or `indeterminate` physical
diagnosis.

## Explicitly deferred physics

The following items have not been invented or implemented as detailed physics:

- real-gas equations of state and variable thermophysical properties;
- condensation or phase-transition limits for a selected fluid;
- geometric passage areas and Mach-number evaluation;
- distributed exchanger, pipe and valve pressure-loss models;
- leakage and mechanical friction;
- piston, linkage, shaft, inertia, bearing and structural calculations;
- conversion from swept volume to detailed cylinder geometry;
- experimental calibration or validation;
- a final sizing objective or optimization policy.

Thermodynamic gas force may later be exposed as `P * piston_area`, but net
mechanical force and structural sizing require additional mechanical data.

## Software status at version 0.1.0

The conservative simulation core, event integration, successive-cycle search,
performance reporting, plotting and controlled analytical tests are present.
The complete project requested in the roadmap is not finished. Initial sizing
and mechanical-boundary interfaces are present, but validated design studies,
Pareto exploration and detailed mechanical response models remain future work.

## Initial sizing and mechanical-interface decisions

The first sizing layer does not select a final physical objective. It composes
bounded SI design variables, one interchangeable scalar objective and explicit
physical constraints. Periodic non-convergence or an unavailable constrained
quantity makes a design point infeasible.

The SLSQP adapter normalizes bounded variables and requires all numerical
objective and constraint scales from the caller. These scales are optimizer
controls, not physical constants. A different optimizer or Pareto explorer may
reuse the same sizing problem.

Multiobjective comparison uses only caller-supplied design points and
minimization-valued objectives. The Pareto filter removes infeasible,
unavailable and dominated points but never chooses a preferred tradeoff or
creates objective weights.

Feasibility exploration uses a seeded Latin hypercube inside explicitly
configured bounds. It is a numerical coverage strategy, not a physical prior.
An empty finite sample does not prove that the physical feasible set is empty.

The current mechanical boundary exposes gas-side `P * area` forces only when
piston face areas are explicitly supplied. It also exposes the work-conjugate
generalized gas torque `sum(P_i * dV_i/dtheta)`. Neither quantity includes
opposing pressures, linkage statics, inertia, friction, bearings or structural
stress. A future mechanical response model must add those effects from new
mechanical inputs rather than reinterpret the thermodynamic loads.

Sizing may optionally constrain the maximum absolute generalized gas torque.
This remains the thermodynamic work-conjugate quantity and is not renamed or
interpreted as a complete shaft torque.

A strict `valid_thermodynamic_model` sizing constraint requires a `valid`
verdict. Because Mach is unavailable without geometric flow areas, such a
constraint will intentionally reject current first-level points as
`indeterminate`. Users must not remove the constraint merely to relabel an
unknown validity criterion as satisfied.

## External cooling-cell load

The first cooling-cell task remains external to the DADA thermodynamic state.
Its ideal energy is the sum of sensible cooling in the initial phase and an
isothermal phase-change energy. Material mass, temperatures, specific heat,
latent energy and transformed fraction are all configuration inputs.

The one-kilogram water example uses explicit approximate values of
`4180 J/(kg K)` and `333500 J/kg`, giving `417100 J`. It excludes container heat
capacity, parasitic ambient gains, contact resistance, gradients and
supercooling. Results must retain the `ideal` qualifier until those effects are
added.

For a fixed machine point, ideal task time is load energy divided by positive
machine cooling power. If the point requires more than the selected human input
limit, it is not rescaled linearly; a new operating point must be simulated.
Pedaling power inferred from task cooling power divided by the evaluated COP is
a local diagnostic and does not assert constant COP under rescaling.

The cooling-cell sensitivity study initially used `263.15 K`; the user then
raised the machine cold-reservoir temperature by 5 K to `268.15 K`. The external water load still
ends near `273.15 K`; these temperatures must not be conflated. The colder
reservoir can increase thermal driving difference but may reduce COP.

Cylinder clearance is parameterized as `V_min / V_swept` with exploratory
bounds `0.005` to `0.10` for each cylinder. `V_C` and `V_H` remain separate
heat-exchanger control volumes and are not included in cylinder clearance.

The harmonic small-cylinder phase may be swept over a full turn solely to
diagnose the dependence of passive valve topology and refrigeration direction
on imposed kinematics. This does not promote the full-turn interval to a
validated mechanical design range, and valve transitions remain pressure-driven.

At the initial cooling-cell exploratory point, a five-point cardinal phase
sweep produced negative cold-reservoir heat-transfer power at every point
(approximately -38 W to -107 W), despite periodic convergence. All observed
topologies were non-nominal and all validity verdicts were invalid. This result
must be retained as a failed operating-point diagnosis: it must not be converted
into positive cooling capacity by sign changes, absolute values or geometric
scaling.

Before modifying any item in this document, record the proposed physical
change, its source or justification, affected equations, compatibility impact
and required regression tests.

## Ideal piecewise-linear kinematics and four-bar context

The home page and the French four-bar-mechanism page were reviewed in full on
2026-09-02. The ideal-cycle GIF is qualitative explanatory material and carries
no quantitative scientific accuracy. Pixel coordinates, frame counts and visual
intermediate positions must not be used as physical input data.

The four-bar page describes a realizable approximation with a near-dwell over
one quarter-turn, a full stroke over the following quarter-turn, and a return
over the remaining half-turn. Its published example geometries are mechanical
candidates, not the thermodynamic solver's source of imposed volume laws.

A piecewise-linear imposed-volume law is implemented as an ideal diagnostic
kinematics. Its slope discontinuities imply unbounded ideal acceleration and
therefore make it unsuitable for inertia or structural-load claims. Passive
valve events remain pressure-driven and are never scheduled at kinematic
breakpoints.

The study defines `Lambda = (V_max - V) / (V_max - V_min)`: zero is maximum
volume and one is minimum volume. The first search initially held both targets
equal at 0.7 and later swept them together, while the general kinematics retains
separate small- and large-cylinder parameters. In the implemented receiver
sequence, `L` is compressed from Lambda 0 to `Lambda_L_target`, and `S` is
expanded from `Lambda_S_target` to Lambda 0. Increasing a common target therefore
increases both adiabatic stroke fractions. It simultaneously reduces the
remaining hot-transfer stroke fraction `1 - Lambda_target`. No general monotonic
temperature-lift or performance claim may be inferred without simulation.

### Corrected angular origin and withdrawn earlier results

The study fixes theta zero at `V_L = V_L_max`, at the start of Phase I. The
first implementation of the piecewise-linear law started at Phase III with
`V_L = V_L_min`. Although this is a cyclic phase rotation, it materially changed
the configured uniform-charge inventory because charging pressure and
temperature are applied at theta zero. All quantitative piecewise-linear results
reported before this correction are withdrawn. They must not be used for design
or physical conclusions. The law now starts with Phase I and follows the study's
I-II-III-IV order; regression tests protect the origin volumes and charge mass.

For the corrected doubled-cylinder/reduced-HX point at common Lambda 0.70, the
theta-zero filling volume is 852 cm3 and the 2 bar, 298.15 K inventory is about
1.991e-3 kg. The corrected periodic result is still non-refrigerating and
invalid, but it is primarily a topology diagnostic rather than a machine
performance result. The `C -> L` valve closes at 31.85 degrees and `H -> S`
opens at 32.35 degrees during the nominal 0--54 degree Phase I, leaving only
about 0.5 degree with both valves closed. During nominal Phase III (180--234
degrees), `H -> S` remains open throughout; `C -> L` opens at 234.02 degrees and
`H -> S` closes at 234.28 degrees. Thus the configured point does not realize
the nominal closed-valve adiabatic phases. Its negative cold power and extreme
temperatures cannot be used to reject the ideal DADA cycle.

The proposed benefit of increasing the disparity between cylinder swept
volumes is retained as a design hypothesis to test, not as an implemented
physical law. Isothermality also depends on finite UA, transferred gas thermal
capacity, hydraulic losses and available time; a larger cylinder-volume ratio
does not by itself establish lower isothermality error or higher efficiency.

The initial duration study uses two equal adiabatic-sector fractions and two
equal transfer-sector fractions, with `f_transfer = 0.5 - f_adiabatic`. Four
equal sectors are the reference point. A diagnostic sweep from
`f_adiabatic = 0.25` down to `0.05` improved signed cold heat-transfer power
from approximately -32.1 W to -19.6 W and reduced mechanical input from 75.7 W
to 51.6 W. It did not produce refrigeration; every sampled topology remained
non-nominal and every validity verdict invalid. This is local evidence for a
duration sensitivity, not proof that longer transfer sectors are globally
optimal.

At `268.15 K`, a common-Lambda sweep from 0.70 to 0.95 with equal sectors had a
local best signed cold power near Lambda 0.90 (-23.8 W), but never produced
refrigeration. The base exploratory point therefore uses equal Lambda targets
of 0.90 as an observed local candidate, not a global optimum. At Lambda 0.90, a
duration sweep had a local best near `f_adiabatic = 0.15` (-20.1 W); extending
the transfer sectors further degraded the result. These interactions prohibit
assuming monotonic benefit from either Lambda or transfer duration.

## Doubled-cylinder, reduced-HX geometry trial

The user requested an exploratory geometry with both swept volumes doubled,
both heat-exchanger volumes divided by five, and cylinder clearance reduced to
one percent of swept volume. The resulting values are: small cylinder 4 to
404 cm3, large cylinder 8 to 808 cm3, and 20 cm3 for each fixed heat exchanger.
The operating point retained 268.15 K cold reservoir temperature, equal Lambda
targets 0.90, and adiabatic-sector fraction 0.15.

The periodic solution converged, but signed cold power degraded from the prior
local candidate's approximately -20.1 W to -26.4 W, while mechanical input rose
from approximately 63.5 W to 117.3 W. Pressure-equalization and isothermality
errors also increased; topology remained non-nominal and validity remained
invalid. This is a coupled geometry result at unchanged UA and CdA. It does not
establish that larger cylinders or smaller HX volumes are independently
detrimental, because their effects were not isolated and hydraulic/thermal
capacity was not rescaled.

An initially reported Lambda sweep from 0.40 to 0.70 did not use the intended
doubled-cylinder/reduced-HX geometry: the sizing-space initial values overrode
the updated base geometry with the preceding 200/400 cm3 swept volumes, 100 cm3
HX volumes and five-percent clearances. Those numerical conclusions are
withdrawn and must not be used. The configuration mismatch was found by
comparing the direct-simulation inventory and extrema with the sensitivity
output. The sizing-space initial geometry was then synchronized to 400/800 cm3
swept volumes, 20 cm3 HX volumes and one-percent clearances before repeating the
study.

## Slower-speed and larger-CdA diagnostic

The corrected doubled-cylinder, 20 cm3 HX and one-percent-clearance geometry
was retained while angular speed was reduced from 10 to 5 rad/s and all four
effective CdA values were doubled. This improves hydraulic capacity relative to
imposed volume rate by a factor of four. The conservative periodic solution
converged with relative mass and energy residuals below 7e-16 and 3e-15.

Maximum temperature fell from about 1275 K to 498 K, maximum pressure from about
90 bar to 13.3 bar, and mechanical input from about 1677 W to 162 W. This
strongly supports the user's hypothesis that the extreme endpoint peaks were
primarily caused by inadequate evacuation relative to piston motion, rather
than one-percent clearance or the piecewise-linear volume law alone.

The point still did not refrigerate: signed cold power was about -27.0 W,
topology remained non-nominal, and isothermality remained invalid. `C -> L`
closed at 9.41 degrees and `H -> S` opened at 10.57 degrees, leaving only about
1.16 degrees of closed-valve Phase I; `H -> S` remained open throughout Phase
III. Hydraulic peak mitigation and passive-valve sequencing are therefore
separate problems.

More aggressive simultaneous changes (10x CdA with 0.1x speed, then 5x CdA
with 0.2x speed) reached the positivity boundary of the conservative state with
both RK45 and Radau. These attempts are numerical failures, not physical
results. A configurable stiff integration method was added, but further work
requires state scaling or a positivity-preserving formulation before these
more widely separated time scales can be evaluated reliably.

## Zero-hysteresis valve interpretation

A requested trial with `delta_P_open = delta_P_close = 0` exposed an ambiguity
in the discrete event formulation at exact pressure equality. With inclusive
published inequalities, a closed valve qualifies to open and an open valve
qualifies to close at the same state. The uniform filling state starts exactly
at this equality, and `solve_ivp` detects a terminal root at the initial angle
even after adopting the tie convention that the valve is closed at equality.

No physical result is reported for this failed trial, and no hidden numerical
pressure epsilon is introduced. A true zero-hysteresis model should be added as
a continuous ideal hydraulic diode: forward orifice flow for strictly positive
pressure difference and zero flow otherwise, without a discrete valve state at
the equality boundary. Until that distinct formulation and its topology
post-processing are implemented and tested, the last convergent explicit-
hysteresis configuration remains the example baseline.

## Cold-to-large hydraulic sensitivity

The sizing-space initial point was synchronized with the current exploratory
configuration before this study: angular speed 5 rad/s, internal-link CdA
values 1e-5 m2, check-valve CdA values 4e-6 m2, equal Lambda targets 0.70,
400/800 cm3 swept volumes, 20 cm3 heat exchangers, and one-percent cylinder
clearances. The sensitivity tool gained optional local bounds so that a
scientifically selected interval can be swept without changing the broad
exploratory design bounds.

Only the `C -> L` effective CdA was varied. At 8e-6, 1.6e-5, and 3.2e-5 m2,
the periodic solver converged and respectively predicted 8.44, 15.85, and
17.52 W signed cold power; 71.89, 51.56, and 46.81 W mechanical input; and
cooling COP values 0.117, 0.307, and 0.374. These are sensitivity results, not
validated designs. All three cycles were non-nominal and invalid under the
configured validity limits.

The residual pressure difference `P_C - P_L` at theta zero fell from 33.23 kPa
to 6.63 kPa and then 1.58 kPa. Maximum temperature fell from 436.9 K to 422.7 K
and then 419.3 K. This supports the hypothesis that the `C -> L` restriction
was a major cause of the Phase-IV pressure split and thermal peak. Its effect
on valve timing was more limited: `C -> L` closing advanced from 2.79 to 1.74
degrees, while `H -> S` opening advanced in parallel from 4.05 to 3.01 degrees.
The interval with both check valves closed remained approximately 1.26 to 1.27
degrees. Enlarging this passage alone therefore does not restore the nominal
Phase-I duration.

The reported pressure-equalization error increased from 0.102 to 0.161 over
the same sweep. This indicator is correctly associated with the permanently
connected internal pairs `S-C` and `L-H`, as defined in Section 4 of the study;
it is not a measure of the external `C -> L` valve pressure difference. The
result indicates that relieving `C -> L` transfers more hydraulic demand to at
least one unchanged internal link. The internal and external passage areas
must consequently be studied together, while remaining independent design
variables.

Relative cycle mass and energy residual magnitudes stayed below 1.8e-15 and
4.9e-15 for the converged points. The result is therefore conservative to
numerical precision even though the configured thermodynamic-validity and
topology criteria fail.

The 4e-6 m2 point reached a non-positive conservative state during the fifth
successive cycle and is reported as `invalid_physical_state`, with no
performance assigned. This conflicts with an earlier manually reported
converged result at nominally the same hydraulic point. The current committed
configuration reproduces the failure, not that earlier result. The earlier
quantitative claim (about -27 W cold power, 162 W mechanical input, and valve
events near 9.41 and 10.57 degrees) is withdrawn until its exact configuration
or numerical cause can be recovered. It must not be used as a baseline.

## 33-liter machine and small-cylinder ratio study

An exploratory large-machine configuration was introduced with 33 L large
swept volume, one revolution per second, one-bar absolute filling pressure,
one-percent cylinder clearances, equal Lambda targets of 0.70, and 0.66 L gas
volume in each heat exchanger. Both UA values are provisionally 25 W/K from
the user's 500 W at 20 K exchanger estimate. This estimate is not a validated
heat-exchanger characteristic.

The prudent hydraulic envelope uses half the proposed geometric free areas:
4.2e-4 m2 for `L-H` and `H -> S`, and 7.5e-4 m2 for `S-C`. The `C -> L` check
valve was provisionally made 20 percent less restrictive than `S-C`, giving
9.0e-4 m2. This 20-percent choice is an explicit exploratory assumption, not a
published physical requirement. Effective CdA remains distinct from actual
free area, so Mach number is still unavailable to the solver.

Small-to-large swept-volume ratios 0.10, 0.20, 0.35, 0.50, 0.70, and 1.00 were
tested at constant filling pressure. The corresponding small swept volumes are
3.3, 6.6, 11.55, 16.5, 23.1, and 33 L. Holding pressure rather than total mass
means that each geometry has its own physically consistent gas inventory.

Ratios 0.10 and 0.20 reached a non-positive conservative state after two and
eight cycles. They have no reported performance and cannot be interpreted as
physical failures without further positivity/stiffness investigation. Ratios
0.35, 0.50, 0.70, and 1.00 converged with relative mass and energy residual
magnitudes below 6.3e-15 and 2.4e-15. Their signed cold powers were -213.6,
-163.1, -145.8, and -181.8 W, respectively; all heat the cold reservoir rather
than refrigerating it. Mechanical inputs were 613, 677, 803, and 1054 W.

Increasing the ratio from 0.35 to 0.70 reduced maximum temperature from 426.5
to 410.2 K and maximum pressure from 9.23 to 5.26 bar, but increased mechanical
input. The thermal trend reversed between 0.70 and 1.00, indicating an
intermediate local optimum rather than monotonic improvement. The equal-volume
point produced nominal valve-event ordering, while all smaller converged
ratios were non-nominal. Nominal topology alone is therefore insufficient for
refrigeration.

The 100 W cold target would require +100 J cold heat per cycle at one
revolution per second. None of these first ratio points approaches that target;
changing the cylinder ratio alone does not repair the cycle. The next study
must inspect the complete pressure, temperature, flow, heat-rate, and valve
chronology near the least-negative region before varying multiple design
parameters.

## Constrained high-COP exploration

A reproducible staged exploration varied `S/L`, `Lambda_S`, and `Lambda_L` at
the 33-liter, 9.091 rpm operating point. Feasibility required at least 100 W
cold power, no more than 150 W modelled mechanical input, periodic convergence,
pressure-equalization error at most 0.05, and both isothermality errors at most
0.05. Non-nominal topology was recorded but not prohibited. Mach number
remained unavailable and was not falsely treated as satisfied.

Latin-hypercube coarse and refinement results were saved incrementally as
`outputs/cop_coarse_exploration.csv`, `outputs/cop_refinement_exploration.csv`,
`outputs/cop_final_refinement.csv`, and `outputs/cop_boundary_extension.csv`.
A local cross-check and cooling-constraint scan were saved as
`outputs/cop_robustness_cross.csv` and
`outputs/cop_cooling_constraint_scan.csv`. Pathologically slow regions were
interrupted rather than assigned a physical result.

The best standard-tolerance point observed was `S/L = 0.84`, `Lambda_S = 0.43`,
and `Lambda_L = 0.44`, producing 111.8 W cold power, 30.4 W mechanical input,
and COP 3.674 with nominal topology. However, it failed immediately with a
non-finite state when integration tolerances were tightened tenfold and the
maximum angular step was halved. It is retained only as an unverified boundary
candidate and must not be reported as the reference result.

The more conservative point `S/L = 0.84`, `Lambda_S = 0.49`, and
`Lambda_L = 0.44` was stable under the tightened verification settings. It
produced 133.0696 W cold power, 39.3206 W mechanical input, and cooling COP
3.384225. Its minimum and maximum gas temperatures were 237.83 and 342.58 K,
maximum pressure was 2.263 bar, pressure-equalization error was 0.00390, and
cold/hot isothermality errors were 0.01798/0.02284. Event ordering was nominal;
relative mass and energy residual magnitudes were below 3.0e-15. The validity
verdict remains indeterminate solely within the implemented checks because
Mach number is unavailable. Real-property and moisture checks are also still
required before physical validation.

## Charging-pressure sensitivity of the high-COP candidate

At fixed 9.091 rpm geometry, Lambda targets, UA values, and CdA values, charge
pressure was tested at 1, 2, and 3 bar absolute. All periodic searches
converged in five or six cycles with nominal topology and relative mass/energy
residual magnitudes below 3.0e-15.

At 1, 2, and 3 bar, cold powers were respectively 133.07, 218.13, and 275.69 W;
mechanical inputs were 39.32, 67.33, and 88.84 W; and cooling COP values were
3.384, 3.240, and 3.103. Power rises sublinearly with pressure because fixed UA
must condition proportionally more gas inventory. The COP decreases gradually,
not catastrophically.

Maximum pressures were 2.263, 4.577, and 6.921 bar absolute. Global gas
temperature ranges were 237.83-342.58, 239.00-341.01, and 239.84-339.89 K.
Despite the narrower global extrema, heat-exchanger isothermality errors rose
with pressure: cold/hot values were 0.01798/0.02284, 0.03191/0.03996, and
0.04339/0.05333. The 3-bar point consequently fails the configured hot
isothermality threshold of 0.05; 1 and 2 bar satisfy all currently calculable
thresholds but remain indeterminate because Mach number is unavailable.

The 2-bar point is the strongest currently acceptable fixed-UA pressure trial.
Increasing pressure while scaling UA proportionally is a separate similarity
study and must not be conflated with this fixed-exchanger sensitivity.

## Reduced-volume similarity checks

All cylinder and heat-exchanger gas volumes of the verified high-COP candidate
were divided by ten while angular speed was multiplied by ten. Charge pressure,
UA values, and CdA values were unchanged. The model reproduced the same
pressure, temperature, topology, power, and COP histories: 3.3 L large swept
volume, 2.772 L small swept volume, 0.066 L per heat exchanger, 90.91 rpm,
133.07 W cold power, 39.32 W mechanical input, and COP 3.384225. Energy per
cycle and gas inventory were divided by ten. This confirms the expected
first-level similarity, not the physical manufacturability of such a compact
exchanger.

A second exact-similarity point targeted 100 W cold power with the same
one-tenth volumes. Speed was 68.318 rpm, both UA values 93.936 W/K, and every
CdA was 0.751486 times the high-COP reference value. It converged with nominal
topology, 100.000 W cold power, 29.549 W mechanical input, and unchanged COP
3.384225 and validity indicators. The required maximum effective area was
1.127e-3 m2 for `C -> L`.

These results show that cylinder size is not the first-level obstacle: the
same gas throughput per unit time can be obtained from smaller swept volumes
at higher speed. The demanding requirement is retaining roughly 94 W/K of UA
and millimeter-scale effective flow areas while limiting each exchanger's gas
volume to 66 cm3. The external hardware volume is not represented. Heat-
exchanger geometry and measured pressure-drop/effectiveness data therefore
become the controlling design problem.

This candidate is stored in `examples/cooling_cell_high_cop_candidate.toml`.
It is a robust local candidate from the sampled domain, not proof of a global
COP optimum.

## Triple-speed heat-exchanger scaling

Starting from the verified high-COP candidate, angular speed, both UA values,
and all four effective CdA values were multiplied by three. With each heat-
exchanger gas volume held at 0.66 L, the dimensionless thermal and hydraulic
time scales were preserved. The periodic solution reproduced the same state
trajectory per crank angle: 878.247 J cold heat per cycle, COP 3.384225,
237.83/342.58 K global temperature extrema, 2.263 bar maximum pressure, and
nominal topology. At 27.273 rpm, powers scaled to 399.21 W cold and 117.96 W
mechanical. This is a numerical similarity result and exceeds the former 350 W
cold target within the 150 W modelled human-input limit.

A separate interpretation representing three complete parallel exchanger
modules multiplied UA, CdA, and gas volume by three. With 1.98 L in each heat
exchanger, the periodic iteration reached a non-positive conservative state
during the second cycle. No performance is reported. Tripling physical heat-
exchanger volume changes compression ratios, gas inventory, and cycle state;
it is not dynamically similar to tripling conductance and hydraulic capacity
alone. The failure may also involve numerical positivity and must not be
treated as proof that parallel modules are physically impossible.

The similarity scaling preserves the current temperature chronology; it does
not make the cylinder transfers more isothermal. Improving that behavior at
triple speed would require UA to grow faster than speed, or a change to phase
duration/kinematics, while hydraulic capacity is checked independently.

## Thermally scaled 33-liter trial beyond the former CdA bound

The user explicitly authorized exceeding the former exploratory upper CdA
bound of 1e-3 m2. The 33-liter configuration was changed to 125 W/K for each
heat exchanger and 1.5e-3 m2 for `C -> L`; the other effective areas remained
4.2e-4 m2 for `L-H` and `H -> S`, and 7.5e-4 m2 for `S-C`. This is a targeted
similarity study and does not silently redefine the older global design space.

At one revolution per second, small-to-large swept-volume ratios 0.50, 0.60,
and 0.70 all converged and produced positive cold power. Their respective cold
powers were 609.8, 660.0, and 686.6 W; mechanical inputs were 941.6, 972.6, and
1015.1 W; and cooling COP values were 0.648, 0.679, and 0.676. Ratio 0.60 had
the best COP among these three discrete points, while ratio 0.70 had the
largest cold power. All three event sequences remained non-nominal and failed
the configured validity criteria.

A direct follow-up at ratio 0.60 and 0.151517 revolution per second (9.091 rpm)
was run rather than assuming powers scale linearly with speed. It converged in
five cycles and produced 1930.3 J cold heat per cycle, 292.5 W cold power,
162.6 W mechanical input, and cooling COP 1.799. Slower operation therefore
improved heat transfer per cycle strongly enough that cold power did not fall
in proportion to speed. The valve sequence became nominal and pressure-
equalization error fell to 0.00335. Cold and hot isothermality errors remained
about 0.0559 and 0.0553, just above the configured 0.05 limit, so the point is
still classified invalid and is not yet a feasible design. It nevertheless
shows that the 100 W cold target and human-power scale are no longer separated
by an order of magnitude.

## Equal-Lambda sensitivity at the low-speed 33-liter point

With swept-volume ratio `S/L = 0.60`, speed 0.151517 revolution per second,
125 W/K heat exchangers, and the beyond-bound 1.5e-3 m2 `C -> L` CdA retained,
equal small- and large-cylinder Lambda targets 0.50, 0.60, 0.70, and 0.80 were
tested. All four periodic searches converged in four or five cycles with mass
and energy residual magnitudes below 1.2e-15 and 2.2e-15.

Cold powers increased monotonically through 150.1, 218.7, 292.5, and 368.2 W,
while mechanical inputs increased more rapidly through 54.5, 95.4, 162.6, and
274.6 W. Cooling COP consequently decreased through 2.756, 2.292, 1.799, and
1.341. Maximum temperatures were 361.6, 382.0, 408.7, and 446.5 K, and minimum
temperatures were 219.3, 208.4, 196.2, and 182.3 K. Increasing Lambda therefore
raises cold power but worsens COP and thermal excursions at this operating
point.

Lambda 0.50 produced a non-nominal extra `H -> S` closing and reopening pair
near 180 degrees. Lambda values 0.60 through 0.80 had nominal event ordering.
At Lambda 0.60, the model predicts 218.7 W cold power using 95.4 W mechanical
input, with pressure-equalization and both isothermality indicators below their
configured thresholds. Its validity verdict remains indeterminate because
Mach number is unavailable without actual geometric areas. The very low gas
temperature extrema and constant-property ideal-gas assumption also require
later real-property and moisture checks; this result is not declared a
validated physical design.

Among these discrete points, Lambda 0.50 maximizes COP but is non-nominal,
whereas Lambda 0.60 is the strongest current candidate combining nominal
topology, less than 150 W modelled mechanical input, more than 100 W cold
power, and satisfied currently computable first-level thresholds. No optimum
between the sampled Lambda values is inferred yet.

## Independent large-cylinder Lambda sensitivity

The user identified that increasing only the large-cylinder Lambda could raise
the gas temperature before hot transfer without imposing the corresponding
additional small-cylinder expansion. Independent `small_lambda_target` and
`large_lambda_target` design parameters were therefore added. They are
mutually exclusive with `common_lambda_target` in one design point, preventing
an ambiguous override while preserving the earlier equal-Lambda interface.

At the 33-liter low-speed point with `S/L = 0.60`, Lambda_S fixed at 0.50, and
all other parameters unchanged, Lambda_L values 0.50 and 0.60 converged. The
0.50/0.50 control point reproduced 150.1 W cold power, 54.5 W mechanical input,
and COP 2.756 with a non-nominal brief `H -> S` closing/reopening event near
180 degrees. The decoupled 0.50/0.60 point produced 183.9 W cold power, 69.4 W
mechanical input, and COP 2.651 with nominal valve ordering. Maximum
temperature decreased from 361.6 to 358.8 K; pressure equalization and both
isothermality indicators remained within their configured limits. Validity is
still indeterminate because Mach number is unavailable.

Lambda_L values 0.70 and 0.80 reached non-positive conservative states before
a periodic result and have no reported performance. An intermediate value
0.65 remained exceptionally stiff and was manually interrupted without a
result after several minutes. It is classified as numerically unresolved, not
as a physical failure. The current evidence therefore supports independent
Lambda control and identifies 0.50/0.60 as a substantially better compromise
than the earlier common 0.60/0.60 point, but it does not establish the physical
feasibility boundary above Lambda_L 0.60.

## Heat-exchanger model boundary

The improvised flattened-tube and copper-mesh performance figures discussed
during exploration are not calibration data and are not retained as a
reference design. A mesh-filled flattened tube remains only one unvalidated
candidate architecture.

The first geometry-based implementation uses identical straight rectangular
gas channels in parallel. Geometry, thermal resistances and operating state
remain separate objects. The thermodynamic cycle continues to use its
published lumped `Q_dot = UA * (T_reservoir - T_HX)` closure; the exchanger
module estimates a physically realizable candidate `UA`, gas volume, pressure
drop and Mach number instead of changing that energy balance.

The overall conductance contains distinct gas-film, wall, external-film and
additional contact/fouling resistances. The external boundary is replaceable.
The currently intended cold-cell architecture uses a water-glycol secondary
loop between the water/ice load and the DADA cold exchanger. Ice therefore
does not form directly on the DADA exchanger. The hot side may initially reject
heat to ambient air, but a water-heating boundary can later use the same
gas-side architecture with a different external resistance and reservoir
dynamics.

Gas viscosity, gas thermal conductivity, wall properties, external
conductance, roughness and minor-loss coefficient are mandatory physical
inputs to a study; no default values are interpreted as DADA hardware data.
The example TOML values are illustrative only.

The current channel correlations assume a steady representative mean state,
uniform distribution, fully developed flow and constant properties. The
transitional Reynolds interval is reported as unavailable rather than
interpolated. A mean-density Darcy pressure drop is also unavailable as a
compressible validity claim unless the caller-selected `delta_P / P` limit is
satisfied. Pulsatile phase lag, conjugate wall storage, axial conduction,
manifold maldistribution and compressible duct choking remain outside this
model.

Because a lumped heat-exchanger control volume can accumulate mass, its two
port flows need not be equal instantaneously. The cycle adapter therefore uses
the maximum absolute adjacent-port flow only as a conservative screening
value. It does not infer a unique internal velocity field. Geometry-derived
hydraulic laws must not silently replace the four independent configured CdA
values until a coupled compressible passage model is explicitly adopted and
validated.

Pump power for the glycol circuit is not thermodynamic piston work. It will be
reported separately and included only in a future complete-system mechanical
input balance.

## Scope decisions after the initial cooling-cell study

The cooling cell is a practical reference case, not the intended limit of the
solver. External applications are represented by an energy demand and allowed
duration; no detailed cooler, domestic-hot-water system or future engine-mode
physics is embedded in the first-level thermodynamic kernel. Engine operation
will require a separate physical specification rather than reversing signs in
the refrigerator implementation.

In the absence of prototype access, exchanger transport properties and
correlations may be selected from published literature for screening. Such
values must retain their source, applicable regime and uncertainty and must
not be described as experimental validation of DADA hardware. Oscillatory-flow
literature shows both heat-transfer enhancement and degradation depending on
geometry and boundary conditions, so no generic pulsation correction factor is
adopted.

Air charged from the surroundings may be humid. The first implementation is a
validity screen: the initial relative humidity defines an ideal-mixture water
mole fraction, and each sampled control-volume state is tested against water or
ice saturation pressure. Once the saturation ratio reaches one, the reported
dry-gas solution is potentially affected by condensation or frost. Condensate
mass, latent heat, drainage, blockage and modified flow properties are not
invented; they require a later multi-species, phase-change model.

For the exploratory high-COP cooling-cell case charged at 298.15 K, 100 kPa
and 50 percent relative humidity, the unchanged-vapour screening trajectory is
already supersaturated at `theta = 0` in S, L and C. Maximum reported ratios
are approximately 64.5 in S, 6.26 in L, 9.48 in C and 1.13 in H. Temperatures
at the first saturated samples in S, L and C are below the water triple point,
so the equilibrium warning is frost rather than freely drainable liquid.
These large post-onset ratios are not predicted condensate quantities: the
constant-vapour trajectory ceases to be physical at the first phase change.

Detailed machine construction, strength, friction, linkage reactions and
actuator sizing remain outside the thermodynamic solver. The mechanical
boundary now permits explicitly supplied absolute pressure on the external
face of each piston. It reports both absolute gas force `P_gas * A` and net
pressure force `(P_gas - P_external) * A`; no atmospheric value is assumed
silently.

## Earlier four-bar optimizer reuse boundary

The separate `dada-engine-4bar-optimizer` repository was audited as prior
implementation work, not adopted as a specification. Its planar closure and
link-attached coordinate frames are reusable concepts. Its prescribed waveform
score, 20-percent precompression filter and fixed angular windows are not
physical laws and will not control passive valve events. Its time-reversal
comparison was intended to examine synchronization of the endpoints of a
complete-cylinder transfer. That quantity may be reported, but its physical
importance remains unproven and it is not a global symmetry requirement.

Most importantly, its E family scores rocker angle rather than the final
Cartesian piston drive point, while its F family treats coupler-point X as
piston displacement without an explicit finite connecting rod. The coupled
solver will instead require an explicit slider/follower constraint and score
the resulting cylinder volume law directly. Assembly branch, continuous crank
rotation and singularity margin will be explicit validity diagnostics. The
detailed findings are in `docs/FOUR_BAR_AUDIT.md`.

For the first low-tech search, both four-bar loops share one physical crank
plate and the same crank pin. The common crank radius and phase are therefore
not independent design variables. Each downstream loop retains independent
coupler/rocker dimensions, fixed rocker-pivot location, output support and
slider geometry. A shared shaft with separate crank throws is retained as a
future comparison family if the single-pin constraint proves too costly.

The first implemented output reference uses a finite connecting rod and an
exact prismatic-slider constraint after each rocker. It therefore distinguishes
the rocker attachment-point path from the actual piston displacement. Slider
position and `dV/dtheta` are analytical; volume normalization retains physical
stroke as a separate result. A Chebyshev approximate straight-line output will
be a later interchangeable family, not an implicit correction to this slider
model.

## Exchanger geometry optimization

Exchanger sizing is a constrained problem, not direct inversion of `UA` and
not an unconstrained search for maximum heat transfer. The initial selectable
objectives are minimum internal gas volume, minimum pressure drop, and minimum
gas-side area. Required `UA`, maximum pressure drop, maximum `delta_P / P`,
maximum Mach, maximum gas volume, optional minimum effectiveness, and
correlation availability remain independent hard constraints.

Internal gas volume is the preferred first objective because it changes the
thermodynamic control-volume inventory. A Pareto set in gas volume, pressure
drop and achieved `UA` is retained so that this preference cannot silently
discard hydraulically superior alternatives.

Channel count is a discrete manufacturing variable. Width, height and length
are locally refined inside explicit manufacturing bounds. Deterministic
feasible grid points are required before local refinement; failure to find one
is reported rather than replaced with a penalized but physically unavailable
solution. The current procedure does not claim global optimality.

A final geometry must be iterated with the periodic cycle because its `UA`,
gas volume and hydraulic behavior alter the operating point from which it was
sized. This mutual-consistency loop is not yet allowed to replace configured
cycle `UA`, exchanger volumes, or CdA silently.

The first such fixed-point loop is now implemented. Each iteration runs the
periodic thermodynamic solver, sizes the cold and hot geometries independently,
replaces their lumped gas volumes and overall conductances, and reruns the
cycle. Convergence is based on both screening flows and all four replaced
quantities, followed by a final periodic verification. Flow under-
relaxation may stabilize the geometry search, but conservative thermodynamic
states are never averaged between simulations.

For flow-capacity screening, each side uses its largest absolute adjacent-port
flow together with the minimum heat-exchanger pressure and maximum heat-
exchanger temperature. Combining separate extrema is deliberately
conservative for density, velocity and Mach, but does not represent a
simultaneous resolved duct state.

The loop now has an explicitly selectable geometric hydraulic extension. The
published compressible-orifice model remains the unmodified default and source
reference. In the extension, each heat-exchanger core is split between its two
ports using a configured inlet fraction; the remainder belongs to the outlet.
Collector loss coefficients are independent explicit inputs. The passive
outlet valve retains its own configured CdA in series with its portion of the
core.

The first closure is quasi-steady mean-density Darcy flow, analytically solved
in the laminar regime and iteratively solved in the turbulent regime. Local
orifice resistance is combined through its low-Mach loss coefficient, while
an isentropic sonic cap remains separate. This is not silently described as a
full Fanno or heated compressible-duct solution. Mach and relative pressure
drop remain mandatory validity indicators.

Both ports must carry the cycle-derived screening flow within the configured
pressure-loss bound before a geometry can be injected into a new cycle. The
reference CdA and geometric flow curves can be compared without calibrating
one to the other.

Passage inertia is not yet integrated as a new state. The implemented decision
diagnostic compares acoustic transit time with cycle period and estimates the
pressure amplitude `(L/A) d(m_dot)/dt`. User-configured thresholds decide
whether a dynamic momentum model is recommended. It must be added only when
those indicators, experiments, or event chronology show the quasi-steady
closure to be inadequate.

An exploratory evaluation used the minimum-volume cold geometry found for the
high-COP candidate under illustrative transport and external-resistance data:
200 parallel channels, 40 mm by 1 mm, 116.86 mm long, 0.935 L internal volume,
125 W/K estimated UA, and about 201 Pa core pressure drop at 0.0752 kg/s. The
former 0.66 L volume bound contained no feasible deterministic seed; the best
otherwise admissible point reached only about 117.8 W/K. Expanding the
exploratory bound to 2 L enabled the 0.935 L candidate. None of these values is
validated hardware data.

Using the reference periodic flow chronology with that candidate gave the
following raw quasi-steady indicators. At the cold inlet, maximum Mach was
0.00812 and estimated inertial pressure fraction was 6.73e-5. At the valved
cold outlet, maximum Mach was 0.02514 and inertial pressure fraction was
0.00229. Acoustic transit time divided by period was 5.38e-5 and the core
pressure-drop fraction was 0.00229. With explicitly exploratory 0.01 limits on
the two inertia indicators, a momentum-state model is not currently
recommended.

The full periodic iteration after geometric hydraulic injection did not reach
an outer-loop result within four minutes and was manually interrupted. It is
classified as numerically unresolved, not physically invalid. The controlled
example instead terminates quickly as `cold_exchanger_infeasible` because its
configured valve cannot pass the requested flow within the pressure-loss
limit. These outcomes show that geometric port capacity is now enforced and
that event/integration performance needs further numerical work before the
high-COP candidate can be compared end to end.

The coupled numerical solve now retains a converged periodic state as the
initial guess for the next geometry. When a changed exchanger volume changes
the gas inventory implied by configured charging pressure, all masses and
internal energies are multiplied by the same factor. This matches the new
inventory while preserving composition, temperature and specific energy.

A configurable hydraulic continuation bridges successive closures.
Intermediate configurations interpolate exchanger volume and `UA` and blend
the old and new mass-flow predictions at the same thermodynamic node states.
This blend is only a numerical homotopy, not reportable physical performance.
The final continuation fraction is exactly one and therefore uses the
unmodified geometric closure. Failure at any intermediate fraction is exposed
as a cycle-evaluation failure.

The first end-to-end rerun with four continuation steps still produced no
outer-loop result after approximately three and a half minutes and was
manually interrupted. Warm starting and homotopy therefore improve the solve
architecture but do not yet establish acceptable runtime for the high-COP
case. The next investigation must time and report each continuation fraction
separately and profile the first slow fraction; it must not infer convergence
from the absence of an integration exception.

Live segment and continuation instrumentation subsequently showed that every
25-percent continuation stage completed. With warm starts, the full geometric
cycle/exchanger loop reached its configured `1e-3` fixed-point tolerance in 10
outer iterations. The converged screening result used 0.934879 L and 125 W/K
per exchanger, predicted 123.711 W cooling, 34.854 W mechanical input and COP
3.54939. Cold/hot screening pressure drops were approximately 193/105 Pa and
Mach numbers 0.0242/0.0124. These remain results of the illustrative steady-
correlation geometry model, not validated hardware values.

Because later outer iterations returned exactly the same hydraulic geometry,
repeating all intermediate homotopy fractions was redundant. The solver now
uses one exact final-model solve in that case; it retains the configured number
of fractions whenever the hydraulic network actually changes.

## Published four-bar reference geometries and coupler outputs

Mechanical exploration must begin from reproducible published geometries before
introducing optimization. The first coupler-output reference is the geometry
labelled F65 on the DADA four-bar page: fixed-pivot distance 100, crank 44.96,
coupler 117.91, rocker 141.73, and coupler-local output coordinates
`(94.4, -88.65)`, with B as the local origin and BC as the positive local axis.
These values are dimensionless ratios until an explicit length scale is chosen.

The published page does not specify the piston slider axis, the slider origin,
or the length and assembly branch of a connecting rod between this coupler point
and a piston. Those quantities must therefore remain explicit configuration
inputs and must not be inferred from the animation. The reference regression
test adds a clearly labelled test-fixture slider only to verify the analytic
coupler-point and slider derivatives; it is not a proposed DADA mechanism.

Coupler-output research is retained alongside rocker-output research. A point
fixed to the coupler is evaluated in its complete moving local frame, including
both translation of B and rotation of BC. It then drives the same explicit
finite connecting-rod and prismatic-slider closure as a rocker output. The old
optimizer's projection of point F onto a selected axis remains useful as a
trajectory diagnostic, but is not itself accepted as piston displacement.

The first complete F65 reference uses the historical opposed construction:
both loops share A, D and the same instantaneous crank pin B; the second loop
is the reflection of the first about AD. Because each coupler frame remains
right-handed, reflection changes both the four-bar assembly branch and the sign
of F's local normal coordinate. Both piston axes are initially constrained
parallel to AD. Each axis is centered between the extrema of its F-point
transverse motion. This centering minimizes the largest transverse connecting-
rod angle for the selected axis direction; it is an explicit installation
choice, not a published F65 dimension.

Piston connecting-rod length is parameterized by the ratio of rod length to
the axial projected stroke of F. At a 100 mm F65 ground distance, that projected
stroke is approximately 83.19 mm. Ratios 3, 5, 8 and 12 produce maximum
normalized waveform differences from the direct axial projection of about
2.78%, 1.66%, 1.04% and 0.69%, respectively. The initial reference ratio is 5,
or approximately 416 mm at this scale. It is a compact exploratory compromise,
not a mechanical optimum or a validated DADA dimension. Thermodynamic studies
must use the exact finite-rod displacement, never the projected approximation.

The historical opposed construction already reverses the two slider motions.
The large-cylinder volume convention must therefore not be reversed a second
time. That double reversal was briefly introduced during initial integration
and made both cylinder volumes move nearly together; it was detected from the
volume chronology and removed before any F65 thermodynamic result was accepted.

The first F65 hydraulic integration exposed the non-differentiability of the
published quasi-steady orifice law at zero pressure difference. An optional
pressure regularization multiplies the original magnitude by
`sqrt(delta_p / (delta_p + delta_p_regularization))`. Its configured default is
exactly zero, which reproduces the published closure without modification. The
F65 exploratory file presently uses 1 Pa, explicitly labelled as numerical.
Results must be subjected to decreasing-regularization sensitivity and must
report use of the affected pressure range before they can be accepted.

Regularization removes the zero-pressure square-root singularity but does not
repair the relay behavior of an inertialess valve that jumps instantly between
zero and full CdA. With F65 and the exploratory large valve areas, zero
hysteresis and tested hysteresis scales from 1 Pa through 100 Pa all reached 100
opening/closing events, sometimes within less than one crank degree. This is a
model-closure limitation, not grounds for suppressing repeated events by angle.
No periodic F65 thermodynamic result is currently accepted. The next physical
choice is between a continuously evaluated ideal hydraulic diode and a valve
lift/momentum state with explicit mechanical parameters.

The rocker-output candidate called E10 on the public DADA page became rank E0
after the optimizer scoring change. Both names refer here to the refined
geometry with ground 100, crank 73, coupler 100, rocker 110.3728, and rocker-
local output coordinates `(97.5310203540, -51.6706401020)`. The software uses
`published_e0_opposed` while documentation preserves the E10 provenance.

E0 uses the same explicit finite-rod, parallel-slider, opposed-pair construction
as F65. At ground distance 100 mm its projected output stroke is approximately
165.22 mm. A rod ratio of 5 therefore gives about 826.11 mm, and its maximum
normalized displacement difference from direct projection is about 0.127%.
The nearby opposed extrema are separated by about 37.75 degrees; this is an
observed feature, not rejected by the old unproven synchronization preference.

Crank direction is now explicit and restricted to `-1` or `+1`. With E0 in the
initially implemented direction, the exploratory machine converged in four
cycles without valve chatter but operated as a non-refrigerator: approximately
-557.85 W cold-side heat absorption and 190.20 W mechanical input at the loose
screening tolerances. Reversing the crank to the intended receiver direction
reintroduced the maximum-event-count failure. The stable opposite-direction
result is useful evidence about numerics and directionality, but it is not a
cooling performance result. It reinforces the need for a continuously evaluated
ideal diode or an explicit valve-dynamics closure before receiver-mode four-bar
performance is accepted.

The reflected F65 pair already produces opposed slider motions. The cylinder
volume convention must therefore have the same sign for both sliders. An early
adapter reversed the large-cylinder convention a second time, causing the two
volumes to move almost together; a regression test now requires the
small-cylinder minimum and large-cylinder maximum to occur within 20 degrees.

The first thermodynamic F65 integration remains numerically unresolved. The
published geometry drives continuously through pressure equality, exposing the
non-differentiable zero-pressure-drop limit of the quasi-steady ideal-orifice
law in both permanently bidirectional links. Radau accumulated tens of
thousands of right-hand-side evaluations over about 1.6 crank degrees. A tiny
explicit valve hysteresis (opening at 1 Pa and closing at 0.5 Pa) removes the
same exact-zero condition at the check valves but cannot regularize S-C and
L-H. No converged F65 thermal performance is claimed. Introducing a configured
low-pressure-drop flow regularization or a viscous transition changes the
hydraulic closure and requires an explicit physical decision before use.
## Continuous ideal-diode valve closure

The selectable `continuous_ideal_diode` closure removes persistent valve
memory. Forward flow follows the configured orifice model whenever upstream
pressure exceeds downstream pressure and is exactly zero otherwise. Zero-
pressure crossings remain non-terminal integration events used only to report
opening and closing chronology. The former `discrete_hysteretic` closure remains
the default for all existing configurations.

With E0 in receiver direction, LSODA, and a 1 Pa orifice regularization, the
periodic state converged in four cycles. It predicted approximately 177.85 W
cooling, 88.11 W mechanical input, and cooling COP 2.019. Global mass and energy
residuals were approximately `8.3e-17 kg` and `8.0e-13 J`. This screening point
is still invalid under its configured criteria because cold isothermality error
was 5.71% against a 5% limit, Mach number was unavailable, and wet-air screening
predicted frost or condensation.

A decreasing regularization study using 0.30, 0.10 and 0.03 Pa converged in four
cycles for every value. Cooling powers were 177.85405, 177.85444 and 177.85357 W;
mechanical powers were 88.09635, 88.09262 and 88.09167 W; COP values were
2.018858, 2.018948 and 2.018960. Integrated performance is therefore insensitive
over this range to much better than 0.01%. The fraction of stored samples with
at least one of the four pressure differences inside the regularization width
fell from about 52.4% to 30.5% and 18.5%. The exact number of reported zero
crossings did not converge monotonically, so fine valve chronology near pressure
equality remains a numerical diagnostic rather than a validated physical lift
history. The E0 exploratory reference uses 0.03 Pa.

F65 was subsequently evaluated with the same receiver direction, continuous
diodes, 0.03 Pa regularization and thermodynamic machine parameters. It also
converged in four cycles and predicted approximately 288.67 W cooling,
225.90 W mechanical input and COP 1.278. Its maximum pressure and temperature
were about 5.98 bar and 448 K, compared with 3.12 bar and 376 K for E0. Both
cases remain invalid screening results, but the comparison supports the
qualitative interpretation that F65 produces more cooling through a more
violent, less efficient thermodynamic excursion.

## First local shared-crank mechanical search

The first thermodynamic search retains one common crank pin, parallel piston
axes, fixed cylinder swept volumes, fixed exchanger and hydraulic hardware, and
the E0 total gas inventory of `0.0395562623435 kg`. Candidate filling pressure
is consequently derived from this fixed mass; it is not reset to 1 bar for each
kinematic waveform. The public `shared_crank_rocker` configuration permits the
two four-bar loops to have independent coupler, rocker, pivot and output-point
coordinates, although this first search deliberately retained reflected loops.

The phase offset was removed from the optimization variables after a +/-3 degree
check changed periodic performance only at numerical-noise scale. For a periodic
autonomous cycle it merely selects the reported angular origin. Local searches
varied dimensionless crank, coupler and rocker lengths and rocker-output angle.
The selected exploratory candidate has crank/ground `0.40`, coupler/ground
`1.00`, rocker/ground `1.17`, and reflected output angles +/-56 degrees. At the
reference numerical tolerances it predicts `100.123 W` cooling, `29.804 W`
mechanical input and cooling COP `3.3594`, versus `177.854 W`, `88.091 W` and
COP `2.0190` for E0.

This candidate is not a global or mechanical optimum. Holding swept volume
fixed while shortening the physical slider stroke implicitly enlarges piston
area. With ground distance 100 mm, its stroke is about 90.49 mm and its implied
small and large piston diameters are about 625 and 681 mm. At 101325 Pa external
pressure the predicted net pressure-force ranges are approximately -8.8 to
22.5 kN and -10.7 to 27.2 kN. These large dimensions and forces must become
explicit constraints before further maximization of COP; otherwise the search
can exploit an unconstrained piston diameter. Mach number remains unavailable,
wet-air screening predicts phase change, and the passively observed valve-event
sequence is non-nominal. The thermodynamic criteria that are available pass,
so the overall validity verdict is `indeterminate`, never `valid`.

The 100 mm ground distance above was only the normalization scale inherited
from the published four-bar proportions, not a required physical frame size.
Uniformly scaling every link preserves the normalized imposed volume laws and
therefore leaves the present thermodynamic solution unchanged. A physical scale
must instead be selected from a piston construction criterion. For candidate
001, imposing a large-piston bore/stroke ratio of 1 gives a ground distance of
about 384 mm, a common stroke of 348 mm, a large bore of 348 mm and a small bore
of 319 mm. This removes the artificial 625--681 mm bores caused by interpreting
the 100 mm normalization as hardware. The selected bore/stroke ratio remains a
configurable engineering choice, not a new physical law or an optimum.

## First-level capacity-speed similarity

For fixed dimensionless kinematics and intensive states, the first-level model
has an exact two-factor similarity. Scaling every control-volume volume and the
gas inventory by capacity factor `s`, cycle frequency by `f`, and every UA and
CdA by `s*f` preserves pressure and temperature histories versus crank angle.
Heat, work and mass transported per cycle scale by `s`; powers and mass-flow
rates scale by `s*f`; COP is unchanged. Mechanical link lengths scale by the
cube root of `s` if piston proportions are to remain geometrically similar.

This is a property of the lumped model, not proof that a manufacturable heat
exchanger has UA proportional to gas volume. Applying it to candidate 001 and
the 20-minute ideal water load requires `s*f = 3.47158`: approximately
`347.58 W` cooling and `103.47 W` thermodynamic mechanical input at unchanged
COP 3.3594. It simultaneously requires about `434 W/K` on each side and CdA
values 3.47158 times larger. Those exchanger and hydraulic requirements must be
realized geometrically before the scaled point can be called feasible.

The user selected the two-capacity variant as the next cooling-cell reference.
`cooling_cell_x2_reference.toml` doubles both cylinder and exchanger volumes,
uses a 1 bar filling pressure, and scales the physical linkage by cube root two.
Its independently integrated periodic cycle predicts 347.585 W cooling,
110.212 W thermodynamic input and COP 3.15378 at 17.897 rpm. The configured
0.8 litre exchanger volumes and 226.0 W/K conductances are targets for geometric
closure, not validated components. Available thermodynamic validity criteria
pass, while Mach remains unavailable, valve topology is non-nominal and wet-air
screening predicts phase change.

Instantaneous generalized gas torque uses absolute cylinder pressures. When
external piston pressure is configured, the mechanical boundary now also
reports generalized net pressure torque by subtracting `P_external*dV/dtheta`
for each piston. This subtraction changes instantaneous loads even though a
constant atmospheric pressure performs zero net work over a closed volume
cycle. Gas-only and net-pressure torque are retained as separate quantities.

Reservoir temperatures are now explicit sizing variables. Increasing local
temperature difference may reduce required UA, but changing reservoir
temperature also changes the thermodynamic lift and therefore the entire cycle.
At the 27.85 rpm, 0.4 litre-HX candidate with 383 W/K, lowering the cold
reservoir from 268.15 to 263.15 K reduced cooling from about 347.6 to 273.7 W
and COP from 3.348 to 3.110. A larger local exchanger difference must therefore
never be credited without rerunning the periodic thermodynamic state.

A speed-UA sensitivity at 268.15 K found a more compact thermal point near
225 W/K and 35.6 rpm, producing about 346.0 W at COP 3.154 and 109.7 W
thermodynamic input. At 150 W/K and 61.3 rpm the model reaches about 350.2 W,
but COP falls to 2.681 and input rises to 130.6 W. Both also require larger
hydraulic CdA with speed. These are Pareto directions, not final selections;
auxiliary and geometric hydraulic constraints decide between them.

## Industrial refrigerator working-fluid search

The domestic refrigerator/freezer case is a new constrained sizing problem;
the human-powered geometry is evaluated only as an off-design seed. Industrial
screening may use helium and elevated absolute charge pressure. Pressure is
treated as a power-density variable, not a direct COP improvement: for the
ideal-gas model an otherwise similar cycle at proportionally scaled UA has the
same COP. Helium changes the physics through its configured `R`, heat capacities
and gamma, while its transport benefit enters only through exchanger geometry.

The first helium freezer screening uses rounded NIST-based calorically-perfect
properties `R=2077.1`, `Cp=5193.0`, and `Cv=3115.9 J/(kg K)`. At reservoir
temperatures 250.15 and 310.15 K, the inherited candidate produces 77.92 W cold,
43.59 W input and COP 1.788 at one bar. Scaling pressure and UA to five bar gives
about 389.6 W and 217.9 W with the same COP. No high-pressure result is declared
valid until the NIST real-gas equation of state or an equivalent source supplies
the compressibility and property-variation checks.

Argon is retained as a monatomic comparison fluid with rounded screening
properties `R=208.13`, `Cp=520.33`, and `Cv=312.20 J/(kg K)`. Helium and argon
have the same first-level intensive ideal cycle when pressure, geometry, speed,
UA and normalized hydraulic time scales are matched. Matching the latter
requires `CdA_target/CdA_reference = sqrt(R_reference/R_target)`. A direct
numerical check gives COP 1.49225 for helium and 1.49229 for argon; the tiny
difference is due to rounded constants and tolerances. Fluid ranking must
therefore include transport properties, Mach margin, containment and real-gas
validity. The ideal solver must not claim an intrinsic helium COP advantage.

Commercial comparison uses an explicitly matched boundary. The approximate
helium COP 1.8 is already in the range of present domestic compressor catalog
data, but the DADA value currently excludes actuator, transmission, bearings,
fans, pumps and controls. It is encouraging but cannot be compared directly to
an appliance label or treated as a commercial-system COP.

Domestic-appliance optimization is constrained by annual electricity rather
than an arbitrary continuous cooling power. Annualization uses explicit on- and
off-cycle powers and an explicit duty fraction. The current 200 kWh/year point,
85 percent transmission efficiency and auxiliary powers are labelled
exploratory; they are not inferred physical properties of an unspecified
cabinet.

The first helium freezer search retains `S/L=0.90`, 2 rpm and 35 W/K per
exchanger as candidate 001. A pressure-volume inverse similarity is now an
explicit solver operation: multiplying all volumes and CdA values by `k` while
dividing filling pressure and pressure thresholds by `k` preserves the ideal
cycle at constant inventory and UA. For `k=0.2`, candidate 001 uses 5 bar,
5.94/6.60 litre swept volumes and 0.132 litre exchanger gas volumes. This
similarity is not evidence of real-gas, exchanger, seal, force or structural
validity at elevated pressure.

## Selected Doty calculation reference

The user selected the Doty fallback as the current working basis. The first
whole-bank scaling implementation and its limitations are documented in
[DOTY_SCREENING.md](DOTY_SCREENING.md). This adds a separate steady screening
calculation, not a change to the conservative cycle equations.

## Four-bar motor demonstrator seed

The user requested four-bar kinematics for the demonstrator. Reuse the existing
shared-crank rocker model and finite slider rods, with automatic motor reversal.
The initial comparison holds hardware and gas inventory fixed across piecewise-linear and
four-bar waveforms at 2, 5 and 10 Hz. The 1 litre seed and exchanger parameters
are exploratory, not selected hardware. No new thermodynamic equation or Doty
closure is introduced. See [MOTOR_FOUR_BAR.md](MOTOR_FOUR_BAR.md).
