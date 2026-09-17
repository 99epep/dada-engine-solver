# Four-stage motor optimization: evidence and current design direction

## Scope and evidence snapshot

Documentation audit: 2026-09-17. Numerical authority is the committed output
artifacts at repository commit `9656fff7199d4c71a129b30ea93aeaa2f3f210fe`,
not earlier prose or a new simulation. In particular, the production 9D report
and history at that commit contain **216 attempts**. The working directory had
newer, uncommitted campaign progress during this audit; it was not overwritten
or silently substituted for this reproducible checkpoint. Read the latest
committed report/history before treating this snapshot as a selection.

This document separates historical constant-property / legacy-gas-model results,
production variable-property gas-model results, controlled diagnostics and
optimization campaigns. All efficiencies are **indicated gas work divided by
external-source heat input**, not shaft efficiency. Wall-to-gas phase heat is a
different boundary. External-air fan losses are excluded by project decision,
not physically zero. Mechanical losses are unknown. All searches provide
best-found evidence under finite budgets, not global optimality.

## Implemented shared four-stage law

[FourStageVolumeKinematics](../src/dada_solver/four_stage_kinematics.py) has seven
free parameters: `t1, t2, t3, a_l, b_l, a_s, b_s`. Forward motor-cycle fraction
runs from T0=0 through T1=t1, T2=t2, T3=t3 to T4=1 (the next T0), with
`0 < t1 < t2 < t3 < 1`. The conventional stage names are low-pressure exchange,
compression, high-pressure exchange and expansion, in that order. They are
kinematic labels, not imposed pressure plateaus, scheduled valve openings or
adiabatic boundary conditions.

| Knot | S normalized swept volume | L normalized swept volume |
| --- | ---: | ---: |
| T0 | b_s | 1 |
| T1 | 1 | b_l |
| T2 | a_s | 0 |
| T3 | 0 | a_l |
| T4 | b_s | 1 |

Each intermediate level independently lies in [0,1]. Zero denotes minimum
**enclosed** volume including clearance; `V = V_min + q * V_swept`.
Interpolation is piecewise linear, velocity is piecewise constant, and velocity
jumps instantaneously at the boundaries (including the cycle wrap). Breakpoints
are supplied to the integrator. No acceleration, inertia, stress or mechanical-
loss penalty is attached to these discontinuities. The campaign minimum stage
duration of 0.02 cycle (7.2 deg) is only a numerical/search bound; the class
itself requires ordered knots. It is not a mechanical-feasibility guarantee.
This is a thermodynamic design probe, mechanically unrealizable as written.

The object adapts forward motor time to study angle; the factory still applies
its usual motor reversal exactly once. Do not add another phase shift or
reversal when injecting it. The later target is a smooth C2 law, followed by
mechanism synthesis and thermodynamic evaluation of the actual mechanism; see
[MOTION_OPTIMALITY.md](MOTION_OPTIMALITY.md).

## Historical legacy-gas-model experiments

All results in this section use the earlier constant-property internal gas film
and hydraulic screening. They retain dynamic wall storage and external-air
coupling. Their historical laminar-only feasibility guard must not be copied
into the current production campaign. Sources are 298.15/598.15 K at 2 Hz.
The fixed-inventory and candidate-specific filling comparisons are distinguished
below; even within one exchanger generation, those are different experiments.

### A. Shared timing: seven-variable search and fine refinement

The [initial experiment](FOUR_STAGE_K2_SEARCH.md) preserves the first 12-candidate
record. The committed continuation and two refinements supersede that initial
performance checkpoint:

| Campaign report | Evaluations in that history | Best indicated efficiency | Indicated power (W) |
| --- | ---: | ---: | ---: |
| [motor_four_stage_k2/report.json](../outputs/motor_four_stage_k2/report.json) | 805 | 21.73617456% | 51.34298375 |
| [motor_four_stage_k2_refine/report.json](../outputs/motor_four_stage_k2_refine/report.json) | 256 | 21.93017374% | 52.37712841 |
| [motor_four_stage_k2_refine_final/report.json](../outputs/motor_four_stage_k2_refine_final/report.json) | 128 | 21.95670623% | 52.68163065 |

These three searches hold the K2 gas inventory fixed. The final fine search uses
half-widths 0.0025 then 0.001 in the seven dimensionless coordinates; its best
record is index 100. The selected shared chronology and intermediate levels are
shown below alongside the later balanced basin for reproducibility. Digits
identify the numerical candidates, not physical measurement precision.

| Parameter | Final shared 7D candidate | Independent 150/30-start candidate |
| --- | ---: | ---: |
| `t1` | 0.5395244395962604 | 0.4639541340225067 |
| `t2` | 0.6990599161470307 | 0.6192082431707531 |
| `t3` | 0.9231798069375672 | 0.8949516786045083 |
| `a_l` | 0.36318997017010785 | 0.3596051765014019 |
| `b_l` | 0.11296032485372168 | 0.11579141985445066 |
| `a_s` | 0.36025860498579815 | 0.35595731584136037 |
| `b_s` | 0.17851659838446757 | 0.14961511405242023 |

| Stage | Final shared 7D duration (deg) | Independent 150/30-start duration (deg) |
| --- | ---: | ---: |
| low pressure exchange | 194.22879825465375 | 167.0234882481024 |
| compression | 57.432771558277295 | 55.8914792933687 |
| high pressure exchange | 80.68316068459315 | 99.26763675615189 |
| expansion | 27.655269502475814 | 37.81739570237701 |

The shared durations are differences of the recorded knots, not four additional
variables. Sources: [motor_four_stage_k2_refine_final/report.json](../outputs/motor_four_stage_k2_refine_final/report.json),
[motor_four_stage_11d/report.json](../outputs/motor_four_stage_11d/report.json) and
[motor_four_stage_from_150_30/report.json](../outputs/motor_four_stage_from_150_30/report.json).

A necessary charge-policy bridge is
[optimize_motor_four_stage_thermo3d.py](../examples/optimize_motor_four_stage_thermo3d.py): with motion frozen, it searches
H_i/H_o tube lengths and S/L swept-volume ratio at fixed total swept volume,
refilling each geometry to 100000 Pa absolute at theta=0. In
[motor_four_stage_thermo3d/report.json](../outputs/motor_four_stage_thermo3d/report.json), the 161-attempt search retains its
control: **22.02886011%**, **51.57574235 W**. That is a reevaluation with the
100 kPa filling convention, not a motion gain over the fixed-mass 21.95670623%
result. Subsequent experiments below use candidate-specific 100 kPa filling.

### B. Independent S/L timing: 11D search

[optimize_motor_four_stage_11d.py](../examples/optimize_motor_four_stage_11d.py) releases the small-cylinder origin,
three timing knots per cylinder and four intermediate levels. Large-cylinder
T0 fixes the common angular gauge. There is no symmetry penalty.
[motor_four_stage_11d/report.json](../outputs/motor_four_stage_11d/report.json) and its
[history](../outputs/motor_four_stage_11d/history.jsonl) contain 167 attempts
(164 converged, 3 interrupted), reaching radius 0.05 of the configured
0.10 / 0.05 / 0.025 schedule. The exact shared seed remains best at
**22.02880543%**, **51.57569136 W**.

This broader search found no improvement. Shared timing appeared an excellent
reduced approximation, but the finite budget and coarse timing resolution do
not prove exact symmetry or exclude a small nearby gain.

### C. Relative event offsets: targeted 4D test

[optimize_motor_four_stage_symmetry4d.py](../examples/optimize_motor_four_stage_symmetry4d.py) freezes the levels, hardware
and paired mean T1/T2/T3 times, varying only the four absolute S-minus-L event
separations. Signed axial probes precede local 4D Sobol sampling.
[motor_four_stage_symmetry4d/report.json](../outputs/motor_four_stage_symmetry4d/report.json) contains 121 attempts; best
index 91 gives **22.05632008%**, **51.40372181 W**.

Its recorded event offsets (deg) are
`(1.3166006207466125, -0.7451270595192909, -0.06738369166851044, 0.18535936251282692)`;
the event-separation RMS is **0.762815538736638 deg**. The gain is about
0.0275 percentage point against the same-model shared control. Thus exact shared
timing is not an exact optimum in this experiment, while the gain from releasing
it is very small. Shared timing remains useful for low-dimensional design work.

### D. Phase heat and hydraulic diagnostics

[analyze_motor_four_stage_phase_losses.py](../examples/analyze_motor_four_stage_phase_losses.py) exports
[motor_four_stage_phase_diagnostics.json](../outputs/motor_four_stage_phase_diagnostics.json) and the corresponding
[CSV](../outputs/motor_four_stage_phase_diagnostics.csv). Its hydraulic quantity
is `integral(sum(abs(dp) * abs(m_dot) / rho_upstream) dt)`, a screening dissipation
proxy, **not** full compressible exergy destruction. The following entries are
rounded directly from that report, not a reintegration:

| Phase | H_i wall-to-gas heat (J) | H_o wall-to-gas heat (J) | Hydraulic proxy (J) |
| --- | ---: | ---: | ---: |
| Low-pressure exchange | 0.322831 | -65.623528 | 1.259471 |
| Compression | 0.612462 | -27.641088 | 0.020198 |
| High-pressure exchange | 100.717600 | -8.470190 | 0.970983 |
| Expansion | 15.415303 | 10.476941 | 0.037740 |

Hydraulic losses are concentrated in the exchange stages; compression and
expansion contributions are much smaller. Thermal coupling, however, remains
active throughout all four stages. Compression and expansion are not isolated
adiabatic periods. The sampled phase quadratures are diagnostics; use the
integrated external-source heat/work outputs for efficiency.

### E. Artificial thermal and hydraulic saturation

[analyze_motor_exchanger_saturation.py](../examples/analyze_motor_exchanger_saturation.py) and
[motor_exchanger_saturation/report.json](../outputs/motor_exchanger_saturation/report.json) isolate sensitivities at frozen
motion and geometry. Thermal factors multiply **both gas-wall and air-wall
conductances** while holding gas volume, wall capacity and hydraulics fixed.
Hydraulic factors divide linear and quadratic loss coefficients while retaining
thermal conductances, volumes and wall capacity (the flow-area cap is unchanged).

The baseline is 22.02876108%. At thermal factor 2, H_i alone gives 22.50630043%,
H_o alone 23.16208915%; at factor 8 both give 24.50875138%. This frozen system
was not thermally saturated, with H_o especially sensitive. Hydraulic benefits
show diminishing returns: H_i factors 4/8 give 22.26275370/22.26451136%, whereas
H_o gives 22.79762300/22.86202198%. H_i reaches a plateau earlier.

These are controlled zero-cost conductance experiments, **not realizable
exchanger resizing** and not optimization campaigns. They do not mean that a
physical exchanger should be enlarged by the same factor.

### F. Physical H_o enlargement and 120/60 versus 150/30 timing

[evaluate_motor_four_stage_120deg_hx.py](../examples/evaluate_motor_four_stage_120deg_hx.py) performs a 2x2 comparison of
retained versus 120/60/120/60 deg timing and 3708 versus 6050 H_o tubes. H_i is
unchanged. H_o length and diameter stay fixed; the enlarged outlet-valve CdA is
scaled by the retained LP duration divided by 120 deg.
[evaluate_motor_four_stage_150_30_hx.py](../examples/evaluate_motor_four_stage_150_30_hx.py) changes timing to
150/30/150/30 deg using those **same two hardware assemblies**, including the
previous enlarged valve. It does not resize H_o for 150 deg.

| Artifact / case | Timing (deg) | H_o tubes | Indicated efficiency | Indicated power (W) |
| --- | --- | ---: | ---: | ---: |
| [motor_four_stage_120deg_hx/report.json](../outputs/motor_four_stage_120deg_hx/report.json) / A | retained shared 7D | 3708 | 22.02887735% | 51.57574909 |
| [motor_four_stage_120deg_hx/report.json](../outputs/motor_four_stage_120deg_hx/report.json) / B | 120/60/120/60 | 3708 | 21.35625764% | 49.65685667 |
| [motor_four_stage_120deg_hx/report.json](../outputs/motor_four_stage_120deg_hx/report.json) / C | retained shared 7D | 6050 | 19.47430173% | 54.11445679 |
| [motor_four_stage_120deg_hx/report.json](../outputs/motor_four_stage_120deg_hx/report.json) / D | 120/60/120/60 | 6050 | 19.43662392% | 54.55389320 |
| [motor_four_stage_150_30_hx/report.json](../outputs/motor_four_stage_150_30_hx/report.json) / E | 150/30/150/30 | 3708 | 21.36449577% | 51.52321495 |
| [motor_four_stage_150_30_hx/report.json](../outputs/motor_four_stage_150_30_hx/report.json) / F | 150/30/150/30 | 6050 | 19.34048344% | 55.35077544 |

Physical enlargement increases internal surface, flow area, gas hold-up,
headers and wall thermal capacity together. In this report H_o total gas volume
changes from 3.3246104237732764e-5 to 5.3061304783787285e-5 m^3, and wall capacity
from 88.30246997753838 to 142.81174308632885 J/K. It does not reproduce a zero-cost
thermal multiplier. The larger H_o reduces the hydraulic difficulty of shorter
LP exchange, but its coupled storage/volume penalties can lower efficiency.
These tests do not isolate the contribution of each penalty individually.

The later 9D report also retains a 4840-tube H_o, 150/30/150/30 control at
20.54130196%, between the 6050-tube and reference-size results. This is a reduced
enlargement, still larger than 3708 tubes; its valve follows the 9D tube-count
scaling policy. Do not silently substitute that policy for the 120-degree valve
rule. See the `controls` field of
[motor_four_stage_hx9d/report.json](../outputs/motor_four_stage_hx9d/report.json).

### G. Counterfactual phase heat ablation

[analyze_motor_phase_heat_ablation.py](../examples/analyze_motor_phase_heat_ablation.py) and
[motor_phase_heat_ablation/report.json](../outputs/motor_phase_heat_ablation/report.json) suppress selected gas-wall heat
paths during compression and/or expansion. External-air/wall coupling remains
active; the wall stores energy and the full periodic wall convergence test is
retained. The two frozen hardware cases use the retained shared chronology.

All tested ablations harmed those frozen configurations. Blocking both exchangers
during compression and expansion gives 16.20401724% instead of 22.02887735% for
the reference hardware, and 14.30346838% instead of 19.47430173% for enlarged H_o.
This is **local causal sensitivity**, not proof that compression/expansion heat
is globally beneficial in every optimized machine. A redesigned machine could
choose less compression rather than relying on cooling to reduce compression
work. No reoptimization in this experiment tests that alternative.

### H. First physical 9D search

[optimize_motor_four_stage_hx9d.py](../examples/optimize_motor_four_stage_hx9d.py) varies
`t1,t2,t3,a_l,b_l,a_s,b_s,n_i,n_o` with physical tube counts, retained lengths,
count-scaled outlet CdA and candidate-specific filling. In
[motor_four_stage_hx9d/report.json](../outputs/motor_four_stage_hx9d/report.json) and its
[history](../outputs/motor_four_stage_hx9d/history.jsonl), 508 attempts comprise
4 controls, 132 broad global Sobol proposals and 372 local proposals; 500
converged and 8 were interrupted. The exact old control remains best at
22.02887735%, 51.57574909 W and 3708/3708 tubes. Best new local index 151 reaches
21.97363023% with 3804/3696 tubes.

This does not establish a true 9D optimum at equal tube counts or at the old
chronology. Good candidates occupy a narrow region: broad Sobol was inefficient,
and even the final local radius 0.03 remained much larger than the 0.001 fine
7D radius. This motivated deliberately local refinement rather than interpreting
an unchanged incumbent as proof of optimal hardware.

### I. Independent climb from 150/30/150/30 only

[refine_motor_four_stage_from_150_30.py](../examples/refine_motor_four_stage_from_150_30.py) uses that chronology as its
only initial search center, with K2 hardware frozen at 3708/3708 tubes. It does
inherit the reference intermediate levels as a seed, but does not use the old
champion timing as a center or its state as a warm start. Seven-variable
incumbent-centered radii decrease from (0.15 timing, 0.06 levels) to
(0.001, 0.0005). All 512 recorded attempts are in
[motor_four_stage_from_150_30/report.json](../outputs/motor_four_stage_from_150_30/report.json) and its history.

The best, index 492, reaches **22.05698563%**, **51.08267642 W**, with the exact
parameters/durations in section A: approximately 167/56/99/38 deg. It did not
return to the old 194/57/81/28 deg basin. A more balanced LP/HP timing basin exists;
the old exchange asymmetry was not uniquely selected by the finite search.
The modest improvement is a historical legacy-model result, not a direct
comparison against the later production physics.

## Production variable-property model and current local 9D campaign

The implemented internal-gas closure uses variable transport properties,
compressible pressure-squared tube hydraulics, instantaneous flow-dependent gas
film, laminar thermal entry and model-domain diagnostics. Stagnant radial
screening remains at zero flow; no empirical pulse multiplier is applied.
See [MICROTUBE_GAS_MODEL.md](MICROTUBE_GAS_MODEL.md) for equations and limitations
and [EXCHANGER_VALIDATION.md](EXCHANGER_VALIDATION.md) for validation evidence.
Static UA is only reference metadata in this mode. External-air modelling is
still screening and remains a major hardware-design uncertainty.

[refine_motor_four_stage_hx9d_variable_gas.py](../examples/refine_motor_four_stage_hx9d_variable_gas.py) starts from the balanced
7D **motion**, reevaluated with production physics; legacy states/performance are
not reused. Its nine coordinates are the seven shared motion parameters and
integer tube counts `n_i,n_o` in [2400,5200]. H_i/H_o lengths remain
0.0463173095241189 / 0.02519308603379876 m and tube ID remains 0.00033 m.
Tube count changes surface, flow area, headers/dead volume and wall capacity;
outlet-valve CdA scales linearly with count. Each candidate is filled at
100000 Pa at theta=0. This campaign keeps cylinder dimensions and speed fixed.

There is **no global Sobol branch**. Local scrambled Sobol proposals perturb the
incumbent, using timing/level/count half-widths
`(0.02,0.01,500)`, `(0.01,0.005,300)`, `(0.005,0.0025,150)`,
`(0.0025,0.00125,80)`, `(0.001,0.0005,40)`, `(0.0005,0.00025,20)`;
48 evaluations are assigned to each radius level.

Feasibility uses the current thermodynamic and microtube **validity verdict**,
plus the campaign guards: 40 W indicated power, 1200000 Pa, 850 K and
0.08 kg/s. Do **not** add the historical blanket `Re < 2300` screen. The model
handles supported laminar states, rejects unsupported transition states, and
can use its supported turbulent closure within its full declared domain.
This is model applicability, not experimental transient validation.

### Best found so far at the committed checkpoint

[motor_four_stage_hx9d_variable_gas/report.json](../outputs/motor_four_stage_hx9d_variable_gas/report.json) and
[history.jsonl](../outputs/motor_four_stage_hx9d_variable_gas/history.jsonl)
contain 216 attempts: 173 converged, 26 invalid-exchanger results and 17
interruptions. The requested phase budget was 25200 s; reported elapsed time
was 25200.096041794983 s. Radius index 4 (zero-based) reached
`(0.001,0.0005,40)`; the full decreasing-radius schedule was **incomplete**.
This is not a final selected machine.

| Production evaluation | Indicated efficiency | Indicated power (W) |
| --- | ---: | ---: |
| Balanced-motion seed, index 0 | 21.95765114% | 52.35746491 |
| Best found so far, index 152 | 22.05791992% | 52.17718204 |

The within-production gain is 0.10026878 percentage point. Do not interpret the
small difference from the legacy balanced candidate as an optimization gain:
changing the gas closure alone moved that seed from 22.05698563% to 21.95765114%.
Old candidates require same-model reevaluation before ranking them here.

| Coordinate | Best found so far |
| --- | ---: |
| `t1` | 0.47924609531431156 |
| `t2` | 0.6064965825648978 |
| `t3` | 0.9024383863589415 |
| `a_l` | 0.3619805269073296 |
| `b_l` | 0.1173749422869344 |
| `a_s` | 0.3554084146606237 |
| `b_s` | 0.1541981502721214 |
| `n_i` | 3222 |
| `n_o` | 3778 |

| Stage | Duration (deg) |
| --- | ---: |
| low pressure exchange | 172.52859431315215 |
| compression | 45.810175410211045 |
| high pressure exchange | 106.53904936585572 |
| expansion | 35.12218091078107 |

The best uses fewer H_i tubes and slightly more H_o tubes than the 3708 reference.
Both exchange stages lengthen relative to the balanced seed, mainly at the
expense of compression. Compression and expansion remain tens of degrees,
well above the 7.2 deg search floor.

Recorded external heat input is 236.54624839296267 W, gas inventory
0.0014086977362127346 kg, maximum gas pressure 500700.7349374512 Pa and maximum
gas temperature 584.310676756102 K. H_i/H_o gas hold-up is
3.5070244882442355e-5 / 3.377818751622287e-5 m^3 and wall capacity is
139.87038566150594 / 89.93169675704962 J/K. These quantities are coupled outputs
of geometry, not independent optimization coordinates.

## Next design study and interpretation limits

The approved next temperature study adapts the **whole machine**, first hardware
and cylinder ratio, then motion, at each source temperature. It retains fixed
total swept volume and candidate-specific atmospheric filling, not fixed mass.
The precise temperature grid, charge rule, power-floor policy and mandatory
Carnot comparisons are in [MOTOR_RESEARCH_OBJECTIVES.md](MOTOR_RESEARCH_OBJECTIVES.md).
This is a plan, not a temperature campaign completed by this documentation audit.

The durable exchanger lessons and remaining external-air uncertainty are in
[EXCHANGER_NEXT_STEPS.md](EXCHANGER_NEXT_STEPS.md). The later path is selected
linear target -> C2 smoothing and limited thermodynamic readjustment -> six-bar
synthesis -> evaluation of actual six-bar motion. Existing free-target six-bar
work does not mean that this new four-stage-to-C2-to-mechanism sequence is done.
