# Motor research objectives

Updated 2026-09-17. This document separates completed numerical exploration
from the agreed next design study and longer-term objectives. It does not
instruct an unrestricted search. The report-based record of completed and
in-progress work is [FOUR_STAGE_OPTIMIZATION.md](FOUR_STAGE_OPTIMIZATION.md).

## 1. Search for optimal motion laws

Maximize efficiency subject to specified useful power, frequency, size and
hardware applicability constraints. Search the small- and large-cylinder laws
independently: timing, dwell, transfer duration and waveform shape need not be
mirrored. The current piecewise-linear law is a reproducible reference, not an
established ideal, a performance ceiling or a required optimum. The desired
cycle chronology is distinct from any particular piston motion law.

Explore smooth periodic laws with explicit stroke, velocity, acceleration and
mechanical feasibility bounds. Report which bounds and losses are modeled;
without them an apparent optimum could require unrealizable motion. Keep
valves passive and measure the resulting chronology rather than schedule it.

Compare symmetric and independently parameterized four-bar mechanisms under
the same thermal, hydraulic, inventory and size assumptions. Mirror symmetry
was a refrigerator construction choice, not a motor requirement. The existing
shared-crank configuration already supports distinct loop parameters; retain
common-crank geometric compatibility where a shared crank is used. Independent
loop parameters do not imply two independent crank speeds. Do not constrain
all candidate waveforms to what the present four-bar seed can reproduce.

Shared timing, independent S/L timing, small relative event offsets and a second
timing basin have now been explored. All remain best-found evidence. The
current linear law has instantaneous velocity jumps, with no acceleration or
mechanical-loss penalty; the 2% stage floor is only a search bound. The agreed
path after the linear machine/temperature exploration is C2 smoothing, limited
thermodynamic reoptimization, six-bar synthesis, then thermodynamic evaluation
of the actual mechanism. See [MOTION_OPTIMALITY.md](MOTION_OPTIMALITY.md) for
that sequence and the distinction between numerical evidence and certification.

## 2. Map Lambdas and cylinder volume ratio

Immediately after the motion-law objective, seek design maps for `Lambda_S`,
`Lambda_L` and `V_swept_S / V_swept_L` as functions of the two reservoir
temperatures and required useful output power. Report clearances separately;
do not conflate swept-volume ratio with maximum enclosed-volume ratio.

Distinguish configured geometric Lambda targets from Lambdas inferred from
actual passive-valve events. Define how each is extracted for a general smooth
law. If repeated or missing events make a single Lambda ambiguous, report that
ambiguity and the event chronology rather than forcing a nominal value.

Temperature and power alone need not determine a unique best design. Condition
maps on frequency, working fluid, charge, exchanger hardware, pressure losses,
mechanical limits and the adopted useful-power boundary. Compare feasible
tradeoffs and report sensitivity to uncertain exchanger properties. Do not
present a sparse local search as a universal law or a global optimum.

Expected deliverables are reproducible configuration families, feasible ranges
and best-found regions, with efficiency, size and uncertainty alongside the
Lambdas and volume ratio. The present 25/325-degree, approximately 100 W,
2–10 Hz demonstrator is an initial application, not the full map domain.

## 3. Agreed first whole-machine temperature map (planned)

This is a **whole-machine design study**, not merely kinematic temperature
sensitivity. Use the production variable-property microtube model and its
validity verdict, not the legacy blanket laminar-only screen. The current local
9D reference is still a best-found checkpoint, not a final selected machine;
finish/select the reference explicitly before building the map. No new
optimization or temperature evaluation was launched for this documentation audit.

### Source grid and continuation order

Cold external-air inlet remains **298.15 K** (25 deg C). The reference hot inlet
is **598.15 K** (325 deg C), giving delta T = 300 K.

| Hot inlet (K) | Hot inlet (deg C) | Source delta T (K) |
| ---: | ---: | ---: |
| 508.15 | 235 | 210 |
| 538.15 | 265 | 240 |
| 568.15 | 295 | 270 |
| 598.15 | 325 | 300 |
| 628.15 | 355 | 330 |
| 658.15 | 385 | 360 |

First evaluate the selected delta-T=300 K machine at all six temperatures
**without reoptimization**. Retain that frozen-design sensitivity separately
from the adapted designs. Then optimize sequentially outward in two branches:

- lower-temperature branch: 300 -> 270 -> 240 -> 210 K delta T;
- higher-temperature branch: 300 -> 330 -> 360 K delta T.

At each next temperature, the previous-temperature champion supplies the starting
center. The higher branch begins from the 300 K reference, not from the end of
the lower branch. Temperature is the imposed condition, not an optimizer escape
variable within each campaign.

### Two design passes per temperature

**Pass 1 — thermodynamic/hardware adaptation:** initially keep the previous
motion fixed; optimize the declared thermal/machine design variables, including
exchanger sizing and S/L swept-volume ratio. Inventory is candidate-specific,
not fixed. Keep geometry-derived conductance, hydraulic behavior, hold-up and
wall capacity coupled; do not optimize them independently of the geometry.

**Pass 2 — motion adaptation:** freeze the new hardware/thermal design and
optimize `t1,t2,t3,a_l,b_l,a_s,b_s`. The filling rule still applies to every
candidate, so a change in theta=0 enclosed volume can change mass even in this
motion pass. Reusing a periodic state is only an initial-guess acceleration.

The two-pass design map is a practical **best-found** approximation. It is not a
fully converged simultaneous global optimum of all hardware and motion variables.
Record active variables, bounds, budgets and unchanged numerical tolerances.

### Fixed scale, variable cylinder ratio

Fix total swept volume to the currently selected repository value:

```text
V_swept_S + V_swept_L = 0.0018217821782178217 m^3
                      = 1.8217821782178217 L
r = V_swept_S / V_swept_L
V_swept_L = V_total / (1 + r)
V_swept_S = V_total * r / (1 + r)
```

Source: [thermo3D definition](../outputs/motor_four_stage_thermo3d/definition.json),
`total_swept_volume_m3`. Optimize `r`, rather than holding the large-cylinder
swept volume fixed, so that ratio changes do not change overall machine scale.
The reference ratio is 0.84. Retain each clearance ratio
`V_min / V_swept = 0.01` unless an explicit later campaign changes it; do not
confuse swept volume with maximum enclosed volume.

### Candidate-specific filling rule

For every candidate, use uniform **100000 Pa absolute at theta=0** and the
configured charge temperature (currently 298.15 K). Calculate mass from the
actual connected cylinder and exchanger gas volumes:

```text
m_total = P_fill * V_total_enclosed(theta=0) / (R * T_charge)
```

Gas mass is not held fixed or independently optimized in this whole-machine
study. Geometry and motion at the origin determine it. This design convention
represents the intended inventory after sufficiently long stopped/leaky
atmospheric equilibration; it is not a simulated leakage model or a requirement
that operating periodic pressure equal 100 kPa at theta=0. Historical fixed-mass
experiments remain controlled comparisons with a different boundary condition.

### Power guard and reporting boundary

The present **40 W minimum indicated-power** constraint is a campaign guard,
not a universal physical requirement or useful shaft-power prediction. If it
starts excluding the high-efficiency region as delta T falls, the next
lower-temperature campaign may deliberately lower it. Record the threshold
explicitly for **every temperature campaign**; never silently change it or
compare different power floors as identical feasibility domains. The purpose
of this first map is to understand design evolution.

Indicated efficiency continues to use external-source heat. External-air fan
power is excluded by explicit decision but is not physically zero; external-air
heat-transfer modelling remains a major screening uncertainty. Useful mechanical
power remains unavailable and mechanical losses unknown.

### Mandatory temperature-report fields

Every temperature report must include:

- actual indicated thermal efficiency `eta`;
- Carnot efficiency `eta_C = 1 - T_cold / T_hot`, using absolute source-inlet K;
- fraction of Carnot `eta / eta_C`;
- indicated power (W) and external heat input (W);
- H_i/H_o geometry, including tube counts, lengths, diameters and headers;
- candidate gas inventory (kg), filling pressure/temperature and volume basis;
- S/L swept-volume ratio, total swept volume (m^3), and clearance ratios;
- the four stage durations and `a_l,b_l,a_s,b_s`;
- physical-model identity, numerical convergence, active constraint margins,
  especially the campaign-specific power floor, and the exact output artifact.

The Carnot fraction is essential for judging low-delta-T performance. It is a
source-temperature reference, not evidence of shaft efficiency or validation
of the exchanger approximation. Keep frozen-design results separate from
adapted results and failed/incomplete evaluations separate from feasible designs.

Julia migration is neither implemented nor an approved architecture decision
for this study. Discuss it separately after the temperature experiment.
