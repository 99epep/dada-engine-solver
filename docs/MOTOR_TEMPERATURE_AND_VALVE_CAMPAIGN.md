# Motor temperature and passive-valve topology campaign

## Purpose

This note records the main lessons from the temperature-map and passive-valve-placement campaigns performed on the DADA engine solver in September 2026.

The study had four successive objectives:

1. map the optimized motor over source temperature difference;
2. locate the region that maximizes indicated efficiency relative to Carnot efficiency;
3. determine whether passive check-valve placement materially affects the cycle;
4. examine how the preferred machine proportions and kinematics evolve with source temperature difference.

The results are exploratory numerical design evidence. They are not a proof of global optimality and they do not yet include mechanical losses, finite valve dynamics, structural constraints, or a realizable smooth mechanism.

Unless stated otherwise, the cold source is fixed at **298.15 K**, the working gas is air, the charge condition is **100 kPa absolute at 298.15 K**, and the total swept volume is fixed at **1.821782 L** with 1% cylinder clearance ratios.

---

## 1. Solver acceleration changed the practical scope of the study

A major enabling result of this campaign was computational rather than thermodynamic.

The introduction of the compiled Numba wall-RHS backend, exact kinematics caching, shared diagnostic replay, and related solver-acceleration work reduced typical candidate evaluation times from roughly order-100-second evaluations to only a few seconds on the same machine.

Representative short-campaign medians were approximately:

- **4.3 s/candidate** for thermo candidates;
- **3.5 s/candidate** for motion candidates.

This is roughly a **25–30× practical speed-up** for this class of campaign.

That improvement changed the optimization strategy itself: instead of avoiding additional sweeps because of cost, it became practical to compare several valve topologies, refine local temperature points, and run independent two-pass searches for each architecture.

---

## 2. Reference DD temperature map

The historical/reference valve arrangement is denoted **DD**, meaning that the passive check valve of each exchanger is placed on its downstream side.

The optimized DD map showed increasing absolute indicated efficiency and power with source temperature difference, while efficiency relative to Carnot reached a broad maximum in the vicinity of 250–260 K.

| ΔT (K) | Indicated efficiency (%) | η / η_Carnot (%) | Power (W) |
|---:|---:|---:|---:|
| 160 | 13.8599 | 39.6869 | 11.107 |
| 200 | 17.5165 | 43.6293 | 21.968 |
| 240 | 19.8259 | 44.4554 | 31.719 |
| 250 | 20.3276 | 44.5702 | 35.214 |
| 256 | 20.6023 | 44.5967 | 37.384 |
| 260 | 20.7994 | 44.6508 | 39.120 |
| 270 | 21.1008 | 44.4016 | 40.841 |
| 300 | 22.0834 | 44.0306 | 52.407 |
| 330 | 22.8184 | 43.4346 | 63.447 |
| 360 | 23.7000 | 43.3282 | 70.352 |

### 2.1 Location of the relative-efficiency maximum

Initial quadratic fits to the coarse map placed the maximum of `η / η_Carnot` in the approximate 247–260 K range.

Additional DD points at 250, 256, and 260 K progressively narrowed the useful region. A full 96+96 refinement at **ΔT = 260 K** produced:

- indicated efficiency: **20.7994%**;
- indicated power: **39.120 W**;
- `η / η_Carnot`: **44.6508%**.

The optimum is broad and shallow, so assigning physical importance to the exact integer-kelvin peak would be misleading. For subsequent work, **260 K** was retained as a convenient reference temperature difference.

This choice is intentionally pragmatic: 260 K is easy to remember, lies essentially at the top of the DD relative-efficiency plateau, and is sufficiently close to the local numerical maximum for comparative architecture work.

---

## 3. Passive-valve topology

### 3.1 Definitions

The nominal circulation is:

`S -> H_i -> L -> H_o -> S`

where:

- `S` = small cylinder;
- `L` = large cylinder;
- `H_i` = heat-input exchanger;
- `H_o` = heat-output exchanger.

For each exchanger the passive check valve may be placed either upstream (`U`) or downstream (`D`) of the exchanger.

The four topologies are therefore:

| Code | H_i valve | H_o valve | Geometric interpretation |
|---|---|---|---|
| DD | downstream | downstream | historical arrangement |
| UD | upstream | downstream | both valves on the **small-cylinder side** |
| DU | downstream | upstream | both valves on the **large-cylinder side** |
| UU | upstream | upstream | one valve on each cylinder side, both before their exchanger in nominal flow |

The terms upstream/downstream refer to the nominal circulation direction, not to a permanent pressure ordering.

---

## 4. Full topology re-optimization at ΔT = 260 K

A frozen replay is useful as a screening tool, but valve placement changes instantaneous pressure differences, mass flow, exchanger Reynolds numbers, and thermal transfer. Therefore the fair comparison is an **independent re-optimization** of each topology.

At 260 K, each alternative topology was optimized with:

- 96 thermo evaluations;
- 96 motion evaluations;
- the same physical bounds and constraints;
- the same total swept volume and charge condition.

The result was unambiguous:

| Topology | Efficiency (%) | η / η_Carnot (%) | Power (W) | Δη vs DD (percentage points) | ΔP vs DD (W) |
|---|---:|---:|---:|---:|---:|
| DD | 20.7994 | 44.6508 | 39.120 | +0.0000 | +0.000 |
| UD | 20.9383 | 44.9489 | 41.438 | +0.1389 | +2.318 |
| DU | 21.7638 | 46.7210 | 42.184 | +0.9644 | +3.065 |
| UU | 22.1326 | 47.5128 | 44.383 | +1.3332 | +5.263 |

### 4.1 Main result

**UU is clearly superior at the 260 K reference point.**

Relative to DD, the fully re-optimized UU solution gains:

- **+1.333 percentage points** of absolute indicated efficiency;
- **+2.862 percentage points** in `η / η_Carnot`;
- **+5.263 W** of indicated power.

DU is also substantially better than DD, while UD provides only a small gain after full optimization.

This establishes that passive-valve placement is not a secondary hydraulic detail. It is a first-order design variable because it changes the pressure history seen by each exchanger and cylinder, which in turn changes mass flow, Reynolds number, exchanger behavior, and the complete thermodynamic cycle.

---

## 5. Does the preferred topology change with temperature?

To answer this without repeating full 96+96 searches at every temperature, a deliberately short topology map was run at:

**160, 200, 240, 270, 300, 330, and 360 K**

For every temperature:

- the existing DD champion was retained as the reference;
- UD, DU, and UU were independently re-optimized;
- each alternative topology received **24 thermo + 24 motion evaluations**.

This screening is not intended to produce final hardware optima. Its purpose is to determine the ordering of the valve architectures.

### 5.1 Result

**UU won at every screened temperature**, both in absolute indicated efficiency and in efficiency relative to Carnot.

The consistency of this result over the full 160–360 K interval is stronger evidence than the 260 K result alone. Within the explored domain, no topology transition was observed.

The practical design conclusion is therefore:

> **UU should be used as the default passive-valve architecture for subsequent thermodynamic optimization, unless later studies introduce new physics that materially alter the ranking.**

---

## 6. Evolution of the optimized UU machine with source temperature difference

The total swept volume was fixed during the campaign. Consequently, this study does **not** determine whether the complete engine should become globally larger or smaller with ΔT.

It does determine how the fixed total swept volume should be divided between the two cylinders, and how the exchanger geometry shifts.

### 6.1 Cylinder swept-volume ratio

| ΔT (K) | Small swept volume (cm³) | Large swept volume (cm³) | S/L swept ratio |
|---:|---:|---:|---:|
| 160 | 891.6 | 930.2 | 0.9584 |
| 200 | 865.5 | 956.3 | 0.9050 |
| 240 | 850.8 | 970.9 | 0.8763 |
| 260 | 843.1 | 978.7 | 0.8614 |
| 270 | 838.5 | 983.3 | 0.8528 |
| 300 | 829.1 | 992.7 | 0.8351 |
| 330 | 825.5 | 996.3 | 0.8285 |
| 360 | 823.1 | 998.7 | 0.8242 |

The trend is exceptionally clear:

- increasing ΔT makes the optimum machine progressively more asymmetric;
- the **large-cylinder share increases**;
- the **small-cylinder share decreases**.

Across the screened range, the optimized S/L swept-volume ratio falls from approximately **0.958 at 160 K** to **0.824 at 360 K**.

Because total swept volume was constrained, this trend should be interpreted as a preferred **volume ratio**, not yet as an absolute machine-size law.

---

## 7. Exchanger-size trends

The microtube inner diameter remained fixed at 0.33 mm. A useful measure of exchanger size is therefore internal tube surface area:

`A_tube = π d N L`

where `N` is tube count and `L` is tube length.

| ΔT (K) | H_i tubes × length | H_i internal area (m²) | H_o tubes × length | H_o internal area (m²) |
|---:|---:|---:|---:|---:|
| 160 | 4746 × 47.6 mm | 0.234 | 4675 × 24.4 mm | 0.118 |
| 200 | 4290 × 45.4 mm | 0.202 | 4623 × 23.5 mm | 0.112 |
| 240 | 4092 × 40.2 mm | 0.170 | 4218 × 24.8 mm | 0.109 |
| 260 | 3834 × 40.3 mm | 0.160 | 4160 × 26.3 mm | 0.113 |
| 270 | 3263 × 42.5 mm | 0.144 | 4028 × 25.6 mm | 0.107 |
| 300 | 3238 × 46.4 mm | 0.156 | 3781 × 25.2 mm | 0.099 |
| 330 | 3223 × 46.9 mm | 0.157 | 3782 × 24.7 mm | 0.097 |
| 360 | 3021 × 40.6 mm | 0.127 | 3595 × 29.3 mm | 0.109 |

### 7.1 Heat-input exchanger H_i

The strongest geometric trend is the reduction in required hot-side exchanger area as ΔT increases.

The optimized H_i internal surface falls from roughly:

- **0.234 m² at 160 K**
- to **0.160 m² at 260 K**
- to **0.127 m² at 360 K**.

This is physically reasonable: as the hot-source temperature rises, the thermal driving force available to H_i increases, so less exchanger area is required for the same qualitative role in the cycle.

### 7.2 Heat-output exchanger H_o

H_o changes much less strongly and remains of order **0.10–0.12 m²** over most of the map.

This asymmetry is also physically reasonable. The cold source remains fixed at 298.15 K, while rejected heat increases as the engine operates at larger temperature difference and higher power.

The hot exchanger therefore shrinks much more strongly with increasing ΔT than the cold exchanger.

### 7.3 Interpretation of tube count and length

Individual tube count and tube length show some non-monotonicity, particularly in the short 24+24 screening points.

Those values should not yet be over-interpreted independently. Several nearby `(N, L)` combinations can provide similar area, pressure drop, hold-up volume, and transfer performance.

For trend analysis, total tube area and exchanger volume are more robust indicators than tube count or tube length alone.

---

## 8. Evolution of the UU kinematics

The implemented four-stage kinematic law uses transition times `t1, t2, t3` and intermediate normalized piston levels `a_l, b_l, a_s, b_s`.

The four phase labels are:

1. low-pressure exchange;
2. compression;
3. high-pressure exchange;
4. expansion.

### 8.1 Phase-duration trends

| ΔT (K) | LP exchange (deg) | Compression (deg) | HP exchange (deg) | Expansion (deg) |
|---:|---:|---:|---:|---:|
| 160 | 180.6 | 36.3 | 117.4 | 25.7 |
| 200 | 178.0 | 46.4 | 101.2 | 34.4 |
| 240 | 179.5 | 38.8 | 105.6 | 36.1 |
| 260 | 177.9 | 54.8 | 92.8 | 34.4 |
| 270 | 173.0 | 50.8 | 100.6 | 35.6 |
| 300 | 174.7 | 47.9 | 104.7 | 32.7 |
| 330 | 169.1 | 50.2 | 104.6 | 36.0 |
| 360 | 169.1 | 52.6 | 100.6 | 37.8 |

The broad trends are:

- **low-pressure exchange becomes shorter** as ΔT increases;
- **compression becomes longer**;
- **expansion tends to become longer**;
- the high-pressure exchange phase tends to shorten, although the short screening produces visible local noise.

The 260 K point is more heavily optimized than the other UU temperature points and therefore deserves more confidence than small point-to-point fluctuations in the screening map.

### 8.2 Intermediate piston levels

The knot values are:

| Knot | Small-cylinder normalized swept level | Large-cylinder normalized swept level |
|---|---:|---:|
| T0 | b_s | 1 |
| T1 | 1 | b_l |
| T2 | a_s | 0 |
| T3 | 0 | a_l |
| T4 | b_s | 1 |

The optimized levels are:

| ΔT (K) | a_l | b_l | a_s | b_s |
|---:|---:|---:|---:|---:|
| 160 | 0.4467 | 0.0964 | 0.4213 | 0.1348 |
| 200 | 0.4098 | 0.1052 | 0.3916 | 0.1505 |
| 240 | 0.3789 | 0.1121 | 0.3703 | 0.1469 |
| 260 | 0.3723 | 0.0996 | 0.3609 | 0.1135 |
| 270 | 0.3666 | 0.1115 | 0.3590 | 0.1456 |
| 300 | 0.3590 | 0.1141 | 0.3504 | 0.1584 |
| 330 | 0.3618 | 0.1191 | 0.3486 | 0.1584 |
| 360 | 0.3482 | 0.1241 | 0.3295 | 0.1570 |

The clearest trend is the parallel reduction of `a_l` and `a_s` with increasing ΔT.

From 160 to 360 K:

- `a_l` falls from about **0.447 to 0.348**;
- `a_s` falls from about **0.421 to 0.330**.

Thus, near the end of high-pressure exchange / beginning of expansion, both cylinders favor progressively lower intermediate normalized volumes as source temperature difference rises.

The trends in `b_l` and especially `b_s` are weaker and noisier.

The fact that several kinematic parameters evolve smoothly with temperature strongly suggests that a lower-dimensional thermodynamic relation may exist between the idealized motion law and the source/gas temperature ratios. Deriving that relation is a promising future analytical task.

---

## 9. Numerical confirmation of the quasi-isobaric exchange-volume relation

Section 6.1.2 of the thermodynamic study derives the approximate quasi-isobaric transfer relation

`dV_r / (-dV_d) ≈ T_r / T_d`

for gas transferred between donor and receiver volumes at approximately equal pressure.

A naive comparison of total cylinder swept volumes with hot- and cold-source temperatures is **not** the correct test of this relation.

The four-stage law uses only part of each cylinder stroke during each exchange phase, so the relevant exchanged volume increments include the intermediate levels `a` and `b`.

### 9.1 Low-pressure exchange

From T0 to T1:

- the small cylinder increases from `b_s` to 1;
- the large cylinder decreases from 1 to `b_l`.

Therefore:

`((1 - b_s) V_S) / ((1 - b_l) V_L) ≈ T_S / T_L`

or equivalently:

`T_L / T_S ≈ (1 - b_l) / ((1 - b_s) (V_S/V_L))`

### 9.2 High-pressure exchange

From T2 to T3:

- the small cylinder decreases from `a_s` to 0;
- the large cylinder increases from 0 to `a_l`.

Therefore:

`(a_l V_L) / (a_s V_S) ≈ T_L / T_S`

or:

`T_L / T_S ≈ a_l / (a_s (V_S/V_L))`

### 9.3 Comparison of the two independent temperature-ratio estimates

| ΔT (K) | T_L/T_S inferred from LP exchange | T_L/T_S inferred from HP exchange |
|---:|---:|---:|
| 160 | 1.090 | 1.106 |
| 200 | 1.164 | 1.156 |
| 240 | 1.188 | 1.168 |
| 260 | 1.179 | 1.198 |
| 270 | 1.219 | 1.197 |
| 300 | 1.260 | 1.226 |
| 330 | 1.263 | 1.252 |
| 360 | 1.261 | 1.282 |

The two independent estimates agree to within only a few percent across the full temperature range.

At the fully optimized **260 K / UU** point:

- LP-exchange inference: **T_L/T_S ≈ 1.179**
- HP-exchange inference: **T_L/T_S ≈ 1.198**

The difference is only about **1.6%**.

This is a particularly important result of the campaign:

> The numerically optimized cylinder ratio and intermediate piston levels independently recover the quasi-isobaric exchange-volume relation derived analytically in the thermodynamic study.

This provides a meaningful internal validation of the thermodynamic interpretation of the four-stage cycle.

It also explains why the optimized total swept-volume ratio must **not** be compared directly with `T_cold / T_hot`. The relevant temperatures are the instantaneous or phase-mean gas temperatures in the donor and receiver volumes, and the exchanged volume increments are generally only fractions of the total cylinder strokes.

---

## 10. What the campaign established

The strongest conclusions are:

1. **The solver is now fast enough for architecture-level optimization.**  
   The acceleration work made multi-temperature and multi-topology campaigns practical on modest hardware.

2. **The DD relative-Carnot optimum is broad and lies near 260 K.**  
   A convenient reference point of ΔT = 260 K was selected.

3. **Valve placement is a first-order thermodynamic design variable.**  
   Moving the check valves changes the entire pressure/flow/heat-transfer history.

4. **UU is the preferred passive-valve topology throughout the explored 160–360 K range.**  
   No topology transition was observed.

5. **At 260 K, full UU re-optimization produces a large gain over DD.**  
   The improvement is +1.333 percentage points of indicated efficiency and +5.263 W.

6. **The preferred cylinder swept-volume ratio decreases smoothly as ΔT rises.**  
   The large cylinder becomes progressively larger relative to the small cylinder.

7. **The hot-side exchanger becomes substantially smaller as ΔT rises.**  
   The cold-side exchanger varies much less strongly.

8. **The preferred kinematic law changes systematically with temperature.**  
   LP exchange shortens, compression lengthens, and the intermediate `a_l` and `a_s` levels decrease.

9. **The optimized kinematics and cylinder ratio numerically reproduce the quasi-isobaric volume/temperature relation.**  
   This is one of the most satisfying theoretical cross-checks produced by the campaign.

---

## 11. Important limitations

These results must be interpreted within the current model scope.

### 11.1 Finite optimization budgets

The 260 K topology comparison used full 96+96 searches.

The topology-vs-temperature map used only 24+24 evaluations per alternative topology. It is suitable for topology ranking and trend identification, but not for claiming fully converged hardware dimensions at every temperature.

### 11.2 Fixed total swept volume

The campaign optimized the **ratio** between cylinder swept volumes while holding their sum fixed.

Therefore it does not yet answer:

- the optimum absolute total displacement;
- optimum machine scale for a target power;
- power density versus efficiency trade-offs.

### 11.3 Fixed operating speed

The current map does not establish how optimum scale changes with rotational frequency.

### 11.4 Idealized valve model

The study concerns ideal passive one-way topology. Valve inertia, lift dynamics, stiffness, leakage, impact, fatigue, and detailed mechanical implementation are intentionally outside the present preliminary solver scope.

### 11.5 Idealized four-stage motion law

The current piecewise-linear law is a thermodynamic design probe, not a mechanically realizable target motion.

The longer-term task remains to determine whether a smooth and mechanically simple motion law can reproduce the main thermodynamic features found here.

---

## 12. Recommended next directions

The campaign suggests the following priorities:

1. adopt **UU** as the default valve topology for subsequent motor studies;
2. repeat selected temperature points with deeper UU optimization only when precise dimensional trends are required;
3. investigate **total swept volume / machine scale** as a free variable;
4. investigate the relation between source temperature, gas temperature ratios, cylinder ratio, and the four-stage kinematic parameters;
5. search for a lower-dimensional or analytic motion law that reproduces the optimized temperature-dependent kinematics;
6. eventually return from the idealized four-stage law to smooth realizable mechanisms;
7. preserve the distinction between thermodynamic topology optimization and later mechanical implementation details.

---

## 13. Data provenance

The principal committed solver artifacts used by this note are:

- `outputs/motor_temperature_map/`
- `outputs/motor_temperature_peak_refinement/dT_260K/`
- `outputs/motor_valve_optimization_260K/`
- `outputs/motor_valve_topology_temperature_screen/`
- `docs/FOUR_STAGE_OPTIMIZATION.md`
- `src/dada_solver/four_stage_kinematics.py`

The valve-topology notation and all numerical values in this document correspond to the solver state used for the September 2026 campaign.

