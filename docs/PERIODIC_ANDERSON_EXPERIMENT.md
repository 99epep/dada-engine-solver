# Periodic Anderson experiment

Completed 2026-09-21 against the post-[Stage-5](SOLVER_ACCELERATION_STAGE5.md)
exact-cache/shared-replay baseline. **The bounded experiment fails the adoption
gate.** All nine cold runs exhaust 100 cycles; the production policy converges in
22. The fastest converged four-bar Anderson run takes 59 cycles, versus 32 for
production. No common setting is selected. The direct solver experiment remains
explicitly opt-in; campaign parsing, campaign identity and defaults are unchanged.

## 1. Exact Anderson formulation

For recent physical map pairs `(x_i, y_i = Phi(x_i))`, project to `z_i, g_i` and
form `f_i = g_i - z_i`. With columns `F`, use the requested constrained least-squares
formulation:

```python
D = F[:, :-1] - F[:, -1, None]
gamma = np.linalg.lstsq(D, -F[:, -1], rcond=None)[0]
alpha = np.r_[gamma, 1 - gamma.sum()]
candidate_z = ((1 - beta) * Z + beta * G) @ alpha
```

Here `memory` counts stored map **pairs** (not differences plus an additional
pair). The implementation is confined to `wall_iteration.py` and an isolated
optional branch of `solve_periodic_wall_motor(..., anderson_acceleration=...)`.
There is no nonlinear optimizer, Jacobian, Newton, Broyden or integrator change.
The ten physical states are periodic unknowns; five quadratures remain excluded
and reset by the ordinary cycle integrator.

## 2. Coordinates and scaling

The original diagnostic SVD construction is now shared source code:
`MassConservingCoordinates` and `tangent_coordinates` in `wall_iteration.py`.
`examples/diagnose_periodic_map.py` imports it; there is no second convention.
For each proposal, the anchor is the latest retained physical endpoint,
`S = diag(abs(anchor))`, and the orthonormal `10 x 9` basis spans the null space
of the mass row multiplied by `S`. Thus `x = anchor + S B z` and
`z = B.T @ ((x - anchor) / scales)`.

Stored map pairs own read-only copies of physical inputs/endpoints and their
normalized residual. Historical states are never renormalized. Mass compatibility
uses `512 * float64_epsilon` relative to the anchor charge (approximately
1.14e-13), allowing accumulated roundoff but not a changed charge. Integrated
endpoints are also checked against the initial charge in the experimental path.
Projection discards only the compatible mass-normal roundoff component when
forming a new initial guess.

## 3. Safeguards and rollback semantics

Controls are strictly validated: memory 2–10 pairs; minimum history between 2 and
memory; finite positive damping <=1; coefficient bound >=1; residual-growth
bound >=1; finite positive relative correction bound. Screen settings are exactly
memory 2/3/4, damping 0.5/0.8/1.0, coefficient L1 <=10, correction <=0.5,
residual growth <=2, minimum history 2.

Rank-deficient or numerically useless least-squares systems discard the oldest
pair and retry down to two pairs. A singular-value ratio below `sqrt(epsilon)`
is treated as unusable. Rank, singular values, effective memory and coefficients
are recorded. Coefficient rejection skips a proposal without failing integration.

The maximum physical relative correction is measured over all ten components.
If needed, the **entire reduced correction** is multiplied by one scalar; no
component is clipped. Preflight requires ten finite positive physical values,
compatible mass and the applied correction bound. Only the next cycle's initial
state changes. The normal Numba/LSODA map, ODE tolerances, exact breakpoints,
exact-angle cache and domain checks are retained.

Only an ordinary complete cycle with the existing maximum normalized periodic
error <=1 certifies convergence. Otherwise a trial is retained only if its error
is <=2 times the previous retained error. Domain/integration failure or excessive
growth rolls back to the last retained physical endpoint, clears Anderson history,
and requires fresh ordinary pairs before another proposal. Ordinary integration
errors still propagate. Interrupted trials return the latest retained complete
endpoint, never a guess or partial trajectory. Every attempted cycle uses budget;
callbacks receive retained completed endpoints only. No proposal is formed after
the last allowed cycle. Anderson is mutually exclusive with adaptive wall
acceleration and suppresses fixed cycle-10 wall extrapolation.

Decisions, accepted/rejected/interrupted trials, reasons, correction bounds,
coefficient norms and timings appear in solver statistics. Completed-cycle history
marks retention explicitly. Final counters distinguish attempts, retained cycles,
proposals, accepted/rejected trials, rollbacks and microtube domain failures.
Default campaign record construction is unchanged; no Anderson TOML section is
accepted or persisted by this stage.

## 4. Focused synthetic tests

The focused suite passes **58 tests**, including the existing periodic-map and
adaptive-wall suites. It covers the known nine-mode linear diagnostic, orthonormal
basis, tangent mass conservation, projection/reconstruction, owned history,
settings validation, synthetic contraction acceleration, duplicate/collinear
residual rank reduction, coefficient rejection, uniform correction scaling,
positivity rejection, incompatible historical charge, excessive and moderate
residual growth, domain/nonfinite rollback, interruption, cycle budget, and
suppression of fixed wall extrapolation.

The deterministic nine-dimensional map has modes 0.95, 0.93, 0.025, 0.01 and
smaller values: Anderson does converge in fewer map calls than plain iteration
there. This establishes algorithm behavior, not superiority to the production
wall-extrapolation policy. A separate scalar test verifies that a 0.95 mode can
require coefficient L1 =39 and is correctly rejected by the experimental cap10.
Test logs are in `outputs/periodic_anderson/`.

## 5. Parameter screen

Reproduction from the source checkout:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 examples/benchmark_periodic_anderson.py \
  --screen --output outputs/periodic_anderson/new_screen
```

The script reuses the frozen Stage-2/5 manifest and the shared frozen candidate
builder, verifies input hashes, and warms the production Numba RHS before timing.
Debian 13.7, Python 3.13.5, NumPy 2.2.4, SciPy 1.15.3; detailed Numba/runtime and
source identities are in `screen/environment.json`. No benchmark or test ran
concurrently with this serial screen. Each entry is **one screening run**, not a
three-repeat timing estimate. Screen times include wrapper construction and the
periodic solve, but exclude final diagnostic replay/reporting, compilation and
artifact serialization. They are not complete candidate-evaluation times.

All policies use the frozen 100-attempt budget. `*` means that budget was exhausted
without convergence; the listed residual is from the last retained physical cycle.
Every completed trial in this screen was retained; attempts equal retained cycles.

| Memory | Damping | Cold cycles / seconds / residual | Four-bar cycles / seconds / residual | Accepted trials cold / four-bar |
|---:|---:|---|---|---|
| 2 | 0.5 | 100* / 18.627 / 8.50484 | 85 / 18.273 / 0.98958 | 8 / 10 |
| 2 | 0.8 | 100* / 18.741 / 6.08140 | 59 / 12.561 / 0.93512 | 4 / 57 |
| 2 | 1.0 | 100* / 18.775 / 5.94881 | 100* / 21.459 / 5.75472 | 2 / 2 |
| 3 | 0.5 | 100* / 18.838 / 6.23187 | 82 / 17.331 / 0.99702 | 3 / 15 |
| 3 | 0.8 | 100* / 18.828 / 7.07812 | 100* / 21.539 / 5.65685 | 3 / 2 |
| 3 | 1.0 | 100* / 18.852 / 5.90661 | 100* / 21.414 / 5.57908 | 2 / 2 |
| 4 | 0.5 | 100* / 19.104 / 6.12220 | 100* / 21.748 / 5.77385 | 2 / 2 |
| 4 | 0.8 | 100* / 19.099 / 5.99291 | 100* / 21.724 / 5.65685 | 2 / 2 |
| 4 | 1.0 | 100* / 18.872 / 5.90661 | 100* / 21.654 / 5.57908 | 2 / 2 |

## 6. Common setting selection

**None selected.** The bounded screen already fails the necessary performance
condition on both target mechanisms. No coefficient bound, tolerance, domain check
or mechanism-specific control was adjusted to manufacture a pass. There is no
second parameter search and no recommendation to use the fastest four-bar setting
as a common policy.

## 7. Comparison with plain fixed point

Plain fixed point also reaches the 100-cycle cap: cold 18.363 s, residual 5.76959;
four-bar 20.991 s, residual 5.52703. Some Anderson settings help four-bar relative
to plain iteration, but none helps cold sufficiently, and plain iteration is not
the adoption target. We do not extend the cycle budget to hide this outcome.

Representative retained-cycle normalized residuals follow. `m2/b0.8` is displayed
because it is the fastest converged four-bar screen entry, **not** a selected
common setting. A dash means the solve had already stopped.

| Case / policy | Cycle 1 | Cycle 2 | Cycle 10 | Cycle 22 | Cycle 32 | Cycle 50 | Cycle 100 |
|---|---:|---:|---:|---:|---:|---:|---:|
| production_cold / plain | 3.442e+05 | 4687 | 3071 | 1315 | 652.6 | 186 | 5.77 |
| production_cold / production | 3.442e+05 | 4687 | 3071 | 0.4322 | — | — | — |
| production_cold / anderson_m2_b0.8 | 3.442e+05 | 4687 | 3242 | 1387 | 688.2 | 196.1 | 6.081 |
| smooth_four_bar / plain | 3.455e+05 | 5064 | 2811 | 1186 | 585.2 | 166.6 | 5.527 |
| smooth_four_bar / production | 3.455e+05 | 5064 | 2811 | 2.46 | 0.2128 | — | — |
| smooth_four_bar / anderson_m2_b0.8 | 3.455e+05 | 5064 | 44.11 | 20.66 | 7.574 | 2.801 | — |

The complete cycle-by-cycle data for **all** policies are in
[`residuals.csv`](../outputs/periodic_anderson/residuals.csv), with detailed
per-case decisions and all 15 trajectory rows under
[`screen/`](../outputs/periodic_anderson/screen/).

## 8. Comparison with current fixed wall extrapolation

The production policy remains substantially better. The table uses the same
`m2/b0.8` comparison only as an illustration; it is not an adopted setting.

| Case | Production cycles | Screen seconds | Anderson cycles | Screen seconds | Rollbacks |
|---|---:|---:|---:|---:|---:|
| Warm | 1 (Stage 5) | not screened | not run | — | — |
| Nearby | 3 (Stage 5) | not screened | not run | — | — |
| Cold | 22 | 3.883 | 100, not converged | 18.741 | 0 |
| Four-bar | 32 | 6.887 | 59 | 12.561 | 0 |
| Six-bar | 27 (Stage 5) | not screened | not run | — | — |

Historical Stage-5 cycle counts are explicitly labelled; their complete-evaluation
timings must not be mixed with this screen's periodic-only measurements.

## 9. Cold/four-bar/six-bar cost result

Cold fails convergence for every Anderson setting. The best four-bar screen entry
is about 1.82 times the production periodic-screen time, rather than at least 25%
faster. Six-bar validation and three serial timing repeats were not launched:
they are downstream of common-setting selection, which failed. This follows the
requested early stop condition rather than treating missing validation as a pass.

Total proposal-decision time is 0.044–0.079 s per screen solve, about
0.5–0.8 ms per decision and below 0.5% of solve time. The problem is the physical
cycle count and rejected proposals, not expensive linear algebra.

## 10. Warm-start result

No selected policy exists, so the expanded warm/nearby campaign-equivalence gate
was not run. The solver checks ordinary convergence before requesting a proposal,
and requires two pairs; a first-cycle convergence therefore cannot trigger
Anderson. This code property is not presented as a measured warm performance gate.
Current defaults and warm-start/candidate identity rules remain unchanged.

## 11. Tighter-periodic accuracy audit

| Case | Stage-5 exact reference cycles | Anderson cycles/time | Audit outcome |
|---|---:|---|---|
| Production cold | 33 | not run | stopped at failed normal screen |
| Smooth four-bar | 52 | not run | stopped at failed normal screen |

No tighter-periodic accuracy claim is made. The tight-manifest runner is available
for reproducibility/future explicitly authorized work, but was not invoked after
failure of the required normal screen.

## 12. Rejected proposals, trials and domain failures

Across all 18 Anderson runs, decision counts are:

- `proposed`: 122.
- `skipped_coefficients`: 1568.
- `skipped_insufficient_history`: 18.

There were **zero rejected integrated trials, zero rollbacks and zero microtube
domain failures** in the physical screen. Skipped coefficient proposals do not
consume a physical trial: the next attempt uses the retained endpoint normally.
Recovery mechanisms are exercised by deterministic regression tests instead.

The restrictive coefficient cap explains the dominant failure mode. For a scalar
contraction `f_next = lambda * f_previous`, a two-pair residual cancellation needs
`alpha = [-lambda/(1-lambda), 1/(1-lambda)]`. At lambda=0.95, this is `[-19,20]`,
with L1=39; at 0.93, L1 is approximately 27.57. These both exceed 10. This algebra
explains why the safe bounded configuration can stagnate near ordinary relaxation;
it does not establish that changing the cap would be physically safe or pass the
scientific-output gate. Larger caps and a different history policy were not tried.

## 13. Final scientific-output comparison

The three converged four-bar Anderson screen runs already fail the original
Stage-3/4/5 aggregate accuracy scales: 1e-6 relative plus 1e-12 absolute, with
scaled final-state acceptance <=1. Work and power have the same relative
difference at fixed speed. Efficiency is indicated gas work divided by heat input.

| Four-bar setting | Maximum scaled state distance | Relative work/power difference | Relative heat-in difference | Relative efficiency difference |
|---|---:|---:|---:|---:|
| anderson_m2_b0.5 | 13.919 | -8.257e-05 | +4.589e-04 | -5.413e-04 |
| anderson_m2_b0.8 | 17.744 | +1.311e-04 | -5.857e-04 | +7.172e-04 |
| anderson_m3_b0.5 | 17.150 | -1.017e-05 | +5.109e-04 | -5.208e-04 |

Their per-cycle absolute relative mass residuals remain <=7.3e-16. Their ordinary
periodic residuals pass <=1, but passing that stopping criterion does **not** prove
agreement with the reference scientific outputs at the separate comparison scale.
Cold has no converged Anderson endpoint to compare.

Full extrema, valve events, Reynolds/Mach, microtube ranges/fractions,
validity/feasibility, constraints and full-report materialization were deliberately
not evaluated for this failed screen. They are not certified equivalent. Raw
physical trajectories, final ten-state vectors, basic heat/work/conservation and
normalized residuals are retained for inspection. No failed criterion was hidden,
no tolerance relaxed, and no converged scientific result is claimed for capped
runs.

## 14. Adoption decision

**Do not adopt this bounded Anderson policy.** Keep the current fixed-wall policy
and Stage-5 exact improvements. Retain the direct opt-in experimental API, tests,
script and negative evidence, but add no campaign parsing or identity controls.
No persistent campaign was launched. This result does not justify claiming that
Anderson removes the need for further convergence research; neither does it
authorize starting Newton, Broyden or Julia in this stage. The requested bounded
experiment stops here.
