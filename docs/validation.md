# Validation and evidence

Tests establish specific numerical or software properties. They do not validate
an experimental DADA machine, a global optimum or all physical assumptions of a
candidate. Always separate convergence, conservation, model applicability and
study feasibility.

## Run checks from the checkout

```sh
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 -m pytest -q
```

Use targeted files during development; run the complete suite after consequential
changes. Optional Numba tests require the `numba` extra. Normal transport tests
use frozen CoolProp oracle data and do not require CoolProp at runtime.

## What to verify

| Layer | Evidence / tests |
|---|---|
| Conservative equations | Ideal-state reconstruction, closed adiabatic `P V^gamma`, donor outflow relation, internal mass/enthalpy cancellation, integrated mass/energy residuals |
| Integration and valves | Event-state continuity, passive directionality, hysteresis where supported, periodic endpoints, refinement and interruption |
| Kinematics/mechanisms | Dense-grid volumes/derivatives, branch/closure diagnostics, family parity and scale conventions; [Research evidence](DADA_ENGINE_RESEARCH_VALIDATION.md) |
| Microtube geometry | Derived pitch/sections, circular-frustum hold-up, legacy compatibility, no double-counting; `tests/test_microtube_circular.py`, `tests/test_microtube_geometry.py` |
| Transport/correlations | Temperature domains, frozen oracle, transition boundaries, Mach/pressure-drop/Kn guards; [microtube evidence](EXCHANGER_VALIDATION.md) |
| External thermal boundary | Signed wall/stream balance, refrigeration performance and legacy air parity; `tests/test_external_stream_v3.py` |
| Compiled execution | RHS parity, actual backend/fallback, exact-angle reuse and replay; [numerical execution](SOLVER_PERFORMANCE.md) |
| Campaign/reports | IDs, cache, recovery, torn tails, Sobol/local resume, constraints, historical inspection and plotting; Research/campaign tests |

An invalid or indeterminate model can converge numerically. A physically valid
cycle can fail an explicit design constraint. A best-found feasible point is not
proof of optimality. A recovered trial-state rejection is not a violation of the
final periodic cycle. Keep raw historical diagnostics separately from final-cycle
constraint evidence.

## Current check record

The complete suite after Bennett/Shah entry and continuous-transition integration
on 2026-10-04 reports **1190 passed, 2 failed, no warnings** (503.98 s).
The focused exchanger/Research/backend check reports **218 passed**.
The two failing test names already failed before this physical revision:

- `test_external_stream_v3.py::test_pre_v3_source_trajectory_and_neutral_roundtrip`:
  the pre-V3 frozen RHS uses historical physics; current maximum absolute
  discrepancy is 1.09678032.
- `test_hybrid_compact_kinematics.py::test_historical_champion_thermodynamic_parity`:
  current evaluation converges in 12 cycles versus seven in the old fixture
  (the pre-entry revision already differed, at three cycles).

Do not describe this as a completely green suite. Historical fixtures remain
unchanged. Four explicit current-physics thermal regression cases are versioned
separately in `tests/data/developing_entry_thermal_reference.json`; tolerances
are unchanged. Bennett has 50 independent HeatLib/Octave reference points;
Shah has Eq.192, Darcy-asymptote, axial-additivity and reverse-flow checks.
Continuity/monotonicity and compiled pointwise tests cover the revised
2300–4000 bridge. Bounded DD13 replay evidence is under
`outputs/research_microtube_developing_entry/continuous_transition_v1/`.

## Physical limits of the evidence

The [microtube reference](MICROTUBE_GAS_MODEL.md) records equation sources and
quasi-steady assumptions; [Doty screening](DOTY_SCREENING.md) preserves absolute
experimental discrepancies. Steady-flow correlations, a uniform wall and lumped
collector K do not certify pulsating flow, spatial distribution or manufactured
hardware. A dilute transport temperature domain does not certify phase or EOS
validity. The ideal-generated property table proves compiled reconstruction,
not real-helium thermodynamics/hydraulics.

Historical timings, refinements and abandoned numerical experiments are condensed
in [numerical execution](SOLVER_PERFORMANCE.md). Historical application results
remain in [the study index](history/README.md), with their original constraints
and source artifacts. Old test counts and intermediate profiler dumps are not
maintained as current status in the user guides.
