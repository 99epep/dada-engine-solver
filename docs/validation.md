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
| Kinematics/mechanisms | Dense-grid volumes/derivatives, branch/closure diagnostics, family parity and scale conventions; `tests/test_research_kinematics_v2.py`, `tests/test_research_v2_thermodynamics.py` |
| Microtube geometry | Derived pitch/sections, circular-frustum hold-up, square/triangular packing, no double-counting; `tests/test_microtube_circular.py`, `tests/test_microtube_geometry.py` |
| Transport/correlations | Temperature domains, frozen oracle, transition boundaries, Mach/pressure-drop/Kn guards; [microtube evidence](EXCHANGER_VALIDATION.md) |
| External thermal boundary | Signed wall/stream balance, refrigeration performance and air/stream parity; `tests/test_external_stream_v3.py` |
| Compiled execution | RHS parity, actual backend/fallback, exact-angle reuse and replay; [numerical execution](SOLVER_PERFORMANCE.md) |
| Campaign/reports | IDs, cache, recovery, torn tails, Sobol/local resume, constraints, historical inspection and plotting; Research/campaign tests |

An invalid or indeterminate model can converge numerically. A physically valid
cycle can fail an explicit design constraint. A best-found feasible point is not
proof of optimality. A recovered trial-state rejection is not a violation of the
final periodic cycle. Keep raw historical diagnostics separately from final-cycle
constraint evidence.

## Software parity and regression tolerances

Software parity establishes equivalence within a measured tolerance, not physical
validation. Check kinematics and instantaneous equations separately from full
periodic integration, with matched inputs and initialization. Tiny representation
differences can alter adaptive integration steps and sampled extrema, so integrated
regressions may require looser measured tolerances than instantaneous comparisons.
Justify tolerances for each quantity and backend; do not relax physical acceptance
constraints merely to recover an old result. The concrete tolerances belong in
the tests, including `tests/test_research_v2_thermodynamics.py`.

An intentional physical-model revision can invalidate historical parity without
being a regression in the revised model. Keep fixtures documenting the old model
unchanged and version current-physics regressions separately.

## Current regression evidence

Current-physics thermal regression cases are versioned separately in
`tests/data/developing_entry_thermal_reference.json` and consumed by
`tests/test_research_v2_thermodynamics.py` and `tests/test_research_sixbar.py`.
These check integrated results and periodic convergence under the declared
model and numerical settings; they do not assert parity with superseded physics.

`tests/test_microtube_developing_entry.py` checks Bennett against independent
HeatLib/Octave reference points, Shah against Eq.192 and its Darcy asymptote,
axial segment additivity and forward/reverse ownership. It also checks continuity
and monotonicity of segment pressure loss across the Reynolds 2300–4000 bridge.
`tests/test_microtube_transition.py` and `tests/test_compiled_transition.py`
cover transition behavior and compiled/reference pointwise parity, including
preservation of domain rejections.

These are equation and software checks. Sources, formulas and domains belong to
the [microtube reference](MICROTUBE_GAS_MODEL.md); experimental comparisons and
their limitations belong to [exchanger validation](EXCHANGER_VALIDATION.md).

## Physical limits of the evidence

The [microtube reference](MICROTUBE_GAS_MODEL.md) records equation sources and
quasi-steady assumptions; [Doty screening](DOTY_SCREENING.md) preserves absolute
experimental discrepancies. Steady-flow correlations, a uniform wall and lumped
collector K do not certify pulsating flow, spatial distribution or manufactured
hardware. A dilute transport temperature domain does not certify phase or EOS
validity. The ideal-generated property table proves compiled reconstruction,
not real-helium thermodynamics/hydraulics.

Current evidence is documented in the references above and in
[numerical execution](SOLVER_PERFORMANCE.md). Earlier application results and
artifacts may remain in Git history or their original sources; they do not define
current validation status.
