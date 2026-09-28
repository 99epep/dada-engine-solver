# Research V3 engineering report

## Implemented boundaries

- `exchangers/external_stream.py`: neutral fluid stream, conservative ten-state
  wall machine, signed external heat and outlet temperature. Historical
  `AirWallMotor` is a compatibility alias; `AirWallExchanger` retains old fields
  and rate keys. `external_microtube.py` reuses production tube geometry,
  material storage, internal gas films and hydraulic links with an explicit
  external wall conductance. No liquid hydraulic correlation was invented.
- `exchangers/wall_cycle.py` and `performance.py`: the same cycle integration
  and periodic convergence support both directions. Performance uses external
  heat with gas-plus-wall storage conservation and all required signed mode
  conditions. Shaft power and mechanical losses remain unavailable.
- `fluids.py`, `state.py`, `dynamics.py`, and charge reconstruction in
  `configuration.py`: conservative rho/u fluid interface, immutable local state,
  and transported energy using state enthalpy. Ideal arithmetic remains
  specialized. Diagnostics now supply volumes for general-fluid reconstruction.
- `tabulated_fluid.py`, `tabulated_backend.py`, `wall_backend.py`, and
  `numerical_primitives.py`: immutable versioned tables, canonical content hash,
  explicit domains/cell masks and compiled interpolation feeding shared
  conservative balances. Python property calls are absent from supported
  compiled table RHS execution. The ideal backend retains its disk cache;
  table JIT caching is currently in-process only.
- `hydraulics.py` and `exchangers/hardware.py`: ideal-gas guards and the
  state-based hydraulic extension protocol. No real-gas choking law is claimed.
- Research schema/basis factories, presets, evaluator and reports: schema 3
  external streams and table definitions, fixed/active stream ownership,
  cooling COP or cooling power objectives, explicit actual backend, fluid-domain
  rejection, neutral stream metrics, current input/cooling powers and existing
  constraint margins. The campaign/persistence engine is reused. The resume
  loader only adds recognition of the new version; old histories are untouched.

## Evidence

All 679 pre-V3 tests remain. New suites and frozen pre-V3 source references
cover instantaneous and trajectory parity, table/domain identity, generic
enthalpy transport, compiled execution, periodic refrigeration, fixed/active
configuration, Sobol/resume and reporting. Complete suite: **728 passed,
0 warnings**; see `README.md` for the command and final run duration.

The exact historical air regression retains its original 2e-11 relative /
1e-11 absolute tolerance and periodic-count checks. The frozen one-cycle
trajectory includes source commit and blob hashes.

Final-source artifacts are under `validated/`:

| Case | Result |
|---|---|
| Declared liquid-scenario refrigerator | Feasible; 11 cycles |
| Cooling / heating / indicated input | 5.36590 / 7.21388 / 1.84798 W |
| Cooling / heating COP | 2.90366 / 3.90365 |
| Absolute cycle conservation residual | approximately 6.1e-13 J, including stored energy |
| Compiled table refrigerator | Feasible; 11 cycles |
| Table relative cooling-power / input / COP errors | 3.50e-8 / 1.11e-7 / 1.46e-7 |
| Neutral-field historical air motor | Feasible; 13 cycles |
| Tiny Sobol run and persisted resume | Two feasible candidates |
| Structured C2 refrigerator setup | Configuration/mechanical preflight only |

The initial high-compression refrigerator startup went below 200 K. The bounded
validation uses explicit 30% clearances and 0.2 Hz instead; the transport limit
was not relaxed. This is not a human-cell optimization or a fair cross-family
optimization comparison.

The local fixed-state benchmark measured 37.68 µs for ideal compiled RHS,
42.40 µs for tabulated compiled RHS and 472.00 µs for Python. The compiled
ratio is **1.125**, with no fallback calls. First calls are separately recorded
(6.023 s and 2.204 s for ideal and table, sharing compiler initialization in one
process). Maximum absolute table/ideal RHS difference was 9.12e-13, maximum
scaled difference 2.60e-13. A separate pre-V3/V3 ideal comparison measured a
1.033 ratio. These timings include the same kinematics and are local measurements,
not universal speed guarantees. Use `examples/benchmark_research_v3_backends.py`
for reproduction; `--ideal-only` also runs against an exported pre-V3 source tree.

The offline HTML was rendered in headless Chromium. The DOM contains both
candidate rows, cooling/input/COP headings, two complete constraint tables and
external fluid/capacity/conductance/heat rows. `validated/report_preview.png`
is the rendered preview. Final saved studies/evaluations and the two-candidate
campaign were rechecked against current runtime identity, including the resume
loader, without extending the campaign.

## Remaining scientific choices and limits

For the human cell, select measured or explicitly uncertain coolant boundary
inputs, motion parameters and physical size, then choose COP or cooling power
under an indicated input bound. Actual human shaft input needs a separate loss
model or measurements. No value for mechanical efficiency was assumed.

For real helium, obtain and validate a single-phase EOS/caloric dataset,
positive energy reference, initialization inversions, interpolation errors and
phase-domain masks. Add compatible flow/choking closures and compiled kernels,
cryogenic transport and heat-transfer evidence. Existing ideal-density duct,
Poiseuille, orifice, sonic and microtube diagnostic equations cannot be labelled
validated real-helium physics. Dilute transport still stops at 200 K.

The table proof uses only analytic ideal-gas data. No helium database,
cryogenic extrapolation, two-phase fluid, pump model or automatic experiment
recommendation was added. Configuration usability and the scientific choices
above remain appropriate topics for the next human review; no approval blocks
use of the completed bounded V3 implementation.
