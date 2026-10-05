# Working fluids: thermodynamics, transport and compiled tables

## Thermodynamics and transport are separate

The production fluid remains `CaloricallyPerfectGas`: P = ρRT, u = CvT,
h = CpT, constant Cp/Cv/gamma and Z = 1. Its original convenience methods and
serialized fields remain unchanged. It now also implements `ThermodynamicFluid`:

- `state_from_rho_u(density, specific_energy)` reconstructs a local state;
- `state_from_rho_t(density, temperature)` supports initial filling;
- `density_from_pt(pressure, temperature)` supports pressure-defined filling.

The conservative solver still stores mass and internal energy, with volume from
kinematics: ρ = m/V, u = U/m. Immutable `FluidState` contains density, u, T, P,
h and optional Z, Cp, Cv and sound speed. Gas transport energy uses the upstream
state's enthalpy. The generic dynamics no longer assumes h(T). Ideal dynamics
retain the original arithmetic; general fluids reconstruct the four local
states once per instantaneous point.

This permits three separate physical levels: calorically perfect ideal gas,
thermally perfect ideal gas with variable caloric properties, and single-phase
real gas with state-dependent EOS/caloric properties. Only the first is a
validated production working fluid here. The table is an execution prototype,
not a new validated physical helium model. The current conservative-state and
table prototype require positive internal energies; a future dataset must use
an appropriate consistent energy reference or explicitly extend that validation.

`DiluteGasTransport` remains distinct. Its temperature-dependent viscosity,
conductivity and correlation Cp do not replace the thermodynamic EOS/caloric
model. Its species-dependent property-temperature domains are air 100–1000 K, He
50–1000 K, and N2/Ar 200–1000 K; see [MICROTUBE_GAS_MODEL.md](MICROTUBE_GAS_MODEL.md).
These are not phase/EOS validity domains. No real-fluid helium EOS, two-phase
model or cryogenic hydraulic validation is introduced.

## Runtime table and scientific identity

`tabulated_fluid.TabulatedFluid` schema 1 stores strictly increasing positive
rho/u axes, a property array in published field order, an explicit boolean
single-phase **cell** mask, T/P validity bounds and provenance. Fields are T, P,
h, Z, Cp, Cv and sound speed, in SI units. Arrays have immutable byte-backed
storage. `to_data` / `from_data` serialize a deterministic JSON representation;
`content_hash` is SHA-256 of its canonical content, including provenance.
`identity` records model/version, hash, interpolation version and domains.
Research schema 3 embeds this artifact in the portable machine basis. Its byte
hash and the canonical fluid hash both participate in study/runtime identity.
A changed table cannot reuse the old candidate ID or exact result cache.

`bilinear_rho_u_v1` performs direct table reconstruction with no iterative
inversion in the RHS and no extrapolation. Axis endpoints are closed. Every
lookup also checks the cell mask and T/P limits. Unavailable states raise
`FluidDomainError`; Research records `invalid_fluid_domain`. At a shared cell
edge the right-hand cell is used; rejecting an invalid adjacent cell is
conservative. The mask is a dataset provider's declared single-phase validity,
not an implemented phase detector. Table stability, thermodynamic consistency
and interpolation error must be established for each future scientific dataset.

`ideal_validation_table` generates data deterministically from the existing
analytic gas. T, P and h are linear/bilinear in these coordinates, making this a
strict reconstruction and transported-energy test. Sound-speed interpolation
has ordinary interpolation error and is not used to replace the exact ideal
hydraulic law in this validation. An optional ideal-reference declaration is
verified against **all** table property nodes before ideal hydraulic use is
allowed; changing P/T/h while retaining that claim is rejected.

A general table can reconstruct rho/u states already. General rho/T or P/T
initialization inversions are deliberately not supplied; only the verified
ideal-reference table provides those operations. A future real-fluid provider
must implement and validate them. Source tables in rho/T coordinates can later
be transformed offline to a rho/u runtime table with consistent energy reference
and a certified domain. V3 does not guess a helium dataset or an inversion.

## Compiled execution

`WallRHS` dispatches explicitly:

| Fluid and capability | Actual backend |
|---|---|
| Calorically perfect gas, supported internal films/links | `numba` (specialized ideal reconstruction) |
| Verified ideal-generated table, supported films/links | `numba_tabulated` |
| Python selected, Numba unavailable or unsupported compiled family | `python`, with reason where applicable |

The ideal kernel reconstructs T/P analytically. The table kernel performs lookup
and interpolation entirely inside Numba, including the enthalpy used for mass
transfer. Both feed the same conservative balance kernel and existing internal
film/flow primitives. No Python property-library call occurs inside a supported
compiled RHS evaluation. Tests replace the Python property function with a
failure to verify this boundary. Unsupported numerical states fall back to the
reference path for the same physical error/domain checks; this never bypasses
table validity. Fallback counts and reasons are recorded.

The ideal disk-cache path retains its content/runtime namespace. The table
prototype intentionally uses in-process JIT only; it does not enable a disk cache
with incomplete transitive dependency invalidation. Runtime records include the
table identity and table-backend source hashes, actual backend, first call,
calls and fallbacks. Research additionally hashes the full solver source tree.

## Hydraulic boundary: why this is not validated real-gas helium

The following existing closures assume ideal density and/or constant gamma:

- `CompressibleOrifice`: ideal isentropic pressure ratio, sonic/choked mass flow;
- `QuasiSteadyCompressibleDuct`: mean ideal density and ideal sonic cap;
- `SeriesDuctOrifice` and historical continuation blends: their component laws;
- `TubeHalfLink`: ideal-density losses, compressible ideal-gas Poiseuille flow,
  pressure-squared form, turbulent loss reconstruction and ideal sonic cap;
- microtube diagnostics/transport: ideal density, dilute-gas sound speed,
  pressure/mean-free-path and Mach/Knudsen assumptions.

`require_ideal_hydraulics` guards historical production flow closures and the
model composition. Only the existing ideal gas or a **verified ideal-reference
table** may use them. This validation does not certify real-helium choking.
`StateHydraulicFlowModel.flow_from_states(first, second, fluid, one_way=...)`
provides the explicit extension boundary. It receives reconstructed states,
including enthalpy and available sound speed. The dynamics can dispatch such
closures while retaining passive valve topology and conservative energy
transport. No production real-fluid flow law or compiled real-fluid hydraulic
kernel is supplied. Unsupported compiled combinations use explicit fallback,
not a relabelled ideal-gas equation.

## Bounded measurements and next scientific work

Run on the local machine, independently of optimization:

```sh
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 examples/benchmark_research_v3_backends.py --calls 10000 --output outputs/my_benchmark.json
```

The same rank01 warm state and crank angle are used for ideal compiled,
tabulated compiled and Python calls. First/JIT calls are separate; steady time
is the median of three fixed-count repetitions. The measured timing and error
are in [the benchmark artifact](../outputs/research_v3/backend_benchmark.json).
The final local measurement was:

| Full RHS call (including the same kinematics) | First/JIT call | Steady call |
|---|---:|---:|
| Ideal compiled | 6.023 s | 37.68 µs |
| Tabulated compiled | 2.204 s | 42.40 µs |
| Python reference | 0.000786 s | 472.00 µs |

The tabulated/ideal steady ratio was **1.125**. First calls are sequential in one
process and share compiler initialization; their times are not independent cold
compiler benchmarks. A separate five-by-30,000-call comparison measured 38.38 µs
for the exported pre-V3 ideal backend and 39.64 µs for V3 (ratio 1.033).
Run-to-run scheduling and thermal/cache variation apply to these small timings.
`--ideal-only --calls 30000` supports reproducing that comparison with each
source tree selected through `PYTHONPATH`. The frozen source commit is recorded
in `tests/data/pre_v3_air_wall.json`; no remote repository is involved.
These are local observations, not universal performance claims. Full-cycle
ideal/table comparison and periodic conservation are separately recorded in
[the acceptance summary](../outputs/research_v3/validated/summary.json).

For each human-cell study: choose the actual external-loop scenarios,
conductance evidence/uncertainty, desired motion, loads, and whether to maximize
COP or cooling power under an indicated input bound; determine how measured
shaft losses will be represented. Those choices belong to the study, not the fluid interface.

Before a cryogenic helium campaign: obtain a traceable single-phase EOS/caloric
dataset; validate positive energy reference, inversion, interpolation error and
phase-domain masks; implement compatible hydraulic closures and their compiled
kernels; supply cryogenic transport data and heat-transfer correlations; validate
conservation and physical reference cases. V3 imports no helium database and
makes no claim of cryogenic accuracy or two-phase capability.

## Historical ideal-gas screening

The early helium/argon comparison used calorically perfect constants, constant UA
and ideal orifices, not the present complete microtube model. Approximate constants
were He `R=2077.1`, `Cp=5193.0`, `Cv=3115.9` and Ar `R=208.13`, `Cp=520.33`,
`Cv=312.20`, in J/(kg K). They are screening inputs, not a real-gas EOS.

For two monatomic ideal gases at identical pressure, temperature, geometry, speed
and UA, scaling orifice CdA by `sqrt(R_reference/R_target)` preserves the
normalized first-level cycle because `Cv/R=3/2`. This does not preserve mass or
prove similarity of real transport, leakage or manufactured components. The
recorded one-bar, 250.15/310.15 K matched check gave He/Ar COP 1.49225/1.49229,
with cooling 176.40/176.44 W; the small difference reflects rounded constants and
tolerances. These historical figures do not certify a modern helium machine.

Higher ideal-gas pressure can scale power at proportionally changed hardware;
it is not by itself a COP improvement. Helium's transport advantages, stronger
adiabatic excursions, containment and nonideality require separate assessment.
Compare indicated cycle COP only with a like-for-like thermal/work boundary;
an appliance's annual electrical label includes a different system boundary.

Retained primary references from the original screening:

- [NIST helium gas properties](https://www.nist.gov/pml/sensor-science/fluid-metrology/database-thermophysical-properties-gases-used-semiconductor-9).
- [NIST Technical Note 1334](https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1334.pdf).
- [NIST helium EOS](https://www.nist.gov/publications/equation-state-thermodynamic-properties-helium).
- [NIST Chemistry WebBook: argon](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440371&Mask=5).

Current transport formulas, oracle data and licensing are documented in the
[microtube reference](MICROTUBE_GAS_MODEL.md); the numbers above do not override them.
