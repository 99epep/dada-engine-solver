# Working fluids: thermodynamics, transport and compiled tables

## Thermodynamics and transport are separate

The production fluid remains `CaloricallyPerfectGas`: P = ρRT, u = CvT,
h = CpT, constant Cp/Cv/gamma and Z = 1. It implements `ThermodynamicFluid`:

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
model. Species-dependent transport domains are documented in
[MICROTUBE_GAS_MODEL.md](MICROTUBE_GAS_MODEL.md); they are not phase/EOS validity
domains. No real-fluid helium EOS, two-phase model or cryogenic hydraulic
validation is introduced.

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
- `SeriesDuctOrifice` and continuation blends: their component laws;
- `TubeHalfLink`: ideal-density losses, compressible ideal-gas Poiseuille flow,
  pressure-squared form, turbulent loss reconstruction and ideal sonic cap;
- microtube diagnostics/transport: ideal density, dilute-gas sound speed,
  pressure/mean-free-path and Mach/Knudsen assumptions.

`require_ideal_hydraulics` guards production flow closures and the
model composition. Only the existing ideal gas or a **verified ideal-reference
table** may use them. This validation does not certify real-helium choking.
`StateHydraulicFlowModel.flow_from_states(first, second, fluid, one_way=...)`
provides the explicit extension boundary. It receives reconstructed states,
including enthalpy and available sound speed. The dynamics can dispatch such
closures while retaining passive valve topology and conservative energy
transport. No production real-fluid flow law or compiled real-fluid hydraulic
kernel is supplied. Unsupported compiled combinations use explicit fallback,
not a relabelled ideal-gas equation.

## Requirements for future real-fluid support

Real-fluid or cryogenic helium support would require:

- a traceable EOS/caloric dataset with a consistent energy reference;
- validated initialization inversions, interpolation errors and phase-domain masks;
- compatible hydraulic closures and, for compiled execution, matching kernels;
- cryogenic transport data and applicable heat-transfer correlations;
- conservation checks and validation against physical reference cases.

The current table interface supplies no helium database, cryogenic accuracy claim
or two-phase capability. See [validation](validation.md) for the distinction
between numerical verification and experimental evidence, and
[solver performance](SOLVER_PERFORMANCE.md) for backend benchmarking methodology.
