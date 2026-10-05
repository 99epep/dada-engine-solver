# Instantaneous single-phase gas model for microtube exchangers

Reference for geometry, internal transport, correlations, domains and diagnostics.
For study setup see [Research configuration](DADA_ENGINE_RESEARCH_REFERENCE.md);
for finite external streams see [the thermal boundary](EXTERNAL_STREAM_THERMAL_MODEL.md).

## Reading map

- [Circular bundle and conical collectors](#circular-bundles-and-conical-collectors-circular_triangular_frustum_v1).
- [Species-dependent transport](#species-dependent-dilute-transport-version-dilute_species_v2).
- [Equation provenance](#equation-provenance-and-applicability-ledger).
- [Transition closure](#gas-transition-endpoint-interpolation-2026-09-29).
- [Diagnostics](#diagnostics-and-interpretation) and [rejected trials](#first-rejected-trial-state-diagnostics).

## Architecture and selection

The original geometry (`MicrotubeBank`), hydraulic links (`TubeHalfLink`),
thermal assembly (`build_exchanger` / `AirWallExchanger`), family assembly
(`MicrotubeExchanger`) and Doty measurement files remain in use. There is no
second exchanger assembly or change to the conservative thermodynamic states,
valves, indicated-work sign, external-source efficiency boundary or periodic
convergence criterion.

The production internal-gas closure is selected explicitly:

```python
from dataclasses import replace
from dada_solver.exchangers.gas_correlations import MicrotubeGasModel

production = replace(exchanger, inputs=replace(
    exchanger.inputs, gas_model=MicrotubeGasModel()))
```

The hardware TOML loader also supports the following frozen campaign input:

```toml
[gas_model]
mode = "variable_properties"
species = "air"
domain_policy = "reject"
thermal_entry = true
```

Supported species: air, nitrogen, argon, helium. Omitting this section preserves
historical constant-property **legacy/screening** behavior and its regressions.
Constant `gas_nusselt`, viscosity and conductivity fields remain for legacy
compatibility/reference metadata; they do not control the selected production
internal film. External-air properties/film are unchanged and still explicit
screening inputs. This change must not be described as a new external-air model.

`HardwareInputs.gas_model` is immutable and serializable through dataclasses.
The hardware TOML contents already participate in campaign identity/snapshots.
Accommodation and optional slip coefficients can be supplied as nested
`gas_model.accommodation` and `gas_model.slip` TOML tables. Unknown accommodation
is represented by omitted values / Python `None`, never by an inferred metal
TMAC. Production defaults reject unsupported states with `MicrotubeDomainError`;
campaigns classify these as `invalid_exchanger`, not an integration failure.
The explicit `domain_policy="report"` allows only available extrapolative
closures with an invalid domain verdict. Transition has its own explicit closure
under either policy; report mode does not relax its underlying domain guards.

## State, flow and thermal coupling

Transport viscosity is evaluated at each link's upstream temperature; direction
changes select the opposite upstream state. Thermal transport uses the lumped
exchanger gas temperature, not a reservoir temperature. Each of the two port
flow estimates owns half the internal surface. Their Nusselt values are averaged
by area; the development length is the full physical tube length, never a new
thermal entrance at the storage node. The wall resistance is then placed in
series. This is a disclosed lumped approximation; there is no axial temperature
field or resolved conjugate wall conduction.

`AirWallMotor.thermal_rates(angle, state)` supplies the actual instantaneous
film to integration and diagnostics. A contextual wall refuses `rates()` without
its flow context, preventing silent substitution of the old static UA. At a
closed valve with zero flow, the blocked port pressure jump is not treated as
an axial tube pressure gradient: the stagnant tube is evaluated at its gas-node
pressure. The raw blocked-port ratio is retained in the context separately.

During zero flow the existing radial heat-storage approximation remains active,
with the Nu=3.66 limiting conductance labelled `stagnant_radial_screening`.
This is **not** experimental validation of a steady forced-convection coefficient
at zero flow. A radial transient temperature field would require additional
states. Heat is not switched off during pauses and no empirical frequency
multiplier is introduced.

## Equation provenance and applicability ledger

The property-temperature domains are air 100–1000 K, nitrogen/argon 200–1000 K,
and helium 50–1000 K. The model assumes dilute gas and smooth circular tubes.
These property domains are not single-phase or EOS validity domains. No universal pressure ceiling
is inferred from dilute-gas sources: nonideality at higher density remains a
separate thermodynamic limitation. The local safeguards Ma <= 0.3 and
`2*abs(p1-p2)/(p1+p2) <= 0.2` are conservative **project screening thresholds**,
not new empirical correlations from Yang. Kn < 0.001 is the default continuum
thermal domain. All sources below apply to steady flow unless explicitly stated.
Unsteady use is quasi-steady and unvalidated for pulse phase response.

### GAS-TRANSPORT-SUTHERLAND

- Equation: `x(T)=x0*(T/T0)^1.5*(T0+S)/(T+S)`, separately for mu and k.
- Meaning: dilute-gas transport temperature dependence, not a fitted DADA loss.
- Geometry/boundary: bulk property, geometry independent; N2 and Ar. Historical air used this law; current air uses the dilute correlation below.
- Re/Pr/Ma/Kn: not property-law fit coordinates. Pressure: dilute-gas limit.
  Temperature: implementation restricted to 200–1000 K; this is an application
  envelope, not a universal accuracy guarantee from the table.
- Type: analytical approximation with tabulated coefficients; steady property.
- Source/location: [COMSOL Sutherland documentation](https://doc.comsol.com/6.3/doc/com.comsol.help.cfd/cfd_ug_fluidflow_high_mach.08.43.html),
  equations 5-9/5-10, tables 5-2/5-3. Coefficients are explicit in code.
- Implementation: `gas_transport.DiluteGasTransport.viscosity/conductivity`.
- Limits: density effects, mixtures other than the stated air approximation,
  and high-temperature chemistry are not resolved.

### GAS-CP-SHOMATE

- Equations: `cp_molar=A+B*t+C*t^2+D*t^3+E/t^2`, `t=T/1000`;
  mass cp divides by molar mass. Ar/He use monatomic `cp=2.5*R`.
- Fluids: [N2 NIST Shomate table](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Mask=1A8F),
  100–500 / 500–2000 K branches; [O2 table](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=11),
  100–700 / 700–2000 K branches. Air transport cp uses a disclosed approximate
  79/21 mole N2/O2 mixture. The solver's calorically perfect cp/cv are unchanged.
- Historical He mu/k used NIST table interpolation (200–1000 K); the current
  dilute Arp / Hands-Arp formulas below replace it.
- Meaning/type: thermochemical reference fits and calculated reference transport
  data, not microtube experiments. Geometry/boundary/Re/Pr/Ma/Kn: not applicable
  to the property fits; pressure: ideal/dilute gas.
- Implementation: `gas_transport.DiluteGasTransport.cp`.
- Limits: temperature range checked; no extrapolation. cp(T) affects Pr and
  diagnostic sound speed, not conservative internal energy or enthalpy transport.

### GAS-POISEUILLE-ISOTHERMAL

- Equation: `m_dot=N*pi*D^4*(pin^2-pout^2)/(256*mu*L*R*T)`.
- Meaning: signed compressible tube flow; smooth circular tubes, ideal gases,
  constant cross-section and isothermal fully developed laminar velocity field.
- Re < 2300; Pr irrelevant to this isothermal momentum relation; low Ma;
  no-slip continuum, or explicit separately justified slip factor below.
  Pressure and temperature must be positive; properties retain their domains.
- Source/type: analytic Navier–Stokes solution; [Ewart et al., CFM 2007](https://hal.science/hal-03361781v1),
  supplied `HAL_microtubes.pdf`, printed p.3 Eq.3, writes the equivalent
  `pi*D^4*delta_p*pm/(128*mu*R*T*L)` with `pm=(pin+pout)/2`.
- Implementation: `gas_correlations.compressible_poiseuille`, called by
  `hardware.TubeHalfLink._gas_flow` with **L/2** for each of the two links.
- Important compatibility result: the old Poiseuille closure evaluated with
  arithmetic mean pressure density is algebraically identical at fixed mu/T.
  Changing notation alone does not add a compressibility correction.
- Header/valve terms remain separate mean-density losses:
  `delta_p_header=K*m_dot^2/(4*rho*A^2)` per half link;
  `delta_p_valve=m_dot^2/(2*rho*CdA^2)` on the outlet. These are the existing
  project component assumptions, not Graur measurements. The compressible
  orifice cap remains. Continuum laminar entrance pressure losses now use the
  cumulative Shah correction below. Axial acceleration, roughness and
  nonisothermal axial fields remain unresolved.

### GAS-KNUDSEN-HS and GAS-SLIP-SECOND-ORDER

- Equations: `lambda=mu/p*sqrt(pi*R*T/2)`, `Kn=lambda/D` in production;
  `S=1+8*A1*Kn_m+16*A2*(P+1)/(P-1)*ln(P)*Kn_m^2`, `P=pin/pout`.
  The last pressure factor tends to 2 at P=1 and is evaluated with log1p.
- Source/type: [Ewart et al.](https://hal.science/hal-03361781v1), printed
  pp.3–5, Eq.1–5, Tables 1–3; analytical slip model fitted to steady isothermal
  N2/Ar/He flow in **fused silica**, not steel. Their experimental Kn uses VHS;
  the hard-sphere alternative is stated in their background theory.
- Table 1 P=5 fits: Aexp/Bexp are N2 11.668/16.626, Ar 13.218/24.274,
  He 10.812/9.156, with `A1=Aexp/8` and
  `A2=Bexp/(16*1.5*ln(5))`. Table 2 Kn ranges are 0.003–0.291,
  0.003–0.302 and 0.009–0.309 respectively; Table 2 second-order TMAC values
  are 0.908, 0.871 and 0.914. Uncertainties are preserved in the validation report.
- Geometry/boundary: circular silica tube; isothermal wall, velocity slip.
  Re: viscous laminar; Pr: not a hydraulic fit coordinate; Ma: low-speed
  continuum approximation. Absolute pressure/T are not reconstructed from
  the Kn-only fit tables. Pressure-ratio-5 regression checks stay at P=5.
- Implementation: `SecondOrderSlip.factor`, `graur_silica_fit`,
  `GasSurfaceAccommodation`, `MicrotubeGasModel.slip_factor`.
- Guards: fitted VHS coefficients cannot be attached directly to the HS
  production transport. No TMAC-to-A1/A2 conversion is invented. An explicitly
  justified HS coefficient pair/provenance may be configured for a gas/surface.
  Below Kn_m=0.001 correction is omitted; otherwise missing coefficients cause
  rejection by default. Above Kn=0.1 the general production continuum domain is
  rejected even though the particular silica fits extend farther. Thermal slip
  and temperature jump remain unavailable; merely supplying thermal accommodation
  does not silently enable a thermal correction.

### GAS-NU-BENNETT-COMBINED-ENTRY

Production continuum laminar heat transfer uses Bennett (2020a), average
constant-wall-temperature circular-tube closure, throughout the supported
laminar range. The exact equations, source oracle and applicability are in
[Combined laminar entry](#combined-laminar-entry-2026-10-04) below. Thermal
conductance always uses the **full physical tube length**. The storage node
at the tube midpoint is not a fresh entrance.

`hydrodynamic_developing` means `L < 0.0565 Re D`; it is informational in the
supported laminar regime. `thermal_developing` retains `L < 0.05 Re Pr D`.
Neither condition alone invalidates a supported laminar state.

Hausen remains a public reference function and the **historical** transition
endpoint (replaced by Bennett in the continuous revision below):
`Nu_bar=3.66+0.0668*Gz/(1+0.04*Gz^(2/3))`, `Gz=Re Pr D/L`.
It assumes a developed velocity profile and constant wall temperature.
It is no longer the production laminar closure. No thermal-jump correction,
axial transient field or empirical pulse multiplier is introduced.

### GAS-TURBULENT-HAALAND-GNIELINSKI

- Equations: smooth Darcy `f=(-1.8*log10(6.9/Re))^-2`;
  `Nu=(f/8)*(Re-1000)*Pr/[1+12.7*sqrt(f/8)*(Pr^(2/3)-1)]`.
- Source/type: Haaland (1983), J. Fluids Engineering 105, pp.89–90,
  [DOI 10.1115/1.3240948](https://doi.org/10.1115/1.3240948);
  Gnielinski (1976), Int. Chem. Eng.16, pp.359–368,
  [permanent record](https://ndlsearch.ndl.go.jp/en/books/R100000136-I1570854174985140864).
  These are established steady engineering correlations, not new microtube fits.
- Implementation: reuse `models._darcy_friction_factor/_nusselt_number` through
  `gas_correlations.darcy_smooth/gnielinski`; circular hydraulically smooth tube,
  single-phase fluid, locally constant wall temperature.
- Domain: existing conservative 4000 <= Re <= 5e6, 0.5 <= Pr <= 2000,
  L/D >= 10, low Ma and relative pressure drop, Kn < 0.001, property scope above.
  Re 2300–4000 now uses the explicit transition bridge described below.
  The fully turbulent closure and its domain are unchanged.
- Pressure/temperature: no additional validated gas-specific absolute range.
  [Yang et al. 2014](https://doi.org/10.1016/j.ijheatmasstransfer.2014.07.017),
  pp.732–740, stainless 750/510/170 micrometre tubes, Re 3000–12000, supplies
  gas evidence and warns about compressibility. Its boundary is imposed heat
  input, not the DADA lumped-wall boundary. No quantitative Yang enhancement
  was implemented from the abstract. Rough-tube extensions remain unavailable.

### GAS-TRANSITION-ENDPOINT-INTERPOLATION (revised 2026-10-04)

The user-authorized transition closure removes the unavailable-correlation gap,
not the laminar Reynolds limit. For `2300 <= Re < 4000`, use
`w=(Re-2300)/1700` and
`Nu=(1-w)*Nu_Bennett(2300,Pr,D/L)+w*Nu_Gnielinski(4000,Pr)`.
When thermal entry is explicitly disabled, the laminar endpoint remains 3.66.
Below this interval Bennett now replaces Hausen; above it smooth Gnielinski is unchanged.
The original transition entry guards are retained.
`correlation_id=bennett_gnielinski_transition_interpolation` identifies the revised bridge.

New primary source: V. Gnielinski, *On heat transfer in tubes*, International
Journal of Heat and Mass Transfer **63** (2013), 134–140,
[DOI 10.1016/j.ijheatmasstransfer.2013.04.015](https://doi.org/10.1016/j.ijheatmasstransfer.2013.04.015).
The publisher's abstract explicitly recommends linear interpolation of Nusselt
values at 2300 and 4000. We apply that principle to the project's existing
endpoint equations; this is not an implementation of every endpoint correction
in that paper. The earlier local 1976 reference alone did not document this
transition construction.

The hydraulic engineering bridge includes the **complete axial segment**
laminar endpoint. For segment `[x1,x2]` with `ell=x2-x1` and core multiplier M:

```
Delta_K_2300 = K_entry(x2/(D*2300)) - K_entry(x1/(D*2300))
f_low = 64/2300 + Delta_K_2300 * D/(M*ell)
f_transition = (1-w)*f_low + w*f_Haaland(4000)
Delta_p_segment = M * f_transition * ell/D * mdot^2/(2*rho_mean*A_flow^2)
```

Header/legacy-valve terms are added separately. Including M in the denominator
of the entrance part keeps that correction independent of the core multiplier,
just as in the laminar closure. The segment excesses still telescope, including
inside transition. The no-slip production path matches both endpoint pressure
losses continuously. The separately retained legacy slip bridge is not evidence
for developing slip flow; thermal slip still rejects such states.

This is an explicit engineering interpolation, not a measured developing-
transition correlation. The full 2300–4000 interval is used. The supported
laminar and turbulent ranges have **no overlap**; Bennett and Shah are evaluated
only at the laminar endpoint, never extrapolated through transition. Continuity
is C0; endpoint derivatives need not match. No unsupported C1/experimental
smoothness claim or empirical pulse correction is made.

Both thermal endpoints must be available: the transition retains the laminar
endpoint's developed-velocity requirement `L/D >= 0.05*2300` and the turbulent
endpoint's `L/D >= 10`. Both thermal endpoints require `0.5 <= Pr <= 500`
(the intersection with Bennett's supported Pr range). Mach, relative pressure drop,
Knudsen/slip and transport-temperature limits are unchanged. Transition is not
an excuse to accept a pressure drop above the configured limit.

`flow_regime` distinguishes laminar, transition and turbulent. To keep Research
constraint semantics compatible, `model_validity` remains `valid` / `invalid`;
`valid` means the declared model domain is respected, not experimental certainty.
Cycle diagnostics add `transition_model_used`, `confidence=transition_uncertainty`,
transition time fraction, absolute-heat fraction and per-passage transition
Reynolds ranges. Per-passage `domains.flow_regime` also reports time and heat
fractions of each thermal regime. Transition usage includes either the thermal
or upstream hydraulic Reynolds state; overall time is the union across passages,
not the sum. Heat is allocated by the existing instantaneous half-film conductance
weights. Zero total absolute heat gives null heat fractions. These are sampled,
trapezoidal physical-time diagnostics, not exact event-duration measurements.

The compiled kernel supports no-slip laminar, transition and turbulent states
using the same friction and Nusselt formulas as Python. Above the laminar branch,
a bracketed bisection solves the existing pressure-loss balance inside Numba to
floating-point precision (the Python reference uses Brent with `xtol=1e-15`).
The ideal and tabulated validation kernels share these primitives. Unsupported
or out-of-domain states still fall back to the authoritative Python path at the
same state, preserving rejection messages and diagnostics. Backend statistics
retain `fallback_calls` and `unsupported_state`. No physical threshold, transition
uncertainty classification or periodic convergence setting changes. See
[the execution benchmark](SOLVER_PERFORMANCE.md#compiled-transition-benchmark).

This remains a quasi-steady approximation in a pulsed machine. Transition onset,
intermittency, inlet disturbances, hysteresis and micro/mini-channel surface effects
are not resolved. A candidate with substantial transition usage needs experimental
validation; neither the interpolation nor its continuity establishes accuracy.
Existing histories are not rewritten. Changed source/runtime identity prevents
silently resuming or reusing results from the previous physical closure.

### GAS-UNSTEADY-DIAGNOSTICS

- Equations: `Wo=(D/2)*sqrt(2*pi*f*rho/mu)`, `St=f*L/abs(u)`;
  `tau_nu=(D/2)^2*rho/mu`, `tau_alpha=tau_nu*Pr`,
  `tau_res=L/abs(u)`, `tau_acoustic=L/sqrt(gamma*R*T)`.
- Meaning/type: dimensionless scaling and characteristic times from momentum,
  thermal diffusion, residence and acoustic propagation; not fitted corrections.
  St and residence time are unavailable at zero flow instead of serialized infinity.
- Geometry/fluid/boundary: circular gas passages, any selected wall closure;
  Re/Pr/Ma/Kn and pressure/T inherit instantaneous model diagnostics. These
  definitions have no experimental validity guarantee for the pulse response.
- Source context: Aubert (1999), supplied `50015150.pdf`, chapter III §3.4.3.2,
  PDF pp.181–186 and annex figures A13–A22. This is experimental **pressure gain
  and phase**, not h(t). PDF p.183 Fig.III.49/50 includes 76.25 micrometre,
  30.13 mm tubing near 113 kPa. Fundamental-frequency Wo alone cannot resolve
  rapid pulse edges or higher harmonics.
- Implementation: `MicrotubeGasModel.diagnose`. No unsteady thermal multiplier.

## Diagnostics and interpretation

`MicrotubeFlowDiagnostics` contains Re, Pr, Ma, local maximum Kn, mean Kn,
Gz, pressure ratio, compressibility indicator, entry flags, continuum/slip
classification, correlation ID, Nu and timescales. Hydraulic extrema use the
actual upstream state; thermal domains use lumped exchanger gas conditions.
`cycle_microtube_diagnostics` reports each inlet/outlet separately with extrema,
trapezoidal physical-time fractions and fractions of **absolute wall-to-gas
heat**. Adaptive solver sample counts are not used as time fractions.
A valid steady closure does not certify pulse response, roughness, flow
maldistribution, axial conjugate heat transfer or the external-air film.

### GAS-STATE-DIMENSIONLESS

- Equations: `rho=p/(R*T)`, `u=abs(m_dot)/(rho*N*pi*D^2/4)`,
  `Re=rho*u*D/mu`, `Pr=cp*mu/k`, `Ma=u/sqrt(gamma*R*T)`,
  `gamma=cp/(cp-R)`, `Gz=Re*Pr*D/L`.
- Meaning/type: ideal-gas equation of state and dimensionless definitions,
  not additional empirical fits; source is the ideal-gas balance already
  implemented in `fluids.CaloricallyPerfectGas` and the cited Graetz/slip models.
- Geometry: circular parallel tubes with equal flow splitting. Fluids, pressure,
  temperature and boundary restrictions inherit the chosen transport/closure;
  the definitions themselves impose no additional Re/Pr/Ma/Kn range.
- Implementation: `MicrotubeGasModel.diagnose`. Outlet/minimum pressure gives
  the conservative local maximum Ma and Kn; mean pressure gives Kn_m for slip.
  No local axial temperature profile is inferred.


## Research V3 external-stream and fluid boundaries

The internal gas film, tube hydraulics and gas-transport domains remain the same.
`ExternalStreamMicrotubeExchanger` connects these components to a declared
finite-capacity stream through explicit wall conductance. Air is one scenario;
liquid labels do not select unvalidated correlations. Wall storage and pause
heat transfer remain active in motor and refrigeration operation. The legacy
`AirWallExchanger` inputs/results remain supported. See
[external-stream equations and sign conventions](EXTERNAL_STREAM_THERMAL_MODEL.md).

Thermodynamic EOS/caloric reconstruction now has a separate conservative-state
interface and a compiled rho/u table validation path; this does not change the
role of transport Cp. The validation table is analytically ideal and may use
the existing ideal-density/constant-gamma hydraulic laws. Those laws explicitly
reject an incompatible real-fluid model. Microtube density, sonic caps,
Poiseuille reconstruction, Mach and rarefaction diagnostics still need scientific
revalidation for real helium. `DiluteGasTransport` now uses the species-dependent temperature domains below;
there is no cryogenic extrapolation or two-phase extension. See
[working-fluid models and compiled backends](WORKING_FLUID_MODELS.md).

### First rejected trial-state diagnostics

Research wall evaluations rejected during integration retain
`diagnostics.first_microtube_failure`. This failure-only snapshot reads the
existing Python exception traceback; it does not replay the cycle or evaluate a
correlation again. The compiled wall RHS already falls back to the Python RHS
at the same unsupported trial state, so the same instrumentation covers both
backends. The exception message and rejection status remain unchanged.

The snapshot records the exact criterion, available flow numbers, passage,
exchanger (`heat_in` / `heat_out` identify configured sides, not operating-mode labels),
pressures, gas temperature and tube geometry. Angle and time refer to the local
cycle's rejected RHS trial, which is not necessarily an accepted solver sample.
Thermal flow is signed in the named passage direction; hydraulic flow is in the
upstream-to-downstream direction. Tube length is the full bank length (the
hydraulic link still uses half that length). If hydraulic closure fails before
producing a flow, unavailable fields remain null; an available preliminary
Reynolds estimate is explicitly named `hydraulic_trial_reynolds`.

The first failure remains attached even after a safe-domain retry or an
acceleration rollback; the top-level status still describes the final outcome.
Missing optional frame locals or properties remain null without discarding the
rest of the snapshot. The integrator callback and model method both named
`derivative` are supported, including the callback without a `self` local.
Snapshot extraction is best effort and cannot replace the scientific exception.
Older artifacts are unchanged and do not acquire these fields retrospectively.

### Transition validation result

Boundary/continuity, retained guards, cycle-weighted usage and explicit Numba
fallback tests pass. Related tests: **146 passed**. Complete suite:
**815 passed, 0 warnings in 295.04 s**.
The single bounded Human Cell evaluation at N=1000, D=0.770 mm, L=0.8 m
reaches Re=2393.30 with `gnielinski_transition_interpolation`, then rejects
relative pressure drop 0.2005927 above the unchanged 0.2 guard. This is not a
converged performance result. See [the artifact and reproduction command](../outputs/research_microtube_transition/README.md).


## Research limit ownership (2026-09-30)

Microtube validity uses each selected gas model's declared Mach limit, including
its hydraulic upstream diagnostics. The historical generic validity value 0.2
no longer adds a second exchanger-domain veto to a model allowing 0.3. An explicit
study `maximum_mach_number` constraint may impose a stricter engineering limit.
No transport, pressure-drop, entry, rarefaction or correlation limit is relaxed.
There is no universal 0.08 kg/s admissible flow: geometry and local state determine
velocity, Reynolds, Mach, pressure drop and heat transfer. Peak absolute mass flow
remains an observable and an optional explicit study constraint. Pressure inequality
is diagnostic-only; isothermality excursion is no longer reported or constrained.
See [current physics decisions](PHYSICS_DECISIONS.md).


## Species-dependent dilute transport, version `dilute_species_v2`

This is a scientific property-law change, not an EOS replacement. New air defaults
are 100–1000 K; helium 50–1000 K; nitrogen and argon remain 200–1000 K.
`minimum_temperature` and `maximum_temperature` may restrict, never expand, these
ranges. Serialized old explicit 200 K lower bounds remain restrictive. Old studies
remain readable, but new code uses the new correlation formulas even at 200–1000 K;
source/backend identity prevents silently resuming a historical runtime. Preserve
old outputs and use a new study/basis/campaign, not a rewritten journal.

### Formulas and provenance

Adapted from CoolProp **v8.0.0**, MIT: [Air.json](https://github.com/CoolProp/CoolProp/blob/v8.0.0/dev/fluids/Air.json)
and [TransportRoutines.cpp](https://github.com/CoolProp/CoolProp/blob/v8.0.0/src/Backends/Helmholtz/TransportRoutines.cpp).
Attribution and the complete MIT license are in [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md).
No residual-density or critical enhancement terms are included.

- Air: [Lemmon & Jacobsen (2004), *Viscosity and Thermal Conductivity Equations
  for Nitrogen, Oxygen, Argon, and Air*, IJT 25, 21–69](https://www.nist.gov/publications/viscosity-and-thermal-conductivity-equations-nitrogen-oxygen-argon-and-air).
  With `x=ln(T/103.3)` and `S=exp(0.431-0.4623*x+0.08406*x²+0.005341*x³-0.00331*x⁴)`,
  `mu=2.66958e-8*sqrt(28.9586*T)/(0.36²*S)` Pa s.
  With `tau=132.6312/T`, `k=0.001308*(mu*1e6)+0.001405*tau^-1.1-0.001036*tau^-0.3` W/(m K).
  The reducing temperature is the EOS reducing value, not the separate critical-state value.
- Helium viscosity: [Arp, McCarty & Friend, NIST TN 1334 revised (1998)](https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1334.pdf),
  dilute limit of CoolProp's `viscosity_helium_hardcoded`.
  For `T<=100 K`, `mu=1e-7*exp(-0.135311743/x+1.00347841+1.20654649*x-0.149564551*x²+0.012520841*x³)`,
  `x=ln(T)`. Above 100 K, `mu=1e-7*196*T^0.71938*exp(12.451/T-295.67/T²-4.1249)` Pa s.
- Helium conductivity: [Hands & Arp (1981), *A Correlation of Thermal Conductivity
  Data for Helium*, Cryogenics 21, 697–703](https://doi.org/10.1016/0011-2275(81)90211-3).
  `k=2.7870034e-3*T^0.7034007057*exp(3.739232544/T-26.20316969/T²+59.82252246/T³-49.26397634/T⁴)` W/(m K).
  The paper's reported data range ends at 830 K; the requested 1000 K DADA ceiling
  follows the implemented CoolProp dilute expression, not new experimental validation above 830 K.
- Air transport Cp retains the NIST Shomate 79/21 molar N2/O2 approximation,
  with branches documented from 100 K. Helium retains `Cp=2.5*2077.1` J/(kg K).
  None of these transport Cp values replaces DADA conservative caloric properties.

The oracle's air EOS is [Lemmon, Jacobsen, Penoncello & Friend (2000), JPCRD 29,
331–385](https://www.nist.gov/publications/thermodynamic-properties-air-and-mixtures-nitrogen-argon-and-oxygen-60-2000-k-pressures).
It is **not** installed as DADA's thermodynamic EOS. `property_temperature_domain`
is not `single_phase_domain`: no condensation, mixture phase equilibrium or
real-gas nonideality check is added. A future `(T,p)` phase/nonideality layer remains
necessary. Tube transport is still quasi-steady in a pulsed machine.

### Frozen oracle and measured changes

`tools/generate_coolprop_transport_reference.py` requires optional CoolProp 8.x
only to regenerate `tests/data/coolprop8_dilute_transport.json`. Runtime and normal
tests do not import CoolProp. The oracle uses `(T,Dmass)` at `1e-10 kg/m³` and
checks density reduction to `1e-11 kg/m³` to bound residual contamination.

On the requested 19 temperature points (air 100–1000 K, He 50–1000 K), maximum
relative errors versus CoolProp 8.0.0 are:

| Species | mu | k | Cp |
|---|---:|---:|---:|
| Air | 7.31e-14 | 3.18e-13 | 0.007061 |
| Helium | 3.17e-13 | 3.67e-13 | 0.00007881 |

Tests allow 2e-11 for directly ported mu/k, 0.8% for approximate air Cp, and 0.01%
for helium Cp using unchanged rounded R. These are oracle agreement tolerances,
not experimental accuracy claims.

Compared with the prior formulas on a 1 K grid from 200 to 1000 K:

| Species | mu relative change | k relative change |
|---|---:|---:|
| Air | +0.1566% to +4.0740% | -0.6811% to +3.0224% |
| Helium | -0.5200% to +0.2476% | -0.5684% to +0.2922% |
| Nitrogen / argon | 0 | 0 |

Transport Cp is unchanged for all four species. Frozen legacy samples are retained
in `tests/data/legacy_dilute_transport_v1.json`. The Doty He variable-flow regression
shifts by about +0.248%; its tolerance and experimental comparison remain unchanged.

### Domain failures and campaign isolation

`TransportDomainError` records species, temperature and allowed bounds, with
`transport_temperature_below_domain`, `transport_temperature_above_domain`, or
`transport_temperature_nonfinite`. Research records `invalid_fluid_domain` and
`diagnostics.transport_failure`, including integration/postprocessing phase.
A postprocessing failure retains convergence evidence but exposes no usable
objective or reusable final state. The next candidate is still evaluated.
Only this typed property-domain error is intercepted; programming errors are not
swallowed in the runner. Compiled temperature guards fall back to Python, which
raises the same typed error. No CoolProp property call occurs inside the RHS.

### Historical machine replays under the changed transport law

`tools/generate_transport_v2_thermal_reference.py` performs four bounded fixed-input
DADA evaluations (no search) and writes a separate versioned test reference.
`tests/data/transport_v2_thermal_reference.json` does not replace historical result
artifacts. On those machines the measured indicated-power / efficiency changes are:

| Historical geometry | Periodic cycles | Indicated-power change | Efficiency change |
|---|---:|---:|---:|
| Fourier C2 | 4 | -0.13685% | -0.09601% |
| Free spline | 3 | -0.11267% | -0.10866% |
| Six-bar | 13 | -0.09318% | -0.06692% |
| Structured C2 15p | 3 | -0.11506% | -0.11116% |

The existing numerical comparison tolerances are retained against these new-law
references; comparing to the old-law numbers with roundoff tolerances would test
a different physical model. Historical constant-transport cases are unchanged.

### Triangular tube banks and internal header volume

Use `pitch_ratio > 1` for triangular (hexagonal-neighbour) packing. The outer
diameter is `Do = inner_diameter_m + 2*wall_thickness_m`; the centre pitch is
`pitch_ratio*Do` (`MicrotubeBank.effective_pitch_m`). Diameter changes therefore
preserve non-overlap without a separately chosen absolute pitch. Ratios equal
to one are rejected, including Research lower bounds. This condition guarantees
geometric separation, not a manufacturing tolerance or structural qualification.

Packing is deterministic: `ceil(sqrt(N))` columns, filled row by row. Odd rows
start half a pitch to the right; vertical separation is `sqrt(3)/2*pitch`.
The rectangular frontal envelope spans the **occupied tube centres**, plus one
outer radius at each edge. An incomplete final row does not introduce fictitious
tubes. No packing optimization is performed. Both internal gas plenums use this
full rectangular face: `V_headers = 2*width*height*header_depth_m`. Total gas
hold-up adds tube bore volume and `additional_internal_volume_m3` exactly once.
No external-fluid interstitial volume is calculated. These are internal header
dimensions; vessel walls and extra fabrication clearances are not modeled.

Compatibility is explicit: old `pitch_m` inputs retain the historical square
packing and `columns*pitch_m` by `rows*pitch_m` envelope. Do not supply both pitch
forms. Existing bases, examples and histories are not migrated implicitly.
Selecting `pitch_ratio` changes header hold-up and scientific study identity;
start a new campaign rather than resuming a square-packing history.

### Circular bundles and conical collectors (`circular_triangular_frustum_v1`)

Supply `collector_half_angle_deg` to select the circular model. It requires
`pitch_ratio > 1` and replaces both `pitch_m` and `header_depth_m` as input fields.
The remaining inputs are tube count, length, bore diameter, wall thickness,
`conduit_area_ratio >= 1` (default 1), and nonnegative
`additional_internal_volume_m3` (default 0). The half-angle is explicit and
strictly between 0 and 90 degrees; no cone-angle optimum is assumed.

With `Do = Di + 2*wall_thickness`, the implemented continuous approximation is:

```text
pitch = pitch_ratio * Do
cell_area = sqrt(3)/2 * pitch^2
bundle_face_area = tube_count * cell_area
Db = sqrt(4*bundle_face_area/pi)
tube_flow_area = tube_count*pi*Di^2/4
conduit_area = conduit_area_ratio*tube_flow_area
Dc = sqrt(4*conduit_area/pi)
height = (Db-Dc)/(2*tan(collector_half_angle))
one_collector_volume = pi*height*(Db^2+Db*Dc+Dc^2)/12
header_gas_volume = 2*one_collector_volume
working_gas_volume = tube_flow_area*tube_length + header_gas_volume + additional_internal_volume
```

`Dc >= Db` is rejected explicitly, without an inverted or degenerate cone.
The triangular-cell area is a continuum envelope approximation, not an exact
placement of a finite number of circles; no discrete edge correction is added.
The model is most representative for large tube counts. `core_width_m` and
`core_height_m` remain bounding-box extents equal to `Db`; their product is **not**
the circular face area. The axial fluid envelope is `tube_length + 2*height`.
No external interstitial fluid volume enters the gas inventory. No tube-sheet,
collector wall storage or conduit length is invented.

For this geometry, both `MicrotubeExchanger` and
`ExternalStreamMicrotubeExchanger` derive the reported valve CdA from
`conduit_area`. A supplied historical outlet CdA is replaced, not clamped.
The ideal diode acts only through the production network's direction logic:
no valve orifice loss is passed to either `TubeHalfLink`. In particular,
merely enlarging a CdA would not suffice, because the old link also adds a
quadratic valve resistance. Both that resistance and the smaller-valve sonic
cap are absent. The **tube-area** compressible cap, tube friction, transport
and gas-model validity guards remain. Upstream/downstream valve placement
continues to determine which port blocks reverse flow. Generic orifices and
legacy finite-CdA links are unchanged.

`header_loss_coefficient` remains the existing lumped loss coefficient referred
to total tube-passage velocity, divided between the two links. It is independent
of cone angle and has not become a distribution, separation or pressure-recovery
correlation for conical manifolds. For legacy external-air thermal models only,
the existing equivalent-passage screen now uses the circular face area and
perimeter; no new empirical correlation is claimed. Declared external-stream
conductances are unchanged.

Reports expose the geometry model, bundle area/diameter, actual pitch, total
bore area, conduit area/diameter/ratio, cone angle/height, both header and total
gas volumes, derived valve area, and loss-model scope. These values also appear
in Research candidate comparisons when present in the recorded hardware data.

**Compatibility:** no file or history is silently migrated. Without a cone-angle
input, square-pitch and rectangular staggered legacy envelopes retain their old
header volumes and finite-CdA hydraulic behavior. Selecting the circular model
changes geometry and valve physics, hence scientific identity; use a new study/
campaign. Exact reconstruction of historical inputs is still possible, subject
to the existing runtime compatibility checks.

## Combined laminar entry (2026-10-04)

Scientific inputs authorized for this revision: Bennett (2020a), DOI
[10.1115/1.4047834](https://doi.org/10.1115/1.4047834), average constant-wall-temperature
branch; [author HeatLib](https://sites.me.ucsb.edu/~bennett/heatlib/conv/index.html),
`LamPipeEntryNuT`, `LamEntryNuL`, `GrtzNuL`, `LevNuL`, `mLevNuL`, `InvGz`.
The official [EES PipeFlow_Laminar documentation](https://fchartsoftware.com/ees/heat_transfer_library/internal_flow/hs1024.htm)
identifies Bennett (2020a) and Shah & London Eq.192 p.98 explicitly. No fitted DADA coefficients or runtime external library are used.

Exact scalar transcription for a circular tube (Darcy fully developed fRe=64):

```
Z = L / (D Re Pr)
A = 0.40377
O = (3.66 - 6.54) / 4.35
n = (3.66 + 41.0) / 13.3
g = 1.10 * (1 + 0.140 / Pr^(2/3))^(3/4)
Nu_Lev = A * (64/Z)^(1/3)
Nu_G = (Nu_Lev^n + (3.66-O)^n)^(1/n) + O
fLRe = 5.312 / (g * sqrt(Pr*Z))
Nu_mLev = A * (fLRe/Z)^(1/3)
F = (Nu_Lev/Nu_mLev) * (1 + 0.565*(Pr*Z)^(1/3))
Nu_mean = Nu_G / tanh(F)
```

Use physical full-tube L, not a half-link length. Re<2300, Z>1e-6;
production retains the pre-existing lower Pr guard 0.5 and caps Pr at 500
(the reported benchmark range). EES states about 2.5% accuracy above Pr=0.5;
the extended Pr>0.1 range with errors up to about 10% is not enabled here.
At Re=0 the existing stagnant Nu=3.66 screening remains. Explicit
`thermal_entry=false` retains the existing 3.66 screening option.

Shah & London (1978), *Laminar Flow Forced Convection in Ducts*, Eq.192 p.98:

```
xplus = x / (D Re)
f_F Re = 3.44/sqrt(xplus)
       + (1.25/(4*xplus) + 16 - 3.44/sqrt(xplus))
         / (1 + 0.00021/xplus^2)
f_D = 4*f_F
K_entry(xplus) = 4*xplus*(f_F Re) - 64*xplus
              = [1.25*xplus^2
                 + 0.00021*(13.76*sqrt(xplus)-64*xplus)]
                / (xplus^2+0.00021)
K_entry(0) = 0
Delta_p_entry[x1,x2] = [K_entry(x2/(D Re))-K_entry(x1/(D Re))]
                      * mdot^2/(2*rho_mean*A_flow^2)
```

The second K expression avoids subtracting large fully developed contributions.
The correction alone uses local mean ideal density `(p1+p2)/(2 R T_upstream)`;
the pressure-squared Poiseuille term is preserved. Header K and diode/legacy
valve terms remain separate. Header K describes manifold/contraction/exit losses,
not the velocity-profile development already counted by Shah. No new multiplier
is introduced. This incompressible correction is used only within the existing
low-Mach/low-relative-pressure-drop continuum assumptions, not as a new
compressible entrance solution.


Each `TubeHalfLink` stores its physical axial half (0 or 1). Forward flow uses
`[0,L/2]` and `[L/2,L]`; reverse flow mirrors these coordinates. Their excess
coefficients telescope to the full-tube coefficient for the same flow and
properties. The lumped storage model can have different instantaneous port
flows and temperatures: each segment uses its own local values, without
introducing another entrance or an axial transient state. Shah tends to
`f_D=64/Re` as `xplus` grows; the finite cumulative excess tends to 1.25 and
becomes negligible relative to the length-proportional Poiseuille loss.
The existing slip-flow hydraulic branch is retained; developing slip flow
is not validated by this no-slip correction and thermal slip remains rejected.

### Regime boundaries and remaining limitations

The user explicitly authorized revising transition to join the new laminar
endpoints continuously, after a bounded DD13 replay exposed switching between
the discontinuous old/new branches. The boundaries remain 2300 and 4000, and
all existing entry/Mach/pressure-drop/slip/transport guards remain. There is no
superposition of validated Reynolds domains: the bridge is marked as transition
uncertainty, with its time/heat fractions retained in reports.

Transition retains `L < 0.05*2300*D` rejection and transition/turbulent flows
retain `L < 10*D` rejection. Bennett is never extrapolated into either regime.
Laminar Pr outside 0.5–500 or inverse Graetz <=1e-6 is unsupported. Mach,
pressure-drop, Knudsen/slip, transport and compressibility guards are unchanged.
The equations are steady constant-property correlations evaluated with local
variable gas properties. They do **not** resolve pulse history or transient
axial velocity/temperature fields. Collector losses remain an independent,
uncalibrated coefficient; any future calibration must exclude the profile
entrance loss now explicitly represented.

Both backends call the same pure numerical primitives. The compiled path
supports the new laminar closure directly; fallback remains reserved for
unsupported states. Exchanger metadata records thermal/hydraulic correlation
IDs, density convention, axial segmentation and header-loss scope.

`tests/data/bennett_heatlib_reference.json` contains 50 independent Octave
outputs from the author's HeatLib at Pr=0.5, 0.7, 1, 10 and 500 over short/long
tubes, with source hashes and the documented unused-branch syntax repair.
HeatLib is not a runtime dependency. Shah tests evaluate Eq.192 independently,
check the Darcy limit, segment additivity and forward/reverse ownership.
