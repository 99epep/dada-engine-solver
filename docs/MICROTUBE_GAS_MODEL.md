# Instantaneous single-phase gas model for microtube exchangers

The current motion/hardware search using this closure is documented in
[FOUR_STAGE_OPTIMIZATION.md](FOUR_STAGE_OPTIMIZATION.md). Its results must be
kept separate from the earlier constant-property campaigns. The planned
[temperature map](MOTOR_RESEARCH_OBJECTIVES.md) adapts the whole machine.

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
closures with an invalid domain verdict. It does not create a transition model.

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

The general implementation guard is 200–1000 K for transport, an ideal dilute
gas, smooth circular tubes and single-phase flow. No universal pressure ceiling
is inferred from dilute-gas sources: nonideality at higher density remains a
separate thermodynamic limitation. The local safeguards Ma <= 0.3 and
`2*abs(p1-p2)/(p1+p2) <= 0.2` are conservative **project screening thresholds**,
not new empirical correlations from Yang. Kn < 0.001 is the default continuum
thermal domain. All sources below apply to steady flow unless explicitly stated.
Unsteady use is quasi-steady and unvalidated for pulse phase response.

### GAS-TRANSPORT-SUTHERLAND

- Equation: `x(T)=x0*(T/T0)^1.5*(T0+S)/(T+S)`, separately for mu and k.
- Meaning: dilute-gas transport temperature dependence, not a fitted DADA loss.
- Geometry/boundary: bulk property, geometry independent; air, N2 and Ar.
- Re/Pr/Ma/Kn: not property-law fit coordinates. Pressure: dilute-gas limit.
  Temperature: implementation restricted to 200–1000 K; this is an application
  envelope, not a universal accuracy guarantee from the table.
- Type: analytical approximation with tabulated coefficients; steady property.
- Source/location: [COMSOL Sutherland documentation](https://doc.comsol.com/6.3/doc/com.comsol.help.cfd/cfd_ug_fluidflow_high_mach.08.43.html),
  equations 5-9/5-10, tables 5-2/5-3. Coefficients are explicit in code.
- Implementation: `gas_transport.DiluteGasTransport.viscosity/conductivity`.
- Limits: density effects, mixtures other than the stated air approximation,
  and high-temperature chemistry are not resolved.

### GAS-CP-SHOMATE and GAS-HE-NIST

- Equations: `cp_molar=A+B*t+C*t^2+D*t^3+E/t^2`, `t=T/1000`;
  mass cp divides by molar mass. Ar/He use monatomic `cp=2.5*R`.
- Fluids: [N2 NIST Shomate table](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Mask=1A8F),
  100–500 / 500–2000 K branches; [O2 table](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=11),
  100–700 / 700–2000 K branches. Air transport cp uses a disclosed approximate
  79/21 mole N2/O2 mixture. The solver's calorically perfect cp/cv are unchanged.
- He mu and k: piecewise-linear interpolation of the
  [NIST helium dilute-gas table](https://www.nist.gov/pml/sensor-science/fluid-metrology/database-thermophysical-properties-gases-used-semiconductor-9),
  200–1000 K rows, columns eta/lambda; provenance Hurly & Moldover (2000).
  Table uncertainty is not the interpolation error.
- Meaning/type: thermochemical reference fits and calculated reference transport
  data, not microtube experiments. Geometry/boundary/Re/Pr/Ma/Kn: not applicable
  to the property fits; pressure: ideal/dilute gas.
- Implementation: `gas_transport.DiluteGasTransport.cp/_helium`.
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
  orifice cap remains. Axial acceleration, hydrodynamic entrance pressure
  losses, roughness and nonisothermal axial fields are unresolved.

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

### GAS-NU-HAUSEN

- Equation: `Nu_bar=3.66+0.0668*Gz/(1+0.04*Gz^(2/3))`,
  `Gz=Re*Pr*D/L`; `h=Nu_bar*k(T)/D`.
- Meaning: mean thermally developing laminar coefficient for a circular tube
  with constant wall temperature and a developed velocity profile.
- Source: Hausen approximation to the Graetz problem; equation reproduced in
  [Rastan et al. research manuscript](https://repository.up.ac.za/bitstream/2263/80636/1/Rastan_Heat_2020.pdf),
  Eq.18 (not its modified heat-flux Eq.19). The constant 3.66 is the
  fully developed asymptote, not a universal microtube value.
- Fluid/ranges: conventional single-phase continuum; Re < 2300; code guard
  0.5 <= Pr <= 2000; low Ma/relative pressure drop as above, Kn < 0.001,
  properties in 200–1000 K scope. These code guards do not assert experimental
  validation across that full parameter box.
- Type: steady analytical-solution approximation, applied quasi-steadily.
  Implementation: `gas_correlations.laminar_entry_nusselt`, `gas_film.MicrotubeGasFilm`.
- Limits: `L<0.05*Re*Pr*D` is diagnosed, not rejected solely for thermal entry.
  Hydrodynamic entry (`L<0.05*Re*D`) remains an explicit unsupported-domain flag.
  No axial gas/metal field, thermal-jump correction or pulsed h multiplier.

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
  Re 2300–4000 is unavailable. A continuous friction bridge is used only to
  bracket the algebraic root; a root in that interval is rejected.
- Pressure/temperature: no additional validated gas-specific absolute range.
  [Yang et al. 2014](https://doi.org/10.1016/j.ijheatmasstransfer.2014.07.017),
  pp.732–740, stainless 750/510/170 micrometre tubes, Re 3000–12000, supplies
  gas evidence and warns about compressibility. Its boundary is imposed heat
  input, not the DADA lumped-wall boundary. No quantitative Yang enhancement
  was implemented from the abstract. Rough-tube extensions remain unavailable.

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
