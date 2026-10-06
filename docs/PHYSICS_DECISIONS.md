# Current physical decisions

This is the current model handoff, consolidated on 2026-10-03. Explicit later
user decisions override it. Earlier decisions and superseded requirements remain
available in Git history; this document is the current source of truth.

## Conservative state and hydraulic topology

Four independent gas control volumes retain their own mass and internal energy:
small cylinder S, large cylinder L, heat-input exchanger H_i, heat-output
exchanger H_o. Temperature and pressure are reconstructed; pressure equality
between connected volumes is never imposed. The quasi-pressure-equalized
reductions in the study are analytical references only.

The eight gas states are `(m_S,U_S,m_L,U_L,m_Hi,U_Hi,m_Ho,U_Ho)`.
Dynamic-wall models add two wall energies. Heat/work quadratures are accounting
variables, not additional physical storage states. Legacy C/cold and H/hot
names remain readable; they identify branches, not necessarily external cold/hot
reservoirs in both operating modes.

Branch circulation is S → H_i → L → H_o → S. Each exchanger may place its
passive valve upstream or downstream; the other half-link is bidirectional.
The actual pressure gradient selects the upstream transported enthalpy.
Valve timing is not prescribed by piston phases or Lambda targets. Reflux on
bidirectional links remains visible. Unsupported/missing event observations
must not be converted into a non-nominal chronology verdict.

Moving a valve upstream or downstream changes the hydraulic network and can
change pressures, flows and heat transfer. Treat placement as a design variable
and reoptimize it with motion and hardware. Valves remain pressure-driven, not
angle-commanded. A topology preferred in one finite motor campaign is not a
universal choice for another machine, fluid or operating regime.

Continuous ideal diodes impose direction without valve memory. The separate
hysteretic event solver preserves state and quadratures at threshold crossings;
there is no instantaneous pressure equalization/reset. The wall wrapper requires
continuous ideal diodes. Real valve lift, inertia, leakage and friction are not
implicitly modeled.

## Kinematics, direction and scale

The solver receives a `KinematicsModel`: cylinder volumes and angular derivatives,
with exact breakpoints where needed. Motion is prescribed, not obtained from
mechanical force balance. Inject study-angle kinematics before the factory's
single motor-direction transformation. Motor uses negative configured angular
speed; integration still advances positive cycle progress. See
[motor conventions](MOTOR_OPERATION.md).

Small and large motion laws need not be symmetric or of the same family.
Shared-crank mechanisms retain their explicit coupling and discrete branches.
Normalized laws do not imply a physical stroke or piston area without scale.
Piecewise-linear laws are design probes, not realizable optima. Mechanism
closure is intrinsic; transmission-quality requirements belong to a study.
Use the [established synthesis method](MECHANISM_SYNTHESIS_SEARCH.md).

## Heat, work and performance

Gas boundary work is positive out of the gas:
`Wdot = P_S * Vdot_S + P_L * Vdot_L`. Each cylinder uses its own pressure.
External heat is positive into the machine. Both operating modes require
`Q_i > 0` and `Q_o < 0`; refrigeration also requires `W < 0`, motor `W > 0`.
Consuming work alone does not establish refrigeration.

For a periodic refrigerator, cooling COP is `Q_i / (-W)` and heating COP is
`(-Q_o) / (-W)`. Motor indicated efficiency is `W / Q_i`. The gas-plus-wall
energy residual retains storage: `delta(U + E_wall) - Q_i - Q_o + W`.
Indicated input/output is distinct from human shaft input/useful output.
Mechanical efficiency is not assumed. External pump/fan losses excluded by a
study are unknown physical losses, not zero losses.
Compare indicated cycle COP only with a like-for-like thermal/work boundary;
an appliance annual electrical label includes a different system boundary.

The reservoir model uses `UA*(T_source-T_gas)`. The wall model resolves internal
film → wall storage → external conductance → finite-capacity stream. Heat
exchange remains active at zero through-flow and during piston pauses.
Wall-to-gas heat is not interchangeable with external-boundary heat.
See [external streams](EXTERNAL_STREAM_THERMAL_MODEL.md) for equations and signs.

## Geometry and exchanger ownership

Microtube geometry determines gas hold-up, tube areas and tube-wall capacity.
Do not independently optimize those derived quantities as unrelated UA, volume
and loss coordinates. Declared external conductance is a separate scenario input
for the generic external-stream family; it is not a universal liquid correlation.

The new circular family uses triangular-cell face area, `pitch_ratio > 1`,
`conduit_area_ratio >= 1` and two coaxial conical-frustum collectors. Only tube
bores, collector interiors and explicit additional internal volume enter gas
hold-up. The external interstitial space does not. Independent valve CdA is
replaced by the conduit area and a lossless directional diode; tube losses and
sonic/domain guards remain. Header K is still a disclosed lumped assumption,
not a loss correlation derived from cone angle.

Old square/staggered rectangular banks retain their stored geometry and valve
loss assumptions. Selecting circular collectors changes the scientific study;
it is not a transparent resume of an old campaign. Formulas, compatibility and
uncertainty are in [the microtube reference](MICROTUBE_GAS_MODEL.md).

At fixed external capacity rate, increasing transfer area has diminishing benefit;
large air flow suppresses outlet-temperature change but does not remove film
resistance. More, shorter tubes at fixed surface can reduce tube friction while
increasing collector hold-up. H_i and H_o geometries may be independent. Multiplying
UA at unchanged volume and hydraulics is a sensitivity experiment, not physical
resizing. Internal hydraulic losses already affect cycle work and must not be
subtracted again as an additional fictitious pump load.

## Fluids and model domains

Production thermodynamics remains calorically perfect ideal gas. Variable
transport Cp, viscosity and conductivity do not replace the EOS or energy law.
The compiled rho/u property table is an architectural proof using ideal-generated
data, not a validated real-helium model. Ideal hydraulic laws require compatible
fluid assumptions; real-fluid choking requires a separate validated closure.

Dilute transport property-temperature domains are air 100–1000 K, He 50–1000 K,
N2/Ar 200–1000 K. User bounds can narrow them, not expand them. Air uses dilute
Lemmon/Jacobsen, He dilute Arp/Hands-Arp; CoolProp is an offline oracle only.
Property availability is not phase or nonideality validity. No condensation,
two-phase working gas or validated cryogenic real-gas cycle is claimed.

Microtube laminar heat transfer uses Bennett (2020a), average constant-wall-
temperature combined entry, with the full physical tube length. Shah & London
(1978), Eq.192 p.98, supplies only the cumulative entrance excess over the
existing compressible Poiseuille term. The two physical half-segments take
successive differences of that excess, mirrored on reversal, never a new
entrance at the gas-storage node. The excess alone uses arithmetic mean
pressure density at upstream temperature; header losses remain separate.
Exact equations and domains: [microtube ledger](MICROTUBE_GAS_MODEL.md#combined-laminar-entry-2026-10-04).

Supported laminar hydrodynamic development is informational, not a rejection.
The user subsequently authorized a continuous transition over all of
Re=2300–4000: Bennett Nu and the complete Shah-corrected segment friction
at 2300 are linearly joined to the existing turbulent endpoints at 4000.
There is no overlap of validated Reynolds domains; this engineering bridge
retains explicit transition uncertainty and the previous entry guards.
It is not a validated developing-transition correlation. Do not weaken Mach, relative
pressure drop, Knudsen/slip or transport guards. Pulsating use remains
quasi-steady; no empirical pulse multiplier or axial transient field is added.

## Limits: four different owners

| Kind | Meaning | Examples |
|---|---|---|
| Model domain | Applicability of the selected equation/model | Transport range, model-specific Mach/pressure-drop guards, real-valued linkage closure |
| Study/design constraint | Requirement selected for a machine | Minimum power, pressure/temperature/flow ceilings, transmission quality |
| Search bound | Region explored, not a validity statement | Parameter bounds, local-region radius |
| Numerical setting | Accuracy/execution policy | ODE tolerances, periodic scales, sampling |

Generic Research does not impose historical 25 W, 1.2 MPa, 850 K, 0.08 kg/s or
66 L campaign guards. Explicit supported study constraints remain available.
Microtube Mach validity uses the chosen exchanger model's domain; a stricter
study Mach ceiling is separate. Cp variation and compressibility deviation
remain approximation checks. Constant Cp and Z=1 yielding zero deviations do
not independently validate real air at every state.

Pressure inequality remains a diagnostic, never a thermodynamic-validity veto.
Isothermality excursion metrics/constraints have been removed. Old configuration
keys load as ignored compatibility inputs. Historical saved verdicts are not
rewritten; changed source identity prevents unsafe resume/cache reuse.

## Inventory, periodic convergence and diagnostics

Explicit inventory remains supported. The optional reference-pressure policy
uses `m = p_ref * max_theta(V_total(theta)) / (R*T_ref)` for the assembled ideal-gas
machine. Both cylinders are evaluated at the same angle; final exchanger volumes
are counted once. This is not the sum of independent piston maxima, nor a
requirement that every law peak at angle zero. See the
[Research reference](DADA_ENGINE_RESEARCH_REFERENCE.md) for the deterministic
maximum method and warm-start compatibility.

Initial guesses may be rescaled compatibly, but integrated mass is never
renormalized between cycles. Only an ordinary complete physical cycle passing
the declared periodic criterion certifies convergence. For hysteretic valves,
the discrete valve state must also recur; matching gas states alone is insufficient.
Continuous ideal diodes have no valve memory to recover. Wall extrapolation or a
safe retry changes an initial guess, not the equations, tolerance or verdict.
A recovered rejected-trial diagnostic stays historical; current constraint
margins describe the final cycle. Sampled screens are not analytic guarantees.

## Scope of ideal-model similarity

For fixed dimensionless motion and the lumped ideal-gas, reservoir-UA/orifice
model, `dada_solver.similarity.scale_capacity_and_speed` scales all volumes and
inventory by `s`, frequency by `f`, and UA and CdA by `s*f`. Corresponding pressure
and temperature histories versus angle are preserved; cycle heat/work scale by
`s`, powers by `s*f`, and COP is unchanged. This mathematical similarity does not
establish how real exchanger conductance, hold-up or hydraulics scale. It is not
the geometry-rebuilding [Research capacity transformation](DADA_ENGINE_RESEARCH_CAPACITY.md).

Within the same ideal model, `scale_volume_at_constant_inventory` trades volume
against pressure while preserving heat, work and COP through matched volume,
CdA and valve-threshold changes. Higher pressure alone is not a COP improvement
or evidence that real hardware retains this similarity.

## Evidence and remaining work

Analytical invariants, mass/energy conservation, refinement and Python/Numba
parity establish numerical properties. They do not calibrate a physical exchanger
or prove a global optimum. Doty whole-bank discrepancies, uniform-wall and
quasi-steady pulse limitations remain explicit. No hidden friction/efficiency
model converts gas work into shaft performance. Distributed exchange, actual
valve losses, mechanics, experimental calibration and real-fluid hydraulics
require separate evidence and decisions.

The application-specific motor demonstrator brief (about 100 W useful output, 2–10 Hz,
25/325 °C sources, 66 L ceiling) describes that application only. Cooling-cell
loads and new Research studies have their own explicit boundaries. See
[validation](validation.md) and [research objectives](MOTOR_RESEARCH_OBJECTIVES.md).
