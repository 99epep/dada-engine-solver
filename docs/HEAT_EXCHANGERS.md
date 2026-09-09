# Heat-exchanger screening and sizing

The exchanger module translates a real gas-passage geometry into quantities
that the thermodynamic sizing problem can constrain. It is a low-order design
tool, not a CFD solver and not an experimental calibration.

Current displays use `H_i` for the heat-in exchanger and `H_o` for the heat-out
exchanger. The historical `cold`/`C` and `hot`/`H` branch identifiers below are
retained in Python and TOML for compatibility; their external reservoirs swap
in motor operation. See [MOTOR_OPERATION.md](MOTOR_OPERATION.md).

## Model boundary

The implemented architecture is a finite set of identical, straight,
rectangular channels in parallel. The model calculates:

- internal gas volume and gas-side area;
- total open flow area and hydraulic diameter;
- mean density, velocity, Reynolds number and Mach number;
- distributed and user-supplied minor pressure losses;
- gas-side heat-transfer coefficient;
- overall `UA` from gas film, wall, external-side and additional resistances;
- reservoir effectiveness `1 - exp(-UA / (m_dot Cp))`;
- feasibility against explicit `UA`, pressure-drop, Mach, gas-volume and
  effectiveness constraints.

All inputs and outputs use SI units. Gas viscosity and thermal conductivity
are explicit inputs at the declared representative state. Wall properties,
external `h*A`, additional contact/fouling resistance, surface roughness and
minor-loss coefficient are also explicit. The program supplies no material,
fluid, manifold or fabrication value silently.

The cold and hot exchangers may use the same internal architecture while
having different external models. A water-glycol loop and ambient air do not
generally have the same external resistance. The present API expresses that
difference through `external_conductance`; future reservoir-side models can
replace it without changing the gas-passage geometry.

## Correlations and validity

For fully developed laminar flow, the implementation uses the classical
rectangular-duct Poiseuille-number and constant-wall-temperature Nusselt-number
polynomials. For turbulent flow it uses the Haaland Darcy friction factor and
the Gnielinski Nusselt correlation. The interval `2300 <= Re < 4000` is not
interpolated: thermal conductance and pressure drop are reported as
`unavailable`. Values outside the declared turbulent Reynolds and Prandtl
ranges are treated the same way.

The entrance-length check is deliberately conservative. The fully developed
correlations are unavailable when the channel is shorter than the estimated
hydrodynamic or thermal entrance length. A future entrance-region model may
recover such designs, but the current model does not extrapolate them.

Pressure drop is calculated with a mean-density Darcy model. It is not a
compressible duct solution. Every sizing study must therefore set both an
absolute pressure-drop limit and a maximum `delta_P / P` limit. Choking,
Fanno-flow effects and pressure-dependent density along a passage are not
resolved here.

The scientific basis and limitations should be reviewed against the original
literature, including:

- V. Gnielinski, *New equations for heat and mass transfer in turbulent pipe
  and channel flow*, International Chemical Engineering 16 (1976), 359-368.
- S. E. Haaland, *Simple and explicit formulas for the friction factor in
  turbulent pipe flow*, Journal of Fluids Engineering 105 (1983), 89-90,
  DOI `10.1115/1.3240948`.
- R. K. Shah and A. L. London, *Laminar Flow Forced Convection in Ducts*,
  Academic Press (1978), for fully developed non-circular duct correlations.
- [Wire-mesh oscillatory-flow experiments and correlations](https://www.sciencedirect.com/science/article/pii/S1110016815000927),
  which are relevant to a future porous-matrix implementation but are not used
  by the rectangular-channel model.
- [Experimental regenerator wire-mesh flow and heat-transfer study](https://cir.nii.ac.jp/crid/1390282679648971648),
  documenting why steady correlations cannot by themselves validate an
  oscillatory porous exchanger.
- [Experimental oscillating-air pin-fin measurements](https://doi.org/10.1155/2013/283830),
  which found a 20--34 percent enhancement without bypass but degradation for
  excessive bypass in its tested apparatus. This is evidence that manifold and
  bypass geometry matter, not a correction factor for DADA.
- [Laminar pulsating-flow analysis in a rectangular channel](https://doi.org/10.1016/j.ijheatmasstransfer.2018.08.109),
  which found a reduction of time-averaged Nusselt number for its boundary
  conditions. Together with the preceding experiment, this rules out assuming
  that pulsation is universally beneficial.
- [Oscillatory shell-and-tube air/water experiments](https://www.osti.gov/servlets/purl/1026487),
  which organize effectiveness using peak-flow Reynolds number and thermal
  penetration depth. These are appropriate future similarity variables for
  DADA exchanger validation.

The literature therefore supports the present steady-flow model as a screening
baseline, but it does not provide one universal multiplier that converts it
into a validated DADA exchanger. Published oscillatory effects can have either
sign and depend on waveform, frequency, penetration depth, bypass and geometry.
The example transport values are ordinary-air order-of-magnitude inputs; the
configured `external_conductance` remains an assumed boundary value until the
complete glycol-side or ambient-side geometry is specified.

## Moisture boundary

Ambient air introduces water into the nominally single-phase working gas. An
optional `[humidity]` section supplies the initial relative humidity as a
fraction from zero to one. The solver preserves the corresponding initial
water mole fraction only for a screening calculation and compares its local
partial pressure with equilibrium saturation pressure throughout S, L, C and
H.

Water and ice saturation pressures use the Murphy--Koop (2005) equations:
[review and parametrizations](https://doi.org/10.1256/qj.04.94). The definition
of dew/frost point is consistent with the
[NIST humid-air reference](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=921756).
Below the water triple point, equilibrium saturation over ice is used.

The possible outcomes are `clear`, `condensation_or_frost_risk`, and
`unavailable`. A risk result marks the dry, calorically-perfect-gas trajectory
as crossing a model boundary. It does not estimate condensed mass, latent
heat, drainage, corrosion, valve adhesion or passage blockage. Those effects
must not be added until water becomes an explicit transported constituent.
For each control volume, the report gives the maximum saturation ratio and the
first saturated sample relative to `theta = 0`, including angle, temperature,
pressure and predicted equilibrium phase. When plotting is enabled, an eighth
panel shows the complete S/L/C/H saturation-ratio histories and the unit
saturation boundary.

## Coupling to the cycle

The thermodynamic model represents each exchanger as a lumped control volume.
Its two instantaneous port flows may differ because the volume accumulates
mass. Consequently, the zero-dimensional state does not define a unique
internal duct velocity.

`conservative_port_flow_operating_point` bridges that gap explicitly for early
screening: it selects the largest absolute flow observed at either adjacent
port and combines it with caller-supplied representative pressure and
temperature. This is conservative for flow capacity, but it is not a derived
pulsatile velocity history.

The geometric result may provide candidate values for exchanger gas volume and
`UA` in a subsequent thermodynamic run. It must not silently replace the four
configured hydraulic `CdA` values. A geometry-based compressible passage law
requires a later coupled implementation and validation.

## Command line

The illustrative configuration
[`examples/exchanger_screening_example.toml`](../examples/exchanger_screening_example.toml)
defines a finite design grid:

```console
dada-exchanger examples/exchanger_screening_example.toml
```

The program reports all feasible geometries or returns a nonzero status when
the supplied grid contains none. Inputs in that file are examples, not
validated DADA exchanger properties.

## Constrained geometry proposal

The grid is also the deterministic starting set for a geometry optimizer:

```console
dada-exchanger examples/exchanger_screening_example.toml \
    --optimize minimum_gas_volume
```

Available objectives are:

- `minimum_gas_volume`;
- `minimum_pressure_drop`;
- `minimum_gas_side_area`.

The `UA`, pressure-drop, pressure-drop-fraction, Mach, gas-volume,
effectiveness and correlation-availability requirements remain hard
constraints. They are not combined into a weighted score.

Channel count remains an integer and is enumerated from the configured list.
For each count with at least one feasible grid point, channel width, height and
length are refined continuously in logarithmic coordinates with SLSQP. The
minimum and maximum grid values are the explicit manufacturing bounds. A
configuration intended for optimization must therefore contain at least two
distinct values for each continuous dimension.

Only feasible grid points seed local refinement. This avoids fabricating a
smooth numerical objective across transitional or otherwise unavailable
correlation regimes. It also means that a coarse grid can miss a narrow
feasible region. The result is the best design found by deterministic
screening and local refinement, not a proof of global optimality.

The API also returns all non-dominated feasible candidates in gas volume,
pressure drop and achieved `UA`. This Pareto set allows the machine-level
sizing process to reconsider a slightly larger exchanger when it offers a
large hydraulic benefit.

For the illustrative file, minimum-volume optimization reduces the best
screened gas volume from `1.60e-4 m3` to approximately `1.17e-4 m3` while
meeting its artificial `20 W/K` target. This verifies the optimization path;
it is not a proposed DADA component.

## Cycle/exchanger fixed point

The coupled command alternates a complete periodic simulation and two
independent constrained geometry searches:

```console
dada-exchanger-coupled examples/exchanger_coupled_example.toml
```

The coupled file references one thermodynamic configuration plus separate
cold- and hot-exchanger files. This permits different external conductances,
materials, channel bounds, constraints and objectives on the glycol and heat-
rejection sides.

At every outer iteration the adapter obtains a conservative port-flow value
for each exchanger. It uses the minimum pressure and maximum temperature of
the corresponding lumped control volume as a conservative density point for
flow-capacity screening. The flow passed to successive searches can be under-
relaxed; the thermodynamic state itself is never blended.

The selected geometries supply new gas volumes and predicted overall `UA`
values. These four values are inserted into a new simulation configuration,
and the periodic cycle is solved again. Convergence requires the maximum
relative change among both screening flows, both volumes and both conductances
to satisfy the configured tolerance. One additional periodic simulation
verifies the final updated configuration before convergence is reported.

Each successful periodic cycle is retained as the initial guess for the next
outer iteration. When a changed exchanger volume changes the inventory implied
by configured charging pressure, all control-volume masses and energies are
scaled by the same factor. This preserves mass fractions, temperatures and
specific energies while matching the new required total inventory exactly.

The transition from one hydraulic closure to the next uses a configurable
numerical continuation. Intermediate solves linearly interpolate exchanger
volume and `UA` and blend the two predicted mass-flow rates at identical node
states. These intermediate closures are not physical results and are never
reported as candidate performance. The final continuation fraction is exactly
one, at which point the configured geometric closure is used without blending.
Failure at any intermediate fraction is reported explicitly rather than being
hidden by relaxation.
When two consecutive optimizations return exactly the same hydraulic network,
the redundant intermediate fractions are skipped and only the final
configuration is solved. This does not weaken continuation when a hydraulic
closure actually changes.

Possible terminal states distinguish convergence, maximum outer iterations,
periodic-cycle failure, cold-exchanger infeasibility and hot-exchanger
infeasibility. The report includes the complete iteration history, both final
geometries and final-cycle performance when available.

## Geometric hydraulic coupling

The coupled loop can now replace the four reference orifice closures with a
geometric quasi-steady model. The reference `CompressibleOrifice` remains the
default for ordinary simulation files and can always be selected for direct
comparison.

Each exchanger core is divided explicitly between its inlet and outlet ports.
The configured fractions must add implicitly to unity: the outlet receives
`1 - inlet_core_resistance_fraction`. Inlet and outlet collector-loss
coefficients are separate required coupled-configuration values. On the
nominal valve outlet, the distributed duct resistance is placed in series with
the independently configured valve `CdA`; the valve remains passive and
one-way.

The geometric closure solves a mean-density Darcy relation for mass flow. In
the laminar regime the rectangular-duct relation reduces to a linear plus
quadratic pressure-loss equation and is solved analytically. Turbulent flow is
solved iteratively. An isentropic mass-flow ceiling prevents the low-Mach
closure from predicting more than the sonic capacity of its smallest local
area. The series conversion of valve `CdA` to a local loss coefficient is a
low-Mach approximation, not a Fanno solution.

Before a geometry enters the periodic simulation, both ports must pass the
required screening mass flow within the smaller of the configured absolute
and relative pressure-drop limits. Failure is reported as exchanger
infeasibility rather than compensated by changing the target `UA`.

The approximation follows the conservative limits documented by the
[NIST comparison of incompressible and isothermal compressible pipe-flow
formulae](https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote356.pdf):
compressibility cannot be neglected freely as Mach approaches one third or
relative pressure loss approaches ten percent. The code therefore reports
Mach and `delta_P/P` independently and does not label results beyond the
configured limits valid. Sonic capacity follows the official
[NASA isentropic-flow relations](https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html).

## Hydraulic-inertia decision

`assess_hydraulic_quasi_steady_validity` calculates four configurable
indicators from a periodic mass-flow history:

- maximum Mach number;
- maximum pressure-drop fraction;
- acoustic transit time divided by cycle period;
- estimated inertial pressure divided by reference pressure, using
  `delta_P_inertia = (L/A) d(m_dot)/dt`.

No acceptance threshold is hard-coded. If either time/inertia threshold is
violated, `inertia_model_recommended` becomes true. Only then is a momentum
state for each passage justified as the next physical level.

`compare_hydraulic_models` evaluates reference `CdA` and geometric closures at
the same pressure ratios and temperature without fitting either model. The
coupled report also retains the initial reference-cycle performance for a
before/after comparison when a geometric iteration has actually run.

The exploratory cooling-cell files are intentionally separate from the
generic demonstration:

- `examples/cooling_cell_exchanger_exploratory.toml`;
- `examples/cooling_cell_exchanger_coupled_exploratory.toml`.

Their air transport properties, wall model and `500 W/K` external conductance
are provisional inputs, not measurements. They are useful for sensitivity and
software verification only.

The 20-minute similarity point requires about `382.9 W/K` per exchanger after
its speed is adjusted to 27.85 rpm. Directly increasing the gas volume to
6.4 litres reduced cooling to about 41.5 W and COP to about one; exchanger gas
volume cannot therefore be traded freely for area. With 0.4 litre per exchanger,
the same first-level cycle instead produces 347.58 W cooling, requires 103.83 W
thermodynamic input and retains COP 3.3476 at a 1 bar filling pressure.

The former lumped external conductance is no longer the only available closure.
`two_sided.py` evaluates an explicitly configured incompressible liquid channel
side, combines gas convection, wall conduction and liquid convection in series,
and reports liquid pressure drop and ideal pump input separately. Liquid density,
heat capacity, viscosity, conductivity, flow and pump efficiency are mandatory;
the solver supplies no invented glycol properties or pump efficiency. The gas-
side geometric pressure loss remains a screening value and must not be added to
the cycle work a second time when it has already been represented by fitted CdA.

## Compact plate-fin candidate for the 20-minute cell

The compact-exchanger literature reports surface densities above
`10000 m2/m3` for passages below 1 mm. Brazed plate-fin construction is an
established gas-to-liquid architecture, with alternating passages and fins:

- [Review of compact and microchannel air-side correlations](https://doi.org/10.1016/j.enconman.2018.06.104)
- [Alfa Laval brazed plate-and-fin architecture](https://www.alfalaval.com/products/heat-transfer/plate-heat-exchangers/plate-and-fin-heat-exchangers/brazed-plate-and-fin-heat-exchangers/)

A straight-fin screening candidate uses 200 parallel gas passages, each 200 mm
wide, 0.30 mm high and 11 mm long. Its calculated core quantities at 0.096 kg/s,
72 kPa and 264 K are:

| Quantity | Screening value |
|---|---:|
| Gas volume | 0.132 L |
| Gas-side area | 0.881 m2 |
| Gas velocity | 8.0 m/s |
| Mach number | 0.026 |
| Core pressure drop including assumed K=1 | 244 Pa |

Alternating 0.75 mm liquid passages provide approximately the same shared area.
Using deliberately conservative 40% propylene-glycol properties tabulated at
-20 degrees Celsius, a 5 K loop rise, 0.1 mm aluminium separating walls and an
exploratory pump efficiency of 30%, the model returns about 227.5 W/K and less
than 1 mW of ideal core pumping power. The fluid data are from the
[DOWFROST HD technical data sheet](https://www.dow.com/content/dam/dcc/documents/en-us/productdatasheet/180/180-01315-01-dowfrost-hd-tds.pdf).

This result establishes plausibility, not hardware validity. The 0.4 L control-
volume allowance leaves about 0.268 L for gas headers, but header distribution
and pressure loss are not yet demonstrated. The 0.30 mm gas gap is vulnerable
to retained condensate, frost, manufacturing variation and fouling. Oscillating
intermittent flow may also change heat transfer and pressure loss relative to
the steady fully-developed correlations. A 0.5 mm option is less vulnerable but
uses roughly 0.35 L of core gas volume for similar UA, leaving almost no header
allowance. Header design and minimum clear passage are therefore the next Pareto
variables; neither candidate is selected yet.

## Required next validation

Before selecting hardware, each promising geometry still requires:

1. a manifold layout and a distribution analysis;
2. reservoir-side geometry and convection data;
3. a pressure-drop test over the relevant bidirectional flow range;
4. a heat-transfer test over the relevant temperatures and mass flows;
5. a pulsatile test at the intended rotational frequency;
6. fitting or replacing the low-order closures using those measurements.

CFD is justified only if collector maldistribution, entrance flow or complex
porous geometry cannot be bounded adequately by correlations and bench tests.
