# Motor demonstrator design brief

Updated 2026-09-09. These are design requirements and screening rules, not a
validated machine or a new exchanger closure. The intended machine is a small
experimental demonstrator, not a commercial product. Commercial refers only to
purchased exchanger candidates. The current reservoir difference is 300 K (25/325 deg C), temporarily raised
from 150 K to guide dimensional optimization. Reduce it progressively after
improving the design; it is not part of the demonstrator name. See
[HIGHER_TEMPERATURE_TRIAL.md](HIGHER_TEMPERATURE_TRIAL.md).

| Quantity | Requirement |
| --- | --- |
| Useful mechanical output | Approximately 100 W; no numerical tolerance specified |
| Cycle frequency | 2–10 Hz |
| Heat input reservoir | 598.15 K (325 degrees Celsius) |
| Heat rejection reservoir | 298.15 K (25 degrees Celsius) |
| Large cylinder maximum enclosed volume | At most 0.066 m^3 (66 litres), including clearance |
| Primary objective | Efficiency subject to power and size requirements |
| Exchanger flow | Intermittent, intended unidirectional flow |
| Component selection | Documented commercial exchanger first; Doty extrapolation fallback |

The cylinder limit is a ceiling, not a design target. It is not a limit on total
machine volume. Report both cylinder volumes, exchanger gas volumes, exchanger
external dimensions and the mechanism envelope where known. Preserve the old
600/300 K, 0.1515 Hz motor example as a historical numerical regression case.

## Power boundary and first screening quantities

The present solver reports indicated gas power, before mechanical and auxiliary
losses. Use `useful_power = indicated_power - mechanical_losses - auxiliary_power`
with an explicit system boundary and without double counting. Unknown losses
remain unknown; a 100 W indicated result cannot certify the useful-output target.

For 100 W useful output, useful work per cycle is 50 J at 2 Hz, 20 J at 5 Hz,
and 10 J at 10 Hz. The reservoir Carnot limit is approximately 50.15 percent.
The following heat duties are illustrative energy balances, not predicted
motor efficiencies. They neglect mechanical and auxiliary losses:

| Assumed efficiency | Heat input for 100 W output | Heat rejection |
| --- | --- | --- |
| 5 percent | 2000 W | 1900 W |
| 10 percent | 1000 W | 900 W |
| 15 percent | 667 W | 567 W |
| 20 percent | 500 W | 400 W |

Actual exchanger sizing also needs simultaneous gas temperatures, pressures,
port flows, pulse durations and allowed pressure loss. Average heat duty alone
cannot determine UA or channel count.

## Commercial evidence status

The [CompAir CWA manufacturer page](https://www.compair.com/en/air-treatment/aftercoolers/)
describes gas/water shell-and-tube aftercoolers with maximum air inlet temperature
200 degrees Celsius and maximum water inlet temperature 90 degrees Celsius.
This makes the family a heat-rejection candidate, but does not qualify it as a
heater supplied by a 175-degree liquid. Both sides require appropriate ratings.

[SWEP's compressed-air application page](https://www.swepgroup.com/applications/manufacturing/compressed-air)
confirms commercial brazed-plate applications for compressed air. An application
page alone does not supply a selected component's gas-side thermal and pressure-
loss maps, internal volume or pulsating-flow validation.

Additional primary documentation reviewed on 2026-09-09:

| Manufacturer | Public evidence | Remaining gaps for this demonstrator |
| --- | --- | --- |
| [Kaeser KAC](https://pr.kaeser.com/en/download.ashx?id=tcm%3A164-37551) | Air-flow capacity tables versus inlet and approach temperatures, dimensions, operating limits; stated pressure drop below 3 psi under catalogue conditions | Full gas pressure-loss map, internal gas volume and pulse response; cooler only |
| [Exergy](https://exergyllc.com/wp-content/uploads/Exergy-Brochure-1.pdf) | Compact shell-and-tube series, heat-transfer areas and pressure ratings; tube-in-tube dimensions and liquid/gas applications | Model-specific gas thermal/hydraulic performance, hold-up and transient behavior |
| [SWEP SSP](https://ssp-wiki.swepgroup.com/Calculations/Singlephase) | Selection and performance calculations with flow, temperature and pressure-drop inputs and dimensional outputs | Actual selected-gas calculation, applicable model limits and pulse validation; an advertised calculation tool is not an acquired experimental map |

Kaeser's smallest listed KAC 50H envelope is approximately 325 x 404 x 348 mm
(12.8 x 15.9 x 13.7 inches), about 46 litres of bounding-box volume. Its air
capacity tables use 80–125 psig operating pressure and its maximum operating
temperature is 350 degrees Fahrenheit (176.7 degrees Celsius). This is useful
stationary gas evidence, but neither a compactness winner nor a heater selection.
A reservoir target is not a bound on instantaneous gas temperature.

Exergy's shell-and-tube series 23 lists a 25 mm shell diameter and transfer areas
0.04–0.21 m^2. This makes it worth investigating for a small prototype, but the
family dimensions and pressure ratings alone do not establish 175-degree gas
performance or suitability. No manufacturer has been contacted.

The practical approach is to fit the demonstrator to an available component's
stationary data, use explicitly qualified extrapolation where needed, and
measure prototype heat transfer and pressure loss. Complete pulsating-flow
manufacturer data are desirable, not a prerequisite for all further progress.

No component has yet been qualified against this motor brief. Qualification
requires a specific model and channel configuration, dimensions, gas hold-up,
fluid and temperature/pressure limits on both sides, and paired heat-transfer
and pressure-drop data with their reference conditions. Do not substitute a
liquid/liquid rating for gas performance.

## Laminar and turbulent operation

For a passage, `Re(t) = rho(t) * u(t) * D_h / mu(t)` or equivalently
`Re(t) = mass_flow(t) * D_h / (mu(t) * flow_area)`. Use the actual passage area,
not an orifice CdA. A single exchanger can span regimes during one pulse.
At fixed mass flow and viscosity, density cancels from this expression;
increasing pressure is not an independent Reynolds-number lever. At fixed
volumetric flow, density matters. Temperature also changes transport properties.

The Reynolds definition follows [NASA](https://www.grc.nasa.gov/WWW/K-12/airplane/reynolds.html).
[SWEP](https://www.swepgroup.com/product/brazed-plate-heat-exchanger/how-does-it-work)
also describes how plate patterns change thermal and hydraulic characteristics.
Do not impose smooth-pipe transition thresholds on corrugated plates. Compare
cycle efficiency, pressure loss and size for each documented geometry. Laminar
parallel microtubes can offer substantial transfer area; turbulence is not an
objective in itself. Cycle frequency alone does not identify the flow regime.

## Authorized Doty fallback

Use the transcribed nitrogen measurements in `examples/data/` as experimental
anchors, preserving their original fluid, pressure, temperatures and uncertainty.
Before predicting a motor component, document prototype geometry and distinguish
replicating parallel passages from changing passage dimensions. Account for
manifold losses and flow distribution. Extrapolated UA and pressure loss must
be identified as estimates, with sensitivity scenarios rather than an invented
statistical confidence interval.

A quasi-steady map evaluated on simultaneous cycle states is the first screening
step. It assumes that steady measurements represent the instantaneous pulse;
that assumption is unvalidated at 2–10 Hz. It introduces no distributed thermal
state. A dynamic wall/gas model remains a separate extension whose benefit and
computational cost must be explained before implementation.

## Selected Doty calculation reference

The user selected the Doty fallback as the current working basis. The first
whole-bank scaling implementation and its limitations are documented in
[DOTY_SCREENING.md](DOTY_SCREENING.md). This adds a separate steady screening
calculation, not a change to the conservative cycle equations.

The four-bar integration and matched speed comparison are documented in
[MOTOR_FOUR_BAR.md](MOTOR_FOUR_BAR.md).

Ordered long-term motion-law and design-map objectives are recorded in
[MOTOR_RESEARCH_OBJECTIVES.md](MOTOR_RESEARCH_OBJECTIVES.md).

The selected external-air boundary and geometry sizing foundation are described
in [AIR_SOURCE_EXCHANGERS.md](AIR_SOURCE_EXCHANGERS.md).

The opt-in air/wall coupling, numerical results and continuation checkpoint
are in [AIR_WALL_COUPLING.md](AIR_WALL_COUPLING.md).

The geometry-to-volume, resistance, wall-capacity, hydraulic and fan connection
is now implemented; see [GEOMETRIC_MOTOR_COUPLING.md](GEOMETRIC_MOTOR_COUPLING.md).

Current doubled-exchanger and independent lower-Lambda trials, with external
aerodynamic losses excluded, are in [DOUBLED_EXCHANGER_TRIAL.md](DOUBLED_EXCHANGER_TRIAL.md).
