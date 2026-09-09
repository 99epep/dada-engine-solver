# Geometry-connected motor exchanger screening

## Implemented connection

`exchangers.hardware.connect_hardware` builds an `AirWallMotor` from independent
heat-input and heat-output tube banks and property sets. The resulting model
uses the actual parametric tube/header gas volumes, tube metal heat capacities,
series thermal resistances and four hydraulic half-passages. The existing motor
kinematics and passive valve rules are retained. The old constant-UA/CdA
exchanger links are replaced, not added on top of the new passage losses.

The thermal network contains inner gas film, cylindrical metal wall and outer
air film. Film conductances are `Nu*k*A/D`; Nusselt numbers and constant
properties are explicit screening inputs. The radial metal resistance is
`ln(d_outer/d_inner)/(2*pi*k_metal*N*L)`, split equally in resistance around
the lumped wall state. Wall capacity is metal volume times density and specific
heat plus explicit extra capacity for tube sheets/collectors. This is a lumped
radial representation, not a resolved transient metal temperature profile.
See [NPTEL cylindrical resistance notes](https://archive.nptel.ac.in/content/storage2/courses/112108149/pdf/M2/Student_Slides_M2.pdf).

External air is modeled as longitudinal flow through the free bundle area.
Its equivalent hydraulic diameter uses the tube outer perimeters and enclosure
perimeter. A caller-supplied Poiseuille number, viscosity and minor-loss K give
an estimated pressure drop. This equivalent laminar passage is not a validated
bundle correlation. Fan electrical power is `pressure_drop * volume_flow /
total_efficiency`; total efficiency includes motor and drive. See the
[DOE fan sourcebook](https://www1.eere.energy.gov/manufacturing/tech_assistance/pdfs/fan_sourcebook.pdf).
A real fan curve, ducts and heater resistance remain to be specified. No fan
heat is silently credited to the declared inlet temperature.

## Hydraulic allocation and limits

Half the tube friction and half the specified total header K are placed on
each side of the exchanger gas-storage node. At equal port flow, their sum is
one full exchanger loss. When gas is stored, the two instantaneous port flows
can differ. Header K is defined relative to total tube-passage velocity, not an
unspecified collector inlet velocity. Additional connection volume remains an
explicit input.

The mean-density laminar law is inverted analytically using a stable quadratic
root. The outlet half-passage includes the specified valve CdA as an additional
quadratic resistance, and an orifice-based flow cap. This is an approximate
series closure, not resolved compressible valve lift or Fanno flow. Reverse
flow is allowed on internal links and blocked by the existing passive valves.
No second full-bank loss is added at the inlet or outlet.

The core loss multiplier is explicit. One means theoretical tube friction; it
is not an empirical fit to Doty's pressure-drop measurements. Fitting a
multiplier to a whole measured bank may absorb its original header losses and
must not then add those same original losses twice. New collector losses and
calibration residuals require separate uncertainty scenarios.

## Example and result boundary

Run `PYTHONPATH=src python3 examples/motor_hardware_screening.py`.
Inputs are in `examples/motor_hardware.toml`, with independent air conditions for
the hot and cold sides. Tube diameter, wall thickness and pitch originate from
Doty's prototype; count, length, square packing and headers describe a redesign.
Illustrative metal properties, Nusselt numbers, external airflow, fan efficiency
and header K are stated explicitly, not presented as manufacturer data.

The runner uses the 1 litre four-bar seed at 2 Hz. Inventory is held equal to
that seed while gas distribution and initial fill pressure adapt to the new
exchanger volume. It checks all ten states for periodicity. Geometry, assumed
properties, resistance components, fan estimates, Reynolds/Mach diagnostics
and convergence history accompany the performance results in JSON.

Indicated power minus fan electricity is reported separately from useful shaft
power. Mechanical losses and any conversion needed to supply fan electricity
are unknown, so useful shaft power remains unavailable. Internal hydraulic
losses already affect the cycle and must not be subtracted again as a separate
pump-power penalty. Negative results must remain negative. A laminar Reynolds
screen does not validate entrance effects, pulsation, constant Nusselt numbers
or uniform flow distribution; all results remain screening estimates.

The geometry-to-cycle chain is implemented. A calibrated Doty counterflow
thermal reference and an optimization over hardware candidates remain further
work; this example is not an optimized or validated 100 W demonstrator.

## Numerical continuation

Wall capacities are derived from the candidate geometry and are not reduced
just to accelerate convergence. The runner allows up to 250 cycles and writes
`outputs/motor_hardware_checkpoint.json` after every cycle. A checkpoint is
an initial guess only, never a converged result. Resume with
`--restart outputs/motor_hardware_checkpoint.json`; hardware inputs must match.
The final report includes periodic status, mass drift and total gas-plus-wall
energy residual. Port-flow CSV and duty summaries retain signed local reflux.

For the first input geometry, each exchanger contains approximately 13.23 ml
of working gas and has wall heat capacity 38.25 J/K. The estimated static
conductance is 3.046 W/K. The gas-film, metal and air-film resistances are about
0.07387, 0.0001769 and 0.25428 K/W respectively. These are consequences of the
stated constant-property/Nusselt inputs, not measured component ratings.
Estimated fan electricity totals about 0.496 W for both streams, excluding
unmodeled duct/heater losses and any fan sizing mismatch.

## First complete coupled result

The four-bar case reaches the ten-state periodic threshold after 196 cycles
starting from the uniform filling state and source-temperature walls. The final
indicated gas power is -39.031 W: this geometry consumes mechanical work at 2 Hz.
Subtracting the estimated 0.496 W of fan electricity gives -39.527 W before
unknown mechanical/conversion losses. No thermal motor efficiency or useful
shaft output is reported for this failed motor candidate.

Tube Reynolds number peaks near 1595 under the stated constant-viscosity
assumption. The representative exchanger-state Mach estimate peaks at 0.219,
above the configured 0.2 screening limit; this is explicitly a failed Mach
screen. Low Reynolds number alone does not make the model valid. External-air
and pulse/entrance correlations also remain unvalidated.

Gas-plus-wall cycle energy residual is approximately -1.9e-11 J; total gas-mass
drift is approximately -1.7e-17 kg. A one-cycle refinement with tenfold tighter
tolerances and half the maximum angular step changes indicated power by about
5.4e-6 W. This checks local numerical sensitivity, not physical validity or an
independently converged refined periodic state.

The result demonstrates the connection, not an optimized design. The seed's
weak external thermal conductance must be addressed jointly with header volume,
internal pressure loss and fan demand. Keep the requested 2–10 Hz range and
100 W useful-output target; do not rescue this candidate by silently returning
to the old slow-frequency case or assigning an unrelated UA.
