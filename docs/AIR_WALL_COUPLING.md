# Conservative air/wall motor coupling

## Implemented opt-in model

`exchangers.air_wall.AirWallMotor` couples the existing four gas control volumes
to two independent wall energies. It reuses the motor geometry and hydraulic
balances. It requires continuous ideal diodes and rejects discrete hysteretic
valves rather than discarding valve memory. The existing configuration loader,
eight-state integrator, reporting and historical results are unchanged.

For each exchanger, let C_air = mass_flow_air * cp_air. A quasi-steady external
stream against one uniform wall has conductance
`G_external = C_air * (1 - exp(-G_air_wall/C_air))`.
Then `Q_air = G_external*(T_air_in-T_wall)`,
`Q_gas = G_gas_wall*(T_wall-T_gas)`, and `dE_wall/dt = Q_air-Q_gas`.
Wall temperature is `E_wall/C_wall`. The gas heat term remains active at zero
working-gas flow. At zero external airflow, forced external heat is zero and
outlet temperature is unavailable; natural convection is not implicitly added.

Total stored energy is the sum of gas internal energies and wall energies.
Its derivative equals external-air heat input minus gas boundary work. The
integrator also tracks gas-side heats separately so wall storage is not mistaken
for external input. Angular integration divides all rates by positive cycle
progress speed. Piecewise-law breakpoints split integration intervals.

Only the two wall energies are new physical states. Five quadratures track heat
and work. This small state count does not guarantee a small runtime penalty:
thermal stiffness and convergence of the walls must be measured. No wall
conductance is fitted from arbitrary mechanical losses or target efficiency.

## Reproducible screening

Run `PYTHONPATH=src python3 examples/motor_air_wall_screening.py`.
The example uses the 1 litre, 2 Hz piecewise-linear reference, 100 W/K on each
side of each wall, 10 J/K per wall, external air at 0.05 kg/s and constant
cp = 1005 J/(kg K). These are explicit numerical assumptions, not selected
hardware or calibrated Doty parameters. The script checks periodicity of all
eight gas states and both walls and reports unavailable fan/useful power.
The maximum cycle count is 50; nonconvergence must remain visible.

The standalone wrapper does not yet export valve event chronology through the
historical event reporting interface. Existing constant-reservoir reports must
not be applied to these trajectories, since they would recompute wrong heats.

## Doty calibration audit and limitation

The three nitrogen rows in the existing CSV have hot-stream temperature drops
about 5.1–5.3 percent larger than cold-stream rises. With the experiment's equal
mass flow and the simplifying equal constant cp assumption, these imply a heat
imbalance of that size; see `outputs/doty_energy_balance_audit.csv`. This is not
proof of a particular loss mechanism. Temperature uncertainty, property changes,
and heat leakage require consideration before fitting both outlet temperatures.

A single uniform wall and a well-mixed working-gas volume cannot in general
reproduce a high-effectiveness counterflow exchanger's axial temperature
profiles. A steady gas/gas UA also does not identify its two resistance shares
or wall capacity. Accordingly, this coupling is an energy-conserving dynamic
screening model, NOT a reproduction or validation of Doty's measurements.

Next calibration work must retain a counterflow steady reference, preserve
reported uncertainties and test an axial discretization if the single-wall
approximation materially affects the ranking. Tube/header hydraulics, finite
external pressure loss and fan efficiency must then be coupled consistently.
No 100 W useful-output design or complete exchanger sizer is established here.

## Numerical results and continuation checkpoint

The initial run used scalar absolute tolerance 1e-9 for both masses and energies.
It reached the 50-cycle limit with scaled periodic error 14.75; that result is
preserved as `motor_air_wall_screening_loose_tolerance.json`. It is not a
converged reference. Reducing wall capacity from the earliest exploratory
100 J/K to 10 J/K was a change of numerical scenario, not a fitted material
property; only the 10 J/K scenario is evaluated here.

The corrected 10 J/K run uses rtol 1e-8 and separate absolute tolerances:
1e-13 kg for masses, 1e-8 J for gas energies, 1e-7 J for wall energies and
1e-8 J for quadratures. All ten states reach the periodic threshold in six
cycles. It gives approximately 23.637 W indicated power, 146.006 W external
heat input and 16.189 percent indicated thermal efficiency. Cycle energy
balance including wall storage closes to about 3.4e-12 J. These are numerical
screening results with assumed conductances, wall capacity and orifice losses;
fan power and useful power remain unavailable. They do not validate the 100 W
demonstrator or replace the four-bar comparison.

The Doty audit additionally applies the paper's apparent-UA reduction using
an explicit approximate nitrogen cp of 1040 J/(kg K). It obtains 3.681, 5.242
and 5.397 W/K against 3.7, 5.3 and 5.4 W/K reported, within the tabulated UA
uncertainties. This is a data-reduction consistency check, not independent
prediction of the measurements. Do not infer that the one-wall dynamic model
reproduces the counterflow experiment.

Remaining work, in order:

1. Calibrate a steady counterflow reference with gas-property and heat-imbalance
   assumptions explicit; compare lumped and spatial thermal descriptions.
2. Link tube, collector and wall geometry to actual gas hold-up, heat capacity,
   calibrated hydraulic loss and thermal resistances; add external fan demand.
3. Recalculate coupled four-bar cases and size hardware toward approximately
   100 W useful output at 2–10 Hz, with the 66 litre ceiling and efficiency priority.
4. Later search independent motion laws and loop geometry, then map Lambdas and
   swept-volume ratio as recorded in MOTOR_RESEARCH_OBJECTIVES.md.

All project text remains English. External source inlet temperatures are
298.15/448.15 K. Preserve intermittent intended forward flow and any computed
local reflux. Neither mirror symmetry nor the piecewise-linear reference is
an imposed optimum. Isothermality remains diagnostic-only.

A refinement cycle starting from the converged state uses tenfold tighter
relative/absolute tolerances and half the maximum angular step. Indicated power
changes by approximately 5.6e-7 relative (0.000056 percent). This is a local
one-cycle numerical check, not a fresh fully converged refined periodic search.
See `outputs/motor_air_wall_refinement.json`.

The geometry-to-volume, resistance, wall-capacity, hydraulic and fan connection
is now implemented; see [GEOMETRIC_MOTOR_COUPLING.md](GEOMETRIC_MOTOR_COUPLING.md).
