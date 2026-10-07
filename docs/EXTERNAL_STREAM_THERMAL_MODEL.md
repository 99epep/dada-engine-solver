# External streams and conservative wall storage

## Physical boundary

The chain is working gas → internal gas film → lumped tube wall and thermal
storage → declared external conductance → finite-capacity external fluid stream.
Arrows indicate thermal coupling, not mandatory heat direction. Passive valves
remain a separate hydraulic subsystem; their equations and placement do not
change with the external fluid label.

`exchangers/external_stream.py` provides `ExternalFluidStream`,
`ExternalStreamWallExchanger`, and `ExternalStreamWallMachine`.
A stream specifies inlet temperature [K], mass flow [kg/s], Cp [J/(kg K)],
wall-to-stream conductance [W/K], and a nonempty fluid/scenario label. A label
such as water/glycol or oil selects **no** empirical correlation. Cp and
conductance are declared scenario inputs, with their own required scientific
justification in a real study.

For stream capacity rate C = mass_flow × Cp and conductance G:

- K = C × (1 − exp(−G/C)), evaluated with `expm1`;
- Q_external = K × (T_inlet − T_wall), positive into the machine;
- T_outlet = T_inlet − Q_external/C;
- Q_gas = G_internal × (T_wall − T_gas);
- dE_wall/dt = Q_external − Q_gas; T_wall = E_wall/C_wall.

Zero flow gives zero external heat and an unavailable outlet temperature. Zero
conductance also gives zero heat. Neither limit discards gas/wall heat exchange
or wall storage during piston pauses. No changes are made to the periodic
convergence test, conservative gas states, or work quadrature.

`ExternalStreamMicrotubeExchanger` retains the production tube bank, finite wall
capacity, internal resistance and optional variable gas film, and the same
`TubeHalfLink` hydraulics. Declared external G is the **total wall-node-to-stream
conductance**, including any external half-wall resistance. It is not an overall
gas-to-stream UA. Geometry/material changes update internal resistance and wall
capacity. They do not automatically invent a new external G correlation.

External pressure drop, pump/fan consumption, liquid-side Nusselt correlations,
fluid phase change, and liquid-loop dynamics are **not modeled** by this family.
Exclusion from the balance does not mean zero physical loss. The air hardware
family computes optional fan estimates; these are not liquid-loop estimates.

## Operating direction and signed performance

The ten stored states are eight gas mass/internal-energy values plus two wall
energies. Five quadratures track external heats, internal gas heats and gas
work. Integration advances positive crank progress; the established kinematics
adapter applies motor direction once. The wrapper accepts either operating
mode with continuous ideal diodes. Hysteretic valve memory is still unsupported
by the wall wrapper; its state dimension has not been silently extended.

Branches remain H_i (heat_in) and H_o (heat_out). In refrigeration, H_i contacts
the cold stream and H_o the hot stream. In a motor the assignment reverses.
`performance.classify_cycle` requires all three signs:

| Mode | Study angular speed | Q_i | Q_o | Gas work W |
|---|---:|---:|---:|---:|
| Refrigeration | positive | positive | negative | negative |
| Motor | negative | positive | negative | positive |

Consuming work alone does not establish refrigeration; invalid heat signs leave
COP unavailable. Cooling power = frequency × Q_i; heating power = −frequency ×
Q_o; indicated input power = −frequency × W. COP_c = Q_i/(−W) and COP_h =
−Q_o/(−W), only for refrigeration. Motor efficiency remains W/Q_i.
The energy residual is Δ(U_gas + E_walls) − Q_i − Q_o + W. The retained storage
term matters before perfect periodicity; COP_h − COP_c approaches 1 as the
periodic storage change vanishes. This is indicated work, not measured human
shaft input. Mechanical losses and human mechanical efficiency remain unknown.

## Research configuration

Use the shell helper in the [Research guide](DADA_ENGINE_RESEARCH.md).

```sh
research init external-stream-refrigeration --output outputs/my_cooling/study.toml
research validate outputs/my_cooling/study.toml
research evaluate outputs/my_cooling/study.toml --output outputs/my_cooling/reference.json --budget 2m
research report outputs/my_cooling/reference.json --html outputs/my_cooling/report.html
```

`external-stream-motor` is also available. Both presets accept `--small` and
`--large` families, including `structured_c2_15p`. The preset is a small thermal
boundary validation fixture, not a human-cell design. The harmonic refrigerator
uses 0.2 Hz, 30% clearance ratios, 278.15/298.15 K streams, declared Cp 4180
J/(kg K), 0.05 kg/s and 150 W/K on each side. The increased clearances keep the
startup compression/expansion inside the existing gas-transport domain. Current species-dependent domains are documented in the
[microtube reference](MICROTUBE_GAS_MODEL.md); this fixture is not a domain definition.

Schema 3 defines parameter ownership, adapters, candidates, history, exact caches,
Sobol and reporting. Its portable `.basis.json` defines exchangers with
`family = "external_stream_wall"` and `inputs.external_stream`. Operation is
selected by the existing signed `configuration.angular_speed`; an active
`operation.frequency_hz` is a positive magnitude and does not change direction.
The existing `configuration.gas` either contains the unchanged ideal-gas fields
or the versioned table described in [working fluids](WORKING_FLUID_MODELS.md).

Each stream contributes these individually fixed/active coordinates:

| Coordinate suffix (prefix `external_stream.heat_in.` or `.heat_out.`) | Unit |
|---|---|
| `inlet_temperature_k` | K |
| `mass_flow_kg_s` | kg/s |
| `cp_j_kg_k` | J/(kg*K) |
| `wall_conductance_w_k` | W/K |

Research coordinates are positive; the lower-level stream model additionally
supports the zero-flow/zero-conductance limits. Fluid labels remain fixed in the
basis and are included in scientific identity. Cp is normally fixed and should
only be active for an intentional scenario study. The existing inventory,
frequency, material, microtube geometry and kinematic parameters coexist in the
same vector. No parameter becomes active automatically.

```toml
[[parameters]]
name = "external_stream.heat_in.mass_flow_kg_s"
unit = "kg/s"
kind = "continuous"
transform = "linear"
initial = 0.05
lower = 0.045
upper = 0.055
```

Replace the corresponding fixed row; do not add a duplicate. The objectives
`maximize_cooling_cop` [1] and `maximize_cooling_power` [W] both work with
`maximum_mechanical_input_power` [W]. Neither fixes the human researcher's later
choice of formulation. This limit bounds indicated input, not actual shaft
power. The fixture does not choose a human-cell objective or power budget.

Reports show the scenario label, inlet/outlet range, mass flow, Cp, capacity
rate, conductance, external cycle heat and average power. Refrigeration tables
make cooling power, indicated input and COP visible. Constraint margins remain
explicit. The cockpit separates stored scientific evidence from deterministic inspection/continuation actions.

## Verification

`AirWallExchanger` describes the explicitly declared air-film model;
`ExternalStreamWallExchanger` describes the generic external stream.
Both expose neutral external heat/outlet diagnostics to
`ExternalStreamWallMachine` and `solve_periodic_wall_machine`.
`tests/test_external_stream_v3.py` verifies their instantaneous thermal parity,
finite-capacity boundaries and signed conservative accounting.
Source/runtime changes prevent incompatible execution resume; stored results are
not rewritten.
