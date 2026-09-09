# Working-Fluid Screening

## Scope

The first-level solver accepts a calorically perfect ideal gas through explicit
`R`, `Cp` and `Cv`. Changing fluid is therefore permitted, but transport
properties used by exchanger sizing remain separate inputs. High-pressure
results must not be labelled physically valid until real-gas and property-
variation checks are available over the simulated pressure-temperature range.

## Helium reference data

NIST reports helium-4 ideal-gas `Cp/R = 2.5` and, at 300 K, thermal conductivity
about `0.1559 W/(m K)` and viscosity about `19.9 uPa s`. The monatomic
calorically-perfect screening values used here are:

```text
R  = 2077.1 J/(kg K)
Cp = 5193.0 J/(kg K)
Cv = 3115.9 J/(kg K)
gamma = 1.6667
```

Primary sources:

- [NIST semiconductor-process gas properties: helium](https://www.nist.gov/pml/sensor-science/fluid-metrology/database-thermophysical-properties-gases-used-semiconductor-9)
- [NIST helium thermophysical-property tables](https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1334.pdf)
- [NIST helium equation of state](https://www.nist.gov/publications/equation-state-thermodynamic-properties-helium)

These rounded constants are screening inputs, not a replacement for the NIST
equation of state.

## Argon reference data

Argon is also monatomic. Rounded calorically-perfect screening values near
ambient temperature are:

```text
R  = 208.13 J/(kg K)
Cp = 520.33 J/(kg K)
Cv = 312.20 J/(kg K)
gamma = 1.6667
```

NIST property data give a thermal conductivity near `0.0177 W/(m K)` and a
dynamic viscosity near `22.7 uPa s` around 300 K and low pressure. These
transport values are not solver constants: they belong in the exchanger and
hydraulic models, with their applicable temperature and pressure range stated.

Primary source:

- [NIST Chemistry WebBook: argon](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440371&Mask=5)

Argon is easier to contain and is generally a more ordinary industrial gas than
helium, but its thermal conductivity is lower than that of air and about nine
times lower than that of helium near 300 K. Its speed of sound is also much
lower. It is consequently a useful monatomic control case, but it is not
expected to beat helium when exchanger compactness and Mach margin dominate.

For two monatomic ideal gases at the same pressure, temperature, geometry,
speed and UA, the first-level thermodynamic solution is identical when the
orifice areas are scaled as

```text
target CdA / reference CdA = sqrt(reference R / target R)
```

This follows because `Cv/R = 3/2` for both gases and because the normalized
orifice mass-relaxation rate is proportional to `CdA * sqrt(R)`. The gas mass
itself is not identical. Relative to the same air reference, matched CdA factors
are approximately `0.3717` for helium and `1.1743` for argon.

## Expected tradeoffs

Helium has much greater thermal conductivity and sound speed than air, which can
reduce exchanger resistance and Mach number. Its larger gamma also creates
stronger adiabatic temperature excursions for a given volume ratio. This can
enable a larger temperature lift but can increase irreversibility or overshoot
if the kinematics are not redesigned.

At the same pressure, helium does not automatically provide greater volumetric
thermal capacity. Increased charge pressure primarily increases power density,
gas inventory, pressure forces and required UA. Within the ideal model, changing
pressure while scaling UA proportionally preserves the intensive cycle and COP.
Pressure must never be presented as a direct COP improvement.

Helium containment is a separate industrial design issue. Leakage, permeation,
seal choice and gas recovery are outside the current thermodynamic model.

## Initial freezer screening

The cooling-cell candidate was tested off-design with a 250.15 K cold reservoir
and 310.15 K hot reservoir, corresponding to explicit 5 K approaches around an
-18 degree Celsius compartment and 32 degree Celsius ambient. Air did not
refrigerate at this point. Helium with approximately matched dimensionless
thermal and hydraulic time scales produced at one bar:

```text
cooling power = 77.92 W
input power = 43.59 W
cooling COP = 1.788
temperature range = 218.6 to 361.1 K
maximum pressure = 2.03 bar absolute
```

Five-bar operation scales this ideal point to approximately 389.6 W cooling and
217.9 W input without materially changing COP. This is evidence that helium is
a promising candidate for the freezer lift, not evidence that the inherited
cooling-cell geometry is optimal. Geometry, volume ratio, UA and speed require a
new constrained search. Real-gas validity must also be evaluated before pressure
is raised further.

## Matched helium-argon check

The current doubled cooling-cell geometry was evaluated at 250.15 K and
310.15 K, one bar charge, 75 W/K on each exchanger and gas-matched CdA values.
This check deliberately isolates ideal-gas similarity; it is not an optimized
freezer candidate:

| Gas | Cooling | Input | COP | Gas mass | Temperature range |
|---|---:|---:|---:|---:|---:|
| helium | 176.40 W | 118.21 W | 1.49225 | 0.01170 kg | 220.53 to 358.74 K |
| argon | 176.44 W | 118.23 W | 1.49229 | 0.11672 kg | 220.53 to 358.75 K |

The tiny difference is caused by rounded property constants and numerical
tolerances. It confirms that fluid selection must include transport,
containment, real-gas validity and hardware scale rather than relying on the
ideal cycle COP alone.

## Commercial COP comparison boundary

A DADA thermodynamic COP must not be compared directly with an appliance energy
label. Compressor catalog COP, complete refrigeration-system COP and annual
cabinet efficiency use different boundaries. For example, current Secop R600a
compressor data at the ASHRAE LBP rating point of approximately -23.3 degrees
Celsius evaporating and 54.4 degrees Celsius condensing include COP values near
1.5, with some operating points and products approaching or exceeding 1.8.

- [Secop SDD compressor-system data](https://www.secop.com/fileadmin/user_upload/technical-literature/leaflets/secop_sdd_system_04-2026_desb101l402.pdf)
- [Secop XV variable-speed compressor data](https://www.secop.com/fileadmin/user_upload/obsolete-compressors/xv72kx_108h7211_r600a_115v_220v_50hz_60hz_04-2020_desd506x202.pdf)

The initial DADA helium point used 250.15 K and 310.15 K reservoirs and reports
gas-cycle mechanical input only. It excludes actuator, bearings, transmission,
fans, pumps and controls. A value near 1.8 is therefore encouraging, but it is
already commercial-compressor territory rather than evidence of superiority.
A defensible comparison requires matched evaporating/cold-boundary,
condensing/hot-boundary and power-boundary definitions.
