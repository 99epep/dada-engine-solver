# Air-source microtube exchanger design

## User-selected boundary

Use external air streams for both demonstrator sources, initially entering at
448.15 K and 298.15 K. These are inlet temperatures, not temperatures imposed
throughout the exchanger. External air cools on the hot side and warms on the
cold side. Its heat-capacity rate and pressure loss must be modeled. The
working gas remains a separate closed circuit inside the tubes; external air
must not enter that circuit. Nitrogen measurements remain experimental anchors,
not exact air-property data.

Report thermal input delivered to the motor separately from electrical heater
input. Report fan consumption and mechanical losses separately from indicated
gas power. For 1 kW transferred and an illustrative 20 K external-air change,
using an explicit approximate cp of 1005 J/(kg K), the required air mass flow
is about 0.050 kg/s. This is an energy-balance illustration, not a selected flow.
An inlet at 25 degrees Celsius is not a gas temperature that can be maintained
throughout a finite-duty cooler without a temperature approach.

## Geometry implementation and optimization boundary

`MicrotubeBank` in `exchangers/microtube_geometry.py` calculates internal tube
volume, surface, wall material volume and two parametric gas plenums. It uses
square-pitch rectangular packing, not a reconstruction of Doty's triangular
prototype. Pitch, header depth and extra connection volume are explicit inputs.
The envelope excludes casing walls, tube sheets, external ducting and insulation.

For circular tubes, `A = N*pi*d*L`, `V = N*pi*d^2*L/4`, and laminar tube pressure
loss at fixed gas state and mass flow is proportional to `L/(N*d^4)` (the same
Poiseuille dependence used in Doty's Eq. 7). At fixed diameter, increasing N by
factor k and reducing L by k keeps tube area and volume unchanged while reducing
tube friction by k squared. Header volume generally grows with bundle face area
for fixed header depth. Thus more short tubes do not necessarily increase tube
volume, but can increase total working-gas hold-up through their collectors.

Header fluid volume is straightforward once a shape is specified. Header loss,
flow distribution and manufacturability are less certain. They require explicit
loss scenarios and geometric checks, not an arbitrary claim of exact prediction.
The implemented tube pressure loss excludes entrance/header losses, compressible
variation along the passage and pulsation. Reynolds number is returned for an
applicability check; no turbulent prediction is inferred from the laminar law.

A future sizing search must vary count, length, diameter, pitch and collectors
jointly, recalculating cycle efficiency with hold-up and pressure loss. Surface
alone does not guarantee unchanged heat transfer: entry effects, external air
flow and axial wall conduction can change when tubes become shorter.

## Heat transfer during flow pauses

The existing cycle already retains exchanger gas mass and energy and applies
`UA*(T_source-T_gas)` at zero through-flow. It does not freeze the gas during
pauses. What it lacks is independent wall storage and finite external airflow.

The next low-order extension should retain the conservative gas balance and add
wall energy per exchanger. For example, wall storage receives external heat
and loses `G_gas_wall*(T_wall-T_gas)` to the working gas. Heat transfer remains
active at zero port flow. External airflow requires a finite-capacity heat
exchange closure; through-flow effectiveness formulas must not erase trapped-
gas heat transfer when working-gas flow is zero. Wall heat capacity is derived
from wall volume and explicit material properties, with tube-sheet/header mass
added separately. Gas/wall and air/wall conductances cannot be uniquely inferred
from a single measured overall UA, even with gas on both sides.

The user authorized accounting for pauses; the purpose and modest state-count
cost of this extension have been explained. An opt-in conservative motor coupling is now implemented; see
[AIR_WALL_COUPLING.md](AIR_WALL_COUPLING.md). It is not yet Doty-calibrated. First reproduce Doty's steady data, then compare wall-storage and
external-resistance sensitivity before increasing spatial resolution. Do not
claim a CPU speedup or a precise added runtime before benchmarking stiffness.

## Verification status

Geometry invariants, the N/L friction tradeoff, signed/zero hydraulic flow and
invalid inputs are tested. The current module is the geometry foundation for a
Doty-style sizer, not a complete calibrated exchanger or an optimized design.
