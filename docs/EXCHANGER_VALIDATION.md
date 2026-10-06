# Microtube gas-model verification and experimental comparison

## Verification versus validation

Analytical and numerical checks establish specific properties of the implemented
model, not experimental accuracy merely because a correlation is implemented.
Tests cover the incompressible limit, signed-flow symmetry, constant-property
Poiseuille equivalence, half-tube allocation, laminar asymptotes, entry effects,
turbulent inversion and transition continuity. Domain checks cover temperature,
compressibility, rarefaction and incompatible accommodation conventions.
Variable-film and wall tests retain heat exchange at rest and check conservation;
compiled/reference parity is a separate software property.

Equations, sources and applicability belong to the
[microtube model](MICROTUBE_GAS_MODEL.md). Finite-stream and wall energy boundaries
are defined in [external streams](EXTERNAL_STREAM_THERMAL_MODEL.md).
Being inside a correlation's domain does not validate pulse response or the
assembled machine. Current test status belongs to [validation](validation.md).

## Graur/Ewart published-regression check

The check reproduces the published pressure-ratio-5 polynomial fits from
Ewart/Graur Tables 1–2. It is not independent validation against raw measurements:
`measured=null` and `reference_kind=published_regression_not_raw_measurement`
are the appropriate report fields. The published coefficient uncertainty and
VHS mean-free-path convention are retained.

Published silica TMAC values do not calibrate a metallic surface. Changing
declared A1/A2 coefficients is a sensitivity exercise, not calibration or an
implicit thermal-jump model. Source equations and surface-transfer limitations
are in the [slip reference](MICROTUBE_GAS_MODEL.md#gas-knudsen-hs-and-gas-slip-second-order).

## Doty hydraulic comparison

Source: [Doty et al. (1991)](https://dotynmr.com/download/pubs/1991_HTE_Doty_HeatExchanger.pdf),
Eq. 7 and Tables 1–2. SI measurements and reported uncertainties are retained in
[the nitrogen data](../tests/data/doty_1991_nitrogen_reference.csv) and
[the helium data](../tests/data/doty_1991_helium_reference.csv).

The table is the **recorded September 2026 comparison**, not a fresh prediction
with the current transport model. It uses 309 equal-flow tubes, 0.33 mm internal
diameter and 127 mm length; temperature is `(T3+T4)/2`, reported pressure is
interpreted as mean tube pressure, and no header loss is added to the tube-only
comparison. Density is ideal-gas; viscosity uses Sutherland nitrogen or the
then-used NIST helium interpolation. Current helium uses Arp/Hands-Arp; see
[transport provenance](MICROTUBE_GAS_MODEL.md#species-dependent-dilute-transport-version-dilute_species_v2).

| Gas | Flow, g/s | Measured tube dP, Pa | Recorded prediction, Pa | Signed relative error |
| --- | ---: | ---: | ---: | ---: |
| N2 | 0.470 | 5100 | 3903 | -23.46% |
| N2 | 0.840 | 9500 | 6819 | -28.22% |
| N2 | 0.930 | 11100 | 7574 | -31.77% |
| He | 0.117 | 2300 | 3319 | +44.30% |
| He | 0.079 | 1500 | 2206 | +47.05% |
| He | 0.213 | 4400 | 5477 | +24.48% |

Signed error is `(predicted - measured) / measured`, calculated before rounding
the displayed pressures. None of these six recorded predictions lies inside the
reported pressure-drop output uncertainty. Nitrogen loss is underpredicted;
helium loss is overpredicted. No opaque multiplier is fitted to hide the
mismatch, and unreported geometric or property-reduction uncertainties are not
invented to obtain agreement.

The independent helium reconstruction is distinct from Doty's own calculated
column. At fixed mean pressure the compressible Poiseuille and mean-density
forms are algebraically equivalent; that identity is verification, not agreement
with experiment. The [Doty screening reference](DOTY_SCREENING.md) describes
relative scaling, extrapolation and the original property convention.

## What is not validated

Measured UA, effectiveness and boundary temperatures exist, but Doty's shell
side and axial thermal field have not been independently reconstructed. The
hydraulic comparison therefore reports no predicted whole-bank UA. Feeding
measured UA back into the model would make a thermal-validation claim circular;
integral validation requires an independent opposite-side/axial thermal model.

Neither the hydraulic comparison nor a valid steady correlation establishes
pulsed-flow accuracy, gas/metal surface transferability, flow distribution or
whole-machine performance. These remain separate experimental questions.

## Reproduction

Run the current equation, domain, geometry and experimental-comparison tests
from a source checkout:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m pytest -q \
  tests/test_doty.py tests/test_microtube_gas_model.py \
  tests/test_microtube_geometry.py tests/test_microtube_circular.py \
  tests/test_microtube_developing_entry.py tests/test_microtube_transition.py
```

The optional `examples/validate_microtube_gas_model.py` script generates Doty
hydraulic comparisons and Graur fit checks using the current implementation,
without replaying a motor campaign. It writes `outputs/microtube_gas_validation.json`;
that path may contain an older report and is overwritten when the script runs.
Its current helium predictions need not reproduce the recorded table above.
The tests exercise the same comparison functions without publishing that report.
See [validation](validation.md) for the complete suite and
[numerical execution](SOLVER_PERFORMANCE.md) for backend verification.
