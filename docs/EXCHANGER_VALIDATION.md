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

Doty's pressure-loss measurements are tube-side only, not combined shell/tube
loss. The independent helium reconstruction in `tests/test_doty.py` uses 309
equal-flow tubes, 0.33 mm internal diameter and 127 mm length, with no header
loss added. Representative temperature is `(T3 + T4) / 2`; density is ideal-gas
helium at the reported pressure, and viscosity is interpolated from NIST values
between 325 and 350 K. This is an explicit reconstruction convention, separate
from the production helium transport closure described in the
[transport reference](MICROTUBE_GAS_MODEL.md#species-dependent-dilute-transport-version-dilute_species_v2).

The tests independently reproduce Eq. 7 and compare with Doty's calculated
column and measured losses. The points are laminar; the independent predictions
remain above the helium measurements under these assumptions. No empirical
multiplier is fitted to force agreement, and unreported geometric or
property-reduction uncertainties are not invented to obtain agreement.

At fixed mean pressure the compressible Poiseuille and mean-density forms are
algebraically equivalent; that identity verifies equations, not experimental
accuracy. The [Doty screening reference](DOTY_SCREENING.md) describes relative
scaling, extrapolation and the declared property convention.

## What is not validated

Measured UA, effectiveness and boundary temperatures exist, but Doty's shell
side and axial thermal field have not been independently reconstructed. The
hydraulic comparison therefore reports no predicted whole-bank UA. Feeding
measured UA back into the model would make a thermal-validation claim circular;
integral validation requires an independent opposite-side/axial thermal model.

Neither the hydraulic comparison nor a valid steady correlation establishes
pulsed-flow or reflux accuracy, gas/metal surface transferability, flow
distribution or whole-machine performance. These remain separate experimental questions.

## Reproduction

Run the current equation, domain, geometry and experimental-comparison tests
from a source checkout:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m pytest -q \
  tests/test_doty.py tests/test_microtube_gas_model.py \
  tests/test_microtube_geometry.py tests/test_microtube_circular.py \
  tests/test_microtube_developing_entry.py tests/test_microtube_transition.py
```

See [validation](validation.md) for the complete suite and
[numerical execution](SOLVER_PERFORMANCE.md) for backend verification.
