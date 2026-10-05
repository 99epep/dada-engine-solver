# Research implementation validation ledger

## Purpose

This historical ledger records the validation methods used while the Research
interface was extracted and extended.

It is not the current test-status page and does not define present physical or
numerical acceptance limits.

For current validation status, known failures and evidence boundaries, see
[`../validation.md`](../validation.md).

For current Research usage, see
[`../DADA_ENGINE_RESEARCH.md`](../DADA_ENGINE_RESEARCH.md).

For detailed source-extraction provenance, see
[`DADA_ENGINE_RESEARCH_MIGRATION_MATRIX.md`](DADA_ENGINE_RESEARCH_MIGRATION_MATRIX.md)
for the retained narrative record. Intermediate machine-readable audit artifacts
and capture scripts were removed after migration; active regression fixtures remain.

## 1. Validation philosophy

### 1.1 Software parity is not physical validation

Numerical parity between an extracted implementation and its historical
implementation establishes software equivalence within a measured tolerance.
It does not validate the underlying physical model or an experimental machine.
Convergence, conservation, model applicability and study feasibility are separate
verdicts; none should substitute for another.

### 1.2 Instantaneous and integrated tolerances

Kinematic and instantaneous model comparisons can be much tighter than full
periodic replay comparisons. Tiny representation differences can change adaptive
integration steps and sampled extrema, amplifying differences in integrated
metrics without changing the intended equations.

Tolerances must be measured and justified for the quantity and backend tested,
not assumed from an instantaneous comparison. No physical acceptance constraint
should be relaxed merely to recover historical parity. Exact tolerances remain
with the regression tests that exercise them.

### 1.3 Historical fixtures and current physics

A historical parity fixture can cease to match after an intentional physical
model revision without implying a regression in the new physical model.
Fixtures intended to document the old model remain unchanged; current-physics
regressions are versioned separately. Present failures and their interpretation
belong on the current validation page, not in this ledger.

## 2. Initial Research extraction

### 2.1 Scientific identity and persistence

Validation covered study scientific identity, candidate identity, exact cache
behaviour, integer decoding, deterministic Sobol continuation, disk resume and
pending-candidate recovery. Runtime/source compatibility is checked separately
from the ability to inspect immutable stored results.

Portable snapshots were exercised without the original study files. A warm start
remains an initial guess followed by periodic validation, not a cached physical
verdict. See the [configuration and persistence reference](../DADA_ENGINE_RESEARCH_REFERENCE.md).

### 2.2 Fixed-pair thermodynamic parity

The initial fixed-mechanism extraction was checked at two levels: kinematic and
instantaneous model equivalence, then full periodic thermodynamic replay.
Comparisons matched machine configuration, mechanisms, volume derivatives,
exchangers, initial states and numerical settings. Statuses and constraint
outcomes were checked alongside physical metrics.

This established wrapper equivalence for the controlled cases. It did not
establish optimality, mechanical efficiency or experimental accuracy.

### 2.3 Offline reporting

Stored records were inspected and compared without a callable physical solver.
Tests separated report inspection from physical reevaluation, preserving the
ability to inspect old results after execution compatibility changes.
Current report and optional replay behaviour is described in the
[cockpit reference](../DADA_ENGINE_RESEARCH_COCKPIT.md).

## 3. Kinematic and mechanism extraction

### 3.1 Dense motion parity

Dense multi-cycle position and derivative comparisons were frozen before
extraction. Checks covered phase, orientation, extrema, available analytic second
derivatives and derivative cross-checks. Mathematical motion representations
were distinguished from the search protocols that had used them.

`tests/test_research_kinematics_v2.py` and `tests/test_research_v2.py` retain this
software evidence. Integrated regressions use looser measured tolerances where
adaptive trajectories are sensitive to tiny representation changes.

### 3.2 Physical mechanisms and artifacts

Physical mechanism artifacts reconstruct through production implementations.
Validation covers complete geometry and settings, loop closure, discrete
assembly branches, serialization round trips and canonical scientific hashes.
Provenance does not redefine scientific mechanism content. Portable campaign
snapshots retain the mechanism inputs needed for resume.

See the [current kinematics reference](../DADA_ENGINE_RESEARCH_KINEMATICS.md),
[primary catalogue](../PRIMARY_FOUR_BAR_FAMILIES.md) and
[six-bar catalogue](../SIX_BAR_MECHANISM_FAMILIES.md).

### 3.3 Mechanical diagnostics

Mechanical comparisons preserve coordinate frames, branch choices, scale
conventions and the method used to obtain extrema. Sampled screens are not
continuous-angle guarantees. Nominal margins are not manufacturing-tolerance
probabilities; different constraints retain their separate units and meanings.

### 3.4 Synthesis handoff boundary

Tests cover primary, downstream, full-mechanism and paired coordinate releases,
fixed branches, independent cylinder mechanisms and multiple retained families.
These are declarative contracts, not evidence that automatic mechanism fitting
or saturation execution exists. The
[synthesis method](../MECHANISM_SYNTHESIS_SEARCH.md) distinguishes discovery
proxies, realizable motion and final thermodynamic assessment.

## 4. External thermal boundaries and property tables

### 4.1 Conservation and signed refrigeration

External-stream work was checked through neutral/legacy parity where applicable,
wall/stream conservation, outlet behaviour, signed refrigeration energy balance
and periodic replay. Fixed/active ownership of stream coordinates was tested
through the same Research campaign interface.

The [external-stream model](../EXTERNAL_STREAM_THERMAL_MODEL.md) defines the
physical boundary and omissions. Signed indicated input is not measured shaft
power, and a numerical refrigerator fixture does not validate cooling hardware.

### 4.2 Compiled property reconstruction

Exact-versus-table checks use a known analytic validation fluid with explicit
instantaneous and integrated tolerances. They exercise conservative-state
reconstruction, transported enthalpy and compiled lookup without Python property
callbacks, as well as periodic replay and property-table identity.

This demonstrates the compiled architecture, not real-fluid accuracy. See
[working-fluid models](../WORKING_FLUID_MODELS.md) for the thermodynamic,
transport and hydraulic boundaries.

### 4.3 Domain rejection

Tests cover table-domain and invalid-cell rejection, immutable table storage and
incompatible hydraulic assumptions. Unavailable properties must produce an
explicit rejection rather than silent extrapolation. A transport temperature
domain is not a phase-domain guarantee.

## 5. Later Research infrastructure

### 5.1 Capacity scaling

Factor-one reconstruction should preserve scientific inputs and results under
matched initialization. Extensive quantities follow the declared capacity
transformation; intensive inputs and constraints are not silently rescaled.
Checks cover active bounds, integer handling, portable mechanism artifacts and
non-mutation of source histories. Capacity scaling creates a new scientific
study; see [capacity scaling](../DADA_ENGINE_RESEARCH_CAPACITY.md).

### 5.2 Journal-first cockpit

Validation covered TTY/non-TTY output, exact candidate selection, offline
sorting/filtering, constraint evidence and report generation without
reevaluation. Solver stderr capture was tested for restoration after both normal
returns and exceptions.

Persistence checks cover journal/recovery ordering, recovery-only completion,
residual recovery after append, torn tails, state-publication failures and legacy
result inspection. UI changes must not change candidate identity or Sobol
continuation. See the [cockpit reference](../DADA_ENGINE_RESEARCH_COCKPIT.md).

### 5.3 Reference-pressure inventory

Checks retain explicit-inventory compatibility while exercising derived charge
across supported motion families, exchanger/dead-volume changes and simultaneous
cylinder volumes. Deterministic reconstruction is independent of the mechanical
screening grid and participates in scientific identity.

Wall/reservoir warm-state handling and capacity-rescale interaction are tested
without silently reusing an incompatible inventory. The current
[reference-pressure policy](../DADA_ENGINE_RESEARCH_REFERENCE.md#filling-at-a-reference-pressure-and-maximum-total-gas-volume)
defines this contract.

### 5.4 Limit ownership

Model domains, study requirements, search bounds and numerical settings have
distinct owners. Historical study constraints must not become universal solver
limits merely because an old candidate used them. See the
[limit-ownership audit](RESEARCH_LIMIT_OWNERSHIP_AUDIT.md).

## 6. Evidence that remains versioned

Historical fixtures under `tests/fixtures/` and `tests/data/` remain available.
The kinematics/Research tests named above retain extraction and wrapper checks;
`tests/test_external_stream_v3.py`, `tests/test_tabulated_fluid_v3.py` and
`tests/test_research_v3.py` retain external-boundary and compiled-property checks.
Campaign, reporting and infrastructure tests retain identity and recovery checks.
The migration matrix and audit inventories preserve detailed source provenance.
No saved campaign report is required to read this ledger.

## 7. What this ledger does not prove

Passing software parity checks does not establish physical validation of the
machine, universality of a correlation, manufacturing robustness or global
optimality. A best feasible candidate is evidence for its declared study only.
Current physical limitations and known failures remain in the current validation
reference; this ledger supplies neither a replacement status nor new limits.

## 8. Current references

- [Validation and evidence](../validation.md): current status and physical limitations.
- [Research guide](../DADA_ENGINE_RESEARCH.md): current usage.
- [Research reference](../DADA_ENGINE_RESEARCH_REFERENCE.md): configuration and persistence.
- [Kinematics](../DADA_ENGINE_RESEARCH_KINEMATICS.md): current family and artifact APIs.
- [Cockpit](../DADA_ENGINE_RESEARCH_COCKPIT.md): reporting and diagnostic presentation.
- [Migration matrix](DADA_ENGINE_RESEARCH_MIGRATION_MATRIX.md): detailed extraction provenance.
