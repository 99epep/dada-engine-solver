# Capacity scaling

## What it does

`research rescale --mode capacity` creates a new portable study from a verified
stored candidate. Capacity scaling changes extensive inputs while retaining the
selected intensive coordinates; it is **not geometric similarity** or a promise
that power will scale exactly. It performs no thermodynamic integration and does
not modify the source. Constraints remain unchanged, without relaxation.

## Command

```sh
research rescale SOURCE \
  --candidate CANDIDATE --mode capacity --factor 5 \
  --output scaled/study.toml
research validate scaled/study.toml
```

`research` is the checkout helper in the [user guide](DADA_ENGINE_RESEARCH.md);
the installed entry is `dada-research`. SOURCE is a Research campaign or compatible
standalone evaluation. Select a candidate by full ID, unambiguous prefix, or the
production selectors `best`/`second`.

The factor must be finite and strictly positive. Output must end in `.toml`;
the study, associated basis and mechanism-artifact destinations must all be new.
The command validates the generated study before publishing it and does not
launch a campaign. Generated execution settings use a 512-attempt default;
this scheduling convention is not a physical criterion. See the
[technical reference](DADA_ENGINE_RESEARCH_REFERENCE.md) for execution semantics.

## Scaling contract

Research coordinate scaling and physical-basis reconstruction are distinct.
The selected candidate supplies the new fixed values and active initial point;
implicit machine defaults are materialized where necessary.

| Research coordinate | Action |
|---|---|
| Total swept volume; explicit total gas inventory | Multiply by `s` |
| Microtube count; declared additional internal gas volume | Multiply by `s` |
| External-stream mass flow; declared external wall conductance | Multiply by `s` |
| Legacy external-air mass flow | Multiply by `s` |
| Frequency, volume/clearance ratios and selected kinematics | Preserve |
| Individual tube dimensions, pitch/packing settings and materials | Preserve |
| Temperatures and specific properties | Preserve |

The rebuilt basis also scales cylinder minimum/maximum volumes, reservoir
exchanger reference volumes, legacy conductance references, additional wall heat
capacity and supported hydraulic CdA quantities. Tube areas, tube-wall capacity
and gas-film conductance are reconstructed from production geometry and models,
not treated as independent target-performance inputs.

Study constraints, mechanical requirements, validity thresholds, screening and
numerical tolerances are preserved. A larger machine may therefore violate a
constraint that the source satisfied; rescaling does not certify feasibility.

## Active parameters and integer counts

Fixed/active ownership and active transforms are retained. Every declaration is
centered on the selected physical candidate:

- An extensive fixed value or active initial value is multiplied by `s`.
- Extensive active bounds are multiplied by `s`.
- Intensive initial values come from the candidate; their bounds stay unchanged.
- Tube counts must remain positive integers. No opportunistic rounding is allowed;
  only absolute roundoff within `1e-9` count is accepted.
- Integer bounds use `ceil(s * lower)` and `floor(s * upper)`.
  Collapsed or invalid intervals are rejected, including equal endpoints.

See the [technical reference](DADA_ENGINE_RESEARCH_REFERENCE.md) for general
parameter declarations and search encoding.

## Geometry is rebuilt, not geometrically similar

The transformation starts from the production machine constructed for the selected
candidate. It creates a new basis and declarations, then validates them through
normal study loading. It does not blindly multiply every old object field.

The selected bank geometry remains authoritative. Legacy rectangular banks retain
their discrete packing; row/column changes can make header volume and envelope
non-proportional to `s`. Circular banks use the continuous approximation described
below. Neither independent UA nor pressure loss nor pressure is adjusted to force
power similarity.

Declared physical mechanism lengths are not homothetically scaled. Retaining
stroke while changing swept volume can change required bore and mechanical loads;
this operation does not certify those loads.

For outlet valves, fixed/source CdA is extensive. A count-ratio policy preserves
the reference CdA/count slope: the new tube count supplies the factor once, not
twice. Circular microtube ideal-diode area is instead derived from conduit area.
Unknown/custom hydraulic closures are refused rather than extrapolated.

## Circular collectors

For `circular_triangular_frustum_v1`, individual tube dimensions, `pitch_ratio`,
`conduit_area_ratio` and collector half-angle stay fixed. With an admissible integer
count change `N -> s*N`, the implemented geometry gives:

- bundle and conduit diameters proportional to `sqrt(s)`;
- collector height proportional to `sqrt(s)`;
- collector volume proportional to `s^(3/2)`;
- conduit/ideal-diode area proportional to `s`.

This follows exactly from the continuous triangular-cell envelope and conical
frustum formulas, apart from floating-point roundoff. It is a consequence of the
chosen geometry, not a scaling error. Tube gas volume and declared additional
volume scale by `s`, but total hold-up need not. Reference-pressure inventory is
therefore re-derived from the final connected volumes rather than forced to scale
linearly. See the [microtube model](MICROTUBE_GAS_MODEL.md) for formulas and domains.

## Warm start and identity

A stored source warm state is scaled as a conservative initial guess only when
its original inventory matches the selected design and charge uses explicit
inventory. Otherwise that stored guess is discarded in favor of uniform
initialization. Uniform initialization remains uniform. Wall-capacity adjustment
and normal periodic convergence still apply; a guess is not a cached solution.

Factor `1` preserves the selected candidate's physical inputs, as checked by the
production regression tests. It still creates a new study with explicit parent
provenance and a new scientific identity, not an alias or cache hit.

## Provenance and portability

The basis records `provenance.capacity_scaling`: version, mode, factor, parent
candidate, source study/definition, source location and basis digest, before/after
changes, bound conventions, unchanged constraints, geometry reconstruction and
warm-start policy. Existing provenance is retained. The generic geometry note
mentions discrete packing; the actual bank family determines whether discrete
packing or the circular continuum formula applies.

Required mechanism artifacts are copied and referenced locally. Move the new TOML
together with its associated basis and artifacts; reconstruction does not require
the original campaign directory. Source paths retained in provenance are evidence,
not live dependencies.

Scaling by `s` can conceptually be followed by `1/s`, but bitwise reversibility is
not promised: integer bounds, ceil/floor, discrete packing, roundoff and geometric
reconstruction can intervene. The original source and recorded `before` values
remain the reference.

## Unsupported cases and related references

Only schema 2/3 (V2/V3) machine studies are supported. Schema 1 remains inspectable
but requires explicit migration before scaling. `local_regions_v1` is rejected:
rescale the global source first, then optionally refine the new study.

Unsupported schemas, custom hydraulic flow models, unknown hydraulic quantities,
non-integral counts, overflow, invalid bounds and existing destinations fail
explicitly. Unknown declarations must pass the normal production schema checks;
there is no fallback scaling rule for unrecognized physical quantities.
These refusals protect the declared scientific transformation.

Result tables, report reconstruction and display corrections belong to the
[cockpit reference](DADA_ENGINE_RESEARCH_COCKPIT.md). Current validation evidence
belongs to [validation](validation.md), not to a particular historical campaign.
