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
| Microtube count and individual length | Coupled total-length transformation below |
| Declared additional internal gas volume | Multiply by `s` |
| External-stream mass flow; declared external wall conductance | Multiply by `s` |
| External-air mass flow | Multiply by `s` |
| Frequency, volume/clearance ratios and selected kinematics | Preserve |
| Inner/outer tube diameters, pitch/packing settings and materials | Preserve |
| Temperatures and specific properties | Preserve |

The rebuilt basis also scales cylinder minimum/maximum volumes, reservoir
exchanger reference volumes, conductance references, additional wall heat
capacity and supported hydraulic CdA quantities. Tube areas, tube-wall capacity
and gas-film conductance are reconstructed from production geometry and models,
not treated as independent target-performance inputs.

Study constraints, mechanical requirements, validity thresholds, screening and
numerical tolerances are preserved. A larger machine may therefore violate a
constraint that the source satisfied; rescaling does not certify feasibility.

## Microtube length and integer counts

For each exchanger independently, capacity multiplies total physical microtube
length, not count alone. For source count `N`, individual length `L` and factor `s`:

```text
target total tube length = s N L
N_ideal = N sqrt(s)
N' = round(N_ideal), nearest integer with ties to even
L' = s N L / N'
N' L' = s N L
```

A quantized count below one is rejected, not clamped. Length must remain finite
and strictly positive. Inner and outer diameters, pitch and material remain
unchanged. Factor `1` retains count and length exactly. The length correction
preserves the total-length invariant even when `N_ideal` is not integral.

Fixed/active ownership and active transforms are retained. Every declaration is
centered on the selected physical candidate:

- Ordinary extensive values and bounds are multiplied by `s`.
- Microtube count initial values and bound endpoints use the same nearest-integer
  quantization of their value times `sqrt(s)`.
- Microtube length initial values use the coupled rule above; length bounds are
  multiplied by the selected candidate's actual length factor `s N / N'`.
- Intensive initial values come from the candidate; their bounds stay unchanged.

Collapsed or invalid active intervals are rejected, including equal endpoints;
no bounds are widened silently. These are independent search bounds around the
transformed candidate, not a coupled constraint on every future count/length pair.
See the [technical reference](DADA_ENGINE_RESEARCH_REFERENCE.md) for general
parameter declarations and search encoding.

## Geometry is rebuilt, not geometrically similar

The transformation starts from the production machine constructed for the selected
candidate. It creates a new basis and declarations, then validates them through
normal study loading. It does not blindly multiply every object field.

The selected bank geometry remains authoritative. Rectangular banks retain
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

For `circular_triangular_frustum_v1`, tube diameters, `pitch_ratio`,
`conduit_area_ratio` and collector half-angle stay fixed. Ignoring integer
quantization, the continuous trend is:

- count and individual tube length proportional to `sqrt(s)`;
- bundle/conduit diameters and collector height proportional to `s^(1/4)`;
- collector volume proportional to `s^(3/4)`;
- conduit/ideal-diode area proportional to `sqrt(s)`.

The actual geometry uses the retained integer count: for `r = N'/N`, diameters
and collector height scale by `sqrt(r)`, collector volume by `r^(3/2)` and conduit
area by `r`. Production reconstructs these quantities; rescale does not force them.
Tube gas volume and declared additional volume scale by `s`, but total hold-up
need not. Reference-pressure inventory is therefore re-derived from the final
connected volumes rather than forced to scale linearly. See the
[microtube model](MICROTUBE_GAS_MODEL.md) for formulas and domains.

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

The basis records `provenance.capacity_scaling` version `2`: mode, factor, parent
candidate, source study/definition, source location and basis digest, before/after
changes, bound conventions, unchanged constraints, geometry reconstruction and
warm-start policy. It identifies ordinary extensive scaling, the coupled microtube
rule, ideal/retained counts, actual length factors and unchanged diameters. The
source/output search domains record normalization of a local source. Existing
provenance is retained without copying local regions into the output search.

Required mechanism artifacts are copied and referenced locally. Move the new TOML
together with its associated basis and artifacts; reconstruction does not require
the original campaign directory. Source paths retained in provenance are evidence,
not live dependencies.

Scaling by `s` can conceptually be followed by `1/s`, but bitwise reversibility is
not promised: integer quantization of counts/bounds, discrete packing, roundoff and geometric
reconstruction can intervene. The original source and recorded `before` values
remain the reference.

## Unsupported cases and related references

Only current schema-3 machine studies are supported. Global and `local_regions_v1`
candidates are accepted: local origin describes how the candidate was found,
not an ambiguity in its physical machine. A local source produces a global
`fixed_global_bounds` Sobol study retaining `seed` and `scramble`. Source regions,
centers, radius, allocation and center scheduling are discarded. The transformed
candidate supplies the new initial point; no parent campaign is required.

Unsupported schemas, custom hydraulic flow models, unknown hydraulic quantities,
counts quantized below one, nonpositive/nonfinite lengths, overflow, invalid bounds and existing destinations fail
explicitly. Unknown declarations must pass the normal production schema checks;
there is no fallback scaling rule for unrecognized physical quantities.
These refusals protect the declared scientific transformation.

Result tables, report reconstruction and display corrections belong to the
[cockpit reference](DADA_ENGINE_RESEARCH_COCKPIT.md). Current validation evidence
belongs to [validation](validation.md).
