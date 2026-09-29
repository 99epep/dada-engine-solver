# Candidate capacity scaling and offline comparison

## Capacity scaling

`rescale` creates a new V2/V3 study and a new portable machine basis from a
verified stored candidate. It accepts a campaign directory or a standalone
Research evaluation artifact, an exact ID or an unambiguous ID prefix, and one
explicit mode, `capacity`. It never integrates or modifies the source.

```sh
research rescale examples/human_cell_stage0/campaign \
  --candidate a16ed2c7 --mode capacity --factor 5 \
  --output outputs/my_human_cell_A5/study.toml
research validate outputs/my_human_cell_A5/study.toml
research evaluate outputs/my_human_cell_A5/study.toml \
  --output outputs/my_human_cell_A5/evaluation.json --budget 3m
```

Here `research` is the shell helper documented in
[DADA_ENGINE_RESEARCH.md](DADA_ENGINE_RESEARCH.md); `dada-research` is the
installed entry point. Output paths must be unused.

The selected candidate becomes the new initial point, including every active
kinematic and intensive coordinate. Existing fixed/active ownership is retained.
Implicit machine defaults are materialized so that edited basis defaults cannot
silently leave the new study at the old capacity.

| Quantity | Capacity action |
| --- | --- |
| Total swept volume, gas inventory | Multiply by `s` |
| Cylinder minimum/maximum and reservoir exchanger reference volumes | Multiply by `s` |
| Tube counts | Multiply by `s`; require integral physical counts |
| Tube flow area, gas-film conductance, wall material capacity | Rebuild from production geometry |
| Additional internal exchanger volume and extra wall capacity | Multiply by `s` |
| External stream mass flow and declared wall conductance | Multiply by `s` |
| Legacy air mass flow | Multiply by `s`; air film remains geometry-derived |
| Hydraulic configuration CdA and fixed exchanger outlet CdA | Multiply by `s` |
| Count-ratio outlet CdA policy | Keep the reference CdA/count slope; scaled tube count supplies the factor once |
| Legacy reservoir conductance references | Multiply by `s` |
| Frequency, kinematics, volume ratio, clearance ratios | Preserve selected candidate values |
| Tube length/diameter/wall thickness/pitch/header depth, materials, temperatures, Cp | Preserve selected candidate values |
| Constraints, validity thresholds and numerical tolerances | Preserve, without relaxation |

This is capacity scaling, **not geometric length similarity**. Physical
mechanism lengths, if declared, remain unchanged: this can imply changed piston
area and loads, which this operation does not certify mechanically. The
production square-pitch rectangular bank packing is retained. Integer row and
column counts mean header volume and envelope are not necessarily exactly `s`
times their previous values. No independent UA, header volume, flow resistance,
or pressure is fitted to force the desired output power.

An active extensive coordinate retains its transform and has initial value and
bounds multiplied by `s`. Integer bounds become `ceil(s*lower)` and
`floor(s*upper)`; collapsed bounds are rejected. The selected tube count must
be integral without rounding (floating-point roundoff within 1e-9 count is
accepted). Intensive bounds remain unchanged. A factor must be finite and
positive. V1 rescaling, custom hydraulic closures, unknown schemas and unsafe
count conversions fail explicitly; they need deliberate migration or a
family-specific implementation first. V1 inspection remains supported.

Source warm states, when present and inventory-compatible, are scaled as
conservative initial guesses, not reused as cached solutions. The regular
wall-capacity adjustment and periodic convergence remain in force. Uniform
initialization stays uniform. A factor of one preserves exactly the same
physical inputs, but receives a new scientific identity because it explicitly
records its parent and transformation.

The basis `provenance.capacity_scaling` contains version 1, mode, factor,
source candidate/study/definition identities, source path and basis digest,
main before/after values, bound policy, unchanged-limit policy and header
packing convention. Existing basis provenance is retained. Associated mechanism
artifacts are copied. No schema bump or history migration is required: this
uses existing V2/V3 provenance and parent-candidate fields. Exact cache and
runtime/source compatibility rules remain unchanged.

Multiplication can conceptually be reversed with `1/s`. Integer bound clipping
and floating-point rounding are explicitly not claimed to be exactly reversible;
the original source and recorded before values remain the reference.

## Valve chronology availability

The ten-state wall integration path (`external_stream_wall` and historical
`air_wall`) uses continuous ideal diodes but does **not** locate or record
pressure-crossing roots. `WallDiagnosticCycle` previously supplied `events=()`
and placeholder closed topologies to generic diagnostics. The replay did not
lose recorded transitions: no transitions had been collected in that integrator.
The separate eight-state integrator does record continuous-diode events through
its root-observation code; that behavior is retained.

Wall-cycle topology is now `unavailable`, with reason
`wall_integrator_does_not_record_valve_events`. Any empty generic sequence is
also unavailable (`no_usable_valve_event_sequence`), rather than automatically
non-nominal. Nonempty recorded sequences retain the existing cyclic nominal /
non-nominal classification. No sampled flow sign is silently promoted to an
accurate valve event, and no thermodynamic or valve equation is changed.

Old records remain immutable. Inspection adds `topology_display` with
`unavailable` for an empty legacy sequence; the original classification and
reasons stay in `diagnostics`, and the display correction records the old
classification. Local reflux remains a separate signed-flow diagnostic.

## Offline table and volume plots

Click visible column headings to toggle ascending/descending order, or use
**Sort by** and **Direction**. Sort keys include COP, cooling power, indicated
input and gas power, pressure/temperature/absolute-flow maxima, cycles, duration,
all active coordinates, and absolute/relative constraint margins. Missing values
stay last in either direction, with candidate ID as a deterministic tie-breaker.
Status and validity filters compose with sorting. The main table exposes
chronology availability, reflux detection and the worst minimum signed flow.

```sh
research report examples/human_cell_stage0/campaign \
  --candidate a16ed2c7 --candidate 256fdb49 --plots volumes \
  --html outputs/my_human_cell_comparison.html
research compare outputs/human_cell_A5/evaluation.json \
  outputs/human_cell_B5/evaluation.json --plots volumes \
  --html outputs/my_scaled_comparison.html
```

`--plots volumes` reconstructs only the selected candidates, samples 721 angles
including both cycle endpoints, and embeds both cylinder volumes in the
standalone HTML. Production factory conventions apply the operation direction
once. The horizontal axis is solver cycle angle [deg]; stored volumes are in
m^3 and the figure displays litres. Candidate colors and solid-small /
dashed-large styles identify the curves. Duplicate selections reuse the same
samples within a report. Historical result artifacts contain no volume arrays,
so reconstruction is needed; no ODE, thermal replay, or optimization is run.
Reports include full selected evidence and lightweight whole-campaign progress.

The structured `plots.volumes` dataset is deliberately separate from thermal
metrics, leaving room for future plot types without implying thermal replay is
available. Only volumes are currently accepted. Standalone V2/V3 evaluations and
campaign snapshots can be reconstructed using current production models, even
when execution resume is incompatible; the runtime warning remains visible.
Reports may regenerate derived HTML, need no network/server, and preserve scientific input
provenance. Cross-study comparison is evidence, not a combined optimization
ranking.

## Bounded Human Cell validation

See [the demonstration and measured results](../outputs/human_cell_stage0/README.md).
No optimization was launched. A×1 and an original-A replay with identical uniform
initialization have exactly equal recorded metrics on this runtime. Small
differences from the historical campaign result arise from its different warm
start and finite convergence tolerance, not changed factor-one physics.
