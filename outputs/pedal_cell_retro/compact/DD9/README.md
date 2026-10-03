# f013 long/thin 100x exchanger study

Place `study.toml` at:

`outputs/pedal_cell_retro/long_thin_f013_100x/study.toml`

It deliberately reuses the existing DD7 machine basis through:

`../compact/DD7/study.basis.json`

## Nominal transformed center

Source candidate: `f013861305dda092291abf0d6af55403127f98684a3c400dc93e17d427261f75`

- heat_in: 68,517 -> **685 tubes**
- heat_out: 81,342 -> **813 tubes**
- initial inner diameter: **1.0 mm**
- wall thickness: **0.10 mm**
- pitch ratio: **1.25**
- conical collector half-angle: **30 deg**
- conduit area ratio: **1.0**
- ideal lossless diode policy: `geometry_conduit_area_v1`
- frequency remains **2.8 Hz**
- total swept volume remains **5 L**

The initial lengths were chosen to retain approximately the source candidate's
gas-side internal tube area after the ~100x count reduction:

- heat_in length: **0.251 m**
- heat_out length: **0.333 m**

Approximate geometry of the initial transformed center, before thermodynamic evaluation:

| quantity | heat_in | heat_out |
|---|---:|---:|
| bundle diameter | 41.2 mm | 44.9 mm |
| conduit diameter | 26.2 mm | 28.5 mm |
| one collector axial height | 13.0 mm | 14.2 mm |
| total exchanger envelope length | 277 mm | 362 mm |
| tube internal heat-transfer area | 0.539 m² | 0.851 m² |
| tube gas volume | 135 cm³ | 213 cm³ |
| both collectors gas volume | 24 cm³ | 31 cm³ |
| total exchanger working-gas volume | 158 cm³ | 243 cm³ |

## Run

```bash
research validate outputs/pedal_cell_retro/long_thin_f013_100x/study.toml

research evaluate \
  outputs/pedal_cell_retro/long_thin_f013_100x/study.toml \
  --output outputs/pedal_cell_retro/long_thin_f013_100x/center.json

research run \
  outputs/pedal_cell_retro/long_thin_f013_100x/study.toml \
  --directory outputs/pedal_cell_retro/long_thin_f013_100x/campaign \
  --budget 6h \
  --max-candidates 1024
```

The center evaluation is useful here because the exchanger geometry is a large
discontinuous design change from f013; it immediately tells whether the chosen
1 mm / ~0.25-0.33 m starting point is within the microtube hydraulic/model domain.
