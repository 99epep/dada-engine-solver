# 5 L / 2.8 Hz — DD compact-microtube COP campaign

Purpose: deliberately push exchanger working-gas volume far below the previous campaign before worrying about useful cooling power.

## Search design

- Frequency fixed at **2.8 Hz**; swept volume fixed at **5 L**.
- Valve topology fixed **DD**.
- `volume.swept_ratio` remains free over **0.5..1.6**; values above 1 are intentionally allowed.
- The 9 `hybrid_compact` kinematic coordinates remain free.
- Per exchanger:
  - tube count: **20,000..200,000**
  - tube length: **6..20 mm**, logarithmic
  - inner diameter: **80..180 µm**, logarithmic
  - wall thickness: **20 µm** fixed
  - pitch: **225 µm** fixed
  - header depth: **1 mm** fixed
- Three previous DD kinematic basins are retained as local Sobol centers; all are recentered on the same compact exchanger geometry.
- Search radius: **0.45**.

The nominal center is:

- H_i: 100,000 × 10 mm × 120 µm → **21.45 mL** working-gas volume, 11.31 cm² total flow area, 0.377 m² internal surface.
- H_o: 120,000 × 12 mm × 120 µm → **28.44 mL**, 13.57 cm², 0.543 m².
- Combined exchanger volume: **49.89 mL = 0.998% of the 5 L swept volume**.

Even the simultaneous upper geometry corner is about **244.1 mL = 4.88%** of swept volume, so this campaign cannot drift back into the previous very-large-exchanger regime.

## Run

Place this directory under:

`outputs/pedal_cell_retro/5L_2p8Hz_DD_compact_microtubes/`

then run:

```bash
PYTHONPATH=src python3 -m dada_solver.research.cli run \
  outputs/pedal_cell_retro/5L_2p8Hz_DD_compact_microtubes/study.toml \
  --directory outputs/pedal_cell_retro/5L_2p8Hz_DD_compact_microtubes/campaign \
  --budget 6h --max-candidates 4096
```

The copied `study.basis.json` retains the air `dilute_species_v2` transport domain down to 100 K. No valve-placement optimization is present in this pass.
