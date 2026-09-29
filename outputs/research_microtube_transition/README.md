# Microtube transition validation

One bounded evaluation, without optimization, uses the current local Human Cell
Stage 2 study and overrides both banks to N=1000, D=0.000770 m, L=0.8 m.
The complete scientific definition and candidate values are embedded in
[D_0.000770.json](D_0.000770.json).

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m dada_solver.research.cli evaluate \
  outputs/human_cell_stage2_microtube/human_cell_stage2_microtube.toml \
  --set microtube.heat_in.tube_count=1000 --set microtube.heat_out.tube_count=1000 \
  --set microtube.heat_in.tube_length_m=0.8 --set microtube.heat_out.tube_length_m=0.8 \
  --set microtube.heat_in.inner_diameter_m=0.00077 --set microtube.heat_out.inner_diameter_m=0.00077 \
  --budget 30s --output /tmp/microtube-transition-reproduction.json
```

The first cycle passes the old unavailable-transition barrier and then rejects
`large_relative_pressure_drop`, without a Reynolds-domain rejection:

- Heat-out, `large_to_hot`, angle 0.3103138843 rad (0.01763856735 s).
- Re=2393.29848, `gnielinski_transition_interpolation`.
- Relative pressure drop=0.2005926643, above the unchanged 0.2 guard.
- Mass flow=0.02698352741 kg/s, Pr=0.710860217, Ma=0.146924741,
  Kn=0.00009091245, gas temperature=303.685025 K.
- Numba requested/selected; 40 explicit Python fallbacks in 412 RHS calls.
- Evaluation duration 8.02 s, including backend preparation, on this machine.

This is not a converged or feasible refrigerator result. It confirms transition
availability while retaining the pressure-drop rejection. Cycle-level transition
fractions are only available for a completed diagnostic trajectory; they are not
inferred from this first-cycle failure.

Validation: 146 directly related Research/exchanger tests passed. The final
complete suite passed **815 tests, 0 warnings in 295.04 s** with
`PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 -m pytest -q`.
