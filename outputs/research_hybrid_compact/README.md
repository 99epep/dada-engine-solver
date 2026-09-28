# Compact hybrid Research replay

The eleventh family is `hybrid_compact`, with nine independently fixed/active
coordinates (five small-side, four large-side). Production implementation:
`src/dada_solver/hybrid_compact_kinematics.py`.

Reproduce from a source checkout, choosing unused output paths:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m dada_solver.research init kinematics --small hybrid_compact --large hybrid_compact --output outputs/my_compact/study.toml
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 -m dada_solver.research evaluate outputs/my_compact/study.toml --output outputs/my_compact/champion.json --budget 2m
```

The paired preset includes the historical candidate 501 machine and preceding
warm state. `champion.json` is the bounded Research evaluation, `champion.html`
is its offline report, and `parity.json` records the observed errors.

Observed parity:

- Scalar and vector volumes/first derivatives: zero difference against two
  frozen pre-extraction cases, each sampled at 4097 angles over three cycles.
- Indicated power: 40.555304211350325 W.
- Indicated efficiency: 0.22525540307940708.
- External heat input: 180.04142700654225 W.
- Periodic convergence: 7 cycles, as in the original record.
- All ten compared physical metrics: zero difference on this runtime; tests
  permit 2e-11 relative / 1e-11 absolute numerical tolerance.

Phase, sign, scalar/vector evaluation, original parameter-domain guards and
sampled monotonicity screens are preserved. The example keeps its historical
search/fitting logic and reexports the production implementation. Analytic
second derivatives remain unavailable, as in the original class. The law is
C2 but no numerical acceleration is substituted for an analytic API.

This is a fixed-input replay, not a new optimization or a physical hardware
validation. Mixed-family and V3 external-stream configuration are supported;
they require their own scientific evaluation.

Full-suite validation: **744 passed, 0 warnings in 293.79 s** with
`PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 -m pytest -q`.
