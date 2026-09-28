# Human Cell Stage 0

Initial Dada-Engine Research files for the first coupled human-powered
refrigeration study.

Files:

- `study.toml` — Research V3 study definition.
- `STUDY_NOTES.md` — rationale, known limitation of `four_stage` ownership,
  and proposed progression.

Expected location in the source checkout:

```text
outputs/human_cell_stage0/
```

`study.toml` intentionally reuses:

```text
outputs/research_v3/refrigeration.basis.json
```

and its recorded SHA-256.  If that validated V3 basis is absent or differs in
the local checkout, regenerate the V3 refrigeration preset/basis first instead
of changing the hash by hand.

Before running a campaign:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src \
  python3 -m dada_solver.research validate outputs/human_cell_stage0/study.toml
```

The study maximizes cooling COP.  There is deliberately no mechanical-input
power target or minimum cooling-power constraint in Stage 0.
