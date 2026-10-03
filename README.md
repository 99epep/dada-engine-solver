# DADA Engine Solver

DADA simulates a thermal machine with conservative gas balances, passive valves,
prescribed piston motion and interchangeable exchangers. **Research** configures
studies, evaluates candidates, runs resumable Sobol searches and compares results
in standalone HTML reports. It supports motor and refrigeration studies,
independent motion families, local multi-centre search and categorical valve
placement.

## Start here

From the local source checkout:

```sh
research() { PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m dada_solver.research "$@"; }
research --help
```

- [Research guide](docs/DADA_ENGINE_RESEARCH.md): study setup, run/resume,
  refine/rescale, reports, curves and mechanism animations.
- [Documentation index](docs/README.md): physics, models, development and historical evidence.
- [Current physical decisions](docs/PHYSICS_DECISIONS.md): conventions, validity,
  model boundaries and deliberate omissions.

The local study, portable basis, source version and result artifacts define a
calculation. Historical examples are reproducibility cases, not universal design
requirements. No useful shaft power is inferred from indicated gas work.

## Install and verify

Python 3.11 or newer is required. The runtime dependencies are NumPy and SciPy.
From a virtual environment, install the checkout and optional test/compiled/plot
support as needed:

```sh
python3 -m pip install -e '.[test,numba,plot]'
PYTHONPATH=src python3 -m pytest -q
```

Numba and plotting dependencies are optional. Research can also run with the
source helper above. The installed command is `dada-research`; this documentation
uses the shorter `research` helper. For a direct simulation TOML, use
`dada-solver configuration.toml`; add `--integration-profile` for solver counters.
The lower-level sizing and campaign interfaces remain available for existing
scripts; see [architecture](docs/PLUGGABLE_MODELS.md).

[Validation](docs/validation.md) separates analytical checks, numerical parity,
correlation evidence and experimental limitations. Passing tests or reaching a
periodic state does not establish that a physical machine is validated.

## Scientific and implementation provenance

The originating [Thermodynamic and Mechanical Study](https://dada-engine.org/Thermodynamic_and_Mechanical_Study)
and its locally recorded revisions are discussed in the
[decision ledger](docs/history/PHYSICS_DECISION_LEDGER.md). Current implementation
choices and later user decisions are summarized in the current physics reference.

> The substantial initial implementation was produced by OpenAI Codex running
> GPT-5.6 Sol (`gpt-5.6-sol`), under human scientific direction, on 2026-09-02
> (initial project version 0.1.0). Later work continues to require independent
> review of equations, conventions, assumptions and results. AI-generated
> scientific software can contain plausible but serious errors.

## License

Copyright (C) 2026 DADA Engine Solver contributors.
Licensed under **GPL-3.0-or-later**; see [LICENSE](LICENSE). There is no warranty,
to the extent permitted by law. Adapted transport formulas and coefficients are
attributed in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Third-party
materials retain their respective terms.
