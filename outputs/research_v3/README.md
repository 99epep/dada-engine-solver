# Research V3 bounded validation artifacts

No large optimization or repository publication was performed.

- `validated/`: final-source evaluations, two-candidate persisted Sobol resume,
  structured-C2 setup validation, scientific inputs, and offline reports.
- `backend_benchmark.json`: fixed-state ideal/table/Python timing, separate first
  calls and steady calls, source/runtime identities and numerical parity errors.
- `pre_v3_ideal_timing.json` / `v3_ideal_timing.json`: five repetitions of 30,000
  ideal compiled calls on the same reference state. The former uses the local
  pre-V3 committed source exported to an isolated temporary directory.
- `acceptance/`: earlier complete acceptance run, retained unchanged. Its source
  identity precedes the final diagnostic-reference adjustment; use `validated/`
  for current execution and resume. Offline inspection remains supported.
- `refrigeration.first.json`: rejected initial high-compression fixture, which
  left the 200 K transport domain. No limit was loosened.
- `refrigeration.clearance.json`: first successful explicit 30% clearance,
  0.2 Hz refrigerator evaluation. Its study input is `refrigeration.toml`.

The declared liquid Cp (4180 J/(kg K)), flow and conductance are scenario inputs,
not a validated water-loop correlation. External hydraulic/pump losses and shaft
losses are not supplied. Table data are generated from the analytic ideal gas;
there is no real-helium or cryogenic validation claim. Different families or
operating modes do not establish a fair optimization comparison.

The complete test command was:

```sh
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONPATH=src python3 -m pytest -q
```

Result: **728 passed, 0 warnings**, 256.47 s on this checkout. The pre-V3 baseline
was 679 passing tests. Existing tests were retained; the frozen pre-V3 source
trajectory is `tests/data/pre_v3_air_wall.json` with source commit/blob hashes.
