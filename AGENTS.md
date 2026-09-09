# Project conventions

All project content must be written in English: documentation, comments,
docstrings, identifiers, messages, example descriptions, reports and figures.
Standard scientific symbols, SI units, proper names and original source URLs
are preserved. The conversation with the user may remain in French.

Read `docs/PHYSICS_DECISIONS.md` before changing the physical model. Current user
decisions override historical project decisions. Motor efficiency is the primary
objective under explicit power and size constraints; ideal chronology remains
the target. Isothermality is diagnostic-only for the current motor work.

The motor is a small experimental demonstrator, not a commercial product.
Keep its name independent of the changeable reservoir temperature difference.
The current motor research envelope is 2–10 Hz, with reservoirs at 298.15 K and
598.15 K (25/325 deg C); reduce the hot-source temperature progressively after
improving dimensions. The earlier 448.15 K hot source is historical. Target approximately 100 W useful mechanical output and limit the
large cylinder maximum enclosed volume (including clearance) to 0.066 m^3.
This is a ceiling, not a target displacement. Report indicated gas power
separately from useful output; mechanical losses remain unknown; the opt-in hardware model estimates fan power.
Prefer a documented commercial exchanger and adapt the machine to it.
If suitable commercial data are unavailable, the user authorizes extrapolation
from Doty's measurements with explicit assumptions and uncertainty scenarios.
The intended exchanger through-flow is intermittent and unidirectional; report
any modeled local reflux without hiding it. Explain the purpose and cost before
adding a dynamic exchanger model.

Run relevant checks from a source checkout using `PYTHONPATH=src python3 -m pytest`.

The piecewise-linear motion law is a reference, not a proven optimum. The motor
search must not inherit mandatory mirror symmetry from the refrigerator seed.
Record the ordered long-term objectives in docs/MOTOR_RESEARCH_OBJECTIVES.md:
first optimize motion laws, then map independent Lambdas and cylinder swept-
volume ratio against reservoir temperatures and useful-power requirements.

External air is the selected hot/cold source medium for the demonstrator.
Treat source temperatures as air inlet conditions, with finite flow. The latest
user decision excludes external-air aerodynamic losses and fan consumption from
the current trial balance; do not confuse exclusion with zero physical losses. Retain working-gas heat exchange during pauses; independent
wall storage is implemented in the opt-in air_wall wrapper, but not calibrated
to Doty or exposed through the historical CLI.

For current progress and reproducible graphs, read docs/PARALLEL_EXCHANGER_TRIAL.md,
docs/HIGHER_TEMPERATURE_TRIAL.md and docs/EXCESS_AIR_TRIAL.md
and the preceding docs/DOUBLED_EXCHANGER_TRIAL.md. External airflow is now sized
from expected peak internal flow with a capacity-rate margin; verify achieved
ratios and retain finite thermal-film resistance.
Use checkpoint intervals and optional bounded wall initial-guess acceleration
without weakening periodic convergence or changing states within a cycle.

Architecture handoff: read docs/PLUGGABLE_MODELS.md before implementing the
global optimization campaign. Preserve generic KinematicsModel and
ExchangerModel boundaries, independent free spline motions, explicit validity
margins, and family-specific parameter ownership. Do not independently optimize
geometrically derived exchanger UA, hold-up or losses. The current task adds
no optimizer, dynamic plugin registry, or assumed mechanical efficiency.
