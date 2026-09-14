# Independent six-bar K2 evaluation

The selected S mechanism is Stage 2F restart 0, secondary branch +1 from
`outputs/small_sixbar_stage2f_r4_freeH_tightaxis.json`. L is the global best from
`outputs/large_sixbar_stageL1.json`, used without an additional mirror or phase.
The full reproducible record is `outputs/motor_champion_sixbar_k2.json`.

With the unchanged K2 thermodynamics (25/325 deg C sources, 2 Hz, fixed gas
inventory of 1.497823285 g, asymmetric microtubes and dynamic walls), the solve
converged in 23 cycles and 38.21 seconds. The last normalized periodic error was
0.967057, below the unchanged threshold of one. The existing bounded wall
initial-guess acceleration was used; convergence was established by an ordinary
physical cycle.

| Quantity | Independent six-bars | Free K2 reference |
| --- | ---: | ---: |
| Indicated thermal efficiency | 20.381817% | 20.561604% |
| Indicated gas power | 51.406522 W | 51.312591 W |
| External-source heat input | 252.217564 W | 249.555386 W |

The efficiency difference is -0.179787 percentage point. Mechanical losses
remain unknown and useful mechanical output is unavailable. These numbers
compare thermodynamics, not mechanism friction, strength or practical packaging.
No physical crank radius has been inferred.

| Gas region | Pressure range, absolute bar | Temperature range, deg C |
| --- | ---: | ---: |
| S | 1.20074–4.98232 | 54.945–212.775 |
| L | 1.22532–4.95407 | 106.524–311.808 |
| Hi (legacy C) | 1.20074–4.96951 | 294.839–309.860 |
| Ho (legacy H) | 1.21413–4.95405 | 41.965–63.230 |

| Passage, positive reference direction | Mass-flow range, g/s |
| --- | ---: |
| S to Hi | -0.90048–10.73758 |
| Hi to L | 0–10.48504 |
| L to Ho | -1.84482–5.88944 |
| Ho to S | 0–5.86615 |

Local reverse flows in the non-valved passages are retained. Maximum tube
Reynolds number is 558.64 and maximum tube Mach number is 0.03867; the existing
thermodynamic validity assessment is valid with no failed/unavailable criteria.
The wall-cycle diagnostic adapter does not reconstruct valve events: its generic
`non_nominal` event-sequence label is not evidence of an observed chronology.
Wall-to-gas heat extrema are separate from the external-source efficiency boundary.

The motion RMS values against the exported target CSV are 0.538276% for S and
0.370029% for L. Differences from synthesis summary rounding also reflect the
CSV endpoint convention and continuous-extremum normalization. The raw position
and derivative regression independently reproduces the synthesis implementation.

Validation: 8 focused six-bar tests pass, covering synthesis equivalence,
analytic derivatives, periodic wrap continuity, volume bounds, candidate
selection, invalid geometry and exactly one motor reversal. The complete suite
passes: 319 tests in 44.04 seconds. The K2 calculation above is the end-to-end
dynamic-wall integration check; no thermodynamic equations were modified.
