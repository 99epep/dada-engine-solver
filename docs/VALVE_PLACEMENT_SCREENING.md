# Passive check-valve placement screening

Completed 2026-09-21 on the frozen 300 K-difference motor champion. This is a
four-case evaluation of one candidate, not an optimization.

## Physical graph and naming

The nominal circulation is `S -> H_i -> L -> H_o -> S`. Legacy solver names use
`C` (state index 2) for `H_i` and `H` (index 3) for `H_o`; `S` and `L` are indices
0 and 1. The four topology codes give the H_i placement first and H_o placement
second: `D` means downstream and `U` upstream.

| Code | H_i branch | H_o branch |
|---|---|---|
| DD | `S <-> H_i -> L` | `L <-> H_o -> S` |
| UD | `S -> H_i <-> L` | `L <-> H_o -> S` |
| DU | `S <-> H_i -> L` | `L -> H_o <-> S` |
| UU | `S -> H_i <-> L` | `L -> H_o <-> S` |

Placement changes the physical hydraulic network: the selected half-link becomes
unidirectional and receives the configured valve CdA, while the other half-link
becomes bidirectional and has no valve CdA. Valve opening remains strictly passive
and pressure-driven. No angle or motion-phase command is used.

`ValveTopology.hot_to_small` continues to report the H_o valve state even when
that valve is physically on `large_to_hot`. Likewise, `cold_to_large` reports the
H_i valve state even when it is on `small_to_cold`. This preserves the public
legacy diagnostic API. Alternate placements are implemented for
`continuous_ideal_diode`; `discrete_hysteretic` accepts DD only.

## Method

The script [`screen_motor_valve_placements_300k.py`](../examples/screen_motor_valve_placements_300k.py)
loads the champion and its ten-state `last_complete_state` from
`outputs/motor_temperature_map/dT_300K/report.json`. It reconstructs the exact
saved geometry, exchanger counts and lengths, motion parameters, 298.15/598.15 K
sources and 100 kPa at 298.15 K charge. That same saved DD endpoint is only an
initial guess for every topology; each case is integrated to the unchanged
periodic criterion.

The production Numba backend, exact-angle cache, LSODA, tolerances, breakpoints,
microtube correlations and full Stage-5 diagnostics are unchanged. Every case
reports `requested_backend = numba`, `actual_backend = numba`, and zero fallback.
Elapsed time includes full candidate diagnostics; DD also includes first JIT
compilation in this serial run, so these one-shot times are descriptive rather
than a backend timing comparison.

## Results

All cases satisfy the 300 K campaign limits: at least 40 W, at most 1.2 MPa,
850 K and 0.08 kg/s, with a valid thermodynamic/microtube verdict.

| Topology | Cycles | Time (s) | Efficiency | Power (W) | Heat in (W) | Heat out (W) | Max pressure (MPa) | Max temperature (K) | Max flow (kg/s) | Max Re | Max Mach | Feasible |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| DD | 1 | 8.168 | 22.0829% | 52.4066 | 237.3172 | -184.9445 | 0.500601 | 584.405 | 0.0100723 | 461.96 | 0.028574 | yes |
| UD | 14 | 5.865 | 21.4509% | 55.4089 | 258.3060 | -202.8156 | 0.490511 | 582.461 | 0.0127062 | 668.28 | 0.028465 | yes |
| DU | 13 | 3.939 | 22.4687% | 54.0578 | 240.5920 | -186.5744 | 0.498162 | 584.057 | 0.0142067 | 750.33 | 0.028560 | yes |
| UU | 13 | 3.902 | 22.7169% | 56.7675 | 249.8906 | -193.1963 | 0.490392 | 585.479 | 0.0141962 | 748.35 | 0.028523 | yes |

The DD replay starts from the saved periodic endpoint and converges after one
ordinary cycle. It gives 22.0829326% and 52.406601 W, versus 22.0833993% and
52.406880 W in the historical report. The small differences (4.67e-6 absolute
efficiency and 2.78e-4 W) arise from advancing the tolerance-level saved endpoint
through another complete adaptive cycle; the default DD algebra and Python/Numba
RHS regression tests remain unchanged.

UU has both the highest raw efficiency and the highest efficiency among feasible
cases in this screen. It is therefore the topology retained for the next 300 K
campaign. This result applies only to the frozen candidate: no geometry, motion,
charge or valve CdA was reoptimized after moving the valves.

The complete machine-readable record is
`outputs/motor_valve_placement_screen_300K.json`.
