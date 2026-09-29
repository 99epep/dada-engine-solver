# Human Cell Stage 2: reference-pressure filling

This is a copy of the local
`outputs/human_cell_stage2_microtube/human_cell_stage2_microtube.toml` seed.
The original study and basis were not modified. The basis is copied byte-for-byte;
its historical inventory is superseded by the declared derived-charge policy.
The copy removes `charge.total_mass_kg`, sets reference pressure 100000 Pa and
reference temperature 293.15 K, and uses uniform initialization. All candidate
geometry, motion, external streams, constraints and convergence settings remain
those of the source seed. This is not the earlier 1000-tube diameter sweep.

Artifacts:

- [Study](study.toml), [portable basis](study.basis.json).
- [Exact evaluation](seed.json), [standalone report](report.html).

Reproduce one bounded evaluation into a new output:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m dada_solver.research.cli evaluate \
  outputs/human_cell_stage2_reference_charge/study.toml \
  --budget 60s --output /tmp/human-cell-reference-charge.json
```

The original evaluation converged after 13 cycles, in 18.09 s on this machine.
It is **converged_infeasible**, and its operating mode is **non_refrigeration**.
No optimization or subsequent candidate adjustment was performed.

| Quantity | Result |
|---|---:|
| Maximum simultaneous total gas volume | 0.007298096764755425 m³ |
| Reference solver angle | 0.012630857211916823 rad (0.72369481 degrees) |
| Derived inventory | 0.008674367589288079 kg |
| Reconstructed reference pressure, m R T / V | 100000.00000000003 Pa |
| Indicated mechanical input | 72.48982878 W |
| Signed cooling power | -1.00842009 W |
| Heat rejected at hot boundary | 71.47487766 W |
| Cooling COP | Unavailable: no positive cooling load |
| Maximum pressure / limit | 184141.04 / 1200000 Pa — satisfied |
| Maximum temperature / limit | 328.20023 / 850 K — satisfied |
| Maximum absolute mass flow / limit | 0.11880496 / 0.08 kg/s — violated |
| Large enclosed volume / limit | 0.00518531 / 0.066 m³ — satisfied |
| Thermodynamic model validity | valid |
| Microtube model domain | satisfied; transition uncertainty remains explicit |

The negative cold-boundary load must not be presented as useful cooling. The
indicated input is not a measured shaft/human power requirement. Existing model
and external-loop limitations are unchanged.

Validation: **256 Research/campaign/interface/coupling tests passed, 0 warnings
in 265.00 s**, including 29 new reference-charge cases. No complete unrelated
suite or optimization campaign was launched for this change.
