# Frozen compact hybrid parity reference

Captured from the local historical implementation before extraction:
`examples/optimize_motor_hybrid_compact_260k.py`, SHA-256
`f9906b877a0b97fddcd13b4dbc5f901eefd8944ba7ecfc3e19ec6e44d0ab2add`.

`trajectories.npz` stores 4097 study angles from -2*pi through 4*pi, scalar
`cylinder_volumes_and_derivatives` results and independently evaluated vector
results. Columns are small volume, large volume, small dV/dtheta, large
dV/dtheta. Case 0 is historical champion 501; case 1 tests symmetric kink
coordinates with the minimum and maximum admitted rounding fractions.

`reference.json` records both parameter sets, the original volume limits,
champion record and complete historical thermal basis. The thermal basis was
constructed using the original Fourier-source machine and fixed source16
inventory, following the historical search script. Its initial wall state is
from the champion's recorded nearest warm-start candidate, not its final state.
This is necessary to compare the periodic count and stored result exactly.

The original implementation generated the frozen arrays; production code did
not generate its own expected values. The current example reexports the
production class and cannot serve as an independent motion oracle anymore.
Tests compare against these frozen data and the previously stored champion
metrics. No historical output/history is changed.
