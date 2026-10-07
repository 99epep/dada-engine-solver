# Frozen compact hybrid parity reference

The independent reference implementation has SHA-256
`f9906b877a0b97fddcd13b4dbc5f901eefd8944ba7ecfc3e19ec6e44d0ab2add`.

`trajectories.npz` stores 4097 study angles from -2*pi through 4*pi, scalar
`cylinder_volumes_and_derivatives` results and independently evaluated vector
results. Columns are small volume, large volume, small dV/dtheta, large
dV/dtheta. Case 0 exercises asymmetric motion; case 1 tests symmetric kink
coordinates with the minimum and maximum admitted rounding fractions.

`reference.json` supplies the parameter sets and volume limits. The frozen arrays
were generated independently of the production class. Tests use these data as a
motion oracle, not as a current thermal regression reference.
