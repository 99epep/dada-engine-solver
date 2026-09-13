# Motor champion four-bar validation

The compact mechanism synthesized by `dada-4bar-synthesis` is connected to the
motor through `SharedCrankCouplerProjectionKinematics`.  Both cylinder mechanisms
share the same crank angle.  Each coupler-point projection is normalized to a
zero-to-one travel coordinate and mapped to the established cylinder minimum and
maximum volumes.

The synthesis export uses the study-angle convention.  It is injected into
`MachineDesign` before operation-direction handling, so the motor transformation
`V(-theta)` and its derivative sign change still occur exactly once in
`build_model`.

`examples/evaluate_motor_champion_four_bar_k2.py` evaluates this mechanism with
the established Stage K2 thermodynamic design.  The run preserves:

* the Stage A5/F6 working-gas inventory and cylinder volume ranges;
* the 25 °C / 325 °C source temperatures and 2 Hz operating point;
* K2's independently optimized H_i and H_o microtube lengths;
* the ten-state dynamic-wall model and external-source efficiency boundary;
* the existing valve, exchanger-correlation, and periodic-convergence models.

Only the free Stage F6 motion law is replaced.  The generated JSON records both
the motion-fit error and the periodic thermodynamic result.  This is a validation
of the projected mechanism, not a claim that its dimensional link sizes, bearing
loads, balancing, or mechanical losses have been resolved.  Because the geometry
is normalized by crank radius, selecting a physical crank radius remains a later
mechanical design decision.
