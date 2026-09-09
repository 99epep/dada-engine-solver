# Motor research objectives

These are user-requested long-term objectives, ordered after the current
exchanger characterization and demonstrator work. They are not completed
optimization results or an instruction to launch an unrestricted search now.

## 1. Search for optimal motion laws

Maximize efficiency subject to specified useful power, frequency, size and
hardware applicability constraints. Search the small- and large-cylinder laws
independently: timing, dwell, transfer duration and waveform shape need not be
mirrored. The current piecewise-linear law is a reproducible reference, not an
established ideal, a performance ceiling or a required optimum. The desired
cycle chronology is distinct from any particular piston motion law.

Explore smooth periodic laws with explicit stroke, velocity, acceleration and
mechanical feasibility bounds. Report which bounds and losses are modeled;
without them an apparent optimum could require unrealizable motion. Keep
valves passive and measure the resulting chronology rather than schedule it.

Compare symmetric and independently parameterized four-bar mechanisms under
the same thermal, hydraulic, inventory and size assumptions. Mirror symmetry
was a refrigerator construction choice, not a motor requirement. The existing
shared-crank configuration already supports distinct loop parameters; retain
common-crank geometric compatibility where a shared crank is used. Independent
loop parameters do not imply two independent crank speeds. Do not constrain
all candidate waveforms to what the present four-bar seed can reproduce.

## 2. Map Lambdas and cylinder volume ratio

Immediately after the motion-law objective, seek design maps for `Lambda_S`,
`Lambda_L` and `V_swept_S / V_swept_L` as functions of the two reservoir
temperatures and required useful output power. Report clearances separately;
do not conflate swept-volume ratio with maximum enclosed-volume ratio.

Distinguish configured geometric Lambda targets from Lambdas inferred from
actual passive-valve events. Define how each is extracted for a general smooth
law. If repeated or missing events make a single Lambda ambiguous, report that
ambiguity and the event chronology rather than forcing a nominal value.

Temperature and power alone need not determine a unique best design. Condition
maps on frequency, working fluid, charge, exchanger hardware, pressure losses,
mechanical limits and the adopted useful-power boundary. Compare feasible
tradeoffs and report sensitivity to uncertain exchanger properties. Do not
present a sparse local search as a universal law or a global optimum.

Expected deliverables are reproducible configuration families, feasible ranges
and best-found regions, with efficiency, size and uncertainty alongside the
Lambdas and volume ratio. The present 25/325-degree, approximately 100 W,
2–10 Hz demonstrator is an initial application, not the full map domain.
