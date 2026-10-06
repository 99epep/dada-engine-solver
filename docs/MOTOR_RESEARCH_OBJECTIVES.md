# Motor research objectives

These are ordered research questions, not permission to launch an unattended
campaign or universal defaults. New studies own their power, size,
source-temperature and frequency requirements explicitly.

## 1. Improve motion laws, then realize them mechanically

Search small/large timing and shape independently under explicit constraints.
Mirror symmetry is optional, and passive valve events follow pressures rather
than prescribed phase labels. Piecewise-linear laws are probes with velocity
jumps; smooth mathematical laws can also contain acceleration artifacts.
A good thermodynamic target is not automatically a buildable mechanism.

Use the established progression: abstract motion → primary family discovery →
downstream fitting → full local mechanism adaptation → paired thermodynamic
replay/optimization → hardware retuning. Preserve several mechanical tradeoffs,
not only the lowest curve-fit error. Position is the primary motion reference;
velocity describes cadence/topology; do not blindly fit sharp acceleration.
See [synthesis method](MECHANISM_SYNTHESIS_SEARCH.md),
[family catalogs](SIX_BAR_MECHANISM_FAMILIES.md) and
[Research family interfaces](DADA_ENGINE_RESEARCH_KINEMATICS.md).

## 2. Map Lambdas and cylinder swept-volume ratio

Study `Lambda_S`, `Lambda_L` and `V_swept_S/V_swept_L` against source temperatures
and required output. Report clearances separately. Distinguish configured
geometric targets from Lambdas measured at actual passive-valve events; missing
or repeated events may make a single Lambda unavailable.

Condition any map on fluid, speed, charge, exchanger geometry, model domains,
mechanical constraints and power boundary. Temperature alone need not select a
unique best design. Keep frozen-design sensitivity separate from reoptimized
hardware/motion, and report budgets, seeds and retained basins. A finite local
search is best-found evidence, not a universal design law.

## New-study boundary

For a new temperature/geometry campaign, define its own filling policy. The
reference-pressure policy uses the maximum **simultaneous total gas volume** of
the assembled machine. Keep hardware-derived hold-up, conductance and capacity
coupled. External stream inputs and excluded auxiliary losses remain explicit.

Slider-crank, four-bar and six-bar are the immediate mechanism classes. Solenoid
actuation and complete mechanical-loss/inertia models remain separate work.
Report indicated gas power independently from unknown useful shaft power.
