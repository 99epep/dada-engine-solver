# Dada Engine Research user guide

Research configures, evaluates, searches and compares machines using the existing
physical solver. Indicated gas power is not automatically useful mechanical
shaft power: mechanical losses remain unknown.

From the repository root, define the shorthand used throughout this guide:

```sh
research() { PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m dada_solver.research "$@"; }
```

The installed command is `dada-research`; use it instead of `research` outside a
source checkout. All paths below are illustrative destinations for new work.
Use `research COMMAND --help` for the complete command options.

## Create a study

Choose one of the current starting templates:

```sh
research init kinematics --small slider_crank --large harmonic \
    --output outputs/my_study/study.toml
research init external-stream-refrigeration \
    --output outputs/my_refrigerator/study.toml
research init external-stream-motor \
    --output outputs/my_motor/study.toml
```

Templates start with fixed coordinates. Edit the generated study for the intended
machine and operating conditions. A fixed parameter has a `value`; to search a
coordinate, use the study editor below or an active declaration with its initial
value and bounds. Do this before `run`, which needs at least one active parameter.

See the [technical reference](DADA_ENGINE_RESEARCH_REFERENCE.md) for parameter
types, policies, objectives and constraints, and the
[kinematics reference](DADA_ENGINE_RESEARCH_KINEMATICS.md) for family selection.

## Edit a study

```sh
pip install "dada-engine-solver[tui]"
research study edit path/to/study.toml
research study release path/to/retuned/study.toml --group mechanisms \
    --radius 0.05 --output path/to/combined/study.toml
research study freeze path/to/study.toml --group frequency \
    --output path/to/fixed-frequency/study.toml
```

Use arrows, Space, Enter and Escape to inspect groups, toggle activity, edit
individual parameters, review and save to a new path. Partial activity is marked
`[~]`. Coordinates stored in mechanical artifacts are available just like explicit
parameters. Reopening SMALL/LARGE mechanisms preserves the current machine and
already active hardware; categories and branches stay fixed. Cancel writes nothing.

The declared parameter bounds describe a numerical search domain. A local radius
is a separate scheduler control: it restricts **all active parameters**, including
hardware, in their normalized declared domains. Batch release defaults to local
Sobol at radius 0.05; use `--search global` to search full declared bounds instead.
Missing domains require explicit bounds; values are never clamped into a domain.
Review shows effective bounds and any necessary recentring or artifact regeneration.
To balance search windows by group, use **Recalibrate bounds for this group**
or `research study bounds recenter path/to/active/study.toml --group mechanisms
--half-width 0.01 --output path/to/narrow/study.toml`. It narrows active numeric
mechanism declarations around their exact initials, leaving exchanger bounds and
the scheduler radius unchanged. The TUI previews old/new intervals and offers
undo; setting the global local radius to 1 is a separate confirmed action.
Toggles stay in the group submenu with the cursor on the same parameter.

Saving validates and copies portable inputs without integrating or changing the
source. Then use `validate` and `run` explicitly; use `evaluate` for a fully fixed
study. Batch commands require no TUI extra. Full options and domain recovery rules
are in the [study editor reference](DADA_ENGINE_RESEARCH_REFERENCE.md#study-editor).

## Validate and evaluate one point

```sh
research validate outputs/my_study/study.toml
research validate outputs/my_study/study.toml --search-window
research evaluate outputs/my_study/study.toml \
    --output outputs/my_study/reference.json --budget 3m \
    --set microtube.heat_in.tube_count=24   # override a parameter (repeatable)
```

The `--set` flag (alias `--set`) can be repeated to override any parameter defined in the study TOML. Each assignment must be `NAME=VALUE`; the CLI validates the name, casts to the correct type (integer, float, or raw string for choices), and uses the overridden value for the evaluation, leaving the original study file unchanged. This is handy for quick sensitivity checks without editing the TOML.

`validate` checks declarations, hashes and production geometry without integrating thermodynamics. `--search-window` shows active parameters by group, their initials, declared domains, and the effective window of each local region. It never overwrites an existing result, so choose a new output filename for each evaluation.

`validate` checks declarations, hashes and production geometry without integrating
thermodynamics. `--search-window` shows active parameters by group, their
initials, declared domains, and the effective window of each local region. It
shows the global normalized radius and clipping, respects log/integer encodings,
and reports eligible choices without treating categories as physical distances.
The grouped table rounds numbers to at most five significant digits without modifying
values. Numeric bound pairs stay on one line, with fewer digits when necessary. Other
long cells wrap within their columns; region identities label each table. The same display is available in `study edit`; add `--json` for structured
output with full precision.
Invalid declarations report the source file, line, parameter or
constraint, and the corresponding TOML excerpt. Coupled geometry failures identify
the affected piston and its related parameters rather than blaming an arbitrary
coordinate. `evaluate` integrates one exact configured point, using active
initial values where applicable; it does not run a search. It never overwrites an
existing result, so choose a new output filename for another evaluation.

## Run, resume and inspect

After declaring active parameters:

```sh
research run outputs/my_study/study.toml --budget 30m --max-candidates 32
research resume outputs/my_study/campaign --budget 2h --max-candidates 128
research status outputs/my_study/campaign
```

`run` creates `campaign` beside the study TOML by default. Use `--directory` to
choose another destination. A nonempty destination is refused; use `resume` to
continue it. Budget and `--max-candidates` bound the current invocation.

`resume` uses the campaign's stored snapshots, not later edits to the original
TOML. It also accepts a TOML path as shorthand for its sibling `campaign` directory.
Resume checks execution compatibility; a missing campaign does not start a new
search. `status` inspects stored results without integrating thermodynamics.

For persistence and compatibility rules, see the
[technical reference](DADA_ENGINE_RESEARCH_REFERENCE.md). Progress display,
interruption recovery and result evidence are covered in the
[cockpit guide](DADA_ENGINE_RESEARCH_COCKPIT.md).

## Reports, curves and animations

```sh
research report outputs/my_study/campaign
research report outputs/my_study/campaign --candidate second
research report outputs/my_study/campaign --candidate best --plots default,flows,mechanisms
research report outputs/my_study/campaign --plots none
research compare outputs/my_study/campaign outputs/my_motor/campaign \
    --html outputs/comparison.html
```

By default, `report` writes `campaign/report.html` and requests positions, heat
exchange, pressures and temperatures for the two best feasible candidates.
`--candidate` also accepts an exact ID or unique prefix; repeat it to select
several candidates. `--plots none` omits curves, while `mechanisms` requests
a whole-machine animated WebP when both cylinders have complete physical
mechanisms. The animation uses production linkage closures, replayed gas
temperatures and effective valve states, with the candidate's actual valve
placements. Its independently scaled visible strokes illustrate motion, not a
physical common-shaft assembly; animation speed is illustrative. The validated
transfer block is uniformly scaled between fixed cylinder inner faces.
The large cylinder's transverse display width is at least 20% of the exported
image height, including margins, so large linkages with short piston strokes
do not collapse the chamber and transfer rows visually. Both cylinders enlarge
by the same factor, preserving their relative transverse dimensions. This changes only the
illustration; piston strokes and scientific volumes remain unchanged.
Pillow with animated WebP support is required; unavailable geometry, saved states or
encoding support are reported without failing the report. The WebP is embedded
in the offline HTML alongside the existing derived curves; no GIF or graph
panels are included in the machine image. Use `--html PATH` for another HTML
destination. `--external-webp` saves standalone machine animations in
`<HTML-stem>.assets/` beside the HTML and uses relative links instead of embedded
image data. Keep that directory with the HTML when moving the report; each WebP
can also be opened independently.

Requested thermodynamic curves are reconstructed from saved states; Research
announces the thermodynamic replay when needed. The HTML is a derived artifact.
Comparison can inspect several campaigns or saved evaluations without inventing
a joint ranking for incompatible studies. Plot choices, evidence interpretation
and animation details belong in the [cockpit guide](DADA_ENGINE_RESEARCH_COCKPIT.md).

## Browse and save candidate graphs

```sh
research graphs path/to/campaign --candidate best
research graphs path/to/campaign --candidate best --candidate second \
    --plots positions,pressures,temperatures --save-directory work/figures
research graphs path/to/campaign-a path/to/campaign-b \
    --candidate CANDIDATE_A_ID --candidate CANDIDATE_B_ID --plots volumes,heat
```

One Matplotlib window uses the native toolbar: Back/Forward traverse graph and
zoom history, Home restores the original view, and Save opens the backend's usual
export dialog. No extra navigation or save buttons are added. PageUp/PageDown
select the previous/next graph directly, including after a zoom has replaced
forward history. Selected candidates are overlaid with candidate IDs and series
labels. Without selectors, the existing report policy prepares the two best
candidates per source. The default graph set includes all report curves except
mechanism animations; request `--plots mechanisms` to add live machine pages. Use `--plots` to request only the curves you need.

`--save-directory` sets the initial path of the native save dialog; it does not
save files automatically. The native dialog handles format and overwrite choices.
Programmatic `save_current(path)` provides a graph-only PNG/SVG/PDF export with
the current zoom and refuses overwrite. `--no-show` supports headless preparation,
but does not itself save a graph. Install the `plot` extra for Matplotlib; GUI
availability depends on its configured backend.

This viewer reuses the report sampling and replay cache. Thermal curves may
require one announced saved-state cycle replay when uncached; they do not run a
new search or periodic reconvergence. Navigation and saving only redraw prepared
curves and never replay or integrate. Unavailable curves show their reason.
Comparison across different studies has no combined scientific ranking; review
the source definitions before interpreting the overlaid curves.

### Live machine drawing and animated SVG export

```sh
research graphs path/to/campaign --candidate best --plots mechanisms
research graphs path/to/campaign --candidate best --candidate second \
    --plots positions,torques,mechanisms --save-directory work/figures
```

The machine is drawn directly with the report's production Matplotlib renderer,
not by decoding its WebP. Each selected candidate gets a separate machine page. Machine pages have no
visible title or annotations and fit the full-cycle envelope with small constant
margins; curve titles return when navigating back to a graph.
Geometry and cycle samples are prepared once; playback, navigation and saving
perform no thermodynamic replay. Space pauses/resumes the live animation.

Use the native Save action and choose SVG to export the **complete animation**:
66 exact renderer frames, 85 ms each, looping over 5.61 seconds. The standalone
SVG uses discrete SMIL visibility; there is no geometry or physics interpolation.
The validated linkages, cylinder outlines and glyphs remain vector paths. The
thermal conduits use native SVG linear gradients aligned to exact exchanger ports.
All assets are incorporated in the file. Browsers with SMIL support animate it;
viewers without support retain a static first-frame fallback. A wiki may sanitize
SVG animation, so check its SVG embedding policy. PNG/PDF save the current frame.

Vector conduits are standard in the live view and report renderer. Small axial
overlaps at the exchanger ports and cylinder walls hide antialiasing seams without
moving the thermal components. Report retains animated WebP output; its frames
are rasterized from vector-strip thermal fills before WebP encoding. SVG export
uses native `linearGradient` fills, preserving clips, valve directions and draw
order. Gas inside the cylinders is a solid vector rectangle.

```sh
research report path/to/campaign --plots mechanisms --show-mechanisms none
research graphs path/to/campaign --plots mechanisms --show-mechanisms small
research graphs path/to/campaign --plots mechanisms --show-mechanisms large
```

`--show-mechanisms` defaults to `both`. The `none` mode hides both linkages and retains the thermal assembly.
The `small` and `large` modes show only the selected mechanism, piston, gas and
cylinder, without exchangers or the other cylinder. Framing follows the visible
assembly; cylinder dimensions, piston motion and component positions are unchanged. `--external-webp` still saves report animations as separate assets.

`save_current("machine.svg")` also exports the complete animation and refuses to
overwrite. SVG animation uses the renderer's fixed full-machine framing, regardless
of a temporary interactive zoom; ordinary curve SVGs remain static exports.

### Gas torque at the cranks

```sh
research graphs path/to/campaign --candidate best --plots torques
research report path/to/campaign --candidate best --plots torques
```

`torques` uses exact replay pressures and production piston volume rates:
`tau_progress = P * dV/dt / abs(omega)`. For a physical mechanism, its signed
crank torque is `tau_progress / (study_direction * crank_direction)`. Thus its
product with signed crank speed equals that cylinder's indicated gas power.
SMALL and LARGE retain their own crank conventions. The total curve is conjugate
to positive cycle progress, not a sum of differently oriented shaft torques.
Abstract motion laws expose generalized progress loads, not physical crank loads.

These are **indicated gas-pressure torques** in N m, matching the existing absolute
pressure work boundary. Piston back-face pressure, inertia, friction, bearings,
auxiliaries and shaft efficiency are not added; this is not total required drive
or useful output torque. No new thermodynamic evaluation is performed beyond the
standard saved-state report replay, and stored candidate metrics are unchanged.

## Refine a promising result

```sh
research refine outputs/my_study/campaign --candidate best --radius 0.20 \
    --output outputs/local/study.toml
research validate outputs/local/study.toml
```

`refine` creates a portable schema-3 study around selected stored candidates;
it does not launch a campaign or copy their results as completed evaluations.
Sources must have the same scientific study identity and an active parameter
space. Repeat `--candidate` to retain several centers. The radius is a fraction
of the original normalized search width, not of each physical parameter value.

Review the generated bounds and then use `run` when ready. Exact scheduling and
radius semantics are in
[Local refinement and explicit initial evaluations](DADA_ENGINE_RESEARCH_REFERENCE.md#local-refinement-and-explicit-initial-evaluations).

## Change capacity

```sh
research rescale outputs/my_study/campaign --candidate best --factor 2 \
    --output outputs/larger/study.toml
research validate outputs/larger/study.toml
```

`rescale` creates a new portable schema-3 study without integrating or modifying
the source. It preserves declared constraints and never overwrites destination
files. Global and local candidates are accepted; source local basins are not reused
in the new study. Review the result before evaluation or search.

The [capacity reference](DADA_ENGINE_RESEARCH_CAPACITY.md) defines which
quantities scale and which remain fixed.

## Further reading

- [Research technical reference](DADA_ENGINE_RESEARCH_REFERENCE.md): schema,
  scientific identity, policies, search, persistence and runtime compatibility.
- [Cockpit and reports](DADA_ENGINE_RESEARCH_COCKPIT.md): progress, recovery,
  result inspection and evidence presentation.
- [Kinematics and mechanisms](DADA_ENGINE_RESEARCH_KINEMATICS.md): motion families,
  physical mechanisms and mechanical diagnostics.
- [Capacity scaling](DADA_ENGINE_RESEARCH_CAPACITY.md): transformation rules.
- [Current validation and limitations](validation.md): what the available
  software and physical evidence establishes.

Machine animations use vector thermal conduits aligned to exchanger ports.
`report --plots mechanisms` still exports animated WebP; `graphs` saves animated
SVG through the native Matplotlib save action. Both commands accept
`--show-mechanisms both|none|small|large`: `none` shows only pistons, cylinders,
gas and the thermal circuit. `small` and `large` omit the thermal circuit and the other cylinder.
Visibility changes framing, not machine geometry.
