# Compact finite-rod synthesis pilot

The local `dada-4bar-synthesis` implementation uses differential evolution over
four-bar ratios and phase, with a linear least-squares projection fit inside each
evaluation. The latter shortcut does not carry over to finite slider rods. The
outer search method does: this pilot evaluates the established
`FourBarSliderAssembly` directly, without thermodynamic integration or another
implementation of the mechanism equations.

Run:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 examples/search_compact_motor_motion.py --iterations 100
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 examples/plot_compact_motor_motion.py
```

The seed is `cooling_cell_x2_reference.toml`. Both sides are searched separately,
retaining their rocker outputs, closure branches, finite rods and crank radius.
Seven coordinates per side vary link lengths, rocker output coordinates, rod
length, slider-axis offset and phase. Length bounds remain near the reference.
The target uses uncentered 0–1 volume and its corresponding analytic derivative;
there is no factor-of-two conversion. The objective is
`sqrt(position_RMS² + 0.1 * derivative_RMS²)`.

Each sampled candidate must retain at least 95% of the reference's:

* stroke / swept mechanism bounding-box diagonal;
* minimum absolute sine of the coupler/rocker transmission angle;
* minimum cosine of the rod/slider-axis angle.

The envelope includes the crank origin, fixed rocker pivot, crank pin, coupler
joint, output point and slider throughout the cycle. It excludes cylinder bore,
clearance, bearings, link thickness and supports. It is an explicit compactness
proxy, not a full machine envelope. Normalization cannot erase the stroke in
these constraints. A 2880-angle check is reported separately from the coarse
search; it is still a sampled screen rather than a proof between samples.

The JSON saves exact assemblies, phase, direction, bounds policy, seed, best
improvements and reference/fit metrics. Side phase offsets are local frames;
assembling one common crank requires rotating each side's complete geometry by
the negative local offset. Do not interpret different offsets as independent
physical crank pins. This pilot does not yet export a combined machine.

The first 100-generation run yields approximately 6.03% / 5.58% position RMS
for S/L, with stroke/envelope approximately 14.3% / 17.8% (reference 13.0%).
Both dense geometric screens pass. These are candidates, not optima or validated
machines. Both fits approach the lower link-ratio bounds, which motivates further
bounded exploration rather than a claim of convergence to the best mechanism.

Next: preserve several trade-offs and restarts, rotate/export the paired geometry
in the shared-crank frame, and inspect real layout/clearances. Collision screening
requires explicit component thickness and axial layering; planar line crossings
alone are not evidence of collision. Only selected candidates should then receive
a K2 thermodynamic evaluation. No thermal efficiency was calculated in this pilot.

## Coupler-point output extension

Select `--family coupler` to attach the finite piston rod to a point fixed in the
coupler frame. The existing `CouplerOutputPoint` and slider closure are reused;
there is no projection approximation. The default remains `rocker` for backward
compatibility. Coupler coordinates are expressed in reference coupler-length
units and independently bounded in [-1.5, 1.5]; they replace the two rocker-output
multipliers. Other bounds, reference-relative feasibility floors and the weighted
position/derivative objective remain identical. Family coordinate bounds differ
by definition, so this is a first family comparison, not an exhaustive ranking.

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 examples/search_compact_motor_motion.py --family coupler --iterations 100 --output outputs/compact_motor_coupler_search.json
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 examples/plot_compact_motor_motion.py
```

The seed-27 run took 31.2 seconds. Position RMS was 1.775% / 2.432% for S/L,
and derivative RMS was 0.0542 / 0.0720 per radian. Dense-screen stroke/envelope
ratios were 0.1340 / 0.1254; both dense screens passed. Transmission floors are
nearly active, and the small-side normal output coordinate approaches its bound.
These observations must remain visible before expanding the search. No axial
thickness constraint is applied and planar rod crossings are allowed. Cylinder
collisions remain an outstanding layout check, not a validated property.

## K2 evaluation of the finite-rod coupler candidate

The local side geometries are rotated into a common crank frame by
`examples/compact_coupler_geometry.py`. A regression checks slider coordinates
and analytic derivatives against the original independent frames. Motor reversal
is still applied once by the thermodynamic factory.

Evaluate with:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 examples/evaluate_motor_champion_four_bar_k2.py --finite-coupler --geometry "$PWD/outputs/compact_motor_coupler_search.json" --output outputs/compact_motor_coupler_k2.json
```

The result converged in 23 cycles, with final normalized periodic error 0.06474.
Indicated thermal efficiency is 0.18497743 and indicated gas power is 43.11694 W,
using 233.093 W external heat input. Maximum pressure is 567075 Pa, peak tube
Reynolds number 708.52 and peak tube Mach number 0.04199. The model validity
verdict is valid. Compared with the free-motion K2 reference (20.5616%, 51.3126 W),
efficiency decreases by 2.064 percentage points and power by approximately 16.0%.

This run preserves K2 exchanger hardware, gas inventory and gas hold-up volumes.
It does not account for additional ducts required by a separated-cylinder layout.
Cylinder ports, exchanger placement and connecting duct volumes remain to be
resolved before treating this result as performance of a packaged machine.
Useful mechanical output remains unknown.

## Fine harmonic K2 benchmark

The 245–255 degree scan at 0.5-degree spacing completed all 21 points.
All evaluations converged. The retained sampled phase is 249 degrees in the
solver convention: indicated thermal efficiency 18.337139%, indicated gas
power 48.346931 W. At 249.5 degrees efficiency is 18.336705%, so the peak
is effectively flat at this resolution; this is not a continuous optimum claim.
The finite-rod coupler candidate gives 18.497743% and 43.116944 W under the
same K2 hardware assumptions: only 0.160604 percentage point more efficiency,
with less indicated power. The harmonic benchmark uses ideal sinusoidal volume
laws, not finite-rod slider-crank motion. Results and the plot are stored in
outputs/motor_harmonic_k2_fine_scan.json and the corresponding PNG/CSV.

## Independent slider-crank comparison

The user selected independently synchronized crank axes and removed cylinder-head
placement from this comparison. `examples/compare_slider_crank_motion.py` fits
centered and offset finite-rod laws to the free target, without thermodynamics.
It uses x=cos(a)+sqrt(l²-(e-sin(a))²), with exact analytic velocity and extrema
located from velocity roots. Volume increasing with x produces the requested
sharper full-volume peak; volume decreasing with x gives the opposite distortion.
The piston-rod extension needed for a physical arrangement does not alter this
law. Plot labels explicitly state volume direction instead of relying on an
ambiguous inverted-mechanism label.

Each cylinder has independent phase, rod/crank ratio (2–10), and optional offset
(up to 0.4 rod length). The fit minimizes position MSE plus 0.1 derivative MSE.
The favorable centered fit has position RMS 4.04% / 3.85% (S/L); adding offset
reduces this to 3.80% / 3.80%. Corresponding offset-fit rod/crank ratios are
2.148 / 2.146 and offset/crank ratios +0.272 / -0.135. Rod angles reach about
36.3 / 31.9 degrees; transverse loads and physical guides are not evaluated.
No efficiency improvement is inferred from geometric agreement. The opposite
volume direction pushes rod length to the upper search bound, approaching a
sinusoid rather than matching the desired distortion.

Artifacts: `outputs/slider_crank_target_comparison.png` and `.json`.
Four focused tests verify derivatives, periodicity, inversion and invalid closure.

## K2 evaluation of the fitted slider-crank laws

`examples/evaluate_slider_crank_k2.py` injects the centered and offset fits into
unchanged K2 wall physics with the established inventory. Its fast analytic
travel limits were checked against the root-normalized fitting implementation.
Both runs converged with valid model verdicts: centered 22 cycles, 14.73425%
indicated efficiency and 30.80964 W; offset 24 cycles, 11.83434% and 22.46272 W.
Both are below the 249-degree harmonic benchmark (18.33714%, 48.34693 W), and
below the 40 W screening power floor. A valid model verdict does not imply a
feasible design. Mechanical losses remain unknown.

The small-minus-large geometric crank phase offsets are 114.64656 degrees
(centered) and 125.51860 degrees (offset), measured relative to each local
positive piston axis. Their opposite directed offsets are 245.35344 and
234.48140 degrees. Physical shaft indexing also depends on the relative
orientation of those local axes. With offset sliders, crank phase and piston
extremum phase are distinct. These phases came from fitting shape, not optimizing
thermal efficiency. Poor performance of these particular fitted geometries does
not establish a family-wide limitation. No plateau-specific causal conclusion
is justified by these two tests.

Full results: `outputs/slider_crank_k2.json`.

## Conventional short-rod K2 scan

`examples/scan_short_rods_k2.py` holds the fitted conventional phase difference
(113.53674 degrees), cylinder volume ranges, gas inventory and K2 exchangers
fixed while shortening both centered rods. Results (rod/crank, indicated
efficiency, indicated power): 10: 17.8723%, 50.7254 W; 6: 17.4997%,
50.9310 W; 4: 16.9452%, 50.7689 W; 3: 16.2997%, 50.1689 W;
2: 14.7492%, 47.6734 W. All four new runs converged (13, 13, 13, 22 cycles).
This is a controlled rod-length comparison, not a search for each ratio's best
phase. Reduced efficiency does not imply a proportional reduction in power,
since external heat input also changes. Full results: outputs/short_rods_k2.json.

## Independent small-cylinder E search

Run `PYTHONPATH=src python3 examples/search_small_rocker_e.py` to search the
small cylinder only. Nine coordinates in crank-radius units vary ground,
coupler and rocker lengths, rocker output radius and angle, finite rod length,
slider perpendicular offset and orientation, and phase. The ground lies on the
positive X axis to remove global rotation redundancy; slider longitudinal origin
is a coordinate gauge and is fixed. Four seeded searches cover both four-bar
closure branches and both volume directions. The positive slider branch is used;
reversing the freely oriented slider axis covers the opposite physical extension.
No shared-pin placement, axial thickness or planar rod-crossing constraint is
imposed. Cylinder-body collisions are not certified.

The 160-generation searches took 151.5 seconds in total. The best densely
validated candidate has position RMS 0.0288901 and derivative RMS 0.0672571 per
radian, with stroke/envelope 0.239884. Lengths are normalized by crank radius,
not fixed dimensional sizing. Four-bar transmission sine is 0.49258 and rod-axis
cosine 0.94937, near the retained floors 0.48733 and 0.94927. Evaluation uses
analytic backend velocities; normalization and validation sample the cycle at
3 degrees during search and 0.25 degrees afterward. This is not a continuous
clearance proof or a global optimum. The velocity shoulder of the free target
is not reproduced. E improves substantially over the previous restricted E fit,
but these searches do not establish whether another E geometry can fit closer.

Outputs: `outputs/small_rocker_e_search.json` and `.png`. Only the small-cylinder
kinematic search was performed; F and thermodynamic integration were not run.

## Extended E search

Eight seeded runs of up to 400 generations reuse previous branch/sign winners.
Each initial population contains 75% nearby perturbations and 25% global samples;
all original bounds and feasibility floors remain unchanged. Reproduction:

```sh
PYTHONPATH=src python3 examples/search_small_rocker_e.py --iterations 400 --restarts 2 --seed 527 --warm-start outputs/small_rocker_e_search.json --output outputs/small_rocker_e_extended.json
PYTHONPATH=src python3 examples/report_extended_e.py
```

The run completed in 991.58 seconds. Best position RMS decreased from 0.0288901
to 0.0231822 and derivative RMS from 0.0672571 to 0.0519224 per radian. The
weighted objective decreased from 0.0358747 to 0.0284079. All eight retained
candidates passed an additional 0.05-degree geometry screen. The best has
stroke/envelope 0.538028, minimum transmission sine 0.701627, and minimum
rod-axis cosine 0.949439. The latter remains close to the unchanged 0.94927
floor. These dimensionless envelope metrics exclude cylinder bodies and do not
certify a packaged layout. The target velocity shoulder remains absent.

No thermodynamic integration or F search was performed. The side remains the
small cylinder only. Outputs are the extended JSON and
`outputs/small_rocker_e_extended_comparison.png`. Twelve targeted mechanism and
search tests passed before the extended run; the reporting step additionally
checked all retained candidate geometries on the finer angular grid.

## Independent small-cylinder F search

The shared search supports `--family F`, replacing the rocker-fixed output E
with a coupler-fixed output F using the existing finite-rod backend. E remains
the default. Warm-start files from the other family are rejected. Reproduction:

```sh
PYTHONPATH=src python3 examples/search_small_rocker_e.py --family F --iterations 200 --seed 827 --output outputs/small_coupler_f_search.json
PYTHONPATH=src python3 examples/report_small_f.py
```

Four branch/sign searches completed in 250.58 seconds. The best F has position
RMS 0.0110050, derivative RMS 0.0325406 per radian, and weighted score 0.0150665.
All four winners passed the additional 0.05-degree geometry screen. The best's
stroke/envelope is 0.124164, close to the 0.12384 floor, and transmission sine
is 0.49049, close to the 0.48733 floor. Rod-axis cosine is 0.98125. Thus F offers
better shape agreement than extended E, but substantially less compact geometry.
The velocity shoulder is partially reproduced, displaced and too high relative
to the target. No thermodynamic performance or collision clearance is inferred.
The E search received more restarts/generations; neither result is a proven
family optimum. Comparison: outputs/small_e_f_comparison.png. Thirteen targeted
tests passed, including the F finite-rod derivative regression.

## F animation and extended search

`examples/animate_small_f.py` renders the original 200-generation F winner in
`outputs/small_coupler_f.gif`. The rigid coupler triangle, output point F, finite
rod and slider retain the same length scale. The guide is schematic and does
not imply a sized cylinder. Playback follows motor-direction traversal.

Four additional 400-generation F runs were started from the original family
winners, using seed 1227 and unchanged bounds/floors. Outputs are
`outputs/small_coupler_f_extended.json` and
`outputs/small_e_f_extended_comparison.png`. The best position RMS improves to
0.00883387 and derivative RMS to 0.0308273 per radian. Stroke/envelope remains
near its floor (0.123943), as does rod-axis cosine (0.949324). The original GIF
is intentionally retained as the original candidate's visualization, not silently
replaced with this extended winner. No thermal evaluation was performed.

## Opposite-side piston check

`examples/check_f_opposite_slider.py` changes only the extended F slider branch
from +1 to -1, preserving loop geometry, phase, guide axis, rod length and the
volume-increases-with-coordinate convention. The lower-side chamber can then
be placed beyond the piston, away from the mechanism, without a rod crossing the
head. Cylinder width/clearance in the drawing remain schematic.

The laws are not equivalent: target position RMS changes from 0.00883387 to
0.0895395 and derivative RMS from 0.0308273 to 0.106814 per radian. Travel changes
from 1.22861 to 1.97542 crank radii, so preserving physical volume limits requires
corresponding piston-area sizing. Slider coordinates follow
x = longitudinal_projection +/- sqrt(rod_length² - transverse_offset²);
the variable square-root term changes sign, rather than providing a constant
translation. A long rod relative to crank radius is not sufficient to neglect
this correction relative to the small resulting piston stroke.

Outputs: `outputs/f_opposite_slider.json`, `f_opposite_slider_comparison.png`,
and `f_lower_piston_layout.png`. No efficiency is inferred from this shape check.
