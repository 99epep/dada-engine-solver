"""Stage F4: two-extrema motion law with four crossed endpoint flatnesses.

This stage keeps the three timing variables from F3 and replaces the shared
low/high flatness pair by four independent S-side endpoint flatnesses.  The
corresponding L-side values are imposed by the mechanical mirror symmetry
requested for this study.

Seven searched quantities
-------------------------

    small_min_deg
    small_max_deg
    large_max_deg
    h_s_low_before
    h_s_low_after
    h_s_high_before
    h_s_high_after

"before" and "after" are defined in the direction of increasing study angle.
The L-side endpoint flatnesses are not independent:

    h_L_low_before  = h_S_low_after
    h_L_low_after   = h_S_low_before
    h_L_high_before = h_S_high_after
    h_L_high_after  = h_S_high_before

For a monotone branch, with u in [0, 1]:

    y = (1 - cos(pi*u)) / 2
    E = y - h0*y*(1-y)^2 + h1*y^2*(1-y)

h0 controls the departure endpoint and h1 the arrival endpoint.  Positive h
flattens the corresponding extremum; negative h sharpens it.  V and dV/dtheta
remain continuous at every extremum.  d2V/dtheta2 may jump, so extrema are
exposed as integration breakpoints.

The large-cylinder minimum remains fixed at 180 deg as the global-phase gauge.
Comparison basis is unchanged from F1/F2/F3: same total gas mass, cylinder
volume limits, 2 Hz operation and dynamic-wall microtube hardware, with no
mechanical-loss model.

Default search: 2**6 = 64 Sobol points plus explicit controls.

Outputs:
    outputs/motor_two_extrema_stageF4.csv
    outputs/motor_two_extrema_stageF4.json
    outputs/motor_two_extrema_stageF4_motion_overlay.png

The CSV is checkpointed after every evaluation.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np
from scipy.stats import qmc

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.kinematics import HarmonicVolumeKinematics

from compare_motor_motion_laws_stage7A5 import (
    ROOT,
    A5_CAMPAIGN,
    _candidate_design_and_mass,
    _evaluate,
    _same_inventory_design,
)


TAU = 2.0 * math.pi
LARGE_MIN_DEG = 180.0
MINIMUM_POWER_W = 40.0

F3_JSON = ROOT / "outputs" / "motor_two_extrema_stageF3.json"
OUTPUT_CSV = ROOT / "outputs" / "motor_two_extrema_stageF5.csv"
OUTPUT_JSON = ROOT / "outputs" / "motor_two_extrema_stageF5.json"
OUTPUT_PLOT = ROOT / "outputs" / "motor_two_extrema_stageF5_motion_overlay.png"

PARAMETERS = (
    "small_min_deg",
    "small_max_deg",
    "large_max_deg",
    "h_s_low_before",
    "h_s_low_after",
    "h_s_high_before",
    "h_s_high_after",
)

# F3 showed that useful timing lives close to the harmonic gauge but its best
# point pushed S_max downward, L_max upward and high-flatness downward.  F4
# keeps a compact timing box while opening the new left/right shape freedom.
#
# The low-side ranges deliberately cross zero to test the hypothesis that one
# side of the minimum wants to flatten while the symmetry-paired side sharpens.

BOUNDS = {
    "small_min_deg": (68.0, 82.0),
    "small_max_deg": (237.0, 250.0),
    "large_max_deg": (356.0, 370.0),

    "h_s_low_before": (-0.8, 0.3),
    "h_s_low_after":  (-0.8, 0.3),

    "h_s_high_before": (-2.50, -1.00),
    "h_s_high_after":  (-1.20, -0.60),
}



@dataclass(frozen=True, slots=True)
class TwoExtremaCrossedFlatnessKinematics:
    """Periodic C1 two-extrema law with crossed S/L endpoint flatness."""

    small_limits: CylinderVolumeLimits
    large_limits: CylinderVolumeLimits
    small_min_deg: float
    small_max_deg: float
    large_max_deg: float
    h_s_low_before: float
    h_s_low_after: float
    h_s_high_before: float
    h_s_high_after: float

    def __post_init__(self) -> None:
        if not 0.0 < self.small_min_deg < LARGE_MIN_DEG:
            raise ValueError("small_min_deg must lie between 0 and 180 deg.")
        if not LARGE_MIN_DEG < self.small_max_deg < self.large_max_deg:
            raise ValueError(
                "Expected 180 < small_max_deg < large_max_deg."
            )
        if not self.large_max_deg < self.small_min_deg + 360.0:
            raise ValueError(
                "large_max_deg must precede the next S minimum."
            )
        for name in (
            "h_s_low_before",
            "h_s_low_after",
            "h_s_high_before",
            "h_s_high_after",
        ):
            value = getattr(self, name)
            if not math.isfinite(value) or not -2.5 <= value <= 1.0:
                raise ValueError(f"{name} must lie in [-2.5, 1].")

    @property
    def small_volume_limits(self):
        return self.small_limits

    @property
    def large_volume_limits(self):
        return self.large_limits

    @property
    def small_physical_stroke(self):
        return None

    @property
    def large_physical_stroke(self):
        return None

    @staticmethod
    def _ease(
        u: float,
        start_flatness: float,
        end_flatness: float,
    ) -> tuple[float, float, float]:
        """Return E, dE/du and d2E/du2."""
        y = 0.5 * (1.0 - math.cos(math.pi * u))
        dy = 0.5 * math.pi * math.sin(math.pi * u)
        d2y = 0.5 * math.pi**2 * math.cos(math.pi * u)

        h0 = start_flatness
        h1 = end_flatness

        a = y * (1.0 - y) ** 2
        b = y * y * (1.0 - y)
        e = y - h0 * a + h1 * b

        de_dy = (
            1.0
            - h0 * (1.0 - 4.0 * y + 3.0 * y * y)
            + h1 * (2.0 * y - 3.0 * y * y)
        )
        d2e_dy2 = (
            h0 * (4.0 - 6.0 * y)
            + h1 * (2.0 - 6.0 * y)
        )

        de = de_dy * dy
        d2e = d2e_dy2 * dy * dy + de_dy * d2y
        return e, de, d2e

    @classmethod
    def _motion(
        cls,
        theta: float,
        limits: CylinderVolumeLimits,
        minimum_deg: float,
        maximum_deg: float,
        low_before: float,
        low_after: float,
        high_before: float,
        high_after: float,
    ) -> tuple[float, float, float]:
        minimum = math.radians(minimum_deg) % TAU
        maximum = math.radians(maximum_deg) % TAU
        angle = theta % TAU

        rise_width = (maximum - minimum) % TAU
        fall_width = TAU - rise_width
        if rise_width <= 0.0 or fall_width <= 0.0:
            raise ValueError("Extrema must be distinct.")

        x = (angle - minimum) % TAU

        if x <= rise_width:
            # Minimum -> maximum: leave the minimum on its "after" side,
            # arrive at the maximum on its "before" side.
            u = x / rise_width
            e, de, d2e = cls._ease(u, low_after, high_before)
            value = limits.minimum + limits.swept * e
            first = limits.swept * de / rise_width
            second = limits.swept * d2e / (rise_width * rise_width)
        else:
            # Maximum -> next minimum: leave the maximum on its "after"
            # side, arrive at the minimum on its "before" side.
            u = (x - rise_width) / fall_width
            e, de, d2e = cls._ease(u, high_after, low_before)
            value = limits.maximum - limits.swept * e
            first = -limits.swept * de / fall_width
            second = -limits.swept * d2e / (fall_width * fall_width)

        return value, first, second

    def _small(self, theta: float) -> tuple[float, float, float]:
        return self._motion(
            theta,
            self.small_limits,
            self.small_min_deg,
            self.small_max_deg,
            self.h_s_low_before,
            self.h_s_low_after,
            self.h_s_high_before,
            self.h_s_high_after,
        )

    def _large(self, theta: float) -> tuple[float, float, float]:
        # Mechanical mirror symmetry, in increasing study-angle direction:
        #   h_L_low_before  = h_S_low_after
        #   h_L_low_after   = h_S_low_before
        #   h_L_high_before = h_S_high_after
        #   h_L_high_after  = h_S_high_before
        return self._motion(
            theta,
            self.large_limits,
            LARGE_MIN_DEG,
            self.large_max_deg,
            self.h_s_low_after,
            self.h_s_low_before,
            self.h_s_high_after,
            self.h_s_high_before,
        )

    def small_cylinder_volume(self, theta: float) -> float:
        return self._small(theta)[0]

    def small_cylinder_volume_derivative(self, theta: float) -> float:
        return self._small(theta)[1]

    def small_cylinder_volume_second_derivative(self, theta: float) -> float:
        return self._small(theta)[2]

    def large_cylinder_volume(self, theta: float) -> float:
        return self._large(theta)[0]

    def large_cylinder_volume_derivative(self, theta: float) -> float:
        return self._large(theta)[1]

    def large_cylinder_volume_second_derivative(self, theta: float) -> float:
        return self._large(theta)[2]

    def cylinder_volumes_and_derivatives(
        self,
        theta: float,
    ) -> tuple[float, float, float, float]:
        small = self._small(theta)
        large = self._large(theta)
        return small[0], large[0], small[1], large[1]

    def breakpoint_angles(self) -> tuple[float, ...]:
        values = (
            math.radians(self.small_min_deg) % TAU,
            math.radians(self.small_max_deg) % TAU,
            math.pi,
            math.radians(self.large_max_deg) % TAU,
        )
        return tuple(sorted(x for x in set(values) if 0.0 < x < TAU))


def _scale(point: np.ndarray) -> dict[str, float]:
    result = {}
    for x, name in zip(point, PARAMETERS, strict=True):
        lo, hi = BOUNDS[name]
        result[name] = lo + float(x) * (hi - lo)
    return result


def _make_kinematics(
    small_limits,
    large_limits,
    values: dict[str, float],
) -> TwoExtremaCrossedFlatnessKinematics:
    return TwoExtremaCrossedFlatnessKinematics(
        small_limits=small_limits,
        large_limits=large_limits,
        **values,
    )


def _sector_durations(values: dict[str, float]) -> dict[str, float]:
    smin = values["small_min_deg"]
    smax = values["small_max_deg"]
    lmax = values["large_max_deg"]
    return {
        "sector_smin_to_lmin_deg": LARGE_MIN_DEG - smin,
        "sector_lmin_to_smax_deg": smax - LARGE_MIN_DEG,
        "sector_smax_to_lmax_deg": lmax - smax,
        "sector_lmax_to_next_smin_deg": smin + 360.0 - lmax,
    }


def _motion_diagnostics(kinematics) -> dict[str, float]:
    angles = np.linspace(0.0, TAU, 4096, endpoint=False)
    ds = np.array([
        kinematics.small_cylinder_volume_derivative(float(x))
        for x in angles
    ])
    dl = np.array([
        kinematics.large_cylinder_volume_derivative(float(x))
        for x in angles
    ])
    d2s = np.array([
        kinematics.small_cylinder_volume_second_derivative(float(x))
        for x in angles
    ])
    d2l = np.array([
        kinematics.large_cylinder_volume_second_derivative(float(x))
        for x in angles
    ])
    return {
        "small_max_abs_dV_dtheta": float(np.max(np.abs(ds))),
        "large_max_abs_dV_dtheta": float(np.max(np.abs(dl))),
        "small_max_abs_d2V_dtheta2": float(np.max(np.abs(d2s))),
        "large_max_abs_d2V_dtheta2": float(np.max(np.abs(d2l))),
    }


def _row(index, label, values, result, kinematics):
    eta = result.get("indicated_thermal_efficiency")
    power = result.get("indicated_power_w")
    validity = result.get("validity", {})
    feasible = (
        result.get("status") == "converged"
        and eta is not None
        and math.isfinite(eta)
        and power is not None
        and power >= MINIMUM_POWER_W
        and validity.get("verdict") == "valid"
    )
    return {
        "index": index,
        "label": label,
        "feasible": feasible,
        **values,
        **_sector_durations(values),
        "indicated_thermal_efficiency": eta,
        "indicated_power_w": power,
        "heat_input_w": result.get("heat_input_w"),
        "heat_out_w": result.get("heat_out_w"),
        "maximum_tube_reynolds": result.get("maximum_tube_reynolds"),
        "maximum_tube_mach_number": result.get("maximum_tube_mach_number"),
        "maximum_absolute_mass_flow_kg_s": result.get(
            "maximum_absolute_mass_flow_kg_s"
        ),
        "cold_isothermality_error": validity.get("cold_isothermality_error"),
        "hot_isothermality_error": validity.get("hot_isothermality_error"),
        "pressure_equalization_error": validity.get(
            "maximum_pressure_equalization_error"
        ),
        **_motion_diagnostics(kinematics),
        "convergence_cycles": result.get("convergence", {}).get(
            "cycles_completed"
        ),
    }


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _load_checkpoint(path: Path) -> dict[int, dict]:
    if not path.exists():
        return {}
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        return {}
    required = set(PARAMETERS) | {"index", "label"}
    if not required.issubset(rows[0]):
        return {}
    return {int(row["index"]): row for row in rows}


def _verify_exact_harmonic(small_limits, large_limits) -> None:
    values = {
        "small_min_deg": 70.0,
        "small_max_deg": 250.0,
        "large_max_deg": 360.0,
        "h_s_low_before": 0.0,
        "h_s_low_after": 0.0,
        "h_s_high_before": 0.0,
        "h_s_high_after": 0.0,
    }
    reduced = _make_kinematics(small_limits, large_limits, values)
    harmonic = HarmonicVolumeKinematics(
        small_limits,
        large_limits,
        math.radians(250.0),
    )

    angles = np.linspace(0.0, TAU, 721, endpoint=False)
    errors = []
    derivative_errors = []
    for theta in angles:
        errors.extend((
            abs(
                reduced.small_cylinder_volume(float(theta))
                - harmonic.small_cylinder_volume(float(theta))
            ),
            abs(
                reduced.large_cylinder_volume(float(theta))
                - harmonic.large_cylinder_volume(float(theta))
            ),
        ))
        derivative_errors.extend((
            abs(
                reduced.small_cylinder_volume_derivative(float(theta))
                - harmonic.small_cylinder_volume_derivative(float(theta))
            ),
            abs(
                reduced.large_cylinder_volume_derivative(float(theta))
                - harmonic.large_cylinder_volume_derivative(float(theta))
            ),
        ))

    if max(errors) > 1e-12 or max(derivative_errors) > 1e-12:
        raise RuntimeError(
            "F4 h=0 baseline does not reproduce the exact harmonic law."
        )

def _verify_monotone_flatness_domain() -> None:
    y = np.linspace(0.0, 1.0, 10001)

    # Actual start/end flatness pairs used by the four branches:
    #
    # S min -> max : low_after  -> high_before
    # S max -> min : high_after -> low_before
    # L min -> max : low_before -> high_after
    # L max -> min : high_before -> low_after
    branch_pairs = (
        ("h_s_low_after",  "h_s_high_before"),
        ("h_s_high_after", "h_s_low_before"),
        ("h_s_low_before", "h_s_high_after"),
        ("h_s_high_before","h_s_low_after"),
    )

    for start_name, end_name in branch_pairs:
        for h0 in BOUNDS[start_name]:
            for h1 in BOUNDS[end_name]:
                derivative = (
                    1.0
                    - h0 * (1.0 - 4.0 * y + 3.0 * y * y)
                    + h1 * (2.0 * y - 3.0 * y * y)
                )
                if float(np.min(derivative)) < -1e-12:
                    raise RuntimeError(
                        f"F5 bounds can create a non-monotone branch: "
                        f"{start_name}={h0}, {end_name}={h1}"
                    )


def _load_f3_champion() -> tuple[dict, dict[str, float]]:
    data = json.loads(F3_JSON.read_text())
    best = data.get("best_feasible")
    if best is None:
        raise RuntimeError("F3 contains no feasible champion.")

    h_low = float(best["low_flatness"])
    h_high = float(best["high_flatness"])
    values = {
        "small_min_deg": float(best["small_min_deg"]),
        "small_max_deg": float(best["small_max_deg"]),
        "large_max_deg": float(best["large_max_deg"]),
        "h_s_low_before": h_low,
        "h_s_low_after": h_low,
        "h_s_high_before": h_high,
        "h_s_high_after": h_high,
    }
    return best, values


def _clamp_h(value: float) -> float:
    return min(1.0, max(-1.5, value))


def _control_candidates(f3_values: dict[str, float]):
    harmonic = {
        "small_min_deg": 70.0,
        "small_max_deg": 250.0,
        "large_max_deg": 360.0,
        "h_s_low_before": 0.0,
        "h_s_low_after": 0.0,
        "h_s_high_before": 0.0,
        "h_s_high_after": 0.0,
    }

    controls = [("exact_harmonic", harmonic)]
    controls.append(("f3_champion_symmetric", dict(f3_values)))

    # Preserve the F3 mean low flatness while explicitly testing both signs
    # of the newly released left/right asymmetry.  With the current F3
    # champion this already puts one low-side h above zero.
    for sign, label in ((1.0, "f3_low_split_A"), (-1.0, "f3_low_split_B")):
        values = dict(f3_values)
        delta = 0.50 * sign
        values["h_s_low_before"] = _clamp_h(
            f3_values["h_s_low_before"] + delta
        )
        values["h_s_low_after"] = _clamp_h(
            f3_values["h_s_low_after"] - delta
        )
        controls.append((label, values))

    # Do the analogous directional test at the high extrema without changing
    # the F3 mean high flatness appreciably.
    for sign, label in ((1.0, "f3_high_split_A"), (-1.0, "f3_high_split_B")):
        values = dict(f3_values)
        delta = 0.50 * sign
        values["h_s_high_before"] = _clamp_h(
            f3_values["h_s_high_before"] + delta
        )
        values["h_s_high_after"] = _clamp_h(
            f3_values["h_s_high_after"] - delta
        )
        controls.append((label, values))

    return controls


def _plot_results(
    small_limits,
    large_limits,
    f3_values,
    rows,
) -> None:
    feasible = [
        row
        for row in rows
        if bool(row["feasible"])
        and row.get("indicated_thermal_efficiency") not in (None, "")
        and math.isfinite(float(row["indicated_thermal_efficiency"]))
    ]
    if not feasible:
        return

    feasible.sort(
        key=lambda row: float(row["indicated_thermal_efficiency"]),
        reverse=True,
    )

    selected = [
        (1, feasible[0]),
        (2, feasible[1]),
        (3, feasible[2]),
        (9, feasible[8]),
    ]

    cache = Path(tempfile.gettempdir()) / "dada_solver_matplotlib"
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib.pyplot as plt

    angle_deg = np.linspace(0.0, 360.0, 1441)
    angle_rad = np.radians(angle_deg)

    harmonic = HarmonicVolumeKinematics(
        small_limits,
        large_limits,
        math.radians(250.0),
    )
    f3_kin = _make_kinematics(small_limits, large_limits, f3_values)

    harmonic_small = np.array([
        harmonic.small_cylinder_volume(float(x)) for x in angle_rad
    ])
    harmonic_large = np.array([
        harmonic.large_cylinder_volume(float(x)) for x in angle_rad
    ])
    f3_small = np.array([
        f3_kin.small_cylinder_volume(float(x)) for x in angle_rad
    ])
    f3_large = np.array([
        f3_kin.large_cylinder_volume(float(x)) for x in angle_rad
    ])

    fig, axes = plt.subplots(2, 1, figsize=(12, 11), sharex=True)

    axes[0].plot(
        angle_deg,
        harmonic_small * 1000.0,
        linestyle="--",
        linewidth=2.0,
        label="Exact harmonic 250°",
    )
    axes[0].plot(
        angle_deg,
        f3_small * 1000.0,
        linestyle=":",
        linewidth=2.0,
        label="F3 champion reproduced in F4",
    )
    axes[1].plot(
        angle_deg,
        harmonic_large * 1000.0,
        linestyle="--",
        linewidth=2.0,
        label="Exact harmonic 250°",
    )
    axes[1].plot(
        angle_deg,
        f3_large * 1000.0,
        linestyle=":",
        linewidth=2.0,
        label="F3 champion reproduced in F4",
    )

    for rank, row in selected:
        values = {name: float(row[name]) for name in PARAMETERS}
        kin = _make_kinematics(small_limits, large_limits, values)
        small = np.array([
            kin.small_cylinder_volume(float(x)) for x in angle_rad
        ])
        large = np.array([
            kin.large_cylinder_volume(float(x)) for x in angle_rad
        ])
        label = (
            f"F4 #{rank}, index {int(row['index'])}, "
            f"eta={100*float(row['indicated_thermal_efficiency']):.4f}%"
        )
        axes[0].plot(angle_deg, small * 1000.0, label=label)
        axes[1].plot(angle_deg, large * 1000.0, label=label)

    axes[0].set_ylabel("Small-cylinder volume (L)")
    axes[1].set_ylabel("Large-cylinder volume (L)")
    axes[1].set_xlabel("Study angle (deg)")

    for axis in axes:
        axis.set_xlim(0.0, 360.0)
        axis.set_xticks(np.arange(0.0, 361.0, 45.0))
        axis.grid(True, alpha=0.3)
        axis.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.15),
            ncol=2,
            frameon=False,
        )

    fig.suptitle(
        "Stage F4 — crossed before/after endpoint flatness",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.97), h_pad=3.0)
    OUTPUT_PLOT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PLOT, dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sobol-power",
        type=int,
        default=7,
        help="2**m Sobol points; default m=6 gives 64.",
    )
    parser.add_argument("--seed", type=int, default=2718)
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse matching rows from a compatible F4 CSV.",
    )
    parser.add_argument("--csv", type=Path, default=OUTPUT_CSV)
    parser.add_argument("--output", type=Path, default=OUTPUT_JSON)
    args = parser.parse_args()

    if args.sobol_power < 1:
        raise ValueError("--sobol-power must be positive.")

    definition = CampaignDefinition(A5_CAMPAIGN)
    base_design, total_mass, _saved, _a5_best, candidate_id = (
        _candidate_design_and_mass(definition)
    )
    small_limits = base_design.configuration.machine_volumes.small_cylinder
    large_limits = base_design.configuration.machine_volumes.large_cylinder

    _verify_exact_harmonic(small_limits, large_limits)
    _verify_monotone_flatness_domain()

    f3_best, f3_values = _load_f3_champion()
    candidates = _control_candidates(f3_values)

    sampler = qmc.Sobol(
        d=len(PARAMETERS),
        scramble=True,
        seed=args.seed,
    )
    for point in sampler.random_base2(args.sobol_power):
        candidates.append(("sobol", _scale(point)))

    print("Stage F4: two extrema with crossed before/after flatness")
    print(f"Fixed gas mass: {total_mass:.12g} kg")
    print(f"L_min fixed at {LARGE_MIN_DEG:.1f} deg")
    print(
        "Symmetry: L low/high before = S low/high after; "
        "L low/high after = S low/high before"
    )
    print(
        f"{len(candidates)} evaluations "
        f"({2**args.sobol_power} Sobol + "
        f"{len(candidates) - 2**args.sobol_power} controls)"
    )
    print(
        "\nF3 champion benchmark: "
        f"eta={100*float(f3_best['indicated_thermal_efficiency']):.6f}% "
        f"P={float(f3_best['indicated_power_w']):.3f} W"
    )
    print()

    checkpoint = _load_checkpoint(args.csv) if args.resume else {}
    rows: list[dict] = []

    for index, (label, values) in enumerate(candidates):
        cached = checkpoint.get(index)
        if cached is not None:
            matches = (
                cached.get("label") == label
                and all(
                    abs(float(cached[name]) - values[name]) < 1e-12
                    for name in PARAMETERS
                )
            )
            if matches:
                row = dict(cached)
                row["index"] = int(row["index"])
                row["feasible"] = row["feasible"].strip().lower() == "true"
                for key in row:
                    if key not in ("label", "feasible"):
                        try:
                            row[key] = float(row[key])
                        except (TypeError, ValueError):
                            pass
                rows.append(row)
                eta = row.get("indicated_thermal_efficiency")
                eta_text = (
                    "n/a"
                    if eta in (None, "")
                    else f"{100*float(eta):8.5f}%"
                )
                print(f"{index:3d} CK {eta_text:>10s} {label}")
                continue

        kinematics = _make_kinematics(small_limits, large_limits, values)
        design = _same_inventory_design(base_design, total_mass, kinematics)
        result, _ = _evaluate(
            f"two_extrema_f4_{index}",
            design,
            definition,
        )
        row = _row(index, label, values, result, kinematics)
        rows.append(row)
        _write_csv(args.csv, rows)

        eta = row["indicated_thermal_efficiency"]
        eta_text = "n/a" if eta is None else f"{100*eta:8.5f}%"
        power = row["indicated_power_w"]
        power_text = "n/a" if power is None else f"{power:7.3f} W"
        flag = "OK" if row["feasible"] else "X"
        sectors = _sector_durations(values)

        print(
            f"{index:3d} {flag:2s} {eta_text:>10s} "
            f"P={power_text:>10s}  "
            f"Smin={values['small_min_deg']:6.2f} "
            f"Smax={values['small_max_deg']:6.2f} "
            f"Lmax={values['large_max_deg']:6.2f} "
            f"hSLb={values['h_s_low_before']:6.3f} "
            f"hSLa={values['h_s_low_after']:6.3f} "
            f"hSHb={values['h_s_high_before']:6.3f} "
            f"hSHa={values['h_s_high_after']:6.3f}  "
            f"sectors="
            f"{sectors['sector_smin_to_lmin_deg']:.1f}/"
            f"{sectors['sector_lmin_to_smax_deg']:.1f}/"
            f"{sectors['sector_smax_to_lmax_deg']:.1f}/"
            f"{sectors['sector_lmax_to_next_smin_deg']:.1f}"
        )

    rows.sort(key=lambda row: int(row["index"]))
    _write_csv(args.csv, rows)

    feasible = [
        row
        for row in rows
        if bool(row["feasible"])
        and row.get("indicated_thermal_efficiency") not in (None, "")
        and math.isfinite(float(row["indicated_thermal_efficiency"]))
    ]
    feasible.sort(
        key=lambda row: float(row["indicated_thermal_efficiency"]),
        reverse=True,
    )

    output = {
        "experiment": (
            "Stage F4 two-extrema study with crossed before/after flatness"
        ),
        "comparison_basis": {
            "same_total_working_gas_mass": True,
            "total_mass_kg": total_mass,
            "same_cylinder_volume_limits": True,
            "same_operation": True,
            "same_microtube_hardware": True,
            "mechanical_losses_modeled": False,
            "stage7A5_candidate_id": candidate_id,
        },
        "kinematics": {
            "large_min_deg_fixed": LARGE_MIN_DEG,
            "flatness_symmetry": {
                "h_L_low_before": "h_S_low_after",
                "h_L_low_after": "h_S_low_before",
                "h_L_high_before": "h_S_high_after",
                "h_L_high_after": "h_S_high_before",
                "before_after_direction": "increasing study angle",
                "positive_h": "flatter corresponding endpoint",
                "negative_h": "sharper corresponding endpoint",
            },
            "search_bounds": BOUNDS,
        },
        "f3_reference": {
            "best_feasible": f3_best,
            "reproduced_f4_parameters": f3_values,
        },
        "sobol": {
            "power": args.sobol_power,
            "sample_count": 2 ** args.sobol_power,
            "seed": args.seed,
            "scramble": True,
        },
        "best_feasible": feasible[0] if feasible else None,
        "ranking": feasible,
        "all_rows": rows,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n")

    _plot_results(small_limits, large_limits, f3_values, rows)

    print("\nTOP FEASIBLE")
    for rank, row in enumerate(feasible[:12], 1):
        print(
            f"{rank:2d} "
            f"eta={100*float(row['indicated_thermal_efficiency']):.6f}% "
            f"P={float(row['indicated_power_w']):.3f} W  "
            f"Smin={float(row['small_min_deg']):.3f} "
            f"Smax={float(row['small_max_deg']):.3f} "
            f"Lmax={float(row['large_max_deg']):.3f} "
            f"hSLb={float(row['h_s_low_before']):.4f} "
            f"hSLa={float(row['h_s_low_after']):.4f} "
            f"hSHb={float(row['h_s_high_before']):.4f} "
            f"hSHa={float(row['h_s_high_after']):.4f}"
        )

    print(f"\nWrote {args.csv.relative_to(ROOT)}")
    print(f"Wrote {args.output.relative_to(ROOT)}")
    if feasible:
        print(f"Wrote {OUTPUT_PLOT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
