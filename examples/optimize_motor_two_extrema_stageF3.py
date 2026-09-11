"""Stage F2: two-extrema motion-law study with shared endpoint flatness.

Each cylinder law has only two extrema:
- one minimum-volume extremum;
- one maximum-volume extremum;
- zero slope at both extrema.

The large-cylinder minimum is fixed at 180 deg to remove the irrelevant
global phase. Five quantities are searched:

    small_min_deg
    small_max_deg
    large_max_deg
    low_flatness
    high_flatness

Unlike the first F2 draft, flatness is NOT attached to a cylinder.

    low_flatness  is shared by the S and L minimum-volume extrema.
    high_flatness is shared by the S and L maximum-volume extrema.

Thus the two half-cycles of one cylinder can have different endpoint shapes,
while corresponding low/high extrema remain symmetric between cylinders.

Shape family
------------
For a monotone branch, with u in [0,1]:

    y = (1 - cos(pi*u)) / 2

and endpoint-specific flatness h0, h1:

    E = y - h0*y*(1-y)^2 + h1*y^2*(1-y)

Properties:
- E(0)=0, E(1)=1;
- dE/du=0 at both ends;
- h0=h1=0 gives the exact half-cosine;
- positive h flattens the corresponding endpoint;
- h=1 cancels the leading quadratic curvature at that endpoint;
- negative h sharpens the corresponding endpoint.

For a rising min->max branch:
    h0 = low_flatness
    h1 = high_flatness

For a falling max->min branch:
    h0 = high_flatness
    h1 = low_flatness

The Stage F2 search box [-1.5, 1] for both h parameters remains monotone.

Comparison basis is unchanged from F1:
- same total working-gas mass;
- same cylinder volume limits;
- same 2 Hz operation;
- same dynamic-wall microtube hardware;
- no mechanical-loss model.

Outputs:
    outputs/motor_two_extrema_stageF2.csv
    outputs/motor_two_extrema_stageF2.json
    outputs/motor_two_extrema_stageF2_motion_overlay.png

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
from optimize_motor_free_local_stageF1 import (
    _make_kinematics as _make_f1_kinematics,
)


TAU = 2.0 * math.pi
LARGE_MIN_DEG = 180.0
MINIMUM_POWER_W = 40.0

F1_JSON = ROOT / "outputs" / "motor_free_local_stageF1.json"
OUTPUT_CSV = ROOT / "outputs" / "motor_two_extrema_stageF3.csv"
OUTPUT_JSON = ROOT / "outputs" / "motor_two_extrema_stageF3.json"
OUTPUT_PLOT = ROOT / "outputs" / "motor_two_extrema_stageF3_motion_overlay.png"

PARAMETERS = (
    "small_min_deg",
    "small_max_deg",
    "large_max_deg",
    "low_flatness",
    "high_flatness",
)

# These ranges guarantee:
#
#   S_min < L_min(180) < S_max < L_max < next S_min
#
# for every sampled point.
BOUNDS = {
    "small_min_deg": (66.5, 82.5),
    "small_max_deg": (240.5, 260.5),
    "large_max_deg": (342.0, 363.5),
    "low_flatness": (-0.75, 0.10),
    "high_flatness": (-1.20, -0.50),
}


@dataclass(frozen=True, slots=True)
class TwoExtremaVolumeKinematics:
    """Periodic C1 law with two extrema and shared low/high flatness."""

    small_limits: CylinderVolumeLimits
    large_limits: CylinderVolumeLimits
    small_min_deg: float
    small_max_deg: float
    large_max_deg: float
    low_flatness: float
    high_flatness: float

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
        for name, value in (
            ("low_flatness", self.low_flatness),
            ("high_flatness", self.high_flatness),
        ):
            if not math.isfinite(value) or not -1.5 <= value <= 1.0:
                raise ValueError(f"{name} must lie in [-1.5, 1].")

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
        """Return E, dE/du and d2E/du2.

        start_flatness acts only on the u=0 endpoint curvature.
        end_flatness acts only on the u=1 endpoint curvature.
        """
        y = 0.5 * (1.0 - math.cos(math.pi * u))
        dy = 0.5 * math.pi * math.sin(math.pi * u)
        d2y = 0.5 * math.pi**2 * math.cos(math.pi * u)

        # E(y) = y - h0*y*(1-y)^2 + h1*y^2*(1-y)
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
        low_flatness: float,
        high_flatness: float,
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
            # Minimum -> maximum.
            u = x / rise_width
            e, de, d2e = cls._ease(
                u,
                low_flatness,
                high_flatness,
            )
            value = limits.minimum + limits.swept * e
            first = limits.swept * de / rise_width
            second = limits.swept * d2e / (rise_width * rise_width)
        else:
            # Maximum -> next minimum.
            u = (x - rise_width) / fall_width
            e, de, d2e = cls._ease(
                u,
                high_flatness,
                low_flatness,
            )
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
            self.low_flatness,
            self.high_flatness,
        )

    def _large(self, theta: float) -> tuple[float, float, float]:
        return self._motion(
            theta,
            self.large_limits,
            LARGE_MIN_DEG,
            self.large_max_deg,
            self.low_flatness,
            self.high_flatness,
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
        # V and dV/dtheta are continuous. d2V/dtheta2 can jump because
        # neighbouring branch durations differ, so expose extrema as clean
        # integration boundaries.
        values = (
            math.radians(self.small_min_deg) % TAU,
            math.radians(self.small_max_deg) % TAU,
            math.pi,
            math.radians(self.large_max_deg) % TAU,
        )
        return tuple(
            sorted(x for x in set(values) if 0.0 < x < TAU)
        )


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
) -> TwoExtremaVolumeKinematics:
    return TwoExtremaVolumeKinematics(
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
        "maximum_tube_mach_number": result.get(
            "maximum_tube_mach_number"
        ),
        "maximum_absolute_mass_flow_kg_s": result.get(
            "maximum_absolute_mass_flow_kg_s"
        ),
        "cold_isothermality_error": validity.get(
            "cold_isothermality_error"
        ),
        "hot_isothermality_error": validity.get(
            "hot_isothermality_error"
        ),
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

    # Silently ignore the interrupted first F2 draft, whose schema used
    # small_flatness / large_flatness instead of low/high flatness.
    required = set(PARAMETERS) | {"index", "label"}
    if not required.issubset(rows[0]):
        return {}

    return {int(row["index"]): row for row in rows}


def _find_extrema_deg(kinematics) -> dict[str, float]:
    degrees = np.linspace(0.0, 360.0, 14401, endpoint=False)
    radians = np.radians(degrees)
    small = np.array([
        kinematics.small_cylinder_volume(float(x)) for x in radians
    ])
    large = np.array([
        kinematics.large_cylinder_volume(float(x)) for x in radians
    ])
    return {
        "small_min_deg": float(degrees[int(np.argmin(small))]),
        "small_max_deg": float(degrees[int(np.argmax(small))]),
        "large_min_deg": float(degrees[int(np.argmin(large))]),
        "large_max_deg": float(degrees[int(np.argmax(large))]),
    }


def _align_f1_extrema_to_large_min(
    extrema: dict[str, float],
) -> dict[str, float]:
    shift = LARGE_MIN_DEG - extrema["large_min_deg"]

    def shifted(value):
        return (value + shift) % 360.0

    smin = shifted(extrema["small_min_deg"])
    smax = shifted(extrema["small_max_deg"])
    lmax = shifted(extrema["large_max_deg"])

    if smax <= LARGE_MIN_DEG:
        smax += 360.0
    if lmax <= smax:
        lmax += 360.0

    return {
        "small_min_deg": smin,
        "small_max_deg": smax,
        "large_max_deg": lmax,
        "low_flatness": 0.0,
        "high_flatness": 0.0,
    }


def _load_f1_reference(small_limits, large_limits):
    data = json.loads(F1_JSON.read_text())
    best = data["best_feasible"]
    if best is None:
        raise RuntimeError("F1 contains no feasible champion.")

    controls = {
        "small_60_deg": float(best["small_60_deg"]),
        "small_90_deg": float(best["small_90_deg"]),
        "large_150_deg": float(best["large_150_deg"]),
        "large_210_deg": float(best["large_210_deg"]),
    }
    kinematics = _make_f1_kinematics(
        small_limits,
        large_limits,
        controls,
    )
    extrema = _find_extrema_deg(kinematics)
    aligned = _align_f1_extrema_to_large_min(extrema)
    return data, best, kinematics, extrema, aligned


def _verify_exact_harmonic(small_limits, large_limits) -> None:
    values = {
        "small_min_deg": 70.0,
        "small_max_deg": 250.0,
        "large_max_deg": 360.0,
        "low_flatness": 0.0,
        "high_flatness": 0.0,
    }
    reduced = _make_kinematics(
        small_limits,
        large_limits,
        values,
    )
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
            "Two-extrema h=0 baseline does not reproduce "
            "the exact harmonic law."
        )


def _verify_monotone_search_box() -> None:
    # dE/dy is affine in h0 and h1. Therefore checking the four corners
    # of the parameter rectangle is sufficient at every y.
    y = np.linspace(0.0, 1.0, 10001)
    for h0 in (-1.5, 1.0):
        for h1 in (-1.5, 1.0):
            derivative = (
                1.0
                - h0 * (1.0 - 4.0*y + 3.0*y*y)
                + h1 * (2.0*y - 3.0*y*y)
            )
            if float(np.min(derivative)) < -1e-12:
                raise RuntimeError(
                    "Configured flatness box can create a non-monotone branch."
                )


def _plot_results(
    small_limits,
    large_limits,
    f1_kinematics,
    f1_extrema,
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
    top = feasible[:3]

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

    # Align F1 globally so that its L minimum is at 180 deg, matching F2.
    shift_deg = LARGE_MIN_DEG - f1_extrema["large_min_deg"]
    f1_arg = np.radians((angle_deg - shift_deg) % 360.0)

    harmonic_small = np.array([
        harmonic.small_cylinder_volume(float(x)) for x in angle_rad
    ])
    harmonic_large = np.array([
        harmonic.large_cylinder_volume(float(x)) for x in angle_rad
    ])
    f1_small = np.array([
        f1_kinematics.small_cylinder_volume(float(x)) for x in f1_arg
    ])
    f1_large = np.array([
        f1_kinematics.large_cylinder_volume(float(x)) for x in f1_arg
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
        f1_small * 1000.0,
        linestyle=":",
        linewidth=2.0,
        label="F1 champion, phase-aligned",
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
        f1_large * 1000.0,
        linestyle=":",
        linewidth=2.0,
        label="F1 champion, phase-aligned",
    )

    for rank, row in enumerate(top, 1):
        values = {name: float(row[name]) for name in PARAMETERS}
        kin = _make_kinematics(
            small_limits,
            large_limits,
            values,
        )
        small = np.array([
            kin.small_cylinder_volume(float(x)) for x in angle_rad
        ])
        large = np.array([
            kin.large_cylinder_volume(float(x)) for x in angle_rad
        ])
        label = (
            f"F2 #{rank}, index {int(row['index'])}, "
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
        "Stage F2 — shared low/high endpoint flatness",
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
        help="Reuse matching rows from a compatible F2 CSV.",
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
    _verify_monotone_search_box()

    _f1_data, f1_best, f1_kinematics, f1_extrema, f1_timing = (
        _load_f1_reference(small_limits, large_limits)
    )

    harmonic_values = {
        "small_min_deg": 70.0,
        "small_max_deg": 250.0,
        "large_max_deg": 360.0,
        "low_flatness": 0.0,
        "high_flatness": 0.0,
    }

    f1_cosine = dict(f1_timing)

    f1_low_flat = dict(f1_timing)
    f1_low_flat["low_flatness"] = 1.0
    f1_low_flat["high_flatness"] = 0.0

    f1_high_flat = dict(f1_timing)
    f1_high_flat["low_flatness"] = 0.0
    f1_high_flat["high_flatness"] = 1.0

    f1_both_flat = dict(f1_timing)
    f1_both_flat["low_flatness"] = 1.0
    f1_both_flat["high_flatness"] = 1.0

    candidates: list[tuple[str, dict[str, float]]] = [
        ("exact_harmonic", harmonic_values),
        ("f1_extrema_cosine", f1_cosine),
        ("f1_extrema_low_flat", f1_low_flat),
        ("f1_extrema_high_flat", f1_high_flat),
        ("f1_extrema_both_flat", f1_both_flat),
    ]

    sampler = qmc.Sobol(
        d=len(PARAMETERS),
        scramble=True,
        seed=args.seed,
    )
    for point in sampler.random_base2(args.sobol_power):
        candidates.append(("sobol", _scale(point)))

    print("Stage F2: two extrema with shared low/high flatness")
    print(f"Fixed gas mass: {total_mass:.12g} kg")
    print(f"L_min fixed at {LARGE_MIN_DEG:.1f} deg")
    print(
        f"{len(candidates)} evaluations "
        f"({2**args.sobol_power} Sobol + 5 controls)"
    )
    print(
        "\nF1 champion benchmark: "
        f"eta={100*float(f1_best['indicated_thermal_efficiency']):.6f}% "
        f"P={float(f1_best['indicated_power_w']):.3f} W"
    )
    print("F1 extrema before phase alignment:")
    for key, value in f1_extrema.items():
        print(f"  {key:18s} {value:9.3f} deg")
    print("F1 extrema timing reduced to F2 gauge:")
    for key in ("small_min_deg", "small_max_deg", "large_max_deg"):
        print(f"  {key:18s} {f1_timing[key]:9.4f}")
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
                row["feasible"] = (
                    row["feasible"].strip().lower() == "true"
                )
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
                print(
                    f"{index:3d} CK {eta_text:>10s} {label}"
                )
                continue

        kinematics = _make_kinematics(
            small_limits,
            large_limits,
            values,
        )
        design = _same_inventory_design(
            base_design,
            total_mass,
            kinematics,
        )
        result, _ = _evaluate(
            f"two_extrema_{index}",
            design,
            definition,
        )
        row = _row(
            index,
            label,
            values,
            result,
            kinematics,
        )
        rows.append(row)

        # Overwrites any incompatible interrupted first-draft F2 CSV.
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
            f"h_low={values['low_flatness']:6.3f} "
            f"h_high={values['high_flatness']:6.3f}  "
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
            "Stage F2 two-extrema study with shared low/high flatness"
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
            "shared_flatness": {
                "low_flatness": (
                    "shared by S_min and L_min"
                ),
                "high_flatness": (
                    "shared by S_max and L_max"
                ),
                "zero": "exact half-cosine endpoint curvature",
                "positive": "flatter corresponding extrema",
                "negative": "sharper corresponding extrema",
                "one": (
                    "leading quadratic curvature cancelled "
                    "at corresponding extrema"
                ),
            },
            "search_bounds": BOUNDS,
        },
        "f1_reference": {
            "best_feasible": f1_best,
            "raw_extrema_degrees": f1_extrema,
            "phase_aligned_f2_timing": f1_timing,
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

    _plot_results(
        small_limits,
        large_limits,
        f1_kinematics,
        f1_extrema,
        rows,
    )

    print("\nTOP FEASIBLE")
    for rank, row in enumerate(feasible[:12], 1):
        print(
            f"{rank:2d} "
            f"eta={100*float(row['indicated_thermal_efficiency']):.6f}% "
            f"P={float(row['indicated_power_w']):.3f} W  "
            f"Smin={float(row['small_min_deg']):.3f} "
            f"Smax={float(row['small_max_deg']):.3f} "
            f"Lmax={float(row['large_max_deg']):.3f} "
            f"h_low={float(row['low_flatness']):.4f} "
            f"h_high={float(row['high_flatness']):.4f}"
        )

    print(f"\nWrote {args.csv.relative_to(ROOT)}")
    print(f"Wrote {args.output.relative_to(ROOT)}")
    if feasible:
        print(f"Wrote {OUTPUT_PLOT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
