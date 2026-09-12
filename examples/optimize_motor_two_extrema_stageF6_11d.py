"""Stage F6-11D: free timings plus eight independent one-sided flatnesses.

This is the symmetry-release successor to Stage F6.

Searched variables (11)
-----------------------
Timing:
    small_min_deg
    small_max_deg
    large_max_deg

The large-cylinder minimum remains fixed at 180 deg as the global-phase gauge.

Independent one-sided endpoint flatnesses:
    h_s_low_before
    h_s_low_after
    h_s_high_before
    h_s_high_after
    h_l_low_before
    h_l_low_after
    h_l_high_before
    h_l_high_after

"before" and "after" are defined in the direction of increasing study angle.

Unlike F4/F5/F6, NO crossed S/L flatness symmetry is imposed.  Each cylinder
therefore owns four independent endpoint shoulder parameters.  Timings remain
free simultaneously so the optimizer can move extrema while changing their
approach/departure shapes.

Motion family
-------------
For each monotone branch, with u in [0, 1]:

    y = (1 - cos(pi*u)) / 2
    E = y - h0*y*(1-y)^2 + h1*y^2*(1-y)

h0 controls departure from the first extremum and h1 arrival at the second.
Positive h flattens an extremum; negative h sharpens it.  V and dV/dtheta are
continuous. d2V/dtheta2 may jump at extrema, so all extrema remain integration
breakpoints.

Thermodynamic basis
-------------------
- Stage F6 champion cylinder-volume limits, gas inventory and 2 Hz operation.
- Stage K2 champion exchanger asymmetry is frozen:
      k_i and k_o are loaded from motor_exchanger_asymmetry_stageK2.json.
- No mechanical losses are introduced.
- External aerodynamic/fan losses remain excluded from the objective.

Search
------
Default: 2**9 = 512 scrambled Sobol points, no controls.
Use --sobol-power 8 for a quicker 256-point pass.
CSV checkpointing and --resume are supported.

Outputs
-------
    outputs/motor_two_extrema_stageF6_11d.csv
    outputs/motor_two_extrema_stageF6_11d.json
    outputs/motor_two_extrema_stageF6_11d_motion_overlay.png
    outputs/motor_two_extrema_stageF6_11d_total_volume.png

Run
---
    PYTHONPATH=src python3 examples/optimize_motor_two_extrema_stageF6_11d.py

Resume
------
    PYTHONPATH=src python3 examples/optimize_motor_two_extrema_stageF6_11d.py --resume
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import argparse
import csv
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np
from scipy.stats import qmc

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.exchangers.microtube import MicrotubeExchanger

from compare_motor_motion_laws_stage7A5 import (
    ROOT,
    A5_CAMPAIGN,
    _candidate_design_and_mass,
    _evaluate,
    _same_inventory_design,
)
from optimize_motor_two_extrema_stageF5 import (
    TAU,
    LARGE_MIN_DEG,
    _motion_diagnostics,
    _sector_durations,
)
from optimize_motor_two_extrema_stageF6 import (
    PARAMETERS as F6_PARAMETERS,
    _make_kinematics as _make_f6_kinematics,
)
from optimize_motor_exchanger_asymmetry_stageK1 import (
    _scaled_exchanger,
)


F6_JSON = ROOT / "outputs" / "motor_two_extrema_stageF6.json"
K2_JSON = ROOT / "outputs" / "motor_exchanger_asymmetry_stageK2.json"

OUTPUT_CSV = ROOT / "outputs" / "motor_two_extrema_stageF6_11d.csv"
OUTPUT_JSON = ROOT / "outputs" / "motor_two_extrema_stageF6_11d.json"
OUTPUT_PLOT = ROOT / "outputs" / "motor_two_extrema_stageF6_11d_motion_overlay.png"
OUTPUT_TOTAL_PLOT = ROOT / "outputs" / "motor_two_extrema_stageF6_11d_total_volume.png"

MINIMUM_POWER_W = 40.0
MAXIMUM_PRESSURE_PA = 1_200_000.0
MAXIMUM_TEMPERATURE_K = 850.0
MAXIMUM_ABSOLUTE_MASS_FLOW_KG_S = 0.08

PARAMETERS = (
    "small_min_deg",
    "small_max_deg",
    "large_max_deg",
    "h_s_low_before",
    "h_s_low_after",
    "h_s_high_before",
    "h_s_high_after",
    "h_l_low_before",
    "h_l_low_after",
    "h_l_high_before",
    "h_l_high_after",
)

# Local-but-open timing box around the F6 basin.  The previous F6 champion was
# approximately Smin=70.23, Smax=234.64, Lmax=369.04 deg.  The timing bounds
# are widened relative to the old local box because releasing four additional
# shoulder coordinates can move the optimum extrema.
#
# All four low-side h values use the same bounds, and all four high-side h
# values use the same bounds.  Thus no hidden S/L or before/after preference is
# encoded in the new independent coordinates.
BOUNDS = {
    "small_min_deg": (65.0, 80.0),
    "small_max_deg": (230.0, 245.0),
    "large_max_deg": (362.0, 374.0),

    "h_s_low_before": (-0.60, 0.35),
    "h_s_low_after": (-0.60, 0.35),
    "h_s_high_before": (-1.80, -0.40),
    "h_s_high_after": (-1.80, -0.40),

    "h_l_low_before": (-0.60, 0.35),
    "h_l_low_after": (-0.60, 0.35),
    "h_l_high_before": (-1.80, -0.40),
    "h_l_high_after": (-1.80, -0.40),
}


@dataclass(frozen=True, slots=True)
class TwoExtremaIndependentFlatnessKinematics:
    """Periodic C1 two-extrema law with eight independent one-sided h values."""

    small_limits: CylinderVolumeLimits
    large_limits: CylinderVolumeLimits

    small_min_deg: float
    small_max_deg: float
    large_max_deg: float

    h_s_low_before: float
    h_s_low_after: float
    h_s_high_before: float
    h_s_high_after: float

    h_l_low_before: float
    h_l_low_after: float
    h_l_high_before: float
    h_l_high_after: float

    def __post_init__(self) -> None:
        if not 0.0 < self.small_min_deg < LARGE_MIN_DEG:
            raise ValueError("small_min_deg must lie between 0 and 180 deg.")
        if not LARGE_MIN_DEG < self.small_max_deg < self.large_max_deg:
            raise ValueError("Expected 180 < small_max_deg < large_max_deg.")
        if not self.large_max_deg < self.small_min_deg + 360.0:
            raise ValueError(
                "large_max_deg must precede the next S minimum."
            )

        for name in PARAMETERS[3:]:
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
            # Minimum -> maximum:
            # leave low on its after side, arrive high on its before side.
            u = x / rise_width
            e, de, d2e = cls._ease(u, low_after, high_before)
            value = limits.minimum + limits.swept * e
            first = limits.swept * de / rise_width
            second = limits.swept * d2e / (rise_width * rise_width)
        else:
            # Maximum -> next minimum:
            # leave high on its after side, arrive low on its before side.
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
        return self._motion(
            theta,
            self.large_limits,
            LARGE_MIN_DEG,
            self.large_max_deg,
            self.h_l_low_before,
            self.h_l_low_after,
            self.h_l_high_before,
            self.h_l_high_after,
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


def _make_kinematics(
    small_limits,
    large_limits,
    values: dict[str, float],
) -> TwoExtremaIndependentFlatnessKinematics:
    return TwoExtremaIndependentFlatnessKinematics(
        small_limits=small_limits,
        large_limits=large_limits,
        **values,
    )


def _scale(point: np.ndarray) -> dict[str, float]:
    result = {}
    for x, name in zip(point, PARAMETERS, strict=True):
        lo, hi = BOUNDS[name]
        result[name] = lo + float(x) * (hi - lo)
    return result


def _load_f6_champion() -> tuple[dict, dict[str, float]]:
    data = json.loads(F6_JSON.read_text())
    best = data.get("best_feasible")
    if best is None:
        raise RuntimeError("F6 contains no feasible champion.")
    values = {name: float(best[name]) for name in F6_PARAMETERS}
    return best, values


def _load_k2_champion() -> dict:
    data = json.loads(K2_JSON.read_text())
    best = data.get("best_feasible")
    if best is None:
        raise RuntimeError("K2 contains no feasible champion.")
    return best


def _independent_reference_values(
    f6_values: dict[str, float],
) -> dict[str, float]:
    """Express the crossed F6 champion exactly in the new 11D coordinates."""
    return {
        "small_min_deg": f6_values["small_min_deg"],
        "small_max_deg": f6_values["small_max_deg"],
        "large_max_deg": f6_values["large_max_deg"],

        "h_s_low_before": f6_values["h_s_low_before"],
        "h_s_low_after": f6_values["h_s_low_after"],
        "h_s_high_before": f6_values["h_s_high_before"],
        "h_s_high_after": f6_values["h_s_high_after"],

        # Historical crossed symmetry:
        # L low before  = S low after
        # L low after   = S low before
        # L high before = S high after
        # L high after  = S high before
        "h_l_low_before": f6_values["h_s_low_after"],
        "h_l_low_after": f6_values["h_s_low_before"],
        "h_l_high_before": f6_values["h_s_high_after"],
        "h_l_high_after": f6_values["h_s_high_before"],
    }


def _verify_f6_reference_equivalence(
    small_limits,
    large_limits,
    f6_values: dict[str, float],
) -> None:
    """The new independent family must contain the old crossed F6 law exactly."""
    crossed = _make_f6_kinematics(
        small_limits,
        large_limits,
        f6_values,
    )
    independent = _make_kinematics(
        small_limits,
        large_limits,
        _independent_reference_values(f6_values),
    )

    angles = np.linspace(0.0, TAU, 1441, endpoint=False)
    errors = []
    derivative_errors = []

    for theta in angles:
        t = float(theta)
        errors.extend(
            (
                abs(
                    crossed.small_cylinder_volume(t)
                    - independent.small_cylinder_volume(t)
                ),
                abs(
                    crossed.large_cylinder_volume(t)
                    - independent.large_cylinder_volume(t)
                ),
            )
        )
        derivative_errors.extend(
            (
                abs(
                    crossed.small_cylinder_volume_derivative(t)
                    - independent.small_cylinder_volume_derivative(t)
                ),
                abs(
                    crossed.large_cylinder_volume_derivative(t)
                    - independent.large_cylinder_volume_derivative(t)
                ),
            )
        )

    if max(errors) > 1e-12 or max(derivative_errors) > 1e-12:
        raise RuntimeError(
            "The independent 11D family does not reproduce the F6 crossed "
            "reference exactly."
        )


def _verify_monotone_flatness_domain() -> None:
    """Prove numerically that every h-box corner gives monotone branches.

    dE/dy is affine in h0 and h1, so the extrema over each rectangular endpoint
    box occur on h-box corners.  A dense y check of all branch box corners is
    therefore sufficient for this bounded study.
    """
    y = np.linspace(0.0, 1.0, 20001)

    branch_pairs = (
        ("h_s_low_after", "h_s_high_before"),
        ("h_s_high_after", "h_s_low_before"),
        ("h_l_low_after", "h_l_high_before"),
        ("h_l_high_after", "h_l_low_before"),
    )

    worst = math.inf
    worst_case = None

    for start_name, end_name in branch_pairs:
        for h0 in BOUNDS[start_name]:
            for h1 in BOUNDS[end_name]:
                derivative = (
                    1.0
                    - h0 * (1.0 - 4.0 * y + 3.0 * y * y)
                    + h1 * (2.0 * y - 3.0 * y * y)
                )
                minimum = float(np.min(derivative))
                if minimum < worst:
                    worst = minimum
                    worst_case = (start_name, h0, end_name, h1)

    if worst < -1e-12:
        raise RuntimeError(
            "11D h bounds can create a non-monotone branch: "
            f"{worst_case}, minimum dE/dy={worst}"
        )

    print(
        "Monotonicity box check passed: "
        f"minimum sampled dE/dy={worst:.6f}"
    )


def _cross_symmetry_diagnostics(values: dict[str, float]) -> dict[str, float]:
    """Measure departure from the old crossed mechanical symmetry."""
    deltas = np.asarray(
        [
            values["h_l_low_before"] - values["h_s_low_after"],
            values["h_l_low_after"] - values["h_s_low_before"],
            values["h_l_high_before"] - values["h_s_high_after"],
            values["h_l_high_after"] - values["h_s_high_before"],
        ],
        dtype=float,
    )
    return {
        "cross_symmetry_rms_h": float(np.sqrt(np.mean(deltas * deltas))),
        "cross_symmetry_max_abs_h": float(np.max(np.abs(deltas))),
        "delta_cross_low_before": float(deltas[0]),
        "delta_cross_low_after": float(deltas[1]),
        "delta_cross_high_before": float(deltas[2]),
        "delta_cross_high_after": float(deltas[3]),
    }


def _diagnostic_maximum(
    diagnostics: dict,
    group_name: str,
) -> float | None:
    group = diagnostics.get(group_name)
    if not isinstance(group, dict):
        return None
    maxima = []
    for item in group.values():
        if isinstance(item, dict) and item.get("maximum") is not None:
            maxima.append(float(item["maximum"]))
    return max(maxima) if maxima else None


def _row(
    index: int,
    values: dict[str, float],
    result: dict,
    kinematics,
) -> dict:
    eta = result.get("indicated_thermal_efficiency")
    power = result.get("indicated_power_w")
    validity = result.get("validity", {})
    diagnostics = result.get("diagnostics", {})

    maximum_pressure = _diagnostic_maximum(
        diagnostics,
        "pressure_extrema",
    )
    maximum_temperature = _diagnostic_maximum(
        diagnostics,
        "temperature_extrema",
    )
    maximum_mass_flow = result.get("maximum_absolute_mass_flow_kg_s")

    constraints = {
        "power": (
            power is not None
            and math.isfinite(float(power))
            and float(power) >= MINIMUM_POWER_W
        ),
        "pressure": (
            maximum_pressure is not None
            and maximum_pressure <= MAXIMUM_PRESSURE_PA
        ),
        "temperature": (
            maximum_temperature is not None
            and maximum_temperature <= MAXIMUM_TEMPERATURE_K
        ),
        "mass_flow": (
            maximum_mass_flow is not None
            and maximum_mass_flow <= MAXIMUM_ABSOLUTE_MASS_FLOW_KG_S
        ),
        "validity": validity.get("verdict") == "valid",
    }

    feasible = (
        result.get("status") == "converged"
        and eta is not None
        and math.isfinite(float(eta))
        and all(constraints.values())
    )

    return {
        "index": index,
        "label": "sobol",
        "feasible": feasible,
        **values,
        **_sector_durations(values),
        **_cross_symmetry_diagnostics(values),
        "indicated_thermal_efficiency": eta,
        "indicated_power_w": power,
        "heat_input_w": result.get("heat_input_w"),
        "heat_out_w": result.get("heat_out_w"),
        "maximum_pressure_pa": maximum_pressure,
        "maximum_temperature_k": maximum_temperature,
        "maximum_tube_reynolds": result.get("maximum_tube_reynolds"),
        "maximum_tube_mach_number": result.get("maximum_tube_mach_number"),
        "maximum_absolute_mass_flow_kg_s": maximum_mass_flow,
        "cold_isothermality_error": validity.get("cold_isothermality_error"),
        "hot_isothermality_error": validity.get("hot_isothermality_error"),
        "pressure_equalization_error": validity.get(
            "maximum_pressure_equalization_error"
        ),
        **_motion_diagnostics(kinematics),
        "convergence_cycles": result.get("convergence", {}).get(
            "cycles_completed"
        ),
        "constraint_power_ok": constraints["power"],
        "constraint_pressure_ok": constraints["pressure"],
        "constraint_temperature_ok": constraints["temperature"],
        "constraint_mass_flow_ok": constraints["mass_flow"],
        "constraint_validity_ok": constraints["validity"],
    }


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
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
    return {int(float(row["index"])): row for row in rows}


def _restore_cached_row(row: dict) -> dict:
    restored = dict(row)
    restored["index"] = int(float(restored["index"]))
    for key in list(restored):
        if key == "label":
            continue
        if key == "feasible" or key.startswith("constraint_"):
            restored[key] = str(restored[key]).strip().lower() == "true"
            continue
        if restored[key] in ("", None):
            restored[key] = None
            continue
        try:
            restored[key] = float(restored[key])
        except (TypeError, ValueError):
            pass
    return restored


def _kinematic_arrays(kinematics, angle_rad: np.ndarray):
    small = np.asarray(
        [
            kinematics.small_cylinder_volume(float(x))
            for x in angle_rad
        ]
    )
    large = np.asarray(
        [
            kinematics.large_cylinder_volume(float(x))
            for x in angle_rad
        ]
    )
    dsmall = np.asarray(
        [
            kinematics.small_cylinder_volume_derivative(float(x))
            for x in angle_rad
        ]
    )
    dlarge = np.asarray(
        [
            kinematics.large_cylinder_volume_derivative(float(x))
            for x in angle_rad
        ]
    )
    return small, large, dsmall, dlarge


def _setup_matplotlib():
    cache = Path(tempfile.gettempdir()) / "dada_solver_matplotlib"
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _plot_results(
    small_limits,
    large_limits,
    rows: list[dict],
    f6_reference_values: dict[str, float],
) -> None:
    feasible = [
        row
        for row in rows
        if bool(row["feasible"])
        and row.get("indicated_thermal_efficiency") is not None
        and math.isfinite(float(row["indicated_thermal_efficiency"]))
    ]
    if not feasible:
        return

    feasible.sort(
        key=lambda row: float(row["indicated_thermal_efficiency"]),
        reverse=True,
    )
    top = feasible[:4]

    reference = _make_kinematics(
        small_limits,
        large_limits,
        f6_reference_values,
    )

    angle_deg = np.linspace(0.0, 360.0, 1441)
    angle_rad = np.radians(angle_deg)

    plt = _setup_matplotlib()

    fig, axes = plt.subplots(2, 1, figsize=(12, 11), sharex=True)

    rs, rl, _, _ = _kinematic_arrays(reference, angle_rad)
    axes[0].plot(
        angle_deg,
        rs * 1000.0,
        linestyle="--",
        linewidth=1.8,
        label="K2 / crossed-F6 motion reference",
    )
    axes[1].plot(
        angle_deg,
        rl * 1000.0,
        linestyle="--",
        linewidth=1.8,
        label="K2 / crossed-F6 motion reference",
    )

    for rank, row in enumerate(top, 1):
        values = {name: float(row[name]) for name in PARAMETERS}
        kin = _make_kinematics(small_limits, large_limits, values)
        small, large, _, _ = _kinematic_arrays(kin, angle_rad)
        label = (
            f"11D #{rank}, index {int(row['index'])}, "
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
        "Stage F6-11D — independent one-sided flatness motion overlay",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.97), h_pad=3.0)
    OUTPUT_PLOT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PLOT, dpi=160)
    plt.close(fig)

    # Total-cylinder-volume diagnostic.
    best = top[0]
    best_values = {name: float(best[name]) for name in PARAMETERS}
    best_kin = _make_kinematics(
        small_limits,
        large_limits,
        best_values,
    )

    curves = [
        ("K2 / crossed-F6 reference", reference, "--", 1.8),
        (
            f"F6-11D best, index {int(best['index'])}",
            best_kin,
            "-",
            2.2,
        ),
    ]

    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

    for label, kin, linestyle, linewidth in curves:
        small, large, dsmall, dlarge = _kinematic_arrays(
            kin,
            angle_rad,
        )
        axes[0].plot(
            angle_deg,
            (small + large) * 1000.0,
            linestyle=linestyle,
            linewidth=linewidth,
            label=label,
        )
        axes[1].plot(
            angle_deg,
            (dsmall + dlarge) * 1000.0,
            linestyle=linestyle,
            linewidth=linewidth,
            label=label,
        )

    axes[0].set_ylabel(r"$V_S + V_L$ (L)")
    axes[1].set_ylabel(r"$d(V_S+V_L)/d\theta$ (L/rad)")
    axes[1].set_xlabel("Study angle (deg)")
    axes[1].axhline(0.0, linewidth=0.8)

    event_angles = (
        best_values["small_min_deg"],
        LARGE_MIN_DEG,
        best_values["small_max_deg"],
        best_values["large_max_deg"],
    )

    for axis in axes:
        for angle in event_angles:
            axis.axvline(
                angle % 360.0,
                linestyle=":",
                linewidth=0.8,
                alpha=0.5,
            )
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
        "Stage F6-11D — total cylinder-volume diagnostic\n"
        "Kinematic indicator only: not actual exchanger mass flow",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.95), h_pad=3.0)
    fig.savefig(OUTPUT_TOTAL_PLOT, dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sobol-power",
        type=int,
        default=9,
        help="2**m Sobol points; default m=9 gives 512.",
    )
    parser.add_argument("--seed", type=int, default=6118)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--csv", type=Path, default=OUTPUT_CSV)
    parser.add_argument("--output", type=Path, default=OUTPUT_JSON)
    args = parser.parse_args()

    if args.sobol_power < 1:
        raise ValueError("--sobol-power must be positive.")

    definition = CampaignDefinition(A5_CAMPAIGN)
    base_design, total_mass, _saved, _a5_best, candidate_id = (
        _candidate_design_and_mass(definition)
    )
    if not isinstance(base_design.heat_in, MicrotubeExchanger):
        raise TypeError("Expected microtube H_i design.")
    if not isinstance(base_design.heat_out, MicrotubeExchanger):
        raise TypeError("Expected microtube H_o design.")

    f6_best, f6_values = _load_f6_champion()
    k2_best = _load_k2_champion()

    small_limits = base_design.configuration.machine_volumes.small_cylinder
    large_limits = base_design.configuration.machine_volumes.large_cylinder

    _verify_monotone_flatness_domain()
    _verify_f6_reference_equivalence(
        small_limits,
        large_limits,
        f6_values,
    )

    reference_values = _independent_reference_values(f6_values)

    # Freeze Stage K2 champion exchanger geometry while releasing only motion.
    base_with_motion = _same_inventory_design(
        base_design,
        total_mass,
        _make_kinematics(
            small_limits,
            large_limits,
            reference_values,
        ),
    )
    heat_in = _scaled_exchanger(
        base_with_motion.heat_in,
        float(k2_best["k_i"]),
    )
    heat_out = _scaled_exchanger(
        base_with_motion.heat_out,
        float(k2_best["k_o"]),
    )

    sampler = qmc.Sobol(
        d=len(PARAMETERS),
        scramble=True,
        seed=args.seed,
    )
    candidates = [
        _scale(point)
        for point in sampler.random_base2(args.sobol_power)
    ]

    print("Stage F6-11D: free timings + eight independent one-sided h values")
    print(f"Large-cylinder minimum fixed at {LARGE_MIN_DEG:.1f} deg")
    print(f"Fixed total gas mass: {total_mass:.12g} kg")
    print(
        "Frozen K2 exchanger asymmetry: "
        f"k_i={float(k2_best['k_i']):.6f}, "
        f"k_o={float(k2_best['k_o']):.6f}"
    )
    print(
        "K2 / crossed-F6 reference: "
        f"eta={100*float(k2_best['indicated_thermal_efficiency']):.6f}% "
        f"P={float(k2_best['indicated_power_w']):.3f} W"
    )
    print(
        f"{len(candidates)} evaluations "
        f"({2**args.sobol_power} Sobol, no controls)"
    )
    print()

    checkpoint = _load_checkpoint(args.csv) if args.resume else {}
    rows: list[dict] = []

    for index, values in enumerate(candidates):
        cached = checkpoint.get(index)
        if cached is not None:
            matches = (
                cached.get("label") == "sobol"
                and all(
                    abs(float(cached[name]) - values[name]) < 1e-12
                    for name in PARAMETERS
                )
            )
            if matches:
                row = _restore_cached_row(cached)
                rows.append(row)
                eta = row.get("indicated_thermal_efficiency")
                eta_text = (
                    "n/a"
                    if eta is None
                    else f"{100*float(eta):8.5f}%"
                )
                print(f"{index:3d} CK {eta_text:>10s}")
                continue

        kinematics = _make_kinematics(
            small_limits,
            large_limits,
            values,
        )
        design = replace(
            _same_inventory_design(
                base_design,
                total_mass,
                kinematics,
            ),
            heat_in=heat_in,
            heat_out=heat_out,
        )

        result, _last_state = _evaluate(
            f"two_extrema_f6_11d_{index}",
            design,
            definition,
        )
        row = _row(
            index,
            values,
            result,
            kinematics,
        )
        rows.append(row)
        _write_csv(args.csv, rows)

        eta = row["indicated_thermal_efficiency"]
        eta_text = (
            "n/a" if eta is None
            else f"{100*float(eta):8.5f}%"
        )
        power = row["indicated_power_w"]
        power_text = (
            "n/a" if power is None
            else f"{float(power):7.3f} W"
        )
        flag = "OK" if row["feasible"] else "X"

        print(
            f"{index:3d} {flag:2s} {eta_text:>10s} "
            f"P={power_text:>10s}  "
            f"Smin={values['small_min_deg']:6.2f} "
            f"Smax={values['small_max_deg']:6.2f} "
            f"Lmax={values['large_max_deg']:6.2f}  "
            f"symRMS={row['cross_symmetry_rms_h']:.3f}"
        )

    rows.sort(key=lambda row: int(row["index"]))
    _write_csv(args.csv, rows)

    feasible = [
        row
        for row in rows
        if bool(row["feasible"])
        and row.get("indicated_thermal_efficiency") is not None
        and math.isfinite(float(row["indicated_thermal_efficiency"]))
    ]
    feasible.sort(
        key=lambda row: float(row["indicated_thermal_efficiency"]),
        reverse=True,
    )

    output = {
        "experiment": (
            "Stage F6-11D: three free extrema timings and eight independent "
            "one-sided flatness parameters"
        ),
        "comparison_basis": {
            "same_total_working_gas_mass": True,
            "total_mass_kg": total_mass,
            "same_cylinder_volume_limits": True,
            "same_operation": True,
            "frequency_hz": 2.0,
            "K2_exchanger_geometry_frozen": True,
            "K2_k_i": float(k2_best["k_i"]),
            "K2_k_o": float(k2_best["k_o"]),
            "mechanical_losses_modeled": False,
            "external_aerodynamic_losses_in_objective": False,
            "stage7A5_candidate_id": candidate_id,
        },
        "kinematics": {
            "large_min_deg_fixed": LARGE_MIN_DEG,
            "parameter_count": len(PARAMETERS),
            "parameters": PARAMETERS,
            "search_bounds": BOUNDS,
            "crossed_flatness_symmetry_imposed": False,
            "continuity": "C1; acceleration may jump at extrema",
            "before_after_direction": "increasing study angle",
        },
        "references": {
            "F6_crossed_champion": f6_best,
            "K2_champion_with_F6_motion": k2_best,
            "F6_champion_expressed_in_11D_coordinates":
                reference_values,
        },
        "constraints": {
            "minimum_indicated_power_w": MINIMUM_POWER_W,
            "maximum_pressure_pa": MAXIMUM_PRESSURE_PA,
            "maximum_temperature_k": MAXIMUM_TEMPERATURE_K,
            "maximum_absolute_mass_flow_kg_s":
                MAXIMUM_ABSOLUTE_MASS_FLOW_KG_S,
            "thermodynamic_validity_required": True,
        },
        "sobol": {
            "power": args.sobol_power,
            "sample_count": 2 ** args.sobol_power,
            "seed": args.seed,
            "scramble": True,
            "control_count": 0,
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
        rows,
        reference_values,
    )

    print("\nTOP FEASIBLE")
    for rank, row in enumerate(feasible[:15], 1):
        print(
            f"{rank:2d} "
            f"eta={100*float(row['indicated_thermal_efficiency']):.6f}% "
            f"P={float(row['indicated_power_w']):.3f} W  "
            f"Smin={float(row['small_min_deg']):.3f} "
            f"Smax={float(row['small_max_deg']):.3f} "
            f"Lmax={float(row['large_max_deg']):.3f}  "
            f"symRMS={float(row['cross_symmetry_rms_h']):.4f}"
        )
        print(
            "    S: "
            f"low(before/after)="
            f"{float(row['h_s_low_before']):+.3f}/"
            f"{float(row['h_s_low_after']):+.3f}  "
            f"high(before/after)="
            f"{float(row['h_s_high_before']):+.3f}/"
            f"{float(row['h_s_high_after']):+.3f}"
        )
        print(
            "    L: "
            f"low(before/after)="
            f"{float(row['h_l_low_before']):+.3f}/"
            f"{float(row['h_l_low_after']):+.3f}  "
            f"high(before/after)="
            f"{float(row['h_l_high_before']):+.3f}/"
            f"{float(row['h_l_high_after']):+.3f}"
        )

    print(f"\nWrote {args.csv.relative_to(ROOT)}")
    print(f"Wrote {args.output.relative_to(ROOT)}")
    if feasible:
        print(f"Wrote {OUTPUT_PLOT.relative_to(ROOT)}")
        print(f"Wrote {OUTPUT_TOTAL_PLOT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
