#!/usr/bin/env python3
"""V4.1 primary panel -> second-dyad synthesis for candidate 3952 LARGE.

This is a corrective wrapper around search_large_downstream_from_primary_panel_3952.

Differences from the old downstream experiment:
- piston POSITION is fitted over the complete cycle, including the HP band;
- the target velocity is linearly bridged through the HP-noise band and used only
  as a very small tie-breaker;
- wrong-direction motion and extra reversals are penalized;
- HP acceleration itself is never fitted;
- the former hard rod/stroke <= 3 limit is disabled (it was arbitrary and would
  reject the historical ~2.08%-RMS control, whose ratio is ~3.31).

Default diversified V4.1 panel:
  1  compact near-champion
  0  scalar V4.1 champion
  4  compact alternative
  12 BP-symmetric counterexample
  9  strongest known BP-symmetry family
  50 historical ~2.08%-RMS primary control

Typical run:
  PYTHONPATH=src python3 examples/search_large_second_dyad_v41_panel_3952.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np

import search_large_downstream_from_primary_panel_3952 as base


ROOT = Path.cwd()
DEFAULT_LIBRARY = ROOT / "outputs" / "large_primary_cadence_3952_v4_1_rescore.json"
DEFAULT_OUTPUT = ROOT / "outputs" / "large_second_dyad_v41_panel_3952.json"

MONOTONICITY_WEIGHT = 0.02
EXTRA_REVERSAL_WEIGHT = 0.02
MONOTONICITY_SPEED_FLOOR = 0.03
REVERSAL_SMOOTHING_DEG = 2.0


def bridge_hp_velocity(y: np.ndarray, valid_mask: np.ndarray) -> np.ndarray:
    """Linearly bridge the one contiguous HP interval (mask False)."""
    y = np.asarray(y, float)
    hp = ~np.asarray(valid_mask, bool)
    if not np.any(hp):
        return y.copy()

    n = len(y)
    starts = [i for i in range(n) if hp[i] and not hp[(i - 1) % n]]
    ends = [i for i in range(n) if hp[i] and not hp[(i + 1) % n]]
    if len(starts) != 1 or len(ends) != 1:
        raise ValueError("expected one contiguous HP-noise band")

    left = (starts[0] - 1) % n
    right = (ends[0] + 1) % n
    span = (right - left) % n

    out = y.copy()
    for k in range(1, span):
        i = (left + k) % n
        f = k / span
        out[i] = (1.0 - f) * y[left] + f * y[right]
    return out


def smooth_periodic(y: np.ndarray, half: int) -> np.ndarray:
    y = np.asarray(y, float)
    if half <= 0:
        return y.copy()
    out = np.zeros_like(y)
    for k in range(-half, half + 1):
        out += np.roll(y, k)
    return out / (2 * half + 1)


def count_zero_crossings(v: np.ndarray) -> int:
    # Infer angular sample spacing from one complete revolution.
    step_deg = 360.0 / len(v)
    half = max(1, int(round(0.5 * REVERSAL_SMOOTHING_DEG / step_deg)))
    s = smooth_periodic(v, half)
    return sum(
        1
        for i in range(len(s))
        if s[i] == 0.0 or s[i] * s[(i + 1) % len(s)] < 0.0
    )


def evaluate_v41(
    x, *, E, Ed, psine, pdata, pv, second_branch,
    tq, tdq, mask, velocity_weight
):
    gx, gy, lef, lgf, ha, hn, rod, axis_offset, axis_angle = map(float, x)

    G = np.array((gx, gy))
    delta = G - E
    dist = np.linalg.norm(delta, axis=1)
    if np.any(dist <= 1e-10):
        raise ValueError("secondary centers")

    u0 = delta / dist[:, None]
    along = (lef * lef - lgf * lgf + dist * dist) / (2 * dist)
    h2 = lef * lef - along * along
    if np.any(h2 <= 1e-10):
        raise ValueError("secondary closure")
    h = np.sqrt(h2)
    n0 = np.column_stack((-u0[:, 1], u0[:, 0]))
    F = E + along[:, None] * u0 + second_branch * h[:, None] * n0

    ef = F - E
    gf = F - G
    det = ef[:, 0] * gf[:, 1] - ef[:, 1] * gf[:, 0]
    if np.any(np.abs(det) <= 1e-10):
        raise ValueError("secondary singular")

    rhs = np.sum(ef * Ed, axis=1)
    Fd = np.column_stack((rhs * gf[:, 1] / det, -rhs * gf[:, 0] / det))
    ssine = float(np.min(np.abs(det) / (lef * lgf)))

    u = ef / lef
    ud = (Fd - Ed) / lef
    nvec = np.column_stack((-u[:, 1], u[:, 0]))
    nd = np.column_stack((-ud[:, 1], ud[:, 0]))
    H = E + ha * lef * u + hn * lef * nvec
    Hd = Ed + ha * lef * ud + hn * lef * nd

    axis = np.array((math.cos(axis_angle), math.sin(axis_angle)))
    normal = np.array((-axis[1], axis[0]))
    origin = axis_offset * normal
    rel = H - origin
    longitudinal = rel @ axis
    transverse = rel @ normal
    longitudinal_d = Hd @ axis
    transverse_d = Hd @ normal

    margin2 = rod * rod - transverse * transverse
    if np.any(margin2 <= 1e-10):
        raise ValueError("rod closure")
    margin = np.sqrt(margin2)

    slider = longitudinal + margin
    slider_d = longitudinal_d - transverse * transverse_d / margin
    stroke = float(np.ptp(slider))
    if stroke <= 1e-10:
        raise ValueError("zero stroke")

    q = 1.0 - (slider - float(np.min(slider))) / stroke
    dq = -slider_d / stroke

    valid = np.asarray(mask, bool)
    hp = ~valid

    # Critical correction: never remove HP from the position objective.
    pos = float(np.sqrt(np.mean((q - tq) ** 2)))
    pos_non_hp = float(np.sqrt(np.mean((q[valid] - tq[valid]) ** 2)))
    pos_hp = (
        float(np.sqrt(np.mean((q[hp] - tq[hp]) ** 2)))
        if np.any(hp) else 0.0
    )

    tdq_bridge = bridge_hp_velocity(tdq, valid)
    der = float(np.sqrt(np.mean((dq - tdq_bridge) ** 2)))

    # Soft monotonicity guard. Ignore target-turnaround neighborhoods by
    # requiring at least 3% of target peak speed before assigning a sign.
    target_scale = max(float(np.max(np.abs(tdq_bridge))), 1e-12)
    active = np.abs(tdq_bridge) >= MONOTONICITY_SPEED_FLOOR * target_scale
    desired = np.sign(tdq_bridge[active])
    wrong = np.maximum(-dq[active] * desired, 0.0)
    dq_scale = max(float(np.sqrt(np.mean(dq * dq))), 1e-12)
    wrong_rms = (
        float(np.sqrt(np.mean(wrong * wrong))) / dq_scale
        if np.any(active) else 0.0
    )

    zeros = count_zero_crossings(dq)
    extra = max(0, zeros - 2)

    score = math.sqrt(
        pos * pos
        + velocity_weight * der * der
        + MONOTONICITY_WEIGHT * wrong_rms * wrong_rms
    ) + EXTRA_REVERSAL_WEIGHT * extra

    rod_cos = float(np.min(margin / rod))
    clearance, inside = base.stage.signed_triangle_clearance_origin(E, F, H)
    line = base.stage.best_line_metrics(H, stroke, axis)

    tc = transverse - float(np.mean(transverse))
    lat_rms = float(np.sqrt(np.mean(tc * tc))) / stroke
    lat_span = float(np.ptp(transverse)) / stroke
    eh = lef * math.hypot(ha, hn)

    full = np.concatenate((pv, np.asarray(x, float)))
    out = {
        "score": score,
        "motion_score": score,
        "position_rms": pos,
        "position_rms_raw": pos,
        "position_rms_full_cycle": pos,
        "position_rms_non_hp": pos_non_hp,
        "position_rms_hp_band": pos_hp,
        "derivative_rms": der,
        "derivative_rms_raw": der,
        "derivative_rms_bridged": der,
        "wrong_sign_velocity_rms": wrong_rms,
        "zero_crossing_count": int(zeros),
        "extra_reversals": int(extra),
        "stroke_over_crank": stroke,
        "minimum_primary_transmission_sine": float(psine),
        "minimum_secondary_transmission_sine": ssine,
        "minimum_rod_axis_cosine": rod_cos,
        "EF_over_crank": lef,
        "GF_over_crank": lgf,
        "EH_over_EF": math.hypot(ha, hn),
        "EH_over_crank": eh,
        "piston_rod_over_crank": rod,
        "piston_rod_over_stroke": rod / stroke,
        "H_axis_lateral_rms_over_stroke": lat_rms,
        "H_axis_lateral_span_over_stroke": lat_span,
        "crank_axis_to_EFH_clearance_over_crank": clearance,
        "crank_axis_inside_EFH_frames": inside,
        "primary": dict(pdata),
        "parameters": dict(zip(base.stage.NAMES, map(float, full))),
        "second_branch": int(second_branch),
    }
    out.update(line)
    return out


def violation_v41(r, args):
    """Historical feasibility envelope, but NO hard rod/stroke ceiling."""
    v = 0.0
    v += max(0.0, args.stroke_floor - r["stroke_over_crank"]) / args.stroke_floor
    v += max(0.0, r["stroke_over_crank"] - args.stroke_ceiling) / args.stroke_ceiling
    v += max(
        0.0,
        args.secondary_sine_floor - r["minimum_secondary_transmission_sine"],
    ) / args.secondary_sine_floor
    v += max(
        0.0,
        args.rod_cos_floor - r["minimum_rod_axis_cosine"],
    ) / max(1.0 - args.rod_cos_floor, 1e-3)
    v += max(0.0, r["EH_over_crank"] - args.maximum_EH) / args.maximum_EH
    v += max(
        0.0,
        args.crank_clearance_floor - r["crank_axis_to_EFH_clearance_over_crank"],
    ) / max(args.crank_clearance_floor, 0.1)
    v += max(
        0.0,
        r["H_axis_lateral_rms_over_stroke"] - args.h_lateral_rms_max,
    ) / args.h_lateral_rms_max
    v += max(
        0.0,
        r["H_axis_lateral_span_over_stroke"] - args.h_lateral_span_max,
    ) / args.h_lateral_span_max
    return float(v)


def ensure_option(name: str, value: str):
    if name not in sys.argv:
        sys.argv.extend([name, value])


def selected_output() -> Path:
    if "--output" in sys.argv:
        i = sys.argv.index("--output")
        return Path(sys.argv[i + 1])
    return DEFAULT_OUTPUT


def main():
    # Install corrected objective into the existing tested search driver.
    base.evaluate = evaluate_v41
    base.violation = violation_v41

    base.PRIMARY_LIBRARY = DEFAULT_LIBRARY
    base.OUTPUT = DEFAULT_OUTPUT

    ensure_option("--primary-library", str(DEFAULT_LIBRARY))
    ensure_option("--output", str(DEFAULT_OUTPUT))
    ensure_option("--primary-ranks", "1,0,4,12,9,50")
    ensure_option("--restarts", "3")
    ensure_option("--generations", "360")
    ensure_option("--population-size", "14")
    ensure_option("--velocity-weight", "0.0005")

    base.main()

    # Correct the inherited V1-era metadata and print the extra diagnostics.
    path = selected_output()
    data = json.loads(path.read_text(encoding="utf-8"))
    data["description"] = (
        "3952 LARGE second-dyad synthesis over diversified V4.1 primaries; "
        "full-cycle position fit including HP; bridged-velocity tie-breaker; "
        "wrong-direction/extra-reversal penalty."
    )
    data["objective"] = {
        "position_rms": "complete cycle, HP INCLUDED",
        "velocity": "target velocity linearly bridged through HP band",
        "velocity_weight": 0.0005,
        "monotonicity_weight": MONOTONICITY_WEIGHT,
        "extra_reversal_weight": EXTRA_REVERSAL_WEIGHT,
        "monotonicity_speed_floor_fraction": MONOTONICITY_SPEED_FLOOR,
        "reversal_smoothing_deg": REVERSAL_SMOOTHING_DEG,
        "target_acceleration": "not fitted",
    }
    data["constraints"]["piston_rod_over_stroke_maximum"] = None
    data["constraints"]["piston_rod_over_stroke_note"] = (
        "reported only; no hard ceiling in V4.1 downstream search"
    )
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("\nV4.1 EXTRA DIAGNOSTICS")
    for fam in data.get("families", []):
        b = fam.get("best")
        if not b:
            continue
        print(
            f"rank {fam['primary_rank']:2d}: "
            f"full={100*b['position_rms_full_cycle']:.5f}% "
            f"nonHP={100*b['position_rms_non_hp']:.5f}% "
            f"HP={100*b['position_rms_hp_band']:.3f}% "
            f"wrong={b['wrong_sign_velocity_rms']:.5f} "
            f"zeros={b['zero_crossing_count']} "
            f"rod/stroke={b['piston_rod_over_stroke']:.3f}"
        )


if __name__ == "__main__":
    main()
