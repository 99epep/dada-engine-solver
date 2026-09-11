#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path

TESTS = Path("tests/test_doty.py")
DOC = Path("docs/DOTY_SCREENING.md")
README = Path("examples/data/README.md")


def replace_between(text: str, start_pattern: str, end_pattern: str, replacement: str, *, label: str) -> str:
    pattern = re.compile(start_pattern + r".*?(?=" + end_pattern + r")", re.MULTILINE | re.DOTALL)
    updated, count = pattern.subn(replacement.rstrip() + "\n\n", text, count=1)
    if count != 1:
        raise RuntimeError(f"Could not uniquely replace {label}; found {count} matches.")
    return updated


def replace_exact(text: str, old: str, new: str, *, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"Could not uniquely replace {label}; found {count} matches.")
    return text.replace(old, new, 1)


def update_tests() -> None:
    text = TESTS.read_text(encoding="utf-8")

    replacement = r'''# ============================================================================
# Hydraulic checks against the helium measurements
#
# MicrotubeBank.laminar_tube_loss() implements the same Poiseuille form as
# Doty Eq. 7:
#
#   dp = 128 * mu * L * mass_flow / (rho * n * pi * d_i**4)
#
# Independent property assumptions used here:
#
# - representative tube temperature = arithmetic mean of measured T3 and T4;
# - helium density = ideal-gas density at the reported pressure and that
#   representative temperature;
# - dilute-gas helium viscosity = linear interpolation of the NIST values
#   21.0 uPa.s at 325 K and 22.1 uPa.s at 350 K.
#
# NIST source:
# https://www.nist.gov/pml/sensor-science/fluid-metrology/database-thermophysical-properties-gases-used-semiconductor-9
#
# The NIST table attributes the helium viscosity data to Hurly and Moldover
# (2000) and reports an estimated uncertainty of 0.1 %. These tests do not
# assume that Doty used exactly the same property-reduction convention.
# ============================================================================

_DOTY_BANK = MicrotubeBank(
    tube_count=309,
    tube_length_m=0.127,
    inner_diameter_m=0.00033,
    wall_thickness_m=0.0001524,
    pitch_m=0.00125,
    # Header depth is not used by laminar_tube_loss(); a positive value is
    # required only to construct the generic MicrotubeBank geometry object.
    header_depth_m=0.001,
)

_R_HE_J_KG_K = 2077.1
_NIST_HE_VISCOSITY_325_PA_S = 21.0e-6
_NIST_HE_VISCOSITY_350_PA_S = 22.1e-6


def _helium_density_from_ideal_gas(pressure_pa, temperature_k):
    """Ideal-gas helium density at the selected representative condition."""
    return pressure_pa / (_R_HE_J_KG_K * temperature_k)


def _helium_viscosity_nist(temperature_k):
    """Interpolate NIST dilute-gas helium viscosity between 325 and 350 K."""
    if not 325.0 <= temperature_k <= 350.0:
        raise ValueError("This local NIST interpolation is restricted to 325-350 K.")
    fraction = (temperature_k - 325.0) / 25.0
    return (
        _NIST_HE_VISCOSITY_325_PA_S
        + fraction
        * (_NIST_HE_VISCOSITY_350_PA_S - _NIST_HE_VISCOSITY_325_PA_S)
    )


def test_helium_viscosity_reference_points():
    assert _helium_viscosity_nist(325.0) == pytest.approx(21.0e-6)
    assert _helium_viscosity_nist(350.0) == pytest.approx(22.1e-6)


@pytest.mark.parametrize(
    "row_index,doty_eq7_pa,independent_dp_pa",
    [
        (0, 3000.0, 3319.0),
        (1, 2000.0, 2206.0),
        (2, 5500.0, 5477.0),
    ],
)
def test_helium_absolute_poiseuille_against_doty(
    row_index, doty_eq7_pa, independent_dp_pa
):
    """Reproduce Doty Eq. 7 independently and compare with measured helium dp."""
    ref = load_references(DATA_HE)[row_index]
    rho = _helium_density_from_ideal_gas(
        ref.pressure_pa, ref.tube_mean_temperature_k
    )
    mu = _helium_viscosity_nist(ref.tube_mean_temperature_k)

    result = _DOTY_BANK.laminar_tube_loss(
        ref.mass_flow_kg_s,
        density_kg_m3=rho,
        viscosity_pa_s=mu,
    )
    dp = result["signed_tube_pressure_drop_pa"]

    assert dp == pytest.approx(independent_dp_pa, rel=0.01)
    assert dp == pytest.approx(doty_eq7_pa, rel=0.12)
    assert dp > ref.tube_pressure_drop_pa
    assert result["reynolds_number"] < 2300


# ============================================================================
# Hydraulic Check 2: relative scaling between helium measurements
#
# Use the 117 mg/s measurement as the single measured anchor. Apply the
# existing dp ~ mass_flow * viscosity / density scaling with the same NIST
# viscosity interpolation and ideal-gas representative density.
# ============================================================================

@pytest.mark.parametrize(
    "row_index,expected_pressure_drop_pa,expected_relative_error_percent",
    [
        (1, 1529.0, 1.9),
        (2, 3796.0, -13.7),
    ],
)
def test_helium_relative_scaling_from_117mg_s_anchor(
    row_index, expected_pressure_drop_pa, expected_relative_error_percent
):
    anchor = load_references(DATA_HE)[0]
    target = load_references(DATA_HE)[row_index]

    density_anchor = _helium_density_from_ideal_gas(
        anchor.pressure_pa, anchor.tube_mean_temperature_k
    )
    density_target = _helium_density_from_ideal_gas(
        target.pressure_pa, target.tube_mean_temperature_k
    )
    density_ratio = density_target / density_anchor

    viscosity_anchor = _helium_viscosity_nist(anchor.tube_mean_temperature_k)
    viscosity_target = _helium_viscosity_nist(target.tube_mean_temperature_k)
    viscosity_ratio = viscosity_target / viscosity_anchor

    result = estimate(
        anchor,
        banks=1,
        flow=target.mass_flow_kg_s,
        viscosity_ratio=viscosity_ratio,
        density_ratio=density_ratio,
    )
    predicted_pa = result["estimated_tube_pressure_drop_pa"]
    measured_pa = target.tube_pressure_drop_pa
    relative_error_percent = (predicted_pa - measured_pa) / measured_pa * 100

    assert predicted_pa == pytest.approx(expected_pressure_drop_pa, rel=0.01)
    assert relative_error_percent == pytest.approx(
        expected_relative_error_percent, abs=0.3
    )
'''

    start = r"^# =+\n# Hydraulic Check 1:.*?\n"
    end = r"^# =+\n# File format and rejection tests"
    text = replace_between(text, start, end, replacement, label="helium hydraulic test sections")
    TESTS.write_text(text, encoding="utf-8")


def update_doc() -> None:
    text = DOC.read_text(encoding="utf-8")

    replacement = r'''### Absolute Poiseuille comparison (Doty Eq. 7)

Doty's Table 2 gives both measured helium tube-side pressure drops and values
calculated from Eq. 7:

| Flow (mg/s) | Measured dP (Pa) | Doty Eq. 7 (Pa) | Doty overprediction |
|---|---:|---:|---:|
| 117 | 2300 | 3000 | +30.4 % |
| 79 | 1500 | 2000 | +33.3 % |
| 213 | 4400 | 5500 | +25.0 % |

Doty explicitly states that the helium pressure drops are somewhat lower than
expected and notes minor tube-side flow maldistribution as a possible
explanation.

The project also performs an independent reproduction of Eq. 7 rather than
using the Table 2 calculated column as an input. The assumptions are:

- 309 tubes, 0.33 mm internal diameter and 127 mm length;
- representative tube temperature `(T3 + T4) / 2`;
- ideal-gas helium density at the reported pressure;
- dilute-gas helium viscosity linearly interpolated from the NIST values
  21.0 uPa.s at 325 K and 22.1 uPa.s at 350 K.

NIST source:
https://www.nist.gov/pml/sensor-science/fluid-metrology/database-thermophysical-properties-gases-used-semiconductor-9

The NIST table attributes these viscosity values to Hurly and Moldover (2000)
and gives an estimated viscosity uncertainty of 0.1 %. Doty's paper does not
state enough detail about the density/viscosity evaluation used for Table 2 to
require exact numerical reproduction of its calculated column.

| Flow (mg/s) | Measured dP (Pa) | Doty Eq. 7 (Pa) | Independent Eq. 7 (Pa) | Independent vs measured | Independent vs Doty |
|---|---:|---:|---:|---:|---:|
| 117 | 2300 | 3000 | 3319 | +44.3 % | +10.6 % |
| 79 | 1500 | 2000 | 2206 | +47.1 % | +10.3 % |
| 213 | 4400 | 5500 | 5477 | +24.5 % | -0.4 % |

Thus `MicrotubeBank.laminar_tube_loss()` is consistent with the algebraic form
of Doty's Eq. 7 and reproduces the published calculated values to about 11 %
under an explicit independent property convention. It does **not** validate
the absolute helium pressure drop experimentally. No empirical multiplier is
introduced to force agreement.

### Relative scaling between helium measurements

Using the 117 mg/s measurement as a single anchor and applying
`dp ~ mass_flow * viscosity / density`, with the same NIST viscosity
interpolation and ideal-gas representative density:

| Flow (mg/s) | Measured dP (Pa) | Predicted from 117 mg/s anchor (Pa) | Relative error |
|---|---:|---:|---:|
| 117 | 2300 | 2300 (anchor) | - |
| 79 | 1500 | 1529 | +1.9 % |
| 213 | 4400 | 3796 | -13.7 % |

The relative scaling is therefore substantially more successful than the
absolute prediction: it tracks the 79 mg/s point closely and remains within
about 14 % at 213 mg/s. This supports use of the laminar scaling as screening
physics while retaining an explicit uncertainty on absolute loss.'''

    start = r"^### Absolute Poiseuille comparison \(Doty Eq\. 7\)\n"
    end = r"^## Reproducible illustrative comparison"
    text = replace_between(text, start, end, replacement, label="DOTY_SCREENING helium validation section")
    DOC.write_text(text, encoding="utf-8")


def update_readme() -> None:
    text = README.read_text(encoding="utf-8")

    new = """**Important:** Doty's own Eq. 7 values overpredict the measured helium pressure
drops by about 25–33 %. An independent project reconstruction using NIST
dilute-gas helium viscosity, ideal-gas density and mean tube temperature gives
approximately 3.319, 2.206 and 5.477 kPa, versus Doty's published 3.0, 2.0 and
5.5 kPa and measured 2.3, 1.5 and 4.4 kPa. The independent reproduction is
within about 11 % of Doty's Eq. 7 values but remains clearly above experiment.
Relative scaling from a measured helium anchor is much closer. No hidden
calibration factor is introduced. See `docs/DOTY_SCREENING.md` for assumptions
and detailed comparisons.
"""

    pattern = re.compile(
        r"\*\*Important:.*?See `docs/DOTY_SCREENING\.md` for detailed comparisons\.\n",
        re.DOTALL,
    )
    text, count = pattern.subn(new, text, count=1)
    if count != 1:
        raise RuntimeError(
            f"Could not uniquely replace helium README interpretation; found {count} matches."
        )
    README.write_text(text, encoding="utf-8")


def main() -> None:
    missing = [str(path) for path in (TESTS, DOC, README) if not path.exists()]
    if missing:
        raise SystemExit(
            "Run this script from the repository root. Missing: " + ", ".join(missing)
        )

    update_tests()
    update_doc()
    update_readme()

    print("Updated:")
    print(f"  {TESTS}")
    print(f"  {DOC}")
    print(f"  {README}")
    print()
    print("No physical model code was changed.")
    print("Now run:")
    print("  PYTHONPATH=src python3 -m pytest tests/test_doty.py")
    print("  PYTHONPATH=src python3 -m pytest")


if __name__ == "__main__":
    main()
