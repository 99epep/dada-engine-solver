from pathlib import Path
import os
import tempfile

import pytest

from dada_solver.exchangers.doty import load_references, scale_bank
from dada_solver.exchangers.microtube_geometry import MicrotubeBank


DATA_N2 = (
    Path(__file__).resolve().parents[1]
    / "examples/data/doty_1991_nitrogen_reference.csv"
)
DATA_HE = (
    Path(__file__).resolve().parents[1]
    / "examples/data/doty_1991_helium_reference.csv"
)


def estimate(reference, banks=1, flow=None, **overrides):
    inputs = dict(
        parallel_banks=banks,
        total_mass_flow_kg_s=(
            reference.mass_flow_kg_s if flow is None else flow
        ),
        ua_flow_exponent=0.0,
        ua_condition_factor=1.0,
        viscosity_ratio=1.0,
        density_ratio=1.0,
        additional_pressure_drop_pa=0.0,
    )
    inputs.update(overrides)
    return scale_bank(reference, **inputs)


# ============================================================================
# Nitrogen reference tests
# ============================================================================

def test_measured_anchors_and_parallel_replication():
    for reference in load_references(DATA_N2):
        result = estimate(
            reference,
            banks=4,
            flow=4 * reference.mass_flow_kg_s,
        )
        assert result["estimated_ua_w_k"] == pytest.approx(
            4 * reference.ua_w_k
        )
        assert result["estimated_tube_pressure_drop_pa"] == pytest.approx(
            reference.tube_pressure_drop_pa
        )
        assert result["sum_bank_envelope_m3"] == pytest.approx(0.000798)


def test_density_viscosity_and_flow_scaling():
    ref = load_references(DATA_N2)[0]
    result = estimate(
        ref,
        flow=2 * ref.mass_flow_kg_s,
        viscosity_ratio=1.5,
        density_ratio=3,
        additional_pressure_drop_pa=100,
    )
    assert result[
        "estimated_tube_and_additional_pressure_drop_pa"
    ] == pytest.approx(5200)
    assert ref.tube_mean_temperature_k == pytest.approx(
        (296.15 + 366.15) / 2
    )


@pytest.mark.parametrize(
    "overrides",
    [
        dict(total_mass_flow_kg_s=0),
        dict(total_mass_flow_kg_s=-1),
        dict(density_ratio=float("nan")),
        dict(parallel_banks=1.5),
        dict(parallel_banks=True),
        dict(ua_flow_exponent=2),
    ],
)
def test_unsupported_inputs_are_explicit(overrides):
    with pytest.raises(ValueError):
        estimate(load_references(DATA_N2)[0], **overrides)


# ============================================================================
# Helium reference data
# ============================================================================

def test_helium_file_loads_and_has_three_rows():
    references = load_references(DATA_HE)
    assert len(references) == 3


def test_helium_reference_values():
    """Verify the SI transcription of Doty Tables 1-2."""
    refs = load_references(DATA_HE)

    assert refs[0].mass_flow_kg_s == pytest.approx(0.000117)
    assert refs[0].ua_w_k == pytest.approx(11.1)
    assert refs[0].tube_pressure_drop_pa == pytest.approx(2300)
    assert refs[0].pressure_pa == pytest.approx(749000)
    assert refs[0].tube_mean_temperature_k == pytest.approx(
        (296.65 + 376.75) / 2
    )

    assert refs[1].mass_flow_kg_s == pytest.approx(0.000079)
    assert refs[1].ua_w_k == pytest.approx(6.9)
    assert refs[1].tube_pressure_drop_pa == pytest.approx(1500)
    assert refs[1].pressure_pa == pytest.approx(756000)
    assert refs[1].tube_mean_temperature_k == pytest.approx(
        (297.15 + 373.65) / 2
    )

    assert refs[2].mass_flow_kg_s == pytest.approx(0.000213)
    assert refs[2].ua_w_k == pytest.approx(15.4)
    assert refs[2].tube_pressure_drop_pa == pytest.approx(4400)
    assert refs[2].pressure_pa == pytest.approx(825000)
    assert refs[2].tube_mean_temperature_k == pytest.approx(
        (296.25 + 376.55) / 2
    )


# ============================================================================
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
        raise ValueError(
            "This local NIST interpolation is restricted to 325-350 K."
        )
    fraction = (temperature_k - 325.0) / 25.0
    return (
        _NIST_HE_VISCOSITY_325_PA_S
        + fraction
        * (
            _NIST_HE_VISCOSITY_350_PA_S
            - _NIST_HE_VISCOSITY_325_PA_S
        )
    )


def test_helium_viscosity_reference_points():
    assert _helium_viscosity_nist(325.0) == pytest.approx(21.0e-6)
    assert _helium_viscosity_nist(350.0) == pytest.approx(22.1e-6)


@pytest.mark.parametrize(
    "row_index,doty_eq7_pa,independent_dp_pa",
    [
        (0, 3000.0, 3318.8793),
        (1, 2000.0, 2205.7510),
        (2, 5500.0, 5477.2102),
    ],
)
def test_helium_absolute_poiseuille_against_doty(
    row_index,
    doty_eq7_pa,
    independent_dp_pa,
):
    """Reproduce Doty Eq. 7 independently and compare with helium data."""
    ref = load_references(DATA_HE)[row_index]
    rho = _helium_density_from_ideal_gas(
        ref.pressure_pa,
        ref.tube_mean_temperature_k,
    )
    mu = _helium_viscosity_nist(ref.tube_mean_temperature_k)

    result = _DOTY_BANK.laminar_tube_loss(
        ref.mass_flow_kg_s,
        density_kg_m3=rho,
        viscosity_pa_s=mu,
    )
    dp = result["signed_tube_pressure_drop_pa"]

    # Regression values from the explicit NIST/ideal-gas assumptions above.
    assert dp == pytest.approx(independent_dp_pa, rel=0.01)

    # Doty's paper does not document enough property-reduction detail to
    # require exact reproduction of the Table 2 calculated column. Under the
    # explicit convention above, the reconstruction remains within about 11 %.
    assert dp == pytest.approx(doty_eq7_pa, rel=0.12)

    # All three independent theoretical predictions remain above experiment.
    assert dp > ref.tube_pressure_drop_pa

    # These particular Doty helium points are well inside the laminar regime.
    assert result["reynolds_number"] < 2300


# ============================================================================
# Relative scaling between helium measurements
#
# Use the 117 mg/s measurement as the single measured anchor. Apply the
# existing dp ~ mass_flow * viscosity / density scaling with the same NIST
# viscosity interpolation and ideal-gas representative density.
# ============================================================================

@pytest.mark.parametrize(
    "row_index,expected_pressure_drop_pa,expected_relative_error_percent",
    [
        (1, 1528.5965, 1.9064),
        (2, 3795.7341, -13.7333),
    ],
)
def test_helium_relative_scaling_from_117mg_s_anchor(
    row_index,
    expected_pressure_drop_pa,
    expected_relative_error_percent,
):
    anchor = load_references(DATA_HE)[0]
    target = load_references(DATA_HE)[row_index]

    density_anchor = _helium_density_from_ideal_gas(
        anchor.pressure_pa,
        anchor.tube_mean_temperature_k,
    )
    density_target = _helium_density_from_ideal_gas(
        target.pressure_pa,
        target.tube_mean_temperature_k,
    )
    density_ratio = density_target / density_anchor

    viscosity_anchor = _helium_viscosity_nist(
        anchor.tube_mean_temperature_k
    )
    viscosity_target = _helium_viscosity_nist(
        target.tube_mean_temperature_k
    )
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
    relative_error_percent = (
        (predicted_pa - measured_pa) / measured_pa * 100
    )

    assert predicted_pa == pytest.approx(
        expected_pressure_drop_pa,
        rel=0.01,
    )
    assert relative_error_percent == pytest.approx(
        expected_relative_error_percent,
        abs=0.3,
    )


# ============================================================================
# File format and rejection tests
# ============================================================================

_CSV_HEADER = (
    "fluid,mass_flow_kg_s,mass_flow_uncertainty_kg_s,"
    "pressure_Pa,pressure_uncertainty_Pa,T1_K,T2_K,T3_K,T4_K,"
    "effectiveness,effectiveness_uncertainty,UA_W_K,UA_uncertainty_W_K,"
    "pressure_drop_Pa,pressure_drop_uncertainty_Pa\n"
)


def _write_temporary_reference(*rows):
    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".csv",
        delete=False,
        newline="",
    ) as stream:
        stream.write(_CSV_HEADER)
        for row in rows:
            stream.write(row)
        return stream.name


def test_mixed_fluid_file_rejected():
    path = _write_temporary_reference(
        (
            "nitrogen,0.000470,0.000020,322000,3000,377.85,304.05,"
            "296.15,366.15,0.857,0.006,3.7,0.4,5100,200\n"
        ),
        (
            "helium,0.000117,0.000007,749000,3000,380.65,301.45,"
            "296.65,376.75,0.954,0.006,11.1,2.5,2300,100\n"
        ),
    )
    try:
        with pytest.raises(ValueError, match="homogeneous fluid"):
            load_references(path)
    finally:
        os.unlink(path)


def test_empty_file_rejected():
    path = _write_temporary_reference()
    try:
        with pytest.raises(ValueError, match="non-empty"):
            load_references(path)
    finally:
        os.unlink(path)


def test_unsupported_fluid_rejected():
    path = _write_temporary_reference(
        (
            "argon,0.000470,0.000020,322000,3000,377.85,304.05,"
            "296.15,366.15,0.857,0.006,3.7,0.4,5100,200\n"
        )
    )
    try:
        with pytest.raises(ValueError, match="Unsupported fluid"):
            load_references(path)
    finally:
        os.unlink(path)
